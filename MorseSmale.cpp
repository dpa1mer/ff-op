#include "mex.hpp"
#include "mexAdapter.hpp"

#include "PersistenceDiagram.h"
#include "TopologicalSimplification.h"
#include "MorseSmaleComplex.h"

#define ML_ERR(errMsg) matlabPtr->feval(u"error", 0, std::vector<matlab::data::Array>({ factory.createScalar(errMsg) }))
#define ML_PRINT(msg) matlabPtr->feval(u"fprintf", 0, std::vector<matlab::data::Array>({ factory.createScalar(msg) }))

using MLArray = matlab::data::TypedArray<double>;
using MLPtr = matlab::data::buffer_ptr_t<double>;
using MLArrayType = matlab::data::ArrayType;
using MLFactory = matlab::data::ArrayFactory;
using MLEngine = matlab::engine::MATLABEngine;

class MexFunction : public matlab::mex::Function {
    MLFactory factory;

    class mystream : public std::streambuf {
    public:
        mystream(std::shared_ptr<MLEngine> & matlabPtr, MLFactory & factory)
        : std::streambuf(), matlabPtr(matlabPtr), factory(factory) {}
    protected:
        virtual std::streamsize xsputn(const char *s, std::streamsize n) { ML_PRINT(std::string(s, n)); return n; }
        virtual int overflow(int c=EOF) { if (c != EOF) { ML_PRINT(std::string(1, c)); } return 1; }
    private:
        std::shared_ptr<MLEngine> matlabPtr;
        MLFactory & factory;
    };
    class scoped_redirect_cout {
    public:
        scoped_redirect_cout(std::shared_ptr<MLEngine> & matlabPtr, MLFactory & factory)
        : mout(matlabPtr, factory) { 
            old_cout = std::cout.rdbuf(); std::cout.rdbuf(&mout);
            old_cerr = std::cerr.rdbuf(); std::cerr.rdbuf(&mout);
        }
        ~scoped_redirect_cout() { std::cout.rdbuf(old_cout); std::cerr.rdbuf(old_cerr); }
    private:
        mystream mout;
        std::streambuf *old_cout;
        std::streambuf *old_cerr;
    };

public:
    void checkArguments(matlab::mex::ArgumentList inputs) {
        std::shared_ptr<MLEngine> matlabPtr;

        if (inputs.size() < 4) ML_ERR("Must provide four inputs.");
        
        if (inputs[0].getType() != MLArrayType::DOUBLE
         || inputs[1].getType() != MLArrayType::DOUBLE
         || inputs[2].getType() != MLArrayType::DOUBLE
         || inputs[3].getType() != MLArrayType::DOUBLE) {
            ML_ERR("Inputs must be of type double.");
        }

        if (inputs[3].getNumberOfElements() != 1) ML_ERR("Minimum persistence must be a scalar.");
    }

    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
        std::shared_ptr<MLEngine> matlabPtr = getEngine();
        scoped_redirect_cout mycout_redirect(matlabPtr, factory);

        ttk::globalDebugLevel_ = 2;

        MLArray simplices = std::move(inputs[0]);
        MLArray verts = std::move(inputs[1]);
        MLArray scalar = std::move(inputs[2]);

        double minPersistence = inputs[3][0];

        std::vector<size_t> dimVerts = verts.getDimensions();
        if (dimVerts.size() != 2) ML_ERR("Inputs must be matrices.");
        size_t D = dimVerts[0];
        size_t nv = dimVerts[1];
        if (D != 3) ML_ERR("Ambient dimension must be 3.");

        std::vector<size_t> dimSimplices = simplices.getDimensions();
        if (dimSimplices.size() != 2) ML_ERR("Inputs must be matrices.");
        size_t d = dimSimplices[1] - 1;
        if (d != 2 && d != 3) ML_ERR("Simplices must be triangles or tetrahedra.");
        size_t nc = dimSimplices[0];

        // order array: every vertex sorted according to the scalar field
        std::vector<ttk::SimplexId> perm(nv);
        std::iota(perm.begin(), perm.end(), 0);
        std::stable_sort(perm.begin(), perm.end(),
            [&scalar] (size_t i1, size_t i2) {return scalar[i1] < scalar[i2];});

        std::vector<ttk::SimplexId> order(nv);
        for (size_t i = 0; i < nv; ++i) {
            order[perm[i]] = i;
        }

        MLPtr vertsPtr = verts.release();
        MLPtr scalarPtr = scalar.release();

        // Copy the data into TTK's format
        std::ostringstream stream;
        stream << "Setting up TTK triangulation: " << nv << " verts, " << nc << " cells...\n";
        ML_PRINT(stream.str());
        ttk::Triangulation triangulation;
        triangulation.setInputPoints(nv, vertsPtr.get(), true);

        std::vector<ttk::LongSimplexId> inputCells;
        for (int i = 0; i < nc; ++i) {
            inputCells.push_back(d + 1); // Size of each simplex
            for (int j = 0; j < d + 1; ++j) {
                // Convert MATLAB 1-offset to 0-offset
                inputCells.push_back(static_cast<ttk::LongSimplexId>(simplices[i][j]) - 1);
            }
        }
        triangulation.setInputCells(nc, inputCells.data());


        // Topological simplification
        // Following example code from https://github.com/topology-tool-kit/ttk/blob/dev/examples/c%2B%2B/main.cpp
        ML_PRINT("Simplifying Morse function...\n");
        using PersistencePair = std::tuple<ttk::SimplexId, ttk::CriticalType, ttk::SimplexId, ttk::CriticalType, double, ttk::SimplexId>;
        std::vector<PersistencePair> diagramOutput;
        std::vector< std::tuple<ttk::dcg::Cell, ttk::dcg::Cell> > dmtPairs;
        ttk::PersistenceDiagram diagram;
        diagram.setupTriangulation(&triangulation);
        diagram.setInputScalars(scalarPtr.get());
        diagram.setInputOffsets(order.data());
        diagram.setDMTPairs(&dmtPairs);
        diagram.setOutputCTDiagram(&diagramOutput);
        diagram.execute<double, ttk::SimplexId>();

        // Select critical point pairs with sufficient persistence
        std::vector<ttk::SimplexId> selectedCriticalPoints;
        for (int i = 0; i < diagramOutput.size(); ++i) {
            if (std::get<4>(diagramOutput[i]) > minPersistence) {
                selectedCriticalPoints.push_back(std::get<0>(diagramOutput[i]));
                selectedCriticalPoints.push_back(std::get<2>(diagramOutput[i]));
            }
        }

        // Modify the scalar function to eliminate spurious critical points
        std::vector<double> simplifiedScalar(nv);
        std::vector<ttk::SimplexId> simplifiedOrder = order;
        ttk::TopologicalSimplification simplification;
        simplification.setupTriangulation(&triangulation);
        simplification.setVertexNumber(nv);
        simplification.setConstraintNumber(selectedCriticalPoints.size());
        simplification.setVertexIdentifierScalarFieldPointer(selectedCriticalPoints.data());
        simplification.setInputScalarFieldPointer(scalarPtr.get());
        simplification.setOutputScalarFieldPointer(simplifiedScalar.data());
        simplification.setInputOffsetScalarFieldPointer(order.data());
        simplification.setOutputOffsetScalarFieldPointer(simplifiedOrder.data());
        simplification.execute<double, ttk::SimplexId>();


        // Compute the Morse-Smale Complex
        // Following example code from https://github.com/topology-tool-kit/ttk/blob/dev/examples/c%2B%2B/main.cpp
        ML_PRINT("Computing MSC...\n");
        ttk::MorseSmaleComplex morseSmaleComplex;

        // critical points
        ttk::SimplexId criticalPoints_numberOfPoints{};
        std::vector<float> criticalPoints_points;
        std::vector<char> criticalPoints_points_cellDimensions;
        std::vector<ttk::SimplexId> criticalPoints_points_cellIds;
        std::vector<char> criticalPoints_points_isOnBoundary;
        std::vector<float> criticalPoints_points_cellScalars;
        std::vector<ttk::SimplexId> criticalPoints_points_PLVertexIdentifiers;
        std::vector<ttk::SimplexId> criticalPoints_points_manifoldSize;

        // 1-separatrices
        ttk::SimplexId separatrices1_numberOfPoints{};
        std::vector<float> separatrices1_points;
        std::vector<char> separatrices1_points_smoothingMask;
        std::vector<char> separatrices1_points_cellDimensions;
        std::vector<ttk::SimplexId> separatrices1_points_cellIds;
        ttk::SimplexId separatrices1_numberOfCells{};
        std::vector<ttk::SimplexId> separatrices1_cells_connectivity;
        std::vector<ttk::SimplexId> separatrices1_cells_sourceIds;
        std::vector<ttk::SimplexId> separatrices1_cells_destinationIds;
        std::vector<ttk::SimplexId> separatrices1_cells_separatrixIds;
        std::vector<char> separatrices1_cells_separatrixTypes;
        std::vector<char> separatrices1_cells_isOnBoundary;
        std::vector<double> separatrices1_cells_separatrixFunctionMaxima;
        std::vector<double> separatrices1_cells_separatrixFunctionMinima;
        std::vector<double> separatrices1_cells_separatrixFunctionDiffs;

        // segmentation
        std::vector<ttk::SimplexId> ascendingSegmentation(nv, -1),
                                    descendingSegmentation(nv, -1),
                                    mscSegmentation(nv, -1);
        
        morseSmaleComplex.setupTriangulation(&triangulation);
        morseSmaleComplex.setInputScalarField(simplifiedScalar.data());
        morseSmaleComplex.setInputOffsets(simplifiedOrder.data());
        morseSmaleComplex.setIterationThreshold(1e8);
        
        morseSmaleComplex.setOutputMorseComplexes(
            ascendingSegmentation.data(),
            descendingSegmentation.data(),
            mscSegmentation.data());
        morseSmaleComplex.setOutputCriticalPoints(
            &criticalPoints_numberOfPoints, &criticalPoints_points,
            &criticalPoints_points_cellDimensions, &criticalPoints_points_cellIds,
            &criticalPoints_points_cellScalars, &criticalPoints_points_isOnBoundary,
            &criticalPoints_points_PLVertexIdentifiers,
            &criticalPoints_points_manifoldSize);
        morseSmaleComplex.setOutputSeparatrices1(
            &separatrices1_numberOfPoints, &separatrices1_points,
            &separatrices1_points_smoothingMask, &separatrices1_points_cellDimensions,
            &separatrices1_points_cellIds, &separatrices1_numberOfCells,
            &separatrices1_cells_connectivity, &separatrices1_cells_sourceIds,
            &separatrices1_cells_destinationIds, &separatrices1_cells_separatrixIds,
            &separatrices1_cells_separatrixTypes,
            &separatrices1_cells_separatrixFunctionMaxima,
            &separatrices1_cells_separatrixFunctionMinima,
            &separatrices1_cells_separatrixFunctionDiffs,
            &separatrices1_cells_isOnBoundary);


        morseSmaleComplex.execute<double, ttk::SimplexId>();

        // Pass results back to MATLAB
        MLArray ascending = factory.createArray<double>({nv, 1});
        MLArray descending = factory.createArray<double>({nv, 1});
        MLArray msc = factory.createArray<double>({nv, 1});
        for (int i = 0; i < nv; ++i) {
            ascending[i] = ascendingSegmentation[i];
            descending[i] = descendingSegmentation[i];
            msc[i] = mscSegmentation[i];
        }

        size_t nCrit = static_cast<size_t>(criticalPoints_numberOfPoints);
        MLArray critDims = factory.createArray<double>({nCrit, 1});
        MLArray critical = factory.createArray<double>({nCrit, 3});
        for (int i = 0; i < nCrit; ++i) {
            critDims[i] = static_cast<double>(criticalPoints_points_cellDimensions[i]);
            for (int j = 0; j < 3; ++j) {
                critical[i][j] = static_cast<double>(criticalPoints_points[3 * i + j]);
            }
        }

        size_t nSepPts = static_cast<size_t>(separatrices1_numberOfPoints);
        MLArray sepPts = factory.createArray<double>({nSepPts, 3});
        for (int i = 0; i < nSepPts; ++i) {
            for (int j = 0; j < 3; ++j) {
                sepPts[i][j] = static_cast<double>(separatrices1_points[3 * i + j]);
            }
        }

        size_t nSepCells = static_cast<size_t>(separatrices1_numberOfCells);
        MLArray sepCells = factory.createArray<double>({nSepCells, 2});
        for (int i = 0; i < nSepCells; ++i) {
            for (int j = 0; j < 2; ++j) {
                sepCells[i][j] = static_cast<double>(separatrices1_cells_connectivity[3 * i + j + 1] + 1);
            }
        }


        outputs[0] = ascending;
        outputs[1] = descending;
        outputs[2] = msc;
        outputs[3] = critical;
        outputs[4] = critDims;
        outputs[5] = sepPts;
        outputs[6] = sepCells;
    }
};