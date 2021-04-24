#include "mex.hpp"
#include "mexAdapter.hpp"

#define TTK_CELL_ARRAY_NEW
#include "PersistenceDiagram.h"
#include "TopologicalSimplification.h"
#include "MorseSmaleComplex.h"
#include "MorseSmaleQuadrangulation.h"

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

        if (inputs.size() < 5) ML_ERR("Must provide five inputs.");
        
        if (inputs[0].getType() != MLArrayType::DOUBLE
         || inputs[1].getType() != MLArrayType::DOUBLE
         || inputs[2].getType() != MLArrayType::DOUBLE
         || inputs[3].getType() != MLArrayType::DOUBLE) {
            ML_ERR("Inputs must be of type double.");
        }

        if (inputs[3].getNumberOfElements() != 1) ML_ERR("Minimum persistence must be a scalar.");
        if (inputs[4].getNumberOfElements() != 1) ML_ERR("Dualization must be a scalar boolean.");
    }

    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
        std::shared_ptr<MLEngine> matlabPtr = getEngine();
        scoped_redirect_cout mycout_redirect(matlabPtr, factory);

        ttk::globalDebugLevel_ = 2;

        MLArray simplices = std::move(inputs[0]);
        MLArray verts = std::move(inputs[1]);
        MLArray scalar = std::move(inputs[2]);

        double minPersistence = inputs[3][0];
        bool dualize = inputs[4][0];

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

        MLPtr vertsPtr = verts.release();
        MLPtr scalarPtr = scalar.release();
        
        // order array: every vertex sorted according to the scalar field
        std::vector<ttk::SimplexId> order(nv);
        ttk::preconditionOrderArray(nv, scalarPtr.get(), order.data());

        // Copy the data into TTK's format
        std::ostringstream stream;
        stream << "Setting up TTK triangulation: " << nv << " verts, " << nc << " cells...\n";
        ML_PRINT(stream.str());
        ttk::Triangulation triangulation;
        triangulation.setInputPoints(nv, vertsPtr.get(), true);

        std::vector<ttk::LongSimplexId> inputCellsConnectivity, inputCellsOffset;
        inputCellsOffset.push_back(0);
        for (int i = 0; i < nc; ++i) {
            for (int j = 0; j < d + 1; ++j) {
                // Convert MATLAB 1-offset to 0-offset
                inputCellsConnectivity.push_back(static_cast<ttk::LongSimplexId>(simplices[i][j]) - 1);
            }
            inputCellsOffset.push_back(inputCellsConnectivity.size());
        }
        triangulation.setInputCells(nc, inputCellsConnectivity.data(), inputCellsOffset.data());


//         // Topological simplification
//         // Following example code from https://github.com/topology-tool-kit/ttk/blob/dev/examples/c%2B%2B/main.cpp
//         ML_PRINT("Simplifying Morse function...\n");
//         std::vector<ttk::PersistencePair> diagramOutput;
//         ttk::PersistenceDiagram diagram;
//         diagram.preconditionTriangulation(&triangulation);
//         diagram.execute(diagramOutput, scalarPtr.get(), order.data(), &triangulation);
// 
//         // Select critical point pairs with sufficient persistence
//         std::vector<ttk::SimplexId> selectedCriticalPoints;
//         for (int i = 0; i < diagramOutput.size(); ++i) {
//             if (diagramOutput[i].persistence > minPersistence) {
//                 selectedCriticalPoints.push_back(diagramOutput[i].birth);
//                 selectedCriticalPoints.push_back(diagramOutput[i].death);
//             }
//         }
// 
//         // Modify the scalar function to eliminate spurious critical points
//         std::vector<double> simplifiedScalar(nv);
//         std::vector<ttk::SimplexId> simplifiedOrder = order;
//         ttk::TopologicalSimplification simplification;
//         simplification.preconditionTriangulation(&triangulation);
//         simplification.execute(
//                 scalarPtr.get(), simplifiedScalar.data(),
//                 selectedCriticalPoints.data(),
//                 order.data(), simplifiedOrder.data(),
//                 selectedCriticalPoints.size(), triangulation);


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
        std::vector<long long> separatrices1_cells_connectivity;
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
        
        morseSmaleComplex.preconditionTriangulation(&triangulation);
        morseSmaleComplex.setInputScalarField(scalarPtr.get());//simplifiedScalar.data());
        morseSmaleComplex.setInputOffsets(order.data());//simplifiedOrder.data());
//         morseSmaleComplex.setIterationThreshold(1e8);
        
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


        morseSmaleComplex.execute<double>(triangulation);

//         // Morse-Smale Quadrangulation
//         ML_PRINT("Extracting quadrangulation...\n");
//         ttk::MorseSmaleQuadrangulation quadrangulator;
//         quadrangulator.preconditionTriangulation(&triangulation);
//         quadrangulator.setCriticalPoints(
//             criticalPoints_numberOfPoints,
//             criticalPoints_points.data(),
//             criticalPoints_points_PLVertexIdentifiers.data(),
//             criticalPoints_points_cellIds.data(),
//             criticalPoints_points_cellDimensions.data());
//         quadrangulator.setSeparatrices(
//             separatrices1_numberOfPoints,
//             separatrices1_points_cellIds.data(),
//             separatrices1_points_cellDimensions.data(),
//             separatrices1_points_smoothingMask.data(),
//             separatrices1_points.data());
//         quadrangulator.setDualQuadrangulation(dualize);
//         quadrangulator.execute(triangulation);

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
        MLArray sepMask = factory.createArray<double>({nSepPts, 1});
        for (int i = 0; i < nSepPts; ++i) {
            sepMask[i] = static_cast<double>(separatrices1_points_smoothingMask[i]);
            for (int j = 0; j < 3; ++j) {
                sepPts[i][j] = static_cast<double>(separatrices1_points[3 * i + j]);
            }
        }

        size_t nSepCells = static_cast<size_t>(separatrices1_numberOfCells);
        MLArray sepCells = factory.createArray<double>({nSepCells, 2});
        for (int i = 0; i < nSepCells; ++i) {
            for (int j = 0; j < 2; ++j) {
                sepCells[i][j] = static_cast<double>(separatrices1_cells_connectivity[2 * i + j] + 1);
            }
        }


        outputs[0] = ascending;
        outputs[1] = descending;
        outputs[2] = msc;
        outputs[3] = critical;
        outputs[4] = critDims;
        outputs[5] = sepPts;
        outputs[6] = sepCells;
        outputs[7] = sepMask;

//         size_t nqv = quadrangulator.outputPoints_.size() / 3;
//         size_t nq = quadrangulator.outputCells_.size() / 5;
//         MLArray quadVerts = factory.createArray<double>({nqv, 3});
//         MLArray quadCells = factory.createArray<double>({nq, 4});
//         for (int i = 0; i < nqv; ++i) {
//             for (int j = 0; j < 3; ++j) {
//                 quadVerts[i][j] = static_cast<double>(quadrangulator.outputPoints_[3 * i + j]);
//             }
//         }
// 
//         for (int i = 0; i < nq; ++i) {
//             for (int j = 0; j < 4; ++j) {
//                 quadCells[i][j] = static_cast<double>(quadrangulator.outputCells_[5 * i + j + 1] + 1);
//             }
//         }
//         outputs[0] = quadVerts;
//         outputs[1] = quadCells;
    }
};