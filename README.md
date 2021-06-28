# ff-op: Frame Field Operators
This is the code for the paper
> Palmer, David, Oded Stein, and Justin Solomon. "Frame Field Operators." SGP 2021, online.

# Dependencies
This code depends on
- [`gptoolbox`](https://github.com/alecjacobson/gptoolbox) and [`distmesh`](http://persson.berkeley.edu/distmesh/) for triangle mesh processing
- [`ARFF`](https://github.com/dpa1mer/arff) for 3D frame field manipulation

# Usage
Generally, you should first load a mesh data structure
```MATLAB
% In 2D using gptoolbox:
[verts, faces] = load_mesh('path/to/mesh.obj');
meshData = ProcessMesh2D(verts, faces);

% In 3D using ARFF:
meshData = ImportMesh('path/to/mesh.mesh');
```
before generating a frame field and its operator:
```MATLAB
% In 2D:
z4 = MBO2D(meshData);
[~, Tij] = Frame2Tensor2D(meshData, z4, ellipticity);
[Op, M] = FFOp2D(meshData, Tij, neumann);

% In 3D:
q = MBO(meshData, RayMBO);
octa = OctaMBO;
q = octa.proj(q);

[Op, M] = FFOp3D(meshData, q, ellipticity, neumann);

```
Various figure-generating experiments from the paper can be found in `experiments/`.
