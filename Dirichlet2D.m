function u = Dirichlet2D(meshData, Op, bc)

OpII = Op(meshData.intIdx, meshData.intIdx);
OpIB = Op(meshData.intIdx, meshData.bdryIdx);

u = zeros(meshData.nv, 1);
u(meshData.bdryIdx) = bc;
u(meshData.intIdx) = -OpII \ (OpIB * bc);

end