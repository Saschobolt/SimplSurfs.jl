octahedron() = PolyhedralMesh([[1, 2, 3], [1, 3, 4], [1, 4, 5], [1, 2, 5], [6, 2, 3], [6, 3, 4], [6, 4, 5], [6, 2, 5]])

octahedron_emb() = begin
    m = PolyhedralMesh{3,Float64}([[1, 2, 3], [1, 3, 4], [1, 4, 5], [1, 2, 5], [6, 2, 3], [6, 3, 4], [6, 4, 5], [6, 2, 5]])
    coordinates!(m, transpose([0 0 -1; 1 1 0; -1 1 0; -1 -1 0; -1 1 0; 0 0 1]))
    return m
end

tetrahedron() = PolyhedralMesh([[1, 2, 3], [1, 2, 4], [1, 3, 4], [2, 3, 4]])

tetrahedron_emb() = begin
    m = PolyhedralMesh{3,Float64}([[1, 2, 3], [1, 2, 4], [1, 3, 4], [2, 3, 4]])
    coordinates!(m, transpose([0 0 0; 1 0 0; 0 1 0; 0 0 1]))
    return m
end

double_tetrahedron() = PolyhedralMesh([[1, 2, 3], [1, 2, 4], [1, 3, 4], [5, 2, 3], [5, 2, 4], [5, 3, 4]])

cube() = begin
    m = PolyhedralMesh{3,Float64}([[1, 2, 3, 4], [1, 2, 6, 5], [2, 3, 7, 6], [3, 4, 8, 7], [1, 4, 8, 5], [5, 6, 7, 8]])
    coordinates!(m, transpose([0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1]))
    return m
end