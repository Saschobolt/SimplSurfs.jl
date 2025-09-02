# rigidity analysis of embedded simplicial surfaces

"""
    coord!(mesh::PolyhedralMesh, coordinate_matrix::AbstractMatrix{T}) where T<:Number

Set the vertex positions of `mesh` to the columns of the matrix `coordinate_matrix`
"""
function coord!(mesh::PolyhedralMesh, coordinate_matrix::AbstractMatrix{T}) where T<:Number
    for v in vertices(mesh)
        coord!(v, coordinate_matrix[:, id(v)])
    end

    return mesh
end
coord!(surf::SimplicialSurface, coordinate_matrix::AbstractMatrix{T}) where T<:Number = coord!(surf.mesh)

"""
    coord(mesh::PolyhedralMesh)

Return the coordinate matrix of `mesh`. The coordinates of the vertices are the columns of the returned matrix.
"""
function coord(mesh::PolyhedralMesh)
    return hcat([coord(v) for v in sort(vertices(mesh), by=id)]...)
end
coord(surf::SimplicialSurface) = coord(surf.mesh)

function Graphs.Graph(mesh::PolyhedralMesh)
    g = Graph()

    for v in vertices(mesh)
        add_vertex!(g)
    end

    for (edge_id, e) in edges(mesh)
        add_edge!(g, edge_id[1], edge_id[2])
    end

    return g
end
Graphs.Graph(surf::SimplicialSurface) = Graph(surf.mesh)

RigidityTheoryTools.Framework(mesh::PolyhedralMesh) = RigidityTheoryTools.Framework(Graph(mesh), coord(mesh), labels(mesh))