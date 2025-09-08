# rigidity analysis of embedded simplicial surfaces

"""
    positions!(mesh::PolyhedralMesh{PositionDim,PositionType}, coordinate_matrix::AbstractMatrix{PositionType}) where {PositionDim,PositionType}


Set the vertex positions of `mesh` to the columns of the matrix `coordinate_matrix`
"""
function positions!(mesh::PolyhedralMesh{PositionDim,PositionType}, coordinate_matrix::AbstractMatrix{T}) where {PositionDim,PositionType,T<:Number}
    for v in vertices(mesh)
        position!(v, Point(coordinate_matrix[:, id(v)]...))
    end

    return mesh
end
positions!(surf::SimplicialSurface, coordinate_matrix::AbstractMatrix{T}) where T<:Number = positions!(surf.mesh, coordinate_matrix)

"""
    coord(mesh::PolyhedralMesh)

Return the coordinate matrix of `mesh`. The coordinates of the vertices are the columns of the returned matrix.
"""
function positions(mesh::PolyhedralMesh)
    return hcat([Vector(position(v)) for v in sort(vertices(mesh), by=id)]...)
end
positions(surf::SimplicialSurface) = positions(surf.mesh)

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

RigidityTheoryTools.Framework(mesh::PolyhedralMesh) = RigidityTheoryTools.Framework(Graph(mesh), positions(mesh), labels(mesh))