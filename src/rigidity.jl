# rigidity analysis of embedded simplicial surfaces

function coordinates!(mesh::PolyhedralMesh{PositionDim,PositionType}, positions::AbstractVector{<:Point{PositionDim}}) where {PositionDim,PositionType}
    for (i, v) in enumerate(vertices(mesh))
        coordinates!(v, positions[i])
    end

    return mesh
end

"""
    coordinates!(mesh::PolyhedralMesh{PositionDim,PositionType}, coordinate_matrix::AbstractMatrix{PositionType}) where {PositionDim,PositionType}


Set the vertex positions of `mesh` to the columns of the matrix `coordinate_matrix`
"""
function coordinates!(mesh::PolyhedralMesh{PositionDim,PositionType}, coordinate_matrix::AbstractMatrix{T}) where {PositionDim,PositionType,T<:Number}
    for v in vertices(mesh)
        coordinates!(v, Point(coordinate_matrix[:, id(v)]...))
    end

    return mesh
end
coordinates!(surf::SimplicialSurface, coordinate_matrix::AbstractMatrix{T}) where T<:Number = coordinates!(surf.mesh, coordinate_matrix)

"""
    coordinates(mesh::PolyhedralMesh)

Return the vertex positions of the vertices of `mesh`.
"""
coordinates(mesh::PolyhedralMesh) = [coordinates(v) for v in vertices(mesh)]
coordintes(surf::SimplicialSurface) = coordinates(mesh(surf))

"""
    coordinate_matrix(mesh::PolyhedralMesh)

Return the coordinate matrix of `mesh`. The coordinates of the vertices are the columns of the returned matrix.
"""
function coordinate_matrix(mesh::PolyhedralMesh)
    return hcat([Vector(coordinates(v)) for v in sort(vertices(mesh), by=id)]...)
end
coordinate_matrix(surf::SimplicialSurface) = coordinate_matrix(surf.mesh)

"""
    graph(mesh::PolyhedralMesh)

Return the primal graph of `mesh`.
"""
function graph(mesh::PolyhedralMesh)
    g = Graphs.Graph()

    for v in vertices(mesh)
        add_vertex!(g)
    end

    for (edge_id, e) in edges(mesh)
        add_edge!(g, edge_id[1], edge_id[2])
    end

    return g
end
graph(surf::SimplicialSurface) = Graph(surf.mesh)

RigidityTheoryTools.Framework(mesh::PolyhedralMesh) = RigidityTheoryTools.Framework(graph(mesh), coordinate_matrix(mesh), labels(mesh))