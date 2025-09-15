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

####################################################################################################
# edge length and vertex distance
####################################################################################################
"""
    distance(v1::Vertex, v2::Vertex)

Compute the Euclidean distance between two vertices `v1` and `v2`.

# Examples
```jldoctest
```
"""
distance(v1::Vertex, v2::Vertex) = norm(coordinates(v1) - coordinates(v2))

"""
    distancesq(v1::Vertex, v2::Vertex)

Compute the squared Euclidean distance between two vertices `v1` and `v2`.

# Examples
```jldoctest
```
"""
distancesq(v1::Vertex, v2::Vertex) = sum((coordinates(v1) .- coordinates(v2)) .^ 2)

"""
    elength(e::PrimalEdge)

Compute the Euclidean length of the primal edge `e`.

#Examples
```jldoctest
```
"""
elength(e::PrimalEdge) = distance(tail(e), head(e))


"""
    elengthsq(e::PrimalEdge)

Compute the squared Euclidean length of the primal edge `e`.
#Examples
```jldoctest
```
"""
elengthsq(e::PrimalEdge) = distancesq(tail(e), head(e))

####################################################################################################
# Coordinate transformation of polyhedral meshes
####################################################################################################
"""
    transform!(mesh::PolyhedralMesh{PositionDim}, A::AbstractMatrix, b::AbstractVector) where PositionDim

Transform the vertex coordinates of `mesh` by an affine transformation `x |-> A*x + b`, where `A` is a `PositionDim x PositionDim` matrix and `b` is a vector of length `PositionDim` and return the transformed mesh.
"""
function transform!(mesh::PolyhedralMesh{PositionDim}, A::AbstractMatrix, b::AbstractVector) where PositionDim
    trafo_matrix = vcat(hcat(A, b), hcat(zeros(1, PositionDim), 1))
    ext_coordinate_matrix = vcat(coordinate_matrix(mesh), ones(1, length(vertices(mesh))))
    new_coordinate_matrix = (trafo_matrix*ext_coordinate_matrix)[1:PositionDim, :]
    coordinates!(mesh, new_coordinate_matrix)
    return mesh
end

transform!(surf::SimplicialSurface, A::AbstractMatrix, b::AbstractVector) = transform!(surf.mesh, A, b)

"""
    transform_iso!(mesh::PolyhedralMesh{PositionDim}, A::AbstractMatrix, b::AbstractVector) where PositionDim

Transform the vertex coordinates of `mesh` by an isometry `x |-> A*x + b`, where `A` is an orthogonal `PositionDim x PositionDim` matrix and `b` is a vector of length `PositionDim` and return the transformed mesh.
"""
function transform_iso!(mesh::PolyhedralMesh{PositionDim}, A::AbstractMatrix, b::AbstractVector) where PositionDim
    @assert isapprox(A' * A, I, atol=1e-10) "Matrix A must be orthogonal."
    return transform!(mesh, A, b)
end

transform_iso!(surf::SimplicialSurface, A::AbstractMatrix, b::AbstractVector) = transform_iso!(surf.mesh, A, b)

"""
    kabsch(P::AbstractMatrix, Q::AbstractMatrix; correct_reflection::Bool=true)

Compute the optimal rigid transformation (isometry + translation) that aligns the point set `P` to the point set `Q` in a least-squares sense using the Kabsch algorithm and return the rmsd. The points are given as columns of the matrices `P` and `Q`. If `correct_reflection` is `true`, the function will only return rotations.

# Examples
```jldoctest
julia> P = [0 0 0; 1 0 0; 0 1 0]';

julia> Q = [0 0 0; 0 0 1; 0 1 0]';

julia> R, t, rmsd = kabsch(P,Q);

julia> isapprox(R*P .+ t, Q, atol = 1e-12)
true

julia> isapprox(rmsd, 0.0, atol = 1e-12)
true
```
"""
function kabsch(P::AbstractMatrix, Q::AbstractMatrix; correct_reflection::Bool=true)
    @assert size(P) == size(Q) "Point sets must have the same number of points and dimension."

    N = size(P, 2) # Number of points
    d = size(P, 1) # dimension of underlying space

    # 1. Center the point sets
    centroid_P = sum(P, dims=2) / N
    centroid_Q = sum(Q, dims=2) / N
    P_centered = P .- centroid_P
    Q_centered = Q .- centroid_Q

    # 2. Compute the covariance matrix
    H = P_centered * Q_centered'

    # 3. Compute SVD
    F = svd(H)
    U, V = F.U, F.V # Note: svd returns V, not Vt

    # 4. Calculate the rotation matrix
    R = V * U'

    # 5. Correct for reflection
    if correct_reflection && det(R) < 0
        V_copy = copy(V) # copy before modifying
        V_copy[:, d] *= -1
        R = V_copy * U'
    end

    # 6. Calculate the translation
    t = centroid_Q - R * centroid_P

    # 7. Calculate the RMSD
    P_aligned = R * P .+ t
    rmsd = sqrt(sum((P_aligned - Q) .^ 2) / N)

    return R, t, rmsd
end
randn
"""
    kabsch(points::AbstractVector{<:Point{PositionDim}}, target::AbstractVector{<:Point{PositionDim}}; correct_reflection::Bool=true) where {PositionDim}

Compute the optimal rigid transformation (rotation + translation) that aligns `points` to `target` in a least-squares sense using the Kabsch algorithm. If `correct_reflection` is `true`, the function will only return rotations.
"""
function kabsch(points::AbstractVector{<:Point{PositionDim}}, target::AbstractVector{<:Point{PositionDim}}; correct_reflection::Bool=true) where {PositionDim}
    # see https://www.wikiwand.com/en/articles/Kabsch_algorithm. We do everything transposed compared to the wiki article.
    @assert length(points) == length(target) "Point sets must have the same number of points."
    n = length(points)

    P = hcat(Vector.(points)...)  # n × PositionDim
    Q = hcat(Vector.(target)...)  # n × PositionDim

    return kabsch(P, Q, correct_reflection=correct_reflection)
end