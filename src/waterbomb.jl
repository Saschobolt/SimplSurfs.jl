using LinearAlgebra
using Statistics

# ===================================================================
# Find the rotation around an axis that best aligns two point clouds 
# in the sense that the dihedral angle between their approximating 
# planes is minimized. 
# ===================================================================

# ===================================================================
# PART 1: Rotation Matrix around a Vector (through the origin)
# ===================================================================

"""
    rotate_around_vector(axis::AbstractVector, θ::Real) -> Matrix{Float64}

Computes the 3x3 rotation matrix for a rotation around a given `axis` vector
by an angle `θ` (in radians) using Rodrigues' rotation formula.
"""
function rotate_around_vector(axis::AbstractVector, θ::Real)
    # Ensure the axis is a unit vector
    a = normalize(axis)

    # Cross-product matrix (skew-symmetric)
    K = [0 -a[3] a[2];
        a[3] 0 -a[1];
        -a[2] a[1] 0]

    I_mat = Matrix{Float64}(I, 3, 3)

    # Rodrigues' formula: R = I + sin(θ)K + (1-cos(θ))K²
    R = I_mat + sin(θ) * K + (1 - cos(θ)) * (K * K)
    return R
end

# ===================================================================
# PART 2: Rotation Matrix around an Affine Line
# ===================================================================

"""
    rotate_around_line(p1::AbstractVector, p2::AbstractVector, θ::Real) -> Tuple{Matrix{Float64}, Matrix{Float64}, Vector{Float64}}

Computes the transformation for a rotation around an affine line defined by
points `p1` and `p2` by an angle `θ` (in radians).

Returns a tuple containing:
1.  `M::Matrix{Float64}`: The 4x4 homogeneous transformation matrix.
2.  `R::Matrix{Float64}`: The 3x3 rotation part of the transformation.
3.  `t::Vector{Float64}`: The 3x1 translation part of the transformation.
"""
function rotate_around_line(p1::AbstractVector, p2::AbstractVector, θ::Real)
    # 1. Define the rotation axis direction
    axis = p2 - p1

    # 2. Get the 3x3 rotation part for an axis through the origin
    R = rotate_around_vector(axis, θ)

    # 3. Create 4x4 homogeneous matrices for the translation-rotation-translation sequence
    T_to_origin = [1 0 0 -p1[1]; 0 1 0 -p1[2]; 0 0 1 -p1[3]; 0 0 0 1]
    T_back = [1 0 0 p1[1]; 0 1 0 p1[2]; 0 0 1 p1[3]; 0 0 0 1]

    R_4x4 = Matrix{Float64}(I, 4, 4)
    R_4x4[1:3, 1:3] = R

    # 4. Combine transformations to get the final 4x4 matrix
    # The sequence is: translate to origin, rotate, then translate back.
    M = T_back * R_4x4 * T_to_origin

    # 5. Extract the translation vector from the final 4x4 matrix
    t = M[1:3, 4]

    # 6. Return all three components
    return (M, R, t)
end

# ===================================================================
# PART 3: Find the Approximating Plane of a Point Cloud
# ===================================================================

"""
    find_plane(points::AbstractMatrix) -> Tuple{Vector{Float64}, Vector{Float64}}

Finds the best-fit plane for a 3xN point cloud using PCA.

Returns a tuple containing:
1.  `normal::Vector{Float64}`: The unit normal vector to the plane.
2.  `centroid::Vector{Float64}`: The centroid of the point cloud.
"""
function find_plane(points::AbstractMatrix)
    @assert size(points, 1) == 3 "Input matrix must be 3xN."

    # 1. Calculate the centroid
    centroid = vec(mean(points, dims=2))

    # 2. Center the data
    points_centered = points .- centroid

    # 3. Compute the covariance matrix
    cov_matrix = points_centered * points_centered'

    # 4. Find the eigenvectors using SVD
    F = svd(cov_matrix)

    # 5. The normal is the eigenvector with the smallest eigenvalue (last column of U)
    normal = F.U[:, 3]

    return (normal, centroid)
end

# ===================================================================
# PART 4: Find the Dihedral Angle between Two Planes
# ===================================================================

"""
    dihedral_angle(plane1, plane2) -> Float64

Calculates the dihedral angle (in radians) between two planes.
Each plane is a tuple `(normal, centroid)`.
"""
function dihedral_angle(plane1::Tuple, plane2::Tuple)
    n1 = normalize(plane1[1])
    n2 = normalize(plane2[1])

    # Use clamp to avoid domain errors from floating point inaccuracies
    dot_product = clamp(dot(n1, n2), -1.0, 1.0)

    return acos(dot_product)
end


# ===================================================================
# PART 5: Putting It All Together
# ===================================================================

"""
    align_planes(P::AbstractMatrix, Q::AbstractMatrix, p1::AbstractVector, p2::AbstractVector)

Calculates the optimal transformation to align the approximating plane of point
cloud `P` with the approximating plane of point cloud `Q` via a rotation
around the affine line defined by `p1` and `p2`.

Returns a tuple containing:
1. `M::Matrix{Float64}`: The 4x4 homogeneous transformation matrix.
2. `R::Matrix{Float64}`: The 3x3 rotation part of the transformation.
3. `t::Vector{Float64}`: The 3x1 translation part of the transformation.
"""
function align_planes(P::AbstractMatrix, Q::AbstractMatrix, p1::AbstractVector, p2::AbstractVector)
    # 1. Find the normal vector for each point cloud's plane
    (n_p, _) = find_plane(P)
    (n_q, _) = find_plane(Q)

    # 2. Define the rotation axis direction vector
    a = normalize(p2 - p1)

    # 3. Project normals onto the plane perpendicular to the axis
    u = n_p - dot(n_p, a) * a # Component of n_p perpendicular to a
    v = n_q - dot(n_q, a) * a # Component of n_q perpendicular to a

    # 4. Calculate the optimal rotation angle θ
    # This is the signed angle between the projected vectors u and v
    θ = atan(dot(cross(u, v), a), dot(u, v))

    # 5. Construct the full 4x4 affine transformation matrix
    M, _, _ = rotate_around_line(p1, p2, θ)

    # 6. Extract the 3x3 rotation and 3x1 translation components
    R = M[1:3, 1:3]
    t = M[1:3, 4]

    return (M, R, t)
end

# ===================================================================
# Helper functions
# ===================================================================

"""
    _permutation_matrix(p::AbstractVector{<:Integer})

Return the permutation matrix corresponding to the permutation `p`.
"""
_permutation_matrix(p::AbstractVector{<:Integer}) = I(length(p))[p, :]

# Reflection of x through point p
"""
    _reflect_through_point(x::AbstractVector{<:Number}, p::AbstractVector{<:Number})

Reflection of `x` through point `p`.
"""
function _reflect_through_point(x::AbstractVector{<:Number}, p::AbstractVector{<:Number})
    return 2 .* p .- x
end

# determine whether the points represented by the columns of M are collinear
"""
    _collinear(M::AbstractMatrix{<:Number})

Return whether the points represented by the columns of `M` are collinear.
"""
function _collinear(M::AbstractMatrix{<:Number})
    p1, p2, p3 = M[:, 1], M[:, 2], M[:, 3]
    v1 = p2 - p1
    v2 = p3 - p1
    return rank(hcat(v1, v2)) ≤ 1
end

# ===================================================================
# WaterbombCell Struct and Functions
# ===================================================================

mutable struct WaterbombCell
    mesh::PolyhedralMesh{3,Float64}
    s_connects::Tuple{Int,Int} # which vertices does s connect. Either (3,5) or (2,4) or (5,7) or (2,6)
    t_connects::Tuple{Int,Int} # which vertices does t connect. Either (3,7) or (4,6)
    dependent_on::Vector{WaterbombCell} # cells that this cell depends on for its coordinates
    dependent_cells::Vector{WaterbombCell} # cells that depend on this cell for their coordinates
end

a(cell::WaterbombCell) = elength(edges(cell.mesh)[(3, 4)])  # edge length of edges (3,4), (6,7)
b(cell::WaterbombCell) = elength(edges(cell.mesh)[(1, 3)]) # edge length of edges (1,3), (1,4), (1,6), (1,7)
d(cell::WaterbombCell) = elength(edges(cell.mesh)[(2, 3)]) # edge length of edges (2,3), (4,5), (5,6), (2,7)
r(cell::WaterbombCell) = distance(vertices(cell.mesh)[2], vertices(cell.mesh)[5]) # distance between vertices 2 and 5
s(cell::WaterbombCell) = distance(vertices(cell.mesh)[cell.s_connects[1]], vertices(cell.mesh)[cell.s_connects[2]]) # distance between vertices s_connects
t(cell::WaterbombCell) = distance(vertices(cell.mesh)[cell.t_connects[1]], vertices(cell.mesh)[cell.t_connects[2]]) # distance between vertices t_connects
coordinate_matrix(cell::WaterbombCell) = coordinate_matrix(cell.mesh)
coordinates!(cell::WaterbombCell, coordinate_matrix::AbstractMatrix{T}) where T<:Number = coordinates!(cell.mesh, coordinate_matrix)
mesh(cell) = cell.mesh

"""
    coordinates(a::Float64=1, b::Float64=sqrt(2), d::Float64=1, r::Float64=1, s::Float64=sqrt(2), t::Float64=2, s_connects::Tuple{Int,Int}=(3, 5), t_connects::Tuple{Int,Int}=(3, 7))

Compute the coordinates of a waterbomb cell in standard position where `s` connects `s_connects` and `t` connects `t_connects`.
"""
function coordinates(a::Real=1, b::Real=sqrt(2), d::Real=1, r::Real=1, s::Real=sqrt(2), t::Real=2, s_connects::Tuple{<:Integer,<:Integer}=(3, 5), t_connects::Tuple{<:Integer,<:Integer}=(3, 7))
    # coordinates of a cell in standard position where s connects (3,5) and t connects (3,7)
    p1 = [0, 0, 0]
    p2 = [r / 2, 0, sqrt(-r^2 + 4) / 2]
    p3 = [-(d^2 - s^2) / r / 2, 0.1e1 / r * sqrt(-(b^4 * r^2 - b^2 * d^2 * r^2 + b^2 * r^4 - b^2 * r^2 * s^2 + d^2 * r^2 * s^2 - 2 * r^2 * b^2 + d^4 - r^2 * d^2 - 2 * d^2 * s^2 - r^2 * s^2 + s^4 + r^2) / (-r^2 + 4)), sqrt(-r^2 + 4) * (-2 * b^2 + d^2 + s^2 - 2) / (r^2 - 4) / 2]
    p4 = [(2 * sqrt(-(b^4 - 2 * b^2 * s^2 + s^4 - 2 * b^2 - 2 * s^2 + 1) * r^2) * sqrt(-r^2 + 4) * sqrt((b^2 * r^4 + (b^4 + (-d^2 - s^2 - 2) * b^2 + (s^2 - 1) * d^2 - s^2 + 1) * r^2 + (d - s)^2 * (d + s)^2) / (r^2 - 4)) * sqrt((a^2 * b^4 + ((-d^2 - s^2 - 2) * a^2 + (d - s)^2 * (d + s)^2) * b^2 + (a^2 + (d^2 - 1) * s^2 - d^2 + 1) * a^2) / (b^4 + (-2 * s^2 - 2) * b^2 + s^4 - 2 * s^2 + 1)) + r * (((-a^2 + 2 * b^2) * s^2 + (a^2 - 2 * d^2) * b^2 + a^2) * r^2 - 2 * ((-b^2 / 2 + d^2 / 2 - 1 / 2) * s^2 + b^4 / 2 + (-d^2 / 2 - 1) * b^2 + a^2 - d^2 / 2 + 1 / 2) * (d + s) * (d - s))) / (2 * s^4 + (-4 * b^2 - 4) * s^2 + 2 * b^4 - 4 * b^2 + 2) / r^2, (((b^2 - s^2 + 1) * r^2 - 2 * d^2 + 2 * s^2) * sqrt(-r^2 + 4) * sqrt(-(b^4 - 2 * b^2 * s^2 + s^4 - 2 * b^2 - 2 * s^2 + 1) * r^2) * sqrt((a^2 * b^4 + ((-d^2 - s^2 - 2) * a^2 + (d - s)^2 * (d + s)^2) * b^2 + (a^2 + (d^2 - 1) * s^2 - d^2 + 1) * a^2) / (b^4 + (-2 * s^2 - 2) * b^2 + s^4 - 2 * s^2 + 1)) + 2 * (1 / 2 + b^4 / 2 + (-d^2 / 2 - s^2 / 2 - 1) * b^2 + (d^2 - 1) * s^2 / 2 + a^2 - d^2 / 2) * (r + 2) * r * (r - 2) * sqrt((b^2 * r^4 + (b^4 + (-d^2 - s^2 - 2) * b^2 + (s^2 - 1) * d^2 - s^2 + 1) * r^2 + (d - s)^2 * (d + s)^2) / (r^2 - 4))) / (b^4 - 2 * b^2 * s^2 + s^4 - 2 * b^2 - 2 * s^2 + 1) / r^2 / (r^2 - 4), ((2 * r^2 - 8) * sqrt(-(b^4 - 2 * b^2 * s^2 + s^4 - 2 * b^2 - 2 * s^2 + 1) * r^2) * sqrt((b^2 * r^4 + (b^4 + (-d^2 - s^2 - 2) * b^2 + (s^2 - 1) * d^2 - s^2 + 1) * r^2 + (d - s)^2 * (d + s)^2) / (r^2 - 4)) * sqrt((a^2 * b^4 + ((-d^2 - s^2 - 2) * a^2 + (d - s)^2 * (d + s)^2) * b^2 + (a^2 + (d^2 - 1) * s^2 - d^2 + 1) * a^2) / (b^4 + (-2 * s^2 - 2) * b^2 + s^4 - 2 * s^2 + 1)) - sqrt(-r^2 + 4) * r * (((a^2 - 2 * d^2 + 2 * s^2) * b^2 - a^2 * s^2 + a^2) * r^2 + 2 * b^6 + (-3 * d^2 - 3 * s^2 - 2) * b^4 + (s^4 + (4 * d^2 - 10) * s^2 + d^4 + 6 * d^2 - 2) * b^2 + (-d^2 + 1) * s^4 + (-d^4 + 2 * a^2 + 4 * d^2 - 3) * s^2 + 2 + d^4 + (-2 * a^2 - 3) * d^2)) / (b^4 - 2 * b^2 * s^2 + s^4 - 2 * b^2 - 2 * s^2 + 1) / r / (r^2 - 4) / 2]
    p5 = [-r / 2, 0, sqrt(-r^2 + 4) / 2]
    dist = t # distance between vertices 3 and 7
    t = sqrt((-2 * sqrt(dist^2 * (d^4 + (-2 * b^2 - 2) * d^2 + b^4 - 2 * b^2 + dist^2 + 1) * (b^2 * r^4 + ((-b^2 + s^2 - 1) * d^2 + b^4 + (-s^2 - 2) * b^2 - s^2 + 1) * r^2 + (d - s)^2 * (d + s)^2)) + (2 * dist^2 + d^4 + (-2 * b^2 - 2) * d^2 + b^4 - 2 * b^2 + 1) * s^2 + dist^2 * (-r^2 * b^2 + (r^2 - 2) * d^2 - r^2)) / (b^4 + (-2 * d^2 - 2) * b^2 + d^4 - 2 * d^2 + 1))
    p6 = [(2 * sqrt(-(b^4 - 2 * b^2 * t^2 + t^4 - 2 * b^2 - 2 * t^2 + 1) * r^2) * sqrt(-r^2 + 4) * sqrt((b^2 * r^4 + (b^4 + (-d^2 - t^2 - 2) * b^2 + (t^2 - 1) * d^2 - t^2 + 1) * r^2 + (d - t)^2 * (d + t)^2) / (r^2 - 4)) * sqrt((a^2 * b^4 + ((-d^2 - t^2 - 2) * a^2 + (d - t)^2 * (d + t)^2) * b^2 + (a^2 + (d^2 - 1) * t^2 - d^2 + 1) * a^2) / (b^4 + (-2 * t^2 - 2) * b^2 + t^4 - 2 * t^2 + 1)) + r * (((-a^2 + 2 * b^2) * t^2 + (a^2 - 2 * d^2) * b^2 + a^2) * r^2 - 2 * ((-b^2 / 2 + d^2 / 2 - 1 // 2) * t^2 + b^4 / 2 + (-d^2 / 2 - 1) * b^2 + a^2 - d^2 / 2 + 1 // 2) * (d + t) * (d - t))) / (2 * t^4 + (-4 * b^2 - 4) * t^2 + 2 * b^4 - 4 * b^2 + 2) / r^2, (-((b^2 - t^2 + 1) * r^2 - 2 * d^2 + 2 * t^2) * sqrt(-r^2 + 4) * sqrt(-(b^4 - 2 * b^2 * t^2 + t^4 - 2 * b^2 - 2 * t^2 + 1) * r^2) * sqrt((a^2 * b^4 + ((-d^2 - t^2 - 2) * a^2 + (d - t)^2 * (d + t)^2) * b^2 + (a^2 + (d^2 - 1) * t^2 - d^2 + 1) * a^2) / (b^4 + (-2 * t^2 - 2) * b^2 + t^4 - 2 * t^2 + 1)) - 2 * (1 // 2 + b^4 / 2 + (-d^2 / 2 - t^2 / 2 - 1) * b^2 + (d^2 - 1) * t^2 / 2 + a^2 - d^2 / 2) * (r + 2) * r * (r - 2) * sqrt((b^2 * r^4 + (b^4 + (-d^2 - t^2 - 2) * b^2 + (t^2 - 1) * d^2 - t^2 + 1) * r^2 + (d - t)^2 * (d + t)^2) / (r^2 - 4))) / (b^4 - 2 * b^2 * t^2 + t^4 - 2 * b^2 - 2 * t^2 + 1) / r^2 / (r^2 - 4), ((2 * r^2 - 8) * sqrt(-(b^4 - 2 * b^2 * t^2 + t^4 - 2 * b^2 - 2 * t^2 + 1) * r^2) * sqrt((b^2 * r^4 + (b^4 + (-d^2 - t^2 - 2) * b^2 + (t^2 - 1) * d^2 - t^2 + 1) * r^2 + (d - t)^2 * (d + t)^2) / (r^2 - 4)) * sqrt((a^2 * b^4 + ((-d^2 - t^2 - 2) * a^2 + (d - t)^2 * (d + t)^2) * b^2 + (a^2 + (d^2 - 1) * t^2 - d^2 + 1) * a^2) / (b^4 + (-2 * t^2 - 2) * b^2 + t^4 - 2 * t^2 + 1)) - sqrt(-r^2 + 4) * r * (((a^2 - 2 * d^2 + 2 * t^2) * b^2 - a^2 * t^2 + a^2) * r^2 + 2 * b^6 + (-3 * d^2 - 3 * t^2 - 2) * b^4 + (t^4 + (4 * d^2 - 10) * t^2 + d^4 + 6 * d^2 - 2) * b^2 + (-d^2 + 1) * t^4 + (-d^4 + 2 * a^2 + 4 * d^2 - 3) * t^2 + 2 + d^4 + (-2 * a^2 - 3) * d^2)) / (b^4 - 2 * b^2 * t^2 + t^4 - 2 * b^2 - 2 * t^2 + 1) / r / (r^2 - 4) / 2]
    p7 = [(-d^2 + t^2) / r / 2, -sqrt((b^2 * r^4 + (b^4 + (-d^2 - t^2 - 2) * b^2 + (t^2 - 1) * d^2 - t^2 + 1) * r^2 + (d - t)^2 * (d + t)^2) / (r^2 - 4)) / r, (2 * b^2 - d^2 - t^2 + 2) * (-r^2 + 4)^(-1 // 2) / 2]
    coords = hcat(p1, p2, p3, p4, p5, p6, p7)

    # transform coords if necessary to get other cases.
    if s_connects == (3, 5) && t_connects == (3, 7)
        return coords
    elseif s_connects == (5, 7) && t_connects == (3, 7)
        # swap vertex pairs (3,7), (4,6)
        p = _permutation_matrix([1, 2, 7, 6, 5, 4, 3])
        return coords * p
    elseif s_connects == (2, 4) && t_connects == (4, 6)
        # swap vertex pairs (3,4), (2,5), (6,7)
        p = _permutation_matrix([1, 5, 4, 3, 2, 7, 6])
        return coords * p
    elseif s_connects == (2, 6) && t_connects == (4, 6)
        # swap vertex pairs (3,6), (2,5), (4,7)
        p = _permutation_matrix([1, 5, 6, 7, 2, 3, 4])
        return coords * p
    else
        error("Invalid combination of s_connects and t_connects.")
    end
end

function WaterbombCell(a::Real=1, b::Real=sqrt(2), d::Real=1, r::Real=1, s::Real=sqrt(2), t::Real=2)
    faces = [
        [1, 2, 3],
        [1, 3, 4],
        [1, 4, 5],
        [1, 5, 6],
        [1, 6, 7],
        [1, 7, 2]
    ]
    m = PolyhedralMesh{3,Float64}(faces; labels=1:7)
    coords = coordinates(a, b, d, r, s, t, (3, 5), (3, 7))
    coordinates!(m, coords)
    return WaterbombCell(m, (3, 5), (3, 7), WaterbombCell[], WaterbombCell[])
end

"""
    shared_edge(cell1::WaterbombCell, cell2::WaterbombCell; atol::Float64=1e-10)

Return the shared edges between `cell1` and `cell2` as a tuple `(e1, e2)`, where `e1` is the edge in `cell1` and `e2` is the edge in `cell2`. If no shared edge exists, return `nothing`.
The function checks for shared edges by comparing the coordinates of the tail and head vertices of each edge in both cells, within a specified tolerance `atol`. The cells must have the same edge lengths `a`, `b`, and `d`.
"""
function shared_edge(cell1::WaterbombCell, cell2::WaterbombCell; atol::Real=1e-10)
    @assert a(cell1) ≈ a(cell2) "Cells must have the same edge length a"
    @assert b(cell1) ≈ b(cell2) "Cells must have the same edge length b"
    @assert d(cell1) ≈ d(cell2) "Cells must have the same edge length d"

    m1 = cell1.mesh
    m2 = cell2.mesh

    for e1 in values(edges(m1))
        for e2 in values(edges(m2))
            if (isapprox(coordinates(tail(e1)), coordinates(tail(e2)), atol=atol) && isapprox(coordinates(head(e1)), coordinates(head(e2)), atol=atol)) || (isapprox(coordinates(head(e1)), coordinates(tail(e2)), atol=atol) && isapprox(coordinates(tail(e1)), coordinates(head(e2)), atol=atol))
                return (e1, e2)
            end
        end
    end
end

"""
    coordinates(cell1::WaterbombCell, cell2::WaterbombCell, r::Float64, s::Float64, t::Float64, rot::Float64; atol::Float64=1e-10)

Compute the vertex coordinates of a new `WaterbombCell` that shares an edge with `cell1` and `cell2`, given the parameters `r`, `s`, and `t`.
If the vertices, the new cell attaches to are collinear, it will be rotated around the axis defined by these vertices by an angle `rot` (in radians).
The new cell will be positioned such that it shares the appropriate edges with `cell1` and `cell2` based on their shared edge. 
If the shared edge of `cell1` and `cell2` is
- edge (2,3) of `cell1`, then the new cell will have 
    - vertex 2 at the vertex position of vertex 4 of `cell2`
    - vertex 7 at the vertex position of vertex 3 of `cell1`
    - vertex 6 at the vertex position of vertex 4 of `cell1`
- edge (5,6) of `cell1`, then the new cell will have 
    - vertex 5 at the vertex position of vertex 7 of `cell2`
    - vertex 4 at the vertex position of vertex 6 of `cell1`
    - vertex 3 at the vertex position of vertex 7 of `cell1`
- edge (4,5) of `cell1`, then the new cell will have 
    - vertex 5 at the vertex position of vertex 3 of `cell2`
    - vertex 6 at the vertex position of vertex 4 of `cell1`
    - vertex 7 at the vertex position of vertex 3 of `cell1`
- edge (2,7) of `cell1`, then the new cell will have 
    - vertex 2 at the vertex position of vertex 6 of `cell2`
    - vertex 3 at the vertex position of vertex 7 of `cell1`
    - vertex 4 at the vertex position of vertex 6 of `cell1`
- edge (3,4) of `cell1`, then the new cell will have 
    - vertex 3 at the vertex position of vertex 5 of `cell2`
    - vertex 2 at the vertex position of vertex 4 of `cell1`
    - vertex 7 at the vertex position of vertex 5 of `cell1`
- edge (6,7) of `cell1`, then the new cell will have 
    - vertex 2 at the vertex position of vertex 4 of `cell2`
    - vertex 7 at the vertex position of vertex 5 of `cell1`
    - vertex 2 at the vertex position of vertex 6 of `cell1`
The function returns a 3x7 matrix where each column represents the coordinates of a vertex of the new cell.
"""
function coordinates(cell1::WaterbombCell, cell2::WaterbombCell, r::Real, s::Real, t::Real, rot::Real=0; atol::Real=1e-10)
    @assert a(cell1) ≈ a(cell2) "Cells must have the same edge length a"
    @assert b(cell1) ≈ b(cell2) "Cells must have the same edge length b"
    @assert d(cell1) ≈ d(cell2) "Cells must have the same edge length d"

    m1 = cell1.mesh
    m2 = cell2.mesh

    a1 = a(cell1)
    b1 = b(cell1)
    d1 = d(cell1)

    # determine, which parameter (r,s,t) of the new cell is determined by the other two cells
    shared_edge_1, shared_edge_2 = shared_edge(cell1, cell2, atol=atol)
    if isnothing(shared_edge_1) || isnothing(shared_edge_2)
        throw(ArgumentError("Cells do not share an edge"))
    end
    v1 = tail(shared_edge_1)
    v2 = head(shared_edge_1)
    edge_id = (id(v1), id(v2)) # edge id in the form (tail_id, head_id). Is sorted, because edges in values(edges(mesh)) have tail_id < head_id. 
    # case 1: shared edge is an edge of length d, i.e. either (2,3), (4,5), (5,6), or (2,7).
    # In this case the new cell's s is determined and it will be placed so that it shares the adjacent edge of length a with cell1 and the adjacent edge of length d with cell2.
    if edge_id == (2, 3) # then corresponding edge of cell2 is (5,6)
        other_vertex_cell1 = vertices(m1)[4]
        other_vertex_cell2 = vertices(m2)[4]
        shared_vertex = vertices(m1)[3]
        s = distance(other_vertex_cell1, other_vertex_cell2)
        s_connects = (2, 6)
        t_connects = (4, 6)
        coords = coordinates(a1, b1, d1, r, s, t, s_connects, t_connects)
        R, trans = kabsch(coords[:, [2, 7, 6]], hcat(coordinates(other_vertex_cell2), coordinates(shared_vertex), coordinates(other_vertex_cell1)))
        coords = R * coords .+ trans
    elseif edge_id == (5, 6) # then corresponding edge of cell2 is (2,3)
        other_vertex_cell1 = vertices(m1)[7]
        other_vertex_cell2 = vertices(m2)[7]
        shared_vertex = vertices(m1)[6]
        s = distance(other_vertex_cell1, other_vertex_cell2)
        s_connects = (3, 5)
        t_connects = (3, 7)
        coords = coordinates(a1, b1, d1, r, s, t, s_connects, t_connects)
        R, trans = kabsch(coords[:, [5, 4, 3]], hcat(coordinates(other_vertex_cell2), coordinates(shared_vertex), coordinates(other_vertex_cell1)))
        coords = R * coords .+ trans
    elseif edge_id == (4, 5) # then corresponding edge of cell2 is (2,7)
        other_vertex_cell1 = vertices(m1)[3]
        other_vertex_cell2 = vertices(m2)[3]
        shared_vertex = vertices(m1)[4]
        s = distance(other_vertex_cell1, other_vertex_cell2)
        s_connects = (5, 7)
        t_connects = (3, 7)
        coords = coordinates(a1, b1, d1, r, s, t, s_connects, t_connects)
        R, trans = kabsch(coords[:, [5, 6, 7]], hcat(coordinates(other_vertex_cell2), coordinates(shared_vertex), coordinates(other_vertex_cell1)))
        coords = R * coords .+ trans
    elseif edge_id == (2, 7) # then corresponding edge of cell2 is (4,5)
        other_vertex_cell1 = vertices(m1)[6]
        other_vertex_cell2 = vertices(m2)[6]
        shared_vertex = vertices(m1)[7]
        s = distance(other_vertex_cell1, other_vertex_cell2)
        s_connects = (2, 4)
        t_connects = (4, 6)
        coords = coordinates(a1, b1, d1, r, s, t, s_connects, t_connects)
        R, trans = kabsch(coords[:, [2, 3, 4]], hcat(coordinates(other_vertex_cell2), coordinates(shared_vertex), coordinates(other_vertex_cell1)))
        coords = R * coords .+ trans
        # case 2: shared edge is an edge of length a, i.e. either (3,4) or (6,7).
        # In this case the new cell's t is determined and it will be placed so that it shares either vertex 4 or vertex 7 with cell1.
    elseif edge_id == (3, 4) # then the corresponding edge of cell2 is (6,7). The new cell shares vertices 4 and 5 with cell1 and vertex 6 (shared vertex between cell1 and cell2) and 5 with cell2.
        other_vertex_cell1 = vertices(m1)[5]
        other_vertex_cell2 = vertices(m2)[5]
        shared_vertex = vertices(m1)[4]
        t = distance(other_vertex_cell1, other_vertex_cell2)
        s_connects = (3, 5)
        t_connects = (3, 7)
        coords = coordinates(a1, b1, d1, r, s, t, s_connects, t_connects)
        R, trans = kabsch(coords[:, [3, 2, 7]], hcat(coordinates(other_vertex_cell2), coordinates(shared_vertex), coordinates(other_vertex_cell1)))
        coords = R * coords .+ trans
        if _collinear(coords[:, [3, 2, 7]]) # if points are collinear, rotate around axis defined by these points by angle rot
            # first align them so that dihedral angle between approximating planes is minimized.
            M, R, trans = align_planes(coords, hcat(coordinate_matrix(m1), coordinate_matrix(m2)), coords[:, 3], coords[:, 2])
            coords = R * coords .+ trans
            # then rotate by the defined angle.
            M, R, t = rotate_around_line(coords[:, 3], coords[:, 2], rot)
            coords = R * coords .+ t
        end
    elseif edge_id == (6, 7) # then the corresponding edge of cell2 is (3,4). The new cell shares vertices 2 and 7 with cell1 and vertex 3 (shared vertex between cell1 and cell2) and 2 with cell2.
        other_vertex_cell1 = vertices(m1)[2]
        other_vertex_cell2 = vertices(m2)[2]
        shared_vertex = vertices(m1)[7]
        t = distance(other_vertex_cell1, other_vertex_cell2)
        s_connects = (2, 4)
        t_connects = (4, 6)
        coords = coordinates(a1, b1, d1, r, s, t, s_connects, t_connects)
        R, trans = kabsch(coords[:, [4, 5, 6]], hcat(coordinates(other_vertex_cell2), coordinates(shared_vertex), coordinates(other_vertex_cell1)))
        coords = R * coords .+ trans
        if _collinear(coords[:, [4, 5, 6]]) # if points are collinear, rotate around axis defined by these points by angle rot
            # first align them so that dihedral angle between approximating planes is minimized.
            M, R, trans = align_planes(coords, hcat(coordinate_matrix(m1), coordinate_matrix(m2)), coords[:, 4], coords[:, 5])
            coords = R * coords .+ trans
            # then rotate by the defined angle.
            M, R, t = rotate_around_line(coords[:, 4], coords[:, 5], rot)
            coords = R * coords .+ t
        end
    end

    return coords
end

function attach_cell!(cell1::WaterbombCell, cell2::WaterbombCell; r::Real=1, s::Real=sqrt(2), t::Real=1.5, rot::Real=0, atol::Real=1e-10)
    # compute coordinates of new cell
    faces = [
        [1, 2, 3],
        [1, 3, 4],
        [1, 4, 5],
        [1, 5, 6],
        [1, 6, 7],
        [1, 7, 2]
    ]
    m = PolyhedralMesh{3,Float64}(faces; labels=1:7)
    coords = coordinates(cell1, cell2, r, s, t, rot; atol=atol)
    coordinates!(m, coords)

    # determine s_connects and t_connects of new cell
    shared_edge_1, shared_edge_2 = shared_edge(cell1, cell2, atol=atol)
    if isnothing(shared_edge_1) || isnothing(shared_edge_2)
        throw(ArgumentError("Cells do not share an edge"))
    end
    v1 = tail(shared_edge_1)
    v2 = head(shared_edge_1)
    edge_id = (id(v1), id(v2)) # edge id in the form (tail_id, head_id). Is sorted, because edges in values(edges(mesh)) have tail_id < head_id. 
    if edge_id == (2, 3) # then corresponding edge of cell2 is (5,6)
        s_connects = (2, 6)
        t_connects = (4, 6)
    elseif edge_id == (5, 6) # then corresponding edge of cell2 is (2,3)
        s_connects = (3, 5)
        t_connects = (3, 7)
    elseif edge_id == (4, 5) # then corresponding edge of cell2 is (2,7)
        s_connects = (5, 7)
        t_connects = (3, 7)
    elseif edge_id == (2, 7) # then corresponding edge of cell2 is (4,5)
        s_connects = (2, 4)
        t_connects = (4, 6)
    elseif edge_id == (3, 4) # then the corresponding edge of cell2 is (6,7)
        s_connects = (3, 5)
        t_connects = (3, 7)
    elseif edge_id == (6, 7) # then the corresponding edge of cell2 is (3,4)
        s_connects = (2, 4)
        t_connects = (4, 6)
    end

    new_cell = WaterbombCell(m, s_connects, t_connects, WaterbombCell[cell1, cell2], WaterbombCell[])
    push!(cell1.dependent_cells, new_cell)
    push!(cell2.dependent_cells, new_cell)
    return new_cell
end

function coordinates(cell1::WaterbombCell; r::Real=1, s::Real=sqrt(2), t::Real=2, rot::Real=0, attach_to_a::Bool=true)
    m = mesh(cell1)
    coords = coordinates(a(cell1), b(cell1), d(cell1), r, s, t, (3, 5), (3, 7))

    if attach_to_a # new cell is attached to edge (3,4) of cell1. Corresponding verts of new cell are (7,6)
        v1 = vertices(m)[3]
        v2 = vertices(m)[4]
        R, t = kabsch(coords[:, [7, 6]], hcat(coordinates(v1), coordinates(v2)))
        coords = R * coords .+ t
        M, R, t = align_planes(coords, coordinate_matrix(cell1), coordinates(v1), coordinates(v2))
        coords = R * coords .+ t
        M, R, t = rotate_around_line(coordinates(v1), coordinates(v2), rot)
        coords = R * coords .+ t
    else # new cell is attached to edge (2,3) of cell1. Corresponding verts of new cell are (6,5)
        v1 = vertices(m)[2]
        v2 = vertices(m)[3]
        R, t = kabsch(coords[:, [6, 5]], hcat(coordinates(v1), coordinates(v2)))
        coords = R * coords .+ t
        M, R, t = align_planes(coords, coordinate_matrix(cell1), coordinates(v1), coordinates(v2))
        coords = R * coords .+ t
        M, R, t = rotate_around_line(coordinates(v1), coordinates(v2), rot)
        coords = R * coords .+ t
    end

    return coords
end

function attach_cell!(cell1::WaterbombCell; r::Real=1, s::Real=sqrt(2), t::Real=2, rot::Real=0, attach_to_a::Bool=true)
    m = cell1.mesh

    new_cell = WaterbombCell(a(cell1), b(cell1), d(cell1), r, s, t)
    coords = coordinates(cell1, r=r, s=s, t=t, rot=rot, attach_to_a=attach_to_a)
    coordinates!(new_cell.mesh, coords)

    push!(cell1.dependent_cells, new_cell)
    push!(new_cell.dependent_on, cell1)

    return new_cell
end

function update_coordinates!(cell::WaterbombCell; r::Real=1, s::Real=sqrt(2), t::Real=1.5, rot::Real=0, atol::Real=1e-10)
    m = cell.mesh

    if length(cell.dependent_on) == 0
        coords = coordinates(a(cell), b(cell), d(cell), r, s, t, cell.s_connects, cell.t_connects)
        coordinates!(m, coords)
    elseif length(cell.dependent_on) == 1
        e1, e2 = shared_edge(cell, cell.dependent_on[1], atol=atol)
        if isnothing(e1) || isnothing(e2)
            throw(ArgumentError("Cells do not share an edge"))
        end
        v1 = tail(e1)
        v2 = head(e1)
        edge_id = (id(v1), id(v2)) # edge id in the form (tail_id, head_id). Is sorted, because edges in values(edges(mesh)) have tail_id < head_id.
        if edge_id in [(3, 4), (6, 7)]
            attach_to_a = true
        elseif edge_id in [(2, 3), (4, 5), (5, 6), (2, 7)]
            attach_to_a = false
        else
            throw(ArgumentError("Cells do not share a valid edge"))
        end

        dependent_cell = cell.dependent_on[1]

        coords = coordinates(dependent_cell, r=r, s=s, t=t, rot=rot, attach_to_a=attach_to_a)
        coordinates!(m, coords)
    elseif length(cell.dependent_on) == 2
        coords = coordinates(cell.dependent_on[1], cell.dependent_on[2], r, s, t, atol=atol)
        coordinates!(m, coords)
    else
        error("A WaterbombCell can only depend on up to two other cells.")
    end

    # update dependent cells
    for dependent_cell in cell.dependent_cells
        update_coordinates!(dependent_cell, r=SimplSurfs.r(dependent_cell), s=SimplSurfs.s(dependent_cell), t=SimplSurfs.t(dependent_cell), atol=atol)
    end

    return cell
end