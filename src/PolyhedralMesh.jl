############################################################################
# types
############################################################################
mutable struct Vertex
    id::Int
    label::Union{Int,String}
    coord::Union{Nothing,Vector{Union{AbstractAlgebra.RingElem,Number}}}
    edge::Any # Edge    # a primal edge with this vertex as tail
    mesh::Any # PolyhedralMesh

    Vertex() = new()
end

function Base.show(io::IO, v::Vertex)
    if id(v) != 0
        print(io, "V$(v.id) ($(v.label))")
    else
        print(io, "null vertex")
    end
end

mutable struct Face
    id::Int
    edge::Any # Edge    # a dual edge with this face as tail
    mesh::Any # PolyhedralMesh

    Face() = new()
end

function Base.show(io::IO, f::Face)
    if id(f) > 0
        print(io, "F $(f.id)")
    elseif id(f) < 0
        print(io, "hole $(f.id)")
    else
        print(io, "null face")
    end
end

abstract type Edge end

mutable struct PrimalEdge <: Edge
    tail::Vertex
    next::PrimalEdge
    rot::Any # DualEdge
    flip::PrimalEdge

    mesh::Any # PolyhedralMesh

    PrimalEdge() = new()
end

mutable struct DualEdge <: Edge
    tail::Face
    next::DualEdge
    rot::PrimalEdge # DualEdge
    flip::DualEdge

    mesh::Any # PolyhedralMesh

    DualEdge() = new()
end

function Base.show(io::IO, e::Edge)
    s_head = string(head(e))
    s_tail = string(tail(e))
    s_left = string(left(e))
    s_right = string(right(e))

    left_padding = "    "
    right_padding = "    "
    line3 = s_left * left_padding * "|" * right_padding * s_right
    center_pos = length(s_left) + length(left_padding)
    pad_head = max(0, center_pos - floor(Int, length(s_head) / 2))
    pad_tail = max(0, center_pos - floor(Int, length(s_tail) / 2))

    println(io)
    println(io, " "^pad_head * s_head)
    println(io, " "^center_pos * "^")
    println(io, line3)
    println(io, " "^pad_tail * s_tail)
end

mutable struct PolyhedralMesh
    vertices::Vector{Vertex}
    primal_edges::Vector{PrimalEdge}
    dual_edges::Vector{DualEdge}
    faces::Vector{Face}
    holes::Vector{Face} # convention: holes are faces with negative id

    PolyhedralMesh() = new([], [], [], [], [])
end

############### Vertex functions
# constructor. Construct a new vertex with id, label, edge, coord and mesh. Vertex will be tail of edge. If tail of edge is not nothing, an error is thrown
function Vertex(id::Int; label::Union{Nothing,Int,String}=nothing, edge::Union{Nothing,Edge}=nothing, coord::Union{Nothing,Vector{Union{AbstractAlgebra.RingElem,Number}}}=nothing, mesh::Union{Nothing,PolyhedralMesh}=nothing)
    v = Vertex()

    v.id = id

    if !isnothing(label)
        v.label = id
    end

    if !isnothing(mesh)
        v.mesh = mesh
    end

    if !isnothing(coord)
        v.coord = coord
    end

    if !isnothing(edge)
        if !isdefined(edge, :tail) # if edge tail is not yet set, set edge tail to the newly created vertex
            tail!(edge, v)
            edge!(v, edge)
        else
            throw(ArgumentError("Tail of edge needs to be nothing."))
        end
    end
    return v
end

id(v::Vertex) = v.id
label(v::Vertex) = v.label
mesh(v::Vertex) = v.mesh::PolyhedralMesh
edge(v::Vertex) = v.edge::PrimalEdge

"""
    edge!(v::Vertex, e::Edge)

Set the edge pointer of `v` to the edge `e`. `e` has to have `v` as its tail.
"""
function edge!(v::Vertex, e::Edge)
    if !(tail(e) === v)
        throw(ArgumentError("tail of edge must be the same object as v ($v), but tail is $(tail(e))"))
    end
    if !(mesh(v) === mesh(e))
        throw(ArgumentError("Meshes of vertex and edge do not match."))
    end
    v.edge = e
    return v
end

function label!(v::Vertex, s::String)
    v.label = s
    return v
end

function label!(v::Vertex, i::Int)
    if id(v) != i
        throw(ArgumentError("If label is integer, vertex id and label must match. Vertex id is $(id(v))."))
    end
    v.label = i
    return v
end

# comparison
Base.:(==)(v1::Vertex, v2::Vertex) = id(v1) == id(v2) && v1.mesh === v2.mesh

############### Face functions
# constructor. Construct a new face with id, edge and mesh. If tail of edge is not nothing, an error is thrown.
function Face(id::Int; edge::Union{Nothing,Edge}=nothing, mesh::Union{Nothing,PolyhedralMesh}=nothing)
    f = Face()

    f.id = id

    if !isnothing(mesh)
        f.mesh = mesh
    end

    if !isnothing(edge)
        if isnothing(tail(edge)) # if edge tail is not yet set, set edge tail to the newly created vertex
            tail!(edge, f)
        else
            throw(ArgumentError("Tail of edge needs to be nothing."))
        end
        edge!(f, edge)
    end
    return f
end

id(f::Face) = f.id
edge(f::Face) = f.edge::Edge
mesh(f::Face) = f.mesh::PolyhedralMesh

"""
    edge!(f::Face, e::Edge)

Set the edge pointer of `f` to the edge `e`. `e` has to have `f` as its tail.
"""
function edge!(f::Face, e::Edge)
    if !(tail(e) === f)
        throw(ArgumentError("Tail of edge must be the same object as f ($f), but tail is $(tail(e))"))
    end
    if !(mesh(f) === mesh(e))
        throw(ArgumentError("Meshes of face and edge do not match."))
    end
    f.edge = e
    return f
end

function vertices(f::Face)
    e = edge(f)

end


############### Edge functions
function PrimalEdge(mesh::PolyhedralMesh)
    e = PrimalEdge()
    e.mesh = mesh

    return e
end

function DualEdge(mesh::PolyhedralMesh)
    e = DualEdge()
    e.mesh = mesh

    return e
end

mesh(e::Edge) = e.mesh::PolyhedralMesh

"""
    tail(e::Edge)

Return the tail of the edge `e`.
"""
tail(e::Edge) = e.tail

"""
    next(e::Edge)

Return the next edge in counterclockwise order around `tail(e)`.
"""
next(e::Edge) = e.next

rot(e::PrimalEdge) = e.rot::DualEdge

rot(e::DualEdge) = e.rot::PrimalEdge

"""
    rot(e::Edge)

Return the edge `e` rotated by 90 degrees counterclockwise. This is the edge that has right as its tail, left as its head, tail as its left and head as its right.
"""
rot(e::Edge) = rot(e)

"""
    flip(e::Edge)

Return the flipped edge of `e`. This is the edge with left and right swapped.
"""
flip(e::Edge) = e.flip

"""
    flip(e::Edge)

Return the reversed edge of `e`. This is the edge with head and tail swapped.
"""
rev(e::Edge) = flip(rot(rot(e)))

"""
    head(e::Edge)

Return the head of the edge `e`.
"""
head(e::Edge) = tail(rot(rot(e)))

"""
    right(e::Edge)

Return the object to the right of the edge `e`.
"""
right(e::Edge) = tail(rot(e))

"""
    left(e::Edge)

Return the object to the left of the edge `e`.
"""
left(e::Edge) = head(rot(e))

"""
    prev(e::Edge)

Return the previous edge in counterclockwise order around `tail(e)`.
"""
prev(e::Edge) = flip(next(flip(e)))

"""
    tail!(e::PrimalEdge, v::Vertex)

Set the tail of `e` to `v` and return the updated edge `e`
"""
function tail!(e::PrimalEdge, v::Vertex)
    if !(e.mesh === v.mesh)
        throw(ArgumentError("Meshes of Edge and face/vertex don't match."))
    end

    e.tail = v
    e.flip.tail = v

    return e
end

"""
    tail!(e::DualEdge, f::Face)

Set the tail of `e` to `f` and return the updated edge `e`
"""
function tail!(e::DualEdge, f::Face)
    if !(e.mesh === f.mesh)
        throw(ArgumentError("Meshes of Edge and face/vertex don't match."))
    end

    e.tail = f
    e.flip.tail = f

    return e
end

"""
    head!(e::PrimalEdge, v::Vertex)

Set the head of `e` to `f` and return the updated edge `e`.
"""
head!(e::PrimalEdge, v::Vertex) = tail!(rot(rot(e)), v)

"""
    head!(e::DualEdge, f::Face)

Set the head of `e` to `f` and return the updated edge `e`.
"""
head!(e::DualEdge, f::Face) = tail!(rot(rot(e)), f)

"""
    left!(e::PrimalEdge, f::Face)

Set left of `e` to `f` and return the updated edge `e`.
"""
left!(e::PrimalEdge, f::Face) = tail!(rot(rot(rot(e))), f)

"""
    left!(e::DualEdge, v::Vertex)

Set left of `e` to `v` and return the updated edge `e`.
"""
left!(e::DualEdge, v::Vertex) = tail!(rot(rot(rot(e))), v)

"""
    right!(e::PrimalEdge, f::Face)

Set right of `e` to `f` and return the updated edge `e`.
"""
right!(e::PrimalEdge, f::Face) = tail!(rot(e), f)

"""
    right!(e::PrimalEdge, v::Vertex)

Set right of `e` to `v` and return the updated edge `e`.
"""
right!(e::PrimalEdge, v::Vertex) = tail!(rot(e), v)

"""
    is_boundary(e::PrimalEdge)

Return whether `e` is a boundary edge, i.e. it has only one adjacent face. By convention, this means the `id` of its left xor right face is <= 0.
"""
is_boundary(e::PrimalEdge) = xor(id(left(e)) <= 0, id(right(e)) <= 0)


"""
    make_edge()

Construct an empty edge with all its pointers to other edges correctly set up.
"""
function make_edge()
    e = PrimalEdge() # Edge e
    e_rot = DualEdge() # e rotated by 90 degrees ccw
    e_rot2 = PrimalEdge() # e rotated by 180 degrees
    e_rot3 = DualEdge() # e rotated by -90 degrees ccw

    e.rot, e_rot.rot, e_rot2.rot, e_rot3.rot = e_rot, e_rot2, e_rot3, e # rotation pointers
    e.next, e_rot.next, e_rot2.next, e_rot3.next = e, e_rot3, e_rot2, e_rot # a valid quad edge respects e.rot.next.rot.next === e

    # flipped edge (left and right swapped)
    e_flip = PrimalEdge()
    e_flip_rot = DualEdge() # e_flip rotated by 90 degrees ccw
    e_flip_rot2 = PrimalEdge() # e_flip rotated by 180 degrees
    e_flip_rot3 = DualEdge() # e_flip rotated by 270 ccw

    # e_flip.tail, e_flip_rot.tail, e_flip_rot2.tail, e_flip_rot3.tail = e.tail, e_rot3.tail, e_rot2.tail, e_rot.tail   # as edges are empty, e.tail, ... are initialized to nothing. Is not set up as pointers. Use tail! to set tail.
    e.flip, e_rot.flip, e_rot2.flip, e_rot3.flip = e_flip, e_flip_rot3, e_flip_rot2, e_flip_rot
    e_flip.rot, e_flip_rot.rot, e_flip_rot2.rot, e_flip_rot3.rot = e_flip_rot, e_flip_rot2, e_flip_rot3, e_flip
    e_flip.next, e_flip_rot.next, e_flip_rot2.next, e_flip_rot3.next = e_flip, e_flip_rot3, e_flip_rot2, e_flip_rot
    e_flip.flip, e_flip_rot.flip, e_flip_rot2.flip, e_flip_rot3.flip = e, e_rot3, e_rot2, e_rot

    return e
end

"""
    make_edge(mesh::PolyhedralMesh)

Construct an empty edge in the polyhedral mesh `mesh` with all its pointers to other edges correctly set up.
"""
function make_edge(mesh::PolyhedralMesh)
    e = make_edge()
    e.mesh = mesh
    e_rot = rot(e)
    e_rot.mesh = mesh
    e_rot2 = rot(rot(e))
    e_rot2.mesh = mesh
    e_rot3 = rot(rot(rot(e)))
    e_rot3.mesh = mesh
    e_flip = flip(e)
    e_flip.mesh = mesh
    e_flip_rot = rot(flip(e))
    e_flip_rot.mesh = mesh
    e_flip_rot2 = rot(rot(flip(e)))
    e_flip_rot2.mesh = mesh
    e_flip_rot3 = rot(rot(rot(flip(e))))
    e_flip_rot3.mesh = mesh

    return e
end

"""
    splice!(a::PrimalEdge, b::PrimalEdge)

Perform the splice operation on the two edges `a` and `b` and return the updated edges. 
If `a` and `b` belong to different vertex rings, this operation joins the two vertex rings by topologically identifying the tails and joining the vertex rings.
If `a` and `b` belong to the same vertex ring (have the same tail), this operation topologically splits the tail vertex into two vertices with corresponding vertex rings.
splice! is an involution, meaning splice!(splice!(a, b)...) === a, b
"""
function splice!(a::PrimalEdge, b::PrimalEdge)
    # if a and b are part of different vertex cycles, splice pastes these vertex cycles in a way that new next of a is original next of b and vice versa. 
    # consistency is ensured by also the next pointers of rot(next(a)) and rot(next(b)).
    # for the flipped case the next pointers of flip(prev(a)) and flip(prev(b)) need to be switched. Further the next pointers of rot(next(flip(prev(a)))) and rot(next(flip(prev(b)))) need to be switched to ensure consistency.
    # splice! is self-inverse, meaning calling splice on two edges of the same cycle results in "splitting" the tail and two disconnected cycles.

    # edges that need their next pointers altered.
    a_next_rot = rot(next(a)) # dual edge corr to a, non flipped
    b_next_rot = rot(next(b)) # dual edge corr to b, non flipped
    alpha = flip(next(a)) # primal edge corr to a, flipped
    beta = flip(next(b)) # primal edge corr to b, flipped
    alpha_next_rot = rot(next(alpha)) # dual edge corr to a, flipped
    beta_next_rot = rot(next(beta)) # dual edge corr to b, flipped

    # their next pointers
    a_next = next(a)
    b_next = next(b)
    a_next_rot_next = next(a_next_rot)
    b_next_rot_next = next(b_next_rot)
    alpha_next = next(alpha)
    beta_next = next(beta)
    alpha_next_rot_next = next(alpha_next_rot)
    beta_next_rot_next = next(beta_next_rot)

    # switch next pointers
    a.next = b_next
    b.next = a_next
    a_next_rot.next = b_next_rot_next
    b_next_rot.next = a_next_rot_next
    alpha.next = beta_next
    beta.next = alpha_next
    alpha_next_rot.next = beta_next_rot_next
    beta_next_rot.next = alpha_next_rot_next

    return a, b
end


############### PolyhedralMesh functions
# PolyhedralMesh constructor from list of face vertex lists.
function PolyhedralMesh(faces::AbstractVector{<:AbstractVector{<:Integer}}; labels::AbstractVector{<:Union{Int,String}}=sort(union(faces...)))
    vert_ids = sort(union(faces...))
    if vert_ids != collect(1:length(vert_ids))
        throw(ArgumentError("Vertex ids need to be consecutive integers starting from 1."))
    end

    mesh = PolyhedralMesh() # empty mesh that is updated
    vertices = [Vertex(id, mesh=mesh, label=labels[id]) for id in vert_ids]
    mesh.vertices = vertices

    null_face = Face(0) # default face that is set when left or right is not known - the "outside"/surrounding space of the mesh
    null_face.mesh = mesh

    edge_dict = Dict{SVector{2,Int},PrimalEdge}() # dict that keeps track of which primary edges are already in the mesh. Convention: All 

    for (i, f) in enumerate(faces)
        face = Face(i, mesh=mesh)
        push!(mesh.faces, face)
        prev_edge = nothing

        for (it, j) in enumerate(vcat(eachindex(f), [1]))
            tail_id = f[j]
            head_id = f[mod1(j + 1, length(f))]
            edge_id = [tail_id, head_id]

            if !haskey(edge_dict, edge_id)
                # case 1: edge is not already processed. Construct new edge. Set tail, head and right face.
                current_edge = make_edge(mesh)
                tail!(current_edge, vertices[tail_id])
                head!(current_edge, vertices[head_id])
                right!(current_edge, face)
                left!(current_edge, null_face)

                # add edge and all its related edges to mesh.edges
                mesh.primal_edges = vcat(mesh.primal_edges, [current_edge,
                    rot(rot(current_edge)),
                    flip(current_edge),
                    rot(rot(flip(current_edge)))
                ])

                mesh.dual_edges = vcat(mesh.dual_edges, [rot(current_edge),
                    rot(rot(rot(current_edge))),
                    rot(flip(current_edge)),
                    rot(rot(rot(flip(current_edge))))
                ])

                # set edge entries of tail vertex and right face, if it isn't already set
                if !isdefined(vertices[tail_id], :edge)
                    edge!(vertices[tail_id], current_edge)
                end
                if !isdefined(face, :edge)
                    edge!(face, rot(current_edge))
                end

            elseif haskey(edge_dict, edge_id) && it != length(f) + 1
                # case 2: edge was already processed in another face. By behaviour defined above, the right side is set to this different face
                current_edge = flip(edge_dict[edge_id]) # flip, so that left face is set and right set is null face

                if right(current_edge) !== null_face # check consistency of face list. If right is not the null face, there is an edge with at least three incident faces.
                    throw(ArgumentError("Face list is not consistent. There is at least one edge with three or more faces incident to it."))
                end

                right!(current_edge, face)
            else
                # case 3: edge is the first edge around the current face f.
                current_edge = edge_dict[edge_id] # right side is set. Either left or right is current face.
                if right(current_edge) !== face
                    current_edge = flip(current_edge) # flip edge, if right is not current face.
                end
            end

            # splice prev(edge) with rot(rot(prev_edge)). 
            # prev_edge is the last edge around the face, thus rot(rot(prev_edge)) has same tail as edge and face as left face. Thus next(rot(rot(prev_edge))) should be edge. This is achieved by splice!(prev(e), rot(rot(prev_edge)))
            if !isnothing(prev_edge) && next(rot(rot(prev_edge))) !== current_edge
                splice!(prev(current_edge), rot(rot(prev_edge)))
            end


            prev_edge = current_edge # update prev_edge. Note that prev_edge always has f to its right.

            # update edge_dict
            edge_dict[edge_id] = current_edge
            edge_dict[reverse(edge_id)] = flip(rot(rot(current_edge))) # flip so that right side of all edge_dict entries are set

            next_id = [id(tail(next(current_edge))), id(head(next(current_edge)))]
            if right(next(current_edge)) !== null_face
                edge_dict[next_id] = next(current_edge) # next(e) has right face set
                edge_dict[reverse(next_id)] = rev(next(current_edge)) # ensure that all entries in edge_dict have right face set.
            else
                edge_dict[next_id] = flip(next(current_edge)) # next(e) has left face set
                edge_dict[reverse(next_id)] = rev(flip(next(current_edge))) # ensure that all entries in edge_dict have right face set.
            end

            prev_id = [id(tail(prev(current_edge))), id(head(prev(current_edge)))]
            if right(prev(current_edge)) !== null_face
                edge_dict[prev_id] = prev(current_edge) # next(e) has right face set
                edge_dict[reverse(prev_id)] = rev(prev(current_edge)) # ensure that all entries in edge_dict have right face set.
            else
                edge_dict[prev_id] = flip(prev(current_edge)) # next(e) has left face set
                edge_dict[reverse(prev_id)] = rev(flip(prev(current_edge))) # ensure that all entries in edge_dict have right face set.
            end
        end
    end

    return mesh
end

vertices(mesh::PolyhedralMesh) = mesh.vertices
primal_edges(mesh::PolyhedralMesh) = mesh.primal_edges
dual_edges(mesh::PolyhedralMesh) = mesh.dual_edges
faces(mesh::PolyhedralMesh) = mesh.faces

function holes(mesh::PolyhedralMesh)

end

