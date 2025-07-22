############################################################################
# types
############################################################################
mutable struct Vertex
    id::Int
    label::Union{Int,String}
    coord::Union{Nothing,Vector{Union{AbstractAlgebra.RingElem,Number}}}
    edge::Any # Edge    # a primal edge with this vertex as tail
    mesh::Any # PolyhedralMesh
end

function Base.show(io::IO, v::Vertex)
    print(io, "V$(v.id) ($(v.label))")
end

mutable struct Face
    id::Int
    edge::Any # Edge    # a dual edge with this face as tail
    mesh::Any # PolyhedralMesh
end

function Base.show(io::IO, f::Face)
    print(io, "F $(f.id) ($(string(f.vertices)))")
end

mutable struct Edge
    tail::Union{Nothing,Vertex,Face} # tail of edge

    next::Union{Nothing,Edge} # next edge in counterclockwise order around tail (if tail is a vertex, this is the edge with tail == this.tail and head = other vertex connected to tail in this.left with right = this.left)
    rot::Union{Nothing,Edge} # edge with tail = this.right, head = this.left, left = this.tail, right = this.head
    flip::Union{Nothing,Edge} # edge with tail = this.tail, head = this.head, left = this.right, right = this.left

    mesh::Any # PolyhedralMesh

    Edge() = new(nothing, nothing, nothing, nothing, nothing)
end

function Base.show(io::IO, e::Edge)
    s_head = head(e) === nothing ? "nothing" : string(head(e))
    s_tail = tail(e) === nothing ? "nothing" : string(tail(e))
    s_left = left(e) === nothing ? "nothing" : string(left(e))
    s_right = right(e) === nothing ? "nothing" : string(right(e))

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
    edges::Vector{Edge}
    faces::Vector{Face}

    PolyhedralMesh() = new([], [], [])
end

############### Vertex functions
# constructor. Construct a new vertex with id, label, edge, coord and mesh. Vertex will be tail of edge. If tail of edge is not nothing, an error is thrown
function Vertex(id::Int; label::Union{Nothing,Int,String}=nothing, edge::Union{Nothing,Edge}=nothing, coord::Union{Nothing,Vector{Union{AbstractAlgebra.RingElem,Number}}}=nothing, mesh::Union{Nothing,PolyhedralMesh}=nothing)
    if isnothing(label)
        label = id
    end

    v = Vertex(id, label, coord, nothing, mesh)

    if !isnothing(edge)
        if isnothing(tail(edge)) # if edge tail is not yet set, set edge tail to the newly created vertex
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
edge(v::Vertex) = v.edge::Edge

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
    f = Face(id, nothing, mesh)

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
edge(f::Face) = f.edge
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


############### Edge functions
function Edge(mesh::PolyhedralMesh)
    e = Edge()
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

"""
    rot(e::Edge)

Return the edge `e` rotated by 90 degrees counterclockwise. This is the edge that has right as its tail, left as its head, tail as its left and head as its right.
"""
rot(e::Edge) = e.rot

"""
    flip(e::Edge)

Return the flipped edge of `e`. This is the edge with left and right swapped.
"""
flip(e::Edge) = e.flip

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
    is_primary(e::Edge)

Return whether `e`is a primary edge, i.e. its head and tail are vertices.
"""
is_primary(e::Edge) = typeof(head(e)) <: Union{Nothing,Vertex} && typeof(tail(e)) <: Union{Nothing,Vertex} && typeof(left(e)) <: Union{Nothing,Face} && typeof(right(e)) <: Union{Nothing,Face}

"""
    id_dual(e::Edge)

Return whether `e`is a dual edge, i.e. its head and tail are faces.
"""
is_dual(e::Edge) = typeof(head(e)) <: Union{Nothing,Face} && typeof(tail(e)) <: Union{Nothing,Face} && typeof(left(e)) <: Union{Nothing,Vertex} && typeof(right(e)) <: Union{Nothing,Vertex}

"""
    tail!(e::Edge, x::Union{Vertex,Face})

Set the tail of `e` to `x` and return the updated edge `e`.
"""
function tail!(e::Edge, x::Union{Vertex,Face})
    # if edge is non empty (not both primary and dual), check that typeof(x) matches edge type.
    if !is_dual(e) && is_primary(e) && !(typeof(x) <: Vertex)
        throw(ArgumentError("Can't set the tail, as edge is primary, but type of new tail is $(typeof(x)). (Should be Vertex)."))
    elseif !is_primary(e) && is_dual(e) && !(typeof(x) <: Face)
        throw(ArgumentError("Can't set the tail, as edge is dual, but type of new tail is $(typeof(x)). (Should be Face)."))
    end

    if !(e.mesh === x.mesh)
        throw(ArgumentError("Meshes of Edge and face/vertex don't match."))
    end

    e.tail = x
    e.flip.tail = x
    return e
end

"""
    head!(e::Edge, x::Union{Vertex,Face})

Set the head of `e` to `x` and return the updated edge `e`.
"""
function head!(e::Edge, x::Union{Vertex,Face})
    if is_primary(e) && !(typeof(x) <: Vertex)
        throw(ArgumentError("Can't set the head, as edge is primary, but type of new head is $(typeof(x)). (Should be Vertex)."))
    elseif is_dual(e) && !(typeof(x) <: Face)
        throw(ArgumentError("Can't set the head, as edge is dual, but type of new head is $(typeof(x)). (Should be Face)."))
    end

    tail!(rot(rot(e)), x)
    return e
end

"""
    head!(e::Edge, x::Union{Vertex,Face})

Set left of `e` to `x` and return the updated edge `e`.
"""
function left!(e::Edge, x::Union{Vertex,Face})
    if is_primary(e) && !(typeof(x) <: Face)
        throw(ArgumentError("Can't set left, as edge is primary, but type of new left is $(typeof(x)). (Should be Face)."))
    elseif is_dual(e) && !(typeof(x) <: Vertex)
        throw(ArgumentError("Can't set left, as edge is dual, but type of new left is $(typeof(x)). (Should be Vertex)."))
    end

    head!(rot(e), x)
    return e
end

"""
    right!(e::Edge, x::Union{Vertex,Face})

Set right of `e` to `x` and return the updated edge `e`.
"""
function right!(e::Edge, x::Union{Vertex,Face})
    if is_primary(e) && !(typeof(x) <: Face)
        throw(ArgumentError("Can't set right, as edge is primary, but type of new right is $(typeof(x)). (Should be Face)."))
    elseif is_dual(e) && !(typeof(x) <: Vertex)
        throw(ArgumentError("Can't set right, as edge is dual, but type of new right is $(typeof(x)). (Should be Vertex)."))
    end

    tail!(rot(e), x)
    return e
end


"""
    make_edge()

Construct an empty edge with all its pointers to other edges correctly set up.
"""
function make_edge()
    e = Edge() # Edge e
    e_rot = Edge() # e rotated by 90 degrees ccw
    e_rot2 = Edge() # e rotated by 180 degrees
    e_rot3 = Edge() # e rotated by -90 degrees ccw

    e.rot, e_rot.rot, e_rot2.rot, e_rot3.rot = e_rot, e_rot2, e_rot3, e # rotation pointers
    e.next, e_rot.next, e_rot2.next, e_rot3.next = e, e_rot3, e_rot2, e_rot # a valid quad edge respects e.rot.next.rot.next === e

    # flipped edge (left and right swapped)
    e_flip = Edge()
    e_flip_rot = Edge() # e_flip rotated by 90 degrees ccw
    e_flip_rot2 = Edge() # e_flip rotated by 180 degrees
    e_flip_rot3 = Edge() # e_flip rotated by 270 ccw

    # e_flip.tail, e_flip_rot.tail, e_flip_rot2.tail, e_flip_rot3.tail = e.tail, e_rot3.tail, e_rot2.tail, e_rot.tail   # as edges are empty, e.tail, ... are initialized to nothing. Is not set up as pointers. Use tail! to set tail.
    e.flip, e_rot.flip, e_rot2.flip, e_rot3.flip = e_flip, e_flip_rot3, e_flip_rot2, e_flip_rot
    e_flip.rot, e_flip_rot.rot, e_flip_rot2.rot, e_flip_rot3.rot = e_flip_rot, e_flip_rot2, e_flip_rot3, e_flip
    e_flip.next, e_flip_rot.next, e_flip_rot2.next, e_flip_rot3.next = e_flip, e_flip_rot3, e_flip_rot2, e_flip_rot
    e_flip.flip, e_flip_rot.flip, e_flip_rot2.flip, e_flip_rot3.flip = e, e_rot3.flip, e_rot2.flip, e_rot.flip

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
    splice!(a::Edge, b::Edge)

Perform the splice operation on the two edges `a` and `b` and return the updated edges. 
If `a` and `b` belong to different vertex rings, this operation joins the two vertex rings by topologically identifying the tails and joining the vertex rings.
If `a` and `b` belong to the same vertex ring (have the same tail), this operation topologically splits the tail vertex into two vertices with corresponding vertex rings.
splice! is an involution, meaning splice!(splice!(a, b)...) === a, b
"""
function splice!(a::Edge, b::Edge)
    # non flipped
    alpha = next(rot(a))
    beta = next(rot(b))

    a_next = next(a)
    b_next = next(b)
    a.next = b_next
    b.next = a_next

    alpha_next = next(alpha)
    beta_next = next(beta)
    alpha.next = beta_next
    beta.next = alpha_next

    # flipped
    a_flip = flip(a)
    b_flip = flip(b)
    alpha_flip = next(rot(flip(a)))
    beta_flip = next(rot(flip(b)))

    a_flip_next = next(flip(a))
    b_flip_next = next(flip(b))
    a_flip.next = b_flip_next
    b_flip.next = a_flip_next

    alpha_flip_next = next(alpha_flip)
    beta_flip_next = next(beta_flip)
    alpha_flip.next = beta_flip_next
    beta_flip.next = alpha_flip_next

    return a, b
end


############### PolyhedralMesh functions
vertices(mesh) = mesh.vertices
edges(mesh) = mesh.edges
faces(mesh) = mesh.faces

"""
    _check_consistency(mesh::PolyhedralMesh)

Check the consistency of the mesh, i.e. for each edge `e`: 
- tail(next(e)) === tail(e) 
- left(edge) === right(next(e))
- tail(flip(e)) === tail(e)
- !is_primary(e) || (!isnothing(tail(e)) && !isnothing(head(e)) && (!isnothing(left(e)) || !isnothing(right(e))))
- !is_dual(e) || (!isnothing(left(e)) && !isnothing(right(e)) && (!isnothing(head(e)) || !isnothing(tail(e))))
"""
function _check_consistency(mesh::PolyhedralMesh)
    for e in edges(mesh)
        if !(tail(next(e)) === tail(e))
            return false
        elseif !(left(edge) === right(next(e)))
            return false
        elseif !(tail(flip(e)) === tail(e))
            return false
        elseif !(!is_primary(e) || (!isnothing(tail(e)) && !isnothing(head(e)) && (!isnothing(left(e)) || !isnothing(right(e)))))
            return false
        elseif !(!is_dual(e) || (!isnothing(left(e)) && !isnothing(right(e)) && (!isnothing(head(e)) || !isnothing(tail(e)))))
            return false
        end
    end

    return true
end