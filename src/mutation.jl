"""
    splice_out!(e::PrimalEdge)

Splice out `e` from the edge ring around its tail.
"""
splice_out!(e::PrimalEdge) = splice!(e, prev(e))

"""
    update_edge_dict!(e::PrimalEdge)

Update the entry for edge `e` in the primal edge dictionary of the mesh.
"""
function update_edge_dict!(e::PrimalEdge)
    edge_id = (id(tail(e)), id(head(e)))
    rev_flag = id(tail(e)) > id(head(e))
    m = mesh(e)
    if rev_flag
        m.primal_edges[(edge_id[2], edge_id[1])] = rev(e)
    else
        m.primal_edges[edge_id] = e
    end
end

"""
    delete_edge_dict_entry!(e::PrimalEdge)

Delete the entry for edge `e` from the primal edge dictionary of the mesh.
"""
function delete_edge_dict_entry!(e::PrimalEdge)
    edge_id = (id(tail(e)), id(head(e)))
    rev_flag = id(tail(e)) > id(head(e))
    m = mesh(e)
    if rev_flag
        delete!(m.primal_edges, (edge_id[2], edge_id[1]))
    else
        delete!(m.primal_edges, edge_id)
    end
end

"""
    edge_turn!(e::PrimalEdge)

Perform an edge turn (edge flip) on the edge `e`.
This replaces the edge `e` (shared by two triangles) with the other diagonal of the quadrilateral formed by their union.
Returns the updated mesh.
"""
function edge_turn!(e::PrimalEdge)
    m = mesh(e)
    @assert !is_boundary(e) "Cannot perform edge turn on boundary edge."
    @assert left(e) !== right(e) "Cannot perform edge turn on edge with identical left and right face."
    @assert length(vertices(left(e))) == 3 "Left face of edge is not a triangle."
    @assert length(vertices(right(e))) == 3 "Right face of edge is not a triangle."
    f1 = left(e) # left face
    tip_left = setdiff(vertices(f1), (tail(e), head(e))) |> first # vertex opposite to e in left face
    f2 = right(e) # right face
    tip_right = setdiff(vertices(f2), (tail(e), head(e))) |> first # vertex opposite to e in right face
    v1 = tail(e) # one endpoint of e
    v2 = head(e) # other endpoint of e
    v1_tip_left = next(e) # edge from v1 to tip_left
    v1_tip_right = prev(e) # edge from v1 to tip_right
    v2_tip_left = prev(rot(rot(e))) # edge from v2 to tip_left
    v2_tip_right = next(rot(rot(e))) # edge from v2 to tip_right

    # update edge reference of vertices in case their edge reference was e. This is done, because e will be removed from the mesh.
    if edge(v1) === e || edge(v1) === flip(e)
        edge!(v1, v1_tip_right)
    end
    if edge(v2) === e || edge(v2) === flip(e)
        edge!(v2, v2_tip_left)
    end

    # update edge reference of faces in case their edge reference was e. This is done, because e will be removed from the mesh.
    if rot(edge(f1)) === e || rot(edge(f1)) === rev(e)
        edge!(f1, rot(v1_tip_left))
    end
    if rot(edge(f2)) === flip(e) || rot(edge(f2)) === flip(rev(e))
        edge!(f2, rot(v2_tip_right))
    end

    # splice out e from the edge rings around v1 and v2
    splice_out!(e) # around v1
    splice_out!(rot(rot(e))) # around v2

    # create new edge from tip_left to tip_right
    new_e = make_edge(m)
    tail!(new_e, tip_left)
    head!(new_e, tip_right)
    right!(new_e, f1)
    left!(new_e, f2)

    # splice new_e into the edge rings around tip_left and tip_right
    splice!(new_e, rot(rot(v1_tip_left))) # around tip_left: v1_tip_left has f1 to its right, so rot(rot(v1_tip_left)) has f1 to its left and its next should be new_e with f1 to its right.
    splice!(rot(rot(new_e)), rot(rot(v2_tip_right))) # around tip_right: v2_tip_right has f2 to its right, so rot(rot(v2_tip_right)) has f2 to its left and its next should be rot(rot(new_e)) with f2 to its right.

    # update left and right faces of the surrounding edges
    left!(v1_tip_right, f1)
    left!(v2_tip_left, f2)

    # update edge dict of mesh. Keys of edge dict are sorted!
    rev_flag = id(v1) > id(v2)
    rev_flag ? delete!(m.primal_edges, (id(v2), id(v1))) : delete!(m.primal_edges, (id(v1), id(v2)))

    for edge in (new_e, v1_tip_left, v1_tip_right, v2_tip_left, v2_tip_right)
        update_edge_dict!(edge)
    end

    return m
end

"""
    edge_turn(e::PrimalEdge)

Perform an edge turn (edge flip) on the edge `e`.
This replaces the edge `e` (shared by two triangles) with the other diagonal of the quadrilateral formed by their union.
Returns a new simplicial surface with the edge turned.

See also `edge_turn!`.
"""
function edge_turn(e::PrimalEdge)
    m_copy = deepcopy(mesh(e))
    e_key = (id(tail(e)), id(head(e)))
    rev_flag = id(tail(e)) > id(head(e))
    if rev_flag
        e_copy = rev(edges(m_copy)[reverse(e_key)])
    else
        e_copy = edges(m_copy)[e_key]
    end
    return edge_turn!(e_copy)
end

"""
    vertex_split!(e1::PrimalEdge, e2::PrimalEdge; new_vertex_coords::Point=nothing, new_vertex_label=nothing)

Perform a vertex split on the two edges `e1` and `e2` that share the same tail vertex.
This adds a new vertex and two new faces to the mesh. `e2` will keep its tail vertex, while `e1` will get the new vertex as its tail. 
"""
function vertex_split!(e1::PrimalEdge, e2::PrimalEdge; new_vertex_coords=nothing, new_vertex_label=nothing)
    # check that e1 and e2 have the same mesh
    m = mesh(e1)
    if mesh(e2) !== m
        throw(ArgumentError("Edges must belong to the same mesh."))
    end

    # check that e1 and e2 share the same tail
    v = tail(e1)
    if tail(e2) !== v
        throw(ArgumentError("Edges must share the same tail vertex."))
    end

    # check that e1 and e2 are part of an edge ring around their tail.
    next_e = next(e1)
    while true
        if next_e === e2
            break
        elseif next_e === flip(e2)
            throw(ArgumentError("Edges must be part of an edge ring around their tail vertex. Try `flip(e2)` instead of `e2`."))
        end
        next_e = next(next_e)
    end

    # splice e1 and e2
    splice!(e1, e2)

    # save edges around v for which the entries in the primal_edges dict of the mesh need to be updated in a later step
    edges_to_update = primal_edges(v)

    # create new vertex and set its edge reference to e1
    w = Vertex(length(vertices(m)) + 1, label=new_vertex_label, edge=e1, position=new_vertex_coords, mesh=m)
    m.vertices[w.id] = w
    # set tail of e1 and all edges around w to w
    for e in primal_edges(w)
        tail!(e, w)
    end

    # add new edge from v to w
    edge_vw = make_edge(m)
    tail!(edge_vw, v)
    head!(edge_vw, w)
    splice!(e1, rot(rot(edge_vw)))
    splice!(edge_vw, e2)

    # add edge from v to head(e1)
    edge_v_head_e1 = make_edge(m)
    tail!(edge_v_head_e1, v)
    head!(edge_v_head_e1, head(e1))
    splice!(edge_v_head_e1, edge_vw)

    # add edge from w to head(e2)
    edge_w_head_e2 = make_edge(m)
    tail!(edge_w_head_e2, w)
    head!(edge_w_head_e2, head(e2))
    splice!(edge_w_head_e2, rot(rot(edge_vw)))

    # add two new faces
    f1 = Face(length(faces(m)) + 1, edge=invrot(edge_vw), mesh=m)
    f2 = Face(length(faces(m)) + 2, edge=rot(edge_vw), mesh=m)

    # update left face of the primal edges around f1
    for e in primal_edges(f1)
        left!(e, f1)
    end

    # update left face of the primal edges around f2
    for e in primal_edges(f2)
        left!(e, f2)
    end

    # update edge dict
    update_edge_dict!(edge_vw)
    update_edge_dict!(edge_v_head_e1)
    update_edge_dict!(edge_w_head_e2)
    for e in edges_to_update
        update_edge_dict!(e)
    end

    return m
end


"""
    vertex_split(e1::PrimalEdge, e2::PrimalEdge; new_vertex_coords::Point=nothing, new_vertex_label=nothing)

Perform a vertex split on the two edges `e1` and `e2` that share the same tail vertex.
This adds a new vertex and two new faces to the mesh. `e2` will keep its tail vertex, while `e1` will get the new vertex as its tail. 
Returns a new mesh with the vertex split performed.
"""
function vertex_split(e1::PrimalEdge, e2::PrimalEdge; new_vertex_coords=nothing, new_vertex_label=nothing)
    m_copy = deepcopy(mesh(e1))
    e1_key = (id(tail(e1)), id(head(e1)))
    e2_key = (id(tail(e2)), id(head(e2)))
    rev_flag1 = id(tail(e1)) > id(head(e1))
    rev_flag2 = id(tail(e2)) > id(head(e2))
    if rev_flag1
        e1_copy = rev(edges(m_copy)[reverse(e1_key)])
    else
        e1_copy = edges(m_copy)[e1_key]
    end
    if rev_flag2
        e2_copy = rev(edges(m_copy)[reverse(e2_key)])
    else
        e2_copy = edges(m_copy)[e2_key]
    end
    return vertex_split!(e1_copy, e2_copy; new_vertex_coords=new_vertex_coords, new_vertex_label=new_vertex_label)
end


function contraction!(e::PrimalEdge)
end

function contraction(e::PrimalEdge)
    m_copy = deepcopy(mesh(e))
    e_key = (id(tail(e)), id(head(e)))
    rev_flag = id(tail(e)) > id(head(e))
    if rev_flag
        e_copy = rev(edges(m_copy)[reverse(e_key)])
    else
        e_copy = edges(m_copy)[e_key]
    end
    return contraction!(e_copy)
end

function subdivision!(f::Face)
end

function subdivision(f::Face)
    m_copy = deepcopy(mesh(f))
    f_copy = faces(m_copy)[id(f)]
    return subdivision!(f_copy)
end

function flatten!(v::Vertex)
end

function flatten(v::Vertex)
    m_copy = deepcopy(mesh(v))
    v_copy = vertices(m_copy)[id(v)]
    return flatten!(v_copy)
end