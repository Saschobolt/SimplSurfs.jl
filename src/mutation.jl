"""
    splice_out!(e::PrimalEdge)

Splice out `e` from the edge ring around its tail.
"""
splice_out!(e::PrimalEdge) = splice!(e, prev(e))

"""
    edge_turn!(surf::SimplicialSurface, e::PrimalEdge)

Perform an edge turn (edge flip) on the edge `e` of the simplicial surface `surf`.
This replaces the edge `e` (shared by two triangles) with the other diagonal of the quadrilateral formed by their union.
Returns the updated simplicial surface.
"""
function edge_turn!(surf::SimplicialSurface, e::PrimalEdge)
    m = mesh(surf)
    @assert !is_boundary(e) "Cannot perform edge turn on boundary edge."
    @assert left(e) !== right(e) "Cannot perform edge turn on edge with identical left and right face."
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
        v1 = tail(edge)
        v2 = head(edge)
        rev_flag = id(v1) > id(v2)
        if rev_flag
            m.primal_edges[(id(v2), id(v1))] = edge
        else
            m.primal_edges[(id(v1), id(v2))] = edge
        end
    end

    return surf
end

"""
    edge_turn(surf::SimplicialSurface, e::PrimalEdge)

Perform an edge turn (edge flip) on the edge `e` of the simplicial surface `surf`.
This replaces the edge `e` (shared by two triangles) with the other diagonal of the quadrilateral formed by their union.
Returns a new simplicial surface with the edge turned.

See also `edge_turn!`.
"""
function edge_turn(surf::SimplicialSurface, e::PrimalEdge)
    surf_copy = deepcopy(surf)
    e_key = (id(tail(e)), id(head(e)))
    rev_flag = id(tail(e)) > id(head(e))
    if rev_flag
        e_copy = rev(edges(surf_copy)[reverse(e_key)])
    else
        e_copy = edges(surf_copy)[e_key]
    end
    return edge_turn!(surf_copy, e_copy)
end