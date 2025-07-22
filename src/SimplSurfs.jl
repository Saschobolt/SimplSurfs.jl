module SimplSurfs

using Graphs
using AbstractAlgebra

export Vertex, Face, Edge, PolyhedralMesh, make_edge, head, head!, tail, tail!, left, left!, right, right!, flip, rot, next, mesh, splice!, edge, is_primary, is_dual, splice!
include("PolyhedralMesh.jl")

end # module SimplSurfs
