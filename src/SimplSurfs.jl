module SimplSurfs

using Graphs
using AbstractAlgebra
using StaticArrays

export Vertex, Face, Edge, PolyhedralMesh, make_edge, head, head!, tail, tail!, left, left!, right, right!, flip, rot, next, mesh, splice!, edge, is_primary, is_dual, splice!, vertices, faces, edges, prev
include("PolyhedralMesh.jl")

end # module SimplSurfs
