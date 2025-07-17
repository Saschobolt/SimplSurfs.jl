module SimplSurfs

using Graphs
using AbstractAlgebra

export Vertex, Face, Edge, PolyhedralMesh, make_edge, head, tail, left, right, flip, rot, next, mesh, splice!
include("PolyhedralMesh.jl")

end # module SimplSurfs
