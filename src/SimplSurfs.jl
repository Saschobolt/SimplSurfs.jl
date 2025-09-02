module SimplSurfs

using Graphs
using AbstractAlgebra
using StaticArrays

export Vertex, Face, Edge, PolyhedralMesh, make_edge, head, head!, tail, tail!, left, left!, right, right!, flip, rot, invrot, next, lnext, rnext, mesh, splice!, edge, is_primary, is_dual, splice!, vertices, faces, edges, dual_edges, prev, is_boundary, holes!, labels
include("PolyhedralMesh.jl")

export SimplicialSurface
include("SimplicialSurface.jl")

include("symmetry.jl")

include("rigidity.jl")

export octahedron, tetrahedron, double_tetrahedron
include("examples.jl")

end # module SimplSurfs
