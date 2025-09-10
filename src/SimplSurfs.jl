module SimplSurfs

using Graphs
using AbstractAlgebra
using StaticArrays
using RigidityTheoryTools
using GeometryBasics
import GLMakie

export Vertex, id, coordinate_matrix, coordinate_matrix, Face, Edge, PolyhedralMesh, labels, make_edge, head, head!, tail, tail!, left, left!, right, right!, flip, rot, invrot, next, lnext, rnext, mesh, splice!, edge, is_primary, is_dual, splice!, vertices, faces, edges, dual_edges, prev, is_boundary, holes!, labels
include("PolyhedralMesh.jl")

export SimplicialSurface
include("SimplicialSurface.jl")

include("symmetry.jl")

export coordinates!, coordinate_matrix, graph
include("rigidity.jl")

export octahedron, octahedron_emb, tetrahedron, double_tetrahedron, cube
include("examples.jl")

include("plotting.jl")

end # module SimplSurfs
