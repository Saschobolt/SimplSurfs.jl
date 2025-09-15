module SimplSurfs

using Graphs
using AbstractAlgebra
using LinearAlgebra
using StaticArrays
using RigidityTheoryTools
using GeometryBasics
import GLMakie

export Vertex, id, faces, neighbors, coordinates, primal_edges, Face, Edge, PolyhedralMesh, labels, make_edge, head, head!, tail, tail!, left, left!, right, right!, flip, rot, invrot, next, lnext, rnext, mesh, splice!, edge, is_primary, is_dual, splice!, vertices, faces, edges, dual_edges, prev, is_boundary, holes!, labels, embedding_dim
include("PolyhedralMesh.jl")

export SimplicialSurface, embedding_dim
include("SimplicialSurface.jl")

include("symmetry.jl")

export coordinates!, coordinate_matrix, graph, distance, distancesq, elength, elengthsq, transform!, transform_iso!, kabsch
include("rigidity.jl")

export octahedron, octahedron_emb, tetrahedron, tetrahedron_emb, double_tetrahedron, cube
include("examples.jl")

include("plotting.jl")

end # module SimplSurfs
