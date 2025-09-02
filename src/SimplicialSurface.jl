mutable struct SimplicialSurface
    mesh::PolyhedralMesh

    function SimplicialSurface(mesh::PolyhedralMesh)
        for f in faces(mesh)
            if length(vertices(f)) != 3
                throw(ArgumentError("All faces in a simplicial surface must be triangles"))
            end
        end
    end
end

function SimplicialSurface(faces::AbstractVector{<:AbstractVector{<:Integer}}; labels::AbstractVector{<:Union{Int,String}}=sort(union(faces...)))
    for f in faces
        if length(f) != 3
            throw(ArgumentError("All faces in a simplicial surface must be triangles."))
        end
    end
    mesh = PolyhedralMesh(faces; labels=labels)
    return SimplicialSurface(mesh)
end

vertices(surf::SimplicialSurface) = vertices(surf.mesh)
labels(surf::SimplicialSurface) = [label(v) for v in vertices(surf)]
edges(surf::SimplicialSurface) = edges(surf.mesh)
faces(surf::SimplicialSurface) = faces(surf.mesh)
holes!(surf::SimplicialSurface) = holes!(surf.mesh)

function Base.show(io::IO, surf::SimplicialSurface)
    print(io, "Simplicial surface with $(length(vertices(surf))) vertices, $(length(edges(surf))) edges, $(length(faces(surf))) faces.")
end
