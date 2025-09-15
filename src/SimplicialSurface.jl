mutable struct SimplicialSurface{PositionDim,PositionType}
    mesh::PolyhedralMesh{PositionDim,PositionType}

    function SimplicialSurface{PositionDim,PositionType}(mesh::PolyhedralMesh{PositionDim,PositionType}) where {PositionDim,PositionType}
        for f in faces(mesh)
            if length(vertices(f)) != 3
                throw(ArgumentError("All faces in a simplicial surface must be triangles"))
            end
        end

        return new{PositionDim,PositionType}(mesh)
    end
end

SimplicialSurface(mesh::PolyhedralMesh{PositionDim,PositionType}) where {PositionDim,PositionType} = SimplicialSurface{PositionDim,PositionType}(mesh)

function SimplicialSurface{PositionDim,PositionType}(faces::AbstractVector{<:AbstractVector{<:Integer}}; labels::AbstractVector{<:Union{Int,String}}=sort(union(faces...)), positions::Union{AbstractVector{Point{PositionDim}},Nothing}=nothing) where {PositionDim,PositionType}
    for f in faces
        if length(f) != 3
            throw(ArgumentError("All faces in a simplicial surface must be triangles."))
        end
    end
    mesh = PolyhedralMesh{PositionDim,PositionType}(faces; labels=labels, positions=positions)
    return SimplicialSurface{PositionDim,PositionType}(mesh)
end

SimplicialSurface(faces::AbstractVector{<:AbstractVector{<:Integer}}, positions::AbstractVector{Point{PositionDim,PositionType}}; labels::AbstractVector{<:Union{Int,String}}=sort(union(faces...))) where {PositionDim,PositionType} = SimplicialSurface{PositionDim,PositionType}(faces, labels=labels, positions=positions)
SimplicialSurface(faces::AbstractVector{<:AbstractVector{<:Integer}}; labels::AbstractVector{<:Union{Int,String}}=sort(union(faces...))) = SimplicialSurface{0,Any}(faces, labels=labels)

mesh(surf::SimplicialSurface) = surf.mesh

vertices(surf::SimplicialSurface) = vertices(surf.mesh)
labels(surf::SimplicialSurface) = [label(v) for v in vertices(surf)]
edges(surf::SimplicialSurface) = edges(surf.mesh)
faces(surf::SimplicialSurface) = faces(surf.mesh)
holes!(surf::SimplicialSurface) = holes!(surf.mesh)
embedding_dim(surf::SimplicialSurface) = embedding_dim(surf.mesh)

function Base.show(io::IO, surf::SimplicialSurface)
    print(io, "Simplicial surface with $(length(vertices(surf))) vertices, $(length(edges(surf))) edges, $(length(faces(surf))) faces.")
end
