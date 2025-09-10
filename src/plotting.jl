function GeometryBasics.mesh(m::PolyhedralMesh; face_colors::AbstractVector{<:GLMakie.Colorant}=[GLMakie.RGBAf(0, 102 / 256, 204 / 256, 1) for f in faces(m)], vertexattributes...)
    ps = coordinates(m)
    fs = [GeometryBasics.NgonFace{length(vertices(f)),Int}(id.(vertices(f))...) for f in faces(m)]
    FT = eltype(fs)
    N = length(fs)
    cs = GeometryBasics.FaceView(face_colors, [FT(i) for i in 1:N])
    ns = GeometryBasics.face_normals(ps, fs)
    return GeometryBasics.mesh(ps, fs, normal=ns, color=cs, vertexattributes...)
end

function plot_mesh(m::PolyhedralMesh; title::String="Polyhedral Mesh", show_labels::Bool=true, kwargs...)
    fig = GLMakie.Figure()
    ax = GLMakie.Axis3(fig[1, 1], aspect=:data, title=title)
    me = GeometryBasics.mesh(m; kwargs...)

    GLMakie.mesh!(ax, me)

    # 2. If requested, add the vertex labels
    if show_labels
        vertex_positions = coordinates(m)
        vertex_labels = [string(label(v)) for v in vertices(m)]

        # Defaults for text labels that can be overridden by kwargs
        text_defaults = (
            fontsize=20,
            color=:black,
            align=(:left, :center),
            offset=(5, 0)
        )

        # Merge defaults with user-provided kwargs. User kwargs take precedence.
        text_kwargs = text_defaults
        GLMakie.text!(ax,
            vertex_positions;
            text=vertex_labels,
            text_kwargs...
        )
    end

    # Return the figure so it can be displayed or saved
    return fig
end