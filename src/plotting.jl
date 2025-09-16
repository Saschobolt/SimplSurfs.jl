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
    # fig = GLMakie.Figure()
    # ax = GLMakie.Axis3(fig[1, 1], aspect=:data, title=title)
    me = GeometryBasics.mesh(m; kwargs...)
    fig, ax, plt = GLMakie.mesh(me)

    # GLMakie.mesh!(ax, me)

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
    return fig, ax, plt
end

function plot_meshes(meshes::Vector{<:PolyhedralMesh};
    colors::AbstractVector{<:GLMakie.Colorant}=GLMakie.RGBf[],
    title::String="Polyhedral Meshes",
    show_labels::Bool=false,
    show_edges::Bool=true,
    edge_color=GLMakie.RGBf(0, 0, 0),
    edge_linewidth=1.5,
    kwargs...)

    # Create the figure and axis ONCE.
    fig = GLMakie.Figure(size=(800, 800))
    ax = GLMakie.Axis3(fig[1, 1], aspect=:data, title=title)

    # --- Color Handling ---
    # Set up default colors if none are provided
    final_colors = if isempty(colors)
        default_palette = [GLMakie.RGBf(0, 48 / 255, 144 / 255), GLMakie.RGBf(239 / 255, 123 / 255, 0), GLMakie.RGBf(243 / 255, 207 / 255, 198 / 255), GLMakie.RGBf(237 / 255, 28 / 255, 36 / 255), GLMakie.RGBf(255 / 255, 219 / 255, 88 / 255)]
        # Cycle through the default palette to get enough colors
        [default_palette[mod1(i, length(default_palette))] for i in 1:length(meshes)]
    else
        # Input validation
        if length(colors) != length(meshes)
            error("The number of colors ($(length(colors))) must match the number of meshes ($(length(meshes))).")
        end
        colors
    end

    # --- Loop and Plot ---
    # Iterate through each mesh and its corresponding color
    for (m, mesh_color) in zip(meshes, final_colors)
        # Plot the mesh onto the EXISTING axis using mesh!
        # We pass the single color for this mesh via the `color` keyword.
        me = GeometryBasics.mesh(m; face_colors=[mesh_color for f in faces(m)])
        GLMakie.mesh!(ax, me; color=mesh_color)

        if show_edges
            GLMakie.wireframe!(ax, me;
                color=edge_color,
                linewidth=edge_linewidth,
            )
        end

        if show_labels
            vertex_positions = coordinates(m)
            vertex_labels = [string(label(v)) for v in vertices(m)]

            text_defaults = (fontsize=20, align=(:left, :center), offset=(5, 0))
            # Correctly merge kwargs to allow user overrides
            text_kwargs = merge(text_defaults)

            GLMakie.text!(ax,
                vertex_positions;
                text=vertex_labels,
                # The text color can also be customized via kwargs
                color=haskey(kwargs, :text_color) ? kwargs[:text_color] : :black,
                text_kwargs...
            )
        end
    end

    return fig
end