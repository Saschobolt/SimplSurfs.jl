using GLMakie
using LinearAlgebra
using GeometryBasics
using Statistics

###################################################################################
# Animation of one Waterbomb Cell
###################################################################################
# It is assumed that all the previously defined structs and functions are
# available in the current scope.

# Helper function to convert the 3xN coordinate matrix to a Vector of Point3f
function points_from_matrix(mat::AbstractMatrix)
    return [Point3f(mat[:, i]) for i in 1:size(mat, 2)]
end

# --- Animation Setup ---
duration = 8 # seconds
framerate = 60
n_frames = duration * framerate
filename = "waterbomb_kinematics.mp4"
r_range = (0.95, 1.05)
s_range = (sqrt(2) - 0.1, sqrt(2) + 0.1)
t_range = (1.96, 2.0)


# --- Initial State ---
initial_cell = WaterbombCell()


# --- GLMakie Scene Setup ---
set_theme!(theme_black())
fig = Figure(size=(800, 800))
ax = Axis3(fig[1, 1],
    aspect=:data,
    title="Waterbomb Cell Kinematics",
    protrusions=(0, 0, 0, 0))

# --- Observables for Dynamic Plotting ---
observable_coords = Observable(points_from_matrix(coordinate_matrix(initial_cell)))

fs = [GeometryBasics.NgonFace{length(vertices(f)),Int}(id.(vertices(f))...) for f in SimplSurfs.faces(initial_cell.mesh)]
observable_mesh = @lift(GeometryBasics.mesh($observable_coords, fs))

mesh!(ax, observable_mesh, color=GLMakie.RGBAf(0.5, 0.75, 1.0, 1.0), shading=FastShading)
wireframe!(ax, observable_mesh, color=:white, linewidth=1.5)

# --- Add Observables and Plots for Dynamic Lines and Labels using @lift ---

# A. Line for 't' distance between vertices (3, 7) - Red
line_37_points = @lift([$observable_coords[3], $observable_coords[7]])
dist_37_text = @lift("t = " * string(round(norm($observable_coords[3] - $observable_coords[7]), digits=2)))
pos_37_text = @lift([($observable_coords[3] + $observable_coords[7]) / 2])

lines!(ax, line_37_points, color=:red, linestyle=:dot, linewidth=4)
text!(ax, pos_37_text, text=dist_37_text, color=:red, fontsize=20, align=(:center, :bottom), space=:data)

# B. Line for 'r' distance between vertices (2, 5) - Green
line_25_points = @lift([$observable_coords[2], $observable_coords[5]])
dist_25_text = @lift("r = " * string(round(norm($observable_coords[2] - $observable_coords[5]), digits=2)))
pos_25_text = @lift([($observable_coords[2] + $observable_coords[5]) / 2])

lines!(ax, line_25_points, color=:lime, linestyle=:dot, linewidth=4)
text!(ax, pos_25_text, text=dist_25_text, color=:lime, fontsize=20, align=(:center, :bottom), space=:data)

# C. Line for 's' distance between vertices (3, 5) - Yellow
line_35_points = @lift([$observable_coords[3], $observable_coords[5]])
dist_35_text = @lift("s = " * string(round(norm($observable_coords[3] - $observable_coords[5]), digits=2)))
pos_35_text = @lift([($observable_coords[3] + $observable_coords[5]) / 2])

lines!(ax, line_35_points, color=:yellow, linestyle=:dot, linewidth=4)
text!(ax, pos_35_text, text=dist_35_text, color=:yellow, fontsize=20, align=(:center, :bottom), space=:data)


# --- Animation Loop ---
record(fig, filename, 1:n_frames; framerate=framerate) do frame
    angle = 2π * (frame - 1) / n_frames

    current_r = mean(r_range) + (r_range[2] - r_range[1]) / 2 * cos(angle)
    current_s = mean(s_range) + (s_range[2] - s_range[1]) / 2 * cos(angle)
    current_t = mean(t_range) + (t_range[2] - t_range[1]) / 2 * cos(angle)

    update_coordinates!(initial_cell; r=current_r, s=current_s, t=current_t)

    observable_coords[] = points_from_matrix(coordinate_matrix(initial_cell))

    # *** CHANGE IS HERE: Use the new azimuth function ***
    ax.azimuth[] = 1.75pi + 0.175 * sin(2pi * frame / n_frames)

    # Re-enable autolimits on every frame to adjust the zoom dynamically
    if frame == 1
        autolimits!(ax)
    end
end

println("Animation saved to $(filename)")


####################################################################################
# Animation of Three Connected Waterbomb Cells
####################################################################################

# Helper function to convert the 3xN coordinate matrix to a Vector of Point3f
function points_from_matrix(mat::AbstractMatrix)
    return [Point3f(mat[:, i]) for i in 1:size(mat, 2)]
end

# --- Animation Setup ---
duration = 9 # seconds
framerate = 60
n_frames = duration * framerate
filename = "three_cell_kinematics.mp4"

# Define the start and finish parameters for the animation of cell2
# *** FIX: Renamed field `end` to `finish` ***
r_params = (start=1.0, finish=0.97)
s_params = (start=sqrt(2), finish=1.39)
t_params = (start=2.0, finish=1.96)

# Define distinct colors for each cell
cell_colors = [
    GLMakie.RGBAf(0.5, 0.75, 1.0, 1.0), # Light Blue for Cell 1
    GLMakie.RGBAf(0.5, 1.0, 0.75, 1.0), # Light Green for Cell 2
    GLMakie.RGBAf(1.0, 0.75, 0.5, 1.0)  # Light Orange for Cell 3
]

# --- Initial State Construction ---
cell1 = WaterbombCell()
cell2 = attach_cell!(cell1)
cell3 = attach_cell!(cell1, cell2)


# --- GLMakie Scene Setup ---
set_theme!(theme_black())
fig = Figure(size=(1000, 800))
ax = Axis3(fig[1, 1],
    aspect=:data,
    title="Connected Waterbomb Cell Kinematics",
    protrusions=(0, 0, 0, 0))

# --- Helper function to plot one cell's details ---
function plot_cell_details!(ax, coords_obs, mesh_color)
    # Derived observable for the mesh geometry
    fs = [GeometryBasics.NgonFace{length(vertices(f)),Int}(id.(vertices(f))...) for f in SimplSurfs.faces(cell1.mesh)]
    mesh_obs = @lift(GeometryBasics.mesh($coords_obs, fs))

    # Plot the mesh and wireframe
    mesh!(ax, mesh_obs, color=mesh_color, shading=true)
    wireframe!(ax, mesh_obs, color=:white, linewidth=1.0)

    # Plot the r, s, t lines and labels
    # r-line (2, 5)
    r_pts = @lift([$coords_obs[2], $coords_obs[5]])
    r_txt = @lift("r = " * string(round(norm($coords_obs[2] - $coords_obs[5]), digits=2)))
    r_pos = @lift([($coords_obs[2] + $coords_obs[5]) / 2])
    lines!(ax, r_pts, color=:lime, linestyle=:dot, linewidth=3)
    text!(ax, r_pos, text=r_txt, color=:lime, fontsize=16, align=(:center, :bottom))

    # s-line (3, 5)
    s_pts = @lift([$coords_obs[3], $coords_obs[5]])
    s_txt = @lift("s = " * string(round(norm($coords_obs[3] - $coords_obs[5]), digits=2)))
    s_pos = @lift([($coords_obs[3] + $coords_obs[5]) / 2])
    lines!(ax, s_pts, color=:yellow, linestyle=:dot, linewidth=3)
    text!(ax, s_pos, text=s_txt, color=:yellow, fontsize=16, align=(:center, :bottom))

    # t-line (3, 7)
    t_pts = @lift([$coords_obs[3], $coords_obs[7]])
    t_txt = @lift("t = " * string(round(norm($coords_obs[3] - $coords_obs[7]), digits=2)))
    t_pos = @lift([($coords_obs[3] + $coords_obs[7]) / 2])
    lines!(ax, t_pts, color=:red, linestyle=:dot, linewidth=3)
    text!(ax, t_pos, text=t_txt, color=:red, fontsize=16, align=(:center, :bottom))
end

# --- Create Observables and Plot Each Cell ---
coords1_obs = Observable(points_from_matrix(coordinate_matrix(cell1)))
coords2_obs = Observable(points_from_matrix(coordinate_matrix(cell2)))
coords3_obs = Observable(points_from_matrix(coordinate_matrix(cell3)))

plot_cell_details!(ax, coords1_obs, cell_colors[1])
plot_cell_details!(ax, coords2_obs, cell_colors[2])
plot_cell_details!(ax, coords3_obs, cell_colors[3])


# --- Animation Loop ---
record(fig, filename, 1:n_frames; framerate=framerate) do frame
    # Linear interpolation from start to finish over the duration
    progress = (frame - 1) / (n_frames - 1)

    # An angle from 0 to 2π ensures the animation starts and ends at the same state.
    angle = 2π * (frame - 1) / n_frames

    # This formula oscillates smoothly between the start and finish values.
    # It finds the midpoint and adds an oscillating value scaled by the amplitude.
    current_r = mean([r_params.start, r_params.finish]) + (r_params.start - r_params.finish) / 2 * cos(angle)
    current_s = mean([s_params.start, s_params.finish]) + (s_params.start - s_params.finish) / 2 * cos(angle)
    current_t = mean([t_params.start, t_params.finish]) + (t_params.start - t_params.finish) / 2 * cos(angle)

    # Update cell2. Because cell3 depends on cell2, the `update_coordinates!`
    # function will automatically update cell3's coordinates as well.
    update_coordinates!(cell2; r=current_r, s=current_s, t=current_t)

    # Notify the observables for the cells that have changed
    coords2_obs[] = points_from_matrix(coordinate_matrix(cell2))
    coords3_obs[] = points_from_matrix(coordinate_matrix(cell3))

    # Update camera
    ax.azimuth[] = 1.8pi + 0.1 * sin(2pi * frame / n_frames)

    # Adjust zoom
    autolimits!(ax)
end

println("Animation saved to $(filename)")