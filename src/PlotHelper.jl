using Plots

""" Helper Function for 2D contour plot for linear functions f(x) = c'x"""
function plot_lin_contour!(c, grid_res; kwargs...)
    f(x) = dot(c,x)
    contour!(grid_res, grid_res, (x,y) -> f([x,y]); kwargs...)
end

""" Helper Function for Plotting an arrow in 2D"""
function plot_vec!(x, label="", x_start= [0.0,0.0]; kwargs...)
    x_end = x
    x_pos = [x_start[1], x_end[1]]
    y_pos = [x_start[2], x_end[2]]
    plot!(x_pos, y_pos, arrow=(1,:black); kwargs...)
end
