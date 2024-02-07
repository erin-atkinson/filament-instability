using Oceananigans
using JLD2
using CairoMakie
using Statistics
using ImageFiltering: imfilter, Kernel.gaussian

@inline function xy_state_data(foldername, n, slice; σ=0, field="v")
    #foldername = "../scratch/filament-instability/$runname"
    filename = "down_front_mean.jld2"
    paramfilename = "parameters.jld2"
    
    # Should sort by iteration
    filename_iters = map(readdir("../scratch/filament-instability/$runname/checkpoints")) do f
        m = match(r"checkpoint_iteration(\d+)\.jld2", f)
        parse(Int, m[1])
    end
    iteration = sort(filename_iters)[n]
    checkpoint_filename = "../scratch/filament-instability/$runname/checkpoints/checkpoint_iteration$iteration.jld2"
    
    sp = jldopen("$foldername/$paramfilename") do file
        file["parameters/simulation"]
    end
    
    file = jldopen(checkpoint_filename)
    
    # Grid
    grid = file["grid"]
    xᶜᵃᵃ = xnodes(Center, grid)
    xᶠᵃᵃ = xnodes(Face, grid)
    yᵃᶜᵃ = ynodes(Center, grid)
    zᵃᵃᶜ = znodes(Center, grid)
    zᵃᵃᶠ = znodes(Face, grid)
    z_text = if length(slice) > 1
        "[$(round(zᵃᵃᶜ[slice[1]]; digits=2)), $(round(zᵃᵃᶜ[slice[end]]; digits=2))]"
    else
        "$(round(zᵃᵃᶜ[slice[1]]; digits=2))"
    end
    
    
    # We just need the v and b
    v = mean(file["$field/data"][4:end-3, 4:end-3, slice]; dims=3)[:, :, 1]
    b = mean(file["b/data"][4:end-3, 4:end-3, slice]; dims=3)[:, :, 1]
    t = file["clock"].time
    close(file)
    
    figtitle = "Ro=$(round(sp.Ro; digits=1)), Ri=$(round(sp.Ri; digits=2)), z = $z_text"
    axtitle = "t = $(round(t; digits=2))"
    
    (v, b) = map([v, b]) do field
        imfilter(field, gaussian((σ, σ)), "circular")
    end
    
    return (; figtitle, axtitle, xs=xᶜᵃᵃ, ys=yᵃᶜᵃ, v, b)
end

@inline function xy_state!(layout_cell; axtitle, xs, ys, v, b, max_val=5, kwargs...)
    # Create a single plot in layout_cell
    axis_kwargs = (;
        xlabel="x",
        ylabel="y",
        title=axtitle,
        limits=(-5, 5, -5, 5))
    
    # The steps between contours
    bstep = 2*1.875
    
    brange = minimum(b):bstep:maximum(b)
    
    ax = Axis(layout_cell; axis_kwargs...)
    ht = heatmap!(ax, xs, ys, v; colormap=:balance, colorrange=(-max_val, max_val))
    contour!(ax, xs, ys, b; color=(:black, 0.3), levels=brange, linewidth=1)
    return ht
end

@inline function xy_states(runname, ns, slice; resolution=(1000, 500), σ=16, field="v")
    n_plots = length(ns)
    plot_datas = [xy_state_data(runname, n, slice; σ, field) for n in ns]
    max_val = maximum(map(plot_data->maximum(abs.(plot_data.v)), plot_datas))
    fig = Figure(; resolution)
    # Make each plot
    hts = map(enumerate(plot_datas)) do (i, plot_data)
        xy_state!(fig[1, i]; plot_data..., max_val)
    end
    # Add a colourbar
    Colorbar(fig[1, n_plots+1], hts[1], label=L"%$field (x, y, z)")
    # Add a figure title
    Label(fig[0, :], plot_datas[1].figtitle)
    # Clean up the figure
    n_plots > 1 && hideydecorations!.(fig.content[2:n_plots])
    colgap!(fig.layout, 15)
    rowgap!(fig.layout, 5)
    return fig
end