using Oceananigans
using JLD2
using CairoMakie
using Statistics
using ImageFiltering: imfilter, Kernel.gaussian

@inline function yz_state_data(runname, n, slice; σ=0, field="v")
    foldername = "../scratch/filament-instability/$runname"
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
    x_text = if length(slice) > 1
        "[$(round(xᶜᵃᵃ[slice[1]]; digits=2)), $(round(xᶜᵃᵃ[slice[end]]; digits=2))]"
    else
        "$(round(xᶜᵃᵃ[slice[1]]; digits=2))"
    end
    
    
    # We just need the v and b
    v = mean(file["$field/data"][slice, 4:end-3, 4:end-3]; dims=1)[1, :, end-length(zᵃᵃᶜ)+1:end]
    b = mean(file["b/data"][slice, 4:end-3, 4:end-3]; dims=1)[1, :, end-length(zᵃᵃᶜ)+1:end]
    t = file["clock"].time
    close(file)
    
    figtitle = "Ro=$(round(sp.Ro; digits=1)), Ri=$(round(sp.Ri; digits=2)), x = $x_text"
    axtitle = "t = $(round(t; digits=2))"
    
    (v, b) = map([v, b]) do field
        imfilter(field, gaussian((σ, 0)), "circular")
    end
    
    return (; figtitle, axtitle, ys=yᵃᶜᵃ, zs=zᵃᵃᶜ, v, b)
end

@inline function yz_state!(layout_cell; axtitle, ys, zs, v, b, max_val=5, kwargs...)
    # Create a single plot in layout_cell
    axis_kwargs = (;
        xlabel="y",
        ylabel="z",
        title=axtitle,
        limits=(-5, 5, -0.1, 0))
    
    # The steps between contours
    bstep = 2*1.875
    
    brange = minimum(b):bstep:maximum(b)
    
    ax = Axis(layout_cell; axis_kwargs...)
    ht = heatmap!(ax, ys, zs, v; colormap=:balance, colorrange=(-max_val, max_val))
    contour!(ax, ys, zs, b; color=(:black, 0.3), levels=brange, linewidth=1)
    return ht
end

@inline function yz_states(runname, ns, slice; resolution=(1000, 500), σ=16, field="v")
    n_plots = length(ns)
    plot_datas = [yz_state_data(runname, n, slice; σ, field) for n in ns]
    max_val = maximum(map(plot_data->maximum(abs.(plot_data.v)), plot_datas))
    fig = Figure(; resolution)
    # Make each plot
    hts = map(enumerate(plot_datas)) do (i, plot_data)
        yz_state!(fig[1, i]; plot_data..., max_val)
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