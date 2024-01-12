using Oceananigans
using JLD2
using CairoMakie
using Statistics
using ImageFiltering: imfilter, Kernel.gaussian

@inline function ψᶜᶜᶜ(uᶠᶜᶜ, wᶜᶜᶠ, Δx, Δzᵃᵃᶜ)
    aᶠᶜᶜ = cumsum(uᶠᶜᶜ .* Δzᵃᵃᶜ; dims=3)
    bᶜᶜᶠ = cumsum(wᶜᶜᶠ .* Δx; dims=1)
    aᶜᶜᶜ = (circshift(aᶠᶜᶜ, (-1, 0)) .+ aᶠᶜᶜ) / 2
    bᶜᶜᶜ = (bᶜᶜᶠ[:, 1:end-1] .+ bᶜᶜᶠ[:, 2:end]) ./ 2
    return -aᶜᶜᶜ .+ bᶜᶜᶜ
end

@inline function xz_states_data(runname, n; σ=0)
    # Read data from files to be plotted
    foldername = "../scratch/filament-instability/$runname"
    filename = "down_front_mean.jld2"
    paramfilename = "parameters.jld2"
    
    frames, grid = jldopen("$foldername/$filename") do file
        keys(file["timeseries/t"]), file["serialized/grid"]
    end
    
    xᶜᵃᵃ = xnodes(Center, grid)
    xᶠᵃᵃ = xnodes(Face, grid)
    zᵃᵃᶜ = znodes(Center, grid)
    zᵃᵃᶠ = znodes(Face, grid)
    
    Δzᵃᵃᶜ = reshape(diff(zᵃᵃᶠ), 1, length(zᵃᵃᶜ))
    Δx = xᶠᵃᵃ[2] - xᶠᵃᵃ[1]
    
    sp = jldopen("$foldername/$paramfilename") do file
        file["parameters/simulation"]
    end
    
    file = jldopen("$foldername/$filename")
    frame = frames[n]
    
    t = file["timeseries/t/$frame"] .- 1
    u = file["timeseries/u_dfm/$frame"][:, 1, :]
    v = file["timeseries/v_dfm/$frame"][:, 1, :]
    w = file["timeseries/w_dfm/$frame"][:, 1, :]
    b = file["timeseries/b_dfm/$frame"][:, 1, :]
    
    # Get the secondary cirulation streamfunction
    ψ = ψᶜᶜᶜ(u, w, Δx, Δzᵃᵃᶜ)
    
    figtitle = "Ro=$(round(sp.Ro; digits=1)), Ri=$(round(sp.Ri; digits=2))"
    axtitle = "t = $(round(t; digits=2))"
    (v, b, ψ) = map([v, b, ψ]) do field
        imfilter(field, gaussian((σ, 0)), "circular")
    end
    close(file)
    return (; figtitle, axtitle, xs=xᶜᵃᵃ, zs=zᵃᵃᶜ, v, b, ψ)
end

@inline function xz_state!(layout_cell; axtitle, xs, zs, v, b, ψ, kwargs...)
    # Create a single plot in layout_cell
    axis_kwargs = (;
        xlabel="x",
        ylabel="z",
        title=axtitle,
        limits=(-5, 5, -0.12, 0))
    
    # The steps between contours
    ψstep = 6e-4
    bstep = 1.875
    
    ψrange = minimum(ψ):ψstep:maximum(ψ)
    brange = minimum(b):bstep:maximum(b)
    
    ax = Axis(layout_cell; axis_kwargs...)
    ht = heatmap!(ax, xs, zs, v; colormap=:balance, colorrange=(-5, 5))
    contour!(ax, xs, zs, ψ; colormap=:BrBG_10, levels=ψrange, alpha=1, linewidth=1)
    contour!(ax, xs, zs, b; color=(:black, 0.5), levels=brange, linewidth=1)
    return ht
end

@inline function xz_states(runname, ns; resolution=(1000, 500))
    n_plots = length(ns)
    plot_datas = xz_states_data.(runname, ns; σ=0)
    
    fig = Figure(; resolution)
    # Make each plot
    hts = map(enumerate(plot_datas)) do (i, plot_data)
        xz_state!(fig[1, i]; plot_data...)
    end
    # Add a colourbar
    Colorbar(fig[1, n_plots+1], hts[1], label=L"\overline{v}")
    # Add a figure title
    Label(fig[0, :], plot_datas[1].figtitle)
    # Clean up the figure
    n_plots > 1 && hideydecorations!.(fig.content[2:n_plots])
    colgap!(fig.layout, 15)
    rowgap!(fig.layout, 5)
    return fig
end