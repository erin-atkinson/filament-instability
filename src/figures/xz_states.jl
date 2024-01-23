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

@inline function xz_state_data(runname, n; σ=0, field="v", Δ=false)
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
    v = file["timeseries/v_dfm/$frame"][:, 1, :] .- file["timeseries/v_dfm/$(frames[101])"][:, 1, :]
    w = file["timeseries/w_dfm/$frame"][:, 1, :]
    b = file["timeseries/b_dfm/$frame"][:, 1, :]
    # Get the secondary cirulation streamfunction
    ψ = ψᶜᶜᶜ(u, w, Δx, Δzᵃᵃᶜ)
    
    plt = if field in ["u", "v", "b"]
        file["timeseries/$(field)_dfm/$frame"][:, 1, :] .- Δ*file["timeseries/$(field)_dfm/$(frames[101])"][:, 1, :]
    elseif field == "w"
        let a = file["timeseries/$(field)_dfm/$frame"][:, 1, :] .- Δ*file["timeseries/$(field)_dfm/$(frames[101])"][:, 1, :]
            (a[:, 1:end-1] .+ a[:, 2:end]) / 2
        end
    elseif field == "w′v′"
        vfile = jldopen("$foldername/down_front.jld2")
        w′v′ = let vwFLUX = (vfile["timeseries/vwFLUX/$frame"][:, 1, 1:length(zᵃᵃᶜ)] .+ vfile["timeseries/vwFLUX/$frame"][:, 1, 2:length(zᵃᵃᶜ)+1]) / 2
            cumsum(Δzᵃᵃᶜ.* vwFLUX; dims=2)
        end
        close(vfile)
        w′v′
    elseif field == "ζ"
        let a = file["timeseries/v_dfm/$frame"][:, 1, :] .- Δ*file["timeseries/v_dfm/$(frames[101])"][:, 1, :]
            (circshift(a, (-1, 0)) .- circshift(a, (1, 0))) / (2Δx)
        end
    elseif field == "ω"
        let u = file["timeseries/u_dfm/$frame"][:, 1, :] .- Δ*file["timeseries/u_dfm/$(frames[101])"][:, 1, :]
            a = file["timeseries/w_dfm/$frame"][:, 1, :] .- Δ*file["timeseries/w_dfm/$(frames[101])"][:, 1, :]
            w = (a[:, 1:end-1] .+ a[:, 2:end]) / 2
            du = (u .- circshift(u, (0, 1))) / (Δzᵃᵃᶜ)
            du[:, 1] .= 0
            -(circshift(w, (-1, 0)) .- circshift(w, (1, 0))) / (2Δx) .+ du
        end
    end
    close(file)
    
    
    figtitle = "Ro=$(round(sp.Ro; digits=1)), Ri=$(round(sp.Ri; digits=2))"
    axtitle = "t = $(round(t; digits=2))"
    (plt, b, ψ) = map([plt, b, ψ]) do field
        imfilter(field, gaussian((σ, 0)), "circular")
    end
    ψ = field in ["ω", "w′v′"] ? nothing : ψ
    return (; figtitle, axtitle, xs=xᶜᵃᵃ, zs=zᵃᵃᶜ, v=plt, b, ψ)
end

@inline function xz_state!(layout_cell, cmax, limits=(-3, 3, -0.12, 0); axtitle, xs, zs, v, b, ψ, kwargs...)
    # Create a single plot in layout_cell
    axis_kwargs = (;
        xlabel="x",
        ylabel="z",
        title=axtitle,
        limits)
    
    # The steps between contours
    ψstep = 6e-4
    bstep = 1.875
    
    ψrange = ψ == nothing ? nothing : minimum(ψ):ψstep:maximum(ψ)
    brange = minimum(b):bstep:maximum(b)
    
    ax = Axis(layout_cell; axis_kwargs...)
    ht = heatmap!(ax, xs, zs, v; colormap=:balance, colorrange=(-cmax, cmax))
    ψ != nothing && contour!(ax, xs, zs, ψ; colormap=:BrBG_10, levels=ψrange, alpha=1, linewidth=1)
    contour!(ax, xs, zs, b; color=(:black, 0.5), levels=brange, linewidth=1)
    return ht
end

@inline function xz_states(runname, ns; resolution=(1000, 500), σ=0, field="v", Δ=false, cmax=nothing, limits=(-3, 3, -0.12, 0))
    n_plots = length(ns)
    plot_datas = xz_state_data.(runname, ns; σ, field, Δ)
    cmax = cmax==nothing ? maximum([maximum(abs.(a.v)) for a in plot_datas]) : cmax
    fig = Figure(; resolution)
    # Make each plot
    hts = map(enumerate(plot_datas)) do (i, plot_data)
        xz_state!(fig[1, i], cmax, limits; plot_data...)
    end
    # Add a colourbar
    Colorbar(fig[1, n_plots+1], hts[1], label=Δ ? L"\overline{%$field} - \overline{%$field}_0" : L"\overline{%$field}")
    # Add a figure title
    Label(fig[0, :], plot_datas[1].figtitle)
    # Clean up the figure
    n_plots > 1 && hideydecorations!.(fig.content[2:n_plots])
    colgap!(fig.layout, 15)
    rowgap!(fig.layout, 5)
    return fig
end