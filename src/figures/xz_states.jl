using Oceananigans
using JLD2
using CairoMakie
using Statistics
using ImageFiltering: imfilter, Kernel.gaussian

include("../filament_state.jl")

@inline function ψᶜᶜᶜ(uᶠᶜᶜ, wᶜᶜᶠ, Δx, Δzᵃᵃᶜ)
    aᶠᶜᶜ = cumsum(uᶠᶜᶜ .* Δzᵃᵃᶜ; dims=3)
    bᶜᶜᶠ = cumsum(wᶜᶜᶠ .* Δx; dims=1)
    aᶜᶜᶜ = (circshift(aᶠᶜᶜ, (-1, 0)) .+ aᶠᶜᶜ) / 2
    bᶜᶜᶜ = (bᶜᶜᶠ[:, 1:end-1] .+ bᶜᶜᶠ[:, 2:end]) ./ 2
    return -aᶜᶜᶜ .+ bᶜᶜᶜ
end

@inline function ref_q(sp)
    # Return a function that describes the reference state for potential vorticity
    (b, v) = get_filament_state(sp; verbose=false)
    @inline ∂xb(x, z) = (b(x + 1e-5, z) - b(x - 1e-5, z)) / 2e-5
    @inline ∂zb(x, z) = (b(x, z + 1e-6) - b(x, z - 1e-6)) / 2e-6
    @inline ∂xv(x, z) = (v(x + 1e-5, z) - v(x - 1e-5, z)) / 2e-5
    @inline ∂zv(x, z) = (v(x, z + 1e-6) - v(x, z - 1e-6)) / 2e-6
    
    @inline q(x, z) = (∂xv(x, z) + sp.f) * ∂zb(x, z) .- ∂zv(x, z) * ∂xb(x, z)
end

@inline function ref_Ri(sp)
    # Return a function that describes the reference state for potential vorticity
    (b, v) = get_filament_state(sp; verbose=false)
    @inline ∂zb(x, z) = (b(x, z + 1e-6) - b(x, z - 1e-6)) / 2e-6
    @inline ∂zv(x, z) = (v(x, z + 1e-6) - v(x, z - 1e-6)) / 2e-6
    
    @inline Ri(x, z) = ∂zb(x, z) / ∂zv(x, z)^2
end

@inline function xz_state_data(foldername, n; σ=0, field="v", Δ=false)
    n₀ = 21
    
    # Read data from files to be plotted
    #foldername = "../scratch/filament-instability/$runname"
    filename = "down_front_mean.jld2"
    paramfilename = "parameters.jld2"
    
    frames, grid = jldopen("$foldername/$filename") do file
        keys(file["timeseries/t"]), file["serialized/grid"]
    end
    
    xᶜᵃᵃ = xnodes(grid, Center())
    xᶠᵃᵃ = xnodes(grid, Face())
    zᵃᵃᶜ = znodes(grid, Center())
    zᵃᵃᶠ = znodes(grid, Face())
    
    Δzᵃᵃᶜ = reshape(diff(zᵃᵃᶠ), 1, length(zᵃᵃᶜ))
    Δx = xᶠᵃᵃ[2] - xᶠᵃᵃ[1]
    
    sp = jldopen("$foldername/$paramfilename") do file
        file["parameters/simulation"]
    end
    
    file = jldopen("$foldername/$filename")
    frame = frames[n]
    
    t = file["timeseries/t/$frame"] .- 1
    u = file["timeseries/u_dfm/$frame"][:, 1, :]
    #v = file["timeseries/v_dfm/$frame"][:, 1, :]
    w = file["timeseries/w_dfm/$frame"][:, 1, :]
    b = file["timeseries/b_dfm/$frame"][:, 1, :]
    # Get the secondary cirulation streamfunction
    ψ = ψᶜᶜᶜ(u, w, Δx, Δzᵃᵃᶜ)
    
    plt = if field in ["u", "v", "b"]
        file["timeseries/$(field)_dfm/$frame"][:, 1, :] .- Δ*file["timeseries/$(field)_dfm/$(frames[n₀])"][:, 1, :]
    elseif field == "w"
        let a = file["timeseries/$(field)_dfm/$frame"][:, 1, :] .- Δ*file["timeseries/$(field)_dfm/$(frames[n₀])"][:, 1, :]
            (a[:, 1:end-1] .+ a[:, 2:end]) / 2
        end
    elseif field == "w′v′"
        vfile = jldopen("$foldername/down_front.jld2")
        w′v′ = let vwFLUX = (vfile["timeseries/vwFLUX/$frame"][:, 1, 1:length(zᵃᵃᶜ)] .+ vfile["timeseries/vwFLUX/$frame"][:, 1, 2:length(zᵃᵃᶜ)+1]) / 2
            cumsum(Δzᵃᵃᶜ.* vwFLUX; dims=2)
        end
        close(vfile)
        w′v′
    elseif field == "w′v′_z"
        vfile = jldopen("$foldername/down_front.jld2")
        vwFLUX = (vfile["timeseries/vwFLUX/$frame"][:, 1, 1:length(zᵃᵃᶜ)] .+ vfile["timeseries/vwFLUX/$frame"][:, 1, 2:length(zᵃᵃᶜ)+1]) / 2
        close(vfile)
        vwFLUX
    elseif field == "wv_z"
        v_z = diff(file["timeseries/v_dfm/$frame"][:, 1, :]; dims=2) ./ Δzᵃᵃᶜ[:, 2:end]
        w = file["timeseries/w_dfm/$frame"][:, 1, :]
        w[:, 2:end-1] .*= v_z
        (w[:, 1:end-1] .+ w[:, 2:end]) ./ 2
    elseif field == "ζ"
        let a = file["timeseries/v_dfm/$frame"][:, 1, :] .- Δ*file["timeseries/v_dfm/$(frames[n₀])"][:, 1, :]
            (circshift(a, (-1, 0)) .- circshift(a, (1, 0))) / (2Δx)
        end
    elseif field == "ω"
        let u = file["timeseries/u_dfm/$frame"][:, 1, :] .- Δ*file["timeseries/u_dfm/$(frames[n₀])"][:, 1, :]
            a = file["timeseries/w_dfm/$frame"][:, 1, :] .- Δ*file["timeseries/w_dfm/$(frames[n₀])"][:, 1, :]
            w = (a[:, 1:end-1] .+ a[:, 2:end]) / 2
            du = (u .- circshift(u, (0, 1))) / (Δzᵃᵃᶜ)
            du[:, 1] .= 0
            -(circshift(w, (-1, 0)) .- circshift(w, (1, 0))) / (2Δx) .+ du
        end
    end
    close(file)
    
    axtitle = L"Ri_{\text{min}}=%$(round(sp.Ri; digits=2))\quad t / 2\pi = %$(round(t / (2π); digits=2))"
    (plt, b, ψ) = map([plt, b, ψ]) do field
        imfilter(field, gaussian((σ, 0)), "circular")
    end
    ψ = field in ["ω", "w′v′", "ζ", "w"] ? nothing : ψ
    q = nothing #field in ["ω", "w′v′"] ? ref_q(sp) : nothing
    Ri = nothing #field in ["ω", "w′v′"] ? ref_Ri(sp) : nothing
    return (; axtitle, xs=xᶜᵃᵃ, zs=zᵃᵃᶜ, v=plt, b, ψ, q, Ri)
end

@inline function xz_state!(layout_cell, cmax, limits=(-3, 3, -0.12, 0); axtitle, xs, zs, v, b, ψ, q, Ri, kwargs...)
    # Create a single plot in layout_cell
    axis_kwargs = (;
        xlabel=L"x",
        ylabel=L"z",
        title=axtitle,
        limits)
    
    # The steps between contours
    ψstep = 12e-4
    bstep = 3.75
    
    ψrange = ψ == nothing ? nothing : -maximum(abs.(ψ)):ψstep:maximum(abs.(ψ))
    brange = minimum(b):bstep:maximum(b)
    
    ax = Axis(layout_cell; axis_kwargs...)
    ht = heatmap!(ax, xs, zs, v; colormap=:balance, colorrange=(-cmax, cmax))
    ψ != nothing && contour!(ax, xs, zs, ψ; colormap=:BrBG_10, levels=ψrange, alpha=1, linewidth=2)
    q != nothing && contour!(ax, xs, zs, q; color=:red, levels=[0], alpha=1, linewidth=1)
    Ri != nothing && contour!(ax, xs, zs, Ri; color=:blue, levels=[0.25, 1], alpha=1, linewidth=1)
    contour!(ax, xs, zs, b; color=(:black, 0.5), levels=brange, linewidth=2)
    #contour!(ax, xs, zs, b; color=(:red, 0.5), levels=brange[end-7:end-5], linewidth=1)
    return ht
end

@inline function xz_states(runname, ns; size=(1000, 500), σ=0, field="v", Δ=false, cmax=nothing, limits=(-3, 3, -0.12, 0))
    n_plots = length(ns)
    plot_datas = xz_state_data.(runname, ns; σ, field, Δ)
    cmax = cmax==nothing ? maximum([maximum(abs.(a.v)) for a in plot_datas]) : cmax
    fig = Figure(; size, backgroundcolor = (:white, 0.0), fontsize=16)
    # Make each plot
    hts = map(enumerate(plot_datas)) do (i, plot_data)
        xz_state!(fig[1, i], cmax, limits; plot_data...)
    end
    # Add a colourbar
    Colorbar(fig[1, n_plots+1], hts[1], label=Δ ? L"\overline{%$field} - \overline{%$field}_0" : L"\overline{%$field}")
    # Clean up the figure
    n_plots > 1 && hideydecorations!.(fig.content[2:n_plots])
    if n_plots > 1
    for i in 1:n_plots-1
            break
        colgap!(fig.layout, i, 40)
    end
    end
    #rowgap!(fig.layout, 5)
    return fig
end