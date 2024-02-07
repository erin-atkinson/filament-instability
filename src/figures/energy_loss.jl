using Oceananigans
using JLD2
using CairoMakie
using Statistics
using ImageFiltering: imfilter, Kernel.gaussian


@inline function energy_loss_data(foldername)
    bfilename = "buoyancy.jld2"
    vfilename = "down_front.jld2"
    filename = "down_front_mean.jld2"
    paramfilename = "parameters.jld2"
    sp = jldopen("$foldername/$paramfilename") do file
        file["parameters/simulation"]
    end
    file = jldopen("$foldername/$filename")
    bfile = jldopen("$foldername/$bfilename")
    vfile = jldopen("$foldername/$vfilename")
    
    frames = keys(file["timeseries/t"])[1:end-1]
    grid = file["serialized/grid"]
    
    xᶜᵃᵃ = xnodes(Center, grid)
    xᶠᵃᵃ = xnodes(Face, grid)
    zᵃᵃᶜ = znodes(Center, grid)
    zᵃᵃᶠ = znodes(Face, grid)
    Δzᵃᵃᶜ = reshape(diff(zᵃᵃᶠ), 1, length(zᵃᵃᶜ))
    Δx = xᶠᵃᵃ[2] - xᶠᵃᵃ[1]
    
    @inline function ∂x(fᶜᵃᵃ)
    return (circshift(fᶜᵃᵃ, (-1, 0)) - circshift(fᶜᵃᵃ, (1, 0))) / (2Δx)
    end
    @inline function ∂z(fᶜᵃᵃ)
        let a = (circshift(fᶜᵃᵃ, (0, -1)) - circshift(fᶜᵃᵃ, (0, 1))) ./ (2Δzᵃᵃᶜ)
            a[:, 1] .= 0
            a[:, end] .= 0
            a
        end
    end
    
    z_omit_fraction=0.1
    # Boundary layer cells
    blc = zᵃᵃᶜ .> -sp.H
    # Central boundary layer cells
    cblc = -sp.H*z_omit_fraction .> zᵃᵃᶜ .> -sp.H * (1-z_omit_fraction)
    axtitle = "Ro=$(round(sp.Ro; digits=1)), Ri=$(round(sp.Ri; digits=2))"
    @inline BFLUX(frame)= -cumsum(Δzᵃᵃᶜ.*(bfile["timeseries/bwFLUX/$frame"][:, 1, 1:length(zᵃᵃᶜ)] .+ bfile["timeseries/bwFLUX/$frame"][:, 1, 2:length(zᵃᵃᶜ)+1]); dims=2)./2
    @inline LSP(frame) = cumsum(Δx.*vfile["timeseries/vuFLUX/$frame"][:, 1, :]; dims=1) .* ∂x(file["timeseries/v_dfm/$frame"][:, 1, :])
    @inline VSP(frame) = let vwFLUX = (vfile["timeseries/vwFLUX/$frame"][:, 1, 1:length(zᵃᵃᶜ)] .+ vfile["timeseries/vwFLUX/$frame"][:, 1, 2:length(zᵃᵃᶜ)+1]) / 2
        cumsum(Δzᵃᵃᶜ.* vwFLUX; dims=2) .* ∂z(file["timeseries/v_dfm/$frame"][:, 1, :])
    end
    ts = [file["timeseries/t/$f"] for f in frames] .- 1
    
    data = map(frames) do frame
        fields = [BFLUX(frame) .* Δzᵃᵃᶜ .* Δx, LSP(frame) .* Δzᵃᵃᶜ .* Δx, VSP(frame) .* Δzᵃᵃᶜ .* Δx]
        sum.([[field[:, cblc] for field in fields]..., [field[:, blc] for field in fields]...])
    end
    close(file)
    close(bfile)
    close(vfile)
    
    BFLUX = map(x->x[1], data)
    LSP = map(x->x[2], data)
    VSP = map(x->x[3], data)
    BFLUX_surface = map(x->x[4], data)
    LSP_surface = map(x->x[5], data)
    VSP_surface = map(x->x[6], data)
    
    return (; ts, BFLUX, LSP, VSP, BFLUX_surface, LSP_surface, VSP_surface, axtitle)
end

@inline function energy_loss!(layout_cell; ts, BFLUX, LSP, VSP, BFLUX_surface, LSP_surface, VSP_surface, axtitle, kwargs...)
    axis_kwargs = (;
        xlabel="t",
        ylabel="Energy loss",
        title=axtitle,
        limits=(0, ts[2601], -0.0015, 0.0085),
        xlabelsize=16,
        ylabelsize=16)
    
    cols = [:blue, :red, :green]
    cols_boundary = [(:blue, 0.3), (:red, 0.3), (:green, 0.3)]
    
    ax = Axis(layout_cell; axis_kwargs...)
    lns = [lines!(ax, ts, x; color=c) for (x, c) in zip((BFLUX, LSP, VSP), cols)]
    [lines!(ax, ts, x; color=c) for (x, c) in zip((BFLUX_surface, LSP_surface, VSP_surface), cols_boundary)]
    return lns
end
@inline function energy_loss(runnames; resolution=(1000, 250))
    n_plots = length(runnames)
    plot_datas = energy_loss_data.(runnames)
    fig = Figure(; resolution, backgroundcolor = (:white, 0))
    lnss = map(enumerate(plot_datas)) do (i, plot_data)
        energy_loss!(fig[1, i]; plot_data...)
    end
    Legend(fig[1, n_plots+1], lnss[1], ["BFLUX", "LSP", "VSP"])
    if n_plots > 1
        hideydecorations!.(fig.content[2:n_plots])
        for ax in fig.content[2:n_plots]
            ax.ygridvisible = true
        end
    end
    colgap!(fig.layout, 15)
    fig
end