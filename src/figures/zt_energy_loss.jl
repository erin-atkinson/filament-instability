using Oceananigans
using JLD2
using CairoMakie
using Statistics
using ImageFiltering: imfilter, Kernel.gaussian


@inline function zt_energy_loss_data(runname)
    foldername = "../scratch/filament-instability/$runname"
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
    
    frames = keys(file["timeseries/t"])[101:2601]
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
    
    z_omit_fraction=0.01
    # Boundary layer cells
    blc = zᵃᵃᶜ .> -sp.H
    # Central boundary layer cells
    cblc = -sp.H*z_omit_fraction .> zᵃᵃᶜ .> -sp.H * (1-z_omit_fraction)
    figtitle = "Ro=$(round(sp.Ro; digits=1)), Ri=$(round(sp.Ri; digits=2))"
    @inline BFLUX(frame)= -cumsum(Δzᵃᵃᶜ.*(bfile["timeseries/bwFLUX/$frame"][:, 1, 1:length(zᵃᵃᶜ)] .+ bfile["timeseries/bwFLUX/$frame"][:, 1, 2:length(zᵃᵃᶜ)+1]); dims=2)./2
    @inline LSP(frame) = cumsum(Δx.*vfile["timeseries/vuFLUX/$frame"][:, 1, :]; dims=1) .* ∂x(file["timeseries/v_dfm/$frame"][:, 1, :])
    @inline VSP(frame) = let vwFLUX = (vfile["timeseries/vwFLUX/$frame"][:, 1, 1:length(zᵃᵃᶜ)] .+ vfile["timeseries/vwFLUX/$frame"][:, 1, 2:length(zᵃᵃᶜ)+1]) / 2
        cumsum(Δzᵃᵃᶜ.* vwFLUX; dims=2) .* ∂z(file["timeseries/v_dfm/$frame"][:, 1, :])
    end
    ts = [file["timeseries/t/$f"] for f in frames] .- 1
    
    data = map(frames) do frame
        fields = [BFLUX(frame) .* Δx, LSP(frame) .* Δx, VSP(frame) .* Δx]
        sum.(fields; dims=1)
    end
    close(file)
    close(bfile)
    close(vfile)
    
    BFLUX = [data[i][1][1, j] for i in 1:length(data), j in 1:size(data[1][1], 2)]
    LSP = [data[i][2][1, j] for i in 1:length(data), j in 1:size(data[1][1], 2)]
    VSP = [data[i][3][1, j] for i in 1:length(data), j in 1:size(data[1][1], 2)]
    
    return (; ts, zs=zᵃᵃᶜ, BFLUX, LSP, VSP, figtitle)
end

@inline function zt_energy_loss!(fig; ts, zs, BFLUX, LSP, VSP, kwargs...)
    
    axis_kwargs = (;
        xlabel="t",
        ylabel="z",
        limits=(ts[1], ts[end], zs[33], zs[end]),
    )
    # Color limit should be the largest value at the start, ignore the end bit
    cmax = max(
        maximum(abs.(BFLUX[1:2000, :])),
        maximum(abs.(LSP[1:2000, :])),
        maximum(abs.(VSP[1:2000, :]))
    )
    axtitles = ["BFLUX", "LSP", "VSP"]
    hts = map(enumerate([BFLUX, LSP, VSP])) do (i, field)
        ax = Axis(fig[1, i]; axis_kwargs..., title=axtitles[i])
        heatmap!(ax, ts, zs, field; colormap=:balance, colorrange=(-cmax, cmax))
    end
    return hts
end
@inline function zt_energy_loss(runname; resolution=(1000, 300))
    plot_data = zt_energy_loss_data(runname)
    fig = Figure(; resolution)
    hts = zt_energy_loss!(fig; plot_data...)
    Colorbar(fig[1, 4], hts[1])
    Label(fig[0, :], plot_data.figtitle)
    hideydecorations!.(fig.content[2:3])
    for ax in fig.content[2:3]
        ax.ygridvisible = true
    end
    colgap!(fig.layout, 15)
    rowgap!(fig.layout, 5)
    fig
end