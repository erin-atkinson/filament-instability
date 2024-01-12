using Oceananigans
using JLD2
using CairoMakie
using Statistics
using ImageFiltering: imfilter, Kernel.gaussian


@inline function xz_state_figure(runname, n; resolution=(1000, 500))
    foldername = "../scratch/filament-instability/$runname"
    filename = "down_front_mean.jld2"
    paramfilename = "parameters.jld2"
    frames, grid = jldopen("$foldername/$filename") do file
        keys(file["timeseries/t"]), file["serialized/grid"]
        end;
    xᶜᵃᵃ = xnodes(Center, grid)
    xᶠᵃᵃ = xnodes(Face, grid)
    zᵃᵃᶜ = znodes(Center, grid)
    zᵃᵃᶠ = znodes(Face, grid)
    function ψᶜᶜᶜ(uᶠᶜᶜ, wᶜᶜᶠ, xᶜᵃᵃ, xᶠᵃᵃ, zᵃᵃᶜ, zᵃᵃᶠ)
        # Integrate
        Δzᵃᵃᶜ = reshape(diff(zᵃᵃᶠ), 1, length(zᵃᵃᶜ))
        Δx = xᶠᵃᵃ[2] - xᶠᵃᵃ[1]
        aᶠᶜᶜ = cumsum(uᶠᶜᶜ .* Δzᵃᵃᶜ; dims=3)

        bᶜᶜᶠ = cumsum(wᶜᶜᶠ .* Δx; dims=1)
        aᶜᶜᶜ = (circshift(aᶠᶜᶜ, (-1, 0)) .+ aᶠᶜᶜ) / 2
        bᶜᶜᶜ = (bᶜᶜᶠ[:, 1:end-1] .+ bᶜᶜᶠ[:, 2:end]) ./ 2
        return -aᶜᶜᶜ .+ bᶜᶜᶜ
    end
    sp = jldopen("$foldername/$paramfilename") do file
        file["parameters/simulation"]
    end
    file = jldopen("$foldername/$filename")


    frame = frames[n]


    ts = [file["timeseries/t/$f"] for f in frames] .- 1
    v = file["timeseries/v_dfm/$frame"][:, 1, :]

    u = file["timeseries/u_dfm/$frame"][:, 1, :]
    w = file["timeseries/w_dfm/$frame"][:, 1, :]

    b = file["timeseries/b_dfm/$frame"][:, 1, :]
    σ=0
    # Get the secondary cirulation streamfunction
    ψ = imfilter(ψᶜᶜᶜ(u, w, xᶜᵃᵃ, xᶠᵃᵃ, zᵃᵃᶜ, zᵃᵃᶠ), gaussian((σ, 0), (4σ+1, 5)), "circular")
    # Vorticity
    
    ζ = imfilter((circshift(v, (-1, 0)) .- circshift(v, (1, 0))) / (xᶜᵃᵃ[3] - xᶜᵃᵃ[1]), gaussian((σ, 0), (4σ+1, 5)), "circular")
    
    title = "Ro=$(round(sp.Ro; digits=1)), Ri=$(round(sp.Ri; digits=2)), t = $(round(ts[n]; digits=2))"

    axis_kwargs = (; xlabel="x", ylabel="z", title, limits=(-5, 5, -0.12, 0))

    fig = Figure(; resolution)
    ax = Axis(fig[1, 1]; axis_kwargs...)

    ht = heatmap!(ax, xᶜᵃᵃ, zᵃᵃᶜ, v; colormap=:balance, colorrange=(-5, 5))
    contour!(ax, xᶜᵃᵃ, zᵃᵃᶜ, ψ; colormap=:BrBG_10, levels=range(-0.024, 0.024, 40), alpha=1, linewidth=1)
    bstep = 1.875
    brange = minimum(b):bstep:maximum(b)
    contour!(ax, xᶜᵃᵃ, zᵃᵃᶜ, b; color=(:black, 0.5), levels=brange, linewidth=1)
    Colorbar(fig[1, 2], ht, label=L"\langle ζ \rangle")
    close(file)
    return fig
end