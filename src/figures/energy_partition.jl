using Oceananigans
using JLD2
using CairoMakie
using Statistics
using ImageFiltering: imfilter, Kernel.gaussian


@inline function energy_partition_data(foldername)
    #foldername = "../scratch/filament-instability/$runname"
    filename = "down_front_mean.jld2"
    paramfilename = "parameters.jld2"
    sp = jldopen("$foldername/$paramfilename") do file
        file["parameters/simulation"]
    end
    file = jldopen("$foldername/$filename")
    
    frames = keys(file["timeseries/t"])[1:end-1]
    grid = file["serialized/grid"]
    
    xᶜᵃᵃ = xnodes(Center, grid)
    xᶠᵃᵃ = xnodes(Face, grid)
    zᵃᵃᶜ = znodes(Center, grid)
    zᵃᵃᶠ = znodes(Face, grid)
    Δzᵃᵃᶜ = reshape(diff(zᵃᵃᶠ), 1, length(zᵃᵃᶜ))
    Δx = xᶠᵃᵃ[2] - xᶠᵃᵃ[1]
    
    z_omit_fraction=0.01
    # Boundary layer cells
    blc = zᵃᵃᶜ .> -sp.H
    # Central boundary layer cells
    cblc = -sp.H*z_omit_fraction .> zᵃᵃᶜ .> -sp.H * (1-z_omit_fraction)
    axtitle = "Ro=$(round(sp.Ro; digits=1)), Ri=$(round(sp.Ri; digits=2))"
    @inline Ev(frame) = 0.5 * file["timeseries/v_dfm/$frame"][:, 1, :] .^ 2
    @inline Eb(frame) = let b = file["timeseries/b_dfm/$frame"][:, 1, :]; cumsum(-[b[i, j] * Δzᵃᵃᶜ[j] * zᵃᵃᶜ[j] for i in 1:length(xᶜᵃᵃ), j in 1:length(zᵃᵃᶜ)]; dims=2) end
    @inline function Euw(frame)
        u = file["timeseries/u_dfm/$frame"][:, 1, :]
        w = file["timeseries/w_dfm/$frame"][:, 1, :]
        u = (circshift(u, (1, 0)) + circshift(u, (-1, 0))) / 2
        w = (w[:, 2:end] + w[:, 1:end-1]) / 2
        (u.^2 + w.^2) / 2
    end
    
    ts = [file["timeseries/t/$f"] for f in frames] .- 1
    
    data = map(frames) do frame
        fields = [Ev(frame) .* Δzᵃᵃᶜ .* Δx, Eb(frame) .* Δzᵃᵃᶜ .* Δx, Euw(frame) .* Δzᵃᵃᶜ .* Δx]
        map(x->sum(x[:, blc]), fields)
    end
    
    close(file)
    
    Ev = map(x->x[1], data)
    Eb = map(x->x[2], data)
    Eb = Eb .- Eb[1]
    Euw = map(x->x[3], data)
    Etotal = Ev .+ Eb .+ Ev
    #Ev_surface = map(x->x[4], data)
    #Eb_surface = map(x->x[5], data)
    #Euw_surface = map(x->x[6], data)
    
    return (; ts, Ev, Eb, Euw, Etotal, axtitle)
end

@inline function energy_partition!(layout_cell; ts, Ev, Eb, Euw, Etotal, axtitle, kwargs...)
    axis_kwargs = (;
        xlabel="t",
        ylabel="Energy",
        title=axtitle,
        limits=(0, ts[end], -0.05, 0.1)
    )
    
    cols = [:blue, :red, :green]
    
    ax = Axis(layout_cell; axis_kwargs...)
    lns = [lines!(ax, ts, x; color=c) for (x, c) in zip((Ev, Eb, Euw), cols)]
    [lns..., lines!(ax, ts, Etotal; color=:black)]
end
@inline function energy_partition(runnames; resolution=(1000, 250))
    n_plots = length(runnames)
    plot_datas = energy_partition_data.(runnames)
    fig = Figure(; resolution)
    lnss = map(enumerate(plot_datas)) do (i, plot_data)
        energy_partition!(fig[1, i]; plot_data...)
    end
    Legend(fig[1, n_plots+1], lnss[1], [L"\frac{1}{2}v^2", L"\int bz \text{d}z", L"\frac{1}{2}(u^2 + w^2)", "Total"])
    if n_plots > 1
        hideydecorations!.(fig.content[2:n_plots])
        for ax in fig.content[2:n_plots]
            ax.ygridvisible = true
        end
    end
    colgap!(fig.layout, 15)
    fig
end