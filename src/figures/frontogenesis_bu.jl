using Oceananigans
using JLD2
using CairoMakie
using Statistics
using ImageFiltering: imfilter, Kernel.gaussian

@inline function average_∇b²(runname, slice; σ=0)
    foldername = "../scratch/filament-instability/$runname"
    filename = "down_front_mean.jld2"
    
    paramfilename = "parameters.jld2"
    sp = jldopen("$foldername/$paramfilename") do file
        file["parameters/simulation"]
    end
    
    file = jldopen("$foldername/$filename")
    
    frames = keys(file["timeseries/t"])[101:2601]
    grid = file["serialized/grid"]
    xᶠᵃᵃ = xnodes(Face, grid)
    Δx = xᶠᵃᵃ[2] - xᶠᵃᵃ[1]
    
    ts = [file["timeseries/t/$frame"] for frame in frames] .- 1
    ∇b² = map(frames) do frame
        b = imfilter(mean(file["timeseries/b_dfm/$frame"][:, 1, slice]; dims=2)[:, 1], gaussian((σ, )))
        mean(((circshift(b, -1) - circshift(b, 1)) ./ Δx).^2)# / mean(b.^2)
    end
    close(file)
    return (; ts, ∇b²)
end

@inline function symmetric_u(runname, slice)
    foldername = "../scratch/filament-instability/$runname"
    filename = "down_front_mean.jld2"
    
    paramfilename = "parameters.jld2"
    sp = jldopen("$foldername/$paramfilename") do file
        file["parameters/simulation"]
    end
    
    file = jldopen("$foldername/$filename")
    
    frames = keys(file["timeseries/t"])[101:2601]
    
    ts = [file["timeseries/t/$frame"] for frame in frames] .- 1
    u_sym = map(frames) do frame
        mean(file["timeseries/u_dfm/$frame"][end÷2:end, 1, slice]) - mean(file["timeseries/u_dfm/$frame"][1:end÷2, 1, slice])
    end
    close(file)
    return (; ts, u_sym)
end

@inline function plot_∇b²!(layout_cell, runnames, slice; σ=0)
    plot_datas = map(rn->average_∇b²(rn, slice; σ=0), runnames)
    
    axis_kwargs = (;
        xlabel="t",
        ylabel=L"|\nabla b|^2"
    )
    
    ax = Axis(layout_cell; axis_kwargs...)
    αs = reverse(range(0.3, 1.0, length(runnames)))
    map(zip(plot_datas, αs)) do (plot_data, α)
        lines!(ax, plot_data.ts, plot_data.∇b²; color=(:red, α))
    end
end

@inline function plot_u_sym!(layout_cell, runnames, slice)
    plot_datas = map(rn->symmetric_u(rn, slice), runnames)
    
    axis_kwargs = (;
        xlabel="t",
        ylabel=L"u_\text{asym}"
    )
    
    ax = Axis(layout_cell; axis_kwargs...)
    αs = reverse(range(0.3, 1.0, length(runnames)))
    map(zip(plot_datas, αs)) do (plot_data, α)
        lines!(ax, plot_data.ts, plot_data.u_sym; color=(:red, α))
    end
end

@inline function frontogenesis_bu(runnames, slice, legendtitle, runlabels; σ=0, resolution=(1000, 500))
    fig = Figure(; resolution)
    plot_∇b²!(fig[1, 1], runnames, slice; σ=0)
    lns = plot_u_sym!(fig[1, 2], runnames, slice)
    Legend(fig[1, 3], lns, runlabels, legendtitle)
    colgap!(fig.layout, 15)
    fig
end