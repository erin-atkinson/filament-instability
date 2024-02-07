using Oceananigans
using JLD2
using CairoMakie
using Statistics
using ImageFiltering: imfilter, Kernel.gaussian
using ZipFile

@inline function ψᶜᶜᶜ(uᶠᶜᶜ, wᶜᶜᶠ, Δx, Δzᵃᵃᶜ)
    aᶠᶜᶜ = cumsum(uᶠᶜᶜ .* Δzᵃᵃᶜ; dims=3)
    bᶜᶜᶠ = cumsum(wᶜᶜᶠ .* Δx; dims=1)
    aᶜᶜᶜ = (circshift(aᶠᶜᶜ, (-1, 0)) .+ aᶠᶜᶜ) / 2
    bᶜᶜᶜ = (bᶜᶜᶠ[:, 1:end-1] .+ bᶜᶜᶠ[:, 2:end]) ./ 2
    return -aᶜᶜᶜ .+ bᶜᶜᶜ
end

@inline function xz_state_data(foldername; σ=0, field="v", Δ=false)
    # Read data from files to be plotted
    #foldername = "../scratch/filament-instability/$runname"
    @info "Reading data from file"
    
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
    ts = [file["timeseries/t/$frame"] for frame in frames] .- 1
    
    @info "Creating streamfunction"
    ψs = field in ["ω", "w′v′"] ? nothing : map(frames) do frame
        u = file["timeseries/u_dfm/$frame"][:, 1, :]
        w = file["timeseries/w_dfm/$frame"][:, 1, :]
        imfilter(ψᶜᶜᶜ(u, w, Δx, Δzᵃᵃᶜ), gaussian((σ, 0)), "circular")[205:end-205, 33:end]
    end
    
    @info "Creating buoyancy"
    bs = map(frames) do frame
        imfilter(file["timeseries/b_dfm/$frame"][:, 1, :], gaussian((σ, 0)), "circular")[205:end-205, 33:end]
    end
    
    @info "Creating heatmap data"
    plot_datas = map(frames) do frame
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
        imfilter(plt, gaussian((σ, 0)), "circular")[205:end-205, 33:end]
    end
    close(file)
    
    figtitle = "Ro=$(round(sp.Ro; digits=1)), Ri=$(round(sp.Ri; digits=2))"
    
    return (; figtitle, ts, xs=xᶜᵃᵃ[205:end-205], zs=zᵃᵃᶜ[33:end], plot_datas, bs, ψs)
end

@inline function xz_state_video(foldername, outputfoldername; resolution=(800, 800), σ=0, field="v", Δ=false, cmax=nothing, t_max=25, axis_kwargs...)
    
    if ispath(outputfoldername)
        @info "$outputfoldername already exists!"
        return nothing
    end
    
    figtitle, ts, xs, zs, plot_datas, bs, ψs = xz_state_data(foldername; σ, field, Δ)
    
    @info "Generating scene"
    fig = Figure(; resolution, backgroundcolor = (:white, 0.0))
    
    n = Observable(101)
    n_max = argmin(abs.(ts .- t_max))
    t = @lift ts[$n]
    title = @lift "$figtitle, t = $(round($t; digits=2))"
    plot_data = @lift plot_datas[$n]
    b = @lift bs[$n]
    ψ = ψs == nothing ? nothing : @lift ψs[$n] 
    
    # The steps between contours
    ψstep = 6e-4
    bstep = 1.875
    ψ_min = minimum(minimum.(ψs))
    ψ_max = maximum(maximum.(ψs))
    b_min = minimum(minimum.(bs))
    b_max = maximum(maximum.(bs))
    ψrange = ψ == nothing ? nothing : ψ_min:ψstep:ψ_max
    brange = b_min:bstep:b_max
    
    # Maximum heatmap value
    cmax = cmax == nothing ? maximum([maximum(abs.(a)) for a in plot_datas]) : cmax
    
    ax = Axis(fig[1, 1]; title, xlabel="x", ylabel="z", limits=(xs[1], xs[end], zs[1], zs[end]),  axis_kwargs...)
    
    ht = heatmap!(ax, xs, zs, plot_data; colormap=:balance, colorrange=(-cmax, cmax))
    ψ != nothing && contour!(ax, xs, zs, ψ; colormap=:BrBG_10, levels=ψrange, alpha=1, linewidth=1)
    
    Colorbar(fig[1, 2], ht, label=Δ ? L"\overline{%$field} - \overline{%$field}_0" : L"\overline{%$field}")
    
    # Record the figure
    @info "Saving frames at $outputfoldername"
    mkpath(outputfoldername)
    
    w = ZipFile.Writer("$outputfoldername.zip")
    for i in 101:n_max
        print("\r$(100 * round((i-101) / n_max; digits=4))%")
        n[] = i
        filename = "$outputfoldername/$(lpad(i, 4, '0')).png"
        save(filename, fig; px_per_unit=2)
        
        zipfile = ZipFile.addfile(w, "$(lpad(i, 4, '0')).png")
        open(r -> write(zipfile, r), filename)
        rm(filename)
    end
    rm(outputfoldername)
    close(w)
    
    return nothing
end
