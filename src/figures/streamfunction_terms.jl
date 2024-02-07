using Oceananigans
using JLD2
using CairoMakie
using Statistics
using ImageFiltering: imfilter, Kernel.gaussian

@inline function ψterms(foldername; σ=0)
    #foldername = "../scratch/filament-instability/$runname"
    filename = "down_front_mean.jld2"
    vfilename = "down_front.jld2"
    ufilename = "across_front.jld2"
    wfilename = "vertical.jld2"
    bfilename = "buoyancy.jld2"
    paramfilename = "parameters.jld2"
    frames, grid = jldopen("$foldername/$filename") do file
        keys(file["timeseries/t"])[101:2601], file["serialized/grid"]
        end;
    xᶜᵃᵃ = xnodes(Center, grid)
    xᶠᵃᵃ = xnodes(Face, grid)
    zᵃᵃᶜ = znodes(Center, grid)
    zᵃᵃᶠ = znodes(Face, grid)
    Δzᵃᵃᶜ = reshape(diff(zᵃᵃᶠ), 1, length(zᵃᵃᶜ))
    Δx = xᶠᵃᵃ[2] - xᶠᵃᵃ[1]
    function ψᶜᶜᶜ(uᶠᶜᶜ, wᶜᶜᶠ, xᶜᵃᵃ, xᶠᵃᵃ, zᵃᵃᶜ, zᵃᵃᶠ)
        # Integrate
        aᶠᶜᶜ = cumsum(uᶠᶜᶜ .* Δzᵃᵃᶜ; dims=2)

        bᶜᶜᶠ = cumsum(wᶜᶜᶠ .* Δx; dims=1)
        aᶜᶜᶜ = (circshift(aᶠᶜᶜ, (-1, 0)) .+ aᶠᶜᶜ) / 2
        bᶜᶜᶜ = (bᶜᶜᶠ[:, 1:end-1] .+ bᶜᶜᶠ[:, 2:end]) ./ 2
        return -aᶜᶜᶜ .+ bᶜᶜᶜ
    end
    function ∂z(a)
        let a = -(circshift(a, (0, 1)) - a) ./ Δzᵃᵃᶜ
            a[:, 1] .= 0
            a
        end
    end
    function ∂x(a)
        (circshift(a, (-1, 0)) - circshift(a, (1, 0))) ./ 2Δx
    end
    ∇²(a) = ∂x(∂x(a)) + ∂z(∂z(a))
    sp = jldopen("$foldername/$paramfilename") do file
        file["parameters/simulation"]
    end
    z_omit_fraction=0.1
    # Boundary layer cells
    blc = zᵃᵃᶜ .> -sp.H
    # Central boundary layer cells
    cblc = -sp.H*z_omit_fraction .> zᵃᵃᶜ .> -sp.H * (1-z_omit_fraction)
    
    axtitle = "Ro=$(round(sp.Ro; digits=1)), Ri=$(round(sp.Ri; digits=2))"
    file = jldopen("$foldername/$filename")
    ufile = jldopen("$foldername/$ufilename")
    vfile = jldopen("$foldername/$vfilename")
    wfile = jldopen("$foldername/$wfilename")
    bfile = jldopen("$foldername/$bfilename")
    ts = [file["timeseries/t/$f"] for f in frames] .- 1
    Δt = ts[2] - ts[1]
    v(i) = file["timeseries/v_dfm/$(frames[i])"][:, 1, :]
    u(i) = file["timeseries/u_dfm/$(frames[i])"][:, 1, :]
    w(i) = file["timeseries/w_dfm/$(frames[i])"][:, 1, :]
    b(i) = file["timeseries/b_dfm/$(frames[i])"][:, 1, :]
    Fx(i) = (let a=ufile["timeseries/uwFLUX/$(frames[i])"]; a[:, 1, 2:end] + a[:, 1, 1:end-1] end) / 2
    w′v′(i) = -cumsum((let a=vfile["timeseries/vwFLUX/$(frames[i])"]; a[:, 1, 2:end] + a[:, 1, 1:end-1] end) .* Δzᵃᵃᶜ / 2; dims=2)
    Fz(i) = wfile["timeseries/wwFLUX/$(frames[i])"][:, 1, :]
    w′b′(i) = -cumsum((let a=bfile["timeseries/bwFLUX/$(frames[i])"]; a[:, 1, 2:end] + a[:, 1, 1:end-1] end) .* Δzᵃᵃᶜ / 2; dims=2)
    ψ(i) = ψᶜᶜᶜ(u(i), w(i), xᶜᵃᵃ, xᶠᵃᵃ, zᵃᵃᶜ, zᵃᵃᶠ)
    J₁₁(i) = ∂z(sp.f * v(1)) .* ∂x(∂z(ψ(i)))
    J₁₂(i) = ∂x(sp.f * v(1)) .* ∂z(∂z(ψ(i)))
    J₂₁(i) = ∂z(b(1)) .* ∂x(∂x(ψ(i)))
    J₂₂(i) = ∂x(b(1)) .* ∂z(∂x(ψ(i)))
    ∂²ψ∂t²(i) = (ψ(i+1) + ψ(i-1) - 2ψ(i)) / Δt^2
    data = map(2:length(frames)-1) do i
       map([∇²(∂²ψ∂t²(i)), ∂z(∂z(sp.f^2 * ψ(i))), -J₁₁(i), J₁₂(i), J₂₁(i), -J₂₂(i), -sp.f*∂z(∂z(w′v′(i))), ∂x(∂z(w′b′(i)))]) do field
            sum((field * Δx .* Δzᵃᵃᶜ)[1:length(xᶜᵃᵃ)÷2, cblc])
        end
    end
    # Non linear terms
    δJ₁₁(i) = ∂z(sp.f * (v(i)-v(1))) .* ∂x(∂z(ψ(i)))
    δJ₁₂(i) = ∂x(sp.f * (v(i)-v(1))) .* ∂z(∂z(ψ(i)))
    δJ₂₁(i) = ∂z(b(i)-b(1)) .* ∂x(∂x(ψ(i)))
    δJ₂₂(i) = ∂x(b(i)-b(1)) .* ∂z(∂x(ψ(i)))
    δatw₁(i) = ∂z(sp.f * ∂z(v(i)) - ∂x(b(i))) .* ∂x(ψ(i))
    δatw₂(i) = ∂x(sp.f * ∂z(v(i)) - ∂x(b(i))) .* ∂z(ψ(i))
    δ∂tFψ(i) = let a(i)=∂x(Fz(i)) - ∂z(Fx(i)); (a(i+1) - a(i)) / Δt end
    data_nl = map(2:length(frames)-1) do i
       map([-δJ₁₁(i), δJ₁₂(i), δJ₂₁(i), -δJ₂₂(i), -δatw₁(i), δatw₂(i), -δ∂tFψ(i)]) do field
            sum((field * Δx .* Δzᵃᵃᶜ)[1:length(xᶜᵃᵃ)÷2, cblc])
        end
    end
    close(file)
    close(ufile)
    close(vfile)
    close(wfile)
    close(bfile)
    ts = ts[2:end-1]
    data = [data[i][j] for j in 1:length(data[1]), i in 1:length(data)]
    data_nl = [data_nl[i][j] for j in 1:length(data_nl[1]), i in 1:length(data_nl)]
    linear_terms = sum(data[2:end, :]; dims=1)[1, :]
    nonlinear_terms = sum(data_nl; dims=1)[1, :]
    ∂²C∂t² = -data[1, :]
    ∂²C∂t², linear_terms, nonlinear_terms = map(x->imfilter(x, gaussian((σ, )), "replicate"), (∂²C∂t², linear_terms, nonlinear_terms))
    return (; ts, ∂²C∂t², linear_terms, nonlinear_terms, axtitle)
end

@inline function ψterms!(layout_cell; ts, ∂²C∂t², linear_terms, nonlinear_terms, axtitle)
    axis_kwargs = (;
        xlabel = "t",
        ylabel = L"\frac{\text{d}^2C}{\text{d}t^2}",
        title = axtitle,
        limits = (0, ts[end], -1, 1)
    )
    ax = Axis(layout_cell; axis_kwargs...)
    ln1 = lines!(ax, ts, ∂²C∂t²; color=:black)
    ln2 = lines!(ax, ts, linear_terms)
    ln3 = lines!(ax, ts, nonlinear_terms)
    ln4 = lines!(ax, ts, linear_terms .+ nonlinear_terms)
    return [ln1, ln2, ln3, ln4]
end

@inline function ψterms_plot(runnames; resolution=(1000, 500),  σ=0)
    n_plots = length(runnames)
    plot_datas = ψterms.(runnames; σ)
    fig = Figure(; resolution, backgroundcolor = (:white, 0.0))
    lnss = map(enumerate(plot_datas)) do (i, plot_data)
        ψterms!(fig[1, i]; plot_data...)
    end
    Legend(fig[1, n_plots+1], lnss[1], ["Actual", "Linear terms", "Non-linear terms", "Sum"])
    if n_plots > 1
        hideydecorations!.(fig.content[2:n_plots])
        for ax in fig.content[2:n_plots]
            ax.ygridvisible = true
        end
    end
    colgap!(fig.layout, 15)
    fig
end