#= 
make_forcing.jl
    Functions to extract the vertical turbulent momentum flux due to SI
=#

@inline function get_vwFLUX(foldername, tmax)
    vfilename = "down_front.jld2"
    filename = "down_front_mean.jld2"
    paramfilename = "parameters.jld2"
    frames, grid = jldopen("$foldername/$filename") do file
        keys(file["timeseries/t"])[101:1101], file["serialized/grid"]
        end;
    xᶜᵃᵃ = xnodes(Center, grid)
    xᶠᵃᵃ = xnodes(Face, grid)
    zᵃᵃᶜ = znodes(Center, grid)
    zᵃᵃᶠ = znodes(Face, grid)
    Δzᵃᵃᶜ = reshape(diff(zᵃᵃᶠ), 1, length(zᵃᵃᶜ))
    Δx = xᶠᵃᵃ[2] - xᶠᵃᵃ[1]

    function ∂z(v)
        let a = (circshift(v, (0, -1)) - circshift(v, (0, 1))) ./ (2Δzᵃᵃᶜ)
            a[:, 1] .= 0
            a[:, end] .= 0
            a
        end
    end
    
    sp = jldopen("$foldername/$paramfilename") do file
        file["parameters/simulation"]
    end
    vfile = jldopen("$foldername/$vfilename")
    file = jldopen("$foldername/$filename")
    function vwFLUX(frame)
        (vfile["timeseries/vwFLUX/$frame"][:, 1, 1:length(zᵃᵃᶜ)] .+ vfile["timeseries/vwFLUX/$frame"][:, 1, 2:length(zᵃᵃᶜ)+1]) / 2
    end
    
    ts = [file["timeseries/t/$f"] for f in frames] .- 1
    
    vwFLUX_fit(t) = t < tmax ? 1 : 0
    
    δvwFLUX = sum(map(zip(frames, ts)) do (frame, t)
        vwFLUX(frame) * vwFLUX_fit(t) * (ts[2] - ts[1])
            end) / sum(VSP_fit.(ts))
    close(file)
    close(vfile)
    return (ts, δVSP, δψ, δvwFLUX)
end

function u_initial_condition(foldername)
    
end