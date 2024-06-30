using Oceananigans
using JLD2

foldername = ARGS[1]
scriptname = ARGS[2]
# Might use a lot of IO
buffer = length(ARGS) > 2 ? ARGS[3] : nothing
# Function to update input fields
function update_inputs!(input_fields, frame, path)
    jldopen(path) do file
        input_fields.u .= file["timeseries/u/$frame"]
        input_fields.v .= file["timeseries/v/$frame"]
        input_fields.w .= file["timeseries/w/$frame"]
        input_fields.b .= file["timeseries/b/$frame"]
        input_fields.φ .= file["timeseries/φ/$frame"]
        map(fill_halo_regions!, input_fields)
    end
    return nothing
end

macro postpostprocess() end
temp_outputs = (; )

function write_grid_times(grid, frames, ts, path)
    jldopen(path, "a") do file
        for (i, frame) in enumerate(frames)
            file["timeseries/t/$frame"] = ts[i]
        end
        # Now copy over grid things so Oceananigans isn't needed
        for k in fieldnames(typeof(grid))
            file["grid/$k"] = getproperty(grid, k)
        end
    end
    return nothing
end

function write_outputs(filename, outputs, iteration)
    jldopen(filename, "a") do file
        for (k, v) in zip(keys(outputs), outputs)
            file["timeseries/$k/$iteration"] = v.data
        end
    end
end
@info "Reading timeseries from file"
u_series = FieldTimeSeries("$foldername/output.jld2", "u"; backend=OnDisk())
v_series = FieldTimeSeries("$foldername/output.jld2", "v"; backend=OnDisk())
w_series = FieldTimeSeries("$foldername/output.jld2", "w"; backend=OnDisk())

b_series = FieldTimeSeries("$foldername/output.jld2", "b"; backend=OnDisk())

φ_series = FieldTimeSeries("$foldername/output.jld2", "φ"; backend=OnDisk())

ν_series = FieldTimeSeries("$foldername/output.jld2", "ν"; backend=OnDisk())

u = u_series[1]
v = v_series[1]
w = w_series[1]

b = b_series[1]

φ = φ_series[1]

ν = ν_series[1]
#=
@info u
@info v
@info w
@info b
@info φ
@info ν
=#
grid, iters, ts = jldopen("$foldername/output.jld2") do file
    iters = keys(file["timeseries/t"])
    file["serialized/grid"], iters, [file["timeseries/t/$iter"] for iter in iters]
end

sp, init_time, run_time = jldopen("$foldername/parameters.jld2") do file
    file["parameters/simulation"], file["parameters/init_time"], file["parameters/run_time"]
end

# Include script file which should define a named tuple of outputs, temp_outputs (which are deleted after)
# and a function called update_outputs!
@info "Including $scriptname.jl"
include("$scriptname.jl")

@info outputs

filename = buffer == nothing ? "$foldername/$scriptname.jld2" : "$buffer/$scriptname.jld2"
temp_filename = buffer == nothing ? "$foldername/temp_$scriptname.jld2" : "$buffer/temp_$scriptname.jld2"

for (i, iter) in enumerate(iters)
    print("Computing $iter\r")
    u .= u_series[i]
    v .= v_series[i]
    w .= w_series[i]

    b .= b_series[i]

    φ .= φ_series[i]

    ν .= ν_series[i]
    update_outputs!(outputs)
    write_outputs(filename, outputs, iter)
    write_outputs(temp_filename, temp_outputs, iter)
end
println()
write_grid_times(grid, iters, ts, filename)

@postpostprocess
rm(temp_filename)

if buffer != nothing
    @info "Moving from $buffer to $foldername"
    mv(filename, "$foldername/$scriptname.jld2"; force=true)
end

@info "Finished!"
