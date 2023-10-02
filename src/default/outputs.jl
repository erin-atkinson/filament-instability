#=
default/outputs.jl
    Additional fields to output, none in this case
=#
using Oceananigans

# This macro should define and create output writers for any desired saved fields
macro additional_outputs!(simulation)
    return quote
        nothing
    end
end