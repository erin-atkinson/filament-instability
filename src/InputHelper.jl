module InputHelper
export CLIOpt, readparams, parseargs, parseargs!, getargs

mutable struct CLIOpt{T}
    aliases::Tuple{Vararg{String}}
    value::T
end
# For some reason this won't work with {Any} specified..? So I need wrap()
CLIOpt(alias::String, value) = CLIOpt((alias,), value)

wrap(x) = (x,)
wrap(x::Tuple) = x

function readparams(filename)
    """
    Read parameters from a file, return a tuple of named tuples
    """
    params_str  = open(filename) do io
        parameter_text = chomp(replace(read(io, String), r"#.*?\n" => ""))
        params_str = split.(split(parameter_text, "\n\n"), "\n")
    end
    return Tuple((;zip(Symbol.(replace.(i, r"=.*" => "")), Base.Meta.parse.(replace.(i, r".*?=" => "")))...) for i in params_str)
end

function getargs(args; kwargs...)
    kwargs = NamedTuple(kwargs)
    cli_opts = map((x) -> CLIOpt{Any}(wrap(x), nothing), kwargs)
    p_args = parseargs!(cli_opts, args)
    kw_args = map((x) -> x.value, cli_opts)
    return (p_args, kw_args)
end

function parseargs!(cli_opts, args)
    """
    Parse all CLI arguments, return a tuple of positional arguments, and modify the input options in cli_opts
    """
    args = Iterators.Stateful(args)
    p_args = []
    for k in args
        if startswith(k, r"--?[A-Z]")
            # It is a flag
            for x in cli_opts
                if k in x.aliases
                    x.value = true
                end
            end
        elseif startswith(k, r"--?[a-z]")
            # It is a key, value pair, get the value from the iterator
            v = popfirst!(args)
            for x in cli_opts
                if k in x.aliases
                    x.value = parse_str(v)
                end
            end
        else
            # It is a positional argument
            append!(p_args, (parse_str(k), ))
        end
    end
    return Tuple(p_args)
end

function parse_str(str)
    try
        v = Base.Meta.parse(str)
        (v isa Symbol || v isa Expr) ? str : v
    catch e
        e isa Base.Meta.ParseError && str
    end
end

end
