using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using Serialization

include(joinpath(@__DIR__, "../src/conventions.jl"))

skip_df_prep = false

if ! @isdefined exp_dir
    exp_dir = "../exp_results/"
end

if (@isdefined _df) && _df !== nothing
    # Assume preamble already ran
    println("_df is defined, assuming preamble already ran, data is not being reloaded.")
    skip_df_prep = true
elseif (@isdefined df_file) && (df_file !== nothing) && (isfile(df_file)) # load date from local file
    println("Loading results from $(df_file)")
    include(joinpath(@__DIR__, "../src/new_eval_utils.jl"))
    _df = deserialize(df_file)
else
    println("Aggregating results from $(exp_dir)")
    include(joinpath(@__DIR__, "../src/new_eval_utils.jl"))
    _df = load_results(joinpath(@__DIR__, exp_dir))
    df_file = joinpath(@__DIR__, "../preloaded", "df_raw.serial")
    mkpath(joinpath(@__DIR__, "../preloaded"))
    serialize(df_file, _df)
end

if ! skip_df_prep
    const sm = skipmissing

    sum_inm(x) = ismissing(x) ? x : sum(x)
    l_inm(x) = ismissing(x) ? x : length(x)
    im_to_0(x) = ismissing(x) ? 0. : x
    max_inm(x) = ismissing(x) ? x : maximum(x)
    _df.total_penalty = sum_inm.(_df.erp .+ _df.erm)

    first_dir_after_results(path) = splitpath(split(path, "results/")[end])[1]
    _df.first_dir = first_dir_after_results.(_df.path)

    _df.infeasible = _df.feasible .== false
    _df.optimized = _df.termination_status .=== "OPTIMAL"

    get_F_idx(x) = ismissing(x) ? x : F_dict[x]

    _df.F_idx = get_F_idx.(_df.F)

    _df.total_n_events = _df.t_xi .|> l_inm

    # + penalty * 0.1 * (sum(hfp) + sum(hfm))
    # + penalty * (sum(ebp[1]) + sum(ebm[1]))
    # + penalty * 10. * (sum(ebp[2:end]) + sum(ebm[2:end]))
    # + penalty * 10. * (sum(hbp) + sum(hbm))

    _df.energy_not_served = im_to_0.(sum_inm.(_df.erp .+ _df.erm))

    _df.penalty_cost = im_to_0.((_df.expected_events_per_full_period ./ _df.total_n_events) .*_df.penalty .* _df.total_penalty)

    _df.no_penalty_cost = _df.full_cost .- _df.penalty_cost

    _df.pathological = ( _df.gci .|> max_inm) .> 5e6
end

df = copy(_df)

if df.pathological[ .! df.infeasible] |> any
    println("Experimental runs with pathological behavior are present.")
end

if df.infeasible |> any
    println("$(count(df.infeasible)) infeasible")
end

# deleteat!(df, findall(df.infeasible))