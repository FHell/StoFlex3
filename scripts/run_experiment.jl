using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

println(@__DIR__)

using LinearAlgebra
println("BLAS Threads: $(LinearAlgebra.BLAS.get_num_threads())")
println("Julia Threads: $(Threads.nthreads())")

LinearAlgebra.BLAS.set_num_threads(1)

include(joinpath(@__DIR__, "../src/load_energy_system_data.jl"))
include(joinpath(@__DIR__, "../src/sp_model.jl"))
include(joinpath(@__DIR__, "../src/prepare_stochastic_programs.jl"))
include(joinpath(@__DIR__, "../src/conventions.jl"))

using CSV
using DelimitedFiles

#-
# TODO dual_objective value does not work for deserialized models, and neither does the
# ScenarioDependentModelAttribute, load model has to be re-optimized (probably quicker)
# TODO Think about how to control the random seed
#-

# ENV["GUROBI_HOME"] = "/home/hellmann/gurobi"
# ENV["GUROBI_JL_USE_GUROBI_JLL"] = "false"

using Gurobi
using Random
Random.seed!(42)


_ARGS = []

if length(ARGS) == 0 # running locally
    using Gurobi
    OUR_OPTIMIZER = Gurobi.Optimizer
    F = F_array[5]
    flex_interval = flex_interval_array[3]
    n_samples = n_samples_array[1]
    data_es = load_industry_park1(time_steps = 1:1000)
    F = 500
    flex_interval = 6
    n_samples = 10
    data_es = load_industry_park1()
end

if length(ARGS) >= 3 # running on cluster
    _ARGS = parse.(Int64, ARGS[1:3])

    # using Gurobi
    using Gurobi
    OUR_OPTIMIZER = Gurobi.Optimizer
    F = F_array[_ARGS[1]]
    flex_interval = flex_interval_array[_ARGS[2]]
    n_samples = n_samples_array[_ARGS[3]]
    data_es = load_industry_park1()
end

flex_pars = (F = F, flex_interval = flex_interval, n_samples = n_samples)
println("Running experiments with $flex_pars.")

if length(ARGS) >= 4
    resultspath = joinpath(@__DIR__, "../exp_results/", ARGS[4])
    mkpath(resultspath)
elseif length(ARGS) > 0
    resultspath = joinpath(@__DIR__, "../exp_results/", "defaults/")
    mkpath(resultspath)
else
    resultspath = joinpath(@__DIR__, "../exp_results/", "local/")
    mkpath(resultspath)
end

println("Saving results to $resultspath")

add_rp(s) = joinpath(resultspath, s)

bkgpath = nothing
cappath = nothing
sppath = add_rp("F_$(F)/flex_$(flex_interval)/n_samples_$(n_samples)/")
mkpath(sppath)

if length(_ARGS) == 0 || all(_ARGS .== 1)
    bkgpath = resultspath
end

if length(_ARGS) == 0 || all(_ARGS[2:end] .== 1)
    cappath = add_rp("F_$(F)/")
end

if length(ARGS) >= 5
    for arg in ARGS[5:end]
        par_name, par_val = split(arg, "=")
        println("Setting $par_name to $par_val")
        data_es.pars[Symbol(par_name)] = parse(Float32, par_val)
    end
end

if data_es.pars[:heat_demand_scale] != 1.0
    println("Rescaling heat demand by $(data_es.pars[:heat_demand_scale])")
    data_es.heatdemand .*= data_es.pars[:heat_demand_scale]
end

## Background model
println("Optimizing the background model...")

# To avoid faile access conflicts when using job arrays
# only the first job saves data

sp_bkg = get_sp_bkg(data_es; savepath = bkgpath)

set_attribute(sp_bkg, MOI.NumberOfThreads(), 4)

optimize!(sp_bkg) # Jump model is optimized without cache

save_no_scen_results(bkgpath, sp_bkg, data_es.pars;
        name = "background",)

## Capacity only model

println("Optimizing the capacity only model...")

sp_cap_only = get_sp_cap_only(data_es, flex_pars.F;
savepath = cappath,
optimizer = OUR_OPTIMIZER,)

set_attribute(sp_cap_only, MOI.NumberOfThreads(), 4)


optimize!(sp_cap_only) # Jump model is optimized without cache

save_no_scen_results(cappath, sp_cap_only, data_es.pars; name = "capacity_only", F=flex_pars.F)

## SP models

# ini_sample_data == re_sample_data
ini_sample_data, ini_scens = generate_sample(data_es; flex_pars...)
re_sample_data, re_scens = generate_sample(data_es; flex_pars...)

sp_ini = get_sp_full(data_es;
    sample_data = ini_sample_data,
    scens = ini_scens,
    savepath = sppath,
    optimizer = OUR_OPTIMIZER,)

sp_re = get_sp_full(data_es;
    sample_data = re_sample_data,
    scens = re_scens,
    savepath = sppath,
    optimizer = OUR_OPTIMIZER,)

set_attribute(sp_ini, MOI.NumberOfThreads(), 4)
set_attribute(sp_re, MOI.NumberOfThreads(), 4)

##



## OFIOR
println("Optimizing the OFIOR model...")

ini_scens_objs, re_scens_objs = optimize_resample!(sp_ini, sp_re)

save_results(sppath, sp_ini, data_es.pars;
    flex_pars...,
    sample_data = ini_sample_data,
    scens = ini_scens,
    scens_objectives = ini_scens_objs,
    name = "OFIOR_ini")
save_results(sppath, sp_re, data_es.pars;
    flex_pars...,
    sample_data = re_sample_data,
    scens = re_scens,
    scens_objectives = re_scens_objs,
    name = "OFIOR")

# OFOR
println("Optimizing the OFOR model...")

fix_investment!(sp_ini, get_investments(sp_cap_only))

ini_scens_objs, re_scens_objs = optimize_resample!(sp_ini, sp_re)

save_results(sppath, sp_ini, data_es.pars;
    flex_pars...,
    sample_data = ini_sample_data,
    scens = ini_scens,
    scens_objectives = ini_scens_objs,
    name = "OFOR_ini")
save_results(sppath, sp_re, data_es.pars;
    flex_pars...,
    sample_data = re_sample_data,
    scens = re_scens,
    scens_objectives = re_scens_objs,
    name = "OFOR")

## OFR
println("Optimizing the OFR model...")

fix_operation!(sp_ini, get_operation(sp_cap_only))

ini_scens_objs, re_scens_objs = optimize_resample!(sp_ini, sp_re)

save_results(sppath, sp_ini, data_es.pars;
    flex_pars...,
    sample_data = ini_sample_data,
    scens = ini_scens,
    scens_objectives = ini_scens_objs,
    name = "OFR_ini")
save_results(sppath, sp_re, data_es.pars;
    flex_pars...,
    sample_data = re_sample_data,
    scens = re_scens,
    scens_objectives = re_scens_objs,
    name = "OFR")

## Stochastic optimization of background model
# This is somewhat more complex due to infeasibilities, hence unpack optimize_resample
println("Optimizing the background investment with scenarios...")

sp_ini = get_sp_full(data_es;
    sample_data = ini_sample_data,
    scens = ini_scens,
    savepath = sppath,
    optimizer = OUR_OPTIMIZER,)

sp_re = get_sp_full(data_es;
    sample_data = re_sample_data,
    scens = re_scens,
    savepath = sppath,
    optimizer = OUR_OPTIMIZER,)

set_attribute(sp_ini, MOI.NumberOfThreads(), 4)
set_attribute(sp_re, MOI.NumberOfThreads(), 4)
    

fix_investment!(sp_ini, get_investments(sp_bkg))


ini_scens_objs, re_scens_objs = optimize_resample!(sp_ini, sp_re)

save_results(sppath, sp_ini, data_es.pars;
    flex_pars...,
    sample_data = ini_sample_data,
    scens = ini_scens,
    scens_objectives = ini_scens_objs,
    name = "stochastic_background_ini")
save_results(sppath, sp_re, data_es.pars;
    flex_pars...,
    sample_data = re_sample_data,
    scens = re_scens,
    scens_objectives = re_scens_objs,
    name = "stochastic_background")

##
