using Serialization
using Gurobi
using StatsBase

include("sp_model.jl")

function get_sp_bkg(data_es;
        silent = true, savepath = nothing, optimizer = Gurobi.Optimizer)
    pv, wind, demand, heatdemand, pars = data_es

    sp_bkg = define_energy_system(pv, wind, demand, heatdemand, pars;
        capacity_flex_constraint = false,
        F_pos = 0.0,
        F_neg = 0.0,
        optimizer)

    silent && set_silent(sp_bkg)

    if !isnothing(savepath)
        serialize(joinpath(savepath, "background_blank.serial"), sp_bkg)
    end

    return sp_bkg
end

function get_sp_cap_only(data_es,
        F;
        silent = true,
        savepath = nothing,
        optimizer = Gurobi.Optimizer)
    pv, wind, demand, heatdemand, pars = data_es

    sp_cap_only = define_energy_system(pv, wind, demand, heatdemand, pars;
        capacity_flex_constraint = true,
        F_pos = F,
        F_neg = -F,
        optimizer)

    silent && set_silent(sp_cap_only)

    if !isnothing(savepath)
        serialize(joinpath(savepath, "capacity_only_blank.serial"), sp_cap_only)
    end

    return sp_cap_only
end

function generate_sample(data_es; F, flex_interval, n_samples)
    pars = data_es.pars
    l_ts = length(data_es.pv)
    # Start no event less than the recovery time from the end of the timeseries
    t_latest = l_ts - pars[:recovery_time] - 1
    # the time between requests can be no shorter than the recovery time
    delta_t = flex_interval - pars[:recovery_time]
    # We need to know the expected number of requests per time window
    expected_events_per_full_period = t_latest / (delta_t + pars[:recovery_time] + 1)

    #TODO: Think about the + 1 here.

    total_flex_reqs = round(Int, expected_events_per_full_period * n_samples)

    F_min = 0.6 * F
    F_mean = 0.5 * (F + F_min)

    scens = non_overlapping_sample(total_flex_reqs,
        F,
        t_latest;
        F_min = 0.6 * F)

    sample_data = (total_flex_reqs = total_flex_reqs,
        flex_interval = flex_interval,
        n_samples = n_samples,
        F = F,
        F_mean = F_mean,
        expected_events_per_full_period = expected_events_per_full_period,
        expected_flex_energy_per_full_period = expected_events_per_full_period * F_mean,
        expected_flex_energy_per_hour = expected_events_per_full_period * F_mean / l_ts)
    return sample_data, scens
end

function get_sp_full(data_es;
        sample_data,
        scens,
        silent = true,
        savepath = nothing,
        optimizer = Gurobi.Optimizer,
        nfsc = false)
    pv, wind, demand, heatdemand, pars = data_es

    if isnothing(scens) || isnothing(sample_data)
        error("Scen or sample_data missing. Call generate_sample(data_es; F, flex_interval, n_samples) and supply the outputs!")
    end

    sp_full = define_energy_system(pv, wind, demand, heatdemand, pars;
        events_per_period = sample_data.expected_events_per_full_period,
        scenarios=scens,
        capacity_flex_constraint = true,
        F_pos = sample_data.F,
        F_neg = -sample_data.F)

    silent && set_silent(sp_full)

    if !isnothing(savepath)
        serialize(joinpath(savepath, "full_blank.serial"), sp_full)
    end

    return sp_full
end
