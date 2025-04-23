using DataFrames
using CSV
using Statistics

"""
Default values of system parameters.
- c_i - price of buying energy from the grid, Euro/kW
- c_o - price of selling energy to the grid, Euro/kW
- asset_lifetime - expected lifetime of a component in years. Needed to bring investment and operational costs to the same timescale
- c_pv, c_wind - price of installing pv and wind components per kW peak, Euro/kWp
- c_storage - price of storage, Euro/kW
- c_heat_storage - TODO: heat layer units
- c_heatpump
- inv_budget - maximum investment budget, Euro
- recovery_time - time after the event, in which the system is allowed to deviate from background optimal schedule, h
- COP - heatpump efficiency coefficient
- heat_losses - losses in the heat storage
- penalty - price of non-delivery of flexibility, Euro
"""

# PyPSA Heat: https://github.com/PyPSA/technology-data/blob/master/outputs/costs_2020.csv
# central water tank charger,efficiency,1.0
# central water tank discharger,efficiency,1.0
# central water tank storage,investment,3.1374,EUR/kWhCapacity,
# central water tank storage,lifetime,40.0,years,
# central air-sourced heat pump,efficiency,3.1,per unit,
# central air-sourced heat pump,investment,1006.7765,EUR/kW_th,
# central air-sourced heat pump,lifetime,25.0,years,
# heat losses per hour = 0.0001

if !@isdefined default_es_pars
    default_es_pars = Dict((:c_i => 0.3,
        :c_o => 0.05,
        :asset_lifetime => 20.0,
        :c_pv => 700.0,
        :c_wind => 2000.0,
        :c_storage => 600.0,
        :c_heat_storage => 3.1 * 0.5, # lifetime factor
        :c_heatpump => 1006.0 * 0.8, # lifetime factor
        :inv_budget => 500000000.0,
        :recovery_time => 4,
        :COP => 3.1,
        :feedincap => 10^7,
        :heat_losses => 0.0001,
        :heat_eff => 1.0,
        :heat_demand_scale => 1.0,
        :sto_ef_ch => 0.95,
        :sto_ef_dis => 0.95,
        :penalty => 10000.0,
        :expected_events_per_full_period => 1,
        :max_sto_flow => 0.8,
        :max_pv => 10^3,
        :max_wind => 10^3))
end

function load_max_boegl(;
    basepath=joinpath(@__DIR__, ".."),
    time_steps=nothing,
    offset=0,
    heat="when2heat",
    scale_factor=1.0)
    pv_data = CSV.read(joinpath(basepath, "timeseries/", "pv_Halle18.csv"),
        DataFrame,
        header=false)
    wind_data = CSV.read(joinpath(basepath, "timeseries/", "wind_Karholz.csv"),
        DataFrame,
        header=false)
    demand_data = CSV.read(joinpath(basepath,
            "timeseries/",
            "demand_Industriepark.csv"),
        DataFrame,
        header=false)

    if isnothing(time_steps)
        timesteps = 1:length(pv_data[:, 1])
    else
        timesteps = time_steps
    end

    if heat == "when2heat"
        heatdemand_data = CSV.read(joinpath(basepath, "timeseries", "heatdemand.csv"),
            DataFrame)
        heatdemand = heatdemand_data[timesteps.+offset, 1] ./ (50.0 * scale_factor)
    else
        heatdemand = zeros(length(timesteps))
    end
    pv = pv_data[timesteps.+offset, 1]
    wind = wind_data[timesteps.+offset, 1]
    demand = demand_data[timesteps.+offset, 1]

    pars = copy(default_es_pars)

    pars[:recovery_time] = 4
    pars[:c_storage] = 300.0
    pars[:c_pv] = 800.0
    pars[:c_wind] = 1150.0
    pars[:c_in] = 0.165
    pars[:c_out] = 0.02
    pars[:inv_budget] = 10000000.0
    pars[:sto_ef_ch] = 0.97
    pars[:sto_ef_dis] = 0.97
    pars[:feedincap] = 1e7
    pars[:max_pv] = 2 * 10^4 # increased by factor of 2
    pars[:max_wind] = 10^4
    pars[:inv_budget] = 10^10
    println("Hourly electricity demand: min = $(minimum(demand)),  mean = $(mean(demand)), max = $(maximum(demand))")
    println("Hourly heat demand: min = $(minimum(heatdemand)),  mean = $(mean(heatdemand)), max = $(maximum(heatdemand))")
    println("Hourly PV production at max. investment: mean = $(mean(pv)*pars[:max_pv]), max = $(maximum(pv)*pars[:max_pv])")
    println("Hourly wind production at max. investment: mean = $(mean(wind)*pars[:max_wind]), max = $(maximum(wind)*pars[:max_wind])")
    return (pv=pv,
        wind=wind,
        demand=demand,
        heatdemand=heatdemand,
        pars=pars)
end

function load_basic_example(timesteps;
    basepath=joinpath(@__DIR__, ".."),
    offset=0,)
    data = CSV.read(joinpath(basepath, "timeseries", "basic_example.csv"), DataFrame)
    heatdemand_data = CSV.read(joinpath(basepath, "timeseries", "heatdemand.csv"),
        DataFrame)

    pv = data[timesteps.+offset, 3]
    wind = data[timesteps.+offset, 4]
    demand = data[timesteps.+offset, 2]
    heatdemand = heatdemand_data[timesteps.+offset, 1]

    data = nothing
    heatdemand_data = nothing # Free the memory

    pars = copy(default_es_pars)

    pars[:recovery_time] = 3
    pars[:c_storage] = 100.0
    pars[:c_pv] = 300.0
    pars[:c_wind] = 550.0
    pars[:penalty] = 1000000.0
    return (pv=pv,
        wind=wind,
        demand=demand,
        heatdemand=heatdemand,
        pars=pars)
end
