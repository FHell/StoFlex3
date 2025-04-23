###############################################################################
# Energy system model without StochasticPrograms.jl                          #
#                                                                             #
# Author: ChatGPT                                                             #
# Date:   2025‑04‑21                                                          #
#                                                                             #
###############################################################################

using JuMP
using StatsBase
using LinearAlgebra

################################################################################
# Data structures & scenario sampling                                          #
################################################################################

mutable struct Scenario
    t_xi::Int
    s_xi::Int
    F_xi::Float64
    probability::Float64
end

function non_overlapping_sample(n, F_max, t_latest; F_min=0.0)
    times = repeat(1:t_latest, n ÷ t_latest)
    if n ÷ t_latest > 0
        @warn "More samples than time points – times will repeat."
    end
    append!(times, StatsBase.sample(1:t_latest, n % t_latest; replace=false))
    scenarios = Scenario[]
    for t in times
        push!(scenarios, Scenario(t, rand([-1, 1]), F_min + rand() * (F_max - F_min), 1 / n))
    end
    return scenarios
end

################################################################################
# Model builder                                                                #
################################################################################

inv_vars = [:u_pv, :u_wind, :u_storage, :u_heat_storage, :u_heatpump]

function define_energy_system(
    pv::AbstractVector{<:Real},
    wind::AbstractVector{<:Real},
    demand::AbstractVector{<:Real},
    heatdemand::AbstractVector{<:Real},
    pars;
    scenarios=[],
    debug_cap=1e9,
    reg_lossy_flows=1e-7,
    events_per_period=0,
    capacity_flex_constraint=true,
    F_pos=1000.0,
    F_neg=-1000.0,
    optimizer=Gurobi.Optimizer)

    number_of_hours = minimum((length(pv), length(wind), length(demand)))
    Nscen = length(scenarios)
    Trec = pars[:recovery_time] + 1
    lifetime_factor = pars[:asset_lifetime] * 365 * 24 / number_of_hours
    if Nscen > 0
        @assert isapprox(sum(s.probability for s in scenarios), 1.0; atol=1e-8)
    end

    model = Model(optimizer)

    #──────────────── First‑stage variables ────────────────#
    @variables(model, begin
        0 <= u_pv <= pars[:max_pv]
        0 <= u_wind <= pars[:max_wind]
        0 <= u_storage <= debug_cap
        0 <= u_heatpump <= debug_cap
        0 <= u_heat_storage <= debug_cap
        0 <= pv_cur[1:number_of_hours] <= debug_cap
        0 <= wind_cur[1:number_of_hours] <= debug_cap
        0 <= gci[1:number_of_hours] <= pars[:feedincap]
        0 <= gco[1:number_of_hours] <= pars[:feedincap]
        0 <= sto_to_bus[1:number_of_hours] <= debug_cap
        0 <= sto_from_bus[1:number_of_hours] <= debug_cap
        0 <= sto_soc[1:number_of_hours] <= debug_cap
        0 <= heat_sto_to_bus[1:number_of_hours] <= debug_cap
        0 <= heat_sto_from_bus[1:number_of_hours] <= debug_cap
        0 <= heat_sto_soc[1:number_of_hours] <= debug_cap
        0 <= flow_energy2heat[1:number_of_hours] <= debug_cap
    end)


    @constraint(model, [t = 1:number_of_hours], pv_cur[t] <= u_pv * pv[t])
    @constraint(model, [t = 1:number_of_hours], wind_cur[t] <= u_wind * wind[t])

    @constraint(model, [t = 1:number_of_hours-1],
        sto_soc[t+1] == sto_soc[t] + sto_from_bus[t] * pars[:sto_ef_ch] - sto_to_bus[t] / pars[:sto_ef_dis])

    @constraint(model, [t = 1:number_of_hours-1],
        heat_sto_soc[t+1] == heat_sto_soc[t] + heat_sto_from_bus[t] * pars[:heat_eff] - heat_sto_to_bus[t] / pars[:heat_eff])

    @constraint(model, [t = 1:number_of_hours], sto_soc[t] <= u_storage)
    @constraint(model, [t = 1:number_of_hours], heat_sto_soc[t] <= u_heat_storage)
    @constraint(model, [t = 1:number_of_hours], sto_from_bus[t] <= pars[:max_sto_flow] * u_storage)
    @constraint(model, [t = 1:number_of_hours], sto_to_bus[t] <= pars[:max_sto_flow] * u_storage)

    @constraint(model, [t = 1:number_of_hours], flow_energy2heat[t] <= (1 / pars[:COP]) * u_heatpump)

    @constraint(model, sto_soc[1] == u_storage / 2)
    @constraint(model, sto_soc[end] + sto_from_bus[end] * pars[:sto_ef_ch] - sto_to_bus[end] / pars[:sto_ef_dis] == sto_soc[1])

    @constraint(model, heat_sto_soc[1] == u_heat_storage / 2)
    @constraint(model, heat_sto_soc[end] + heat_sto_from_bus[end] * pars[:heat_eff] - heat_sto_to_bus[end] / pars[:heat_eff] == heat_sto_soc[1])

    @constraint(model, pars[:c_pv] * u_pv +
                       pars[:c_wind] * u_wind +
                       pars[:c_storage] * u_storage +
                       pars[:c_heat_storage] * u_heat_storage +
                       pars[:c_heatpump] * u_heatpump
                       <=
                       pars[:inv_budget])

    @constraint(model, [t = 1:number_of_hours],
        gci[t] - gco[t] + u_pv * pv[t] - pv_cur[t] + u_wind * wind[t] - wind_cur[t] - demand[t] + sto_to_bus[t] - sto_from_bus[t] - flow_energy2heat[t] == 0)

    @constraint(model, [t = 1:number_of_hours],
        -heatdemand[t] + heat_sto_to_bus[t] - heat_sto_from_bus[t] + pars[:COP] * flow_energy2heat[t] - pars[:heat_losses] * heat_sto_soc[t] == 0)

    if capacity_flex_constraint
        # For both electical storage and heat storage we require that the separate
        # storage and the spare flow are large enough to satisfy the request.
        # as we can't take the minimum of storage and flow we have to write the
        # constraint four times, for each possible combination of flow and storage
        # being binding:

        # min(s_e, f_e) + min(s_h, f_h) >= F
        # s_e + s_h >= F
        # s_e + f_h >= F
        # f_e + s_h >= F
        # f_e + f_h >= F

        sto_ef_dis = pars[:sto_ef_dis]
        COP = pars[:COP]
        heat_eff = pars[:heat_eff]

        @constraint(model, [t = 1:(number_of_hours-1)],
            pv_cur[t] + wind_cur[t] +
            sto_soc[t+1] * sto_ef_dis +
            flow_energy2heat[t] >= F_pos) # this binds if flow_energy2heat[t] is smaller than heat_sto_soc[t+1]/COP
        @constraint(model, [t = 1:(number_of_hours-1)],
            pv_cur[t] + wind_cur[t] +
            sto_soc[t+1] * sto_ef_dis
            + heat_sto_soc[t+1] * heat_eff / COP >= F_pos) # ... and the other way around. This way we avoid taking the minimum.
        @constraint(model, [t = 1:(number_of_hours-1)],
            pv_cur[t] + wind_cur[t] +
            (pars[:max_sto_flow] * u_storage - (sto_to_bus[t] - sto_from_bus[t])) +
            flow_energy2heat[t] >= F_pos) # this binds if flow_energy2heat[t] is smaller than heat_sto_soc[t+1]/COP
        @constraint(model, [t = 1:(number_of_hours-1)],
            pv_cur[t] + wind_cur[t] +
            (pars[:max_sto_flow] * u_storage - (sto_to_bus[t] - sto_from_bus[t])) +
            heat_sto_soc[t+1] * heat_eff / COP >= F_pos) # ... and the other way around. This way we avoid taking the minimum.

        @constraint(model, [t = 1:(number_of_hours-1)],
            pv_cur[t] - pv[t] * u_pv + wind_cur[t] - wind[t] * u_wind +
            (sto_soc[t+1] - u_storage) / sto_ef_dis +
            (flow_energy2heat[t] - u_heatpump / COP) <= F_neg)
        @constraint(model, [t = 1:(number_of_hours-1)],
            pv_cur[t] - pv[t] * u_pv + wind_cur[t] - wind[t] * u_wind +
            (sto_soc[t+1] - u_storage) / sto_ef_dis +
            (heat_sto_soc[t+1] - u_heat_storage) / (COP * heat_eff) <= F_neg)
        @constraint(model, [t = 1:(number_of_hours-1)],
            pv_cur[t] - pv[t] * u_pv + wind_cur[t] - wind[t] * u_wind -
            (pars[:max_sto_flow] * u_storage - (sto_from_bus[t] - sto_to_bus[t])) +
            (flow_energy2heat[t] - u_heatpump / COP) <= F_neg)
        @constraint(model, [t = 1:(number_of_hours-1)],
            pv_cur[t] - pv[t] * u_pv + wind_cur[t] - wind[t] * u_wind -
            (pars[:max_sto_flow] * u_storage - (sto_from_bus[t] - sto_to_bus[t])) +
            (heat_sto_soc[t+1] - u_heat_storage) / (COP * heat_eff) <= F_neg)
    end


    #──────────────── Second‑stage variables ───────────────#
    if Nscen > 0

        @variables(model, begin
            0 <= pv_cur2[s=1:Nscen, t=1:Trec] <= debug_cap
            0 <= wind_cur2[s=1:Nscen, t=1:Trec] <= debug_cap
            0 <= gci2[s=1:Nscen, t=1:Trec] <= debug_cap
            0 <= gco2[s=1:Nscen, t=1:Trec] <= debug_cap
            0 <= sto_to_bus2[s=1:Nscen, t=1:Trec] <= debug_cap
            0 <= sto_from_bus2[s=1:Nscen, t=1:Trec] <= debug_cap
            0 <= sto_soc2[s=1:Nscen, t=1:Trec] <= debug_cap
            0 <= heat_sto_to_bus2[s=1:Nscen, t=1:Trec] <= debug_cap
            0 <= heat_sto_from_bus2[s=1:Nscen, t=1:Trec] <= debug_cap
            0 <= heat_sto_soc2[s=1:Nscen, t=1:Trec] <= debug_cap
            0 <= flow_energy2heat2[s=1:Nscen, t=1:Trec] <= debug_cap
            # 0 <= ebp[s=1:Nscen, t=1:pars[:recovery_time]] <= debug_cap
            # 0 <= ebm[s=1:Nscen, t=1:pars[:recovery_time]] <= debug_cap
            # 0 <= hbp[s=1:Nscen, t=1:(4+pars[:recovery_time])] <= debug_cap
            # 0 <= hbm[s=1:Nscen, t=1:(4+pars[:recovery_time])] <= debug_cap
            0 <= erp[s=1:Nscen] <= debug_cap
            0 <= erm[s=1:Nscen] <= debug_cap
            0 <= hfp[s=1:Nscen] <= debug_cap
            0 <= hfm[s=1:Nscen] <= debug_cap
        end)

        t_a = []
        tf_a = []
        p_a = []
        for (s, sc) in enumerate(scenarios)
            t_xi, s_xi, F_xi, prob = sc.t_xi, sc.s_xi, sc.F_xi, sc.probability
            t_final = t_xi + pars[:recovery_time]
            push!(t_a, t_xi)
            push!(tf_a, t_final)
            push!(p_a, prob)

            @constraints(model, begin
                [t = 1:Trec], pv_cur2[s, t] <= u_pv * pv[t+t_xi-1]
                [t = 1:Trec], wind_cur2[s, t] <= u_wind * wind[t+t_xi-1]
                pv_cur2[s, Trec] == pv_cur[t_final]
                wind_cur2[s, Trec] == wind_cur[t_final]
                gci2[s, Trec] == gci[t_final]
                gco2[s, Trec] == gco[t_final]
                gci2[s, 1] == gci[t_xi]
                gco2[s, 1] == gco[t_xi]
                [t = 2:Trec], gci2[s, t] <= pars[:feedincap]
                [t = 1:(Trec-1)], sto_soc2[s, t+1] == sto_soc2[s, t] + sto_from_bus2[s, t] * pars[:sto_ef_ch] - sto_to_bus2[s, t] / pars[:sto_ef_dis]
                [t = 1:(Trec-1)], heat_sto_soc2[s, t+1] == heat_sto_soc2[s, t] + heat_sto_from_bus2[s, t] * pars[:heat_eff] - heat_sto_to_bus2[s, t] / pars[:heat_eff]
                [t = 1:Trec], sto_soc2[s, t] <= u_storage
                [t = 1:Trec], heat_sto_soc2[s, t] <= u_heat_storage
                [t = 1:Trec], sto_from_bus2[s, t] <= pars[:max_sto_flow] * u_storage
                [t = 1:Trec], sto_to_bus2[s, t] <= pars[:max_sto_flow] * u_storage
                [t = 1:Trec], flow_energy2heat2[s, t] <= (1 / pars[:COP]) * u_heatpump
                sto_soc2[s, 1] == sto_soc[t_xi]
                sto_soc2[s, Trec] == sto_soc[t_final]
                sto_from_bus2[s, Trec] == sto_from_bus[t_final]
                sto_to_bus2[s, Trec] == sto_to_bus[t_final]
                heat_sto_soc2[s, 1] == heat_sto_soc[t_xi]
                heat_sto_soc2[s, Trec] == heat_sto_soc[t_final] + pars[:COP] * (hfp[s] - hfm[s])
            end)

            @constraint(model, erp[s] - erm[s] + gci2[s, 1] - gco2[s, 1] + u_pv * pv[t_xi] - pv_cur2[s, 1] +
                               u_wind * wind[t_xi] - wind_cur2[s, 1] - demand[t_xi] + sto_to_bus2[s, 1] -
                               sto_from_bus2[s, 1] - flow_energy2heat2[s, 1] - F_xi * s_xi == 0)
            @constraint(model, [τ = 2:Trec], gci2[s, τ] - gco2[s, τ] +
                                             u_pv * pv[t_xi+τ-1] - pv_cur2[s, τ] +
                                             u_wind * wind[t_xi+τ-1] - wind_cur2[s, τ]
                                             -
                                             demand[t_xi+τ-1]
                                             +
                                             sto_to_bus2[s, τ] - sto_from_bus2[s, τ]
                                             -
                                             flow_energy2heat2[s, τ] == 0)
            @constraint(model, [τ = 1:Trec], -heatdemand[t_xi+τ-1] + heat_sto_to_bus2[s, τ] -
                                                 heat_sto_from_bus2[s, τ] + pars[:COP] * flow_energy2heat2[s, τ] - pars[:heat_losses] * heat_sto_soc2[s, τ] == 0)

        end
        @expression(model, scenario_cost[s=1:Nscen],
            pars[:c_i] * (sum(gci2[s, :]) - sum(gci[t_a[s]:tf_a[s]])) -
            pars[:c_o] * (sum(gco2[s, :]) - sum(gco[t_a[s]:tf_a[s]])) +
            pars[:c_i] * hfm[s] - pars[:c_o] * hfp[s] + 0.001 * (hfm[s] + hfp[s]) +
            pars[:penalty] * (erp[s] + erm[s]) -
            reg_lossy_flows * (sum(heat_sto_from_bus[t_a[s]:tf_a[s]]) + sum(heat_sto_to_bus[t_a[s]:tf_a[s]]) + sum(sto_from_bus[t_a[s]:tf_a[s]]) + sum(sto_to_bus[t_a[s]:tf_a[s]])) +
            reg_lossy_flows * (sum(heat_sto_from_bus2[s, :]) + sum(heat_sto_to_bus2[s, :]) + sum(sto_from_bus2[s, :]) + sum(sto_to_bus2[s, :])))
        @expression(model, recourse_expr[s=1:Nscen], p_a[s] * events_per_period * scenario_cost[s])

    end

    @expression(model, investment_cost,
        (u_pv * pars[:c_pv] + u_wind * pars[:c_wind] + u_storage * pars[:c_storage] + u_heat_storage * pars[:c_heat_storage] + u_heatpump * pars[:c_heatpump]) / lifetime_factor)

    @expression(model, base_operation_cost,
        pars[:c_i] * sum(gci) - pars[:c_o] * sum(gco) + reg_lossy_flows * (sum(heat_sto_from_bus) + sum(heat_sto_to_bus) + sum(sto_from_bus) + sum(sto_to_bus)))

    if Nscen > 0
        @objective(model, Min, investment_cost + base_operation_cost + sum(recourse_expr))
    else
        @objective(model, Min, investment_cost + base_operation_cost)
    end

    return model
end

################################################################################
# Helper utilities (investment & operation extraction + fix/unfix)             #
################################################################################

function get_investments(model::Model)
    return NamedTuple((var => value(model[var]) for var in inv_vars))
end

function get_operation(model::Model)
    return (
        pv_cur=value.(model[:pv_cur]),
        wind_cur=value.(model[:wind_cur]),
        gci=value.(model[:gci]),
        gco=value.(model[:gco]),
        sto_to_bus=value.(model[:sto_to_bus]),
        sto_from_bus=value.(model[:sto_from_bus]),
        sto_soc=value.(model[:sto_soc]),
        heat_sto_to_bus=value.(model[:heat_sto_to_bus]),
        heat_sto_from_bus=value.(model[:heat_sto_from_bus]),
        heat_sto_soc=value.(model[:heat_sto_soc]),
        flow_energy2heat=value.(model[:flow_energy2heat]))
end

function fix_investment!(model::Model, investments; verbose=false)
    for (sym, val) in pairs(investments)
        fix(model[sym], val; force=true)
    end
    verbose && println("The investments are fixed.")
    return nothing
end

function unfix_investment!(model::Model)
    error("Unfixing variables does not restore their bounds")
    for sym in inv_vars
        unfix(model[sym])
    end
    println("The investments are released")
    return nothing
end

function fix_operation!(model::Model, operation; verbose=false)
    for (sym, vec) in pairs(operation)
        var_array = model[sym]
        @assert length(var_array) == length(vec) "Dimension mismatch fixing $(sym)"
        for i in eachindex(vec)
            fix(var_array[i], max(0.0, vec[i]); force=true)
        end
    end
    verbose && println("Operational schedule is fixed")
    return nothing
end

function unfix_operation!(model::Model)
    error("Unfixing variables does not restore their bounds")
    op_vars = [:gci, :gco, :sto_soc, :sto_to_bus, :sto_from_bus,
        :heat_sto_soc, :heat_sto_to_bus, :heat_sto_from_bus,
        :flow_energy2heat, :pv_cur, :wind_cur]
    for sym in op_vars
        for v in model[sym]
            unfix(v)
        end
    end
    println("Operational schedule is released")
    return nothing
end


function get_scenario_objectives(sp)
    if termination_status(sp) != MOI.OPTIMAL
        println("Solution not optimal. Returning empty objectives.")
        return (;)
    end

    second_stage_costs = value.(sp[:scenario_cost])

    erp = value.(sp[:erp])
    erm = value.(sp[:erm])
    hfp = value.(sp[:hfp])
    hfm = value.(sp[:hfm])

    return (; erp, erm, hfp, hfm, second_stage_costs)
end


function optimize_resample!(sp_ini, sp_re)
    println("Optimizing sp_ini...")
    optimize!(sp_ini)

    if termination_status(sp_ini) != MOI.OPTIMAL
        return (feasible = false, ini_feasible = false,) , (feasible = false, ini_feasible = false,)
    end

    ini_scens_objs = (feasible = true, ini_feasible = true, get_scenario_objectives(sp_ini)...)

    inv = get_investments(sp_ini)
    op = get_operation(sp_ini)

    fix_investment!(sp_re, inv)
    fix_operation!(sp_re, op)

    println("Optimizing sp_re...")
    optimize!(sp_re)

    if termination_status(sp_ini) != MOI.OPTIMAL
        return ini_scens_objs , (feasible = false, ini_feasible = true,)
    end

    re_scens_objs = (feasible = true, ini_feasible = true, get_scenario_objectives(sp_re)...)

    return ini_scens_objs, re_scens_objs
end

function save_no_scen_results(path, sp::Model, pars;
    name, F=0,)

    if isnothing(path)
        println("skipping save")
        return nothing
    end

    println("Saving the first-stage model $(name) at $(path).")

    status = @show termination_status(sp)

    # Some basic asserts

    if name == "capacity_only"
        @assert F != 0
    end
    
    serialize(joinpath(path, name * "_optimal.serial"), sp)

    # Scalar data

    if status != MOI.OPTIMAL
        serialize(joinpath(path, name * "_all_data.serial"), (
            feasible = false,
            termination_status=string(status),
            name=name,
            path=path,
            F=F, 
            pars...)
        )
    end


    model_costs = (full_cost=objective_value(sp),
        investment_cost=value(sp[:investment_cost]),
        operation_cost=value(sp[:base_operation_cost]))

    inv = get_investments(sp)

    scalar_data = (;
        feasible = true,
        termination_status=string(status),
        name=name,
        path=path,
        F=F, 
        pars...,
        inv..., 
        model_costs...,)

    # Timestep data

    timestep_data = get_operation(sp)

    # All data

    model_data = (; scalar_data..., timestep_data...)

    # To get a 1xN DataFrame, wrap every value of the NamedTuple in an Array
    # CSV.write(path * name * "_scalar_data.csv", DataFrame(map(x -> [x], scalar_data)))

    # CSV.write(path * name * "_timestep_data.csv", DataFrame(timestep_data))
    # CSV does not store the type information, so serialize all data as well
    serialize(joinpath(path, name * "_all_data.serial"), model_data)

    return nothing
end


function save_results(path, sp, pars;
    name,
    scens,
    scens_objectives,
    F=0,
    flex_interval=0,
    n_samples=0,
    sample_data=(;))
    
    serialize(joinpath(path, name * "_optimal.serial"), sp)

    status=string(termination_status(sp))

    if scens_objectives.ini_feasible == false
        status *= " (ini false)"
    end

    scens_data = (;
    t_xi=scens .|> x -> x.t_xi,
    s_xi=scens .|> x -> x.s_xi,
    F_xi=scens .|> x -> x.F_xi)

    always_data = (;
    scens_objectives..., 
    name=name,
    path=path,
    modeltypename=typeof(sp).name.name,
    F=F,
    flex_interval=flex_interval,
    n_samples=n_samples,
    termination_status = status,
    pars..., 
    sample_data..., 
    scens_data...
    )

    if scens_objectives.feasible == false
        serialize(joinpath(path, name * "_all_data.serial"),  always_data)
        return nothing
    end

    model_costs = (full_cost=objective_value(sp),
        investment_cost=value(sp[:investment_cost]),
        operation_cost=value(sp[:base_operation_cost]))

    inv = get_investments(sp)

    scalar_data = (; inv..., model_costs...,)

    # Timestep data

    timestep_data = get_operation(sp)

    # All data

    model_data = (; always_data..., scalar_data..., timestep_data...)

    # To get a 1xN DataFrame, wrap every value of the NamedTuple in an Array
    # CSV.write(path * name * "_scalar_data.csv", DataFrame(map(x -> [x], scalar_data)))
    # CSV.write(path * name * "_scenario_data.csv", DataFrame(scens_data))
    # CSV.write(path * name * "_timestep_data.csv", DataFrame(timestep_data))
    # CSV does not store the type information, so serialize all data as well
    serialize(joinpath(path, name * "_all_data.serial"), model_data)

    return nothing
end
