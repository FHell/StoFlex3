
df_file = "preloaded/df_raw.serial"

# point the variable df_file to the raw aggregated experimental
# results in order to skip loading from the remote.

include(joinpath(@__DIR__, "load_data_preamble.jl"))

##

using LinearAlgebra

include("../src/load_energy_system_data.jl")
include("../src/prepare_stochastic_programs.jl")
include("../src/sp_model.jl")

data_es = load_max_boegl()

##

# experiment = "no_solar"
# experiment = "large_heat"

if ! @isdefined experiment
    experiment = "defaults"
end

if ! @isdefined dont_display
    dont_display = false
end


##
using CairoMakie
using LaTeXStrings

figure_dir = joinpath(@__DIR__, "../plots/", experiment)
mkpath(figure_dir)


##

df_exp = filter_kwd(df; first_dir = experiment)


CB = get_CB(df)

flex_interval = 12
n_samples = 6
F=500
df = filter_kwd(df_exp; flex_interval, name = "OFIOR", n_samples, F)


##

model_idx = 1

##

if isfile("preloaded/sp_rec_schedule_$experiment.serial")
    println("Loading preloaded stochastic model")
    sp_optimal = deserialize("preloaded/sp_rec_schedule_$experiment.serial")
else
    println("Loading stochastic model from HPC ../remote_results/")
    sp_optimal = load_sp(joinpath(df.path[model_idx], "$(df.name[model_idx])_optimal.serial"))
    serialize("preloaded/sp_rec_schedule_$experiment.serial", sp_optimal)
end

##


function get_optimized_outcome_model(sp, s)
    t_xi = scenarios(sp)[s].data.t_xi
    s_xi = scenarios(sp)[s].data.s_xi
    F_xi = scenarios(sp)[s].data.F_xi

    scen = @scenario t_xi=t_xi s_xi=s_xi F_xi=F_xi probability=1.0
    om = outcome_model(sp,
        optimal_decision(sp),
        scen;
        optimizer = subproblem_optimizer(sp))
    optimize!(om)
    om
end

function plot_recovery_window_deviation_raw(sp; s = 1, om = get_optimized_outcome_model(sp, s), vars = ["flow_energy2heat", "gci", "gco", "sto_soc", "sto_to_bus", "sto_from_bus", "heat_sto_soc", "pv_cur", "wind_cur"])

    t_xi = scenarios(sp)[s].data.t_xi
    s_xi = scenarios(sp)[s].data.s_xi
    F_xi = scenarios(sp)[s].data.F_xi
    println("time=$(t_xi), sign=$s_xi, size=$F_xi")

    recovery_time = sp.stages[2].parameters[:recovery_time]

    recovery_window = (t_xi):(t_xi + recovery_time)

    fig = Figure(;
    figure_padding = (5, 5, 10, 10),
    backgroundcolor = :snow2,
    size = (1200, 400 * length(vars)),)

    for (i, var) in enumerate(vars)
        ax = Axis(fig[i, 1];
        ylabel = var)

        lines!(ax, recovery_window,
            value.(sp[1, Symbol(var)])[recovery_window],
            label = "1st stage", color = :black)
        lines!(ax, recovery_window,
            value.(om[Symbol(var * "2")]),
            label = "2nd stage", color = :red)
    end

    fig
end

function app2(sym)
    Symbol(String(sym)*"2")
end

##

function plot_recovery_agg(sp, s, om)
    update_theme!(fontsize =30)

    df_1 = DataFrame()
    df_2 = DataFrame()

    t_xi = scenarios(sp)[s].data.t_xi
    s_xi = scenarios(sp)[s].data.s_xi
    F_xi = scenarios(sp)[s].data.F_xi
    println("time=$(t_xi), sign=$s_xi, size=$F_xi")

    recovery_time = sp.stages[2].parameters[:recovery_time]

    recovery_window = (t_xi):(t_xi + recovery_time)

    df_1.t = recovery_window
    df_2.t = recovery_window

    var_in = [:gci, :sto_from_bus, :sto_soc, :flow_energy2heat,]
    var_out = [:gco, :sto_to_bus, nothing, nothing]

    recs = ["Grid Connection\nkW",
    "Battery Charging\nkW",
    "Battery SOC\nkWh",
    "Heatpump\nkW"]

    for (vi, vo, vr) in zip(var_in, var_out, recs)
        df_1[!, vr] = value.(sp[1, vi])[recovery_window]
        !isnothing(vo) && (df_1[!, vr] .-= value.(sp[1, vo])[recovery_window])

        df_2[!, vr] = value.(om[app2(vi)])
        !isnothing(vo) && (df_2[!, vr] .-= value.(om[app2(vo)]))
    end

    fig = Figure(;
    figure_padding = (5, 5, 10, 10),
    size = (800, 300 * length(recs)),)

    for (i, rec) in enumerate(recs)
        ax = Axis(fig[i, 1];
        ylabel = rec)

        lines!(ax, df_1.t,
            df_1[!, rec],
            label = "1st Stage", color = :black, linewidth=5)
        lines!(ax, df_2.t,
            df_2[!, rec],
            label = "2nd Stage", color = :red, linestyle=:dash, linewidth=5)
    end
    fig.content 

    reset_limits!.(fig.content)
    yspace = maximum(tight_yticklabel_spacing!, fig.content)
    fig.content .|> x -> x.yticklabelspace = yspace
    fig.content[end].xlabel = "Time (h)"

    fig
end


##

# Find scenarios that are excessively expensive:

##
mkpath(joinpath(figure_dir,"publication"))

##
s = 10

om = get_optimized_outcome_model(sp_optimal, s)

fig = plot_recovery_agg(sp_optimal, s, om)

save(joinpath(figure_dir, "publication", "recovery_schedule_$s.png"), fig)

##
s = 11

om = get_optimized_outcome_model(sp_optimal, s)

fig = plot_recovery_agg(sp_optimal, s, om)

save(joinpath(figure_dir, "publication", "recovery_schedule_$s.png"), fig)

##
s = 12

om = get_optimized_outcome_model(sp_optimal, s)

fig = plot_recovery_agg(sp_optimal, s, om)

save(joinpath(figure_dir, "publication", "recovery_schedule_$s.png"), fig)

##
s = 13

om = get_optimized_outcome_model(sp_optimal, s)

fig = plot_recovery_agg(sp_optimal, s, om)

save(joinpath(figure_dir, "publication", "recovery_schedule_$s.png"), fig)

##
s = 14

om = get_optimized_outcome_model(sp_optimal, s)

fig = plot_recovery_agg(sp_optimal, s, om)

save(joinpath(figure_dir, "publication", "recovery_schedule_$s.png"), fig)

##


# for s in 10:14
#     om = get_optimized_outcome_model(sp_optimal, s)

#     fig = plot_recovery_agg(sp_optimal, s, om)
    
#     save(joinpath(figure_dir, "publication", "recovery_schedule_$s.png"), fig)
# end


##

# fig = plot_recovery_window_deviation_resim(sp_optimal; s, om)

##

expensive_idx = findmax(df.total_costs[model_idx])[2]

s = expensive_idx

om = get_optimized_outcome_model(sp_optimal, s)

fig = plot_recovery_agg(sp_optimal, s, om)

##

save(joinpath(figure_dir, "recovery_schedule_$s.png"), fig)

##

# We find the scenario where the optimizer could get the cost as low as possible.

cheap_idx = findmin(df.total_costs[model_idx])[2]

s = cheap_idx

om = get_optimized_outcome_model(sp_optimal, s)

fig = plot_recovery_agg(sp_optimal, s, om)

##

save(joinpath(figure_dir, "recovery_schedule_cheap.png"), fig)
