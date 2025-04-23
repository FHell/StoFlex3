
df_file = "preloaded/df_raw.serial"

# point the variable df_file to the raw aggregated experimental
# results in order to skip loading from the remote.

include(joinpath(@__DIR__, "load_data_preamble.jl"))

##

using LinearAlgebra

include("../src/load_energy_system_data.jl")
include("../src/prepare_stochastic_programs.jl")
include("../src/StochasticFlexibility.jl")

data_es = load_max_boegl()

##

experiment = "no_solar"
experiment = "large_heat"


if length(ARGS) == 0
    experiment = "defaults"
else
    experiment = ARGS[1]
end


if ! @isdefined dont_display
    dont_display = false
end

##

df_exp = filter_kwd(df; first_dir = experiment)


CB = get_CB(df)

flex_interval = 4
n_samples = 4
df = filter_kwd(df_exp; flex_interval, name = "OFIOR", n_samples, F=500)


##

model_idx = 1

##

# Find scenarios that are excessively expensive:
expensive_idx = findall(df.total_costs[model_idx] .> CB * 1.5)

##

# Requires that the remote is mounted, and probably also that the data was aggregated
# from the same place at which it is now being analyzed.

sp_optimal = load_sp(joinpath(df.path[model_idx], "$(df.name[model_idx])_optimal.serial"))

# optimize!(sp_optimal)

##

sp_data = deserialize(joinpath(df.path[model_idx], "$(df.name[model_idx])_all_data.serial"))


##

# No time occurs twice:
length(df.t_xi[model_idx]) == length(df.t_xi[model_idx] |> unique)

problem_times = df.t_xi[model_idx][expensive_idx]

##

using CairoMakie
using LaTeXStrings

##

function df_results_figure(df, row;
    plot_window = nothing,
    vars = ["gci", "gco", "sto_soc", "heat_sto_soc", "flow_energy2heat"])
    if isnothing(plot_window)
        plot_window = 1:length(df[row, "gci"])
    end

    fig = Figure(;
    figure_padding = (5, 5, 10, 10),
    backgroundcolor = :snow2,
    size = (800, 200*length(vars)),)


    for (i, var) in enumerate(vars)

        ax = Axis(fig[i, 1];
        ylabel = var)

        lines!(ax, plot_window, df[row, var][plot_window])
    end

    fig
end

##

plot_window = (-10:10) .+ problem_times[1]

fig = df_results_figure(df, 1; plot_window, vars=["flow_energy2heat", "heat_sto_from_bus", "heat_sto_to_bus", "sto_from_bus", "sto_to_bus"])

##

function plot_recovery_window_deviation(sp; s = 1, vars = ["flow_energy2heat", "gci", "gco", "sto_soc"])
    fig = Figure(;
    figure_padding = (5, 5, 10, 10),
    backgroundcolor = :snow2,
    size = (800, 200 * length(vars)),)

    t_xi = scenarios(sp)[s].data.t_xi
    println(t_xi)
    recovery_time = sp.stages[2].parameters[:recovery_time]

    recovery_window = (t_xi):(t_xi + recovery_time)

    for (i, var) in enumerate(vars)
        ax = Axis(fig[i, 1];
        ylabel = var)

        lines!(ax, recovery_window,
            value.(sp[1, Symbol(var)])[recovery_window],
            label = "1st stage", color = :black)
        lines!(ax, recovery_window,
            value.(sp[2, Symbol(var * "2")], s),
            label = "2nd stage", color = :red)
    end

    fig
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

function plot_recovery_window_deviation_resim(sp; s = 1, om = get_optimized_outcome_model(sp, s), vars = ["flow_energy2heat", "gci", "gco", "sto_soc", "sto_to_bus", "sto_from_bus", "heat_sto_soc", "pv_cur", "wind_cur"])

    t_xi = scenarios(sp)[s].data.t_xi
    s_xi = scenarios(sp)[s].data.s_xi
    F_xi = scenarios(sp)[s].data.F_xi
    println("time=$(t_xi), sign=$s_xi, size=$F_xi")

    recovery_time = sp.stages[2].parameters[:recovery_time]

    recovery_window = (t_xi):(t_xi + recovery_time)

    fig = Figure(;
    figure_padding = (5, 5, 10, 10),
    backgroundcolor = :snow2,
    size = (800, 200 * length(vars)),)

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



##

fig = plot_recovery_window_deviation_resim(sp_optimal; s=1)

##

exp_oms = [(s=s, om=get_optimized_outcome_model(sp_optimal, s)) for s in expensive_idx]

##

fig = plot_recovery_window_deviation_resim(sp_optimal; exp_oms[1]...)

save("penalty.png", fig)
