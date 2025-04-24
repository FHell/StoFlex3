
df_file = "preloaded/df_raw.serial"

# point the variable df_file to the raw aggregated experimental
# results in order to skip loading from the remote.
df = nothing
include(joinpath(@__DIR__, "load_data_preamble.jl"))

##

using LinearAlgebra
using CairoMakie
using LaTeXStrings

include("../src/load_energy_system_data.jl")
include("../src/prepare_stochastic_programs.jl")
include("../src/sp_model.jl")
include("../src/new_eval_utils.jl")

data_es = load_industry_park1()

##


function plot_split_flows(df, row; in_var=:gci, out_var=:gco)

    fig = Figure(;
    figure_padding = (5, 5, 10, 10),
    backgroundcolor = :snow2,
    size = (800, 800),)

    ax1 = Axis(fig[1, 1])

    lines!(ax1, df[row, in_var])
    lines!(ax1, -1 .* df[row, out_var])

    ax2 = Axis(fig[2, 1])

    lines!(ax2, df[row, in_var] .- df[row, out_var])

    fig
end

function plot_raw_results(df, row;
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

function plot_recovery_window_deviation_raw(sp; scens, recovery_time, s = 1, vars = ["flow_energy2heat", "gci", "gco", "sto_soc", "sto_to_bus", "sto_from_bus", "heat_sto_soc", "pv_cur", "wind_cur"])

    t_xi = scens.t_xi[s]
    s_xi = scens.s_xi[s]
    F_xi = scens.F_xi[s]
    println("time=$(t_xi), sign=$s_xi, size=$F_xi")


    recovery_window = (t_xi):(t_xi + recovery_time)

    fig = Figure(;
    figure_padding = (5, 5, 10, 10),
    backgroundcolor = :snow2,
    size = (1200, 400 * length(vars)),)

    for (i, var) in enumerate(vars)
        ax = Axis(fig[i, 1];
        ylabel = var)

        lines!(ax, recovery_window,
            value.(sp[Symbol(var)])[recovery_window],
            label = "1st stage", color = :black)
        lines!(ax, recovery_window,
            value.(sp[Symbol(var * "2")][s, :]),
            label = "2nd stage", color = :red)
    end

    fig
end

function app2(sym)
    Symbol(String(sym)*"2")
end

##

model_idx = 450

plot_split_flows(df, model_idx)

##

plot_raw_results(df, model_idx)

##

sp, scens, recovery_time = load_model(df, model_idx)

##

optimize!(sp)
@show termination_status(sp)

##

plot_recovery_window_deviation_raw(sp; scens, recovery_time, s = 1)

