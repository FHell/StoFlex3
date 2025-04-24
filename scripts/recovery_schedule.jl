
df_file = "preloaded/df_raw.serial"

# point the variable df_file to the raw aggregated experimental
# results in order to skip loading from the remote.

include(joinpath(@__DIR__, "load_data_preamble.jl"))

##

using LinearAlgebra

include("../src/load_energy_system_data.jl")
include("../src/prepare_stochastic_programs.jl")
include("../src/sp_model.jl")
include("../src/new_eval_utils.jl")

data_es = load_industry_park1()

##

# experiment = "no_solar"
# experiment = "large_heat"

if !@isdefined experiment
    experiment = "defaults"
end

if !@isdefined dont_display
    dont_display = false
end


##
using CairoMakie
using LaTeXStrings

figure_dir = joinpath(@__DIR__, "../plots/", experiment)
mkpath(figure_dir)


##

df_exp = filter_kwd(df; first_dir=experiment)


CB = get_CB(df)

flex_interval = 12
n_samples = 12
F = 500
df = filter_kwd(df_exp; flex_interval, name="OFIOR", n_samples, F)


##

model_idx = 1

##

sp, scens, recovery_time = load_model(df, 1)

##

optimize!(sp)

##

function app2(sym)
    Symbol(String(sym) * "2")
end


function plot_recovery_agg(sp; scens, recovery_time, s=1)
    update_theme!(fontsize=30)

    df_1 = DataFrame()
    df_2 = DataFrame()

    t_xi = scens.t_xi[s]
    s_xi = scens.s_xi[s]
    F_xi = scens.F_xi[s]
    println("time=$(t_xi), sign=$s_xi, size=$F_xi")

    recovery_window = (t_xi):(t_xi+recovery_time)

    df_1.t = recovery_window
    df_2.t = recovery_window

    var_in = [:gci, :sto_from_bus, :sto_soc, :flow_energy2heat,]
    var_out = [:gco, :sto_to_bus, nothing, nothing]

    recs = ["Grid Connection\nkW",
        "Battery Charging\nkW",
        "Battery SOC\nkWh",
        "Heatpump\nkW"]

    for (vi, vo, vr) in zip(var_in, var_out, recs)
        df_1[!, vr] = value.(sp[vi][recovery_window])
        !isnothing(vo) && (df_1[!, vr] .-= value.(sp[vo][recovery_window]))

        df_2[!, vr] = value.(sp[app2(vi)][s, :])
        !isnothing(vo) && (df_2[!, vr] .-= value.(sp[app2(vo)][s, :]))
    end

    fig = Figure(;
        figure_padding=(5, 5, 10, 10),
        size=(800, 300 * length(recs)),)

    for (i, rec) in enumerate(recs)
        ax = Axis(fig[i, 1];
            ylabel=rec)

        lines!(ax, df_1.t,
            df_1[!, rec],
            label="1st Stage", color=:black, linewidth=5)
        lines!(ax, df_2.t,
            df_2[!, rec],
            label="2nd Stage", color=:red, linestyle=:dash, linewidth=5)
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
mkpath(joinpath(figure_dir, "publication"))

##
for s in [10:14..., 100:140...]

    fig = plot_recovery_agg(sp; scens, recovery_time, s)

    save(joinpath(figure_dir, "publication", "recovery_schedule_$s.png"), fig)

end


##

_, s = findmax(df.second_stage_costs[model_idx])
fig = plot_recovery_agg(sp; scens, recovery_time, s)
save(joinpath(figure_dir, "publication", "recovery_schedule_$(s)_expensive.png"), fig)

_, s = findmin(df.second_stage_costs[model_idx])
fig = plot_recovery_agg(sp; scens, recovery_time, s)
save(joinpath(figure_dir, "publication", "recovery_schedule_$(s)_cheap.png"), fig)

##
