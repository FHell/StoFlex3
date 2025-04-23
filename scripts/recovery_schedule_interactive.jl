
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

data_es = load_max_boegl()

##


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

function plot_recovery_agg(sp, s, om;
    var_in = [:gci, :sto_from_bus, :sto_soc, :flow_energy2heat, :heat_sto_soc],
    var_out = [:gco, :sto_to_bus, nothing, nothing, nothing],
    recs = ["Grid Connection\nkW",
    "Battery Charging\nkW",
    "Battery SOC\nkWh",
    "Heatpump\nkW",
    "Heat SOC\nkWh",
    ],
    )
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
#    fig.content 

    reset_limits!.(fig.content)
    yspace = maximum(tight_yticklabel_spacing!, fig.content)
    fig.content .|> x -> x.yticklabelspace = yspace
    fig.content[end].xlabel = "Time (h)"

    fig, df_1, df_2
end


##

model_idx = 450

recovery_time = df.recovery_time[model_idx]

sp = deserialize(joinpath(df.path[model_idx], "$(df.name[model_idx])_optimal.serial"))
scens = (t_xi = df.t_xi[model_idx], s_xi = df.s_xi[model_idx], F_xi = df.F_xi[model_idx])
set_optimizer(sp, Gurobi.Optimizer)

##

optimize!(sp)
@show termination_status(sp)

##

plot_recovery_window_deviation_raw(sp; scens, recovery_time, s = 1)


##
findfirst(df.infeasible)

##



##

hbps = [sum(value.(sp_ini[2, :hbp], s)) for s in 1:length(scenarios(sp))]
hbms = [sum(value.(sp_ini[2, :hbm], s)) for s in 1:length(scenarios(sp))]
ebps = [sum(value.(sp_ini[2, :ebp], s)) for s in 1:length(scenarios(sp))]
ebms = [sum(value.(sp_ini[2, :ebm], s)) for s in 1:length(scenarios(sp))]

hfps = [value(sp_ini[2, :hfp], s) for s in 1:length(scenarios(sp))]
hfms = [value(sp_ini[2, :hfm], s) for s in 1:length(scenarios(sp))]
erps = [value(sp_ini[2, :erp], s) for s in 1:length(scenarios(sp))]
erms = [value(sp_ini[2, :erm], s) for s in 1:length(scenarios(sp))]

##

using Serialization

for s in 1:length(scenarios(sp))
    println("$s of $(length(scenarios(sp)))")
    om = get_optimized_outcome_model(sp, 1)
    if OPTIMAL != termination_status(om)
        serialize("om_$s.serial", om)
    end

end


##

for s in 1:length(scenarios(sp))
    println("$s of $(length(scenarios(sp)))")
    global om = get_optimized_outcome_model(sp, 1)
    if OPTIMAL != termination_status(om)
        break
    end
end

##

optimize!(sp)
@show typeof(sp)
@show termination_status(sp)

dp = DEP(sp) # Deterministic Equivalent Model
optimize!(dp) # JuMP optimization
@show typeof(dp)
@show termination_status(dp)

##

conflict = compute_conflict!(dp)

##

s = t.hfp[model_idx] |> argmax

experiment = t.first_dir[model_idx]


t.t_xi[model_idx][s]

##

a_hfp  = zeros(length(df.name))
a_hsoc = zeros(length(df.name))

for i = 1:length(df.name)
    println(i)
    sp = load_sp(joinpath(df.path[i], "$(df.name[i])_optimal.serial"))
    hfp = 0.
    hsoc = 0.
    if ! ismissing(df.n_samples[i])
        for s in 1:length(scenarios(sp))
            hfp += sum(value.(sp[2, :hfp], s)[1:3])
            hsoc += sum(value.(sp[2, :hfp], s)[4])
        end
    end
    a_hfp[i] = hfp
    a_hsoc[i] = hsoc
end
    
##

om = get_optimized_outcome_model(sp, s)

##

println(t.F_xi[model_idx][s])
println(t.t_xi[model_idx][s])
println(t.s_xi[model_idx][s])
println(t.ebp[model_idx][s])
println(t.ebm[model_idx][s])
println(t.erp[model_idx][s])
println(t.erm[model_idx][s])
println(t.hfp[model_idx][s])
println(t.hfm[model_idx][s])

##

fig, df_1, df_2 = plot_recovery_agg(sp, s, om)
fig

##

hfp     = 0.
hsoc = 0.
for s in 1:length(scenarios(sp))
    hfp += sum(value.(sp[2, :hfp], s)[1:3])
    hsoc += sum(value.(sp[2, :hfp], s)[4])
end
hfp, hsoc

##

set_optimizer(sp, Gurobi.Optimizer)
optimize!(sp)
@show typeof(sp)
@show termination_status(sp)

dp = DEP(sp) # Deterministic Equivalent Model
optimize!(dp) # JuMP optimization
@show typeof(dp)
@show termination_status(dp)

set_optimizer(dp, Gurobi.Optimizer)
optimize!(dp) # JuMP optimization
@show typeof(dp)
@show termination_status(dp)

##

for s in 1:7504
    println(s)
    global om = get_optimized_outcome_model(sp_optimal, s)
    if termination_status(om) != OPTIMAL
        break
    end
end

##

dp = DEP(sp_optimal)
optimize!(dp)

##

using Gurobi
set_optimizer(dp, Gurobi.Optimizer)

##


##

m_ebp = df.ebp[model_idx]

sidx = sortperm(m_ebp) |> reverse

##
s = sidx[1]
om = get_optimized_outcome_model(sp_optimal, s)

##

println(df.F_xi[model_index][s])
println(df.t_xi[model_index][s])
println(df.s_xi[model_index][s])
println(df.ebp[model_index][s])
println(df.ebm[model_index][s])
println(df.hbp[model_index][s])
println(df.hbm[model_index][s])

##

fig = plot_recovery_agg(sp_optimal, 1, om)
fig

##

save("interactive.png", fig)
