df_file = "preloaded/df_raw.serial"

# point the variable df_file to the raw aggregated experimental
# results in order to skip loading from the remote.

include(joinpath(@__DIR__, "load_data_preamble.jl"))

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

df = filter_kwd(df; first_dir = experiment)

CB = get_CB(df)

df.cost_above_baseline = df.no_penalty_cost .- CB
df.cab_pkWh = df.cost_above_baseline ./ df.expected_flex_energy_per_full_period
df.cab_pMWh = df.cost_above_baseline ./ df.expected_flex_energy_per_full_period .* 1e3

sort!(df, :F)

##

using CairoMakie
using LaTeXStrings

figure_dir = joinpath(@__DIR__, "../plots", experiment)
mkpath(figure_dir)

##

df_rF = filter_kwd(df, n_samples = 12, name = "OFIOR")

hm_vars = [:cab_pMWh :operation_cost :investment_cost;
    :u_storage :u_heat_storage :u_heatpump]
hm_factor = [1.0 1e-6 1e-6; 1e-4 1e-4 1e-4]
hm_labels = ["Cost above Baseline (€ per MWh)" "Operation Cost (M€)" "Total investment (M€)";
    "Elec. Storage (MWh)" "Heat Storage (MWh)" "Heat Pump (MW)"]

F_idx = Dict(n => i for (i, n) in enumerate(F_array))
int_idx = Dict(n => i for (i, n) in enumerate(flex_interval_array))

f_loc = map(x -> F_idx[x], df_rF.F)
i_loc = map(x -> int_idx[x], df_rF.flex_interval);

##
update_theme!(fontsize = 30)

fig = Figure(size = (800 .* size(hm_vars)[2], 500 .* size(hm_vars)[1]))

for i in 1:size(hm_vars)[1]
    for j in 1:size(hm_vars)[2]
        ax = Axis(fig[i, 2 * j - 1],
            xticklabelsvisible = (i == size(hm_vars)[1]),
            #            xticklabelrotation=-π / 2,
            yticklabelsvisible = (j == 1),
            xlabel = L"F_G,\, (kW)",
            xlabelvisible = (i == size(hm_vars)[1]),
            title = hm_labels[i, j],
            xticks = (1:5, string.(F_array)),
            yticks = (1:3, string.(flex_interval_array)))
        hm = heatmap!(ax,
            f_loc,
            i_loc,
            df_rF[!, hm_vars[i, j]] .* hm_factor[i, j],
            colormap = :haline)
        Colorbar(fig[i, 2 * j], hm)
    end
end

Label(fig[:, 0], "Average interval between flexibility requests (h)", rotation = π / 2)

save(joinpath(figure_dir, "rate_F_plot.png"), fig)

fig

##
