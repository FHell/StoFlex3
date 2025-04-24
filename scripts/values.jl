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

# select!(df,
#     [:name,
#         :F,
#         :flex_interval,
#         :n_samples,
#         :cab_pkWh,
#         :cost_above_baseline,
#         :total_penalty,
#         :penalty_cost,
#         :optimized,
#         :expected_flex_energy_per_full_period,
#         :F_idx,
#         :operation_cost
#     ])

sort!(df, :F)

##


function values_figures(dfs_cab, dfs_MWph)
    tick_strings = string.(F_array)
    tick_strings[2] = ""

    fig = Figure(;
        figure_padding = (5, 5, 10, 10),
        backgroundcolor = :snow2,
        size = (800, 800),)

    ax1 = Axis(fig[1, 1];
        xlabel = "Maximum size of Flexibility Requests",
        xticks = F_array, # (1:5, string.(F_array)),
        xlabelvisible = false,
        xticklabelsvisible = false,
        xticksvisible = false,
        ylabel = "Value of stochastic optimization (€)")

    for df_m in dfs_cab
        scatterlines!(ax1,
            df_m.F,
            df_m.value,
            label = "$(df_m.plot_name[1])")
    end

    ax2 = Axis(fig[2, 1];
        xlabel = "Maximum size of Flexibility Requests",
        xticks = (F_array, tick_strings), # (1:5, string.(F_array)),
        ylabel = "Value of stochastic optimization per MWh (€)")

    for df_m in dfs_MWph
        scatterlines!(ax2,
            df_m.F,
            df_m.value_pMWh,
            label = "$(df_m.plot_name[1])")
    end

    Legend(fig[:, 2], ax2)
    fig
end

##
update_theme!(fontsize = 20)

for flex_interval in flex_interval_array
    n_samples = 12

    local df_OFIOR = filter_kwd(df; flex_interval, name = "OFIOR", n_samples)

    local df_OFOR = filter_kwd(df; flex_interval, name = "OFOR", n_samples)

    local df_OFR = filter_kwd(df; flex_interval, name = "OFR", n_samples)

    df_val_sto_operation = DataFrame()
    df_val_sto_operation.value = df_OFR.cost_above_baseline .- df_OFOR.cost_above_baseline
    df_val_sto_operation.value_pMWh = df_OFR.cab_pMWh .- df_OFOR.cab_pMWh
    df_val_sto_operation.plot_name .= "for Operations" # L"$C(S|I_F, O_F) - C(S|I_F)$"

    df_val_sto_invest = DataFrame()
    df_val_sto_invest.value = df_OFOR.cost_above_baseline .- df_OFIOR.cost_above_baseline
    df_val_sto_invest.value_pMWh = df_OFOR.cab_pMWh .- df_OFIOR.cab_pMWh
    df_val_sto_invest.plot_name .= "for Investments" # L"$C(S|I_F) - C(S)$"

    df_val_sto = DataFrame()
    df_val_sto.value = df_OFR.cost_above_baseline .- df_OFIOR.cost_above_baseline
    df_val_sto.value_pMWh = df_OFR.cab_pMWh .- df_OFIOR.cab_pMWh
    df_val_sto.plot_name .= "Total" # L"$C(S|I_F, O_F) - C(S)$"

    dfs_val = [df_val_sto, df_val_sto_invest, df_val_sto_operation]

    for df in dfs_val
        df.F .= df_OFIOR.F
    end

    local fig = values_figures(dfs_val,
    dfs_val)

    save(joinpath(figure_dir, "stochastic_optimization_added_value_interv-$flex_interval.png"), fig)

end

##
flex_interval = 6
n_samples = 12

df_OFIOR = filter_kwd(df; flex_interval, name = "OFIOR", n_samples)

df_OFOR = filter_kwd(df; flex_interval, name = "OFOR", n_samples)

df_OFR = filter_kwd(df; flex_interval, name = "OFR", n_samples)

df_val_sto_operation = DataFrame()
df_val_sto_operation.value = df_OFR.cost_above_baseline .- df_OFOR.cost_above_baseline
df_val_sto_operation.value_pMWh = df_OFR.cab_pMWh .- df_OFOR.cab_pMWh
df_val_sto_operation.plot_name .= "for Operations" # L"$C(S|I_F, O_F) - C(S|I_F)$"

df_val_sto_invest = DataFrame()
df_val_sto_invest.value = df_OFOR.cost_above_baseline .- df_OFIOR.cost_above_baseline
df_val_sto_invest.value_pMWh = df_OFOR.cab_pMWh .- df_OFIOR.cab_pMWh
df_val_sto_invest.plot_name .= "for Investments" # L"$C(S|I_F) - C(S)$"

df_val_sto = DataFrame()
df_val_sto.value = df_OFR.cost_above_baseline .- df_OFIOR.cost_above_baseline
df_val_sto.value_pMWh = df_OFR.cab_pMWh .- df_OFIOR.cab_pMWh
df_val_sto.plot_name .= "Total" # L"$C(S|I_F, O_F) - C(S)$"

dfs_val = [df_val_sto, df_val_sto_invest, df_val_sto_operation]

for df in dfs_val
    df.F .= df_OFIOR.F
end


##

function values_publication(dfs_MWph)
    update_theme!(fontsize = 35)

    tick_strings = string.(F_array)
    tick_strings[2] = ""

    fig = Figure(;
        figure_padding = (5, 5, 10, 10),
        size = (1000, 600),
        )

    ax2 = Axis(fig[1, 1];
        xlabel = "Maximum size of Flexibility Requests",
        xticks = (F_array, tick_strings), # (1:5, string.(F_array)),
        ylabel = "Stochastic Value per MWh (€)")

    for df_m in dfs_MWph
        scatterlines!(ax2,
            df_m.F,
            df_m.value_pMWh,
            label = "$(df_m.plot_name[1])",
            markersize = 20,
            )
    end

    Legend(fig[1, 2], ax2)
    fig
end


fig = values_publication(dfs_val)

##

mkpath(joinpath(figure_dir, "publication"))

save(joinpath(figure_dir, "publication", "stochastic_optimization_added_value.png"), fig)
