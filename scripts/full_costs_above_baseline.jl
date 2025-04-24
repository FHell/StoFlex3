df_file = "preloaded/df_raw.serial"

# point the variable df_file to the raw aggregated experimental
# results in order to skip loading from the remote.


include(joinpath(@__DIR__, "load_data_preamble.jl"))


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

println("Experiment set to $(experiment)")

df = filter_kwd(df; first_dir=experiment)

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

# df.pathological = (df.gci .|> maximum) .> 5e6
# count(df.pathological)

# pathological indicates that the optimizer is abusing the discretization of the stochastic
# program. These experiments should not be evaluated.



##

using CairoMakie
using LaTeXStrings

figure_dir = joinpath(@__DIR__, "../plots", experiment)
mkpath(figure_dir)

println("Saving to $(figure_dir)")

##

function cost_above_baseline_fig(dfs_cab, dfs_MWph)
    # update_theme!(fontsize = 30)

    tick_strings = string.(F_array)
    tick_strings[2] = ""

    fig = Figure(;
        figure_padding=(5, 5, 10, 10),
        backgroundcolor=:snow2,
        size=(800, 1200),)

    ax1 = Axis(fig[1, 1];
        xlabel="Maximum size of Flexibility Requests",
        xticks=F_array, # (1:5, string.(F_array)),
        xlabelvisible=false,
        xticklabelsvisible=false,
        xticksvisible=false,
        ylabel="Cost above baseline (k€)")

    for df_m in dfs_cab
        scatterlines!(ax1,
            df_m.F,
            df_m.cost_above_baseline ./ 1e3,
            color=m_colors[df_m.name[1]],
            label="$(m_names[df_m.name[1]])")
    end

    ax2 = Axis(fig[2, 1];
        xlabel="Maximum size of Flexibility Requests",
        xticks=(F_array, tick_strings), # (1:5, string.(F_array)),
        ylabel="Cost above baseline per MWh (€)")

    for df_m in dfs_MWph
        scatterlines!(ax2,
            df_m.F,
            df_m.cab_pMWh,
            color=m_colors[df_m.name[1]],
            label="$(m_names[df_m.name[1]])")
    end

    ax3 = Axis(fig[3, 1];
        xlabel="Maximum size of Flexibility Requests",
        xticks=(F_array, tick_strings), # (1:5, string.(F_array)),
        ylabel="Requested Energy not served (kWh)",
    )

    tp = 0.0

    for df_m in dfs_MWph
        tp += sum(df_m.energy_not_served)
        scatterlines!(ax3,
            df_m.F,
            df_m.energy_not_served ./ df_m.expected_flex_energy_per_full_period,
            color=m_colors[df_m.name[1]],
            label="$(m_names[df_m.name[1]])")
    end

    if tp < 1
        ylims!(ax3, -1.0, 1.0)
    end

    ax4 = Axis(fig[4, 1];
        xlabel="Maximum size of Flexibility Requests",
        xticks=(F_array, tick_strings), # (1:5, string.(F_array)),
        ylabel="Final Heat Mismatch (kWh)",
    )

    for df_m in dfs_MWph
        scatterlines!(ax4,
            df_m.F,
            sum_inm.(df_m.hfp),
            color=m_colors[df_m.name[1]],
            label="$(m_names[df_m.name[1]])")
        scatterlines!(ax4,
            df_m.F,
            -1.0 .* sum_inm.(df_m.hfm),
            color=m_colors[df_m.name[1]],
            label="$(m_names[df_m.name[1]])")
    end

    # If there is no penalty in the system, the automatic y limits are not nice,
    # Set them manually here:

    Legend(fig[:, 2], ax1)
    fig
end

##


for suffix in ["", "_ini"]
    for flex_interval in flex_interval_array
        for n_samples in [6, 12]
            local df_OFIOR = filter_kwd(df; flex_interval, name="OFIOR" * suffix, n_samples)

            local df_OFOR = filter_kwd(df; flex_interval, name="OFOR" * suffix, n_samples)

            local df_OFR = filter_kwd(df; flex_interval, name="OFR" * suffix, n_samples)

            # experiment == "defaults" && (df_OFR.cab_pkWh[3] = missing) # There is a non-zero penalty here.

            local df_CapOnly = filter_kwd(df; name="capacity_only")

            local df_SB = filter_kwd(df;
                name="stochastic_background",
                optimized=true,
                flex_interval,
                n_samples)

            # df_OFIOR_n100 = filter_kwd(df; flex_interval, name="OFIOR", n_samples=100)


            local fig = cost_above_baseline_fig([df_OFIOR, df_OFOR, df_OFR, df_SB, df_CapOnly],
                [df_OFIOR, df_OFOR, df_OFR, df_SB])
            save(joinpath(figure_dir, "cost_above_baseline_n_samples-$(n_samples)_flex_interval-$(flex_interval)$suffix.png"), fig)

        end
    end
end

##

n_samples = 12
suffix = ""
flex_interval = 6

df_OFIOR = filter_kwd(df; flex_interval, name="OFIOR" * suffix, n_samples)

df_OFOR = filter_kwd(df; flex_interval, name="OFOR" * suffix, n_samples)

df_OFR = filter_kwd(df; flex_interval, name="OFR" * suffix, n_samples)

# experiment == "defaults" && (df_OFR.cab_pkWh[3] = missing) # There is a non-zero penalty here.

df_CapOnly = filter_kwd(df; name="capacity_only")

df_SB = filter_kwd(df;
    name="stochastic_background",
    optimized=true,
    flex_interval,
    n_samples)

# df_OFIOR_n100 = filter_kwd(df; flex_interval, name="OFIOR", n_samples=100)

##


function cost_above_baseline_publication(dfs_cab, dfs_MWph)

    update_theme!(fontsize=35)

    tick_strings = string.(F_array)
    tick_strings[2] = ""

    fig = Figure(;
        size=(1600, 600),
        figure_padding=(5, 5, 10, 30),
    )

    ax1 = Axis(fig[1, 1];
        xlabel="Flex Capacity",
        xticks=(F_array, tick_strings), # (1:5, string.(F_array)),
        ylabel="Cost above baseline (k€)")

    for df_m in dfs_cab
        scatterlines!(ax1,
            df_m.F,
            df_m.cost_above_baseline ./ 1e3,
            color=m_colors[df_m.name[1]],
            label="$(m_names[df_m.name[1]])",
            markersize=20)
    end

    ax2 = Axis(fig[1, 2];
        xlabel="Flex Capacity",
        xticks=(F_array, tick_strings), # (1:5, string.(F_array)),
        ylabel="Cost above baseline per MWh (€)")

    for df_m in dfs_MWph
        scatterlines!(ax2,
            df_m.F,
            df_m.cab_pMWh,
            color=m_colors[df_m.name[1]],
            label="$(m_names[df_m.name[1]])",
            markersize=20)
    end

    Legend(fig[1, 3], ax1)
    fig
end


fig = cost_above_baseline_publication([df_OFIOR, df_OFOR, df_OFR, df_SB, df_CapOnly],
    [df_OFIOR, df_OFOR, df_OFR, df_SB])


##

mkpath(joinpath(figure_dir, "publication"))

save(joinpath(figure_dir, "publication", "cost_above_baseline.png"), fig)

##


n_samples = 12
suffix = "_ini"
flex_interval = 6

df_OFIOR = filter_kwd(df; flex_interval, name="OFIOR" * suffix, n_samples)

df_OFOR = filter_kwd(df; flex_interval, name="OFOR" * suffix, n_samples)

df_OFR = filter_kwd(df; flex_interval, name="OFR" * suffix, n_samples)

# experiment == "defaults" && (df_OFR.cab_pkWh[3] = missing) # There is a non-zero penalty here.

df_CapOnly = filter_kwd(df; name="capacity_only")

df_SB = filter_kwd(df;
    name="stochastic_background",
    optimized=true,
    flex_interval,
    n_samples)

fig = cost_above_baseline_publication([df_OFIOR, df_OFOR, df_OFR, df_SB, df_CapOnly],
[df_OFIOR, df_OFOR, df_OFR, df_SB])

save(joinpath(figure_dir, "publication", "cost_above_baseline_no_resample.png"), fig)
