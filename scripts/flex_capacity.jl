df_file = "preloaded/df_raw.serial"

# point the variable df_file to the raw aggregated experimental
# results in order to skip loading from the remote.

include(joinpath(@__DIR__, "load_data_preamble.jl"))

include("../src/load_energy_system_data.jl")
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

df = filter_kwd(df; first_dir = experiment)

CB = get_CB(df)

##

##

using CairoMakie
using LaTeXStrings

figure_dir = joinpath(@__DIR__, "../plots", experiment)
mkpath(figure_dir)

##

function capacity_plot(df_plt; timespan = 1:100, idx = 1)
    update_theme!(fontsize =30)
    sectors_names = ["PV", "Wind", "Elec. Storage", "Heat Sector"]
    pos_flex = [:pv_flex_pos, :wind_flex_pos, :el_flex_pos, :heat_flex_pos]
    neg_flex = [:pv_flex_neg, :wind_flex_neg, :el_flex_neg, :heat_flex_neg]

    colors = Makie.wong_colors()

    sec_colors = [colors[i] for i in 1:length(sectors_names)]

    fig = Figure(;
        figure_padding = (5, 5, 10, 10),
        size = (1000, 800),)

    ax1 = Axis(fig[1, 1];
        xlabel = "Time (h)",
        xticklabelsvisible = false,
        xlabelvisible = false,
        xticksvisible = false,)

    acc = zeros(timespan)

    for (c, r, n) in zip(sec_colors, pos_flex, sectors_names)
        band!(ax1,
            timespan,
            acc[timespan],
            acc[timespan] .+ 1e-3 .* df_plt[idx, r][timespan];
            color = c,
            label = n)
        acc[timespan] .+= 1e-3 .* df_plt[idx, r][timespan]
    end

    ax2 = Axis(fig[2, 1];
        xlabel = "Time (h)")

    acc = zeros(timespan)

    for (c, r, n) in zip(sec_colors, neg_flex, sectors_names)
        band!(ax2,
            timespan,
            acc[timespan],
            acc[timespan] .+ 1e-3 .* df_plt[idx, r][timespan];
            color = c,
            label = n)
        acc[timespan] .+= 1e-3 .* df_plt[idx, r][timespan]
    end

    Label(fig[:, 0],
        "Negative and Positive Capacity (MW)",
        fontsize = 30,
        rotation = pi / 2)

    Legend(fig[:, 2], ax1)

    fig
end

##

df_plt = filter_kwd(df; n_samples = 6, flex_interval = 6, F = 5000, name = "OFIOR")

calculate_flex_capacity!(df_plt, data_es)

##

timespan = 2500:2600
idx = 1

fig = capacity_plot(df_plt; timespan, idx)

##

save(joinpath(figure_dir, "publication", "flex_capacity.png"), fig)
