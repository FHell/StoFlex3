df_file = "preloaded/df_raw.serial"

# point the variable df_file to the raw aggregated experimental
# results in order to skip loading from the remote.

include(joinpath(@__DIR__, "load_data_preamble.jl"))

##

if ! @isdefined experiment
    experiment = "defaults"
end

if ! @isdefined dont_display
    dont_display = false
end


##

df = filter_kwd(df; first_dir = experiment)

CB = get_CB(df)

df.cost_above_baseline = df.full_cost .- CB .- df.penalty_cost
df.cab_pkWh = df.cost_above_baseline ./ df.expected_flex_energy_per_full_period
df.cab_pMWh = df.cost_above_baseline ./ df.expected_flex_energy_per_full_period .* 1e3

##

sort!(df, :n_samples)

##

using CairoMakie
using LaTeXStrings

figure_dir = joinpath(@__DIR__, "../plots/", experiment)
mkpath(figure_dir)
mkpath(joinpath(figure_dir, "publication"))

##

models_to_plot = ["OFIOR", "OFOR", "OFR"]

include("../src/conventions.jl")

##
update_theme!(fontsize = 30,
markersize = 15)

fig = Figure(;
    figure_padding = (5, 5, 10, 10),
    backgroundcolor = :snow2,
    size = (1200, 1600),)

ax1 = Axis(fig[1, 1];
    xlabel = "Number of Samples",
    ylabel = "Cost above baseline",
    title = "Infrequent small requests")

for name in models_to_plot
    local df_pl = filter_kwd(df; F = 250, flex_interval = 24, name = name)
    scatterlines!(ax1,
        df_pl.n_samples,
        df_pl.cost_above_baseline,
        color = m_colors[name],
        label = "$(m_names[name])")
end

ax2 = Axis(fig[2, 1];
    xlabel = "Number of Samples",
    ylabel = "Cost above baseline",
    title = "Infrequent large requests")

for name in models_to_plot
    local df_pl = filter_kwd(df; F = 5000, flex_interval = 24, name = name)
    scatterlines!(ax2,
        df_pl.n_samples,
        df_pl.cost_above_baseline,
        color = m_colors[name],
        label = "$(m_names[name])")
end

ax3 = Axis(fig[3, 1];
    xlabel = "Number of Samples",
    ylabel = "Cost above baseline",
    title = "Frequent large requests")

for name in models_to_plot
    local df_pl = filter_kwd(df; F = 5000, flex_interval = 6, name = name)
    scatterlines!(ax3,
        df_pl.n_samples,
        df_pl.cost_above_baseline,
        color = m_colors[name],
        label = "$(m_names[name])")
end

Legend(fig[:, 2], ax1)
dont_display || display(fig)
save(joinpath(figure_dir, "Convergence_cost_above_baseline.png"), fig)

##
colors = Makie.wong_colors()

us = [:u_pv, :u_wind, :u_storage, :u_heat_storage, :u_heatpump]
u_names = ["PV", "Wind", "Elec. Storage", "Heat Storage", "Heat Pump"]
u_colors = [colors[i] for i in 1:5];

m_name = "OFIOR"

fig = Figure(;
    figure_padding = (5, 5, 10, 10),
    size = (1200, 1600),)

ax1 = Axis(fig[1, 1];
    limits = (nothing, (0.9, 1.1)),
    xlabelvisible = false,
    xticksvisible = false,
    xticklabelsvisible = false,
    title = "Infrequent small requests (250kW, once per day)")

df_pl = filter_kwd(df; F = 250, flex_interval = 24, name = m_name)

for (u, u_name, u_color) in zip(us, u_names, u_colors)
    scatterlines!(ax1,
        df_pl.n_samples,
        df_pl[!, u] ./ df_pl[end, u],
        color = u_color,
        label = "$(u_name)")
end

ax2 = Axis(fig[2, 1];
    limits = (nothing, (0.9, 1.1)),
    xlabelvisible = false,
    xticksvisible = false,
    xticklabelsvisible = false,
    title = "Infrequent large requests (5000kW, one per day)")

df_pl = filter_kwd(df; F = 5000, flex_interval = 24, name = m_name)

for (u, u_name, u_color) in zip(us, u_names, u_colors)
    scatterlines!(ax2,
        df_pl.n_samples,
        df_pl[!, u] ./ df_pl[end, u],
        color = u_color,
        label = "$(u_name)")
end

ax3 = Axis(fig[3, 1];
    limits = (nothing, (0.9, 1.1)),
    xlabel = "Number of Samples",
    xticks = df_pl.n_samples,
    title = "Frequent large requests (5000kW, four per day)")

df_pl = filter_kwd(df; F = 5000, flex_interval = 6, name = m_name)

for (u, u_name, u_color) in zip(us, u_names, u_colors)
    scatterlines!(ax3,
        df_pl.n_samples,
        df_pl[!, u] ./ df_pl[end, u],
        color = u_color,
        label = "$(u_name)")
end

Legend(fig[:, 2], ax1)
Label(fig[1:3, 0], "Relative change in investment",  rotation = pi/2, fontsize = 40)


dont_display || display(fig)
save(joinpath(figure_dir, "publication", "Convergence_investment_decisions.png"), fig)

##


fig = Figure(;
    figure_padding = (5, 5, 10, 10),
    size = (1200, 1600),)

ax1 = Axis(fig[1, 1];
    xlabelvisible = false,
    xticksvisible = false,
    xticklabelsvisible = false,
#    ylabel = "Flex cost per MWh",
    title = "Infrequent small requests (250kW, one per day)")

for name in models_to_plot
    local df_pl = filter_kwd(df; F = 250, flex_interval = 24, name = name)
    scatterlines!(ax1,
        df_pl.n_samples,
        df_pl.cab_pMWh,
        color = m_colors[name],
        label = "$(m_names[name])")
end

ax2 = Axis(fig[2, 1];
    xlabelvisible = false,
    xticksvisible = false,
    xticklabelsvisible = false,
#    ylabel = "Flex cost per MWh",
    title = "Infrequent large requests (5000kW, one per day)")

for name in models_to_plot
    local df_pl = filter_kwd(df; F = 5000, flex_interval = 24, name = name)
    scatterlines!(ax2,
        df_pl.n_samples,
        df_pl.cab_pMWh,
        color = m_colors[name],
        label = "$(m_names[name])")
end

ax3 = Axis(fig[3, 1];
    xlabel = "Number of Samples",
    xticks = df_pl.n_samples,
#    ylabel = "Flex cost per MWh ",
    title = "Frequent large requests (5000kW, four per day)")

for name in models_to_plot
    local df_pl = filter_kwd(df; F = 5000, flex_interval = 6, name = name)
    scatterlines!(ax3,
        df_pl.n_samples,
        df_pl.cab_pMWh,
        color = m_colors[name],
        label = "$(m_names[name])")
end

Legend(fig[:, 2], ax1)
Label(fig[1:3, 0], "Flexibility cost per MWh",  rotation = pi/2, fontsize = 40)

dont_display || display(fig)
save(joinpath(figure_dir, "publication", "Convergence_cab_pMWh.png"), fig)



##
