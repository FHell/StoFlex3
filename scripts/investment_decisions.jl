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
df = _df
df = filter_kwd(df, first_dir = experiment)

select!(df,
    [:name,
        :F,
        :flex_interval,
        :n_samples,
        :u_heat_storage,
        :u_heatpump,
        :u_pv,
        :u_storage,
        :u_wind,
        :c_heat_storage,
        :c_heatpump,
        :c_pv,
        :c_storage,
        :c_wind,
    ])

## request_rate_plot(df; model="OFIOR_ini")

##

us = [:u_pv, :u_wind, :u_storage, :u_heat_storage, :u_heatpump]
cs = [:c_pv, :c_wind, :c_storage, :c_heat_storage, :c_heatpump]
invs = [:i_pv, :i_wind, :i_storage, :i_heat_storage, :i_heatpump]
r_invs = [:r_pv, :r_wind, :r_storage, :r_heat_storage, :r_heatpump]
invs_names = ["PV", "Wind", "Elec. Storage", "Heat Storage", "Heat Pump"]

@assert df[1, :name] == "background"

for i in 1:5
    df[!, invs[i]] = df[!, us[i]] .* df[!, cs[i]]
    df[!, r_invs[i]] = df[!, invs[i]] ./ df[1, invs[i]]
end

##

cases = (small = (F = 250, flex_interval = 6, c_name = "Small Requests"),
    large = (F = 5000, flex_interval = 6, c_name = "Large Requests"),
    rare = (F = 5000, flex_interval = 24, c_name = "Rare Requests"))

##

function bar_plot_selector(r; cases, n_samples)
    Fs = [c.F for c in cases]
    FIs = [(c.F, c.flex_interval) for c in cases]
    if r.name == "capacity_only" && r.F in Fs
        return true
    elseif r.name == "background"
        return true
    elseif r.name == "OFIOR" && (r.F, r.flex_interval) in FIs && r.n_samples === n_samples
        return true
    else
        return false
    end
end

df_bp = filter(r -> bar_plot_selector(r; n_samples = 6, cases), df)

select!(df_bp, [:name, :F, :flex_interval, invs..., r_invs...])

##

function cat_exp(r)
    s = ""
    if (r.name == "background")
        return "No Flex Modelling"
    elseif (r.name == "capacity_only")
        if r.F == 250
            return "Capacity Constraint\n (250kW)"
        elseif r.F == 5000
            return "Capacity Constraint\n (5000kW)"
        end
    elseif (r.name == "OFIOR")
        if r.F == 250
            return "Full model\n (250kW)"
        elseif r.F == 5000 && r.flex_interval == 6
            return "Full model\n (5000kW)"
        elseif r.F == 5000 && r.flex_interval == 24
            return "Full model\n (5000kW, infrequent)"
        end
    end
end

df_bp.exp_name = cat_exp.(eachrow(df_bp))

##

df_bp.exp = [1.0, 2.5, 5.0, 3.5, 6.0, 7.0]

df_bar = DataFrames.stack(df_bp, invs)

inv_idx(sym) = findfirst(String.(invs) .== (sym))

df_bar.grp = inv_idx.(df_bar.variable)

##

using CairoMakie

figure_dir = joinpath(@__DIR__, "../plots/", experiment)
mkpath(figure_dir)

##

colors = Makie.wong_colors()

inv_colors = [colors[i] for i in 1:length(invs)];

##

fig = Figure(;
    figure_padding = (5, 5, 10, 10),
    backgroundcolor = :snow2,
    size = (800, 800),)

ax = Axis(fig[1, 1]; xticks = (df_bp.exp, df_bp.exp_name),
    title = "Investment Decisions", xticklabelrotation = -π / 2)

barplot!(ax, df_bar.exp, df_bar.value,
    stack = df_bar.grp,
    color = inv_colors[df_bar.grp],
    label = invs_names[df_bar.grp])

elements = [PolyElement(polycolor = inv_colors[i]) for i in 1:length(invs)]
labels = invs_names
title = "Investments"

Legend(fig[1, 2], elements, labels, title)

dont_display || display(fig)

##

save(joinpath(figure_dir, "investment_decisions.png"), fig)

##

##

update_theme!(fontsize =35)

fig = Figure(;
    figure_padding = (5, 5, 10, 10),
    size = (1000, 600),)

exp_names = [L"$I_B$", L"$I_F$ (small)", L"$I_S$ (small)", L"$I_S$ (rare)", L"$I_F$ (large)", L"$I_S$ (large)"]

ax = Axis(fig[1, 1];
    xticks = (df_bp.exp, exp_names),
    ylabel = "Relative Investment",
    xticklabelrotation = -π / 2)

for i in 3:5
    scatter!(ax,
        df_bp.exp,
        df_bp[!, r_invs[i]],
        color = inv_colors[i],
        marker = :hline,
        markersize = 60)
end

elements = [PolyElement(polycolor = inv_colors[i]) for i in 3:length(invs)]
labels = invs_names[3:5]
title = "Investments"

Legend(fig[1, 2], elements, labels, title)

dont_display || display(fig)

##

mkpath(joinpath(figure_dir, "publication"))

save(joinpath(figure_dir, "publication", "relative_investment_decisions.png"), fig)

##
