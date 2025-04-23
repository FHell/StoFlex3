df_file = "preloaded/df_raw.serial"

# point the variable df_file to the raw aggregated experimental
# results in order to skip loading from the remote.

include(joinpath(@__DIR__, "load_data_preamble.jl"))

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

df = filter_kwd(df; first_dir = experiment)

##

using CairoMakie
using LaTeXStrings

##

abs.(df.gco[3] .- df.gci[3]) |> sum

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
##

fig = plot_split_flows(df, 2)

##

save("gci-gco.png", fig)
##
# sp = load_sp("/home/micha/cluster/coen/micha/stochasticflexibility2/exp_results/F_500/flex_48/n_samples_10/OFR_optimal.serial")
