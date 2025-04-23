df_file = joinpath(@__DIR__, "../preloaded", "df_raw.serial")
include(joinpath(@__DIR__, "load_data_preamble.jl"))

dont_display = true

experiment = ""

for exp in (df.first_dir |> unique)
    global experiment = exp
    include(joinpath(@__DIR__, "convergence.jl"))
    include(joinpath(@__DIR__, "flex_capacity.jl"))
    include(joinpath(@__DIR__, "full_costs_above_baseline.jl"))
    include(joinpath(@__DIR__, "investment_decisions.jl"))
    include(joinpath(@__DIR__, "rate_F_plots.jl"))
    include(joinpath(@__DIR__, "values.jl"))
    include(joinpath(@__DIR__, "recovery_schedule.jl"))
end