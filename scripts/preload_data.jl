# This file preloads all data needed for generating the plots into the preload directory.
# This can then be uploaded to Zenodo to ensure reproducibility without access to the HPC
# system.

##

df_file = nothing
_df = nothing

include(joinpath(@__DIR__, "load_data_preamble.jl"))

##

flex_interval = 12
n_samples = 6
F=500
first_dir = "defaults"
df = filter_kwd(df; flex_interval, name = "OFIOR", n_samples, F, first_dir)

##

model_idx = 1

##

sp_optimal = load_sp(joinpath(df.path[model_idx], "$(df.name[model_idx])_optimal.serial"))
serialize("preloaded/sp_rec_schedule_$first_dir.serial", sp_optimal)
