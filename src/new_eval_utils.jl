using DataFrames
using CSV
using Serialization
using Statistics

function load_results_csv(simpath)
    error("Loading of tables with different formats not implemented.")
    df = DataFrame()
    for (path, _dirs, files) in walkdir(simpath)
        for filename in files
            if endswith(filename, ".csv")
                append!(df, DataFrame(CSV.File(joinpath(path, filename))),
                    cols = :union)
            end
        end
    end
    return df
end

function load_results(simpath)
    df = DataFrame()
    for (path, _dirs, files) in walkdir(simpath)
        for filename in files
            if endswith(filename, "all_data.serial")
                df_file = DataFrame(map(x -> [x], deserialize(joinpath(path, filename))))
                df_file.path .= path
                append!(df, df_file;
                    cols = :union,)
            end
        end
    end
    return df
end

function get_CB(df)
    b = subset(df,
        :name => name -> name .== "background")
    if nrow(b) > 1
        println("More than one background run. Taking mean.")
    end
    Statistics.mean(b.full_cost)
end

function filter_kwd(df; kwargs...)
    function select_this(r)
        for (key, val) in kwargs
            if !(r[key] === val) # Do not select missing values
                return false
            end
        end
        return true
    end
    filter(select_this, df)
end

function calculate_flex_capacity!(df, data_es)
    df[!, :pv_flex_pos] .= df.pv_cur
    df[!, :pv_flex_neg] .= df.pv_cur
    df[!, :wind_flex_pos] .= df.pv_cur
    df[!, :wind_flex_neg] .= df.pv_cur
    df[!, :el_flex_pos] .= df.pv_cur
    df[!, :el_flex_neg] .= df.pv_cur
    df[!, :heat_flex_pos] .= df.pv_cur
    df[!, :heat_flex_neg] .= df.pv_cur
    for dr in eachrow(df)
        if !ismissing(dr.pv_cur)
            fill_flex_capacity(dr, data_es)
        end
    end
end

function fill_flex_capacity(dr::DataFrameRow, data_es)
    dr.pv_flex_pos = dr.pv_cur
    dr.pv_flex_neg = dr.pv_cur - data_es.pv .* dr.u_pv
    dr.wind_flex_pos = dr.wind_cur
    dr.wind_flex_neg = dr.wind_cur - data_es.wind .* dr.u_wind
    dr.el_flex_pos = min.(dr.sto_soc[2:end] .* dr.sto_ef_dis,
        (dr.max_sto_flow * dr.u_storage .-
         (dr.sto_to_bus[1:(end - 1)] .- dr.sto_from_bus[1:(end - 1)])))
    dr.el_flex_neg = max.((dr.sto_soc[2:end] .- dr.u_storage) ./ dr.sto_ef_dis,
        0.0 .- ((dr.max_sto_flow * dr.u_storage) .-
         (dr.sto_from_bus[1:(end - 1)] .- dr.sto_to_bus[1:(end - 1)])))
    dr.heat_flex_pos = min.(dr.heat_sto_soc[2:end] .* (dr.heat_eff / dr.COP),
        dr.flow_energy2heat[1:(end - 1)])
    dr.heat_flex_neg = max.((dr.heat_sto_soc[2:end] .- dr.u_heat_storage) ./
                            (dr.COP * dr.heat_eff),
        dr.flow_energy2heat[1:(end - 1)] .- (dr.u_heatpump / dr.COP))
end
