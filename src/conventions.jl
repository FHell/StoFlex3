using LaTeXStrings

if ! @isdefined F_array
    const F_array = [250, 500, 1000, 2500, 5000]
    F_dict = Dict(zip(F_array, 1:5))
    F_dict[0] = 6
    const flex_interval_array = [6, 12, 24]
    const n_samples_array = [4, 6, 8, 10, 12]
end

m_colors = Dict("OFIOR" => :black,
    "OFOR" => :blue,
    "OFR" => :orange,
    "capacity_only" => :teal,
    "stochastic_background" => :green,
    "OFIOR_ini" => :black,
    "OFOR_ini" => :blue,
    "OFR_ini" => :orange,
    )

m_names = Dict("OFIOR" => L"C(S)",
    "OFOR" => L"C(S|I_F)",
    "OFR" => L"C(S|I_F, O_F)",
    "capacity_only" => L"C_F",
    "stochastic_background" => L"C(S|I_B)",
    "OFIOR_ini" => L"C'(S)",
    "OFOR_ini" => L"C'(S|I_F)",
    "OFR_ini" => L"C'(S|I_F, O_F)",
)
