using Plots, Plots.Measures, LaTeXStrings

include("GAME_dyn.jl")

#--- GROUP-BASED ---#
# params
begin
    n̄, m̄ = 5, 3  # average size and membership
    n_min, n_max, m_min, m_max = 5, 5, 3, 3 # min/max size and membership
    ns, ms = [n_min:n_max;], [m_min:m_max;]  # size and membership vectors
    īs = [1:5;]
end;

# I^star vs δ plot

ν = 3.
# δs = [0.0925:0.0001:0.0974;]  # ν = 1
# δs = [0.205:0.0005:0.22;]  # ν = 2
δs = [0.395:0.001:0.43;]  # ν = 3

# ϵ = 0.06  # ν = 1
ϵ = 0.8  # ν = 2, 3

res_store = []
for δ in δs
    res = load("/results/GAME/complex_threshold_res_δ_$(δ)_ν_$(ν-1)_ϵ_$(ϵ).jld")  # results for ν were stored as "ν-1"
    res = [res["$i"][3][end] for i in īs]
    push!(res_store, res)
end
I = zeros(length(īs),length(δs))
for i in īs
    I[i,:] = [sum(res_store[j][i]) for j in 1:length(δs)]
end

I_mc, σ_I_mc = load("/results/GAME/MC/MC_complex_threshold_res_ν_$(ν-1)_ϵ_$(ϵ).jld", "I", "σ_I")  # results for ν were stored as "ν-1"

# palette_ = [palette(:Oranges)[[4]];reverse(palette(:Blues)[[3;5;7;9]]);]  # ν = 1
# palette_ = [palette(:Blues)[[9]];palette(:Oranges)[[4]];reverse(palette(:Blues)[[3;5;7]]);]  # ν = 2
palette_ = [reverse(palette(:Blues)[[7;9]]);palette(:Oranges)[[4]];reverse(palette(:Blues)[[3;5]]);]  # ν = 3

# δs_mc = [0.093:0.001:0.097;]  # ν = 1
# δs_mc = [0.205:0.0025:0.2175;]  # ν = 2
δs_mc = [0.395:0.005:0.425;]  # ν = 3

plot(δs, [I'[:,i] for i in 1:size(I,1)], xlabel = L"\delta", width = 3.5, ls = [:solid :solid :dot :dot :dot],
    ylabel = L"I^{\star}", legend_title = L"\bar{i}", label = [īs;]', title = L"\nu\ =\ %$(Int(ν))",
    palette=palette_, legend = (0.96,0.6), fg_legend = :transparent)
scatter!(δs_mc, I_mc, yerr = σ_I_mc, ms = 6., labels =:none, color =:black)



#--- NODE-BASED ---#
# # params
# @everywhere begin
#     ī = 1  # threshold index
#     n̄, m̄ = 5, 3  # average size and membership
#     n_min, n_max, m_min, m_max = 5, 5, 3, 3 # min/max size and membership
#     k_max = m_max*(n_max-1)  # maximum possible pairwise degree
#     ns, ms = [n_min:n_max;], [m_min:m_max;]  # size and membership vectors
#     P̄ᵢ = [[zeros(k_max-n+2, k_max-n+2) for i = 0:n] for n in ns]  # initialize PGF P̄ᵢ
#     P̄ₛ = [[zeros(k_max-n+2, k_max-n+2) for i = 0:n] for n in ns]  # initialize PGF P̄ₛ
#     P̃ᵢ = [[zeros(k_max+1, k_max+1) for l = 0:m] for m in ms]  # initialize PGF P̃ᵢ
#     P̃ₛ = [[zeros(k_max+1, k_max+1) for l = 0:m] for m in ms]  # initialize PGF P̃ₛ
#     ϵ = 0.8  # initial fraction of active nodes
#     t_max = 100
#     δ, ξ = 0.095, 1.
#     ν = 1.
# end;

# ī = 1
# p = ī, ns, ms, k_max, P̄ᵢ, P̄ₛ, P̃ᵢ, P̃ₛ;
# res = run_GAME_node(p, ī, ns, ms, n̄, m̄; ϵ = ϵ, t_max = t_max);
# C, S, I = [res.u[t].x[1] for t in 1:t_max], [res.u[t].x[2] for t in 1:t_max], [res.u[t].x[3] for t in 1:t_max];