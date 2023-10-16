using Plots, Plots.Measures, LaTeXStrings, JLD

default(legendfont = ("Computer modern", 14),
        tickfont = ("Computer modern", 12),
        guidefontsize = 18, ms = 5, legend_title_font_pointsize = 16,
        linewidth=1, framestyle=:axis, grid=:none,
        background_color_legend = AGray32(1.,0.6),
        bottom_margin = -5mm, left_margin = -2.6mm, right_margin = 5.2mm)
gr(size=(400,400))

include("A_GAME_dyn.jl")

#--- GROUP-BASED ---#
# params
begin
    n̄, m̄ = 4, 3 # average size and membership
    n_min, n_max, m_min, m_max = 2, 8, 3, 3  # min/max size and membership
    nₘ = 8  # max allowed size
    k_max = m_max*(nₘ-1)  # max possible pairwise degree
    ns₀, ms = [n_min:n_max;], [m_min:m_max;]  # size and membership vectors
    t_max = 500
    ξ = 1.  # recovery rate
    īs = [1:8;]
end;


function λ(n, i, δ; ν = ν)
    # complex (linear from threshold ν ≥ 1)
    if i < ν
        return 0.
    else
        return δ*i
    end

    ## simple (ν = 1)
    # return δ*i
end

# params
begin
    ϵ = 0.525  # initial fraction of active nodes
    δ = 0.3  # infection rate
    ν = 2.
end;

# I^star vs ī plot
γₛ, γᵢ = 0.5, 0.
ηs = [0.:0.25:1.;]
I = zeros(length(īs),length(ηs))
for i in 1:length(ηs)
    η = ηs[i]
    res = load("/results/A_GAME/complex_threshold_res_δ_$(δ)_ν_$(ν-1)_η_$(η)_γS_$(γₛ)_γI_$(γᵢ).jld")  # results for ν were stored as "ν-1"
    I_η = [sum(res["$i"][3][end]) for i in 1:8]
    I[:,i] = I_η
end
plot([I[:,i] for i in 1:size(I,2)], labels =:none, width = 3., palette = palette(:Blues)[5:9], yticks = round.([0.0:0.2:0.6;],digits=2),
        xticks = īs, xlabel = L"\bar{i}", ylabel = L"I^\star", title = L"\nu = %$(Int(ν)),\ \delta = %$δ,\ \gamma = %$γₛ");
scatter!([I[:,i] for i in 1:size(I,2)], labels = [ηs;]', legend_title = L"\eta",
        palette = palette(:Blues)[5:9], legend=:bottomright, fg_legend = :transparent)

# average degree plot
η = 1.
res = map(ī -> run_A_GAME_group([ī, nₘ, ms, γₛ, γᵢ, η, δ], ī, nₘ, ns₀, ms, n̄, m̄; ϵ = ϵ, t_max = t_max), īs);
C_η = [[res[i].u[t].x[1] for t in 1:length(res[1].u)] for i in īs]
k_η = m̄.*[[sum([sum(C_η[j][t][i,:]) for i in 1:size(C_η[j][t],1)] .* [i*(i-1) for i in 1:size(C_η[j][t],1)]) ./ sum([sum(C_η[j][t][i,:])*i for i in 1:size(C_η[j][t],1)]) for t in 1:length(C_η[j])] for j in īs]
plot(1:length(k_η[1]), k_η, xscale=:log10, legend=:none, width = 3., labels= [is;]', legend_title = L"\bar{i}", ylims = [9.74,12.02],
    xlabel = L"\textrm{time}", ylabel = L"\textrm{average\ degree}", fg_legend = :transparent, xticks = [1,10,100,1000], yticks = [10,11,12],
    palette = [palette(:Oranges)[4:8];:black], title = L"\nu = %$(Int(ν)),\ \delta = %$δ,\ \gamma = %$γₛ,\ \eta = %$η")

# η vs γ heatmap
ϵ, ν, δ, ī = 0.05, 1., 0.2, 1
ηs = [0.:0.25:1.;]
γs = [0.:0.25:1.;]

res_hm = load("/results/A_GAME/heatmap_δ_$(δ)_ν_$(ν-1).jld", "heatmap")  # results for ν were stored as "ν-1"
heatmap(γs, ηs, res_hm[2,:,:]', ylabel = L"η", xlabel = L"γ", xticks = [0:0.2:1;], tick_direction=:out,
        title = L"\nu = %$(Int(ν)),\ \delta = %$δ,\ \bar{i}=%$(ī)", titlefontsize = 13, c =:viridis, right_margin = 5.4mm);
plot!((γs[indxs1]+γs[indxs1.+1])/2, ηs, color=:black, ls=:solid, label=:none, xlims = [-0.01,1.01], ylims=[-0.0,1.0], width = 1.5)

plot(γs[2:41], [res_hm[2,2:41,i] for i in [1:10:51;]], ylabel = L"I^\star", xlabel = L"γ", title = L"\bar{i} = %$(ī),\ \nu = %$(Int(ν)),\ \delta = %$δ",
    width = 3., legend_title = L"\eta", fg_legend = :transparent, legend=:outerright,
    palette = [:purple;palette(:Greens)[[4;6]];palette(:Oranges)[[4;6;8]]],
    labels = ηs[1:10:51]', xlims = [0.,0.8]);
hline!([res_hm[1,1,1];], color=:black, ls=:dash, width=3, label=:none)



#--- NODE-BASED ---#
# # params
# begin
#     ī = 1  # threshold index
#     n̄, m̄ = 2, 3  # average size and membership
#     n_min, n_max, m_min, m_max = 2, 2, 3, 3 # min/max size and membership
#     nₘ = 4  # max allowed size
#     k_max = m_max*(nₘ-1)  # maximum possible pairwise degree
#     ns₀, ms = [n_min:n_max;], [m_min:m_max;]  # size and membership vectors
#     P̄ᵢ = [[zeros(k_max-n+2, k_max-n+2) for i = 0:n] for n in 1:nₘ]  # initialize PGF P̄ᵢ
#     P̄ₛ = [[zeros(k_max-n+2, k_max-n+2) for i = 0:n] for n in 1:nₘ]  # initialize PGF P̄ₛ
#     P̃ᵢ = [[zeros(k_max+1, k_max+1) for l = 0:m] for m in ms]  # initialize PGF P̃ᵢ
#     P̃ₛ = [[zeros(k_max+1, k_max+1) for l = 0:m] for m in ms]  # initialize PGF P̃ₛ
#     ϵ = 0.01  # initial fraction of active nodes
#     t_max = 100
#     δ, ξ = 0.2, 1.
#     γₛ, γᵢ = 0.05, 0.
#     η = 0.2
#     ν = 1.
# end;

# ī = 1
# p = ī, nₘ, ms, k_max, P̄ᵢ, P̄ₛ, P̃ᵢ, P̃ₛ, γ, η;
# res = run_A_GAME_node(p, ī, nₘ, ns₀, ms, n̄, m̄; ϵ = ϵ, t_max = t_max);
# C, S, I = [res.u[t].x[1] for t in 1:t_max], [res.u[t].x[2] for t in 1:t_max], [res.u[t].x[3] for t in 1:t_max];