using FFTW, Distributions, DifferentialEquations

include("build.jl")


h(n̄) = truncated(Poisson(n̄); lower = n_min, upper = n_max)  # size distribution
g(m̄) = truncated(Poisson(m̄); lower = m_min, upper = m_max)  # membership distribution

function initialize(ī, ns, ms, n̄, m̄, ϵ)
    C = zeros(Float64, length(ns), ns[end]+1)  # Cₙ,ᵢ
    S = zeros(Float64, length(ms), ms[end]+1)  # Sₘ,ₗ
    I = zeros(Float64, length(ms), ms[end]+1)  # Iₘ,ₗ

    for n_ in eachindex(ns)
        n = ns[n_]
        for i_ in 1:(n+1)
            i = i_ - 1
            C[n_,i_] = pdf.(h(n̄),n)*binomial(n,i)*ϵ^i*(1-ϵ)^(n-i)
        end
    end
    q = sum([n*pdf.(h(n̄),n)*(1-sum([binomial(n-1,i)*ϵ^i*(1-ϵ)^(n-1-i) for i in 0:(ī-1)])) for n in n_min:n_max]) / sum([n*pdf.(h(n̄),n) for n in n_min:n_max])  # active group probability
    for m_ in eachindex(ms)
        m = ms[m_]
        for l_ in 1:(m+1)
            l = l_ - 1
            S[m_,l_] = (1-ϵ)*pdf.(g(m̄),m)*binomial(m,l)*q^l*(1-q)^(m-l)
            I[m_,l_] = ϵ*pdf.(g(m̄),m)*binomial(m,l)*q^l*(1-q)^(m-l)
        end
    end

    return ArrayPartition(C, S, I)
end

# function initialize_ī_to_∞(ns, ms, n̄, m̄, ϵ)
#     C = zeros(Float64, length(ns), ns[end]+1)  # Cₙ,ᵢ
#     S = zeros(Float64, length(ms))  # Sₘ
#     I = zeros(Float64, length(ms))  # Iₘ

#     for n_ in eachindex(ns)
#         n = ns[n_]
#         for i_ in 1:(n+1)
#             i = i_ - 1
#             C[n_,i_] = pdf.(h(n̄),n)*binomial(n,i)*ϵ^i*(1-ϵ)^(n-i)
#         end
#     end
#     for m_ in eachindex(ms)
#         m = ms[m_]
#         S[m_] = (1-ϵ)*pdf.(g(m̄),m)
#         I[m_] = ϵ*pdf.(g(m̄),m)
#     end

#     return ArrayPartition(C, S, I)
# end

function GAME_node!(du, u, p, t)
    ī, ns, ms, k_max, P̄ᵢ, P̄ₛ, P̃ᵢ, P̃ₛ = p
    
    C, S, I = u.x[1], u.x[2], u.x[3]
    dC, dS, dI = du.x[1], du.x[2], du.x[3]

    # joint probability distributions
    compute_coeffs_P̄!(Ēᵢ, Ēₛ, P̄ᵢ, P̄ₛ, C, I, S, ī, ns, ms, k_max)
    compute_coeffs_P̃!(Ẽᵢ, Ẽₛ, P̃ᵢ, P̃ₛ, C, ī, ns, ms, k_max)

    # mean-fields
    θₛ, θᵢ, ϕₛ, ϕᵢ = mean_fields(C, ī, ns, k_max)

    ## dynamic equations
    # cliques
    for n_ in eachindex(ns)
        n = ns[n_]
        for i in 0:n
            dC[n_,i+1] = - (ᾱ(n, i, n_, P̄ᵢ, k_max)*i + β̄(n, i, n_, P̄ₛ, k_max)*(n-i))*C[n_,i+1]
            i > 0 && ( dC[n_,i+1] += β̄(n, i-1, n_, P̄ₛ, k_max)*(n-i+1)*C[n_,i] )
            i < n && ( dC[n_,i+1] += ᾱ(n, i+1, n_, P̄ᵢ, k_max)*(i+1)*C[n_,i+2] )
        end
    end
    # nodes
    for m_ in eachindex(ms)
        m = ms[m_]
        for l_ in 1:(m+1)
            l = l_ - 1
            α̃ᵐₗ = α̃(m_, l_, P̃ᵢ, k_max)
            β̃ᵐₗ  = β̃(m_, l_, P̃ₛ, k_max)

            dS[m_,l_] = α̃ᵐₗ*I[m_,l_] - (β̃ᵐₗ + (m-l)*θₛ + l*ϕₛ)*S[m_,l_]
            dI[m_,l_] = β̃ᵐₗ*S[m_,l_] - (α̃ᵐₗ + (m-l)*θᵢ + l*ϕᵢ)*I[m_,l_]
            l > 0 && ( dS[m_,l_] += (m-l+1)*θₛ*S[m_,l_-1]; dI[m_,l_] += (m-l+1)*θᵢ*I[m_,l_-1] )
            l < m && ( dS[m_,l_] += (l+1)*ϕₛ*S[m_,l_+1];   dI[m_,l_] += (l+1)*ϕᵢ*I[m_,l_+1] )
        end
    end

    if (maximum(abs.(dI)) < 1e-7) || (maximum(abs.(dC)) < 1e-7) || (sum(I) < 1e-7)
        fill!(dC,0.)
        fill!(dS,0.)
        fill!(dI,0.)
    end
end

function GAME_group!(du, u, p, t)
    ī, ns, ms, δ = p
    
    C, S, I = u.x[1], u.x[2], u.x[3]
    dC, dS, dI = du.x[1], du.x[2], du.x[3]

    # mean-fields
    θₛ, θᵢ, ϕₛ, ϕᵢ = mean_fields(C, S, I, ī, ns, ms, δ)

    ## dynamic equations
    # cliques
    for n_ in eachindex(ns)
        n = ns[n_]
        for i in 0:n
            dC[n_,i+1] = - (ᾱ(C, I, n, i, ī, ns, ms)*i + β̄(C, S, n, i, ī, ns, ms, δ)*(n-i))*C[n_,i+1]
            i > 0 && ( dC[n_,i+1] += β̄(C, S, n, i-1, ī, ns, ms, δ)*(n-i+1)*C[n_,i] )
            i < n && ( dC[n_,i+1] += ᾱ(C, I, n, i+1, ī, ns, ms)*(i+1)*C[n_,i+2] )
        end
    end
    # nodes
    for m_ in eachindex(ms)
        m = ms[m_]
        for l_ in 1:(m+1)
            l = l_ - 1
            α̃ᵐₗ = α̃(C, m, l, ī, ns)
            β̃ᵐₗ  = β̃(C, m, l, ī, ns, δ)

            dS[m_,l_] = α̃ᵐₗ*I[m_,l_] - (β̃ᵐₗ + (m-l)*θₛ + l*ϕₛ)*S[m_,l_]
            dI[m_,l_] = β̃ᵐₗ*S[m_,l_] - (α̃ᵐₗ + (m-l)*θᵢ + l*ϕᵢ)*I[m_,l_]
            l > 0 && ( dS[m_,l_] += (m-l+1)*θₛ*S[m_,l_-1]; dI[m_,l_] += (m-l+1)*θᵢ*I[m_,l_-1] )
            l < m && ( dS[m_,l_] += (l+1)*ϕₛ*S[m_,l_+1];   dI[m_,l_] += (l+1)*ϕᵢ*I[m_,l_+1] )
        end
    end

    if (maximum(abs.(dI)) < 1e-7) || (maximum(abs.(dC)) < 1e-7) || (sum(I) < 1e-7)
        fill!(dC,0.)
        fill!(dS,0.)
        fill!(dI,0.)
    end
end

# function GAME_group_ī_to_∞!(du, u, p, t)
#     ns, ms, δ = p
    
#     C, S, I = u.x[1], u.x[2], u.x[3]
#     dC, dS, dI = du.x[1], du.x[2], du.x[3]

#     ## dynamic equations
#     # cliques
#     for n_ in eachindex(ns)
#         n = ns[n_]
#         for i in 0:n
#             dC[n_,i+1] = - (ᾱ(C, I, n, i, ns, ms)*i + β̄(C, S, n, i, ns, ms, δ)*(n-i))*C[n_,i+1]
#             i > 0 && ( dC[n_,i+1] += β̄(C, S, n, i-1, ns, ms, δ)*(n-i+1)*C[n_,i] )
#             i < n && ( dC[n_,i+1] += ᾱ(C, I, n, i+1, ns, ms)*(i+1)*C[n_,i+2] )
#         end
#     end
#     # nodes
#     for m_ in eachindex(ms)
#         m = ms[m_]
#         dS[m_] = α̃(C, m, ns)*I[m_] - β̃(C, m, ns, δ)*S[m_]
#         dI[m_] = -dS[m_]
#     end

#     if (maximum(abs.(dI)) < 1e-7) || (maximum(abs.(dC)) < 1e-7) || (sum(I) < 1e-7)
#         fill!(dC,0.)
#         fill!(dS,0.)
#         fill!(dI,0.)
#     end
# end

function run_GAME_node(p, ī, ns, ms, n̄, m̄; ϵ = 0.001, t_max = 1000)
    u₀ = initialize(ī, ns, ms, n̄, m̄, ϵ)
    tspan = (1, t_max)

    prob = ODEProblem(GAME_node!, u₀, tspan, p)
    return solve(prob, Tsit5(), saveat = 1, reltol=1e-8, abstol=1e-8)
end

function run_GAME_group(p, ī, ns, ms, n̄, m̄; ϵ = 0.001, t_max = 1000)
    u₀ = initialize(ī, ns, ms, n̄, m̄, ϵ)
    tspan = (1, t_max)

    prob = ODEProblem(GAME_group!, u₀, tspan, p)
    return solve(prob, Tsit5(), saveat = 1, reltol=1e-8, abstol=1e-8)
end

# function run_GAME_group_ī_to_∞(p, ns, ms, n̄, m̄; ϵ = 0.001, t_max = 1000)
#     u₀ = initialize_ī_to_∞(ns, ms, n̄, m̄, ϵ)
#     tspan = (1, t_max)

#     prob = ODEProblem(GAME_group_ī_to_∞!, u₀, tspan, p)
#     return solve(prob, Tsit5(), saveat = 1, reltol=1e-8, abstol=1e-8)
# end
