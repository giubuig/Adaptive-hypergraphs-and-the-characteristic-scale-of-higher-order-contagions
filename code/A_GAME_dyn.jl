using FFTW, Distributions, DifferentialEquations

include("build.jl")


h(n̄) = truncated(Poisson(n̄); lower = n_min, upper = n_max)  # size distribution
g(m̄) = truncated(Poisson(m̄); lower = m_min, upper = m_max)  # membership distribution

function initialize(ī, nₘ, ns₀, ms, n̄, m̄, ϵ)
    C = zeros(Float64, nₘ, nₘ+1)  # Cₙ,ᵢ
    S = zeros(Float64, length(ms), ms[end]+1)  # Sₘ,ₗ
    I = zeros(Float64, length(ms), ms[end]+1)  # Iₘ,ₗ

    for n_ in eachindex(ns₀)
        n = ns₀[n_]
        for i_ in 1:(n+1)
            i = i_ - 1
            C[n,i_] = pdf.(h(n̄),n)*binomial(n,i)*ϵ^i*(1-ϵ)^(n-i)
        end
    end
    q = sum([n*pdf.(h(n̄),n)*(1-sum([binomial(n-1,i)*ϵ^i*(1-ϵ)^(n-1-i) for i in 0:(ī-1)])) for n in ns₀[1]:ns₀[end]]) / sum([n*pdf.(h(n̄),n) for n in ns₀[1]:ns₀[end]])  # active group probability
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

function A_GAME_node!(du, u, p, t)
    ī, nₘ, ms, k_max, P̄ᵢ, P̄ₛ, P̃ᵢ, P̃ₛ, γ, η, δ = p
    ns = [1:nₘ;]

    C, S, I = u.x[1], u.x[2], u.x[3]
    dC, dS, dI = du.x[1], du.x[2], du.x[3]

    # joint probability distributions
    compute_coeffs_P̄!(Ēᵢ, Ēₛ, P̄ᵢ, P̄ₛ, C, I, S, ī, ns, ms, k_max)
    compute_coeffs_P̃!(Ẽᵢ, Ẽₛ, P̃ᵢ, P̃ₛ, C, ī, ns, ms, k_max)

    # mean-fields
    θₛ, θᵢ, ϕₛ, ϕᵢ = mean_fields(C, ī, ns, k_max, δ)

    # dynamics
    ## auxiliary quantities
    Cꜜ₋ = sum([sum([C[n,i+1] for i in 0:minimum([ī-1,n])]) for n in 1:(nₘ-1)])
    C₋ = 1 - sum([C[nₘ,i+1] for i in 0:nₘ])

    ī < nₘ ? ΓI_ī = sum([ī*(ī+1)*C[n,ī+2] for n in (ī+1):nₘ])                : ΓI_ī = 0.
    ī < nₘ ? ΩI_ī = sum([ī*C[n,ī+1] for n in ī:(nₘ-1)])                      : ΩI_ī = 0.
    ī < nₘ ? ΩIꜛ = sum([sum([i*C[n,i+1] for n in i:nₘ]) for i in (ī+1):nₘ])  : ΩIꜛ = 0.
    ΩI = sum([sum([i*C[n,i+1] for n in i:(nₘ-1)]) for i in 1:(nₘ-1)])

    ī < nₘ ? ΩS_ī = sum([(n-ī+1)*C[n,ī] for n in ī:(nₘ-1)])                  : ΩS_ī = 0.
    ī < nₘ ? ΩSꜛ = sum([sum([(n-i)*C[n,i+1] for n in i:nₘ]) for i in ī:nₘ])  : ΩSꜛ = 0.
    ΩSꜜ = sum([sum([(n-i)*C[n,i+1] for n in (i+1):(nₘ-1)]) for i in 0:(ī-1)])
    ΩS = sum([sum([(n-i)*C[n,i+1] for n in (i+1):(nₘ-1)]) for i in 0:(nₘ-1)])

    ## equations
    ### cliques
    for n_ in eachindex(ns)
        n = ns[n_]
        for i in 0:n
            # spreading
            dC[n_,i+1] = - (ᾱ(n, i, n_, P̄ᵢ, k_max)*i + β̄(n, i, n_, P̄ₛ, k_max, δ)*(n-i))*C[n_,i+1]
            i > 0 && ( dC[n_,i+1] += β̄(n, i-1, n_, P̄ₛ, k_max, δ)*(n-i+1)*C[n_,i] )
            i < n && ( dC[n_,i+1] += ᾱ(n, i+1, n_, P̄ᵢ, k_max, δ)*(i+1)*C[n_,i+2] )

            # rewiring
            ## nodes leaving
            if i > ī - 1
                dC[n_,i+1] += - γₛ*(n-i)*C[n_,i+1] # (n,i) -> (n-1,i)
                n < nₘ && ( dC[n_,i+1] += γₛ*(n+1-i)*C[n_+1,i+1] + γᵢ*(i+1)*C[n_+1,i+2] ) # (n+1,i) & (n+1,i+1) -> (n,i)
                i > ī && ( dC[n_,i+1] += - γᵢ*i*C[n_,i+1] ) # (n,i) -> (n-1,i-1)
            end
            ## nodes joining
            if Cꜜ₋ > 0.
                if n < nₘ # n -> n+1
                    i < ī ? dC[n_,i+1] += - (η/Cꜜ₋ + (1-η)/C₋)*(γᵢ*ΩIꜛ + γₛ*ΩSꜛ)*C[n_,i+1]  :  dC[n_,i+1] += - (1-η)/C₋*(γᵢ*ΩIꜛ + γₛ*ΩSꜛ)*C[n_,i+1]
                end
                if n > 1 # n-1 -> n
                    dC[n_,i+1] += γₛ*(1-η)/C₋*ΩSꜛ*C[n_-1,i+1]
                    i > 0 && ( dC[n_,i+1] += γᵢ*(1-η)/C₋*ΩIꜛ*C[n_-1,i] )
                    if i < ī
                        dC[n_,i+1] += γₛ*η/Cꜜ₋*ΩSꜛ*C[n_-1,i+1]
                        i > 0 && ( dC[n_,i+1] += γᵢ*η/Cꜜ₋*ΩIꜛ*C[n_-1,i] )
                    elseif i == ī
                        dC[n_,i+1] += γᵢ*η/Cꜜ₋*ΩIꜛ*C[n_-1,i]
                    end
                end
            elseif C₋ > 0.
                if n < nₘ # n -> n+1
                    dC[n_,i+1] += - (1-η)/C₋*(γᵢ*ΩIꜛ + γₛ*ΩSꜛ)*C[n_,i+1]
                end
                if n > 1 # n-1 -> n
                    dC[n_,i+1] += γₛ*(1-η)/C₋*ΩSꜛ*C[n_-1,i+1]
                    i > 0 && ( dC[n_,i+1] += γᵢ*(1-η)/C₋*ΩIꜛ*C[n_-1,i] )
                end
            end
        end
    end
    ### nodes
    for m_ in eachindex(ms)
        m = ms[m_]
        for l_ in 1:(m+1)
            l = l_ - 1
            # spreading
            α̃ᵐₗ = α̃(m_, l_, P̃ᵢ, k_max)
            β̃ᵐₗ  = β̃(m_, l_, P̃ₛ, k_max, δ)

            dS[m_,l_] = α̃ᵐₗ*I[m_,l_] - (β̃ᵐₗ + (m-l)*θₛ + l*ϕₛ)*S[m_,l_]
            dI[m_,l_] = β̃ᵐₗ*S[m_,l_] - (α̃ᵐₗ + (m-l)*θᵢ + l*ϕᵢ)*I[m_,l_]
            l > 0 && ( dS[m_,l_] += (m-l+1)*θₛ*S[m_,l_-1]; dI[m_,l_] += (m-l+1)*θᵢ*I[m_,l_-1] )
            l < m && ( dS[m_,l_] += (l+1)*ϕₛ*S[m_,l_+1];   dI[m_,l_] += (l+1)*ϕᵢ*I[m_,l_+1] )

            # rewiring
            if Cꜜ₋ > 0.
                dS[m_,l_] += - γₛ*(η + (1-η)*Cꜜ₋/C₋)*l*S[m_,l_]
                dI[m_,l_] += - γᵢ*(η + (1-η)*Cꜜ₋/C₋)*l*I[m_,l_]
                l < m && ( dS[m_,l_] += γₛ*(η + (1-η)*Cꜜ₋/C₋)*(l+1)*S[m_,l_+1] ; dI[m_,l_] += γᵢ*(η + (1-η)*Cꜜ₋/C₋)*(l+1)*I[m_,l_+1] )
            end
            if ΩSꜜ > 0.
                dS[m_,l_] += - γᵢ*ΩIꜛ*(η/ΩSꜜ + (1-η)/ΩS)*ΩS_ī*(m-l)*S[m_,l_]
                l > 0 && ( dS[m_,l_] += γᵢ*ΩIꜛ*(η/ΩSꜜ + (1-η)/ΩS)*ΩS_ī*(m-l+1)*S[m_,l_-1] )
            elseif ΩS > 0.
                dS[m_,l_] += - γᵢ*ΩIꜛ*(1-η)/ΩS*ΩS_ī*(m-l)*S[m_,l_]
                l > 0 && ( dS[m_,l_] += γᵢ*ΩIꜛ*(1-η)/ΩS*ΩS_ī*(m-l+1)*S[m_,l_-1] )
            end
            if ΩI > 0.
                dI[m_,l_] += - γᵢ*ΩIꜛ*(1-η)/ΩI*ΩI_ī*(m-l)*I[m_,l_]
                l > 0 && ( dI[m_,l_] += γᵢ*ΩIꜛ*(1-η)/ΩI*ΩI_ī*(m-l+1)*I[m_,l_-1] )
            end
            if ΩIꜛ > 0.
                dI[m_,l_] += - γᵢ*ΓI_ī/ΩIꜛ*l*I[m_,l_]
                l < m && ( dI[m_,l_] += γᵢ*ΓI_ī/ΩIꜛ*(l+1)*I[m_,l_+1] )
            end
        end
    end

    if (maximum(abs.(dI)) < 1e-7) || (maximum(abs.(dC)) < 1e-7) || (sum(I) < 1e-7)
        fill!(dC,0.)
        fill!(dS,0.)
        fill!(dI,0.)
    end
end

function A_GAME_group!(du, u, p, t)
    ī, nₘ, ms, γₛ, γᵢ, η, δ = p
    ns = [1:nₘ;]

    C, S, I = u.x[1], u.x[2], u.x[3]
    dC, dS, dI = du.x[1], du.x[2], du.x[3]

    # mean-fields
    θₛ, θᵢ, ϕₛ, ϕᵢ = mean_fields(C, S, I, ī, ns, ms, δ)
    
    # dynamics
    ## auxuliary quantities
    Cꜜ₋ = sum([sum([C[n,i+1] for i in 0:minimum([ī-1,n])]) for n in 1:(nₘ-1)])
    C₋ = 1 - sum([C[nₘ,i+1] for i in 0:nₘ])

    ī < nₘ ? ΓI_ī = sum([ī*(ī+1)*C[n,ī+2] for n in (ī+1):nₘ])                : ΓI_ī = 0.
    ī < nₘ ? ΩI_ī = sum([ī*C[n,ī+1] for n in ī:(nₘ-1)])                      : ΩI_ī = 0.
    ī < nₘ ? ΩIꜛ = sum([sum([i*C[n,i+1] for n in i:nₘ]) for i in (ī+1):nₘ])  : ΩIꜛ = 0.
    ΩI = sum([sum([i*C[n,i+1] for n in i:(nₘ-1)]) for i in 1:(nₘ-1)])

    ī < nₘ ? ΩS_ī = sum([(n-ī+1)*C[n,ī] for n in ī:(nₘ-1)])                  : ΩS_ī = 0.
    ī < nₘ ? ΩSꜛ = sum([sum([(n-i)*C[n,i+1] for n in i:nₘ]) for i in ī:nₘ])  : ΩSꜛ = 0.
    ΩSꜜ = sum([sum([(n-i)*C[n,i+1] for n in (i+1):(nₘ-1)]) for i in 0:(ī-1)])
    ΩS = sum([sum([(n-i)*C[n,i+1] for n in (i+1):(nₘ-1)]) for i in 0:(nₘ-1)])
    
    ## equations
    ### cliques
    for n_ in eachindex(ns)
        n = ns[n_]
        for i in 0:n
            # spreading
            dC[n_,i+1] = - (ᾱ(C, I, n, i, ī, ns, ms)*i + β̄(C, S, n, i, ī, ns, ms, δ)*(n-i))*C[n_,i+1]
            i > 0 && ( dC[n_,i+1] += β̄(C, S, n, i-1, ī, ns, ms, δ)*(n-i+1)*C[n_,i] )
            i < n && ( dC[n_,i+1] += ᾱ(C, I, n, i+1, ī, ns, ms)*(i+1)*C[n_,i+2] )

            # rewiring
            ## nodes leaving
            if i > ī - 1
                dC[n_,i+1] += - γₛ*(n-i)*C[n_,i+1] # (n,i) -> (n-1,i)
                n < nₘ && ( dC[n_,i+1] += γₛ*(n+1-i)*C[n_+1,i+1] + γᵢ*(i+1)*C[n_+1,i+2] ) # (n+1,i) & (n+1,i+1) -> (n,i)
                i > ī && ( dC[n_,i+1] += - γᵢ*i*C[n_,i+1] ) # (n,i) -> (n-1,i-1)
            end
            ## nodes joining
            if Cꜜ₋ > 0.
                if n < nₘ # n -> n+1
                    i < ī ? dC[n_,i+1] += - (η/Cꜜ₋ + (1-η)/C₋)*(γᵢ*ΩIꜛ + γₛ*ΩSꜛ)*C[n_,i+1]  :  dC[n_,i+1] += - (1-η)/C₋*(γᵢ*ΩIꜛ + γₛ*ΩSꜛ)*C[n_,i+1]
                end
                if n > 1 # n-1 -> n
                    dC[n_,i+1] += γₛ*(1-η)/C₋*ΩSꜛ*C[n_-1,i+1]
                    i > 0 && ( dC[n_,i+1] += γᵢ*(1-η)/C₋*ΩIꜛ*C[n_-1,i] )
                    if i < ī
                        dC[n_,i+1] += γₛ*η/Cꜜ₋*ΩSꜛ*C[n_-1,i+1]
                        i > 0 && ( dC[n_,i+1] += γᵢ*η/Cꜜ₋*ΩIꜛ*C[n_-1,i] )
                    elseif i == ī
                        dC[n_,i+1] += γᵢ*η/Cꜜ₋*ΩIꜛ*C[n_-1,i]
                    end
                end
            elseif C₋ > 0.
                if n < nₘ # n -> n+1
                    dC[n_,i+1] += - (1-η)/C₋*(γᵢ*ΩIꜛ + γₛ*ΩSꜛ)*C[n_,i+1]
                end
                if n > 1 # n-1 -> n
                    dC[n_,i+1] += γₛ*(1-η)/C₋*ΩSꜛ*C[n_-1,i+1]
                    i > 0 && ( dC[n_,i+1] += γᵢ*(1-η)/C₋*ΩIꜛ*C[n_-1,i] )
                end
            end
        end
    end
    ### nodes
    for m_ in eachindex(ms)
        m = ms[m_]
        for l_ in 1:(m+1)
            l = l_ - 1
            # spreading
            α̃ᵐₗ = α̃(C, m, l, ī, ns)
            β̃ᵐₗ  = β̃(C, m, l, ī, ns, δ)

            dS[m_,l_] = α̃ᵐₗ*I[m_,l_] - (β̃ᵐₗ + (m-l)*θₛ + l*ϕₛ)*S[m_,l_]
            dI[m_,l_] = β̃ᵐₗ*S[m_,l_] - (α̃ᵐₗ + (m-l)*θᵢ + l*ϕᵢ)*I[m_,l_]
            l > 0 && ( dS[m_,l_] += (m-l+1)*θₛ*S[m_,l_-1]; dI[m_,l_] += (m-l+1)*θᵢ*I[m_,l_-1] )
            l < m && ( dS[m_,l_] += (l+1)*ϕₛ*S[m_,l_+1];   dI[m_,l_] += (l+1)*ϕᵢ*I[m_,l_+1] )

            # rewiring
            if Cꜜ₋ > 0.
                dS[m_,l_] += - γₛ*(η + (1-η)*Cꜜ₋/C₋)*l*S[m_,l_]
                dI[m_,l_] += - γᵢ*(η + (1-η)*Cꜜ₋/C₋)*l*I[m_,l_]
                l < m && ( dS[m_,l_] += γₛ*(η + (1-η)*Cꜜ₋/C₋)*(l+1)*S[m_,l_+1] ; dI[m_,l_] += γᵢ*(η + (1-η)*Cꜜ₋/C₋)*(l+1)*I[m_,l_+1] )
            end
            if ΩSꜜ > 0.
                dS[m_,l_] += - γᵢ*ΩIꜛ*(η/ΩSꜜ + (1-η)/ΩS)*ΩS_ī*(m-l)*S[m_,l_]
                l > 0 && ( dS[m_,l_] += γᵢ*ΩIꜛ*(η/ΩSꜜ + (1-η)/ΩS)*ΩS_ī*(m-l+1)*S[m_,l_-1] )
            elseif ΩS > 0.
                dS[m_,l_] += - γᵢ*ΩIꜛ*(1-η)/ΩS*ΩS_ī*(m-l)*S[m_,l_]
                l > 0 && ( dS[m_,l_] += γᵢ*ΩIꜛ*(1-η)/ΩS*ΩS_ī*(m-l+1)*S[m_,l_-1] )
            end
            if ΩI > 0.
                dI[m_,l_] += - γᵢ*ΩIꜛ*(1-η)/ΩI*ΩI_ī*(m-l)*I[m_,l_]
                l > 0 && ( dI[m_,l_] += γᵢ*ΩIꜛ*(1-η)/ΩI*ΩI_ī*(m-l+1)*I[m_,l_-1] )
            end
            if ΩIꜛ > 0.
                dI[m_,l_] += - γᵢ*ΓI_ī/ΩIꜛ*l*I[m_,l_]
                l < m && ( dI[m_,l_] += γᵢ*ΓI_ī/ΩIꜛ*(l+1)*I[m_,l_+1] )
            end
        end
    end
            
    if (maximum(abs.(dI)) < 1e-7) || (maximum(abs.(dC)) < 1e-7) || (sum(I) < 1e-7)
        fill!(dC,0.)
        fill!(dS,0.)
        fill!(dI,0.)
    end
end

function run_A_GAME_node(p, ī, nₘ, ns₀, ms, n̄, m̄; ϵ = 0.001, t_max = 1000)
    u₀ = initialize(ī, nₘ, ns₀, ms, n̄, m̄, ϵ)
    tspan = (1, t_max)

    prob = ODEProblem(A_GAME_node!, u₀, tspan, p)
    return solve(prob, Tsit5(), saveat = 1, reltol=1e-8, abstol=1e-8)
end

function run_A_GAME_group(p, ī, nₘ, ns₀, ms, n̄, m̄; ϵ = 0.001, t_max = 1000)
    u₀ = initialize(ī, nₘ, ns₀, ms, n̄, m̄, ϵ)
    tspan = (1, t_max)

    prob = ODEProblem(A_GAME_group!, u₀, tspan, p)
    return solve(prob, Tsit5(), saveat = 1, reltol=1e-8, abstol=1e-8)
end