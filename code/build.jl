include("DFT.jl")

#--- NODE-BASED ---#

function Kᵢꜛ(x, y, C, ī, ns)
    z, Z = 0., 0.
    for n_ in eachindex(ns)
        n = ns[n_]
        if ī < n
            z += sum([i*C[n_,i+1]*x^(n-1)*y^(i-1) for i in (ī+1):n])
            Z += sum([i*C[n_,i+1] for i in (ī+1):n])
        end
    end
    Z > 0 ? z /= Z : z = 0.
    return z
end

function Kᵢꜜ(x, y, C, ī, ns)
    z, Z = 0., 0.
    for n_ in eachindex(ns)
        n = ns[n_]
        z += sum([i*C[n_,i+1]*x^(n-1)*y^(i-1) for i in 1:minimum([ī,n])])
        Z += sum([i*C[n_,i+1] for i in 1:minimum([ī,n])])
    end
    Z > 0 ? z /= Z : z = 0.
    return z
end

function Kₛꜛ(x, y, C, ī, ns)
    z, Z = 0., 0.
    for n_ in eachindex(ns)
        n = ns[n_]
        if ī < n
            z += sum([(n-i)*C[n_,i+1]*x^(n-1)*y^i for i in ī:(n-1)])
            Z += sum([(n-i)*C[n_,i+1] for i in ī:(n-1)])
        end
    end
    Z > 0 ? z /= Z : z = 0.
    return z
end

function Kₛꜜ(x, y, C, ī, ns)
    z, Z = 0., 0.
    for n_ in eachindex(ns)
        n = ns[n_]
        z += sum([(n-i)*C[n_,i+1]*x^(n-1)*y^i for i in 0:minimum([ī-1,n-1])])
        Z += sum([(n-i)*C[n_,i+1] for i in 0:minimum([ī-1,n-1])])
    end
    Z > 0 ? z /= Z : z = 0.
    return z
end

function Gᵢ(x, y, I, ī, i, ms)
    z, Z = 0., 0.
    if i < ī + 1
        for m_ in eachindex(ms)
            m = ms[m_]
            z += sum([(m-l)*I[m_,l+1]*x^(m-l-1)*y^l for l in 0:(m-1)])
            Z += sum([(m-l)*I[m_,l+1] for l in 0:(m-1)])
        end
    else
        for m_ in eachindex(ms)
            m = ms[m_]
            z += sum([l*I[m_,l+1]*x^(m-l)*y^(l-1) for l in 1:m])
            Z += sum([l*I[m_,l+1] for l in 1:m])
        end
    end
    Z > 0 ? z /= Z : z = 0.
    return z
end

function Gₛ(x, y, S, ī, i, ms)
    z, Z = 0., 0.
    if i < ī
        for m_ in eachindex(ms)
            m = ms[m_]
            z += sum([(m-l)*S[m_,l+1]*x^(m-l-1)*y^l for l in 0:(m-1)])
            Z += sum([(m-l)*S[m_,l+1] for l in 0:(m-1)])
        end
    else
        for m_ in eachindex(ms)
            m = ms[m_]
            z += sum([l*S[m_,l+1]*x^(m-l)*y^(l-1) for l in 1:m])
            Z += sum([l*S[m_,l+1] for l in 1:m])
        end
    end
    Z > 0 ? z /= Z : z = 0.
    return z
end

function Ēᵢ(x, y, C, I, ī, i, ns, ms)
    return Gᵢ(Kᵢꜜ(x, y, C, ī, ns), Kᵢꜛ(x, y, C, ī, ns), I, ī, i, ms)
end

function Ēₛ(x, y, C, S, ī, i, ns, ms)
    return Gₛ(Kₛꜜ(x, y, C, ī, ns), Kₛꜛ(x, y, C, ī, ns), S, ī, i, ms)
end

function Ẽᵢ(x, y, C, ī, m, l, ns)
    return (Kᵢꜜ(x, y, C, ī, ns))^(m-l) * (Kᵢꜛ(x, y, C, ī, ns))^l
end

function Ẽₛ(x, y, C, ī, m, l, ns)
    return (Kₛꜜ(x, y, C, ī, ns))^(m-l) * (Kₛꜛ(x, y, C, ī, ns))^l
end

function α(k, l; ξ = ξ)
    return ξ
end

function β(k, l, δ)
    # complex (linear from threshold ν ≥ 1)
    if i < ν
        return 0.
    else
        return δ*l
    end

    ## simple (ν = 1)
    # return δ*l
end

function ᾱ(n::Int, i::Int, n_::Int, P̄ᵢ::Vector{Vector{Matrix{Float64}}}, k_max::Int)
    # return sum([sum([α(n-1+r,i-1+j) * P̄ᵢ[n_][i+1][r+1,j+1] for j in 0:r]) for r in 0:(k_max-n+1)])
    return ξ
end

function β̄(n::Int, i::Int, n_::Int, P̄ₛ::Vector{Vector{Matrix{Float64}}}, k_max::Int, δ)
    return sum([sum([β(n-1+r,i+j,δ) * P̄ₛ[n_][i+1][r+1,j+1] for j in 0:r]) for r in 0:(k_max-n+1)])
end

function α̃(m_::Int, l_::Int, P̃ᵢ::Vector{Vector{Matrix{Float64}}}, k_max::Int)
    # return sum([sum([α(k,ℓ) * P̃ᵢ[m_][l_][k+1,ℓ+1] for ℓ in 0:k]) for k in 0:k_max])
    return ξ
end

function β̃(m_::Int, l_::Int, P̃ₛ::Vector{Vector{Matrix{Float64}}}, k_max::Int, δ)
    return sum([sum([β(k,ℓ,δ) * P̃ₛ[m_][l_][k+1,ℓ+1] for ℓ in 0:k]) for k in 0:k_max])
end

function mean_fields(C::Matrix{Float64}, ī::Int, ns::Vector{Int}, k_max::Int, δ)
    mf, Z = zeros(4), zeros(4)
    if ī < ns[end]
        for n_ in eachindex(ns)
            n = ns[n_]

            Z[1] += sum([(n-i)*C[n_,i+1] for i in 0:minimum([ī-1,n-1])])
            Z[2] += sum([i*C[n_,i+1] for i in 1:minimum([ī,n])])
            if ī < n
                mf[1] += (n-ī+1)*C[n_,ī]*(n-ī)*β̄(n,ī-1,n_,P̄ₛ,k_max,δ)
                mf[2] += ī*C[n_,ī+1]*(n-ī)*β̄(n,ī,n_,P̄ₛ,k_max,δ)
                mf[3] += (n-ī)*C[n_,ī+1]*ī*ᾱ(n,ī,n_,P̄ᵢ,k_max)
                mf[4] += (ī+1)*C[n_,ī+2]*ī*ᾱ(n,ī+1,n_,P̄ᵢ,k_max)
                Z[3] += sum([(n-i)*C[n_,i+1] for i in ī:(n-1)])
                Z[4] += sum([i*C[n_,i+1] for i in (ī+1):n])
            end
        end
        for i in 1:4
            Z[i] > 0 ? mf[i] /= Z[i] : mf[i] = 0.
        end
    end
    return mf
end

function compute_coeffs_P̄!(Ēᵢ, Ēₛ, P̄ᵢ, P̄ₛ, C, I, S, ī, ns, ms, k_max)
    for n_ in eachindex(ns)
        n = ns[n_]
        fill!(P̄ᵢ[n_][1], 0.)
        P̄ₛ[n_][1] = compute_coefficients2(Ēₛ, k_max-n+1, k_max-n+1, C, S, ī, 0, ns, ms)
        for i_ in 2:n
            i = i_ - 1
            P̄ᵢ[n_][i_] = compute_coefficients2(Ēᵢ, k_max-n+1, k_max-n+1, C, I, ī, i, ns, ms)
            P̄ₛ[n_][i_] = compute_coefficients2(Ēₛ, k_max-n+1, k_max-n+1, C, S, ī, i, ns, ms)
        end
        fill!(P̄ₛ[n_][n+1], 0.)
        P̄ᵢ[n_][n+1] = compute_coefficients2(Ēᵢ, k_max-n+1, k_max-n+1, C, I, ī, n, ns, ms)
    end
end

function compute_coeffs_P̃!(Ẽᵢ, Ẽₛ, P̃ᵢ, P̃ₛ, C, ī, ns, ms, k_max)
    for m_ in eachindex(ms)
        m = ms[m_]
        for l_ in 1:(m+1)
            l = l_ - 1
            P̃ᵢ[m_][l_] = compute_coefficients2(Ẽᵢ, k_max, k_max, C, ī, m, l, ns)
            P̃ₛ[m_][l_] = compute_coefficients2(Ẽₛ, k_max, k_max, C, ī, m, l, ns)
        end
    end
end




#--- GROUP-BASED ---#

function μ(n, i; ξ = ξ)
    return ξ
end

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

function μ̄ꜛ(C, ī, ns)
    x, y = 0., 0.
    if ī < ns[end]
        for n_ in eachindex(ns)
            n = ns[n_]
            x += sum([i * C[n_,i+1] * μ(n,i) for i in (ī+1):n])
            y += sum([i * C[n_,i+1] for i in (ī+1):n])
        end
    end
    y > 0 ? x /= y : x = 0.
    return x
end

function μ̄ꜜ(C, ī, ns)
    x, y = 0., 0.
    for n_ in eachindex(ns)
        n = ns[n_]
        x += sum([i * C[n_,i+1] * μ(n,i) for i in 1:minimum([ī,n])])
        y += sum([i * C[n_,i+1] for i in 1:minimum([ī,n])])
    end
    y > 0 ? x /= y : x = 0.
    return x
end

function λ̄ꜛ(C, ī, ns, δ)
    x, y = 0., 0.
    if ī < ns[end]
        for n_ in eachindex(ns)
            n = ns[n_]
            if ī < n
                x += sum([(n-i) * C[n_,i+1] * λ(n,i,δ) for i in ī:(n-1)])
                y += sum([(n-i) * C[n_,i+1] for i in ī:(n-1)])
            end
        end
    end
    y > 0 ? x /= y : x = 0.
    return x
end

function λ̄ꜜ(C, ī, ns, δ)
    x, y = 0., 0.
    for n_ in eachindex(ns)
        n = ns[n_]
        x += sum([(n-i) * C[n_,i+1] * λ(n,i,δ) for i in 0:minimum([ī-1,n-1])])
        y += sum([(n-i) * C[n_,i+1] for i in 0:minimum([ī-1,n-1])])
    end
    y > 0 ? x /= y : x = 0.
    return x
end

function μ̄(C, I, i, ī, ns, ms)
    d, u = μ̄ꜜ(C, ī, ns), μ̄ꜛ(C, ī, ns)
    x, y = 0., 0.
    if i > ī
        for m_ in eachindex(ms)
            m = ms[m_]
            x += sum([l * I[m_,l+1] * ((m-l)*d + (l-1)*u) for l in 1:m])
            y += sum([l * I[m_,l+1] for l in 1:m])
        end
    else
        for m_ in eachindex(ms)
            m = ms[m_]
            x += sum([(m-l) * I[m_,l+1] * ((m-l-1)*d + l*u) for l in 0:(m-1)])
            y += sum([(m-l) * I[m_,l+1] for l in 0:(m-1)])
        end
    end
    y > 0 ? x /= y : x = 0.
    return x
end

function λ̄(C, S, i, ī, ns, ms, δ)
    d, u = λ̄ꜜ(C, ī, ns, δ), λ̄ꜛ(C, ī, ns, δ)
    x, y = 0., 0.
    if i > ī - 1
        for m_ in eachindex(ms)
            m = ms[m_]
            x += sum([l * S[m_,l+1] * ((m-l)*d + (l-1)*u) for l in 1:m])
            y += sum([l * S[m_,l+1] for l in 1:m])
        end
    else
        for m_ in eachindex(ms)
            m = ms[m_]
            x += sum([(m-l) * S[m_,l+1] * ((m-l-1)*d + l*u) for l in 0:(m-1)])
            y += sum([(m-l) * S[m_,l+1] for l in 0:(m-1)])
        end
    end
    y > 0 ? x /= y : x = 0.
    return x
end

function ᾱ(C::Matrix{Float64}, I::Matrix{Float64}, n::Int, i::Int, ī::Int, ns::Vector{Int}, ms::Vector{Int}; ξ = ξ)
    # return μ(n, i) + μ̄(C, I, i, ī, ns, ms)
    return ξ
end

function β̄(C::Matrix{Float64}, S::Matrix{Float64}, n::Int, i::Int, ī::Int, ns::Vector{Int}, ms::Vector{Int}, δ)
    return λ(n, i, δ) + λ̄(C, S, i, ī, ns, ms, δ)
end

function α̃(C::Matrix{Float64}, m::Int, l::Int, ī::Int, ns::Vector{Int}; ξ = ξ)
    # return (m-l)*μ̄ꜜ(C, ī, ns) + l*μ̄ꜛ(C, ī, ns)
    return ξ
end

function β̃(C::Matrix{Float64}, m::Int, l::Int, ī::Int, ns::Vector{Int}, δ)
    return (m-l)*λ̄ꜜ(C, ī, ns, δ) + l*λ̄ꜛ(C, ī, ns, δ)
end

function mean_fields(C::Matrix{Float64}, S::Matrix{Float64}, I::Matrix{Float64}, ī::Int, ns::Vector{Int}, ms::Vector{Int}, δ)
    mf, Z = zeros(4), zeros(4)
    if ī < ns[end]
        for n_ in eachindex(ns)
            n = ns[n_]

            Z[1] += sum([(n-i)*C[n_,i+1] for i in 0:minimum([ī-1,n-1])])
            Z[2] += sum([i*C[n_,i+1] for i in 1:minimum([ī,n])])
            if ī < n
                mf[1] += (n-ī+1)*C[n_,ī]*(n-ī)*β̄(C,S,n,ī-1,ī,ns,ms, δ)
                mf[2] += ī*C[n_,ī+1]*(n-ī)*β̄(C,S,n,ī,ī,ns,ms, δ)
                mf[3] += (n-ī)*C[n_,ī+1]*ī*ᾱ(C,I,n,ī,ī,ns,ms)
                mf[4] += (ī+1)*C[n_,ī+2]*ī*ᾱ(C,I,n,ī+1,ī,ns,ms)
                Z[3] += sum([(n-i)*C[n_,i+1] for i in ī:(n-1)])
                Z[4] += sum([i*C[n_,i+1] for i in (ī+1):n])
            end
        end
        for i in 1:4
            Z[i] > 0 ? mf[i] /= Z[i] : mf[i] = 0.
        end
    end
    return mf
end



### ī → ∞ (equivalent to not track group activity)
# function μ̄ꜜ(C, ns)
#     x, y = 0., 0.
#     for n_ in eachindex(ns)
#         n = ns[n_]
#         x += sum([i * C[n_,i+1] * μ(n,i) for i in 1:n])
#         y += sum([i * C[n_,i+1] for i in 1:n])
#     end
#     y > 0 ? x /= y : x = 0.
#     return x
# end

# function λ̄ꜜ(C, ns, δ)
#     x, y = 0., 0.
#     for n_ in eachindex(ns)
#         n = ns[n_]
#         x += sum([(n-i) * C[n_,i+1] * λ(n,i, δ) for i in 0:(n-1)])
#         y += sum([(n-i) * C[n_,i+1] for i in 0:(n-1)])
#     end
#     y > 0 ? x /= y : x = 0.
#     return x
# end

# function μ̄(C, I, ns, ms)
#     x, y = 0., 0.
#     for m_ in eachindex(ms)
#         m = ms[m_]
#         x += m * I[m_] * (m-1)
#         y += m * I[m_]
#     end
#     y > 0 ? x = μ̄ꜜ(C, ns) * x/y : x = 0.
#     return x
# end

# function λ̄(C, S, ns, ms, δ)
#     x, y = 0., 0.
#     for m_ in eachindex(ms)
#         m = ms[m_]
#         x += m * S[m_] * (m-1)
#         y += m * S[m_]
#     end
#     y > 0 ? x = λ̄ꜜ(C, ns, δ)* x/y : x = 0.
#     return x
# end

# function ᾱ(C::Matrix{Float64}, I::Vector{Float64}, n::Int, i::Int, ns::Vector{Int}, ms::Vector{Int}; ξ = ξ)
#     # return μ(n, i) + μ̄(C, I, ns, ms)
#     return ξ
# end

# function β̄(C::Matrix{Float64}, S::Vector{Float64}, n::Int, i::Int, ns::Vector{Int}, ms::Vector{Int}, δ)
#     return λ(n, i, δ) + λ̄(C, S, ns, ms, δ)
# end

# function α̃(C::Matrix{Float64}, m::Int, ns::Vector{Int}; ξ = ξ)
#     # return m*μ̄ꜜ(C, ns)
#     return ξ
# end

# function β̃(C::Matrix{Float64}, m::Int, ns::Vector{Int}, δ)
#     return m*λ̄ꜜ(C, ns, δ)
# end