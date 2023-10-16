inflate2(f, xs, ys, args...) = [f(x, y, args...) for x in xs, y in ys]

function compute_coefficients2(G::Function, n_max::Int64, m_max::Int64, args...)
    n, m = 0:n_max, 0:m_max
    cₙ, cₘ = exp.(2*pi*im*n/(n_max+1)), exp.(2*pi*im*m/(m_max+1))
    pₙₘ = abs.(fft(fftshift(inflate2(G, cₙ, cₘ, args...)))/((n_max+1)*(m_max+1)))
    return pₙₘ
end