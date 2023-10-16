using Random, StatsBase, ProgressMeter

function CM_uniform_hypergraph(N::Int, n::Int, m::Vector{Int})
    vec = Vector{Int}()
    [append!(vec,repeat([i],m[i])) for i in 1:N]
    shuffle!(vec)
    M = length(vec)

    edges = Vector{Vector{Int}}()
    if isinteger(M/n)
        j = 1
        while j < M-n+2
            s = vec[j:(j+n-1)]
            while length(s) != length(unique(s))
                vec[j:end] = shuffle!(vec[j:end])
                s = vec[j:(j+n-1)]
            end
            push!(edges, s)
            j += n
        end
        return edges
    else
        return "non-integer number of edges"
    end
end

function montecarlo_qs(N, n, m, δ, ξ, ν, ϵ, max_timesteps, Δt, average_window, M, update_p)
    edges = CM_uniform_hypergraph(N, n, m)

    L = length(edges)

    σ = zeros(Int,N)
    σ[1:Int(N*ϵ)] .= 1
    shuffle!(σ)

    stored_states = zeros(Int64, N, M)
    for i in 1:M
      stored_states[:,i] = shuffle(σ)
    end
  
    σ_ = copy(σ)
  
    ρ_evo = zeros(Float64, max_timesteps)
    ρ_evo[1] = mean(σ)
  
    ω = ones(N)
  
    @inbounds @showprogress 1 "time" for t in 2:max_timesteps
      fill!(ω,1)
      @simd for e in 1:L
        nodes = edges[e]
        inf = sum(σ[nodes])
        λ_e = λ(n̄, inf, δ; ν = ν)
        for i in nodes
            ω[i] *= 1 - Δt*λ_e
        end
      end
  
      @simd for i in 1:N
        r = rand()
        if σ[i] == 0
          if r < 1 - ω[i]
            σ_[i] = 1
          end
        else
          if r < ξ*Δt
            σ_[i] = 0
          end
        end
      end
  
      ρ_evo[t] = mean(σ_)
  
      if sum(σ_) == 0
        σ .= stored_states[:,rand(1:M)]
      else
        σ .= σ_
        if rand() < update_p
          rand_q = rand(1:M)
          for i in 1:N
            stored_states[i, rand_q] = σ[i]
          end
        end
      end
    end
  
    avg_rho = mean(ρ_evo[(max_timesteps - average_window):end])
    std_rho = std(ρ_evo[(max_timesteps - average_window):end])
    return avg_rho, std_rho, ρ_evo
end