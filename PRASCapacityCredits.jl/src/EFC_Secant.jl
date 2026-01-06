struct EFC_Secant{M} <: CapacityValuationMethod{M}

    capacity_max::Int
    capacity_gap::Int
    p_value::Float64
    regions::Vector{Tuple{String,Float64}}
    verbose::Bool

    function EFC_Secant{M}(
        capacity_max::Int, regions::Vector{Pair{String,Float64}};
        capacity_gap::Int=1, p_value::Float64=0.05, verbose::Bool=false) where M

        @assert capacity_max > 0
        @assert capacity_gap > 0
        @assert 0 < p_value < 1
        @assert sum(x.second for x in regions) ≈ 1.0

        return new{M}(capacity_max, capacity_gap, p_value, Tuple.(regions), verbose)

    end

end

function EFC_Secant{M}(
    capacity_max::Int, region::String; kwargs...
) where M
    return EFC_Secant{M}(capacity_max, [region=>1.0]; kwargs...)
end

function assess(sys_baseline::S, sys_augmented::S,
                params::EFC_Secant{M}, simulationspec::SequentialMonteCarlo
) where {N, L, T, P, S <: SystemModel{N,L,T,P}, M <: ReliabilityMetric}

    _, powerunit, _ = unitsymbol(sys_baseline)

    regionnames = sys_baseline.regions.names
    regionnames != sys_augmented.regions.names &&
        error("Systems provided do not have matching regions")

    # Add firm capacity generators to the relevant regions
    efc_gens, sys_variable, sys_target =
        add_firmcapacity(sys_baseline, sys_augmented, params.regions)

    target_metric = M(first(assess(sys_target, simulationspec, Shortfall())))

    capacities = Int[]
    metrics = typeof(target_metric)[]

    # Initial points for Secant method
    # Start at 0 MW and search for a point that brackets the target metric
    c_prev = 0
    update_firmcapacity!(sys_variable, efc_gens, c_prev)
    m_prev = M(first(assess(sys_variable, simulationspec, Shortfall())))
    f_prev = m_prev - target_metric

    # Try using capacity_max as the second point first
    bracket_found = false
    c_curr = params.capacity_max
    update_firmcapacity!(sys_variable, efc_gens, c_curr)
    m_curr = M(first(assess(sys_variable, simulationspec, Shortfall())))
    f_curr = m_curr - target_metric

    if f_prev * f_curr <= 0
        bracket_found = true
    else
        # Scan capacities in steps of capacity_gap to find a bracketing pair
        for c in params.capacity_gap:params.capacity_gap:params.capacity_max
            c_curr = c
            update_firmcapacity!(sys_variable, efc_gens, c_curr)
            m_curr = M(first(assess(sys_variable, simulationspec, Shortfall())))
            f_curr = m_curr - target_metric

            if f_prev * f_curr <= 0
                bracket_found = true
                break
            end

            # Move forward: current point becomes previous for next iteration
            c_prev = c_curr
            m_prev = m_curr
            f_prev = f_curr
        end
    end

    if !bracket_found
        error("EFC_Secant initialization failed to bracket target metric within [0, capacity_max]. " *
              "Consider increasing capacity_max or adjusting capacity_gap.")
    end

    # Record the bracketing points as initial values for the Secant method
    empty!(capacities)
    empty!(metrics)
    push!(capacities, c_prev)
    push!(metrics, m_prev)
    push!(capacities, c_curr)
    push!(metrics, m_curr)

    iter = 0
    max_iter = 100 
    
    final_val = c_curr

    while iter < max_iter
        iter += 1

        params.verbose && println(
            "Iteration $iter: c_prev=$c_prev, c_curr=$c_curr, m_prev=$(val(m_prev)), m_curr=$(val(m_curr)), target=$(val(target_metric))"
        )

        # g(c) = Metric(c) - Target
        g_prev = val(m_prev) - val(target_metric)
        g_curr = val(m_curr) - val(target_metric)

        if abs(g_curr - g_prev) < 1e-9
            params.verbose && @info "Denominator too small in Secant method, stopping."
            final_val = c_curr
            break
        end

        # Secant update
        c_next_float = c_curr - g_curr * (c_curr - c_prev) / (g_curr - g_prev)
        c_next = clamp(round(Int, c_next_float), 0, params.capacity_max)

        # Check stopping criteria: Capacity gap
        if abs(c_next - c_curr) <= params.capacity_gap
             params.verbose && @info "Capacity change within tolerance ($(params.capacity_gap)), stopping."
             final_val = c_next

             # Evaluate final point if different
             if c_next != c_curr
                update_firmcapacity!(sys_variable, efc_gens, c_next)
                m_next = M(first(assess(sys_variable, simulationspec, Shortfall())))
                push!(capacities, c_next)
                push!(metrics, m_next)
             end
             break
        end

        # Evaluate new point
        update_firmcapacity!(sys_variable, efc_gens, c_next)
        m_next = M(first(assess(sys_variable, simulationspec, Shortfall())))
        push!(capacities, c_next)
        push!(metrics, m_next)

        # Setup for next iteration
        c_prev = c_curr
        m_prev = m_curr
        c_curr = c_next
        m_curr = m_next
        final_val = c_curr
        
    end

    if iter >= max_iter && params.verbose
        @warn "EFC_Secant did not converge within maximum iterations ($(max_iter)); using last capacity value $(final_val)."
    end
    return CapacityCreditResult{typeof(params), typeof(target_metric), P}(
        target_metric, final_val, final_val, capacities, metrics)

end
