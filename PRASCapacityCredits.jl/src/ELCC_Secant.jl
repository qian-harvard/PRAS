struct ELCC_Secant{M} <: CapacityValuationMethod{M}

    capacity_max::Int
    capacity_gap::Int
    p_value::Float64
    regions::Vector{Tuple{String,Float64}}
    verbose::Bool

    function ELCC_Secant{M}(
        capacity_max::Int, regions::Vector{Pair{String,Float64}};
        capacity_gap::Int=1, p_value::Float64=0.05, verbose::Bool=false) where M

        @assert capacity_max > 0
        @assert capacity_gap > 0
        @assert 0 < p_value < 1
        @assert sum(x.second for x in regions) ≈ 1.0

        return new{M}(capacity_max, capacity_gap, p_value, Tuple.(regions), verbose)

    end

end

function ELCC_Secant{M}(
    capacity_max::Int, region::String; kwargs...
) where M
    return ELCC_Secant{M}(capacity_max, [region=>1.0]; kwargs...)
end

function assess(sys_baseline::S, sys_augmented::S,
                params::ELCC_Secant{M}, simulationspec::SequentialMonteCarlo
) where {N, L, T, P, S <: SystemModel{N,L,T,P}, M <: ReliabilityMetric}

    _, powerunit, _ = unitsymbol(sys_baseline)

    regionnames = sys_baseline.regions.names
    regionnames != sys_augmented.regions.names &&
        error("Systems provided do not have matching regions")

    target_metric = M(first(assess(sys_baseline, simulationspec, Shortfall())))

    capacities = Int[]
    metrics = typeof(target_metric)[]

    elcc_regions, base_load, sys_variable =
        copy_load(sys_augmented, params.regions)

    # Initial points for Secant method
    # Point 0: 0 MW
    c_prev = 0
    update_load!(sys_variable, elcc_regions, base_load, c_prev)
    m_prev = M(first(assess(sys_variable, simulationspec, Shortfall())))
    push!(capacities, c_prev)
    push!(metrics, m_prev)

    # Point 1: Max Capacity
    c_curr = params.capacity_max
    update_load!(sys_variable, elcc_regions, base_load, c_curr)
    m_curr = M(first(assess(sys_variable, simulationspec, Shortfall())))
    push!(capacities, c_curr)
    push!(metrics, m_curr)

    # Ensure the initial Secant points (c_prev, c_curr) bracket the target metric.
    # If they do not, step through capacity in increments of capacity_gap to find
    # a pair of points with a sign change in (metric - target_metric).
    f_prev = m_prev - target_metric
    f_curr = m_curr - target_metric

    if f_prev * f_curr > 0
        # Existing endpoints do not bracket the target; search for a bracketing pair.
        empty!(capacities)
        empty!(metrics)

        # Recompute at 0 MW to reset the first point.
        c_prev = 0
        update_load!(sys_variable, elcc_regions, base_load, c_prev)
        m_prev = M(first(assess(sys_variable, simulationspec, Shortfall())))
        f_prev = m_prev - target_metric
        push!(capacities, c_prev)
        push!(metrics, m_prev)

        step = max(params.capacity_gap, 1)
        c_curr = step
        found_bracket = false

        while c_curr <= params.capacity_max
            update_load!(sys_variable, elcc_regions, base_load, c_curr)
            m_curr = M(first(assess(sys_variable, simulationspec, Shortfall())))
            f_curr = m_curr - target_metric

            push!(capacities, c_curr)
            push!(metrics, m_curr)

            if f_prev * f_curr <= 0
                # Found bracketing interval [c_prev, c_curr].
                found_bracket = true
                break
            end

            # Advance to next interval.
            c_prev = c_curr
            m_prev = m_curr
            f_prev = f_curr
            c_curr += step
        end

        # If no bracketing pair was found, we keep the last two points encountered.
        # This preserves existing behaviour while having explored the space more.
    end
    iter = 0
    max_iter = 100 

    final_val = c_curr

    while iter < max_iter
        iter += 1

        params.verbose && println(
            "Iteration $iter: c_prev=$c_prev, c_curr=$c_curr, m_prev=$(val(m_prev)), m_curr=$(val(m_curr)), target=$(val(target_metric))"
        )

        # g(c) = Metric(c) - Target, i.e., the difference between the metric and the target
        metric_diff_prev = val(m_prev) - val(target_metric)
        metric_diff_curr = val(m_curr) - val(target_metric)

        if abs(metric_diff_curr - metric_diff_prev) < 1e-9
            params.verbose && @info "Denominator too small in Secant method, stopping."
            final_val = c_curr
            break
        end

        # Secant update
        c_next_float = c_curr - metric_diff_curr * (c_curr - c_prev) / (metric_diff_curr - metric_diff_prev)
        c_next = round(Int, c_next_float)

        # Clamp to feasible capacity range [0, capacity_max]
        c_next_clamped = min(max(c_next, 0), params.capacity_max)
        if params.verbose && c_next_clamped != c_next
            @info "Secant update out of bounds (got $c_next, clamped to $c_next_clamped within [0, $(params.capacity_max)])"
        end
        c_next = c_next_clamped

        # Check stopping criteria: Capacity gap
        if abs(c_next - c_curr) <= params.capacity_gap
             params.verbose && @info "Capacity change within tolerance ($(params.capacity_gap)), stopping."
             final_val = c_next
            
             # Evaluate final point if different
             if c_next != c_curr
                update_load!(sys_variable, elcc_regions, base_load, c_next)
                m_next = M(first(assess(sys_variable, simulationspec, Shortfall())))
                push!(capacities, c_next)
                push!(metrics, m_next)
             end
             break
        end

        # Evaluate new point
        update_load!(sys_variable, elcc_regions, base_load, c_next)
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
    
    return CapacityCreditResult{typeof(params), typeof(target_metric), P}(
        target_metric, final_val, final_val, capacities, metrics)

end
