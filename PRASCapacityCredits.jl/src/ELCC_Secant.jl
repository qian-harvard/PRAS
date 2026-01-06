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
