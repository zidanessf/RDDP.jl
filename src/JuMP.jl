struct StateInfo
    in::JuMP.VariableInfo
    out::JuMP.VariableInfo
    initial_value::Float64
    kwargs
end
struct State{T}
    # The incoming state variable.
    in::T
    # The outgoing state variable.
    out::T
end
struct UncertainInfo
    var::JuMP.VariableInfo
    kwargs
end
struct Uncertain{T}
    # The outcome.
    outcome::T
end

function JuMP.build_variable(
    _error::Function,
    info::JuMP.VariableInfo,
    ::Type{State};
    initial_value = NaN,
    is_binary = false,
    kwargs...,
)
    if isnan(initial_value)
        _error(
            "When creating a state variable, you must set the " *
            "`initial_value` keyword to the value of the state variable at" *
            " the root node.",
        )
    end
    return StateInfo(
        JuMP.VariableInfo(
            false,
            NaN,  # lower bound
            false,
            NaN,  # upper bound
            false,
            NaN,  # fixed value
            false,
            NaN,  # start value
            false,
            false, # binary and integer
        ),
        info,
        initial_value,
        kwargs,
    )
end
function JuMP.build_variable(
    _error::Function,
    info::JuMP.VariableInfo,
    ::Type{Uncertain};
    kwargs...,
)
    return UncertainInfo(
        info,
        kwargs,
    )
end

function JuMP.add_variable(subproblem::JuMP.Model, state_info::StateInfo, name::String)
    state = State(
        JuMP.add_variable(subproblem, JuMP.ScalarVariable(state_info.in), name * "_in"),
        JuMP.add_variable(subproblem, JuMP.ScalarVariable(state_info.out), name * "_out"),
    )
    # for iter in eachindex(state.in)
    # @info(has_lower_bound(state.in))
    push!(subproblem[:state],state)
    push!(subproblem[:initial_value],state_info.initial_value)
    return state
end
function JuMP.add_variable(subproblem::JuMP.Model, state_info::UncertainInfo, name::String)
    uncertain = JuMP.add_variable(subproblem, JuMP.ScalarVariable(state_info.var), name)
    # if has_lower_bound(uncertain_copy)
    #     delete_lower_bound(uncertain_copy)
    # end
    # if has_upper_bound(uncertain_copy)
    #     delete_upper_bound(uncertain_copy)
    # end
    # for iter in eachindex(state.in)
    push!(subproblem[:uncertain],uncertain)
    return uncertain
end

JuMP.variable_type(model::JuMP.Model, ::Type{State}) = State
JuMP.variable_type(model::JuMP.Model, ::Type{Uncertain}) = Uncertain

function JuMP.value(state::State{JuMP.VariableRef})
    return State(JuMP.value(state.in), JuMP.value(state.out))
end

# Overload for broadcast syntax such as `JuMP.value.([state_1, state_2])`.
Broadcast.broadcastable(state::State{JuMP.VariableRef}) = Ref(state)
Broadcast.broadcastable(unc::Uncertain{JuMP.VariableRef}) = Ref(unc)

# ==============================================================================