include(joinpath(@__DIR__, "chapter14_common.jl"))

function setit!(state, value)
    set_kick!(state.model, "Q_KICK"; vkick=value)
    return state
end

state = chapter14_state()
setit!(state, 1e-6)

q_hkick, q_vkick = element_by_name(state.model, "Q_KICK").transport_map_params
@printf("Q_KICK hkick = %.3e, vkick = %.3e\n", q_hkick, q_vkick)

@assert isapprox(q_vkick, 1e-6)
