# Chapter 11, Exercise 2: Move a bend while preserving lattice length.

using SciBmad
using Printf

d = Drift(name="D", L=1.0)
b = SBend(name="B", L=1.0, angle=0.10)
q = Quadrupole(name="Q", L=1.0, Kn1=0.20)

lat = Beamline(
    [d, b, q];
    species_ref=Species("muon"),
    pc_ref=1e9,
)

function bend_positions(lat)
    d_ele, b_ele, q_ele = lat.line
    return (
        bend_upstream=b_ele.s,
        bend_downstream=q_ele.s,
        lattice_length=q_ele.s + q_ele.L,
        drift_length=d_ele.L,
        bend_length=b_ele.L,
    )
end

move_controller = Controller(
    (d, :L) => (ele; dL_drift, dL_bend) -> ele.L + dL_drift,
    (b, :L) => (ele; dL_drift, dL_bend) -> ele.L + dL_bend;
    vars=(; dL_drift=0.0, dL_bend=0.0),
)

old_shift = Ref(0.0)

function set_bend_upstream_shift!(controller, shift, old_shift)
    step = shift - old_shift[]
    set!(controller; dL_drift=step, dL_bend=-step)
    old_shift[] = shift
end

before = bend_positions(lat)
set_bend_upstream_shift!(move_controller, 0.20, old_shift)
after = bend_positions(lat)

println("Before move: ", before)
println("After move:  ", after)

@assert isapprox(after.bend_upstream, before.bend_upstream + 0.20)
@assert isapprox(after.bend_downstream, before.bend_downstream)
@assert isapprox(after.lattice_length, before.lattice_length)
@assert isapprox(after.drift_length, before.drift_length + 0.20)
@assert isapprox(after.bend_length, before.bend_length - 0.20)

@printf(
    "B upstream moved by %.3f m; B downstream and lattice end stayed fixed.\n",
    after.bend_upstream - before.bend_upstream,
)
