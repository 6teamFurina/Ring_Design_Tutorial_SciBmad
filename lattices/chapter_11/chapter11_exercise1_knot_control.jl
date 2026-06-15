# Chapter 11, Exercise 1: Control from knot points.

using SciBmad
using Printf

struct NaturalCubicSpline
    x::Vector{Float64}
    y::Vector{Float64}
    second_derivative::Vector{Float64}
end

function NaturalCubicSpline(x, y)
    xs = Float64.(x)
    ys = Float64.(y)
    length(xs) == length(ys) || error("x and y must have the same length.")
    length(xs) >= 3 || error("At least three knot points are required.")
    all(diff(xs) .> 0) || error("Knot positions must be strictly increasing.")

    n = length(xs)
    h = diff(xs)
    matrix = zeros(n - 2, n - 2)
    rhs = zeros(n - 2)

    for i in 1:(n - 2)
        matrix[i, i] = 2 * (h[i] + h[i + 1])
        i > 1 && (matrix[i, i - 1] = h[i])
        i < n - 2 && (matrix[i, i + 1] = h[i + 1])
        rhs[i] = 6 * ((ys[i + 2] - ys[i + 1]) / h[i + 1] -
                      (ys[i + 1] - ys[i]) / h[i])
    end

    second_derivative = zeros(n)
    second_derivative[2:(n - 1)] = matrix \ rhs
    return NaturalCubicSpline(xs, ys, second_derivative)
end

function (spline::NaturalCubicSpline)(x)
    spline.x[1] <= x <= spline.x[end] ||
        error("Control value $x is outside [$(spline.x[1]), $(spline.x[end])].")

    i = min(searchsortedlast(spline.x, x), length(spline.x) - 1)
    h = spline.x[i + 1] - spline.x[i]
    left = (spline.x[i + 1] - x) / h
    right = (x - spline.x[i]) / h

    return left * spline.y[i] + right * spline.y[i + 1] +
           ((left^3 - left) * spline.second_derivative[i] +
            (right^3 - right) * spline.second_derivative[i + 1]) * h^2 / 6
end

hh_knots = [-0.04, -0.02, 0.0, 0.02, 0.04]
Kn1_spline = NaturalCubicSpline(hh_knots, 0.7 .* hh_knots)
x_offset_spline = NaturalCubicSpline(hh_knots, 0.1 .* hh_knots)

q = Quadrupole(name="Q", L=1.0)
knot_controller = Controller(
    (q, :Kn1) => (ele; hh) -> Kn1_spline(hh),
    (q, :x_offset) => (ele; hh) -> x_offset_spline(hh);
    vars=(; hh=0.0,),
)

println("Knot-point controller checks:")
for hh in [-0.035, -0.011, 0.013, 0.037]
    knot_controller.hh = hh
    @printf(
        "  hh = % .3f  Kn1 = % .6f  x_offset = % .6f m\n",
        hh,
        q.Kn1,
        q.x_offset,
    )
    @assert isapprox(q.Kn1, 0.7 * hh; atol=1e-14)
    @assert isapprox(q.x_offset, 0.1 * hh; atol=1e-14)
end

outside_interval_errors = try
    knot_controller.hh = 0.05
    false
catch
    true
end
@assert outside_interval_errors

