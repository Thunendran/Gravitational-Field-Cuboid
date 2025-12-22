# utils_Julia/safe_math.jl

# Special logarithm: log(a) if a > 0, else 0
@inline function sp_log(argument::Real)
    a = Float64(argument)
    return a <= 0.0 ? 0.0 : log(a)
end

# Special arctan: atan(num/den) if both nonzero, else 0
@inline function sp_arctan(numerator::Real, denominator::Real)
    num = Float64(numerator)
    den = Float64(denominator)
    return (num == 0.0 || den == 0.0) ? 0.0 : atan(num / den)
end
