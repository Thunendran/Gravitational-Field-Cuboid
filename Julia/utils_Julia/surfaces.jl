# utils_Julia/surfaces.jl

include("safe_math.jl")

# ============================================
# COMMON RADIUS + LOG + ATAN UTILS
# ============================================

@inline function _R(a::Real, b::Real, c::Real)
    A = Float64(a)
    B = Float64(b)
    C = Float64(c)
    return sqrt(A*A + B*B + C*C)
end

@inline function _log_term(a::Real, b::Real, c::Real, X::Real)
    R = _R(a, b, c)
    return sp_log(Float64(X) + R), R
end

@inline function _atan_term(a::Real, b::Real, c::Real, X::Real)
    R = _R(a, b, c)
    return sp_arctan(Float64(b) * Float64(c), Float64(X) * R), R
end


# ============================================
# FACES — 1:1 with Python
# ============================================

function Front(X1, X2, X3, X4, X5, X6)
    X1 = Float64(X1); X2 = Float64(X2); X3 = Float64(X3)
    X4 = Float64(X4); X5 = Float64(X5); X6 = Float64(X6)

    # Front face uses terms involving (-X2 + R).
    log1, _ = _log_term(X2, X4, X6, -X2)
    log2, _ = _log_term(X2, X4, X5, -X2)
    log3, _ = _log_term(X2, X3, X5, -X2)
    log4, _ = _log_term(X2, X3, X6, -X2)

    atan1, _ = _atan_term(X2, X4, X6, X2)
    atan2, _ = _atan_term(X2, X4, X5, X2)
    atan3, _ = _atan_term(X2, X3, X5, X2)
    atan4, _ = _atan_term(X2, X3, X6, X2)

    F_L = -(
        X6 * X4 * log1 +
        X4 * X5 * log2 +
        X5 * X3 * log3 +
        X3 * X6 * log4
    )

    F_T = -(X2 * X2 / 2.0) * (
        atan1 + atan2 + atan3 + atan4
    )

    return F_L + F_T
end


function Back(X1, X2, X3, X4, X5, X6)
    X1 = Float64(X1); X2 = Float64(X2); X3 = Float64(X3)
    X4 = Float64(X4); X5 = Float64(X5); X6 = Float64(X6)

    # logs
    log1, _ = _log_term(X1, X4, X5, X1)
    log2, _ = _log_term(X1, X4, X6, X1)
    log3, _ = _log_term(X1, X3, X5, X1)
    log4, _ = _log_term(X1, X3, X6, X1)

    # atans
    atan1, _ = _atan_term(X1, X4, X5, X1)
    atan2, _ = _atan_term(X1, X4, X6, X1)
    atan3, _ = _atan_term(X1, X3, X5, X1)
    atan4, _ = _atan_term(X1, X3, X6, X1)

    B_L = (
        X4 * X5 * log1 +
        X6 * X4 * log2 +
        X5 * X3 * log3 +
        X3 * X6 * log4
    )

    B_T = -(X1 * X1 / 2.0) * (
        atan1 + atan2 + atan3 + atan4
    )

    return B_L + B_T
end


function Right(X1, X2, X3, X4, X5, X6)
    X1 = Float64(X1); X2 = Float64(X2); X3 = Float64(X3)
    X4 = Float64(X4); X5 = Float64(X5); X6 = Float64(X6)

    log1, _ = _log_term(X1, X3, X6, X3)
    log2, _ = _log_term(X1, X3, X5, X3)
    log3, _ = _log_term(X2, X3, X5, X3)
    log4, _ = _log_term(X2, X3, X6, X3)

    atan1, _ = _atan_term(X3, X1, X6, X3)
    atan2, _ = _atan_term(X3, X1, X5, X3)
    atan3, _ = _atan_term(X3, X2, X5, X3)
    atan4, _ = _atan_term(X3, X2, X6, X3)

    R_L = (
        X6 * X1 * log1 +
        X1 * X5 * log2 +
        X5 * X2 * log3 +
        X2 * X6 * log4
    )

    R_T = -(X3 * X3 / 2.0) * (
        atan1 + atan2 + atan3 + atan4
    )

    return R_L + R_T
end


function Left(X1, X2, X3, X4, X5, X6)
    X1 = Float64(X1); X2 = Float64(X2); X3 = Float64(X3)
    X4 = Float64(X4); X5 = Float64(X5); X6 = Float64(X6)

    log1, _ = _log_term(X1, X4, X6, -X4)
    log2, _ = _log_term(X1, X4, X5, -X4)
    log3, _ = _log_term(X2, X4, X5, -X4)
    log4, _ = _log_term(X2, X4, X6, -X4)

    atan1, _ = _atan_term(X4, X1, X6, X4)
    atan2, _ = _atan_term(X4, X1, X5, X4)
    atan3, _ = _atan_term(X4, X2, X5, X4)
    atan4, _ = _atan_term(X4, X2, X6, X4)

    L_L = -(
        X6 * X1 * log1 +
        X1 * X5 * log2 +
        X5 * X2 * log3 +
        X2 * X6 * log4
    )

    L_T = -(X4 * X4 / 2.0) * (
        atan1 + atan2 + atan3 + atan4
    )

    return L_L + L_T
end


function Top(X1, X2, X3, X4, X5, X6)
    X1 = Float64(X1); X2 = Float64(X2); X3 = Float64(X3)
    X4 = Float64(X4); X5 = Float64(X5); X6 = Float64(X6)

    log1, _ = _log_term(X1, X4, X6, -X6)
    log2, _ = _log_term(X2, X4, X6, -X6)
    log3, _ = _log_term(X2, X3, X6, -X6)
    log4, _ = _log_term(X1, X3, X6, -X6)

    atan1, _ = _atan_term(X6, X1, X4, X6)
    atan2, _ = _atan_term(X6, X4, X2, X6)
    atan3, _ = _atan_term(X6, X2, X3, X6)
    atan4, _ = _atan_term(X6, X3, X1, X6)

    T_L = -(
        X1 * X4 * log1 +
        X4 * X2 * log2 +
        X2 * X3 * log3 +
        X3 * X1 * log4
    )

    T_T = -(X6 * X6 / 2.0) * (
        atan1 + atan2 + atan3 + atan4
    )

    return T_L + T_T
end


function Bottom(X1, X2, X3, X4, X5, X6)
    X1 = Float64(X1); X2 = Float64(X2); X3 = Float64(X3)
    X4 = Float64(X4); X5 = Float64(X5); X6 = Float64(X6)

    log1, _ = _log_term(X1, X4, X5, X5)
    log2, _ = _log_term(X2, X4, X5, X5)
    log3, _ = _log_term(X2, X3, X5, X5)
    log4, _ = _log_term(X1, X3, X5, X5)

    atan1, _ = _atan_term(X5, X1, X4, X5)
    atan2, _ = _atan_term(X5, X2, X4, X5)
    atan3, _ = _atan_term(X5, X2, X3, X5)
    atan4, _ = _atan_term(X5, X3, X1, X5)

    B_L = (
        X1 * X4 * log1 +
        X4 * X2 * log2 +
        X2 * X3 * log3 +
        X3 * X1 * log4
    )

    B_T = -(X5 * X5 / 2.0) * (
        atan1 + atan2 + atan3 + atan4
    )

    return B_L + B_T
end
