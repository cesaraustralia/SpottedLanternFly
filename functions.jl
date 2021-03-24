# define Briere struct
struct BriereModel{A, B, TMN, TMX} <: RateModel
    a::A
    b::B
    T_min::TMN
    T_max::TMX
end

@inline function rate(m::BriereModel, x)
    (x < m.T_min) |( x > m.T_max) ? zero(x) : 
        @fastmath m.a * x * (x - m.T_min) * (m.T_max - x)^(1 / m.b)
end

struct LinearModel{A, B} <: RateModel
    a::A
    b::B
end

@inline function rate(m::LinearModel, x)
      y = @fastmath m.a + m.b * x
      y < 0 ? 0 : y
end


struct DegreeDayModel{LC, DD} <: RateModel
    lc::LC
    dd::DD
end

@inline function rate(m::DegreeDayModel, x)
    x = x < m.lc ? m.lc : x
    d = @fastmath  x - m.lc
    d/m.dd
end


struct MultiRate{T<:Tuple} <: RateModel 
    models::T
end

@inline GrowthMaps.rate(m::MultiRate, x) = SVector(map(m -> rate(m, x), m.models))


# Hazard models
struct HazModel{A, B, C, D, E} <: RateModel
    a::A
    b::B
    c::C
    d::D
    e::E
end

@inline function rate(m::HazModel, t)
    exposure = 1 # units 1 day
    # t = t <  5 ?  5 : t
    # t = t > 40 ? 40 : t 
    exp(m.a + m.b*t + m.c*t^2 + m.d*t^3 + m.e*t^4 + log(exposure))
end

struct EggHazModel{A, B, C, D} <: RateModel
    a::A
    b::B
    c::C
    d::D
end

@inline function rate(m::EggHazModel, t)
    exposure = 1 # units 1 day 
    # h = t < -5 ? exp(m.a + m.b*t + log(exposure)) : 0
    # t > 12 ? exp(m.c + m.d*t + log(exposure)) : h
    h = t < -5 ? exp(m.a + m.b*t + log(exposure)) : exp(m.c + m.d*t + log(exposure))
    h > 0 ? h : 0

end