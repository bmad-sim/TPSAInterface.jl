module TPSAInterface

using LinearAlgebra
export AbstractTPSAInit,
       InitGTPSA

import Base: copy!

# Traits

"""
    abstract type TPSBehavior end
    struct IsTPS <: TPSBehavior end
    struct IsNotTPS <: TPSBehavior end

Holy traits group.
""" TPSBehavior, IsTPS, IsNotTPS

abstract type TPSBehavior end
struct IsTPS <: TPSBehavior end
struct IsNotTPS <: TPSBehavior end

"""
    abstract type TPSTypeBehavior end
    struct IsTPSType <: TPSTypeBehavior end
    struct IsNotTPSType <: TPSTypeBehavior end

Holy traits group.
""" TPSTypeBehavior, IsTPSType, IsNotTPSType

abstract type TPSTypeBehavior end
struct IsTPSType <: TPSTypeBehavior end
struct IsNotTPSType <: TPSTypeBehavior end

#

include("inits.jl")

"""
    is_tps(t) 

Returns `TPSAInterface.IsTPS() if `t` is a TPS, else `TPSAInterface.IsNotTPS()."
"""
is_tps(t) = IsNotTPS() 

#

"""
    is_tps_type(t)

Returns `TPSAInterface.IsTPSType() if `t` is a TPS type, else `TPSAInterface.IsNotTPSType()."
"""
is_tps_type(t) = IsNotTPSType()

#

"""
    numtype(t::Number)
    numtype(::Type{T}) where {T<:Number}

Returns the number type of the monomial coefficients of the TPS. `t` may be a TPS or TPS `Type`."
""" numtype

numtype(t::Number) = typeof(t)
numtype(::Type{T}) where {T<:Number} = T

# Constructors using TPSA initializer

"""
    init_tps(::Type, ::AbstractTPSAInit)
    init_tps_type(::Type, ::AbstractTPSAInit) 

Dummy functions that need to be overridden.
""" init_tps

init_tps(::Type, ::AbstractTPSAInit) = error("Please specify a TPSA initializer.")
init_tps_type(::Type, ::AbstractTPSAInit) = error("Please specify a TPSA initializer.")

#

"""
    mono(::Type{T}, init::AbstractTPSAInit, idx::Integer) where {T<:Number}
    mono(::Type{T}, init::AbstractTPSAInit, mono::AbstractArray{<:Integer}) where {T<:Number}
    mono(init::AbstractTPSAInit, idx_or_mono)  # Default to Float64
"""

function mono(::Type{T}, init::AbstractTPSAInit, idx::Integer) where {T<:Number}
    t = init_tps(T, init)
    seti!(t, 1, idx)
    return t
end

function mono(::Type{T}, init::AbstractTPSAInit, mono::AbstractArray{<:Integer}) where {T<:Number}
    t = init_tps(T, init)
    setm!(t, 1, mono)
    return t
end

mono(init::AbstractTPSAInit, idx_or_mono) = mono(Float64, init, idx_or_mono) # Default to Float64

#

"""
    getinit(t::Number)
    getinit(::Type{T}) where {T<:Number}
""" getinit

getinit(t::Number) = error("$t is not a TPS")
getinit(::Type{T}) where {T<:Number} = error("$T is not a TPS type!")

#

"""
    ndiffs(t)

Returns the number of variables + parameters in the TPSA. `t` may be a TPS, TPS `Type`, or a `AbstractTPSAInit`.
"""
ndiffs(t) = ndiffs(getinit(t))

#

"""
    maxord(t)

Returns the maximum truncation order of the TPSA. `t` may be a TPS, TPS `Type`, or a `AbstractTPSAInit`.
"""
maxord(t) = maxord(getinit(t))

#

"""
    nmonos(t)

Returns the number of monomial coefficients in the TPSA. `t` may be a TPS, TPS `Type`, or a `AbstractTPSAInit`.
"""
nmonos(t) = nmonos(getinit(t))

#

"""
    scalar(t::Number)

Returns the zeroth order (scalar) part of the TPS `t`.
"""
scalar(t::Number) = t

#

"""
    norm_tps(t::Number)

Returns a norm of the passed TPS `t` including all monomial coefficients."
"""
norm_tps(t::Number) = norm(t) 

#

"""
    clear!(t)

Sets all monomial coefficients in `t` equal to 0.
"""
clear!(t) = error("Not implemented!") 

#

"""
    geti(t, i::Integer)

Gets the `i`-th monomial coefficient in the TPS `t`, where `i=0` specifies the 
scalar part, and the monomials are sorted by order.
"""
geti(t, i::Integer) = error("Not implemented!") 

#

"""
    getm(t, mono) 

Gets the coefficient of the monomial with orders `mono`.
"""
getm(t, mono) = error("Not implemented!")

#

"""
    seti!(t, v, i::Integer)

Sets the `i`-th monomial coefficient in the TPS `t` to `v`, where `i=0` specifies the 
scalar part, and the monomials are sorted by order.
"""
seti!(t, v, i::Integer) = error("Not implemented!")

#

"""
    setm!(t, v, mono) 

Sets the coefficient of the monomial with orders `mono` to `v`.
"""
setm!(t, v, mono) = error("Not implemented!")

#

"""
    copy!(t, t1)

Sets the entire TPS `t` equal to `t1`, where `t1` may be another TPS or a `Number`. Promotion 
is supported; e.g. if `t` has `numtype` `ComplexF64`, and `t1` has `numtype` `Float64`), 
calling `copy!(t, t1)` is allowed
"""
copy!(t, t1) = error("Not implemented!")

# Arithmetic operators

"""
TPSAInterface functions:

    add!(t, a, b)         # Sets the TPS `t` equal to `a - b`
    sub!(t, a, b)         # Sets the TPS `t` equal to `a - b`
    mul!(t, a, b)         # Sets the TPS `t` equal to `a * b`
    div!(t, a, b)         # Sets the TPS `t` equal to `a / b`
    pow!(t, a, b)         # Sets the TPS `t` equal to `a ^ b`
    real!(t, t1)          # Sets the TPS `t` equal to `real(t1)`
    imag!(t, t1)          # Sets the TPS `t` equal to `imag(t1)`
    getord!(t, t1, ord)   # Sets `t` equal to the homogenous polynomial of order `ord` in `t1`.
    getord(t1, ord)       # Defaults for out of place
    cutord!(t, t1, ord)   # Sets `t` equal to `t1` with monomials at order `ord` and above removed. 
                          #   Or, if `ord` is negative, will remove monomials with orders at and below `abs(ord)`.
    cutord(t1, ord)       # Defaults for out of place
    deriv!(t, t1, i)      # Derivative wrt the i-th differential. Sets `t` equal to the derivative
                          #   of `t1` with respect to the `i`-th differential.
    inv!(m, m1)           # Inverts the map m1 and puts the result into m
                          #   Sets the map `m` equal to the inverse of the map `m1`, ignoring any scalar part.
    compose!(m, m2, m1)   # Composes the TPS or vector of TPSs m2 ∘ m1 and puts the result into m
                          #   Sets `m` equal to `m2 ∘ m1` where `m1` and/or `m2` may be TPS types or arrays of TPS types.
    fgrad!(g, F, h)       # computes F dot grad h
    liebra!(G, F, H)      # Computes Lie bracket for Hamiltonian vector fields
""" add!, sub!, div!, pow!, real!, imag!, getord!, cutord!, deriv!, inv!, compose!, liebra!

add!(t, a, b) = error("Not implemented!")
sub!(t, a, b) = error("Not implemented!")
mul!(t, a, b) = error("Not implemented!")
div!(t, a, b) = error("Not implemented!")
pow!(t, a, b) = error("Not implemented!")
real!(t, t1) = copy!(t, real(t1))
imag!(t, t1) = copy!(t, imag(t1))
getord!(t, t1, ord) = error("Not implemented!")
getord(t1, ord) = (t = zero(t1); getord!(t, t1, ord); return t)
cutord!(t, t1, ord) = error("Not implemented!") 
cutord(t1, ord) = (t = zero(t1); cutord!(t, t1, ord); return t)
deriv!(t, t1, i) = error("Not implemented!")
inv!(m, m1) = error("Not implemented!")
compose!(m, m2, m1) = error("Not implemented!")
### TO-DO: implement defaults of these (maybe in NNF) ###
fgrad!(g, F, h) = error("Not implemented!") # computes F dot grad h
liebra!(G, F, H) = error("Not implemented!") # Computes Lie bracket for Hamiltonian vector fields

"""
    cycle!(t, i::Integer; 
        mono::Union{AbstractArray{<:Integer},Nothing}=nothing, 
        val::Union{Ref{<:Number},Nothing}=nothing)

Cycles through the nonzero monomial coefficients in the TPS `t` starting at `i`. Set 
`i` equal to `-1` to start at the first nonzero monomial coefficient. Returns the next 
nonzero monomial coefficient index and will, if provided, set `m_` and `v_` to the monomial 
orders and coefficient of the next nonzero monomial coefficient, respectively. 
"""
cycle!(
    t, 
    i::Integer; 
    mono::Union{AbstractArray{<:Integer},Nothing}=nothing, 
    val::Union{Ref{<:Number},Nothing}=nothing
) = error("Not implemented!")

end
