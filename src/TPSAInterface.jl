module TPSAInterface
using LinearAlgebra

# Traits
abstract type TPSBehavior end
struct IsTPS <: TPSBehavior end
struct IsNotTPS <: TPSBehavior end

abstract type TPSTypeBehavior end
struct IsTPSType <: TPSTypeBehavior end
struct IsNotTPSType <: TPSTypeBehavior end

"Returns `TPSAInterface.IsTPS() if `t` is a TPS, else `TPSAInterface.IsNotTPS()."
is_tps(t) = IsNotTPS() 
"Returns `TPSAInterface.IsTPSType() if `t` is a TPS type, else `TPSAInterface.IsNotTPSType()."
is_tps_type(t) = IsNotTPSType()

"Returns the number type of the monomial coefficients of the TPS. `t` may be a TPS or TPS `Type`."
numtype(t::Number) = typeof(t)
numtype(::Type{T}) where {T<:Number} = T

"Returns the number of variables in the TPSA. `t` may be a TPS or TPS `Type`."
nvars(t) = 0 
"Returns the number of parameters in the TPSA. `t` may be a TPS or TPS `Type`."
nparams(t) = 0 
"Returns the number of variables + parameters in the TPSA. `t` may be a TPS or TPS `Type`."
ndiffs(t) = 0 
"Returns the maximum truncation order of the TPSA. `t` may be a TPS or TPS `Type`."
maxord(t) = 0 
"Returns the number of monomial coefficients in the TPSA. `t` may be a TPS or TPS `Type`."
numcoefs(t) = 1 

"Returns the zeroth order (scalar) part of the TPS `t`."
scalar(t::Number) = t

"Returns a norm of the passed TPS `t` including all monomial coefficients."
norm_tps(t::Number) = norm(t) 

"Sets all monomial coefficients in `t` equal to 0."
clear!(t) = error("Not implemented!") 

"""
    geti(t, i::Integer)

Gets the `i`-th monomial coefficient in the TPS `t`, where `i=0` specifies the 
scalar part, and the monomials are sorted by order.
"""
geti(t, i::Integer) = error("Not implemented!") 

"""
    getm(t, mono::AbstractArray{<:Integer}) 

Gets the coefficient of the monomial with orders `mono`.
"""
getm(t, mono::AbstractArray{<:Integer}) = error("Not implemented!")

"""
    seti!(t, v, i::Integer)

Sets the `i`-th monomial coefficient in the TPS `t` to `v`, where `i=0` specifies the 
scalar part, and the monomials are sorted by order.
"""
seti!(t, v, i::Integer) = error("Not implemented!")

"""
    setm!(t, v, mono::AbstractArray{<:Integer}) 

Sets the coefficient of the monomial with orders `mono` to `v`.
"""
setm!(t, v, mono::AbstractArray{<:Integer}) = error("Not implemented!")

"Sets the entire TPS `t` equal to `t1`, where `t1` may be another TPS or a `Number`."
copy!(t, t1) = error("Not implemented!")

# Arithmetic operators
"Sets the TPS `t` equal to `a + b`"
add!(t, a, b) = error("Not implemented!")
"Sets the TPS `t` equal to `a - b`"
sub!(t, a, b) = error("Not implemented!")
"Sets the TPS `t` equal to `a * b`"
mul!(t, a, b) = error("Not implemented!")
"Sets the TPS `t` equal to `a / b`"
div!(t, a, b) = error("Not implemented!")
"Sets the TPS `t` equal to `a ^ b`"
pow!(t, a, b) = error("Not implemented!")

"Sets `t` equal to the homogenous polynomial of order `ord` in `t1`."
getord!(t, t1, ord) = error("Not implemented!")

"Sets `t` equal to `t1` with monomials at order `ord` and above removed. Or, if `ord` 
is negative, will remove monomials with orders at and below `abs(ord)`."
cutord!(t, t1, ord) = error("Not implemented!") 

# Defaults for out of place:
getord(t1, ord) = (t = zero(t1); getord!(t, t1, ord); return t)
cutord(t1, ord) = (t = zero(t1); cutord!(t, t1, ord); return t)

# Derivative wrt the i-th differential
"Sets `t` equal to the derivative of `t1` with respect to the `i`-th differential."
deriv!(t, t1, i) = error("Not implemented!")

# Inverts the map m1 and puts the result into m
"Sets the map `m` equal to the inverse of the map `m1`, ignoring any scalar part."
inv!(m, m1) = error("Not implemented!")

# Composes the TPS or vector of TPSs m2 ∘ m1 and puts the result into m
"Sets `m` equal to `m2 ∘ m1` where `m1` and/or `m2` may be TPS types or arrays of TPS types."
compose!(m, m2, m1) = error("Not implemented!")

# TO-DO: implement defaults of these (maybe in NNF)
fgrad!(g, F, h) = error("Not implemented!") # computes F dot grad h
liebra!(G, F, H) = error("Not implemented!") # Computes Lie bracket for Hamiltonian vector fields

end
