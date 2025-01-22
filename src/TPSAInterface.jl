module TPSAInterface

# Traits
abstract type TPSBehavior end
struct IsTPS <: TPSBehavior end
struct IsNotTPS <: TPSBehavior end

abstract type TPSTypeBehavior end
struct IsTPSType <: TPSTypeBehavior end
struct IsNotTPSType <: TPSTypeBehavior end

is_tps(t) = IsNotTPS() 
is_tps_type(t) = IsNotTPSType()

# Returns number type of monomial coefficients
numtype(t::Number) = typeof(t)
numtype(::Type{T}) where {T<:Number} = T

nvars(t) = 0 # Returns number of variables in TPSA
nparams(t) = 0 # Returns number of parameters in TPSA
ndiffs(t) = 0 # Returns number of variables + number of parameters in TPSA
maxord(t) = 0 # Returns the maximum truncation order of the TPSA
numcoefs(t) = 1 # Returns the number of coefficients in the TPSA (1 for regular numbers)

scalar(t::Number) = t # Returns the zeroth order (scalar) part of a passed TPS
norm_tps(t::Number) = norm(t) # Returns a TPS norm including all monomial coefficients

clear!(t) = error("Not implemented!") # Clears the TPS

# Get/set monomials in the TPS t
geti(t, i::Integer) = error("Not implemented!") # 0-based indexing for TPS so that i=0 is scalar part
getm(t, mono::AbstractArray{<:Integer}) = error("Not implemented!") # Monomial as array of orders
seti!(t, v, i::Integer) = error("Not implemented!")
setm!(t, v, mono::AbstractArray{<:Integer}) = error("Not implemented!")

# Sets the entire TPS t equal to t1 (where t1 is a TPS or a Number)
copy!(t, t1) = error("Not implemented!")

# Arithmetic operators
add!(t, a, b) = error("Not implemented!")
sub!(t, a, b) = error("Not implemented!")
mul!(t, a, b) = error("Not implemented!")
div!(t, a, b) = error("Not implemented!")
pow!(t, a, b) = error("Not implemented!")

getord!(t, t1, order) = error("Not implemented!") # Get order as TPS
cutord!(t, t1, order) = error("Not implemented!") # Cuts order from TPS

# Defaults for out of place:
getord(t1, order) = (t = zero(t1); getord!(t, t1, order); return t)
cutord(t1, order) = (t = zero(t1); cutord!(t, t1, order); return t)

# Derivative wrt the i-th differential
deriv!(t, t1, i) = error("Not implemented!")

# Inverts the map m1 and puts the result into m
inv!(m, m1) = error("Not implemented!")

# Composes the TPS or vector of TPSs m2 ∘ m1 and puts the result into m
compose!(m, m2, m1) = error("Not implemented!")

end
