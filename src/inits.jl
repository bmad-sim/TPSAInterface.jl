# TPSA Backends

"""
    abstract type AbstractTPSAInit

Abstract type for Taylor series initialization. 
This is analogous to the AutoGTPSA, AutoForwardDiff, etc types in ADTypes.jl. 
The reason why a TPS type may not be suitable for initialization is because the dimensionality 
may be stored at runtime as opposed to statically (e.g. TPS{GTPSA.Dynamic})
"""
abstract type AbstractTPSAInit end

#

"""
    InitGTPSA{D,DD}

Struct used to define a TPSA using the [GTPSA.jl](https://github.com/bmad-sim/GTPSA.jl) backend.

Defined by [`TPSAInterface.jl`](https://github.com/bmad-sim/TPSAInterface.jl).

# Constructors

    InitGTPSA{D,DD}(; dynamic_descriptor=nothing)

## Fields
  - For static `Descriptor` resolution:
    + `D` is the `Descriptor`
    + `DD` is `Nothing` and `dynamic_descriptor` is `nothing`

  - For dynamic `Descriptor` resolution:
    + `D == GTPSA.Dynamic`
    + `DD == Descriptor` and `dynamic_descriptor` is said `Descriptor`
"""
Base.@kwdef struct InitGTPSA{D,DD} <: AbstractTPSAInit
    dynamic_descriptor::DD = nothing
end