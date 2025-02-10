module TPSAInterfaceGTPSAExt
import TPSAInterface as TI
using TPSAInterface: DefGTPSA
using GTPSA: GTPSA, Descriptor, TPS

# =================================== #
# Static Descriptor Resolution:
TI.init_tps(::Type{T}, ::DefGTPSA{D,Nothing}) where {T,D} = TPS{T,D}()
TI.init_tps_type(::Type{T}, ::DefGTPSA{D,Nothing}) where {T,D} = TPS{T,D}
TI.getdef(::TPS{T,D}) where {T,D} = DefGTPSA{D,Nothing}(; dynamic_descriptor=nothing)
TI.getdef(::Type{TPS{T,D}}) where {T,D} = DefGTPSA{D,Nothing}(; dynamic_descriptor=nothing)

TI.nvars(::DefGTPSA{D,Nothing}) where {D}   = Int(GTPSA.numvars(D))
TI.nparams(::DefGTPSA{D,Nothing}) where {D} = Int(GTPSA.numparams(D))
TI.ndiffs(::DefGTPSA{D,Nothing}) where {D}  = Int(GTPSA.numnn(D))
TI.maxord(::DefGTPSA{D,Nothing}) where {D}  = Int(unsafe_load(D.desc).mo)
TI.nmonos(::DefGTPSA{D,Nothing}) where {D}  = Int(unsafe_load(unsafe_load(D.desc).ord2idx, maxord(D)))

# =================================== #
# Dynamic Descriptor Resolution:
function TI.init_tps(::Type{T}, def::DefGTPSA{GTPSA.Dynamic,Descriptor}) where {T}
  return TPS{T,GTPSA.Dynamic}(use=def.dynamic_descriptor)
end
TI.init_tps_type(::Type{T}, ::DefGTPSA{GTPSA.Dynamic,Descriptor}) where {T} = TPS{T,GTPSA.Dynamic}

function TI.getdef(t::TPS{T,GTPSA.Dynamic}) where {T} 
  return DefGTPSA{GTPSA.Dynamic,Descriptor}(; dynamic_descriptor=GTPSA.getdesc(t))
end

function TI.getdef(::Type{TPS{T,GTPSA.Dynamic}}) where {T}
  return DefGTPSA{GTPSA.Dynamic,Descriptor}(; dynamic_descriptor=GTPSA.desc_current)
end

TI.nvars(def::DefGTPSA{GTPSA.Dynamic,Descriptor})   = Int(GTPSA.numvars(def.dynamic_descriptor))
TI.nparams(def::DefGTPSA{GTPSA.Dynamic,Descriptor}) = Int(GTPSA.numparams(def.dynamic_descriptor))
TI.ndiffs(def::DefGTPSA{GTPSA.Dynamic,Descriptor})  = Int(GTPSA.numnn(def.dynamic_descriptor))
TI.maxord(def::DefGTPSA{GTPSA.Dynamic,Descriptor})  = Int(unsafe_load(def.dynamic_descriptor.desc).mo)
TI.nmonos(def::DefGTPSA{GTPSA.Dynamic,Descriptor})  = Int(unsafe_load(unsafe_load(dynamic_descriptor.desc).ord2idx, maxord(def)))
# =================================== #

TI.is_tps(::TPS) = TI.IsTPS()
TI.is_tps_type(::Type{<:TPS}) = TI.IsTPSType()

TI.numtype(::TPS{T}) where {T} = T
TI.numtype(::Type{<:TPS{T}}) where {T} = T

TI.scalar(t::TPS) = GTPSA.scalar(t)

TI.norm_tps(t::TPS) = GTPSA.normTPS(t)

TI.clear!(t::TPS) = GTPSA.clear!(t)

TI.geti(t::TPS, i::Integer) = t[i]
TI.getm(t::TPS, mono::AbstractArray{<:Integer}) = t[mono]
TI.seti!(t::TPS, v, i::Integer) = (t[i] = v; return v)
TI.setm!(t::TPS, v, mono::AbstractArray{<:Integer}) = (t[mono] = v; return v)

TI.copy!(t::TPS, t1) = GTPSA.setTPS!(t, t1, change=true)

TI.add!(t::TPS, a, b) = (GTPSA.add!(t, a, b); return t)
TI.sub!(t::TPS, a, b) = (GTPSA.sub!(t, a, b); return t)
TI.mul!(t::TPS, a, b) = (GTPSA.mul!(t, a, b); return t)
TI.div!(t::TPS, a, b) = (GTPSA.div!(t, a, b); return t)
TI.pow!(t::TPS, a, b) = (GTPSA.pow!(t, a, b); return t)

TI.getord!(t::TPS, t1::TPS, ord) = GTPSA.getord!(t, t1, ord)
TI.cutord!(t::TPS, t1::TPS, ord) = GTPSA.cutord!(t, t1, ord)

TI.deriv!(t::TPS, t1::TPS, i) = GTPSA.deriv!(t, t1, i)

function TI.inv!(
  m::Union{AbstractArray{TPS{T,D}},TPS{T,D}}, 
  m1::Union{AbstractArray{TPS{T,D1}},TPS{T,D1}}
) where {T<:Union{Float64,ComplexF64},D,D1} 
  return GTPSA.inv!(m, m1)
end

function TI.compose!(
  m::Union{AbstractArray{TPS{T,D}},TPS{T,D}}, 
  m2::Union{AbstractArray{TPS{T,D2}},TPS{T,D2}},
  m1::Union{AbstractArray{TPS{T,D1}},TPS{T,D1}}
) where {T<:Union{Float64,ComplexF64},D,D1,D2}
  return GTPSA.compose!(m, m2, m1)
end

TI.fgrad!(g::TPS{T}, F::AbstractArray{TPS{T,D}}, h::TPS{T}) where {T,D} = GTPSA.fgrad!(g, F, h)
function TI.liebra!(
  G::AbstractArray{TPS{T,DG}}, 
  F::AbstractArray{TPS{T,DF}}, 
  H::AbstractArray{TPS{T,DH}}
) where {T,DG,DF,DH}
  @assert !(G === F) && !(G === H) "Aliasing any source arguments with the destination in lb! is not allowed"
  return GTPSA.liebra!(numvars(F), F, H, G)
end

end