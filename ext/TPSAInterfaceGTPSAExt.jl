module TPSAInterfaceGTPSAExt
import TPSAInterface as TI
using GTPSA: GTPSA, TPS

TI.is_tps(::TPS) = TI.IsTPS()
TI.is_tps_type(::Type{<:TPS}) = TI.IsTPSType()

TI.numtype(::TPS{T}) where {T} = T
TI.numtype(::Type{<:TPS{T}}) where {T} = T

TI.nvars(t::TPS) = GTPSA.numvars(t)
TI.nparams(t::TPS) = GTPSA.numparms(t)
TI.ndiffs(t::TPS) = GTPSA.numnn(t)
TI.maxord(t::TPS) = unsafe_load(GTPSA.getdesc(t).desc).mo
TI.numcoefs(t::TPS) = GTPSA.numcoefs(t)

# Static Descriptor resolution:
TI.nvars(::Type{TPS{T,D}}) where {T,D} = GTPSA.numvars(D)
TI.nparams(::Type{TPS{T,D}}) where {T,D} = GTPSA.numparams(D)
TI.ndiffs(::Type{TPS{T,D}}) where {T,D} = GTPSA.numnn(D)
TI.maxord(::Type{TPS{T,D}}) where {T,D} = unsafe_load(D.desc).mo
TI.numcoefs(::Type{TPS{T,D}}) where {T,D} = unsafe_load(unsafe_load(D.desc).ord2idx, maxord(D))

# Dynamic Descriptor resolution:
TI.nvars(::Type{TPS{T,GTPSA.Dynamic}}) where {T} = TI.nvars(TPS{T,GTPSA.desc_current})
TI.nparams(::Type{TPS{T,GTPSA.Dynamic}}) where {T} = TI.nparams(TPS{T,GTPSA.desc_current})
TI.ndiffs(::Type{TPS{T,GTPSA.Dynamic}}) where {T} = TI.ndiffs(TPS{T,GTPSA.desc_current})
TI.maxord(::Type{TPS{T,GTPSA.Dynamic}}) where {T} = TI.maxord(TPS{T,GTPSA.desc_current})
TI.numcoefs(::Type{TPS{T,GTPSA.Dynamic}}) where {T} = TI.numcoefs(TPS{T,GTPSA.desc_current})

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
TI.liebra!(
  G::AbstractArray{TPS{T,DG}}, 
  F::AbstractArray{TPS{T,DF}}, 
  H::AbstractArray{TPS{T,DH}}
) where {T,DG,DF,DH}
  @assert !(G === F) && !(G === H) "Aliasing any source arguments with the destination in lb! is not allowed"
  return GTPSA.liebra!(numvars(F), F, H, G)
end

end