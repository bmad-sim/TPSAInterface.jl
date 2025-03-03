module TPSAInterfaceGTPSAExt
import TPSAInterface as TI
using TPSAInterface: InitGTPSA
using GTPSA: GTPSA, Descriptor, TPS

# =================================== #
# Static Descriptor Resolution:
TI.init_tps(::Type{T}, ::InitGTPSA{D,Nothing}) where {T,D} = TPS{T,D}()
TI.init_tps_type(::Type{T}, ::InitGTPSA{D,Nothing}) where {T,D} = TPS{T,D}
TI.getinit(::TPS{T,D}) where {T,D} = InitGTPSA{D,Nothing}(; dynamic_descriptor=nothing)
TI.getinit(::Type{TPS{T,D}}) where {T,D} = InitGTPSA{D,Nothing}(; dynamic_descriptor=nothing)

TI.ndiffs(::InitGTPSA{D,Nothing}) where {D}  = Int(GTPSA.numnn(D))
TI.maxord(::InitGTPSA{D,Nothing}) where {D}  = Int(unsafe_load(D.desc).mo)
TI.nmonos(init::InitGTPSA{D,Nothing}) where {D}  = Int(unsafe_load(unsafe_load(D.desc).ord2idx, TI.maxord(init)))

# =================================== #
# Dynamic Descriptor Resolution:
function TI.init_tps(::Type{T}, init::InitGTPSA{GTPSA.Dynamic,Descriptor}) where {T}
  return TPS{T,GTPSA.Dynamic}(use=init.dynamic_descriptor)
end
TI.init_tps_type(::Type{T}, ::InitGTPSA{GTPSA.Dynamic,Descriptor}) where {T} = TPS{T,GTPSA.Dynamic}

function TI.getinit(t::TPS{T,GTPSA.Dynamic}) where {T} 
  return InitGTPSA{GTPSA.Dynamic,Descriptor}(; dynamic_descriptor=GTPSA.getdesc(t))
end

function TI.getinit(::Type{TPS{T,GTPSA.Dynamic}}) where {T}
  return InitGTPSA{GTPSA.Dynamic,Descriptor}(; dynamic_descriptor=GTPSA.desc_current)
end

TI.ndiffs(init::InitGTPSA{GTPSA.Dynamic,Descriptor})  = Int(GTPSA.numnn(init.dynamic_descriptor))
TI.maxord(init::InitGTPSA{GTPSA.Dynamic,Descriptor})  = Int(unsafe_load(init.dynamic_descriptor.desc).mo)
TI.nmonos(init::InitGTPSA{GTPSA.Dynamic,Descriptor})  = Int(unsafe_load(unsafe_load(dynamic_descriptor.desc).ord2idx, TI.maxord(init)))
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

TI.real!(t::TPS, t1) = (GTPSA.real!(t, t1); return t)
TI.imag!(t::TPS, t1) = (GTPSA.imag!(t, t1); return t)

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
  return GTPSA.liebra!(length(G), F, H, G)
end

function TI.cycle!(
  t::TPS, 
  i::Integer; 
  mono::Union{AbstractArray{<:Integer},Nothing}=nothing, 
  val::Union{Ref{<:Number},Nothing}=nothing
)
  if !isnothing(mono) 
    if eltype(mono) == UInt8
      m_ = mono
    else
      m_ = UInt8.(mono)
    end
    n = length(m_)
  else
    m_ = C_NULL
    n = 0
  end

  if !isnothing(val)
    if eltype(val) == TI.numtype(t)
      v_ = val
    else
      v_ = Ref{TI.numtype(t)}(val)
    end
  else
    v_ = C_NULL
  end

  out_i = GTPSA.cycle!(t, i, n, m_, v_)

  if !isnothing(val)
    val[] = v_[]
  end
  
  if !isnothing(mono)
    mono .= m_
  end

  return out_i
end

end