module AlpsTools

export FlowField

include("mesh/mesh.jl")

struct FlowField{T<:Real}
    u::AbstractArray{T, 3}
    v::AbstractArray{T, 3}
    w::AbstractArray{T, 3}
end

FlowField{T}(dims::Dims{3}) where {T} = FlowField(zeros(T, dims), zeros(T, dims), zeros(T, dims))

end # module scripts
