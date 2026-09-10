export AbstractMesh1D
export OneSideClusteredMesh, TwoSideClusteredMesh
export MeshData

abstract type AbstractMesh1D end

"""
    OneSideClusteredMesh{T}
One-side clustered mesh
---
    OneSideClusteredMesh{T}(
    n_intervals::Integer,
    beta::Real;
    n_uniform_left::Integer = 0,
    n_uniform_right::Integer = 0)
Create a OneSideClusteredMesh with stretching coefficient `beta`, optionally `n_uniform_left` 
and/or `n_uniform_right` uniform cells can be added to the left and right sides of the mesh,
respectively.
The mesh is more clustered as the coefficient `beta` increases.
"""
struct OneSideClusteredMesh{T<:Real} <: AbstractMesh1D
    "Total number of intervals"
    n_intervals::Int
    "Number of uniform intervals on the left"
    n_uniform_left::Int
    "Number of uniform intervals on the right"
    n_uniform_right::Int
    "Stretching parameter"
    beta::Float64
    "Node locations of length `n_intervals` + 1"
    x_node::Vector{T}
    "Center locations of length `n_intervals` + 2 with the first and last points being 0 and 1"
    x_center::Vector{T}
end

"""
    TwoSideClusteredMesh{T}
Two-side clustered mesh
---
    TwoSideClusteredMesh{T}(
    n_intervals::Integer,
    beta::Real,
    gamma::Real = 1.0)
Create a TwoSideClusteredMesh with stretching coefficient `beta`, the ratio of the last to first
intervals being `gamma`
"""
struct TwoSideClusteredMesh{T<:Real} <: AbstractMesh1D
    "Total number of intervals"
    n_intervals::Int
    "Stretching parameter"
    beta::Float64
    "Ratio between the sizes of the last and first cell"
    gamma::Float64
    "Node locations of length `n_intervals` + 1"
    x_node::Vector{T}
    "Center locations of length `n_intervals` + 2 with the first and last points being 0 and 1"
    x_center::Vector{T}
end

function convert_node_mesh_to_centers(x_node::Vector{Float64})
    x_center = Vector{eltype(x_node)}(undef, size(x_node, 1) + 1)
    x_center[2:end-1] .= (x_node[1:end-1] .+ x_node[2:end]) ./ 2
    x_center[1] = 0
    x_center[end] = 1
    return x_center
end

"""
Generate a mesh hyperbolically stretched on two ends, return (`n_intervals`+1) node locations from 0 to 1
# `beta` Stretching coefficient, larger -> more clustered
# `gamma` Ratio between the sizes of the last and first cell
"""
function generate_cluster_tanh_two_sided(n_intervals::Int64, beta::Float64, gamma::Float64)
    sgamma = sqrt(gamma)
    xi = Array(range(-0.5, 0.5, n_intervals + 1))
    y = tanh.(beta .* xi) ./ 2 ./ tanh.(beta / 2) .+ 0.5
    z = y ./ ((1 - sgamma) .* y .+ sgamma)
    z[1] = 0
    z[end] = 1
    return z
end

"""
Generate a mesh hyperbolically stretched on one end, return (`n_intervals`+1) node locations from 0 to 1
# `beta` Stretching coefficient, larger -> more clustered
"""
function generate_cluster_tanh_one_sided(n_intervals::Int64, beta::Float64)
    xi = Array(range(0, 1, n_intervals + 1))
    z = 1 .+ tanh.(beta .* (xi .- 1)) ./ tanh(beta)
    return z
end

function TwoSideClusteredMesh{T}(
    n_intervals::Integer,
    beta::Real,
    gamma::Real=1.0,
) where {T<:Real}
    x_node = generate_cluster_tanh_two_sided(Int64(n_intervals), beta, gamma)
    x_center = convert_node_mesh_to_centers(x_node)
    return TwoSideClusteredMesh{T}(
        n_intervals,
        beta,
        gamma,
        Vector{T}(x_node),
        Vector{T}(x_center),
    )
end

function OneSideClusteredMesh{T}(
    n_intervals::Integer,
    beta::Real;
    n_uniform_left::Integer=0,
    n_uniform_right::Integer=0
) where {T<:Real}
    (n_uniform_left + n_uniform_right < n_intervals) ? nothing :
    error(
        "Number of uniform cells $n_uniform_left + $n_uniform_right more than total cells $n_intervals",
    )
    x_stretched_node = generate_cluster_tanh_one_sided(
        Int64(n_intervals - n_uniform_left - n_uniform_right),
        beta,
    )
    x_left_node = [
        x_stretched_node[1] + (x_stretched_node[2] - x_stretched_node[1]) * i for
        i in -n_uniform_left:-1
    ]
    x_right_node = [
        x_stretched_node[end] + (x_stretched_node[end] - x_stretched_node[end-1]) * i
        for i in 1:n_uniform_right
    ]
    x = vcat(x_left_node, x_stretched_node, x_right_node)
    x_node = (x .- x[1]) ./ (x[end] - x[1])
    x_center = convert_node_mesh_to_centers(x_node)
    return OneSideClusteredMesh{T}(
        n_intervals,
        beta,
        Vector{T}(x_node),
        Vector{T}(x_center),
    )
end

struct MeshData{T<:Real}
    zz::Vector{T}
    zw::Vector{T}
    dz::Vector{T}
    dzw::Vector{T}
    kx0::T
    ky0::T
    Lz::T
end

"Return cell center, cell nodes, cell spacings"
function MeshData(global_mesh::AbstractMesh1D, kx0, ky0, Lz)
    global_zz = global_mesh.x_center
    global_zw = similar(global_zz)
    global_zw[1:end-1] .= global_mesh.x_node[1:end]
    global_zw[end] = global_zw[end-1]
    global_dz = similar(global_zz)
    global_dz[1:end-1] .= global_zz[2:end] .- global_zz[1:end-1]
    global_dz[end] = 0
    global_dzw = similar(global_zz)
    global_dzw[1:end-1] .= global_zw[2:end] .- global_zw[1:end-1]
    global_dzw[end] = 0

    T = eltype(global_zz)
    return MeshData{T}(
        global_zz,
        global_zw,
        global_dz,
        global_dzw,
        kx0,
        ky0,
        Lz
    )
end