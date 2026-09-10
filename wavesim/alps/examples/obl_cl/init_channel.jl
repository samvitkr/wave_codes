using ArgParse
using TOML
using HDF5
using AlpsTools
using Logging
import Random
using FFTW
using LinearAlgebra

function set_commandline_opts()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--config", "-c"
        help = "Specify a configuration file"
        arg_type = String
        default = "param.toml"
        "--section"
        help = "Specify which section to read the configuration from"
        arg_type = String
        default = "default_solver"
        "--grid"
        help = "Specify a grid input file"
        arg_type = String
        default = "grid.h5"
        "--output", "-o"
        help = "Specify the output file"
        arg_type = String
        default = "restart.h5"
        "--overwrite"
        help = "Overwrite existing files"
        action = :store_true
    end

    return s
end

"Parabolic profile"
function add_parabolic_profile!(flow::FlowField{T}, mesh::MeshData{T}; B=10.0) where {T}
    z = mesh.zz
    for k in axes(z, 1)
      flow.u[:, :, k] .+= B * z[k]^2
    end
end

"Simple log mean profile"
log_mean_profile(z; B=5.0, kappa=0.4) = log.(z) ./ kappa .+ B

"""
    A composite mean profile from the visous layer to the log layer
---
See the equation listed in https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/introducing-jfm-notebooks/2C28AB44477FE45A913441D8F9AF29CF"
"""
function composite_mean_profile(z; B=4.2, kappa=0.38)
    c1 = 12.0
    c2 = (1 / kappa) * log(c1) + B
    return ((1 / kappa) .* log.(c1 .+ z) .+ B) ./ sqrt.(1 .+ (z / c2) .^ (-2))
end

function add_log_mean_profile!(flow::FlowField{T}, mesh::MeshData{T}, configs::AbstractDict) where {T}
    BC_top::String = configs["BC"]["top"]["type"]
    BC_bottom::String = configs["BC"]["bottom"]["type"]
    pressure_grad::Vector{Float64} = configs["BodyForce"]["PressureGradient"]["values"]
    @assert length(pressure_grad) == 2 "Mean profile: Pressure gradient must be a vector of length 2."
    pressure_grad_mag = sqrt(sum(pressure_grad .^ 2))
    nu_inv::Float64 = configs["Re"]

    if BC_top == "NoSlipWall" && BC_bottom == "NoSlipWall" # both boundaries are walls
        println("Mean profile: BC_top wall && BC_bottom wall")
        z_distance = map(z -> min(z, 1 - z), mesh.zz) .* mesh.Lz
        u_tau = sqrt(pressure_grad_mag * mesh.Lz / 2)
    elseif BC_top == "GradientWall" && BC_bottom == "NoSlipWall" # bottom is a wall, top is stress free
        println("Mean profile: BC_top free slip && BC_bottom wall")
        z_distance = map(z -> z, mesh.zz)
        u_tau = sqrt(pressure_grad_mag * mesh.Lz) .* mesh.Lz
    else
        error("Mean profile: Unsupported boundary condition combination.")
    end
    println("Mean profile: u_tau = ", u_tau)

    # scale the distance from the wall using the viscous units
    z_distance_viscous = z_distance .* (u_tau * nu_inv)

    # compute the mean profile
    u_plus_mean = composite_mean_profile(z_distance_viscous)
    u_plus_x = u_plus_mean .* (pressure_grad[1] / pressure_grad_mag)
    u_plus_y = u_plus_mean .* (pressure_grad[2] / pressure_grad_mag)
    flow.u .+= reshape(u_plus_x .* u_tau, 1, 1, size(flow.u, 3))
    flow.v .+= reshape(u_plus_y .* u_tau, 1, 1, size(flow.u, 3))

    return flow
end

function smooth_z!(a::AbstractArray{T,3}) where {T}
    tmpp = similar(a)
    for i in 1:3
        @inbounds tmpp[:, :, 2:end-1] =
            a[:, :, 2:end-1] .* 0.5 +
            a[:, :, 1:end-2] * 0.25 +
            a[:, :, 3:end] * 0.25
        @inbounds a[:, :, 2:end-1] .= tmpp[:, :, 2:end-1]
    end
    return a
end

function add_disturbance_velocity!(flow::FlowField{T}, mesh::MeshData{T}, configs::AbstractDict,
    u_mag::Real, v_mag::Real, w_mag::Real) where {T}

    seed = 91384773
    Random.seed!(seed)

    kx_cutoff = size(flow.u, 1) ÷ 3 ÷ 2
    ky_cutoff = size(flow.u, 2) ÷ 3 ÷ 2
    tmp = similar(flow.u)
    uk = zeros(Complex{T}, size(tmp, 1) ÷ 2 + 1, size(tmp, 2), size(tmp, 3))
    plan = plan_rfft(tmp, (1, 2))
    i_plan = plan_irfft(uk, size(tmp, 1), (1, 2))
    for (var_name, u, mag) in (("u", flow.u, u_mag), ("v", flow.v, v_mag), ("w", flow.w, w_mag))
        println("Fluctuations: Generating random fluctuations and smoothing $var_name")
        Random.randn!(tmp)
        tmp .*= mag

        # spectral low-pass filter
        mul!(uk, plan, tmp)
        uk[kx_cutoff:end, :, :] .= 0
        uk[:, ky_cutoff:end+2-ky_cutoff, :] .= 0
        mul!(tmp, i_plan, uk)

        # set near-boundary values to zero
        tmp[:, :, 1:8] .= 0
        tmp[:, :, end-7:end] .= 0

        smooth_z!(tmp)

        u .+= tmp
    end
end

function write_restart(filename::String, mesh::MeshData{T}, flow::FlowField{T}, time::Float64=0.0) where {T<:Real}
    h5open(filename, "w") do h5f
        h5f["time"] = time
        h5f["u"] = flow.u
        h5f["v"] = flow.v
        h5f["w"] = flow.w
        h5f["pp"] = zeros(T, size(flow.u)...)

        h5f["zz"] = mesh.zz
        h5f["zw"] = mesh.zw
        h5f["dz"] = mesh.dz
        h5f["dzw"] = mesh.dzw
        h5f["pex"] = mesh.kx0
        h5f["pey"] = mesh.ky0
        h5f["hbar"] = mesh.Lz
    end
    nothing
end

function main()
    cli = set_commandline_opts()
    cmd_options = parse_args(cli)
    file_options = TOML.parsefile(pop!(cmd_options, "config"))
    option_section = pop!(cmd_options, "section")

    # extract the section where the parameters are stored
    param_section = file_options[option_section]

    output_file = cmd_options["output"]
    if (ispath(output_file))
        filepath = islink(output_file) ? "$output_file ($(realpath(output_file)))" : output_file
        if (!cmd_options["overwrite"])
            @error "Output file $filepath already exists. Use --overwrite to overwrite."
            return
        end
        @warn "Overwriting $filepath"
    end

    # read the mesh
    filename = cmd_options["grid"]
    mesh = h5open(filename, "r") do grid_file
        zz = grid_file["zz"][:]
        zw = grid_file["zw"][:]
        dz = grid_file["dz"][:]
        dzw = grid_file["dzw"][:]
        kx0 = grid_file["pex"][]
        ky0 = grid_file["pey"][]
        Lz = grid_file["hbar"][]

        return MeshData{eltype(zz)}(zz, zw, dz, dzw, kx0, ky0, Lz)
    end

    # read the parameters
    nx::Integer = param_section["grid_size"][1]
    ny::Integer = param_section["grid_size"][2]
    nz::Integer = param_section["grid_size"][3]
    if (nz != size(mesh.zz, 1))
        error("The number of points in the grid file does not match the number of points in the configuration file.")
    end

    RealT = eltype(mesh.zz)
    kx0 = RealT(param_section["domain_size"][1])
    ky0 = RealT(param_section["domain_size"][2])
    Lz = RealT(param_section["domain_size"][3])
    if !isapprox(kx0, mesh.kx0) || !isapprox(ky0, mesh.ky0) || !isapprox(Lz, mesh.Lz)
        @warn "The mesh sizes from configuration file do not match the grid file."
        @warn "The mesh sizes from configuration file are used: kx0 = $kx0, ky0 = $ky0, Lz = $Lz"
        mesh = MeshData{RealT}(mesh.zz, mesh.zw, mesh.dz, mesh.dzw, kx0, ky0, Lz)
    end

    flow = FlowField{RealT}((nx, ny, nz))

    add_disturbance_velocity!(flow, mesh, param_section, 10.0, 20.0, 20.0)
    add_parabolic_profile!(flow, mesh)

    write_restart(output_file, mesh, flow, 0.0)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
