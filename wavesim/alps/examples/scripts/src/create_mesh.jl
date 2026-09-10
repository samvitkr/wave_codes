using ArgParse
using TOML
using HDF5
using AlpsTools

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
        "--output", "-o"
        help = "Specify the output file"
        arg_type = String
        default = "grid.h5"
        "--f32"
        help = "Create in Float32 precision"
        action = :store_true
        "--overwrite"
        help = "Overwrite an existing file"
        action = :store_true
    end

    return s
end

function main()
    cli = set_commandline_opts()
    cmd_options = parse_args(cli)
    file_options = TOML.parsefile(pop!(cmd_options, "config"))
    option_section = pop!(cmd_options, "section")

    # extract the section where the parameters are stored
    param_section = file_options[option_section]

    filename = cmd_options["output"]
    if (ispath(filename))
        filepath = islink(filename) ? "$filename ($(realpath(filename)))" : filename
        if (!cmd_options["overwrite"])
            @error "Output file $filepath already exists. Use --overwrite to overwrite."
            return
        end
        @warn "Overwriting existing file: $filepath"
    end

    # create the mesh
    RealT = cmd_options["f32"] ? Float32 : Float64
    mesh = TwoSideClusteredMesh{RealT}(
        param_section["grid_size"][3] - 2,
        param_section["grid_beta"],
        param_section["grid_gamma"],
    )
    kx0::RealT = param_section["domain_size"][1]
    ky0::RealT = param_section["domain_size"][2]
    Lz::RealT = param_section["domain_size"][3]
    mesh_data = MeshData(mesh, kx0, ky0, Lz)

    # write out the mesh
    filename = cmd_options["output"]
    h5open(filename, "w") do h5f
        h5f["zz"] = mesh_data.zz
        h5f["zw"] = mesh_data.zw
        h5f["dz"] = mesh_data.dz
        h5f["dzw"] = mesh_data.dzw
        h5f["pex"] = mesh_data.kx0
        h5f["pey"] = mesh_data.ky0
        h5f["hbar"] = mesh_data.Lz
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end