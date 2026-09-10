Solver
=======

A solver, from bottom up, typically consists of the mesh class, flow field class, and solver class. For example, the following diagram shows the relations among {class}`alps::solver::ChannelFlowSolver`, {class}`alps::solver::FlowField`, and {class}`alps::solver::Mesh`:

```{graphviz}

digraph solverRelation {
    node [shape=record];

    ChannelFlowSolver [label="alps::solver::ChannelFlowSolver"];
    FlowField [label="{alps::solver::FlowField|vec_u: Vector3Field\l pp: HaloView\l...}"];
    Mesh [label="{alps::solver::Mesh|grid: SpectralGrid \l zz: HaloView \l zw: HaloView \l ...}"];

    ChannelFlowSolver -> FlowField [label="flow_field", arrowtail=odiamond, dir=back];
    FlowField -> Mesh [label="mesh",  arrowtail=odiamond, dir=back];
}
```

## `Mesh` and `CurvilinearMesh`
The {class}`alps::solver::Mesh` class is responsible for storing the mesh information, as well as methods for loading and storing the mesh from/to files. The mesh information includes the grid coordinates, the grid spacing, and the domain decomposition ({type}`alps::Grid`). {class}`alps::solver::CurvilinearMesh` is derived from {class}`alps::solver::Mesh`. In addition to the data shared with a rectangular mesh, the transformation data, such as the bottom geometry, the Jacobian, and bottom geometry movement, are also stored.

## `FlowField` and `FlowOverWaveField`
The {class}`alps::solver::FlowField` class stores the flow variables, such as the velocity field, pressure field, and boundary conditions. It also contains a reference to a {class}`alps::solver::Mesh` object. The {class}`alps::solver::FlowOverWaveField` is similar to {class}`alps::solver::FlowField`, but contains the reference to a {class}`alps::solver::CurvilinearMesh` object.

The code snippet below shows how to construct the flow field and mesh objects. A more complete code can be found in the directory `examples/channel`.
```cpp
using namespace alps;

Grid const grid(PencilPlan(comm, {nx, ny, nz}), pex, pey);
solver::Mesh mesh(grid, 1.0); // 1.0 is the height of the domain

solver::FlowField flow(mesh);
flow.set_top_bc(std::make_unique<solver::NoSlipWall>());
flow.set_bottom_bc(std::make_unique<solver::NoSlipWall>());

auto solver_options = solver::ChannelFlowSolverOptions{};
solver::ChannelFlowSolver solver(flow, solver_options);
```