# Features

## Boundary conditions

The code implements a set of boundary conditions for the top and bottom boundaries. The BCs can be specified in the configuration file at keys `BC.top` and `BC.bottom`.

The following BCs are available:

- __No-slip wall__: Uniform velocities at the boundaries. The parameters for the BC are included the {class}`alps::solver::NoSlipWall` class.
```{eval-rst}
.. doxygenstruct:: alps::solver::NoSlipWall
```

- __Fixed gradient wall__: Uniform gradients $\partial u/\partial z$ and $\partial v/\partial z$ at the boundaries. The component $w$ is 0. The parameters for the BC are included the {class}`alps::solver::GradientWall` class.
```{eval-rst}
.. doxygenstruct:: alps::solver::GradientWall
```

- __Fixed tangential stress__: Specified tangential stresses at the boundaries. The velocity component $w$ is 0. THe parameters for the BC are included the {class}`alps::solver::TangentialStressWall` class.
```{eval-rst}
.. doxygenstruct:: alps::solver::TangentialStressWall
```

- __Monochromatic wave__: The boundary velocity and $\eta_t$ are specified as a monochromatic wave. The BC is implemented by the {class}`alps::solver::MonochromaticWaveWall` class.
```{eval-rst}
.. doxygenstruct:: alps::solver::MonochromaticWaveWall
```

## SGS models

The code implements three subgrid-scale (SGS) models: the Smagorinsky model, the Germano-Lilly dynamic Smagorinsky model, and the anisotropic minimum dissipation model. The SGS model is specified in the configuration file under a table named `LES.SGS`, e.g.:
```toml
["parent table".LES.SGS]
model = "cs"
```
The `LES.SGS` table is usually under another table, e.g. `default_solver` (the full name will be `[default_solver.LES.SGS]`). Options for the model can be listed in the same table.

In the code, the options for these SGS models are stored in {class}`alps::solver::ConstantSmagorinskyOptions`, {class}`alps::solver::DynamicSmagorinskyOptions`, and {class}`alps::solver::AnisotropicMinimumDissipationOptions` classes, respectively. The available options for each model are detailed in the documentation of the classes:

```{eval-rst}
.. doxygenstruct:: alps::solver::ConstantSmagorinskyOptions

.. doxygenstruct:: alps::solver::DynamicSmagorinskyOptions

.. doxygenstruct:: alps::solver::AnisotropicMinimumDissipationOptions
```