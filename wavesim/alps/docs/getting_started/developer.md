# Developer guide

### Development environment

The development environment provides a set of tools and configurations to help contributors develop the code. Although for running simulations, the libraries (particularly MPI) provided by the machine should generally preferred and therefore [system-specific setup](linux.md) is more relevant. The development environment is designed to be self-contained and to provide a consistent environment across different systems, including local Linux machines. The tools in the environment also work with Visual Studio Code plugins to enable code navigation and completion.

To set up the development environment, follow the steps below:

1. Install Pixi (skip this step if you already have Pixi installed):
    ```bash
    curl -fsSL https://pixi.sh/install.sh | bash
    ```
   or follow the instructions on the [Pixi website](https://pixi.sh).

2. Enter the project source directory, ensure that `pixi.toml` and `pixi.lock` are present, and run:
    ```bash
    pixi install -a --frozen
    ```

3. Activate the development environment:
    ```bash
    pixi shell -e cuda # environment with CUDA support
    ```
    or
    ```bash
    pixi shell # default environment for CPU-only builds
    ```

### Visual Studio Code setup

1. To use Visual Studio Code with the development environment, install the following extensions:
    - [clangd](https://marketplace.visualstudio.com/items?itemName=llvm-vs-code-extensions.vscode-clangd)
    - [CMake Tools](https://marketplace.visualstudio.com/items?itemName=ms-vscode.cmake-tools)
    - [C/C++](https://marketplace.visualstudio.com/items?itemName=ms-vscode.cpptools)
    - [CMake (optional)](https://marketplace.visualstudio.com/items?itemName=twxs.cmake)

2. Set up `CMakeUserPresets.json` with configure presets that use the installed development environment. An example is provided at `scripts/CMakeUserPresets.json`. If you do not have an existing `CMakeUserPresets.json`, you can copy the example to the project root directory.
    ```
    cp scripts/CMakeUserPresets.json .
    ```

3. Open the project directory in Visual Studio Code.

4. Select the CMake pane in the left sidebar. If under `Configure` you see `No Configure Preset Selected`, click on `Select Configure Preset` to the right of it and choose a preset from the dropdown list, e.g. `local-dev-cuda`.

5. Click `Configure` to configure the project.

6. Click `Build` to build the project.