#!/bin/sh
find src test benchmark examples  -name "CMakeLists.txt" -print -exec cmake-format -i {} \;
cmake-format -i cmake/NVCCWrapperUtils.cmake
cmake-format -i cmake/StaticAnalyzers.cmake
cmake-format -i cmake/Formatter.cmake
cmake-format -i cmake/Dependencies.cmake
cmake-format -i cmake/CompilerSelection.cmake
cmake-format -i cmake/CompileOptions.cmake
cmake-format -i cmake/CompilerWarnings.cmake
find cmake/deps -name "*.cmake" -print -exec cmake-format -i {} \;
