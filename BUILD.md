# Build instructions for the libary
You will need CMake and another build system like Make or Ninja.

### Normal build
If you just want to build everything, just do the following
```
mkdir -p build
cd build/
cmake ..
make
```

You will now find two executables ```mincostmatching_test``` (for running the tests) and ```mincostmatching_example``` (for running the example).
Apart from that, the actual library will be built as well and is called ```libmincostmatching.a```.


### Use as a library in CMake
Just clone the repository (possibly as a submodule), add the project as a subdirectory and 
include the src directory. Additionally, link your project against ```mincostmatching```.

If you don't want to build tests and the example, add the following lines before the add_subdirectory
```
set(MCPM_BUILD_TESTS CACHE BOOL OFF)
set(MCPM_BUILD_EXAMPLE CACHE BOOL OFF)
```
