# HIPC CUDA Lab Example Solutions
Example solutions for the University of York High Performance Computing (HIPC) CUDA course (Unit 7).


## Project Organization
### Hello, world
- `src/hello.cu`: "hello, world" CUDA version
- `src/hello_multi.cu`: parallel version of "hello, world" 

### Vect_add
- `src/vect_add.c`: base version of vectorization addition.
- `src/vect_add_single.cu`: vect_add with single cuda thread (<<<1,1>>>)
- `src/vect_add_multithread.cu`: vect_add with multiple cuda threads
- `src/vect_add_multigrid.cu`: vect_add with multiple grids

### Monte Carlo Simulation
- `src/monte_carlo.c`: C base version of Monte Carlo.
- `src/monte_carlo.cu`: Monte Carlo implemented in CUDA

### Performance
- `src/memory_coalescing.cu`: exploring memory coalescing in CUDA

### Misc
- `src/queryDevice.cu`: query device and environment properties


## Preparation
Make sure you have `nvidia-driver` and `cuda-toolkit` installed.

To test the environment, run:

`$ nvidia-smi`

and 

`$ nvcc --version`

both should return information without an error.


## Usage
In a new terminal, `cd` to the project root folder, and compile:

`$ make`

Then run the relevant program in `bin/`.


## Run Jobs on Viking
An example is provided in `cuda_viking.job`, to use it, login to Viking and in the terminal:

`$ sbatch cuda_viking.job`
