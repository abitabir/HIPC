all:
	mkdir -p bin
	nvcc -o bin/hello src/hello.cu 
	nvcc -o bin/hello_multi src/hello_multi.cu
	nvcc -o bin/queryDevice src/queryDevice.cu 
	nvcc -o bin/vect_add_single src/vect_add_single.cu 
	nvcc -o bin/vect_add_multithread src/vect_add_multithread.cu 
	nvcc -o bin/vect_add_multigrid src/vect_add_multigrid.cu 
	nvcc -o bin/memory_coalescing src/memory_coalescing.cu
	nvcc -o bin/monte_carlo src/monte_carlo.cu -lcurand
	nvcc -o bin/shared_memory src/shared_memory.cu 
	gcc -o bin/vect_add_cpu src/vect_add.c
	gcc -o bin/monte_carlo_cpu src/monte_carlo.c -lm

.PHONY: clean
clean:
	rm -rf bin/*