rm a.out;
nvcc -gencode arch=compute_50,code=sm_50 -Wno-deprecated-gpu-targets -O3 $1 -o a.out;
./a.out;
