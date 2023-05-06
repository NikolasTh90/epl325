/* How to Run
Compile Using:
  gcc -Werror -Wall -O3 -lm n-body_std.c
Run Using: ./a.out [NumberOfInterations NumberOfBodies]
	./a.out 10000 200
	gnuplot plot3D.sh
For gprof:
	gcc -Werror -Wall -lm -pg n-body_std.c
	./a.out 10000 200
	gprof ./a.out > analysis.txt
	gprof ./a.out | ./gprof2dot.py | dot -Tpng -o gprof_output.png
For perf:
	 perf record -g -- ./a.out
	 perf script | c++filt | ./gprof2dot.py -f perf | dot -Tpng -o perf_output.png

Code Ref:https://rosettacode.org/wiki/N-body_problem#C
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "support.h"
// Chnage this value to reflect your ID number
#define ID 1030496
typedef struct
{
	double x, y, z;
} vector;

int bodies, timeSteps;
int SimulationTime = 0;
double *masses, GravConstant;
vector *positions, *velocities, *accelerations;
vector *cuda_positions, *cuda_velocities, *cuda_accelerations;
double *cuda_masses, *cuda_GravConstant;
vector *cuda_local_accelerations;



	/*
	vector addVectors(vector a,vector b){
		vector c = {a.x+b.x,a.y+b.y,a.z+b.z};

		return c;
	}

	vector scaleVector(double b,vector a){
		vector c = {b*a.x,b*a.y,b*a.z};

		return c;
	}

	vector subtractVectors(vector a,vector b){
		vector c = {a.x-b.x,a.y-b.y,a.z-b.z};

		return c;
	}

	double mod(vector a){
		return sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
	}
	*/
	double
	getRND()
{
	return (100.0 - rand() % 200) / 50.0;
}
void initiateSystemRND(int bodies)
{
	int i;
	srand(ID);
	masses = (double *)malloc(bodies * sizeof(double));
	positions = (vector *)malloc(bodies * sizeof(vector));
	velocities = (vector *)malloc(bodies * sizeof(vector));
	accelerations = (vector *)malloc(bodies * sizeof(vector));
	GravConstant = 0.01;
	for (i = 0; i < bodies; i++)
	{
		masses[i] = 0.4; //(rand()%100)/100.0;//0.4 Good Example
		positions[i].x = getRND();
		positions[i].y = getRND();
		positions[i].z = getRND();
		velocities[i].x = getRND() / 5.0; // 0.0;
		velocities[i].y = getRND() / 5.0; // 0.0;
		velocities[i].z = getRND() / 5.0; // 0.0;
	}

	
	// Allocate device memory
	cudaMalloc((void **)&cuda_masses, bodies * sizeof(double));
	cudaMalloc((void **)&cuda_positions, bodies * sizeof(vector));
	cudaMalloc((void **)&cuda_accelerations, bodies * sizeof(vector));
	cudaMalloc((void **)&cuda_local_accelerations, bodies*bodies * sizeof(vector));

	// Transfer data from host to device
	cudaMemcpyAsync(cuda_masses, masses, bodies * sizeof(double), cudaMemcpyHostToDevice);
}
/*
void initiateSystem(char* fileName){
	int i;
	FILE* fp = fopen(fileName,"r");
	fscanf(fp,"%lf%d",&GravConstant,&bodies);

	masses = (double*)malloc(bodies*sizeof(double));
	positions = (vector*)malloc(bodies*sizeof(vector));
	velocities = (vector*)malloc(bodies*sizeof(vector));
	accelerations = (vector*)malloc(bodies*sizeof(vector));

	for(i=0;i<bodies;i++){
		fscanf(fp,"%lf",&masses[i]);
		fscanf(fp,"%lf%lf%lf",&positions[i].x,&positions[i].y,&positions[i].z);
		fscanf(fp,"%lf%lf%lf",&velocities[i].x,&velocities[i].y,&velocities[i].z);
	}
	fclose(fp);
}
*/
void resolveCollisions()
{
	int i, j;
	double dx, dy, dz, md;

	for (i = 0; i < bodies - 1; i++)
		for (j = i + 1; j < bodies; j++)
		{
			md = masses[i] + masses[j];
			dx = fabs(positions[i].x - positions[j].x);
			dy = fabs(positions[i].y - positions[j].y);
			dz = fabs(positions[i].z - positions[j].z);
			// if(positions[i].x==positions[j].x && positions[i].y==positions[j].y && positions[i].z==positions[j].z){
			if (dx < md && dy < md && dz < md)
			{
// Swap Velocities
#ifdef DEBUG
				fprintf(stderr, "T=%d;%lf:%lf:%lf<->%lf:%lf:%lf", SimulationTime, positions[i].x, positions[i].y, positions[i].z, positions[j].x, positions[j].y, positions[j].z);
				fprintf(stderr, "[md:%lf::%lf:%lf:%lf]", md, dx, dy, dz);
				fprintf(stderr, "\tCollision(%d):%d<->%d\n", SimulationTime, i, j);
#endif
				vector temp = velocities[i];
				velocities[i] = velocities[j];
				velocities[j] = temp;
			}
		}
}

#include <cuda.h>
#include <cuda_runtime.h>

// CUDA kernel to compute accelerations for a subset of bodies
__global__ void computeAccelerationsKernel(int bodies, double GravConstant, vector *positions, double *masses, vector *cuda_local_accelerations)
{
	int i = blockIdx.x;
	// int threadID = threadIdx.x;
	//  Calculate the Unique Thread ID
	int threadUniqueID = threadIdx.x + blockIdx.x * blockDim.x;
	// int size =  dInfo[0];
	int j = threadIdx.x;
	
	if(i < bodies) {
			if (i != j)
			{
				vector sij = {positions[i].x - positions[j].x, positions[i].y - positions[j].y, positions[i].z - positions[j].z};
				double mod = sqrt(sij.x * sij.x + sij.y * sij.y + sij.z * sij.z);
				double mod3 = mod * mod * mod;
				double s = GravConstant * masses[j] / mod3;
				vector S = {s * sij.x, s * sij.y, s * sij.z};
				cuda_local_accelerations[threadUniqueID].x, = S.x;
				cuda_local_accelerations[threadUniqueID].y, = S.y;
				cuda_local_accelerations[threadUniqueID].z, = S.z;
			}

			
	}
		
}

__global__ void reduction_accelerations(int bodies, vector* cuda_local_accelerations, vector* cuda_accelerations){
	int threadUniqueID = threadIdx.x + blockIdx.x * blockDim.x;
	if (threadUniqueID < bodies)
    {
		cuda_accelerations[threadUniqueID].x = 0;
		cuda_accelerations[threadUniqueID].y = 0;
		cuda_accelerations[threadUniqueID].z = 0;
		for(int i = 0; i < bodies; i++){
			cuda_accelerations[threadUniqueID].x += cuda_local_accelerations[i+(threadUniqueID*bodies)].x;
			cuda_accelerations[threadUniqueID].y += cuda_local_accelerations[i+(threadUniqueID*bodies)].y;
			cuda_accelerations[threadUniqueID].z += cuda_local_accelerations[i+(threadUniqueID*bodies)].z;
		}
	}

}

void computeAccelerations_cuda()
{
	
	cudaMemcpy(cuda_positions, positions, bodies * sizeof(vector), cudaMemcpyHostToDevice);

	
	// Launch kernels on multiple streams
	int blockSize = 224;
	int gridSize = (bodies + blockSize - 1) / blockSize;
	
	computeAccelerationsKernel<<<gridSize, blockSize>>>(bodies, GravConstant, cuda_positions, cuda_masses, cuda_local_accelerations );
	

	// Wait for all streams to finish
	cudaDeviceSynchronize();
	blockSize = bodies/16;
	gridSize = (bodies + blockSize - 1) / blockSize;
	reduction_accelerations<<<gridSize, blockSize>>>(bodies,cuda_local_accelerations, cuda_accelerations);

	// Transfer data back to host
	cudaMemcpy(accelerations, cuda_accelerations, bodies * sizeof(vector), cudaMemcpyDeviceToHost);

	
	
}



void computeAccelerations()
{
	int i, j;
	for (i = 0; i < bodies; i++)
	{
		accelerations[i].x = 0;
		accelerations[i].y = 0;
		accelerations[i].z = 0;
		for (j = 0; j < bodies; j++)
		{
			if (i != j)
			{
				// accelerations[i] = addVectors(accelerations[i],scaleVector(GravConstant*masses[j]/pow(mod(subtractVectors(positions[i],positions[j])),3),subtractVectors(positions[j],positions[i])));
				vector sij = {positions[i].x - positions[j].x, positions[i].y - positions[j].y, positions[i].z - positions[j].z};
				vector sji = {positions[j].x - positions[i].x, positions[j].y - positions[i].y, positions[j].z - positions[i].z};
				double mod = sqrt(sij.x * sij.x + sij.y * sij.y + sij.z * sij.z);
				double mod3 = mod * mod * mod;
				double s = GravConstant * masses[j] / mod3;
				vector S = {s * sji.x, s * sji.y, s * sji.z};
				accelerations[i].x += S.x;
				accelerations[i].y += S.y;
				accelerations[i].z += S.z;
			}
		}
	}
}

void computeVelocities()
{
	int i;
	for (i = 0; i < bodies; i++)
	{
		// velocities[i] = addVectors(velocities[i],accelerations[i]);
		vector ac = {velocities[i].x + accelerations[i].x, velocities[i].y + accelerations[i].y, velocities[i].z + accelerations[i].z};
		velocities[i] = ac;
	}
}

void computePositions()
{
	int i;
	for (i = 0; i < bodies; i++)
	{
		// positions[i] = addVectors(positions[i],addVectors(velocities[i],scaleVector(0.5,accelerations[i])));
		vector sc = {0.5 * accelerations[i].x, 0.5 * accelerations[i].y, 0.5 * accelerations[i].z};
		vector ac = {velocities[i].x + sc.x, velocities[i].y + sc.y, velocities[i].z + sc.z};
		vector bc = {positions[i].x + ac.x, positions[i].y + ac.y, positions[i].z + ac.z};
		positions[i] = bc;
	}
}


void nBody_CUDA()
{
	SimulationTime++;
	computeAccelerations_cuda();
	computePositions();
	computeVelocities();
	resolveCollisions();
}


void simulate()
{
	SimulationTime++;
	computeAccelerations();
	computePositions();
	computeVelocities();
	resolveCollisions();
}

void printBodiesInfo(FILE *lfp, FILE *dfp)
{
	int j;
	for (j = bodies - 10; j < bodies; j++)
		fprintf(lfp, "Body%d %f\t: %lf\t%f\t%lf\t|\t%lf\t%lf\t%lf\n", j + 1, masses[j], positions[j].x, positions[j].y, positions[j].z, velocities[j].x, velocities[j].y, velocities[j].z);
	fprintf(lfp, "-------------------------------------------------------------------------------------------\n");
	for (j = bodies - 10; j < bodies; j++)
		fprintf(stdout, "Body%d %f\t: %lf\t%f\t%lf\t|\t%lf\t%lf\t%lf\n", j + 1, masses[j], positions[j].x, positions[j].y, positions[j].z, velocities[j].x, velocities[j].y, velocities[j].z);
	fprintf(stdout, "-------------------------------------------------------------------------------------------\n");
}

int main(int argc, char *argv[])
{
	int i;
	FILE *lfp = fopen("./outputs/logfile.txt", "w");
	FILE *dfp = fopen("./outputs/data.dat", "w");
	if (lfp == NULL || dfp == NULL)
	{
		printf("Please create the ./outputs directory\n");
		return -1;
	}
	if (argc == 3)
	{
		timeSteps = atoi(argv[1]);
		bodies = atoi(argv[2]);
	}
	else
	{
		printf("%%*** RUNNING WITH DEFAULT VALUES ***\n");
		timeSteps = 10000;
		bodies = 200;
	}
	initiateSystemRND(bodies);
	// initiateSystem("input.txt");
	fprintf(stdout, "Running With %d Bodies for %d timeSteps. Initial state:\n", bodies, timeSteps);
	fprintf(stderr, "Running With %d Bodies for %d timeSteps. Initial state:\n", bodies, timeSteps);
	fprintf(lfp, "Running With %d Bodies for %d timeSteps. Initial state:\n", bodies, timeSteps);
	fprintf(lfp, "Body   \t\t\t:\t\tx\t\ty\t\t\tz\t\t|\t\tvx\t\t\tvy\t\t\tvz\t\t\n");
	printBodiesInfo(lfp, dfp);
	startTime(0);
	for (i = 0; i < timeSteps; i++)
	{
		simulate();
#ifdef DEBUG
		int j;
		// printf("\nCycle %d\n",i+1);
		for (j = 0; j < bodies; j++)
			fprintf(dfp, "%d\t%d\t%lf\t%lf\t%lf\n", i, j, positions[j].x, positions[j].y, positions[j].z);
#endif
	}
	stopTime(0);
	fprintf(lfp, "\nLast Step = %d\n", i);
	printBodiesInfo(lfp, dfp);
	printf("\nSimulation Time:");
	elapsedTime(0);


	initiateSystemRND(bodies);
	// initiateSystem("input.txt");
	fprintf(stdout, "Running With %d Bodies for %d timeSteps. Initial state:\n", bodies, timeSteps);
	fprintf(stderr, "Running With %d Bodies for %d timeSteps. Initial state:\n", bodies, timeSteps);
	fprintf(lfp, "Running With %d Bodies for %d timeSteps. Initial state:\n", bodies, timeSteps);
	fprintf(lfp, "Body   \t\t\t:\t\tx\t\ty\t\t\tz\t\t|\t\tvx\t\t\tvy\t\t\tvz\t\t\n");
	printBodiesInfo(lfp, dfp);
	startTime(0);
	for (i = 0; i < timeSteps; i++)
	{
		nBody_CUDA();
#ifdef DEBUG
		int j;
		// printf("\nCycle %d\n",i+1);
		for (j = 0; j < bodies; j++)
			fprintf(dfp, "%d\t%d\t%lf\t%lf\t%lf\n", i, j, positions[j].x, positions[j].y, positions[j].z);
#endif
	}
	stopTime(0);
	fprintf(lfp, "\nLast Step = %d\n", i);
	printBodiesInfo(lfp, dfp);
	printf("\nSimulation Time:");
	elapsedTime(0);





	fclose(lfp);
	fclose(dfp);


	cudaFree(d_masses);
	cudaFree(d_positions);
	cudaFree(d_accelerations);
	return 0;
}