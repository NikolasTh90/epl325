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
#include <omp.h>
#include <string.h>
// Chnage this value to reflect your ID number
#define ID 1030496
typedef struct
{
	double x, y, z;
} vector;

int bodies, timeSteps;
int SimulationTime = 0;
int threads = 1;
double *masses, GravConstant;
vector *positions, *velocities, *accelerations;
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
double getRND()
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




void resolveCollisions(char* exec_type)
{
	int i, j;
	double dx, dy, dz, md;
	// int velocity_swaps[bodies-1][bodies];

	// if (strcmp(exec_type, "serial") == 0){
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
			return;
	// }




	// # pragma omp parallel private(i,j,dx,dy,dz,md) shared(bodies,masses,positions,velocity_swaps,velocities,threads, exec_type) default(none)
	// {
	// if (strcmp(exec_type, "static") == 0){
	// 	# pragma omp for schedule(static, bodies/threads)
	// 	for(i=0;i<bodies-1;i++){
	// 		int step = 0;
	// 		for(j=i+1;j<bodies;j++){
	// 			md = masses[i]+masses[j];
	// 			dx = fabs(positions[i].x-positions[j].x);
	// 			dy = fabs(positions[i].y-positions[j].y);
	// 			dz = fabs(positions[i].z-positions[j].z);
	// 			if(dx<md && dy<md && dz<md){
	// 				//Swap Velocities
					
	// 				// Store the swap
	// 				velocity_swaps[i][step++] = j;
	// 			}
	// 		}
	// 		velocity_swaps[i][step] = -1;
	// 	}
	// }
	// if(strcmp(exec_type, "dynamic")==0){
	// 	# pragma omp for schedule(dynamic)
	// 	for(i=0;i<bodies-1;i++){
	// 		int step = 0;
	// 		for(j=i+1;j<bodies;j++){
	// 			md = masses[i]+masses[j];
	// 			dx = fabs(positions[i].x-positions[j].x);
	// 			dy = fabs(positions[i].y-positions[j].y);
	// 			dz = fabs(positions[i].z-positions[j].z);
	// 			if(dx<md && dy<md && dz<md){
	// 				//Swap Velocities
					
	// 				// Store the swap
	// 				velocity_swaps[i][step++] = j;
	// 			}
	// 		}
	// 		velocity_swaps[i][step] = -1;
	// 	}
	// }
	// if(strcmp(exec_type, "guided")==0){
	// 	# pragma omp for schedule(guided)
	// 	for(i=0;i<bodies-1;i++){
	// 		int step = 0;
	// 		for(j=i+1;j<bodies;j++){
	// 			md = masses[i]+masses[j];
	// 			dx = fabs(positions[i].x-positions[j].x);
	// 			dy = fabs(positions[i].y-positions[j].y);
	// 			dz = fabs(positions[i].z-positions[j].z);
	// 			if(dx<md && dy<md && dz<md){
	// 				//Swap Velocities
					
	// 				// Store the swap
	// 				velocity_swaps[i][step++] = j;
	// 			}
	// 		}
	// 		velocity_swaps[i][step] = -1;
	// 	}
	// }

	// # pragma omp barrier
	// 	# pragma omp single
	// 	{	
	// 		for (i=0; i<bodies-1; i++){
	// 			j = 0;
	// 			while(velocity_swaps[i][j] != -1){
	// 				vector temp = velocities[i];
	// 				velocities[i] = velocities[velocity_swaps[i][j]];
	// 				velocities[velocity_swaps[i][j]] = temp;
	// 				j++;
	// 			}
	// 		}
	// 	}


		


	
	// }
}

void computeAccelerations(char* exec_type)
{
	int i,j;
	if (strcmp(exec_type, "static") == 0){
		#pragma omp parallel private(i, j) default(none) shared(accelerations, positions, masses, GravConstant, threads, bodies)
		{
		#pragma omp for schedule(static,bodies/threads) 
		for (i = 0; i < bodies; i++)
		{
			accelerations[i].x = 0;
			accelerations[i].y = 0;
			accelerations[i].z = 0;
		// #pragma omp parallel for schedule(static)
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
	}
	else if (strcmp(exec_type, "dynamic") == 0) {
		#pragma omp parallel private(i, j) default(none) shared(accelerations, positions, masses, GravConstant, bodies, threads)
		{
			int chunksize = (threads > 1) ? (int)(bodies / threads * log(threads)) : bodies;
		#pragma omp for schedule(dynamic, chunksize)
		for (i = 0; i < bodies; i++)
		{
			accelerations[i].x = 0;
			accelerations[i].y = 0;
			accelerations[i].z = 0;
		// #pragma omp parallel for schedule(static)
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
	
	}

	else if (strcmp(exec_type, "guided") == 0){
		#pragma omp parallel private(i, j) default(none) shared(accelerations, positions, masses, GravConstant, bodies, threads)
		{
		#pragma omp for schedule(guided, bodies/threads) 
		for (i = 0; i < bodies; i++)
		{
			accelerations[i].x = 0;
			accelerations[i].y = 0;
			accelerations[i].z = 0;
		// #pragma omp parallel for schedule(static)
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

	}
	else {

		for (i = 0; i < bodies; i++)
		{
			accelerations[i].x = 0;
			accelerations[i].y = 0;
			accelerations[i].z = 0;
		// #pragma omp parallel for schedule(static)
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
}

void computeVelocities(char* exec_type)
{
	int i;
	// if(strcmp(exec_type, "static") == 0)
	// {
	// 	#pragma omp parallel private(i) default(none) shared(accelerations, velocities, threads, bodies)
	// 	{
	// 		#pragma omp for schedule(static, bodies/threads)
	// 			for (i = 0; i < bodies; i++)
	// 			{
	// 				// velocities[i] = addVectors(velocities[i],accelerations[i]);
	// 				vector ac = {velocities[i].x + accelerations[i].x, velocities[i].y + accelerations[i].y, velocities[i].z + accelerations[i].z};
	// 				velocities[i] = ac;
	// 			}
	// 	}
	// }
	// else if (strcmp(exec_type,"dynamic") == 0)
	// {
	// 	#pragma omp parallel private(i) default(none) shared(accelerations, velocities, threads, bodies)
	// 	{
	// 		#pragma omp for schedule(dynamic)
	// 			for (i = 0; i < bodies; i++)
	// 			{
	// 				// velocities[i] = addVectors(velocities[i],accelerations[i]);
	// 				vector ac = {velocities[i].x + accelerations[i].x, velocities[i].y + accelerations[i].y, velocities[i].z + accelerations[i].z};
	// 				velocities[i] = ac;
	// 			}
	// 	}
	// }
	// else if (strcmp(exec_type,"guided") == 0)
	// {
	// 	#pragma omp parallel private(i) default(none) shared(accelerations, velocities, threads, bodies)
	// 	{
	// 		#pragma omp for schedule(guided)
			 
	// 			for (i = 0; i < bodies; i++)
	// 			{
	// 				// velocities[i] = addVectors(velocities[i],accelerations[i]);
	// 				vector ac = {velocities[i].x + accelerations[i].x, velocities[i].y + accelerations[i].y, velocities[i].z + accelerations[i].z};
	// 				velocities[i] = ac;
	// 			}
			
	// 	}
	// }
	// else
	// {
		for (i = 0; i < bodies; i++)
		{
			// velocities[i] = addVectors(velocities[i],accelerations[i]);
			vector ac = {velocities[i].x + accelerations[i].x, velocities[i].y + accelerations[i].y, velocities[i].z + accelerations[i].z};
			velocities[i] = ac;
		}
			
	// }


}

void computePositions(char* exec_type)
{
	int i;
	// if (strcmp(exec_type, "static") == 0)
	// {
	// 	#pragma omp parallel private(i) default(none) shared(accelerations, velocities, positions, threads, bodies)
	// 	{
	// 		#pragma omp for schedule(static, bodies/threads)
			
	// 			for (i = 0; i < bodies; i++)
	// 			{
	// 				// positions[i] = addVectors(positions[i],addVectors(velocities[i],scaleVector(0.5,accelerations[i])));
	// 				vector sc = {0.5 * accelerations[i].x, 0.5 * accelerations[i].y, 0.5 * accelerations[i].z};
	// 				vector ac = {velocities[i].x + sc.x, velocities[i].y + sc.y, velocities[i].z + sc.z};
	// 				vector bc = {positions[i].x + ac.x, positions[i].y + ac.y, positions[i].z + ac.z};
	// 				positions[i] = bc;
	// 			}
				
	// 	}	
	// }
	// else if (strcmp(exec_type, "dynamic") == 0)
	// {
	// 	#pragma omp parallel private(i) default(none) shared(accelerations, velocities, positions, threads, bodies)
	// 	{
	// 		#pragma omp for schedule(dynamic)
			
	// 			for (i = 0; i < bodies; i++)
	// 			{
	// 				// positions[i] = addVectors(positions[i],addVectors(velocities[i],scaleVector(0.5,accelerations[i])));
	// 				vector sc = {0.5 * accelerations[i].x, 0.5 * accelerations[i].y, 0.5 * accelerations[i].z};
	// 				vector ac = {velocities[i].x + sc.x, velocities[i].y + sc.y, velocities[i].z + sc.z};
	// 				vector bc = {positions[i].x + ac.x, positions[i].y + ac.y, positions[i].z + ac.z};
	// 				positions[i] = bc;
	// 			}
				
	// 	}	
	// }
	// else if (strcmp(exec_type, "guided") == 0)
	// {
	// 	#pragma omp parallel private(i) default(none) shared(accelerations, velocities, positions, threads, bodies)
	// 	{
	// 		#pragma omp for schedule(guided)
			
	// 			for (i = 0; i < bodies; i++)
	// 			{
	// 				// positions[i] = addVectors(positions[i],addVectors(velocities[i],scaleVector(0.5,accelerations[i])));
	// 				vector sc = {0.5 * accelerations[i].x, 0.5 * accelerations[i].y, 0.5 * accelerations[i].z};
	// 				vector ac = {velocities[i].x + sc.x, velocities[i].y + sc.y, velocities[i].z + sc.z};
	// 				vector bc = {positions[i].x + ac.x, positions[i].y + ac.y, positions[i].z + ac.z};
	// 				positions[i] = bc;
	// 			}
				
	// 	}	
	// }
	// else
	for (i = 0; i < bodies; i++)
	{
		// positions[i] = addVectors(positions[i],addVectors(velocities[i],scaleVector(0.5,accelerations[i])));
		vector sc = {0.5 * accelerations[i].x, 0.5 * accelerations[i].y, 0.5 * accelerations[i].z};
		vector ac = {velocities[i].x + sc.x, velocities[i].y + sc.y, velocities[i].z + sc.z};
		vector bc = {positions[i].x + ac.x, positions[i].y + ac.y, positions[i].z + ac.z};
		positions[i] = bc;
	}
}

void simulate()
{
	SimulationTime++;
	char* execution_type = "serial";
	computeAccelerations(execution_type);
	computePositions(execution_type);
	computeVelocities(execution_type);
	resolveCollisions(execution_type);
}

void n_body_omp_static(int threads)
{
	SimulationTime++;
	omp_set_num_threads(threads);
	char* execution_type = "static";

	computeAccelerations(execution_type);
	computePositions(execution_type);
	computeVelocities(execution_type);
	resolveCollisions(execution_type);
}

void n_body_omp_dynamic(int threads)
{
	SimulationTime++;
	omp_set_num_threads(threads);
	char* execution_type = "dynamic";

	computeAccelerations(execution_type);
	computePositions(execution_type);
	computeVelocities(execution_type);
	resolveCollisions(execution_type);
}

void n_body_omp_guided(int threads)
{
	SimulationTime++;
	omp_set_num_threads(threads);
	char* execution_type = "guided";

	computeAccelerations(execution_type);
	computePositions(execution_type);
	computeVelocities(execution_type);
	resolveCollisions(execution_type);
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



int myMain(int argc, char *argv[], char* exec_type)
{
	int i;
	char* log_file = calloc(1, strlen("logfile_") + strlen(exec_type) + 1);
	strcpy(log_file, "logfile_");
	strcat(log_file, exec_type);
	char* data_file = calloc(1, strlen("data_") + strlen(exec_type) + 1);
	strcpy(data_file, "data_");
	strcat(data_file, exec_type);		

	FILE *lfp = fopen(log_file, "w");
	FILE *dfp = fopen(data_file, "w");
	if (lfp == NULL || dfp == NULL)
	{
		printf("Please create the ./outputs directory\n");
		return -1;
	}

	if(argc == 3){
		timeSteps = atoi(argv[1]);
		bodies = atoi(argv[2]);
		threads = 1;
		printf("%%*** RUNNING WITH VALUES %d timesteps %d bodies and %d thread(s) ***\n",timeSteps, bodies, threads);

	}
	else
	if (argc == 4)
	{
		timeSteps = atoi(argv[1]);
		bodies = atoi(argv[2]);
		threads = atoi(argv[3]);
		printf("%%*** RUNNING WITH VALUES %d timesteps %d bodies and %d thread(s) ***\n",timeSteps, bodies, threads);

	}
	else
	if (argc == 2){
		timeSteps = 10000;
		bodies = 200;
		threads = atoi(argv[1]);
		printf("%%*** RUNNING WITH DEFAULT VALUES %d timesteps %d bodies and %d thread(s) ***\n",timeSteps, bodies, threads);

		
	}

	else
	{
		timeSteps = 10000;
		bodies = 200;
		threads = 1;
		printf("%%*** RUNNING WITH DEFAULT VALUES %d timesteps %d bodies and %d thread(s) ***\n",timeSteps, bodies, threads);

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
			if(strcmp(exec_type, "static") == 0)
				n_body_omp_static(threads);
			else	
				if(strcmp(exec_type, "dynamic") == 0)
					n_body_omp_dynamic(threads);
			else		
				if(strcmp(exec_type, "guided") == 0)
					n_body_omp_guided(threads);		
			else
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
	fclose(lfp);
	fclose(dfp);
	return 0;
}

int main(int argc, char *argv[])
{
	printf("Static\n");	
	fflush(stdout);
	myMain(argc, argv, "static");

	printf("Dynamic\n");
	fflush(stdout);
	myMain(argc, argv, "dynamic");

	
	printf("Guided\n");
	fflush(stdout);
	myMain(argc, argv, "guided");

	
	// printf("Serial\n");
	// fflush(stdout);
	// myMain(argc, argv, "serial");

	return 0;
}