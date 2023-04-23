/* How to Run
Compile Using:
  gcc -fopenmp -Werror -Wall -o3 n-body_std.c -lm
Run Using: ./a.out [NumberOfInterations NumberOfBodies NoOfThreads static/dynamic/guided]
	./a.out 10000 200 10 static
	gnuplot plot3D.sh
For gprof:
	gcc -Werror -Wall -lm -pg n-body_std.c
	./a.out 10000 200
	gprof ./a.out > analysis.txt
	gprof ./a.out | ./gprof2dot.py | dot -Tpng -o gprof_output.png
For perf:
	 perf record -g -- ./a.out 10000 200 4
	 perf script | c++filt | ./gprof2dot.py -f perf | dot -Tpng -o perf_output.png

Code Ref:https://rosettacode.org/wiki/N-body_problem#C
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "support.h"
#include <omp.h>
#include <string.h>
#include <immintrin.h>

#define VECTOR_SIZE 8 // Number of elements per SIMD vector

// Chnage this value to reflect your ID number
#define ID 1028696
#define END -1
typedef struct
{
	double x, y, z;
} vector;

int bodies, timeSteps;
int SimulationTime = 0;
double *masses, GravConstant;
double *positions_x, *positions_y, *positions_z;
double *velocities_x, *velocities_y, *velocities_z;
vector *accelerations;
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

	// Align memory to a 64-byte boundary
    masses = (double *)_mm_malloc(bodies * sizeof(double), 64);
    positions_x = (double *)_mm_malloc(bodies * sizeof(double), 64);
    positions_y = (double *)_mm_malloc(bodies * sizeof(double), 64);
    positions_z = (double *)_mm_malloc(bodies * sizeof(double), 64);
    velocities_x = (double *)_mm_malloc(bodies * sizeof(double), 64);
    velocities_y = (double *)_mm_malloc(bodies * sizeof(double), 64);
    velocities_z = (double *)_mm_malloc(bodies * sizeof(double), 64);
	accelerations = (vector *)_mm_malloc(bodies * sizeof(vector), 64);

	GravConstant = 0.01;
	for (i = 0; i < bodies; i++)
	{
		masses[i] = 0.4; //(rand()%100)/100.0;//0.4 Good Example
		positions_x[i] = getRND();
        positions_y[i] = getRND();
        positions_z[i] = getRND();
        velocities_x[i] = getRND() / 5.0; // 0.0;
        velocities_y[i] = getRND() / 5.0; // 0.0;
        velocities_z[i] = getRND() / 5.0; // 0.0;
	}
}

void resolveCollisionsDynamic(int noOfThreads)
{
	int i, j, changes;
	double dx, dy, dz, md;

	// Matrix to hold the changes to be made in the form of indexes - velocityPositions[bodies-1][bodies]
	int **velocityPositions = (int **)malloc((bodies - 1) * sizeof(int *));
	for (i = 0; i < (bodies - 1); i++)
	{
		velocityPositions[i] = (int *)malloc(bodies * sizeof(int));
	}
// Set the private variables to be used, so that threads don't mix up. Also seth number of threads to be used.
#pragma omp parallel private(i, j, dx, dy, dz, md) num_threads(noOfThreads) shared(bodies, masses, noOfThreads, changes, positions_x,positions_y,positions_z)
	{
// Set the dynamic parallelization of the for loop. Add reduction clause to changes so that the
// changes made to that variable is visible to each thread
#pragma omp for schedule(dynamic) reduction(+ \
											: changes)
		for (i = 0; i < bodies - 1; i++)
		{
			changes = 0;
			for (j = i + 1; j < bodies; j++)
			{
				md = masses[i] + masses[j];
				dx = fabs(positions_x[i] - positions_x[j]);
				dy = fabs(positions_y[i] - positions_y[j]);
				dz = fabs(positions_z[i] - positions_z[j]);
				// if(positions_x[i]==positions_x[j] && positions_y[i]==positions_y[j] && positions_z[i]==positions_z[j]){
				if (dx < md && dy < md && dz < md)
				{
// Swap Velocities
#ifdef DEBUG
					fprintf(stderr, "T=%d;%lf:%lf:%lf<->%lf:%lf:%lf", SimulationTime, positions_x[i], positions_y[i], positions_z[i], positions_x[j], positions_y[j], positions_z[j]);
					fprintf(stderr, "[md:%lf::%lf:%lf:%lf]", md, dx, dy, dz);
					fprintf(stderr, "\tCollision(%d):%d<->%d\n", SimulationTime, i, j);
#endif
					velocityPositions[i][changes] = j; // Save index for later swapping
					changes++;						   // Increment number of changes for this body
				}
			}
			velocityPositions[i][changes] = END; // Set the end position of changes for this body
		}
	}
#pragma omp barrier // Wait until all threads are done
#pragma omp single	// Get a thread to do the swaps
	{
		for (i = 0; i < bodies - 1; i++)
		{
			for (j = 0; velocityPositions[i][j] != -1; j++)
			{
                double vel_x = velocities_x[i]; double vel_y = velocities_y[i]; double vel_z = velocities_z[i];
                velocities_x[i] = velocities_x[velocityPositions[i][j]];
                velocities_y[i] = velocities_y[velocityPositions[i][j]];
                velocities_z[i] = velocities_z[velocityPositions[i][j]];
                velocities_x[velocityPositions[i][j]] = vel_x;
                velocities_y[velocityPositions[i][j]] = vel_y;
                velocities_z[velocityPositions[i][j]] = vel_z;
			}
		}
	}
	for (i = 0; i < (bodies - 1); i++)
	{
		free(velocityPositions[i]);
	}
	free(velocityPositions);
}

void computeAccelerationsDynamic(int noOfThreads)
{
	int i, j;
// Share the workload dynamically between the threads. Privatize the i and j iterators.
#pragma omp parallel for schedule(dynamic) num_threads(noOfThreads) shared(accelerations, positions_x,positions_y,positions_z, masses) private(i, j)
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
				vector sij = {positions_x[i] - positions_x[j], positions_y[i] - positions_y[j], positions_z[i] - positions_z[j]};
				vector sji = {positions_x[j] - positions_x[i], positions_y[j] - positions_y[i], positions_z[j] - positions_z[i]};
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

void resolveCollisionsStatic(int noOfThreads){
	int i,j,changes;
	double dx,dy,dz,md;

    //Matrix to hold the changes to be made in the form of indexes - velocityPositions[bodies-1][bodies]
    int** velocityPositions = (int**) malloc((bodies-1) * sizeof(int*));
    for (i = 0; i < (bodies-1); i++) {
        velocityPositions[i] = (int*) malloc(bodies * sizeof(int));
    }
    //Set the private variables to be used, so that threads don't mix up. Also seth number of threads to be used.
	#pragma omp parallel private(i,j,dx,dy,dz,md) num_threads(noOfThreads) shared(bodies,masses,noOfThreads,changes,positions_x,positions_y,positions_z)
	{
        //Set the static parallelization of the for loop. Add reduction clause to changes so that the
        //changes made to that variable is visible to each thread
		#pragma omp for schedule(static, VECTOR_SIZE) reduction(+: changes)
		for(i=0;i<bodies-1;i++){
			changes = 0;
			for(j=i+1;j<bodies;j++){
				md = masses[i]+masses[j];
				dx = fabs(positions_x[i]-positions_x[j]);
				dy = fabs(positions_y[i]-positions_y[j]);
				dz = fabs(positions_z[i]-positions_z[j]);
				//if(positions_x[i]==positions_x[j] && positions_y[i]==positions_y[j] && positions_z[i]==positions_z[j]){
				if(dx<md && dy<md && dz<md){
					//Swap Velocities
					#ifdef DEBUG
						fprintf(stderr,"T=%d;%lf:%lf:%lf<->%lf:%lf:%lf",SimulationTime,positions_x[i],positions_y[i],positions_z[i],positions_x[j],positions_y[j],positions_z[j]);
						fprintf(stderr,"[md:%lf::%lf:%lf:%lf]",md,dx,dy,dz);
						fprintf(stderr,"\tCollision(%d):%d<->%d\n",SimulationTime,i,j);
					#endif
					velocityPositions[i][changes] = j;  //Save index for later swapping
					changes++;                          //Increment number of changes for this body
				}
			}
			velocityPositions[i][changes] = END;        //Set the end position of changes for this body
		}
	}
	#pragma omp barrier //Wait until all threads are done
	#pragma omp single //Get a thread to do the swaps
	{
		for(i=0;i<bodies-1;i++){
			for(j=0;velocityPositions[i][j]!=-1;j++){
				double vel_x = velocities_x[i]; double vel_y = velocities_y[i]; double vel_z = velocities_z[i];
                velocities_x[i] = velocities_x[velocityPositions[i][j]];
                velocities_y[i] = velocities_y[velocityPositions[i][j]];
                velocities_z[i] = velocities_z[velocityPositions[i][j]];
                velocities_x[velocityPositions[i][j]] = vel_x;
                velocities_y[velocityPositions[i][j]] = vel_y;
                velocities_z[velocityPositions[i][j]] = vel_z;
			}
		}
	}
    for (i = 0; i < (bodies-1); i++) {
        free(velocityPositions[i]);
    }
    free(velocityPositions);
}

void computeAccelerationsStatic(int noOfThreads)
{
    #ifdef DEBUG
    printf("ComputeAccelerationsStatic\n");
    fflush(stdout);
    #endif
	int i, j, k;
    #pragma omp parallel for schedule(static, bodies/noOfThreads) num_threads(noOfThreads) shared(accelerations, positions_x,positions_y,positions_z, masses) private(i, j, k)
	for (i = 0; i < bodies; i ++)  
	{
        #ifdef DEBUG
        printf("i=%d, Thread %d: Initializing acceleration vectors\n",i,omp_get_thread_num());
        fflush(stdout);
        #endif
		// Initialize acceleration vectors for the current batch of bodies
		__m512d ax = _mm512_setzero_pd();
        __m512d ay = _mm512_setzero_pd();
		__m512d az = _mm512_setzero_pd();

        #ifdef DEBUG
        printf("i=%d, Thread %d: Loading positions & masses for i\n",i,omp_get_thread_num());
        fflush(stdout);
        #endif
        //Load position vector for 8 doubles starting from i
        __m512d pix = _mm512_set1_pd(positions_x[i]);
        __m512d piy = _mm512_set1_pd(positions_y[i]);
        __m512d piz = _mm512_set1_pd(positions_z[i]);

        #ifdef DEBUG
        printf("i=%d, Thread %d: ComputeAccelerations for each body\n",i,omp_get_thread_num());
        fflush(stdout);
        #endif
                
		// Compute acceleration for each body in the batch
		for (j = 0; j < bodies; j+=VECTOR_SIZE)
		{
            #ifdef DEBUG
            printf("i=%d, Thread %d: Loading positions & masses for j\n",i,omp_get_thread_num());
            fflush(stdout);
            #endif
			// Load positions and masses for the current batch of bodies
			__m512d pjx = _mm512_load_pd(&positions_x[j]); //load 8 consecutive doubles starting from the address of positions_x[j]
            #ifdef DEBUG
            printf("i=%d, Thread %d: Loaded position x for j=%d\n",i,omp_get_thread_num(),j);
            fflush(stdout);
            #endif
			__m512d pjy = _mm512_load_pd(&positions_y[j]); //load 8 consecutive doubles starting from the address of positions_y[j]
            #ifdef DEBUG
            printf("i=%d, Thread %d: Loaded position y for j=%d\n",i,omp_get_thread_num(),j);
            fflush(stdout);
            #endif
			__m512d pjz = _mm512_load_pd(&positions_z[j]); //load 8 consecutive doubles starting from the address of positions_z[j]
            #ifdef DEBUG
            printf("i=%d, Thread %d: Loaded position z for j=%d\n",i,omp_get_thread_num(),j);
            fflush(stdout);
            #endif
			__m512d m = _mm512_load_pd(&masses[j]); //load 8 consecutive doubles starting from the address of masses[j]
            #ifdef DEBUG
            printf("i=%d, Thread %d: Loaded mass for j=%d\n",i,omp_get_thread_num(),j);
            fflush(stdout);
            #endif

            // Create a mask for the current batch of bodies, where the elements are 0 if j equals i and 1 otherwise (for if(i!=j))
            __mmask8 mask = (1 << VECTOR_SIZE) - 1;
            mask &= ~(1 << i);


            #ifdef DEBUG
            printf("i=%d, Thread %d: Compute distances and inverse distances\n",i,omp_get_thread_num());
            fflush(stdout);
            #endif
			// Compute distances for sij (Subtract 8 position doubles at the same time, for x,y and z)
			__m512d sij_x = _mm512_sub_pd(pix, pjx);
			__m512d sij_y = _mm512_sub_pd(piy, pjy);
			__m512d sij_z = _mm512_sub_pd(piz, pjz);

            // Compute distances for sji (Subtract 8 position doubles at the same time, for x,y and z)
			__m512d sji_x = _mm512_sub_pd(pjx, pix);
			__m512d sji_y = _mm512_sub_pd(pjy, piy);
			__m512d sji_z = _mm512_sub_pd(pjz, piz);

            //Get mod3
			__m512d inside_sqrt = _mm512_add_pd(_mm512_add_pd(_mm512_mul_pd(sij_x, sij_x), _mm512_mul_pd(sij_y, sij_y)), _mm512_mul_pd(sij_z, sij_z));
            __m512d mod = _mm512_sqrt_pd(inside_sqrt);
			__m512d mod3 = _mm512_mul_pd(mod, _mm512_mul_pd(mod, mod));

            #ifdef DEBUG
            printf("i=%d, Thread %d: Compute accelerations for the current batch of bodies\n",i,omp_get_thread_num());
            fflush(stdout);
            #endif
			// Compute accelerations for the current batch of bodies
			__m512d s = _mm512_mul_pd(_mm512_set1_pd(GravConstant), _mm512_div_pd(m, mod3));
            
            //Calculate S
            __m512d S_x = _mm512_mul_pd(s,sji_x);
            __m512d S_y = _mm512_mul_pd(s,sji_y);
            __m512d S_z = _mm512_mul_pd(s,sji_z);

            // Reduce S using the mask to only store the result for the j's that are not equal to i
            S_x = _mm512_maskz_mov_pd(mask, S_x);
            S_y = _mm512_maskz_mov_pd(mask, S_y);
            S_z = _mm512_maskz_mov_pd(mask, S_z);

			ax = _mm512_add_pd(ax, S_x);
			ay = _mm512_add_pd(ay, S_y);
			az = _mm512_add_pd(az, S_z);
		}

        #ifdef DEBUG
        printf("i=%d, Thread %d: Store accelerations for the current batch of bodies\n",i,omp_get_thread_num());
        fflush(stdout);
        #endif
		// Store accelerations for the current batch of bodies
		for (k = 0; k < VECTOR_SIZE; k++)
		{
			int index = i + k;
			if (index < bodies)
			{
				accelerations[index].x += ax[k];
				accelerations[index].y += ay[k];
				accelerations[index].z += az[k];
			}
		}
	}
}

void resolveCollisions()
{
	int i, j;
	double dx, dy, dz, md;

	for (i = 0; i < bodies - 1; i++)
		for (j = i + 1; j < bodies; j++)
		{
			md = masses[i] + masses[j];
			dx = fabs(positions_x[i] - positions_x[j]);
			dy = fabs(positions_y[i] - positions_y[j]);
			dz = fabs(positions_z[i] - positions_z[j]);
			// if(positions_x[i]==positions_x[j] && positions_y[i]==positions_y[j] && positions_z[i]==positions_z[j]){
			if (dx < md && dy < md && dz < md)
			{
				// Swap Velocities
				/*#ifdef DEBUG
								fprintf(stderr, "T=%d;%lf:%lf:%lf<->%lf:%lf:%lf", SimulationTime, positions_x[i], positions_y[i], positions_z[i], positions_x[j], positions_y[j], positions_z[j]);
								fprintf(stderr, "[md:%lf::%lf:%lf:%lf]", md, dx, dy, dz);
								fprintf(stderr, "\tCollision(%d):%d<->%d\n", SimulationTime, i, j) GravConstant;
				#endif*/
                double vel_x = velocities_x[i]; double vel_y = velocities_y[i]; double vel_z = velocities_z[i];
                velocities_x[i] = velocities_x[j];
                velocities_y[i] = velocities_y[j];
                velocities_z[i] = velocities_z[j];
                velocities_x[j] = vel_x;
                velocities_y[j] = vel_y;
                velocities_z[j] = vel_z;
			}
		}
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
				vector sij = {positions_x[i] - positions_x[j], positions_y[i] - positions_y[j], positions_z[i] - positions_z[j]};
				vector sji = {positions_x[j] - positions_x[i], positions_y[j] - positions_y[i], positions_z[j] - positions_z[i]};
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
		vector ac = {velocities_x[i] + accelerations[i].x, velocities_y[i] + accelerations[i].y, velocities_z[i] + accelerations[i].z};
		velocities_x[i] = ac.x; velocities_y[i] = ac.y; velocities_z[i] = ac.z;
	}
}

void computePositions()
{
	int i;
	for (i = 0; i < bodies; i++)
	{
		// positions[i] = addVectors(positions[i],addVectors(velocities[i],scaleVector(0.5,accelerations[i])));
		vector sc = {0.5 * accelerations[i].x, 0.5 * accelerations[i].y, 0.5 * accelerations[i].z};
		vector ac = {velocities_x[i] + sc.x, velocities_y[i] + sc.y, velocities_z[i] + sc.z};
		vector bc = {positions_x[i] + ac.x, positions_y[i] + ac.y, positions_z[i] + ac.z};
        positions_x[i] = bc.x; positions_y[i] = bc.y; positions_z[i] = bc.z;
	}
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
		fprintf(lfp, "Body%d %f\t: %lf\t%f\t%lf\t|\t%lf\t%lf\t%lf\n", j + 1, masses[j], positions_x[j], positions_y[j], positions_z[j], velocities_x[j], velocities_y[j], velocities_z[j]);
	fprintf(lfp, "-------------------------------------------------------------------------------------------\n");
	for (j = bodies - 10; j < bodies; j++)
		fprintf(stdout, "Body%d %f\t: %lf\t%f\t%lf\t|\t%lf\t%lf\t%lf\n", j + 1, masses[j], positions_x[j], positions_y[j], positions_z[j], velocities_x[j], velocities_y[j], velocities_z[j]);
	fprintf(stdout, "-------------------------------------------------------------------------------------------\n");
}

void n_body_omp_static(int noOfThreads)
{
	SimulationTime++;
	// Functions for static parallelization
	computeAccelerationsStatic(noOfThreads);
	computePositions();
	computeVelocities();
	resolveCollisionsStatic(noOfThreads);
}

void n_body_omp_dynamic(int noOfThreads)
{
	SimulationTime++;
	// Functions for dynamic parallelization
	computeAccelerationsDynamic(noOfThreads);
	computePositions();
	computeVelocities();
	resolveCollisionsDynamic(noOfThreads);
}

int main(int argc, char *argv[])
{
	int i;
	int noOfThreads; // Default thread value

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
		if (timeSteps <= 0 || bodies <= 0)
		{
			printf("Number of particles and timesteps must be positive integers\n");
			return 1;
		}
	}
	else if (argc == 5)
	{
		timeSteps = atoi(argv[1]);
		bodies = atoi(argv[2]);
		if (timeSteps <= 0 || bodies <= 0)
		{
			printf("Number of particles and timesteps must be positive integers\n");
			return 1;
		}
		int num = atoi(argv[3]);
		// Make sure the 3rd argument is a number, for the threads
		if (num == 0 && argv[3][0] != '0')
		{
			printf("Error %s is not a number\n", argv[3]);
			return 1;
		}
		else
		{
			noOfThreads = atoi(argv[3]);
		}
		// #ifdef DEBUG
		printf("<=======THREAD COUNT = %d========>\n", noOfThreads);
		// #endif
		const char *schedule = argv[4];
		if (strcmp(schedule, "static") != 0 && strcmp(schedule, "dynamic") != 0 && strcmp(schedule, "guided") != 0)
		{
			printf("Schedule must be 'static', 'dynamic', or 'guided'\n");
			return 1;
		}
	}
	else
	{
		printf("%%*** RUNNING WITH DEFAULT VALUES ***\n");
		timeSteps = 10000;
		bodies = 200;
	}

	initiateSystemRND(bodies);
	// initiateSystem("input.txt");
	if (argc == 5)
	{
		fprintf(stdout, "Running With %d Bodies for %d timeSteps using %d threads. Initial state:\n", bodies, timeSteps, noOfThreads);
		fprintf(stderr, "Running With %d Bodies for %d timeSteps using %d threads. Initial state:\n", bodies, timeSteps, noOfThreads);
	}
	else
	{
		fprintf(stdout, "Running With %d Bodies for %d timeSteps. Initial state:\n", bodies, timeSteps);
		fprintf(stderr, "Running With %d Bodies for %d timeSteps. Initial state:\n", bodies, timeSteps);
	}
	fprintf(lfp, "Running With %d Bodies for %d timeSteps. Initial state:\n", bodies, timeSteps);
	fprintf(lfp, "Body   \t\t\t:\t\tx\t\ty\t\t\tz\t\t|\t\tvx\t\t\tvy\t\t\tvz\t\t\n");
	printBodiesInfo(lfp, dfp);
	startTime(0);
	if (argc != 5)
	{
		for (i = 0; i < timeSteps; i++)
		{
			simulate();
#ifdef DEBUG
			int j;
			// printf("\nCycle %d\n",i+1);
			for (j = 0; j < bodies; j++)
				fprintf(dfp, "%d\t%d\t%lf\t%lf\t%lf\n", i, j, positions_x[j], positions_y[j], positions_z[j]);
#endif
		}
	}
	else
	{
		if (strcmp(argv[4], "static") == 0)
		{
			printf("Using static parallelization\n");
			for (i = 0; i < timeSteps; i++)
			{
				n_body_omp_static(noOfThreads);
			}
		}
		else
		{
			if (strcmp(argv[4], "dynamic") == 0)
			{
				printf("Using dynamic parallelization\n");
				for (i = 0; i < timeSteps; i++)
				{
					n_body_omp_dynamic(noOfThreads);
				}
			}
			else
			{
				printf("Invalid arguments in command line");
			}
		}
	}
	// n_body_omp_static(noOfThreads);	//Change to the static var numbers functon
	stopTime(0);
	// fprintf(lfp,"\nLast Step = %d\n",i);
	printBodiesInfo(lfp, dfp);
	printf("\nSimulation Time:");
	elapsedTime(0);
	fclose(lfp);
	fclose(dfp);
	free(positions_x);
    free(positions_y);
    free(positions_z);
	free(velocities_x);
    free(velocities_y);
    free(velocities_z);
	free(accelerations);
	free(masses);
	return 0;
}