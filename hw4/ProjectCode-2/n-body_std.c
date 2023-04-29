#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "support.h"
#include <immintrin.h>
#include <omp.h>
// Chnage this value to reflect your ID number
#define ID 1038399
typedef struct
{
    double x, y, z;
} vector;

int bodies, timeSteps;
int SimulationTime = 0;
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
    masses = (double *)aligned_alloc(64, bodies * sizeof(double)); // allocate aligned memory
    positions = (vector *)aligned_alloc(64, bodies * sizeof(vector));
    velocities = (vector *)aligned_alloc(64, bodies * sizeof(vector));
    accelerations = (vector *)aligned_alloc(64, bodies * sizeof(vector));
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

// void resolveCollisionsParallel(int threads)
// {
//     int i, j;
//     double dx, dy, dz, md;
//     int arr[bodies-1][bodies]; // Create an array to hold the swaps -> maximum swaps = bodies-1 * bodies
// #pragma omp parallel shared(velocities, positions, arr) private(i, j, md, dx, dy, dz) firstprivate(bodies, masses, threads) default(none)
//     {
// #pragma omp for schedule(static, bodies / threads)
//         for (i = 0; i < bodies - 1; i++)
//         {
//             int swap = 0;
//             for (j = i + 1; j < bodies; j++)
//             {
//                 md = masses[i] + masses[j];
//                 dx = fabs(positions[i].x - positions[j].x);
//                 dy = fabs(positions[i].y - positions[j].y);
//                 dz = fabs(positions[i].z - positions[j].z);
//                 // if(positions[i].x==positions[j].x && positions[i].y==positions[j].y && positions[i].z==positions[j].z){
//                 if (dx < md && dy < md && dz < md)
//                 {
// // Swap Velocities
// #ifdef DEBUG
//                     fprintf(stderr, "T=%d;%lf:%lf:%lf<->%lf:%lf:%lf", SimulationTime, positions[i].x, positions[i].y, positions[i].z, positions[j].x, positions[j].y, positions[j].z);
//                     fprintf(stderr, "[md:%lf::%lf:%lf:%lf]", md, dx, dy, dz);
//                     fprintf(stderr, "\tCollision(%d):%d<->%d\n", SimulationTime, i, j);
// #endif
//                     // Find in parallel all swaps
//                     // this means that all swaps of i can be found in arr[i] row with all j
//                     arr[i][swap] = j; // save the swap i,j in format arr[i][swap] = j
//                     swap++;
//                 }
//             }
//             arr[i][swap] = bodies + 1; // when swaps finish for row i save in next position bodies+1 to indicate that swap finished
//         }
//     }
//     // swaps will executed only from one thread to avoid wrong swaps among the many threads
//     for (i = 0; i < bodies-1; i++) // for all rows
//     {
//         j = 0;
//         while (arr[i][j] != bodies + 1) // Swap i,j until bodies + 1 number found in array
//         {
//             vector temp = velocities[i];
//             velocities[i] = velocities[arr[i][j]];
//             velocities[arr[i][j]] = temp;
//             j++;
//         }
//     }
// }

// void resolveCollisionsParallel2(int threads)
// {
//     int i, j;
//     double dx, dy, dz, md;
//     int arr[bodies-1][bodies];
// #pragma omp parallel shared(velocities, positions, arr) private(i, j, md, dx, dy, dz) firstprivate(bodies, masses, threads) default(none)
//     {
// #pragma omp for schedule(dynamic)
//         for (i = 0; i < bodies - 1; i++)
//         {
//             int swap = 0;
//             for (j = i + 1; j < bodies; j++)
//             {
//                 md = masses[i] + masses[j];
//                 dx = fabs(positions[i].x - positions[j].x);
//                 dy = fabs(positions[i].y - positions[j].y);
//                 dz = fabs(positions[i].z - positions[j].z);
//                 // if(positions[i].x==positions[j].x && positions[i].y==positions[j].y && positions[i].z==positions[j].z){
//                 if (dx < md && dy < md && dz < md)
//                 {
// // Swap Velocities
// #ifdef DEBUG
//                     fprintf(stderr, "T=%d;%lf:%lf:%lf<->%lf:%lf:%lf", SimulationTime, positions[i].x, positions[i].y, positions[i].z, positions[j].x, positions[j].y, positions[j].z);
//                     fprintf(stderr, "[md:%lf::%lf:%lf:%lf]", md, dx, dy, dz);
//                     fprintf(stderr, "\tCollision(%d):%d<->%d\n", SimulationTime, i, j);
// #endif
//                     arr[i][swap] = j;
//                     swap++;
//                 }
//             }
//             arr[i][swap] = bodies + 1;
//         }
//     }
//     for (i = 0; i < bodies-1; i++)
//     {
//         j = 0;
//         while (arr[i][j] != bodies + 1)
//         {
//             vector temp = velocities[i];
//             velocities[i] = velocities[arr[i][j]];
//             velocities[arr[i][j]] = temp;
//             j++;
//         }
//     }
// }

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

void computeAccelerations_parallel(int threads)
{
    int i, j;
    __m512d gc = _mm512_set1_pd(GravConstant); // Move GravConstant to gc
    __m512d s_ij = _mm512_setzero_pd();        //  Case i==j

#pragma omp parallel shared(accelerations, positions) private(i, j) firstprivate(bodies, gc, masses, s_ij, threads) default(none)
    {

#pragma omp for schedule(static, bodies / threads)
        for (i = 0; i < bodies; i++)
        {
            __m512d i_if = _mm512_set1_pd(i); // Save i indices for if checking in mask
            // Init acceleration to zero
            __m512d acc_x = _mm512_setzero_pd();
            __m512d acc_y = _mm512_setzero_pd();
            __m512d acc_z = _mm512_setzero_pd();
            // Get position x y z for body i. Save the same pos 8 times
            __m512d pos_x = _mm512_set1_pd(positions[i].x);
            __m512d pos_y = _mm512_set1_pd(positions[i].y);
            __m512d pos_z = _mm512_set1_pd(positions[i].z);

            for (j = 0; j < bodies; j += 8)
            {
                __m512d pos_xj = _mm512_setr_pd(positions[j].x, positions[j + 1].x, positions[j + 2].x, positions[j + 3].x, positions[j + 4].x, positions[j + 5].x, positions[j + 6].x, positions[j + 7].x);
                __m512d pos_yj = _mm512_setr_pd(positions[j].y, positions[j + 1].y, positions[j + 2].y, positions[j + 3].y, positions[j + 4].y, positions[j + 5].y, positions[j + 6].y, positions[j + 7].y);
                __m512d pos_zj = _mm512_setr_pd(positions[j].z, positions[j + 1].z, positions[j + 2].z, positions[j + 3].z, positions[j + 4].z, positions[j + 5].z, positions[j + 6].z, positions[j + 7].z);

                __m512d sij_x = _mm512_sub_pd(pos_x, pos_xj);
                __m512d sij_y = _mm512_sub_pd(pos_y, pos_yj);
                __m512d sij_z = _mm512_sub_pd(pos_z, pos_zj);

                __m512d sji_x = _mm512_sub_pd(pos_xj, pos_x);
                __m512d sji_y = _mm512_sub_pd(pos_yj, pos_y);
                __m512d sji_z = _mm512_sub_pd(pos_zj, pos_z);

                __m512d mod = _mm512_sqrt_pd(_mm512_fmadd_pd(sij_x, sij_x, _mm512_fmadd_pd(sij_y, sij_y, _mm512_mul_pd(sij_z, sij_z))));
                __m512d mod3 = _mm512_mul_pd(mod, _mm512_mul_pd(mod, mod));

                __m512d j_if = _mm512_setr_pd(j, j + 1, j + 2, j + 3, j + 4, j + 5, j + 6, j + 7);
                __mmask8 mask = _mm512_cmp_pd_mask(i_if, j_if, _CMP_EQ_OQ);
                __m512d s_inotj = _mm512_div_pd(_mm512_mul_pd(gc, _mm512_load_pd(&masses[j])), mod3);
                __m512d s = _mm512_mask_blend_pd(mask, s_inotj, s_ij);

                __m512d Sx = _mm512_mul_pd(s, sji_x);
                __m512d Sy = _mm512_mul_pd(s, sji_y);
                __m512d Sz = _mm512_mul_pd(s, sji_z);

                acc_x = _mm512_add_pd(acc_x, Sx);
                acc_y = _mm512_add_pd(acc_y, Sy);
                acc_z = _mm512_add_pd(acc_z, Sz);
            }
            accelerations[i].x = _mm512_reduce_add_pd(acc_x);
            accelerations[i].y = _mm512_reduce_add_pd(acc_y);
            accelerations[i].z = _mm512_reduce_add_pd(acc_z);
        }
    }
}

void computeAccelerations_parallel2(int threads)
{
    int i, j;
    __m512d gc = _mm512_set1_pd(GravConstant); // Move GravConstant to gc
    __m512d s_ij = _mm512_setzero_pd();        //  Case i==j

#pragma omp parallel shared(accelerations, positions) private(i, j) firstprivate(bodies, gc, masses, s_ij, threads) default(none)
    {

#pragma omp for schedule(dynamic)
        for (i = 0; i < bodies; i++)
        {
            __m512d i_if = _mm512_set1_pd(i); // Save i indices for if checking in mask
            // Init acceleration to zero
            __m512d acc_x = _mm512_setzero_pd();
            __m512d acc_y = _mm512_setzero_pd();
            __m512d acc_z = _mm512_setzero_pd();
            // Get position x y z for body i. Save the same pos 8 times
            __m512d pos_x = _mm512_set1_pd(positions[i].x);
            __m512d pos_y = _mm512_set1_pd(positions[i].y);
            __m512d pos_z = _mm512_set1_pd(positions[i].z);

            for (j = 0; j < bodies; j += 8)
            {
                __m512d pos_xj = _mm512_setr_pd(positions[j].x, positions[j + 1].x, positions[j + 2].x, positions[j + 3].x, positions[j + 4].x, positions[j + 5].x, positions[j + 6].x, positions[j + 7].x);
                __m512d pos_yj = _mm512_setr_pd(positions[j].y, positions[j + 1].y, positions[j + 2].y, positions[j + 3].y, positions[j + 4].y, positions[j + 5].y, positions[j + 6].y, positions[j + 7].y);
                __m512d pos_zj = _mm512_setr_pd(positions[j].z, positions[j + 1].z, positions[j + 2].z, positions[j + 3].z, positions[j + 4].z, positions[j + 5].z, positions[j + 6].z, positions[j + 7].z);

                __m512d sij_x = _mm512_sub_pd(pos_x, pos_xj);
                __m512d sij_y = _mm512_sub_pd(pos_y, pos_yj);
                __m512d sij_z = _mm512_sub_pd(pos_z, pos_zj);

                __m512d sji_x = _mm512_sub_pd(pos_xj, pos_x);
                __m512d sji_y = _mm512_sub_pd(pos_yj, pos_y);
                __m512d sji_z = _mm512_sub_pd(pos_zj, pos_z);

                __m512d mod = _mm512_sqrt_pd(_mm512_fmadd_pd(sij_x, sij_x, _mm512_fmadd_pd(sij_y, sij_y, _mm512_mul_pd(sij_z, sij_z))));
                __m512d mod3 = _mm512_mul_pd(mod, _mm512_mul_pd(mod, mod));

                __m512d j_if = _mm512_setr_pd(j, j + 1, j + 2, j + 3, j + 4, j + 5, j + 6, j + 7);
                __mmask8 mask = _mm512_cmp_pd_mask(i_if, j_if, _CMP_EQ_OQ);
                __m512d s_inotj = _mm512_div_pd(_mm512_mul_pd(gc, _mm512_load_pd(&masses[j])), mod3);
                __m512d s = _mm512_mask_blend_pd(mask, s_inotj, s_ij);

                __m512d Sx = _mm512_mul_pd(s, sji_x);
                __m512d Sy = _mm512_mul_pd(s, sji_y);
                __m512d Sz = _mm512_mul_pd(s, sji_z);

                acc_x = _mm512_add_pd(acc_x, Sx);
                acc_y = _mm512_add_pd(acc_y, Sy);
                acc_z = _mm512_add_pd(acc_z, Sz);
            }
            accelerations[i].x = _mm512_reduce_add_pd(acc_x);
            accelerations[i].y = _mm512_reduce_add_pd(acc_y);
            accelerations[i].z = _mm512_reduce_add_pd(acc_z);
        }
    }
}

void computeAccelerations_parallel3(int threads)
{
    int i, j;
    __m512d gc = _mm512_set1_pd(GravConstant); // Move GravConstant to gc
    __m512d s_ij = _mm512_setzero_pd();        //  Case i==j

#pragma omp parallel shared(accelerations, positions) private(i, j) firstprivate(bodies, gc, masses, s_ij, threads) default(none)
    {

#pragma omp for schedule(guided)
        for (i = 0; i < bodies; i++)
        {
            __m512d i_if = _mm512_set1_pd(i); // Save i indices for if checking in mask
            // Init acceleration to zero
            __m512d acc_x = _mm512_setzero_pd();
            __m512d acc_y = _mm512_setzero_pd();
            __m512d acc_z = _mm512_setzero_pd();
            // Get position x y z for body i. Save the same pos 8 times
            __m512d pos_x = _mm512_set1_pd(positions[i].x);
            __m512d pos_y = _mm512_set1_pd(positions[i].y);
            __m512d pos_z = _mm512_set1_pd(positions[i].z);

            for (j = 0; j < bodies; j += 8)
            {
                __m512d pos_xj = _mm512_setr_pd(positions[j].x, positions[j + 1].x, positions[j + 2].x, positions[j + 3].x, positions[j + 4].x, positions[j + 5].x, positions[j + 6].x, positions[j + 7].x);
                __m512d pos_yj = _mm512_setr_pd(positions[j].y, positions[j + 1].y, positions[j + 2].y, positions[j + 3].y, positions[j + 4].y, positions[j + 5].y, positions[j + 6].y, positions[j + 7].y);
                __m512d pos_zj = _mm512_setr_pd(positions[j].z, positions[j + 1].z, positions[j + 2].z, positions[j + 3].z, positions[j + 4].z, positions[j + 5].z, positions[j + 6].z, positions[j + 7].z);

                __m512d sij_x = _mm512_sub_pd(pos_x, pos_xj);
                __m512d sij_y = _mm512_sub_pd(pos_y, pos_yj);
                __m512d sij_z = _mm512_sub_pd(pos_z, pos_zj);

                __m512d sji_x = _mm512_sub_pd(pos_xj, pos_x);
                __m512d sji_y = _mm512_sub_pd(pos_yj, pos_y);
                __m512d sji_z = _mm512_sub_pd(pos_zj, pos_z);

                __m512d mod = _mm512_sqrt_pd(_mm512_fmadd_pd(sij_x, sij_x, _mm512_fmadd_pd(sij_y, sij_y, _mm512_mul_pd(sij_z, sij_z))));
                __m512d mod3 = _mm512_mul_pd(mod, _mm512_mul_pd(mod, mod));

                __m512d j_if = _mm512_setr_pd(j, j + 1, j + 2, j + 3, j + 4, j + 5, j + 6, j + 7);
                __mmask8 mask = _mm512_cmp_pd_mask(i_if, j_if, _CMP_EQ_OQ);
                __m512d s_inotj = _mm512_div_pd(_mm512_mul_pd(gc, _mm512_load_pd(&masses[j])), mod3);
                __m512d s = _mm512_mask_blend_pd(mask, s_inotj, s_ij);

                __m512d Sx = _mm512_mul_pd(s, sji_x);
                __m512d Sy = _mm512_mul_pd(s, sji_y);
                __m512d Sz = _mm512_mul_pd(s, sji_z);

                acc_x = _mm512_add_pd(acc_x, Sx);
                acc_y = _mm512_add_pd(acc_y, Sy);
                acc_z = _mm512_add_pd(acc_z, Sz);
            }
            accelerations[i].x = _mm512_reduce_add_pd(acc_x);
            accelerations[i].y = _mm512_reduce_add_pd(acc_y);
            accelerations[i].z = _mm512_reduce_add_pd(acc_z);
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

void n_body_omp_static(int threads)
{
    SimulationTime++;
    omp_set_num_threads(threads);
    computeAccelerations_parallel(threads); // omp compute accelleration with static schedule
    computePositions();
    computeVelocities();
    resolveCollisions(); // omp resolve colision with static schedule
}

void n_body_omp_dynamic(int threads)
{
    SimulationTime++;
    omp_set_num_threads(threads);
    computeAccelerations_parallel2(threads); // compute acceleration with dynamic scheduling
    computePositions();
    computeVelocities();
    resolveCollisions(); // omp resolve colision with dynamic schedule
}

void n_body_omp_guided(int threads)
{
    SimulationTime++;
    omp_set_num_threads(threads);
    computeAccelerations_parallel3(threads); // compute acceleration with dynamic scheduling
    computePositions();
    computeVelocities();
    resolveCollisions(); // omp resolve colision with dynamic schedule
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
    fclose(lfp);
    fclose(dfp);

    lfp = fopen("./outputs/logfile.txt", "w");
    dfp = fopen("./outputs/data.dat", "w");
    if (lfp == NULL || dfp == NULL)
    {
        printf("Please create the ./outputs directory\n");
        return -1;
    }
    return 0;
}
