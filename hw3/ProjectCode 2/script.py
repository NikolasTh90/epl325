def create_sbatch(o_flag=0, number_of_threads=40):
    with open('sbatch', 'w') as sbatch:
        sbatch.write('#!/bin/bash\n')
        sbatch.write('#SBATCH --partition=CSUG\n')
        sbatch.write('### Partition\n')
        sbatch.write('#SBATCH --nodes=1\n')
        sbatch.write('### Number of Nodes\n')
        sbatch.write('#SBATCH --time=0:01:00\n')
        sbatch.write('### WallTime\n')
        sbatch.write('#SBATCH --ntasks=1\n')
        sbatch.write('### Number of tasks\n')
        sbatch.write('#SBATCH --ntasks-per-node=1 ### Number of tasks (MPI processes)\n')
        sbatch.write('#SBATCH --cpus-per-task=' + str(number_of_threads) + ' ### Number of threads per task (OMP threads)\n')
        sbatch.write('#SBATCH --job-name=coun3\n')
        sbatch.write('### Job name\n')
        sbatch.write('#SBATCH --output=output.txt ### Output file\n')
        sbatch.write('lscpu > lscpu_output.txt\n')
        sbatch.write('gcc -lpthread -fopenmp -Werror -Wall -lm -O' + str(o_flag) + ' n-body_std.c\n')
        sbatch.write('./a.out ' + str(number_of_threads) + '\n')
    
def finished():
    counter=0
    with open('output.txt', 'r') as output:
        for line in output.readlines():
            if line.__contains__('Simulation Time'):
                counter+=1
    return True if counter==3 else False

def submit_sbatch():
    import subprocess
    output = subprocess.check_output(['sbatch', 'sbatch'])

def get_simulation_times():
    times = list()
    with open('output.txt', 'r') as output:
        for line in output.readlines():
            if line.__contains__('Simulation Time'):
                times.append(line.split(' ')[-1][:-1])
    return times

from time import *
times = list()
for o_flag in (0,3):
    for threads in range(1,41):
        f = open('output.txt', 'w+')
        f.close()
        print('Creating O' + str(o_flag) + ' with ' + str(threads) + ' threads')
        create_sbatch(o_flag=o_flag, number_of_threads=threads)
        submit_sbatch()


        while not finished():
            pass   
        times.append([str(o_flag) + ',' + str(threads)] + get_simulation_times())

# create_sbatch()
# submit_sbatch()
# sleep(60)
# times = list()
# hello = 123
# times.append(['O3,' + str(hello)] + get_simulation_times())
# times.append(['O0,2'] + get_simulation_times())
# print(times)

with open('times.csv', 'w') as timesfile:
    timesfile.write('\"O flag\", threads, static, dynamic\n')
    for sbatch_try in times:
        for record in sbatch_try:
            timesfile.write(record + ',')
        timesfile.write('\n')