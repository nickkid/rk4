# -*- coding: utf-8 -*-
import numpy as np
import subprocess
import os
import time
from multiprocessing import Pool
'''
0. global constant definitions
'''
rk4_path = './rk4'
experiment_data_path = ['/mnt/d/lab_data_in_paper1.csv',
'/mnt/d/lab_data_in_paper2.csv']
number_of_process = 8
number_of_iteration = 100000
MAX_OF_SAVE = 50
concentration_of_NO = [1.53e-6, 6.6e-6]
concentration_of_sGC = [0.47e-6, 0.47e-6]
concentration_of_6C = 0
concentration_of_5C = 0
length_of_optical_cell = 1
start_time = 0
end_time = [130, 11]
eps = 1e-3
h1 = 1e-6
hmin = 1e-40
'''
1. set lower bound and upper bound for k's
'''
number_of_k = 5
k = np.zeros(number_of_k)
lower_bound_of_k = np.zeros(number_of_k)
upper_bound_of_k = np.zeros(number_of_k)

lower_bound_of_k[0] = 1.55e8
upper_bound_of_k[0] = 1e10

lower_bound_of_k[1] = 1e-6
upper_bound_of_k[1] = 1e2

lower_bound_of_k[2] = 1e3
upper_bound_of_k[2] = 1e6

lower_bound_of_k[3] = 1e-8
upper_bound_of_k[3] = 1e2

lower_bound_of_k[4] = 1e-8
upper_bound_of_k[4] = 1e2

sample_of_k = np.zeros((number_of_iteration, number_of_k))
std_summary = np.zeros((number_of_iteration,9))
log_interval = 10000

'''
input_content: string
input_filename: string
result_filename: string

create input file for rk4 executable and call rk4 with argument input_filename, result_filename.

'''
def run_simulation(input_content, experiment_data_path, index, program_start_time):
    '''
    with open(input_filename, 'w') as f:
        f.write(input_content)  
    '''
    res = subprocess.check_output([rk4_path,
    str(input_content[0]),
    str(input_content[1]),
    str(input_content[2]),
    str(input_content[3]),
    str(input_content[4]),
    str(input_content[5]),
    str(input_content[6]),
    str(input_content[7]),
    str(input_content[8]),
    str(input_content[9]),
    str(input_content[10]),
    str(input_content[11]),
    str(input_content[12]),
    str(input_content[13]),
    str(input_content[14]),
    experiment_data_path])
    simulation_result = np.genfromtxt([res], delimiter=',')
    '''
    result formatted as:
    "%e,%e\n":std, time_of_computation
    '''
    group_order = index + 1
    if (group_order % log_interval == 0):
        current_time = time.time()
        print('After {0:.3} seconds, one case in group {1:} finished.'.format(current_time - program_start_time ,group_order))
    '''
    if os.path.exists(input_filename):
        os.remove(input_filename)
    if os.path.exists(result_filename):
        os.remove(result_filename)
    '''
    return [index, simulation_result]

def save_result(result_filename, std_summary):
    np.savetxt(result_filename, std_summary, delimiter=',',
    header='k1,k2,k3,k4,k5,std_1,std_2,time_of_computation_1,time_of_computation_2',
    comments='')
'''
2. read in experiment data
'''
#experiment_data = np.genfromtxt(experiment_data_path, skip_header=1, delimiter=',')

'''
3.1 random sampling based on uniform distribution
'''
#with Pool(processes=number_of_process) as pool:
'''
The 'with' statement would not work. Why???
'''

for save_order in range(MAX_OF_SAVE):
    pool = Pool(number_of_process)
    results = [[],[]]
    program_start_time = time.time()
    for i in range(number_of_iteration):
        sample = np.random.uniform(lower_bound_of_k, upper_bound_of_k, size=5)
        sample_of_k [i,:] = sample
        for j in range(2):
            '''
            input_content = '{:.2}\n{:.2}\n{:.2}\n{:.2}\n{:.2}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n'.format(
                                                                                            sample[0],
                                                                                            sample[1],
                                                                                            sample[2],
                                                                                            sample[3],
                                                                                            sample[4], 
                                                                                            length_of_optical_cell,
                                                                                            concentration_of_NO[j],
                                                                                            concentration_of_sGC[j],
                                                                                            concentration_of_6C,
                                                                                            concentration_of_5C,
                                                                                            start_time, end_time,
                                                                                            eps,
                                                                                            h1,
                                                                                            hmin)   
            '''
            input_content = (sample[0],
                            sample[1],
                            sample[2],
                            sample[3],
                            sample[4], 
                            length_of_optical_cell,
                            concentration_of_NO[j],
                            concentration_of_sGC[j],
                            concentration_of_6C,
                            concentration_of_5C,
                            start_time, end_time[j],
                            eps,
                            h1,
                            hmin)                                                                           
            #input_filename = 'input{0:}.{1:}.txt'.format(i, j)
            #result_filename = 'result{0:}.{1:}.csv'.format(i, j)
            res = pool.apply_async(run_simulation, (input_content, experiment_data_path[j], i, program_start_time,))#notice the comma at the end of the tuple in the second parameter!!!!!!!
            results[j].append(res)
    pool.close()
    pool.join()

    for i in range(2):
        for res in results[i]:
            entry = res.get()
            if i == 0:
                std_summary[entry[0], 0:5] = sample_of_k[entry[0],:]
            std_summary[entry[0], [5+i, 7+i]] = entry[1]
    std_summary_filename = 'std_summary{}.csv'.format(save_order)
    log_filename = 'log{}.txt'.format(save_order)
    save_result(std_summary_filename, std_summary)
    program_end_time = time.time()
    with open(log_filename, 'w') as f:
        f.write('The program takes {:.3f} seconds to complete.'.format(program_end_time - program_start_time))
'''
3.4 read the simulation result
'''

'''
3.5 compare simulation result and experiment data, generating error
'''

'''
3.6 save the error in the array
'''

'''
4. write result into the output file
'''