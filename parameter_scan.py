# -*- coding: utf-8 -*-
import numpy as np
import subprocess
from multiprocessing import Pool
'''
0. global constant definitions
'''
rk4_path = './rk4'
experiment_data_path = '/home/u24658ly/Dropbox/lab/lab_data_in_paper1.csv'
number_of_process = 8
number_of_iteration = 10000
concentration_of_NO = 1.53e-6
concentration_of_sGC = 0.47e-6
concentration_of_6C = 0
concentration_of_5C = 0
length_of_optical_cell = 1
start_time = 0
end_time = 1000
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

'''
input_content: string
input_filename: string
result_filename: string

create input file for rk4 executable and call rk4 with argument input_filename, result_filename.

'''
def run_simulation(input_content, input_filename, result_filename):
    with open(input_filename, 'w') as f:
        f.write(input_content)   
    subprocess.run([rk4_path, input_filename, result_filename])

'''
2. read in experiment data
'''
#experiment_data = np.genfromtxt(experiment_data_path, skip_header=1, delimiter=',')

'''
3.1 random sampling based on uniform distribution
'''
print('Program starts...\n')
#with Pool(processes=number_of_process) as pool:
'''
The 'with' statement would not work. Why???
'''
pool = Pool(number_of_process)
for i in range(16):
    sample = np.random.uniform(lower_bound_of_k, upper_bound_of_k, size=5)
    input_content = '{:.2}\n{:.2}\n{:.2}\n{:.2}\n{:.2}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n'.format(
                                                                                    sample[0],
                                                                                    sample[1],
                                                                                    sample[2],
                                                                                    sample[3],
                                                                                    sample[4], 
                                                                                    length_of_optical_cell,
                                                                                    concentration_of_NO,
                                                                                    concentration_of_sGC,
                                                                                    concentration_of_6C,
                                                                                    concentration_of_5C,
                                                                                    start_time, end_time,
                                                                                    eps,
                                                                                    h1,
                                                                                    hmin)
    input_filename = 'input{}.txt'.format(i)
    result_filename = 'result{}.csv'.format(i)
    pool.apply_async(run_simulation, (input_content, input_filename, result_filename,))#notice the comma at the end of the tuple in the second parameter!!!!!!!
pool.close()
pool.join()

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