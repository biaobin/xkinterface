# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 06:27:51 2023

@author: lixiangk
"""

from interface import *

def post_job(Ipeak, Qtot):
    direc = 'beam_%.0fA_%.1fnC' % (Ipeak, Qtot); print(direc)
    os.chdir(direc)
    
    os.system('cp ../post.py ./post.py')
    #job = QsubJob(command = 'python')
    #job.create(jobName = 'tmp', inputName = '$Py/postGenesis_qsub.py', submit = True)
    
    job = CondorJob(command = 'python')
    job.create(jobName = 'temp.sh', inputName = '../post.py', submit = True)
    
    os.chdir('..')
    
    return


### Scan 
currents = [50, 100, 150, 200]
currents = [112]
charges = np.linspace(1.5, 5, 8)
charges = [2.]
seeds = np.arange(100)+1
#seeds = [1]

combi = np.array([[v1, v2, v3] for v1 in currents for v2 in charges for v3 in seeds])
#combi = [[78, 1], [100, 1.5]]
#combi = [[112, 2]]

version = 4
for v in combi:
    
    Ipeak, Qtot, seed = v[:] # A, nC, iprime
    
    if version == 3:
        inputName = 'run_job3.py'
        cmd1 = 'module purge && module add python/3.9 mpi/openmpi-x86_64 phdf5/1.14.4 szip/2.1.1'
    elif version == 4:
        inputName = 'run_job4_.py'
        cmd1 = 'module purge && module add python/3.9 mpi/openmpi-x86_64 phdf5/1.14.4-2 szip/2.1.1'
        
    jobName = 'beam_%.0fA_%.1fnC_%d_quiet2' % (Ipeak, Qtot, seed)
    
    ### simulation
    job = CondorJob(command = 'python')
    job.create(jobName+'.sh', inputName, args = v, cmd1 = cmd1)
    
#exit()

#%%
