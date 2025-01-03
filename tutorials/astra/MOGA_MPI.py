# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 09:44:03 2024

@author: lixiangk
modified  30 oct 
"""

# Problem definition with constraints; defined callback to print/save results after each iteration
# Use multiprocessing/threading for parallelization, with each process/threading running the evaluation function individually

import numpy as np
from time import sleep
import matplotlib.pyplot as plt

from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.core.problem import Problem
from pymoo.optimize import minimize
from pymoo.visualization.scatter import Scatter
from pymoo.core.callback import Callback
from pymoo.operators.sampling.rnd import FloatRandomSampling
from pymoo.core.population import Population
from pymoo.core.evaluator import Evaluator
from pymoo.problems.static import StaticProblem
from pymoo.termination.max_gen import MaximumGenerationTermination

import dill

import os
import timeit

import sys
use_mpi = 1

if len(sys.argv)>1:
    use_mpi = sys.argv[1]
print("Use MPI: ", use_mpi)

if use_mpi:
    from mpi4py import MPI
    from mpi4py.futures import MPIPoolExecutor as Pool
    
    # Get the parallel enviroment
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    n_process = size
    print("rank/size: ", rank, ' / ', size)
    
else:
    # Decide the number of processes automatically
    import multiprocessing
    from multiprocessing import Pool
    #from multiprocessing.pool import ThreadPool as Pool
    n_process = multiprocessing.cpu_count()
    print('# of processes: ', n_process)
   
class PrintResultsCallback(Callback):
    def __init__(self):
        super().__init__()
        self.episode = 0

    def notify(self, algorithm):
        self.episode += 1

        if use_mpi:
            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()
            if rank == 0:
                self.notify_0(algorithm)
            # Synchronize all processes
            comm.Barrier()
        else:
            self.notify_0(algorithm)
            
    def notify_0(self, algorithm):
        
        print(f"Generation {self.episode}:")
        
        # save algorithm
        with open(f"checkpoint_{self.episode:d}", "wb") as f:
            dill.dump(algorithm, f)
               
        #print("Objective values:")
        X = algorithm.pop.get("X")
        F = algorithm.pop.get("F")
        G = algorithm.pop.get("G")
        
        #A = np.hstack((X, F))
        A = np.hstack((X, F, G))
        
        fname = f'episode_{self.episode:d}.dat'
        with open(fname, 'w') as f_handle:
            np.savetxt(f_handle, A, fmt = '%14.6E')

        plot = Scatter()
        plot.add(F, edgecolor="red", facecolor="none")
        plot.show()
        figname = str.format('episode_%d.png' % self.episode)
        plot.fig.savefig(figname)
        plt.close(plot.fig)
        
# Define the callback
callback = PrintResultsCallback()


### Backup
import shutil
from datetime import datetime

def _backup_(flags, target_dir, source_dir = '.'):

    # Define the desired format
    ts = datetime.now().strftime("__%d_%m_%Y__%H_%M_%S__")
    target_dir = target_dir+ts
    
    if not os.path.exists(target_dir):
        os.mkdir(target_dir)
        print('Create folder: ', target_dir)
    else:
        print('Folder exits: ', target_dir)
    
    for filename in os.listdir(source_dir):
        for flag in flags:
            if flag in filename:
                # Construct full file path
                full_file_path = os.path.join(source_dir, filename)
                # Copy the file to the target directory
                shutil.copy(full_file_path, target_dir)
    return

### Define test goal/constraint function
# def _my_eval(x):
#     f1 = np.sum(np.power(x, 2))
#     f2 = 1/f1
    
#     f1 = (10-x[1])**2+2*x[0]
#     f2 = -5*x[0]/(10+x[1])
#     #sleep(0.0001)
#     return np.array([f1, f2])
# 
# n_var = 2
# n_obj = 2
# n_ieq = 0
# xl = np.array([1, 1])
# xu = np.array([10, 10])
# n_gen = 5
# pop_size = 128
###

### Define user probelm
from my_eval import *
def _my_eval(x):
    x1 = [x[0], x[1], 45, x[2], 32, x[3], x[4]]
    #return [1, 1]
    return obj_MicroBeam(x1)

n_var = 5
n_obj = 2
n_ieq = 0
xl = np.array([ 6,  0.6, -12, -36, 365])
xu = np.array([20,  1.2,  0, -18, 370])
n_gen = 20
pop_size = 128

load_episode = 0  # this will become the first episode in the output of the resumed run
load_checkpoint = 0 # this output will start from load_checkpoint+1
###

# Define the problem
class MyProblem(Problem):

    def __init__(self, **kwargs):
        super().__init__(n_var=n_var,
                         n_obj=n_obj,
                         n_ieq_constr=n_ieq,
                         xl=xl,
                         xu=xu)
        
    def _evaluate(self, X, out, *args, **kwargs):

        n_var = self.n_var
        n_obj = self.n_obj
        n_ieq = self.n_ieq_constr
        
        # calculate the function values in a parallelized manner and wait until done
        if not use_mpi:
            X_list = X.tolist()
            n_process = multiprocessing.cpu_count()
            with Pool(n_process) as pool:
                results = list(pool.map(_my_eval, X_list))
                
        else: # Use MPI to evaluate the functions

            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()
            size = comm.Get_size()
            n_process = size
            #print("rank/size: ", rank, ' / ', size)
            
            X_array = np.array(X)
            pop_size = len(X_array)
            
            if rank == 0:
                X_array = X_array

                chunk_size = pop_size // n_process

                requests = []
                for i in np.arange(1, size):
                    ss = chunk_size*i
                    x_chunk = X_array[ss:ss+chunk_size]
                    req = comm.isend(x_chunk, dest = i, tag = i)
                    requests.append(req)
                    # req.wait()
                    print(f'Process {i} sent x', x_chunk)

                MPI.Request.waitall(requests)
                x = X_array[0:chunk_size]

            else:
                X_array = None
                req = comm.irecv(source = 0, tag = rank)
                x_chunk = req.wait()
                
                print(f'Process {rank} received x', x_chunk)


            # Synchronize all processes at this point
            comm.Barrier()

            result = np.array([_my_eval(x) for x in x_chunk])
            print(rank, ': ', result)
            
            if rank == 0:

                results = np.vstack((result))
                for i in range(1, size):
                    result = comm.recv(source = i, tag = i)
                    results = np.vstack((results, result))
                
            else:
                #status = MPI.Status()
                comm.send(result, dest = 0, tag = rank)
                
            # Synchronize all processes at this point
            comm.Barrier()

            if rank > 0:
                results = np.empty((pop_size, n_obj+n_ieq), dtype=float)  # Allocate buffer for receiving broadcast
            
            # Broadcast the objective values to all processes
            comm.Bcast(results, root = 0)
            
        out["F"] = np.array(results)[:,0:n_obj]
        if n_ieq>0:
            out["G"] = np.array(results)[:,n_obj:n_obj+n_ieq]
        
# Define the problem
if __name__ == "__main__":

    time1 = timeit.default_timer()
    
    if (use_mpi) or (not use_mpi):
        
        # Define the problem
        problem = MyProblem()
        
        # Define the initial population
        pop = FloatRandomSampling() # This is the default initial sampling method for NSGA2
        if load_episode and not load_checkpoint:
            fname = f'episode_{load_episode:d}.dat'
            ALL = np.loadtxt('..'+os.sep+fname)
            X = ALL[:,0:n_var]
            F = ALL[:,n_var:n_var+n_obj]
            G = ALL[:,n_var+n_obj:n_var+n_obj+n_ieq]
            
            pop = Population.new("X", X)
            pop = Evaluator().eval(StaticProblem(problem, F=F, G=G), pop)
            
        # Define the algorithm    
        algorithm = NSGA2(pop_size=pop_size, sampling=pop)
        if load_checkpoint:
            fname = f"checkpoint_{load_checkpoint:d}"
            with open('..'+os.sep+fname, 'rb') as f:
                algorithm = dill.load(f)
                print("Loaded Checkpoint:", algorithm)
        algorithm.termination = MaximumGenerationTermination(n_gen)
        
        # Always backup old results if existing
        #_backup_(['episode', 'checkpoint'], 'backup')
        
        res = minimize(problem,
                        algorithm,
                        ("n_gen", n_gen),
                        verbose = True,
                        callback = callback,
                        seed=1)
        
        plot = Scatter()
        plot.add(res.F, edgecolor="red", facecolor="none")
        plot.show()

        print('Optimization time:', res.exec_time)

    time2 = timeit.default_timer()
    print('time elapsed: ', time2-time1)
    
