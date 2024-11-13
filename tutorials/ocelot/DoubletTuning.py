import sys

from scipy.interpolate import interp1d, interp2d
import matplotlib as mpl
import matplotlib.pyplot as plt
#plt.switch_backend('agg')

import timeit, math
import shutil, os

# Get the directory of the currently running script
scriptdir = os.path.dirname(os.path.abspath(__file__))

# Add this directory to sys.path if it's not already there
if scriptdir not in sys.path:
    sys.path.insert(0, scriptdir)

# Now you can import the module as usual
from ObjOcelot import *

#%% Use or not use MPI
rank = 0
isParallel = False
try:
    from mpi4py import MPI
    isParallel = True
except Exception as err:
    print(err)
    
#%% Go the the working folder
workdir = r'\\afs\ifh.de\group\pitz\data\lixiangk\sim3\2024\Biaobin'
os.chdir(workdir)

if __name__ == "__main__":
    
    if True:
        print(sys.argv)
        
        grad1 = 3 # gradient of the first quad
        if len(sys.argv)>1:
            grad1 = np.float(sys.argv[1])
        grad2 = -grad1*1.3
        
        grad3 =  grad1*0
        if len(sys.argv)>2:
            grad2 = np.float(sys.argv[2])
        if len(sys.argv)>3:
            grad3 = np.float(sys.argv[3])
        print(grad1, grad2, grad3)
        
        ### Parameter Setup here
        z0 = 4.5 # position of initial beam
        
        Run0 = 1
        fname0 = 'ast.0450.002'
        print(fname0)
        
        quadNames = ['High1.Q1', 'High1.Q2', 'High1.Q3']
        Zstop = (quads[quadNames[-1]]+quads[quadNames[-2]])/2
        #Zstop = pitz['High1.Scr2']
        
        Zfinal = pitz['High1.Scr3']
        ###
        
        ### Control of scan, not used with the current tuning strategy
        Q2_scanned, Q2_found = 0, 1
        Q3_scanned, Q3_found = 0, 0
        ###
        
        zquad = [quads[name] for name in quadNames] # positions of quads to be tuned
        #zmidd = int((zquad[1]+zquad[2])*5)/10. # output position between the second and third quads
        
        # Parameter scan
        var2 = np.array([0.25, 0.5, 0.75, 0.9, 1.0, 1.1, 1.25, 1.50])*grad2
        combi = np.array([[grad1, v2, 0] for v2 in var2])
        
        args = (quadNames,)
        Distribution = '..'+os.sep+fname0
        kwargs = {'Run':Run0+1, 'Distribution':Distribution, 'Zstart':z0,
                  'Zstop':Zstop, 'step':5, 'Screen':[zquad[-1]]}
        
        def func_wrapper(x):
            r = 0
            r = ObjOcelot(x, zquad[:-1], **kwargs, obj_fun = get_obj2)
            return [r]

        # parallel mapping
        if not Q2_scanned:
            if isParallel:
                
                m = int(math.ceil(float(len(combi))/size))
                x_chunk = combi[rank*m:(rank+1)*m]
                r_chunk = list(map(func_wrapper, x_chunk))

                # collect results
                r = comm.allreduce(r_chunk)
                
            else:
                if 1:
                    res = []
                    #grad2 = -grad1
                    y = [grad1, grad2, 0]
                    result = func_wrapper(y)
                    res.append(result[0])
                    
                    nTrials, maxTrials = 0, 5
                    while nTrials<maxTrials:
                        nTrials += 1
                        if len(res)>=2:
                            data = np.array(res)
                            data = data[np.abs(data[:,1]).argsort()[::]]
    
                            func = interp1d(1.0-data[:,7]/data[:,8], data[:,1], fill_value = 'extrapolate')
                            #func = interp1d(data[:,7]-data[:,8], data[:,1], fill_value = 'extrapolate') # rms size match
                            #func = interp1d(data[:,10]-data[:,11], data[:,1], fill_value = 'extrapolate') # beta match
                            tmp = func([0])[0]                            
                
                        else:
                            if result[0][3]*np.sign(grad1)>0:
                                tmp = grad2*0.75
                            else:
                                tmp = grad2*1.25
                    
                        direc = ''
                        for i, xi in enumerate(y):
                            direc += str.format('Q%d-%.2fT_m-' %  (i+1, xi))
                        direc = direc[:-1]
                        shutil.rmtree(direc)
                        
                        if np.abs(1-tmp/grad2)<1e-3:
                            print(('Grad2 = %.4f T/m' % grad2))
                            break
                        else:
                            grad2 = tmp
                            
                        y = [grad1, grad2, 0]
                        result = func_wrapper(y)
                        res.append(result[0])
                    
                else:
                    r = []
                    for x in combi:
                        result = func_wrapper(x)
                        r.append(result)
            
            fig, ax = plt.subplots()
            ax.plot(data[:,1], data[:,7], '-*')
            ax.plot(data[:,1], data[:,8], '-o')
            ax.grid()
            ax.set_xlabel(r'$G$ (T/m)')
            ax.set_ylabel(r'rms size (mm)')
            fig.savefig('rms-size-vs-grad2@%.2fT_m.eps' % grad1)

        # else:
        #     grad1 = None #comm.bcast(grad1, root = 0)
        #     x = None #comm.bcast(x, root = 0)
            
        if isParallel:
            grad1 = comm.bcast(grad1, root = 0)
            grad2 = comm.bcast(grad2, root = 0)
            x = comm.bcast(y, root = 0)
            
        
    if Zfinal > 0:
        #shutil.rmtree(direc)
    
        Zstop = Zfinal
        def func_wrapper(x):
            r = ObjOcelot(x, zquad, Distribution = Distribution, Zstart = z0,
                        Zstop = Zstop, obj_fun = get_obj2)
            return [r]

        y = [grad1, grad2, grad3]; 
        result = func_wrapper(y)

    print('Setting found: ', y)
    # make plot
    # make_plot(y, z0 = z0)
