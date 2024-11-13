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

    if isParallel:
        comm = MPI.COMM_WORLD
        size = comm.Get_size()
        rank = comm.Get_rank()
    
    if True:
        print(sys.argv)
        
        grad0 = -1.4
        if len(sys.argv)>1:
            grad0 = np.float(sys.argv[1])
            
        grad0 = -1.6
        
        grad_ = [grad0, -grad0*1.657, grad0*1]
        quadNames_ = ['High1.Q1', 'High1.Q2', 'High1.Q3']
        
        grad_ = []
        quadNames_ = []
        
        grad2 = -grad0*0.7
        grad3 = grad0*0.7
        
        grad2 = 5
        grad3 = -3
        
        # if len(sys.argv)>2:
        #     grad1 = np.float(sys.argv[2])
        # if len(sys.argv)>3:
        #     grad2 = np.float(sys.argv[3])
        # if len(sys.argv)>4:
        #     grad3 = np.float(sys.argv[4])
            
        print(grad_+[grad2, grad3])
        
        
        ### Parameter Setup here
        z0 = 0.443 # position of initial beam
        
        Run0 = 1
        fname0 = 'ast.2931.002'
        print(fname0)
        
        quadNames = ['High1.Q3', 'High1.Q4', 'High1.Q6', 'High1.Q7'] # 1.4
        quadNames = quadNames_+['High1.Q6', 'High1.Q7'] # 1.657
        #quadNames = ['High1.Q1', 'High1.Q2', 'High1.Q3', 'High1.Q4'] # 1.618
        
        quadNames = quadNames_+['High4.Q1', 'High4.Q2'] # 1.657
        
        Zstop = 1
        Zstop = scrns['HIGH4.window']
        ###
        
        
        # ### Back tracking
        # z0 = 0 # position of initial beam
        
        # Run0 = 1
        # fname0 = 'temp.matched'
        # print(fname0)

        # quadNames = ['Back.Q1', 'Back.Q2', 'Back.Q3']
        # Zstop = 2.5
        # #Zstop = scrns['HIGH3.SCR1']
        
        # Zfinal = 5
        # ###
        
        ### Control of scan, not used with the current tuning strategy
        Q2_scanned, Q2_found = 0, 1
        Q3_scanned, Q3_found = 0, 0
        ###
        
        zquad = [quads[name] for name in quadNames] # positions of quads to be tuned
        #zmidd = int((zquad[1]+zquad[2])*5)/10. # output position between the second and third quads
        
        # Parameter scan
        var2 = np.array([0.25, 0.5, 0.75, 0.9, 1.0, 1.1, 1.25, 1.50])*grad2
        combi = np.array([grad_+[v2, 0] for v2 in var2])
                
        args = (quadNames,)
        Distribution = '..'+os.sep+fname0
        kwargs = {'Run':Run0+1, 'Distribution':Distribution, 'Zstart':z0,
                  'Zstop':zquad[-1]+1e-6, 'step':5, 'Screen':[zquad[-1]]}
        
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
                    y = grad_+[grad2, 0]
                    
                    iq = len(grad_)
                    ig = len(y)
                    ix = len(y)+4
                    iy = len(y)+5
                    
                    result = func_wrapper(y)
                    res.append(result[0])
                    
                    nTrials, maxTrials = 0, 5
                    while nTrials<maxTrials:
                        nTrials += 1
                        if len(res)>=2:
                            data = np.array(res)
                            data = data[np.abs(data[:,ig]).argsort()[::]]
    
                            #func = interp1d(1.0-data[:,8]/data[:,9], data[:,2], fill_value = 'extrapolate')
                            func = interp1d(data[:,ix]-data[:,iy], data[:,iq], fill_value = 'extrapolate') # rms size match
                            #func = interp1d(data[:,10]-data[:,11], data[:,1], fill_value = 'extrapolate') # beta match
                            tmp = func([0])[0]                            
                
                        else:
                            
                            if result[0][ig]*np.sign(grad2)>0:
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
                            
                        y = grad_+[grad2, 0]
                        result = func_wrapper(y)
                        res.append(result[0])
                    
                else:
                    r = []
                    for x in combi:
                        result = func_wrapper(x)
                        r.append(result)
            
            fig, ax = plt.subplots()
            ax.plot(data[:,iq], data[:,ix], '-*')
            ax.plot(data[:,iq], data[:,iy], '-o')
            ax.grid()
            ax.set_xlabel(r'$G$ (T/m)')
            ax.set_ylabel(r'rms size (mm)')
            fig.savefig('rms-size-vs-grad2@%.2fT_m.eps' % grad2)

        # else:
        #     grad1 = None #comm.bcast(grad1, root = 0)
        #     x = None #comm.bcast(x, root = 0)
            
        if isParallel:
            grad1 = comm.bcast(grad1, root = 0)
            grad2 = comm.bcast(grad2, root = 0)
            x = comm.bcast(y, root = 0)
        
        # Find Q3
        fname1 = fname0
        Distribution = '..'+os.sep+fname1
        
        # Parameter scan
        var3 = np.linspace(0.5, 2.0, 8)*grad3
        var3 = np.array([0.25, 0.5, 0.75, 0.9, 1.0, 1.1, 1.25, 1.50])*grad3
        combi = np.array([grad_+[grad2, v3] for v3 in var3])
        
        kwargs = {'Run':Run0+2, 'Distribution':Distribution, 'Zstart':z0,
                  'Zstop':Zstop, 'step':2, 'Screen':[Zstop],
                  'zmin':zquad[-2]-z0}
        
        def func_wrapper(x):
            r = 0
            r = ObjOcelot(x, zquad, **kwargs, obj_fun = get_obj3)
            return [r]
        
        if not Q3_scanned:
            if isParallel:
                m = int(math.ceil(float(len(combi)) / size))
                x_chunk = combi[rank*m:(rank+1)*m]
                r_chunk = list(map(func_wrapper, x_chunk))
    
                r = comm.allreduce(r_chunk)
            else:
                if 0:
                    def obj_func(x):
                        y = grad_+[grad2, 0]+[x[0], 0]
                        result = func_wrapper(y); print(result)
                        obj = result[0][3]**2
                        return obj
                    from scipy.optimize import minimize
                    grad2 = -grad1
                    init = np.array([[grad2], [grad2*1.5]])
                    res = minimize(obj_func, [grad2], method='Nelder-Mead', 
                                   tol=1e-6, options = {'maxiter':100, 'initial_simplex':init})
                elif 1:
                    res = []
                    #grad3 = grad1
                    y = grad_+[grad2, grad3]
                    iq = len(y)-1
                    
                    result = func_wrapper(y)
                    res.append(result[0])
                    
                    while 1:
                        if len(res)>=2:
                            data = np.array(res)
                            
                            select = (data[:,ig]<1e4); data = data[select]
                            data = data[np.abs(data[:,ig]).argsort()[::]]
    
                            #func = interp1d(data[:,3], data[:,2], fill_value = 'extrapolate') # rms size match
                            func = interp1d(data[:,ix]-data[:,iy], data[:,iq], fill_value = 'extrapolate') # rms size match
                            #func = interp1d(data[:,10]-data[:,11], data[:,2], fill_value = 'extrapolate') # beta match
                            tmp = func([0])[0]                            
                
                        else:
                            if result[0][i]*np.sign(grad3)>0:
                                tmp = grad3*1.25
                            else:
                                tmp = grad3*0.75
                    
                        if np.abs(1-tmp/grad3)<5e-3:
                                print(('Grad3 = %.6f T/m' % grad3))
                                break
                        else:
                            direc = CreateFolderName(y, flag = 'transport')
                            shutil.rmtree(direc)
                            grad3 = tmp
                            
                        y = grad_+[grad2, grad3]
                        result = func_wrapper(y)
                        res.append(result[0])
                    
                else:
                    r = []
                    for x in combi:
                        result = func_wrapper(x)
                        r.append(result)
                        
            fig, ax = plt.subplots()
            ax.plot(data[:,iq], data[:,ix], '-*')
            ax.plot(data[:,iq], data[:,iy], '-o')
            ax.grid()
            ax.set_xlabel(r'$G$ (T/m)')
            ax.set_ylabel(r'rms size (mm)')
            fig.savefig('rms-size-vs-grad3@%.2fT_m.eps' % grad3)

            with open('results.dat','a') as f_handle:
                np.savetxt(f_handle,np.atleast_2d(y),fmt='%14.6E')
    
    #plt.close('all')
    direc = CreateFolderName(y, flag = 'transport')
    
    y = [-2.000000E+00,   2.800000E+00,  -1.332976E+00,   5.023331E-01]
    y = [-1.700000E+00,   2.380000E+00,  -1.050372E+00,   3.258825E-01]
    #y = [-1.600000E+00,   2.240000E+00,  -9.668652E-01,   2.856148E-01]
    y = [-1.650000E+00,   2.310000E+00,  -1.007444E+00,   3.038712E-01]
    
    
    #y = [-1.000000E+00,   1.657000E+00,  -1.109959E+00,   3.982906E-01]
    #Zfinal = 28.09
    if Zfinal > Zstop:
        #shutil.rmtree(direc)
    
        Zstop = Zfinal
        def func_wrapper(x):
            r = ObjOcelot(x, zquad, Distribution = Distribution, Zstart = z0,
                        Zstop = Zstop, obj_fun = get_obj2, save = True, Distribution1 = 'ast.2809.002')
            return [r]

        y = grad_+[grad2, grad3]
        result = func_wrapper(y)

    print('Setting found: ', y)
    # make plot
    # make_plot(y, z0 = z0)
