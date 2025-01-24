import sys
sys.path.append('/afs/ifh.de/group/pitz/data/lixiangk/work/sync/python')
#sys.path.append('\\afs\ifh.de\group\pitz\data\lixiangk\work\sync\python')

from xkinterface.tutorials.genesis13.Astra2GenesisSlices import *

def CreateDir(path):
    if not os.path.exists(path):
        os.mkdir(path)
    return

def run_genesis4(x, npart = 8192*4, lambda0 = 100e-6, **kwargs):
    # x: curpeak, curlen, Ntail
    
    curpeak, curlen, Nslice, Nbins, seed = x[:]
    
    partfile, distfile, beamfile = None, None, None
    if len(kwargs)>0:
        if 'partfile' in kwargs.keys():
            partfile = kwargs['partfile']
            
    Nu, lam_u, B = 120, 3e-2, 1.2799839
    lam_s = lambda0
    #K = undulator_parameter(B,lam_u); print ('undulator parameter: ', K)
    K = 3.49
    
    P0 = 17.05 # MeV/c
    Ek = momentum2kinetic(P0) # MeV
    gamma0 = M2G(P0)
    
    lambda0 = 100e-6
    lam_u = 30e-3
    delz = lam_u/2
    
    rootname = 'g4.%d' % (seed)
    
    setup = Namelist('setup', 
                 rootname = rootname,
                 lattice = '..'+os.sep+'gen4lat.lat',
                 beamline = 'THzBL',
                 gamma0 = gamma0,
                 lambda0 = lambda0,
                 delz = delz,
                 seed = seed,
                 npart = npart,
                 nbins = Nbins,
                 one4one = False,
                 shotnoise = True)
    
    time = Namelist('time',
                s0 = 0,
                slen = lambda0*Nslice,
                sample = 1,
                time = True)
    
    profile_gauss = Namelist('profile_gauss', 
                          label = 'beampro',
                          c0 = 112.0,
                          s0 = 0,
                          sig = 2.2136e-3)
    
    # profile_step = Namelist('profile_step', 
    #                      label = 'beampro',
    #                      c0 = 150.0,
    #                      s_start = 0e-3,
    #                      s_end = 5e-3)
    
    beam = Namelist('beam',
                gamma = gamma0,
                delgam = gamma0*0.5e-2,
                current = '@beampro',
                ex = 4e-6,
                ey = 4e-6,
                betax = 20,
                betay = 0.75,
                alphax = 10.92,
                alphay = 3.25,
                bunch = 0,
                emod = 0)
    
    # importdist = Namelist('importdistribution',
    #                   file = 'test_dist.h5',
    #                   charge = 1e-9,
    #                   slicewidth = 0.01)
    
    # Note that if `importbeam` is used, the namelist `beam` should be deleted
    # e.g., partfile = 'test.0.par.h5', which provides slice-wise particle distributions
    importbeam = Namelist('importbeam',
                          file = partfile,
                          time = True)
    
    efield = Namelist('efield',
                      longrange = 1,
                      rmax = 0.01,
                      nz = 5,
                      nphi = 1,
                      ngrid = 100)
    
    field = Namelist('field',
                     power = 0,
                     phase = 0,
                     waist_size = 0.3,
                     waist_pos = 0,
                     dgrid = 40e-3,
                     ngrid = 501)
    
    track = Namelist('track',
                 output_step = 1,
                 field_dump_step = 0,
                 beam_dump_step = 0)
    
    g4 = Genesis4(setup, time, importbeam, efield, field, track)
    #print(g4.output)
    
    fname = rootname
    g4.write(fname+'.in')
    
    #module add python/3.9 phdf5/1.14.4-2 szip/2.1.1 mpi/openmpi-x86_64
    #cmd = '/afs/ifh.de/group/pitz/data/lixiangk/work/apps/genesis/4.6.6/bin/genesis4 '+fname+'.in 2>&1 | tee '+fname+'.log'
    #os.system(cmd)
    
    if seed>2:
        #os.remove(partfile)
        pass
    
    #os.chdir('..')

### Scan 
if __name__=="main":
    currents = np.linspace(50, 200, 4)
    currents = [112]
    charges = np.linspace(1.5, 5, 8)
    charges = [2.0]
    seeds = np.arange(20)+1
    #seeds = [1]
    
    combi = np.array([[v1, v2, v3] for v1 in currents for v2 in charges for v3 in seeds])
    
    ### Or read input arguments from shell
    if len(sys.argv)>3:
        current = float(sys.argv[1])
        charge = float(sys.argv[2])
        seed = int(sys.argv[3])
    else:
        print('No enough input arguments! Exit!')
        sys.exit()
    combi = [[current, charge, seed]]
    print("Input arguments:, ", combi)
    
    #
    npart = 4096*8
    lambda0 = 100e-6
    for v in combi:
        Ipeak, Qtot, seed = v[:] # A, nC, iprime
        
        direc = 'beam_%.0fA_%.1fnC_shotnoise_v4_gaussian' % (Ipeak, Qtot); print(direc)
        CreateDir(direc)
        os.chdir(direc)
        
        fname = '..'+os.sep+'368A.2809.002.1000' # 2 nC
        fname = '..'+os.sep+'oce_2.5nC.2809.001' # 2.5 nC
        fname = '..'+os.sep+'beam_112A_2.0nC.ini' # 2.0 nC, perfect Gaussian profile
        curpeak, curlen, Nslice, Nbins, outputName = Astra2Genesis4Slices(fname, 
                                                      outputName = 'scan', 
                                                      seed = seed, 
                                                      npart = npart,
                                                      lambda0 = lambda0,
                                                      Qscale = 1, 
                                                      zscale = 1,
                                                      bunch = 0, 
                                                      useHammersley = 0)
        x = [curpeak, curlen, Nslice, Nbins, seed]
        #x = [112, 0.00185, 187, 16, seed]
        
        outputName = 'scan'+str.format('.%d.out.par.h5' % seed)
        
        run_genesis4(x, lambda0 = lambda0, partfile = outputName, direc = direc, npart = npart)
        
        # Use the distribution from Astra as input
        # distfile = '..'+os.sep+'368A.2809.002.1000@gen2'
        # x = [112, 0.00185, 180, 16, seed]
        # run_genesis3(x, lambda0 = lambda0, distfile = distfile, direc = direc, npart = npart)
        
        os.chdir('..')
        
    #sys.exit()
