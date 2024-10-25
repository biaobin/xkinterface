import sys
sys.path.append('/afs/ifh.de/group/pitz/data/lixiangk/work/sync/python')
#sys.path.append('\\afs\ifh.de\group\pitz\data\lixiangk\work\sync\python')

from Astra2GenesisSlices import *

def CreateDir(path):
    if not os.path.exists(path):
        os.mkdir(path)
    return

def run_genesis3(x, npart = 8192*4, lambda0 = 100e-6, **kwargs):
    # x: curpeak, curlen, Ntail
    
    curpeak, curlen, Nslice, Nbins, seed = x[:]
    
    partfile, distfile, beamfile = None, None, None
    if len(kwargs)>0:
        if 'partfile' in kwargs.keys():
            partfile = kwargs['partfile']
        if 'beamfile' in kwargs.keys():
            beamfile = kwargs['beamfile']
        if 'distfile' in kwargs.keys():
            distfile = kwargs['distfile']
            
    Nu, lam_u, B = 120, 3e-2, 1.2799839
    lam_s = lambda0
    #K = undulator_parameter(B,lam_u); print ('undulator parameter: ', K)
    K = 3.49
    
    #Q = Qtot # Coulumm

    #FWHM = Q*1e-9/Ipeak # second
    
    sigma_z = curlen
    Q = curpeak*np.sqrt(2*np.pi)*sigma_z/g_c
    #curpeak = Q*1e-9/FWHM; print ('peak current: ', curpeak)
    
    P0 = 17.05 # MeV/c
    #P0 = 22 # MeV/c
    Ek = momentum2kinetic(P0) # MeV
    gamma = kinetic2gamma(Ek); # unitless
    delgam = gamma*0.5e-2
    
    #nslice = Nu+int(np.ceil(FWHM*g_c/lam_s*2)); print ('# of slice: ', nslice-Nu)
    #nslice = 160
    nslice = Nslice
    ntail = 0
    
    emit_x, emit_y = 4e-6, 4e-6
    sigma_x, sigma_y = 1.55e-3, 0.3e-3
    alpha_x, alpha_y = 10.93, 3.25
    
    delz = 0.5
        
    gen = Genesis2()
    
    # set undulator parameters
    gen.set(aw0   = K/np.sqrt(2.),
            awd   = K/np.sqrt(2.),
            nwig  = Nu,
            xlamd = lam_u,
            iertyp= 0,
            delaw = 0,
            iseed = -1)
    
    # set electron beam parameters
    gen.set(gamma0  = gamma,
            delgam  = delgam,
            curpeak = curpeak,
            rxbeam  = sigma_x,
            rybeam  = sigma_y,
            emitx   = emit_x,
            emity   = emit_y,
            alphax  = alpha_x,
            alphay  = alpha_y,
            npart   = npart)
    
    # set particle-loading parameters
    gen.set(ipseed = 1,
            nbins  = Nbins)
    
    # set mesh paremeters
    # Changed dgrid to 0.04->0.02, 31.07.2024
    # Changed nptr from 40 to 100, 02.08.2024
    gen.set(nptr = 40,
            dgrid= 0.04, 
            nscr = 1, 
            nscz = 5)
    
    # set time-dependent parameters/grids
    gen.set(itdp   = 1,
            nslice = nslice,
            ntail  = ntail,
            iotail = 1,
            curlen = sigma_z)
    
    # set radiation parameters
    # Increased ncar to 251->501, 31.07.2024
    # Changed zrayl from 0.03 to 0.3 m, 02.08.2024
    gen.set(xlamds = lam_s, 
            zrayl  = 0.3,
            ncar = 501,
            prad0 = 0) 
    
    # set focusing parameters
    gen.set(quadf = 0, quadd = 0)
    
    # set simulation control parameters
    gen.set(delz  = delz,
            zstop = Nu*lam_u)
    gen.set(zstop = 3.6)
    # set input and ouput control parameters
    gen.set(ffspec = 0,
            ipradi = 0,
            lout = [1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1])

    gen.set(shotnoise = 0)
    
    ipseed = seed
    fname = 'pithz.%d' % (ipseed)
    
    gen.set(outputfile = fname+'.out', maginfile = '../LCLS-I.lat')
    
    if beamfile != None:
        gen.set(beamfile = beamfile)
    if distfile != None:
        gen.set(distfile = distfile)
    if partfile != None:
        gen.set(partfile = partfile)
        
    gen.write(fname+'.in')
    
    #gen.qsub(fname+'.sh', fname+'.in')
    #gen.submit(fname+'.sh', fname+'.in')
    
    #cmd = 'condor_submit '+fname+'.submit'
    #cmd = 'echo '+fname+'.in | genesis 2>&1 | tee '+fname+'.log'
    
    #module add python/3.9 phdf5/1.14.4 szip/2.1.1 mpi/openmpi-x86_64
    cmd = '/afs/ifh.de/group/pitz/data/lixiangk/work/apps/Genesis-1.3-Version3_Alma9/genesis '+fname+'.in 2>&1 | tee '+fname+'.log'
    os.system(cmd)
    
    if seed>2:
        os.remove(partfile)
        pass
    
    #os.chdir('..')

### Scan 
currents = np.linspace(50, 200, 4)
currents = [112]
charges = np.linspace(1.5, 5, 8)
charges = [2]
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
    exit()
combi = [[current, charge, seed]]
print("Input arguments:, ", combi)

#
npart = 4096*8
lambda0 = 100e-6
for v in combi:
    Ipeak, Qtot, seed = v[:] # A, nC, iprime
    
    direc = 'beam_%.0fA_%.1fnC_shotnoise_ffspec0' % (Ipeak, Qtot); print(direc)
    CreateDir(direc)
    os.chdir(direc)
    
    fname = '..'+os.sep+'368A.2809.002.1000' # 2nC
    #fname = '.'+os.sep+'366A.2809.003' # 1nC
    #fname = '..'+os.sep+'ast_1.0nC.2809.1500' # 1nC
    #fname = '..'+os.sep+'ast_1.5nC.2809.1500' # 1.5nC
    #fname = '..'+os.sep+'ast_%.1fnC.2809.1500' % Qtot 
    #fname = '..'+os.sep+'ast.2745.013' # 2.5nC
    #fname = '..'+os.sep+'ast.2745.019' # 4 nC, 17 MeV/c
    #fname = '..'+os.sep+'ast.2745.020' # 4 nC, 22 MeV/c
    
    curpeak, curlen, Nslice, Nbins, outputName = Astra2Genesis3Slices(fname, 
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
    
    run_genesis3(x, lambda0 = lambda0, partfile = outputName, direc = direc, npart = npart)
    
    # Use the distribution from Astra as input
    # distfile = '..'+os.sep+'368A.2809.002.1000@gen2'
    # x = [112, 0.00185, 180, 16, seed]
    # run_genesis3(x, lambda0 = lambda0, distfile = distfile, direc = direc, npart = npart)
    
    os.chdir('..')
    
#exit()