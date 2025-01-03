import sys
#if sys.platform == "linux" or sys.platform == "linux2":
#    sys.path.append(r'/afs/ifh.de/group/pitz/data/lixiangk/work/sync/python')
#elif sys.platform == "win32":
#    sys.path.append(r'\\afs\ifh.de\group\pitz\data\lixiangk\work\sync\python')
#elif sys.platform == "darwin":
#    print("OS X")
#else:
#    print("Unknown platform!")

#from PITZ import *
from interface import *


def obj_PhotoInjector(x, *args, farm = False, submit = False, get_folder = False, **kwargs):
    
    '''
    Simulation from photocathode to EMSY1 at 5.277 m 
    Variables that could be optimized: laser BSA diameter and FWHM length, gun and booster gradients and phases, solenoid current
    Parameters:
      x: an array or list of the variables to be optimized, from beginning: FWHM, BSA, MaxE_gun, phi_gun, MaxE_booster, phi_booster, Imain
      kwargs: can be used to define other parameters, such as bunch charge, e.g., "Q_total = -1" will define a bunch charge of 1 nC.
              The variable here can be anything defined in Astra manual.
    Returns:
      if: the quantity that is envisoned as the "energy" of the sample
    '''

    Ipart = 250000
    Run = 1
    Q_total = -2 # nC
    
    ### Laser
    FWHM = x[0]*1e-3 # ns in Astra
    #FWHM = 250e-6 # for Pharos shortest pulse
    BSA = x[1]

    # Uniform
    sigma_x = sigma_y = BSA/4.
    C_sigma_x = C_sigma_y = 0
    # Or Gaussian truncated
    # sigma_x = sigma_y = 0.96 
    # C_sigma_x = C_sigma_y = BSA/2./sigma_x

    ### Gun and booster
    phi_gun, phi_booster = x[3], x[5]

    # Define Egun and Ebooster by the beam momentum
    MaxE_gun = get_MaxE_gun(phi_gun, 6.3)
    MaxE_booster = get_MaxE_booster(MaxE_gun, phi_gun, phi_booster, 17)
    # Or define them by user
    # MaxE_gun, MaxE_booster = x[2], x[4]

    ### Solenoid
    Imain = x[6]
    MaxB = I2B(Imain)
    
    ### Path to field maps, by default `path_to_interface/PITZsim/field_maps'
    # Users can also define their own path
    #field_maps = rootdir+os.sep+'sync'+os.sep+'field-maps'
    
    Distribution = 'BSA-%.1fmm-FWHM-%.1fps.ini' % (BSA, FWHM*1e3)

    Zstart, Zstop = 0, 5.28
    
    if len(kwargs)>0:
        keys = kwargs.keys()
        if 'Ipart' in keys:
            Ipart = kwargs['Ipart']
        if 'Run' in keys:
            Run = kwargs['Run']
        if 'Distribution' in keys:
            Distribution = kwargs['Distribution']
        if 'Q_total' in keys:
            Q_total = kwargs['Q_total']
        if 'Zstart' in keys:
            Zstart = kwargs['Zstart']
        if 'Zstop' in keys:
            Zstop = kwargs['Zstop']

    generator = Generator1(FNAME = Distribution, IPart = Ipart, Species = 'electrons', Q_total = Q_total,
                           Ref_Ekin = 0.0e-6, LE = 0.55e-3*1, dist_pz = 'i',
                           Dist_z = 'g', sig_clock = FWHM/2.355, Cathode = True,
                           Dist_x = '2', sig_x = sigma_x, Dist_px = 'g', Nemit_x = 0,
                           Dist_y = '2', sig_y = sigma_y, Dist_py = 'g', Nemit_y = 0,
                           C_sig_x = C_sigma_x, C_sig_y = C_sigma_y) # Gaussian by default
    #generator.set(Ref_Ekin = 1e-6, LE = 0, dist_pz = 'g', Nemit_x = 0.17, Nemit_y = 0.18)
    #generator.set(Dist_x = 'r', sig_x = BSA/4.0, Dist_y = 'r', sig_y = BSA/4.0) # Transverse uniform

    #Rt = 2.5e-3
    #generator.set(Dist_z = 'p', Lt = FWHM, Rt = Rt) # Temporal flattop
    
    newrun = Module('Newrun', Run = Run, Head = 'PITZ beam line simulation', Distribution = Distribution, CathodeS = True,
                    Auto_Phase = True, Track_All = True, check_ref_part = False, Lprompt = False, Max_step=200000)
    #newrun.set(Xoff = 0.5, Yoff = 0.5)
    #newrun.set(Qbunch = Q_total)
    #newrun.set(Run = 1, XYrms = sigma_x)
    
    charge = Module('Charge', LSPCH = True, Lmirror = True, Nrad = 40, Nlong_in = 50, N_min = 50, Max_scale = 0.05, Max_count = 20)
    charge.set(N_min= 100)
    #charge.set(L2D_3D = True, Z_TRANS = 5, NXF = 32, NYF = 32, NZF = 32, min_grid_trans = 0.03e-3)
    #charge.set(LSPCH3D = True, NXF = 32, NYF = 32, NZF = 32)
    charge.set(LSPCH = False)
    
    cavity = Module('Cavity', LEfield = True, File_Efield = 
                    [field_maps+os.sep+gun_profile, field_maps+os.sep+'CDS14_15mm.txt'],
                    MaxE = [MaxE_gun, MaxE_booster], C_pos = [0., 2.675], Nue = [1.3, 1.3], Phi = [phi_gun, phi_booster])
    
    soleno = Module('Solenoid', LBfield = True, File_Bfield = [field_maps+os.sep+'gunsolenoidsPITZ.txt'],
                    MaxB = [-MaxB], S_pos = [0.], S_xrot = [0*1e-3], S_yrot = [0])
    
    #Zstart, Zstop = 0, 5.28
    Zemit = int((Zstop-Zstart)/0.01)
    output = Module('Output', Zstart = Zstart, Zstop = Zstop, Zemit = Zemit, Zphase = 1, RefS = True, EmitS = True,
                    PhaseS = True, TrackS = False, LandFS = True, C_EmitS = True, LPROJECT_EMIT = True,
                    LOCAL_EMIT = False, Screen = [5.28])
    #output.set(Zstop = 0.5, Zemit = 50) # emission
    
    apertu = Module('Aperture', LApert = True, File_Aperture = [field_maps+os.sep+'app.txt'])
    
    #Qsol = 0 #-0.0034*x[-2]/365.0
    #Qsol_zrot = 0 # -0.4363 #x[-1]/180.0*np.pi
    #quadru = Module('Quadrupole', LQuad = True, Q_type = ['../Qsol.data'], Q_pos = [0], Q_grad = [Qsol], Q_zrot = [Qsol_zrot])
    
    astra = Astra()
    astra.add_modules([newrun, charge, cavity, soleno, output])
    
    direc = str.format('L-%.2fps-D-%.2fmm-E1-%.2fMV_m-phi1-%.2fdeg-E2-%.2fMV_m-phi2-%.2fdeg-I-%.2fA' %
                       (FWHM*1e3, BSA, MaxE_gun, phi_gun, MaxE_booster, phi_booster, Imain))
    #direc = str.format('D-%.2fmm-I-%.0fA' % (x[1], x[-1]))

    ### If asking for folder name only, return the folder name for post-processing
    if get_folder:
        return direc
    
    ### For farm jobs, create the batch files for job submission
    if farm:
        CreateFolder(direc)
        #os.system('mkdir -p '+direc)

        #os.chdir(direc)
        job_name = direc+'-%03d' % Run
        #job_name = 'myjob.'+('BSA-%.2fmm-I-%.2fA' % (x[0], x[3]))
        gen_name = 'gen' #+`Run`
        ast_name = 'ast' #+`Run`
    
        generator.write(direc+os.sep+gen_name+'.in')
        astra.write(direc+os.sep+ast_name+'.in')
        #astra.qsub(job_name, ast_name, gen_name, direc)
        astra.submit(job_name, ast_name, gen_name, direc, submit = submit)

        return
    
    ### If not farming, continue
    cwd = os.getcwd()
    
    CreateFolder(direc)
    #os.system('mkdir -p '+direc)
    os.chdir(direc)
    
    generator.write('gen.in')
    astra.write('ast.in')
    
    os.system('generator gen.in 2>&1 | tee gen.log')
    os.system('astra ast.in 2>&1 | tee ast.log')
    
    z0 = Zstop
    try:
        fname = 'ast.%04d.001' % (z0*10)
        if not os.path.isfile(fname):
            fname = 'ast.%04d.001' % (z0*100)
            if not os.path.isfile(fname):
                fname = 'ast.%04d.001' % (z0*1000)
                
        dist = np.loadtxt(fname)
            
        dist[1:,2] += dist[0,2]; dist[1:,5] += dist[0,5]
        diag = BeamDiagnostics(dist = dist)

        if diag.loss_cathode/diag.NoP < 0.05:
            f1 = diag.nemit_x # longitudinal core emittance 
            f2 = diag.std_Ek2 # Bunch length or high-order energy spread
            f2 = np.abs(diag.cor_Ekin*diag.std_z+108)/108
        else:
            f1, f2 = 999, 999

        with open('../results.dat', 'a') as f_handle:
            np.savetxt(f_handle, np.atleast_2d(list(x)+list(diag.x)), fmt = '%14.6E')

    except Exception as err:
        print(err)
        f1, f2 = 999, 999
    os.chdir(cwd)
    return np.array([f1, f2])

import difflib

def find_closest_folder(target_folder_name, search_directory = '.'):
    # Get a list of folders in the specified directory
    folders = [f for f in os.listdir(search_directory) if os.path.isdir(os.path.join(search_directory, f))]
    
    # Use difflib to find the closest match
    closest_match = difflib.get_close_matches(target_folder_name, folders, n=1)
    
    # Return the closest match or None if no match is found
    return closest_match[0] if closest_match else None

def post_PhotoInjector(x, *args, Zstop = 5.28, **kwargs):
    
    '''
    Parameters:
      x: an array or list of the variables to be optimized
    Returns:
      energy: the quantity that is envisoned as the "energy" of the sample
    '''
    
    cwd = os.getcwd()

    direc = obj_PhotoInjector(x, *args, get_folder = True, **kwargs)
    direc = find_closest_folder(direc);print(direc)
    
    os.chdir(direc)
    
    z0 = Zstop
    try:
        fname = 'ast.%04d.001' % (z0*10)
        if not os.path.isfile(fname):
            fname = 'ast.%04d.001' % (z0*100)
            if not os.path.isfile(fname):
                fname = 'ast.%04d.001' % (z0*1000)
                
        dist = np.loadtxt(fname)
        dist[1:,2] += dist[0,2]
        dist[1:,5] += dist[0,5]
        
        diag = BeamDiagnostics(dist = dist)
        
        with open('../episode_18_post_.dat', 'a') as f_handle:
            np.savetxt(f_handle, np.atleast_2d(list(x)+list(diag.x)), fmt = '%14.6E')

    except Exception as err:
        print(err)

    os.chdir(cwd)
    return cwd

if __name__ == "__main__":
    workdir = 'test7'
    os.chdir(workdir)
    
    it = 10
    filename = f'episode_{it:d}.dat'
    data = np.loadtxt(filename)

    n_var = 5
    for i, x in enumerate(data[:10]):
        print(i, x)
        x1 = [x[0], x[1], 57.55, x[2], 13.38, x[3], x[4]]
        post_PhotoInjector(x1)
