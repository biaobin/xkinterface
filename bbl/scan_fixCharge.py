from xkinterface.interface.Namelists import Namelist
from xkinterface.interface.Genesis13 import Genesis4
from math import *
import numpy as np
import os
import shutil

def gen_impzin(profile='parabolic',Lbuncht=6e-12,tol_charge=2e-9, Np=1e6, folderName='.'):    
    
    if profile=='parabolic':
        # sigz=c*sigt
        # Q=13.3422*Ipeak*sigz*1e-9
        sigz = Lbuncht/6*c
        Ipeak = tol_charge/(13.3422*sigz*1e-9)  #check, *1e-9?
        distribution_type = 46        
        
    elif profile=='flattop':
        Ipeak=tol_charge/(0.9*Lbuncht)
        distribution_type=49
        sigz = Lbuncht*c/2  #half bunch length [m]
        
    else:
        print("Error...")
        sys.exit(0)
    
    Ekin = E0-0.511e6
    
    control = Namelist('control',
        default_order = 2,
        steps=5,
        lsc=0,
        zwake=0,
        csr=0,
        tsc=0,
        trwake=0,
        core_num_T = 2,
        core_num_L = 2,
        meshx = 32,
        meshy = 32,
        meshz = 64,
        kinetic_energy =Ekin,
        freq_rf_scale = 1.3e9,
        slice_bin=100,
        sample_out=1e5, 
        pipe_radius=20e-3,
        error = 1)
    
    beam = Namelist('beam',
        mass = 0.511001e6,
        charge = -1.0,
        distribution_type = distribution_type,
        Np = Np,
        total_charge = tol_charge,
        emit_nx=enx, beta_x=betax, alpha_x=alphax,
        emit_ny=eny, beta_y=betay, alpha_y=alphay,
        sigz=sigz,sigE=sigE)           
                       
    g4 = Genesis4(control, beam)
    #print(g4.output)
    
    g4.output = g4.output.lower()
    g4.write(folderName+'/lte.impz')
    
    with open(folderName+'/lte.impz',"a",encoding="utf-8") as f:
        f.write("&lattice \n")
        f.write(f"match: match2twiss, betax={betax}, alphax={alphax},betay={betay},alphay={alphay} \n")
        
        if one4one == True:
            f.write("w1: watch, filename_id=1001, sample_freq=1, coord_conv='GenesisV4',coord_info=1 \n")
        else:
            f.write("w1: watch, filename_id=1001, sample_freq=1, coord_conv='Astra',coord_info=1 \n")
            
        f.write("line: line=(match,w1) \n")
        f.write("&end")
    
    if profile=="flattop":
        gen_flattopProfile(folderName=folderName,plot=False)
        
import matplotlib.pyplot as plt
def gen_flattopProfile(folderName='.',plot=False):
    
    clight = 2.998e8
    a = 0.1*20e-12 *clight
    b = 0.8*20e-12*clight
    c = 0.1*20e-12 *clight

    I = 120
    dx = b/200

    x1 = np.arange(0,a,dx)
    x2 = np.arange(a,a+b,dx)
    x3 = np.arange(a+b,a+b+c+dx,dx)

    fx1 = I/a*x1
    fx2 = I *np.ones(len(x2))
    fx3 = -I/c*x3 + I*(1+(a+b)/c)

    # concatenate x and fx
    x = np.concatenate((x1, x2, x3))
    fx = np.concatenate((fx1, fx2, fx3))

    Fx1 = 0.5*I/a*x1**2
    Fx2 = 0.5*I*a +(x2-a)*I

    x3_shift = x3 - (a + b)
    Fx3 = 0.5*I*a + I*b + (-0.5 *I/c *x3_shift**2 +I *x3_shift)
    # Fx3 = 0.5*I*a +b*I -0.5*I/c*x3**2 +I*(1+(a+b)/c)*x3

    Fx = np.concatenate((Fx1, Fx2, Fx3))

    # normalize
    x = x/np.max(x)
    Fx = Fx/np.max(Fx)
    
    if plot==True:
        plt.figure()
        plt.plot(x,fx,'-.')
    
        plt.figure()
        plt.plot(x,Fx,'-.')

    #save data
    tmp = np.vstack((x,fx,Fx)).T
    # np.savetxt("zprofile.txt",tmp)
    length = len(x)

    with open(folderName+"/zprofile.in", "w") as f:
        f.write(f"{length}\n")
        np.savetxt(f, tmp, fmt="%.6e", delimiter="\t")

def gen_gen4in(Lbuncht,folderName='.',slsc="ON",llsc="ON"):

    # c=2.998e8
    gam0=E0/0.511e6
    gambet=np.sqrt(gam0**2-1)
    
    # lamda0=99.93e-6
    # lamdau = 30e-3
    delz=lamdau/2
    
    Lbunch=c*Lbuncht
    
    
    setup = Namelist('setup', 
                     rootname = 'gen4',
                     lattice = 'gen4.lat',
                     beamline = 'THzBL',
                     gamma0 = gam0,
                     lambda0 = lamda0,
                     delz = delz,
                     seed = 2,
                     one4one = one4one,
                     npart = npart,
                     nbins = 4,
                     shotnoise = True)
    
    time = Namelist('time',
                    s0 = 0,
                    slen = Lbunch+100*lamda0,
                    sample = 1,
                    time = True)
    
    if one4one == True:    
        importbeam = Namelist('importbeam',
                         file = '../01_impz/gen4_sliced_one4one.h5',
                         time=True)
    else:
        importbeam = Namelist('importbeam',
                         file = f'../01_impz/scan.seed{seed}.nper{nperlambda}.fitdeg{degree}.out.par.h5',
                         time=True)
    
    field = Namelist('field',
                     power = 0,
                     waist_size=1e-3,
                     phase = 0,
                     dgrid = 20e-3,
                     ngrid = 501)
    
    if llsc=="OFF":
        longrange=0
    else:
        longrange=1
        
    if slsc=="OFF":
        nz=0
    else:
        nz=5
    efield = Namelist('efield',
                      longrange=longrange,
                      nz=nz,
                      rmax = 0.01,
                      nphi = 1,
                      ngrid = 100)
    
    track = Namelist('track',
                     output_step = 1,
                     field_dump_step = 0,
                     beam_dump_step = 0)
    
    g4 = Genesis4(setup, time, importbeam, field, efield, track)
    #print(g4.output)
    
    g4.output = g4.output.lower()
    g4.write(folderName+'/gen4.in')
    
    # generate the lattice
    with open(folderName+'/gen4.lat',"w",encoding="utf-8") as f:
        f.write(f"lcls:undulator={{aw={aw:.4f},lambdau={lamdau:.2e},nwig={nwig},helical={helical},kx={kx},ky={ky},ax=0,ay=0,gradx=0,grady=0}};\n")
        f.write("screen: marker={dumpfield = 0, dumpbeam = 0, stop = 1}; \n")
        f.write("THzBL: line={lcls, screen}; \n")
        
def gen_one(Q,lamdas,path="."):   
    
    with open(path+'/one',"w",encoding="utf-8") as f:
        f.write("#!/bin/bash \n\n")
        f.write("basedir=$(pwd) \n")
        
        f.write("source /etc/profile.d/modules.sh \n")
        f.write("module purge \n")
        f.write("module add mpi/mpich-x86_64 \n\n")
        
        f.write("cd $basedir \n")
        
        f.write("cd ./01_impz \n")
        f.write("genimpactzin lte.impz line \n")
        f.write("mpirun -np 4 ImpactZ.exe \n")
        
        if one4one == True:
            f.write(f"impz2sliceh5 gen4_dist.1001 {Q} {lamdas} 1 \n")
            #delete the .1001 file
            f.write("rm -f gen4_dist.1001 \n")
        else:
            f.write(f"xk_astra2gen4 fort.1001.001 {seed} {npart} {degree} {nperlambda} {lamda0} \n")
            f.write("rm -f fort.1001.001 \n")  
            
        f.write("\n\n")
        f.write("cd ../02_gen4 \n")
        f.write("module purge && module add mpi/openmpi-x86_64 hdf5/1.14.4 szip/2.1.1 \n")
        f.write("mpirun -np 4 genesis4 gen4.in \n")
    
    os.chmod(path+"/one", 0o755)

#------------------------------------------------------------------------
#------------------------------------------------------------------------
c=2.998e8

one4one = False
seed = 2
npart = 8192
degree = 3
nperlambda = 10

# profile='parabolic'
profile='flattop'

lamda0 = 100e-6
E0=17e6
gam0=E0/0.511e6
enx=5e-6;  eny=10e-6
sigE=10e3  #slice energy spread eV


lamdau=30e-3
nwig=113
helical=False
# aw=2.4678
aw= np.sqrt( 2*gam0**2*lamda0/lamdau -1)

if helical==False:
    unduu="Planar"
    Ku=np.sqrt(2)*aw
    K1 = 2*np.pi**2*Ku**2/(gam0**2*lamdau**2)    

    kx=0; ky=1
    #betax=8;  betay=0.4
    #alphax=4.55; alphay=2.5
    
    betax=8  
    alphax=4.55 
    betay=1/np.sqrt(K1)
    alphay=0

else:
    unduu="helical"
    Ku=aw
    K1 = 2*np.pi**2*Ku**2/(gam0**2*lamdau**2)    

    kx=0.5; ky=0.5
    # betax=0.4;  betay=0.4
    # alphax=0; alphay=0

    betax=1/np.sqrt(K1); betay=1/np.sqrt(K1)
    alphax=0;            alphay=0


tol_chargel = [1e-9, 2e-9, 3e-9]    
Lbunchtl = np.arange(4.0e-12, 21e-12, 1e-12)  # [s]

slsc="ON"
llsc="ON"
Np=1e6

# rootdir=f'/mnt/f/simu_2025/202510_idealMachineTHzFEL'
# rootdir=f'/mnt/f/simu_2025/202510_idealMachineTHzFEL/04_scan_scripts_impz_gen4'
# rootdir='/mnt/f/simu_2025/202510_idealMachineTHzFEL/scan_scripts'
# rootdir = '/mnt/f/simu_2025/202510_idealMachineTHzFEL/202511_scan_currentProfile_parabolic_flattop/00_ScanPythonScripts'
# rootdir = '/mnt/f/simu_2025/202510_idealMachineTHzFEL/202511_scan_currentProfile_parabolic_flattop/00_ScanPythonScripts'

# rootdir='/lustre/fs25/group/pitz/biaobin/202510_idealMachine/PlanarUndulator'
rootdir='/lustre/fs25/group/pitz/biaobin/202511_idealMachine'
#-------------------------------------------------------------------------

print(f"total job is {len(Lbunchtl)}*{len(tol_chargel)}=",len(Lbunchtl)*len(tol_chargel))
os.chdir(rootdir)

# import sys
# sys.exit()


#mkdir the folder
rootdir2=f'{unduu}_{profile}_slsc{slsc}_llsc{llsc}_one4one{one4one}_fixCharge_1125'
os.makedirs(rootdir2, exist_ok=True)
os.chdir(rootdir2)

for tol_charge in tol_chargel:
    for Lbuncht in Lbunchtl:

        #mkdir a new folder
        folderName = f"Q{tol_charge*1e9:.1f}nC_Lbuncht{Lbuncht*1e12:.2f}ps"
        os.makedirs(folderName, exist_ok=True)        
        os.chdir(folderName)
        
        #for impz
        os.makedirs("01_impz", exist_ok=True)
        gen_impzin(profile=profile,Lbuncht=Lbuncht,tol_charge=tol_charge, Np=Np, folderName='01_impz')
        
        #for gen4
        os.makedirs("02_gen4", exist_ok=True)
        # shutil.copy2(f"{rootdir}/runFile/gen4.lat", "02_gen4")        
        gen_gen4in(Lbuncht,folderName="02_gen4",slsc=slsc,llsc=llsc)
        
        # Now, for one & submit.sh
        gen_one(tol_charge,lamdas=lamda0,path=".",)        
        shutil.copy2(f"{rootdir}/runFile/submit.sh", ".")
                
        # #ready to submit all the jobs
        os.system("condor_submit submit.sh")
        os.chdir(rootdir+'/'+rootdir2)
        
os.chdir(rootdir)        
        
        
        



