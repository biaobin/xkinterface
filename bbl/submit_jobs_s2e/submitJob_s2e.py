import os
import numpy as np
import matplotlib.pyplot as plt
import h5py
from xkinterface.interface.PITZsim import PITZsim as pitz
import re
import sys
from impzpy.impactz_parser import impactz_parser as impz
from imptpy.impactt_parser import impactt_parser as impt
import time
from pitz_simu.submit_jobs_s2e import gentwobeam 

script_dir = os.path.dirname(os.path.abspath(__file__))
           
#=================================================================================================
def submitJob():
    #basedir = f"/afs/ifh.de/group/pitz/data/biaobin/sim2/2025/{beamprof}_bc{bc_stat}_{charge_q}nC"
    #basedir = f"/mnt/d/NextCloud/subjects/2025/01_twoBeam/{beamprof}_bc{bc_stat}_{charge_q}nC"
    if beamprof =="twobeam":
        basedir = f"/lustre/fs22/group/pitz/biaobin/2025/{beamprof}_bc{bc_stat}_{charge_q}nC_delay{delay:.1f}sigt_Q2_{charge_q2:.0f}pC_sigt2_{sigt2/1e-15:.0f}fs"                   
    else:
        basedir = f"/lustre/fs22/group/pitz/biaobin/2025/{beamprof}_bc{bc_stat}_{charge_q}nC"
    
    os.makedirs(basedir,exist_ok=True)
    os.chdir(basedir)
    
    charge= charge_q*1e-9
    one4one = 0    #if false, then resample will be called
    
    if beamprof=="gauss": 
        if len(seedl) > 1:
            # seed scan cases
            cores=[2,2]
            impt_run = 0
            impz_run = 0
            gen4_run = 1
        else:
            cores= [4,4]
            # whether to run impt, impz, or gen4 ?
            impt_run = 1
            impz_run = 1
            gen4_run = 1
    
        #mesh_impt= [64,64,64]
        #mesh_impz= [64,64,64]
        #Np=20*mesh_impt[0]*mesh_impt[1]*mesh_impt[2]

        mesh_impt= [32,32,64]
        mesh_impz= [32,32,64]
        Np=10*mesh_impt[0]*mesh_impt[1]*mesh_impt[2]
        
    elif beamprof == "flattop":
        if len(seedl) > 1:
            # seed scan cases
            cores=[2,2]
            impt_run = 0
            impz_run = 0
            gen4_run = 1
        else:
            cores= [4,4]
            # whether to run impt, impz, or gen4 ?
            impt_run = 1
            impz_run = 1
            gen4_run = 1
        #for study setting, to speed up the simulation 
        #mesh_impt= [32,32,64]
        #mesh_impz= [32,32,64]
        #Np=10*mesh_impt[0]*mesh_impt[1]*mesh_impt[2]
        mesh_impt= [64,64,64]
        mesh_impz= [64,64,64]
        Np=20*mesh_impt[0]*mesh_impt[1]*mesh_impt[2]
    
    elif beamprof == "twobeam":
        if len(seedl) > 1:
            # seed scan cases
            cores=[2,2]
            impt_run = 0
            impz_run = 0
            gen4_run = 1
        else:
            cores= [2,2]
            # whether to run impt, impz, or gen4 ?
            impt_run = 1
            impz_run = 1
            gen4_run = 1
        #for study setting, to speed up the simulation 
        mesh_impt= [32,32,64]
        mesh_impz= [32,32,64]
        if Np == None:
            Np=10*mesh_impt[0]*mesh_impt[1]*mesh_impt[2]                
    else:
        print("Error, beam prof not available!")
        sys.exit()
    
    mesh = mesh_impz
    
    dt= 2e-12
    ###Egun= 57.55
    ###phi1= 0
    Energy_gun = 5.81+0.511
    Energy= 17.05 #MeV
    
    # parameters for genesis simulation, resampling parameter
    npart = 4096*8
    
    if bc_stat == "OFF": 
        degree = 1
        nperlambda = 5
    else:
        degree = 3
        nperlambda = 10
    lambda0 = 100e-6  #for debug
        
    #=================================================================================================
    #=================================================================================================
    
    print(f"cores={cores}, Np={Np}")
    print("mesh_impt=",mesh_impt,"mesh_impz=",mesh_impz)
    print(f"profile={beamprof}, bc_stat={bc_stat}, q={charge_q}nC,phi1={phi1},phi2={phi2}, Isol={Isol}, seed={seed}")

    #% touch the folder
    # folderName = str.format('cores_%03d_impt_%02d%02d%04d_impz_%02d%02d%04d_Np_%.2e_phi2_%.1fdeg' 
    #                         % (cores[0]*cores[1],mesh_impt[0],mesh_impt[1],mesh_impt[2],
    #                            mesh_impz[0],mesh_impz[1],mesh_impz[2],Np,phi2) )
    #folderName = str.format('cores_%03d_impt_%02d%02d%04d_impz_%02d%02d%04d_Np_%.2e_phi2_%.1fdeg_phi1_%.0fdeg' 
    #                        % (16,mesh_impt[0],mesh_impt[1],mesh_impt[2],
    #                           mesh_impz[0],mesh_impz[1],mesh_impz[2],Np,phi2,phi1) )
    folderName = f"Np_{Np:.2e}_phi1_{phi1:.0f}deg_phi2_{phi2:.0f}deg"
    
    os.makedirs(folderName, exist_ok=True)

    #for solenoid scan
    if len(Isoll) > 1:
        imptFold = f"00_impt_{Isol}A"
    else:
        imptFold = "00_impt" 
    os.makedirs(folderName+f"/{imptFold}", exist_ok=True)
    os.makedirs(folderName+"/01_impz", exist_ok=True)
    os.makedirs(folderName+"/02_gen4", exist_ok=True)
    os.chdir(basedir)
        
    #% impt
    if impt_run == 1:
        os.chdir(basedir)
        #copy rfdata files for impt
        commd = f"cp {script_dir}/ini_simu/00_impt/rfdata* "
        commd += f"{script_dir}/ini_simu/00_impt/one "
        commd += f"{script_dir}/ini_simu/submit_condor.sh "
        commd += f"./{folderName}/{imptFold}"
        os.system(commd)
        
        wpath  = folderName+f"/{imptFold}/"
        os.chdir(wpath)
        
        fname= script_dir+'/ini_simu/00_impt/lte.impt'
        line = 'line'

        Egun = pitz.get_MaxE_gun(phi_gun=phi1, EG=Energy_gun) *1e6
        Eboost = pitz.get_MaxE_booster(MaxE_gun=Egun/1e6, phi_gun=phi1, phi_booster=phi2, Ef =Energy) *1e6
        # print("Eboost, phi2=",Eboost,phi2)
        if np.isnan(Eboost) or np.isnan(Egun):
            print(f"Error, Egun or Eboost is Nan. Stop. Egun={Egun:.2e} V/m,Eboo={Eboost:.2e} V/m")
            sys.exit()
        
        #update lte.impt with the given values
        lte = impt(fname,line)
        
        lte.control["CORE_NUM_T"] = cores[0]
        lte.control["CORE_NUM_L"] = cores[1]
        
        lte.control["MESHX"] = mesh_impt[0]
        lte.control["MESHY"] = mesh_impt[1]
        lte.control["MESHZ"] = mesh_impt[2]
        
        lte.control["DT"] = dt
       
        # pay attention to overwrite everything different from the template 
        if beamprof == "gauss":    
            BSA     = 3.5e-3
            sigx    = 0.83e-3
            sigy    = sigx

            #sigt    = 2.972399e-12
            #c_sig_t = 3.0
            # longer beam for 2nC
            sigt    = 4e-12
            c_sig_t = 4.0
    
            lte.beam["NP"] = Np
            lte.beam["TOTAL_CHARGE"] = charge
            lte.beam["DISTRIBUTION_TYPE"]="2"
    
            lte.beam["SIGX"]   = sigx
            lte.beam["SIGPX"] = 0.0
            lte.beam["C_SIG_X"] = BSA/2/sigx
    
            lte.beam["SIGY"]  = sigy
            lte.beam["SIGPY"] = 0.0
            lte.beam["C_SIG_Y"] = BSA/2/sigy
    
            lte.beam["SIGZ"] = 2.998e8*sigt
            lte.beam["SIGPZ"] = 0.0
            lte.beam["C_SIG_Z"] = c_sig_t
            lte.beam["DZ"] = -c_sig_t*sigt*2.998e8   #shift the beam behind the wall      
            lte.control["INI_T"] = -c_sig_t*sigt     #put the beam center to on-crest phase            
    
            lte.control["TEMISSION"] = sigt*2*c_sig_t                   
            lte.control["NEMISSION"] = int((sigt*2*c_sig_t)/0.02e-12)  #0.02ps/step     
            
        elif beamprof == "flattop":
            BSA=3.5e-3
            sigx=1.36e-3
            ft=16e-12
            rs=6e-12
            sigy=sigx
    
            lte.beam["NP"] = Np
            lte.beam["TOTAL_CHARGE"] = charge
            lte.beam["DISTRIBUTION_TYPE"]="221"
    
            lte.beam["SIGX"]   = sigx
            lte.beam["XSCALE"] = BSA/2/sigx
            lte.beam["SIGPX"] = 0.0
    
            lte.beam["SIGY"]  = sigy
            lte.beam["SIGPY"] = 0.0
            lte.beam["YSCALE"] = BSA/2/sigy
    
            lte.beam["SIGZ"] = ft*2.998e8  #flattop full width
            lte.beam["ZSCALE"] = rs/2*2.998e8 #zrise -> 2*zscale
            lte.beam["SIGPZ"] = 0.0
            lte.beam["DZ"] = 0.0           #overwrite the shift for gauss profile
            
            lte.control["INI_T"] = -1.0*(ft+2*rs)/2
            lte.control["TEMISSION"] = ft+2*rs                  
            lte.control["NEMISSION"] = int((ft+2*rs)/0.02e-12)  #0.02ps/step
        
        elif beamprof == "twobeam":    
            # generate initial two beam dist
            Q1= charge
            Q2= charge_q2*1e-12
            
            Np1=Np
            Np2=int(Q2/Q1*Np1)
            
            dt1 = 0
            dt2 = -sigt1*delay      #-sigt1*5
            
            Ekin=1 #eV
            
            param1={"Np":Np1,
                    "sigx":0.83e-3,
                    "sigy":0.83e-3,
                    "sigt":sigt1,
                    "c_sig_x":2.108434,
                    "c_sig_y":2.108434,
                    "c_sig_t":3.0,
                    "dt":dt1
                    }
            
            param2={"Np":Np2,
                    "sigx":0.83e-3,
                    "sigy":0.83e-3,
                    "sigt":sigt2,
                    "c_sig_x":2.108434,
                    "c_sig_y":2.108434,
                    "c_sig_t":3.0,
                    "dt":dt2
                    }
            
            inidis = gentwobeam.twobeam(param1, param2, Ekin)
            inidis.plotphase()
            inidis.dump_partcl_impt()
    
            # update settings for two beam scheme
            lte.control["TEMISSION"] = inidis.length_t                    
            lte.control["NEMISSION"] = int(inidis.length_t/0.02e-12)  #0.02ps/step
            lte.control["INI_T"] = -1.0*param1["c_sig_t"]*param1["sigt"]  # => head bunch center@design phase                   
    
            lte.beam["DISTRIBUTION_TYPE"] = 16
            lte.beam["NP"] = inidis.nptol
            lte.beam["TOTAL_CHARGE"] = Q1+Q2
            lte.beam["SIGZ"] = inidis.sigz_ct
            lte.beam["DZ"] = 0 #already shifted
    
        else:
            print("Error, wrong beam-profile, stop.")
            sys.exit(0)

        #scan the solenoid
        lte.lattice["GUN51SOL"]["SCALEB"]=Bsol
        
        phi_gun = 134.31+phi1
        phi_cds = 116.19+phi2

        lte.lattice["GUN51SOL"]["EMAX"] = Egun
        lte.lattice["GUN51SOL"]["PHASE"] = phi_gun
        
        lte.lattice["COUP1"]["PHASE"]= phi_cds
        lte.lattice["COUP1"]["EMAX"] = -Eboost
        
        lte.lattice["CELL"]["PHASE"] = phi_cds
        lte.lattice["CELL"]["EMAX"]  = Eboost
        
        lte.lattice["COUP2"]["PHASE"]= phi_cds
        lte.lattice["COUP2"]["EMAX"] = -Eboost
        try: 
            lte.lattice["W1"]["SAMPLE_FREQ"]=int(lte.beam["NP"]/10000)
            lte.lattice["W2"]["SAMPLE_FREQ"]=int(lte.beam["NP"]/10000)
        except:
            pass
        
        #generate ImpactT.in
        lte.back2lte()  #back2lte must be called before write
        lte.write_impacttin()
        print(f"ImpactT.in is updated.")
        os.chdir(basedir)
    
    #% impz
    #--------------------
    if impz_run == 1:
        os.chdir(basedir)
        commd = f"cp {script_dir}/ini_simu/01_impz/one "
        commd += f"{script_dir}/ini_simu/submit_condor.sh "
        commd += f"./{folderName}/01_impz"
        os.system(commd)
    
        wpath = folderName+'/01_impz/'
        os.chdir(wpath)
        if bc_stat == "ON": 
            fname = script_dir+"/ini_simu/01_impz/lte.impz"
        else:
            if beamprof == "flattop":
                fname = script_dir+"/ini_simu/01_impz/lte_bcoff_flattop.impz"
            else:
                fname = script_dir+"/ini_simu/01_impz/lte_bcoff.impz"
        line = 'line'
        lte = impz(fname,line)
        
        lte.control["CORE_NUM_T"] = cores[0]
        lte.control["CORE_NUM_L"] = cores[1]
        
        lte.control["MESHX"] = mesh_impz[0]
        lte.control["MESHY"] = mesh_impz[1]
        lte.control["MESHZ"] = mesh_impz[2]
        
        # may suffer particle loss, Np will be less than Np
        lte.beam["NP"] = Np
        lte.beam["TOTAL_CHARGE"] = charge
        
        #update lte.impz with the given values
        lte.back2lte()  #back2lte must be called before write
        lte.write_impactzin()
        print(f"ImpactZ.in is updated.")
        os.chdir(basedir)
    
    #% genesis-4
    if gen4_run == 1:
        os.chdir(basedir)
        #copy gen4.in and gen4.lte
        commd = str.format('cp %s/ini_simu/02_gen4/gen4.lte ./%s/02_gen4' % (script_dir,folderName))
        os.system(commd)
        
        wpath  = folderName+'/02_gen4/'
        os.chdir(wpath)
    
    #% submit the job
    os.chdir(basedir)
    # #copy submit.sh file
    # commd = str.format('cp %s/ini_simu/submit_condor.sh ./%s' % (script_dir,folderName))
    # os.system(commd)
    
    os.chdir(folderName)
    
    # touch one 
    content = "#!/bin/bash \n"
    content += "basedir=$(pwd) \n\n"
    
    #source the mpich module for impt&impz
    content += "source /etc/profile.d/modules.sh \n"
    content += "module purge \n"
    content += "module add mpi/mpich-x86_64 \n\n"
    
    #% now for impt
    if beamprof =="twobeam":
        #only for study
        #exe = "ImpactT_movie.exe"
    
        exe = "ImpactT.exe"
            
    else:
        exe = "ImpactT.exe"
    
    if impt_run == 1:
        content +="cd $basedir \n"
        content +=f"cd ./{imptFold} \n"
        content +="mpirun -np "+str(cores[0]*cores[1]) +f" {exe} \n"
        content += "mv fort.50 ../01_impz/particle.in \n\n"
    
    #%now for impz
    if impz_run == 1:
        content +="cd $basedir \n"
        content +="cd ./01_impz \n"
        ###content +="rm -f fort.* gen4_dist.* \n"
        content +="mpirun -np "+str(cores[0]*cores[1]) +" ImpactZ.exe \n"
        
        # move outputs to 02_gen4
        content += "mv gen4_dist.1031 ../02_gen4 \n"
        content += "mv fort.1032.001 ../02_gen4 \n\n"     
    
    #% now for genesis
    if gen4_run == 1:
        content +="cd $basedir \n"
        content +="cd ./02_gen4 \n"
        
        #get Nslice and beam_P0 from the beam distribution
        if one4one == 1:
            # for one4one
            content += "out=$(impz2sliceh5 gen4_dist.1031 " +f"{charge} \n" 
            content += "echo $out \n"
            content += 'read Nslice beam_P0 <<< "$out" \n'
            
        else:
            #not one4one, call xiangkun's module to transform astra dist => gen4 dist "
            content += "out=$(xk_astra2gen4 fort.1032.001 "+f"{seed} {npart} {degree} {nperlambda} {lambda0}) \n"
            content += "echo $out \n"
            content += 'read Nslice beam_P0 <<< "$out" \n'
        
        #update genesis input
        #2025-02-01, let's fix it to 17.05 MeV/c
        beam_P0 = 17.05
        # content += "xk_gen_gen4in "+f"{seed} {npart} $Nslice {lambda0} $beam_P0 {one4one} {nperlambda} {degree}\n\n"
        content += "xk_gen_gen4in "+f"{seed} {npart} $Nslice {lambda0} {beam_P0} {one4one} {nperlambda} {degree}\n\n"
    
        if one4one == 1:
            filein = "gen4.in"
        else:
            filein = "gen4.seed%d.nper%d.fitdeg%d.in" % (seed,nperlambda,degree)
        
        # Stop here, do not run genesis
        content +="module purge && module add python/3.9 mpi/openmpi-x86_64 phdf5/1.14.4 szip/2.1.1 \n"
        content +="mpirun -np " +str(cores[0]*cores[1]) +f" genesis4 {filein} \n"
    
        #delete the ini. beam dist. for seed scan loops, to save space
        if len(seedl)>1:    
            content +="rm -f " +str.format('scan.seed%d.nper%d.fitdeg%d.out.par.h5 \n' % (seed,nperlambda,degree))
    
    onefile = f"one_{seed}"
    with open(onefile, "w", encoding="utf-8") as file:
        file.write(content)
    os.system(f"chmod 755 {onefile}")
    
    # touch submit.sh
    if server == "condor":
        content  = "#Executable  = /usr/bin/bash\n"
        content += f"Executable = ./{onefile}\n"
        content += "request_cpus="+str(cores[0]*cores[1])+"\n"
        content += "\n"
        content += f"Log         = {onefile}.log\n"
        content += f"Output      = {onefile}.out\n"
        content += f"Error       = {onefile}.err\n"
        content += "\n"
        content += "getenv = True\n"
        content += "\n"
        content += "# as long as jobs run at DESY Zeuthen batch, please disable file transfer feature\n"
        content += "should_transfer_files = no\n"
        content += "# request 2GB RAM\n"
        content += "request_memory = 8GB\n"
        content += "\n"
        content += "Queue 1\n"
        
        file_path = f"submit_condor_{seed}.sh"
    
    else:
        node= int( np.ceil(cores[0]*cores[1]/32) )
        
        #pax server
        content = ""
        content += "#!/bin/bash\n\n"
        content += "#set job name:\n"
        content += "#SBATCH -J s2e\n\n"
        
        content += "#128 tasks on 4 nodes:\n"
        content += f"#SBATCH -N {node} \n\n"
        
        content += "#one process for each physical core:\n"
        content += "#SBATCH -c 2\n\n"
        
        content += "#run on the broadwell partition (pax11):\n"
        content += "#SBATCH -p broadwell\n\n"
        
        content += "#runtime of 20 minutes:\n"
        content += "#SBATCH -t 47:59:59\n\n"
        
        content += "#copy environment variables from submit environment:\n"
        content += "#SBATCH --get-user-env\n\n"
       
        content += "#send mail on all occasions:\n"
        content += "#SBATCH --mail-type=ALL\n\n"                
        content += f"./{onefile}\n"
        
        file_path = f"submit_pax_{seed}.sh"
    
    with open(file_path, "w", encoding="utf-8") as file:
        file.write(content)
    
    if server == "condor":
        commd = f"condor_submit submit_condor_{seed}.sh"
    else:
        commd = f"sbatch submit_pax_{seed}.sh"
        
    os.system(commd) 
    os.chdir(basedir)
    
    if len(seedl) > 1:
        time.sleep(5)

#===========================================================================
beamprofilel = ["gauss","flattop","twobeam"]
serverl = ["condor","pax"]

gauss  =0
flattop=1
twobeam=0

#%% gauss profile study
#---------------------------
#---------------------------
if gauss==1:
    server = serverl[0]
    beamprof = beamprofilel[0]

    #bc_statl = ["ON", "OFF"]
    #chargel =  [1, 1.5, 1.75, 2]
    #phi2l = np.arange(-16, -40 -1, -2)   #phi2 cannot be smaller than -40deg
    #Isoll = [368.1265306122449]

    bc_statl = ["ON"]
    chargel =  [2]
    phi2l = [-26]   #phi2 cannot be smaller than -40deg
    Isoll = [368.1265306122449]

    phi1l = [5]
    #phi1l = [-30,-35,-20,-15,-10,5,10,15,20]

    #1. run with seed=5
    seedl = [5]
    #2. run with seed scan, after 1 is finished
    #seedl = np.arange(6,25,1) 

#% flattop profile study
#---------------------------
#---------------------------
if flattop==1:
    server = serverl[0]
    beamprof = beamprofilel[1]
    
    #bc_statl = ["ON", "OFF"]
    #chargel =  [1, 1.5, 2]
    #phi2l = np.arange(0, -26-1, -2)
    ##seedl = [5]
    #seedl = np.arange(6,25,1)
    #Isoll = [367]  #opt found

    #test
    bc_statl = ["ON","OFF"]
    chargel =  [4]
    #phi2l = np.arange(-20,-4,+2)
    phi2l = np.arange(-32,-4,+2)
    seedl= [5]
    # scan the solenoid 
    Isoll = [367] 
    phi1l = [0]

    #phi1l = [-30,-20,-10,10,20,30]

#% settings for two beam study
#----------------------------
#---------------------------
if twobeam == 1:
    server = serverl[0]
    beamprof = beamprofilel[2]
    
    chargel =  [1]
    bc_statl =["ON"]
    # phi2l = np.arange(-20,-40-1,-2)
    phi2l = [-32]
    seedl = [5]
    Isoll = [368.1265306122449]

    phi1l = [0]
        
    #delay = 8  #2sig
    sigt1 = 2.9724e-12
    #sigt2 = 0.25e-12
    #charge_q2 = 60 #pC

    delay = 6  #2sig
    sigt2 = 2.9724e-12
    charge_q2 = 1000
    Np=5000

#%% loops
for Isol in Isoll:
    for bc_stat in bc_statl:
        for charge_q in chargel:
            for phi1 in phi1l:
                for phi2 in phi2l:
                    for seed in seedl:
                        Bsol = np.abs(pitz.I2B(Isol))
                        submitJob()
