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
script_dir = os.path.dirname(os.path.abspath(__file__))

phi2l = np.arange(-16, -40 -1, -2)
# phi2l = np.arange(-18, -38 -1, -4)
# phi2l= [-28]

seedl = np.arange(6,26,1)
# chargel=[1, 1.5, 2] #nC
chargel=[1.75]
bc_statl = ["ON","OFF"]

for bc_stat in bc_statl:
    for charge_q in chargel:
        for phi2 in phi2l:
            for seed in seedl:
                #%% inputs
                #=================================================================================================
                #=================================================================================================
                
                # basedir=r'/afs/ifh.de/group/pitz/data/biaobin/sim1/2025/00_sase_chicane/Q_1_1.5_2nC/bcON_2nC'
                #basedir=r'/afs/ifh.de/group/pitz/data/biaobin/sim1/2025/00_sase_chicane/Q_1_1.5_2nC/bcON_1.5nC'
                #basedir=r'/afs/ifh.de/group/pitz/data/biaobin/sim1/2025/00_sase_chicane/Q_1_1.5_2nC/bcON_2nC'
                
                # charge_q = 1     #nC
                # bc_stat  = "ON"
                
                basedir = f'/afs/ifh.de/group/pitz/data/biaobin/sim2/2025/bc{bc_stat}_{charge_q}nC'
                os.chdir(basedir)
                
                # whether to run impt, impz, or gen4 ?
                impt_run = 0
                impz_run = 0
                gen4_run = 1
                
                charge= charge_q*1e-9
                one4one = 0    #if false, then resample will be called
                
                cores= [1,1]
                # cores=[1,1]  #for runing python only
                
                mesh_impt= [64,64,64]
                mesh_impz= [64,64,64]
                
                mesh = mesh_impz
                Np=20*mesh[0]*mesh[1]*mesh[2]
                
                dt= 2e-12
                Egun= 57.55
                phi1= 0
                Energy= 17.05 #MeV
                
                # parameters for genesis simulation, resampling parameter
                #seed = 5
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
           
                print("cores=",cores)
                print("mesh_impt=",mesh_impt,"mesh_impz=",mesh_impz)
                print("phi2=",phi2)
                print("Np=",Np)
                #% touch the folder
                # folderName = str.format('cores_%03d_impt_%02d%02d%04d_impz_%02d%02d%04d_Np_%.2e_phi2_%.1fdeg' 
                #                         % (cores[0]*cores[1],mesh_impt[0],mesh_impt[1],mesh_impt[2],
                #                            mesh_impz[0],mesh_impz[1],mesh_impz[2],Np,phi2) )
                folderName = str.format('cores_%03d_impt_%02d%02d%04d_impz_%02d%02d%04d_Np_%.2e_phi2_%.1fdeg' 
                                        % (16,mesh_impt[0],mesh_impt[1],mesh_impt[2],
                                           mesh_impz[0],mesh_impz[1],mesh_impz[2],Np,phi2) )
                
                os.makedirs(folderName, exist_ok=True)
                os.makedirs(folderName+"/00_impt", exist_ok=True)
                os.makedirs(folderName+"/01_impz", exist_ok=True)
                os.makedirs(folderName+"/02_gen4", exist_ok=True)
                os.chdir(basedir)
                    
                #%% impt
                if impt_run == 1:
                    os.chdir(basedir)
                    #copy rfdata files for impt
                    commd = str.format('cp %s/ini_simu/00_impt/rfdata* ./%s/00_impt' % (script_dir,folderName))
                    os.system(commd)
                    
                    wpath  = folderName+'/00_impt/'
                    os.chdir(wpath)
                    
                    fname= script_dir+'/ini_simu/00_impt/lte.impt'
                    line = 'line'
                    
                    Eboost = pitz.get_MaxE_booster(MaxE_gun=Egun, phi_gun=phi1, phi_booster=phi2, Ef =Energy) *1e6
                    # print("Eboost, phi2=",Eboost,phi2)
                    
                    #update lte.impt with the given values
                    lte = impt(fname,line)
                    
                    lte.control["CORE_NUM_T"] = cores[0]
                    lte.control["CORE_NUM_L"] = cores[1]
                    
                    lte.control["MESHX"] = mesh_impt[0]
                    lte.control["MESHY"] = mesh_impt[1]
                    lte.control["MESHZ"] = mesh_impt[2]
                    
                    lte.control["DT"] = dt
                    lte.beam["NP"] = Np
                    lte.beam["TOTAL_CHARGE"] = charge
                    
                    # phi_gun = 134.31
                    phi_cds = 116.19+phi2
                    
                    lte.lattice["COUP1"]["PHASE"]= phi_cds
                    lte.lattice["COUP1"]["EMAX"] = -Eboost
                    
                    lte.lattice["CELL"]["PHASE"] = phi_cds
                    lte.lattice["CELL"]["EMAX"]  = Eboost
                    
                    lte.lattice["COUP2"]["PHASE"]= phi_cds
                    lte.lattice["COUP2"]["EMAX"] = -Eboost
                    
                    #generate ImpactT.in
                    lte.write_impacttin()
                    print(f"ImpactT.in is updated.")
                    os.chdir(basedir)
                
                #%% impz
                if impz_run == 1:
                    os.chdir(basedir)
                    wpath = folderName+'/01_impz/'
                    os.chdir(wpath)
                    if bc_stat == "ON": 
                        fname = script_dir+"/ini_simu/01_impz/lte.impz"
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
                    lte.write_impactzin()
                    print(f"ImpactZ.in is updated.")
                    os.chdir(basedir)
                
                #%% genesis-4
                if gen4_run == 1:
                    os.chdir(basedir)
                    #copy gen4.in and gen4.lte
                    commd = str.format('cp %s/ini_simu/02_gen4/gen4.lte ./%s/02_gen4' % (script_dir,folderName))
                    os.system(commd)
                    
                    wpath  = folderName+'/02_gen4/'
                    os.chdir(wpath)
                
                #%% submit the job
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
                
                #%%% now for impt
                if impt_run == 1:
                    content +="cd $basedir \n"
                    content +="cd ./00_impt \n"
                    content +="mpirun -np "+str(cores[0]*cores[1]) +" ImpactT.exe \n"
                    content += "mv fort.50 ../01_impz/particle.in \n\n"
                
                #%%% now for impz
                if impz_run == 1:
                    content +="cd $basedir \n"
                    content +="cd ./01_impz \n"
                    ###content +="rm -f fort.* gen4_dist.* \n"
                    content +="mpirun -np "+str(cores[0]*cores[1]) +" ImpactZ.exe \n"
                    
                    # move outputs to 02_gen4
                    content += "mv gen4_dist.1031 ../02_gen4 \n"
                    content += "mv fort.1032.001 ../02_gen4 \n\n"     
                
                #%%% now for genesis
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
                    #2025-02-04, let's fix it to 17.05 MeV/c
                    beam_P0 = 17.05
                    content += "xk_gen_gen4in "+f"{seed} {npart} $Nslice {lambda0} {beam_P0} {one4one} {nperlambda} {degree}\n\n"
                    
                    if one4one == 1:
                        filein = "gen4.in"
                    else:
                        filein = "gen4.seed%d.nper%d.fitdeg%d.in" % (seed,nperlambda,degree)
                        
                    content +="module purge && module add python/3.9 mpi/openmpi-x86_64 phdf5/1.14.4 szip/2.1.1 \n"
                    content +="mpirun -np " +str(cores[0]*cores[1]) +f" genesis4 {filein} \n"
    
                    #delete the ini. beam dist. to save space
                    content +="rm -f " +str.format('scan.seed%d.nper%d.fitdeg%d.out.par.h5 \n' % (seed,nperlambda,degree))
    
                onefile = f"one_{seed}"
                with open(onefile, "w", encoding="utf-8") as file:
                    file.write(content)
                os.system(f"chmod 755 {onefile}")
                
                # touch submit.sh
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
                content += "request_memory = 16GB\n"
                content += "\n"
                content += "Queue 1\n"
                
                filename = f"submit_condor_{seed}.sh"
                with open(filename, "w", encoding="utf-8") as file:
                    file.write(content)
                
                commd = f"condor_submit submit_condor_{seed}.sh"
                os.system(commd) 
                os.chdir(basedir)
                
                time.sleep(60)
