import numpy as np
from scipy.special import erfinv
import matplotlib.pyplot as plt

clight = 2.998e8
mass = 0.511001e6

class twobeam:
    def __init__(self,param1,param2,Ekin=1):
        
        self.param1 = param1
        self.param2 = param2
        self.Ekin = Ekin #eV 
        self.nptol, self.phase = self.com2bunch()

    def dump_partcl_impt(self,path='./'):    
        # nptol, phase = self.com2bunch()
        
        # dump to impt 16 distribution
        if not path.endswith('/'):
            path += '/'
        with open(path+"partcl.data","w") as f:
            f.write(f"{self.nptol}\n")
    
            np.savetxt(f, self.phase, fmt="%16.8e")        

    def gen_gauss3d(self,Np=5000, sigx=0.83e-3,sigy=0.83e-3,sigt=2.9724e-12, \
                    c_sig_x=2.108434,c_sig_y=2.108434,c_sig_t=3.0,dt=0):
        x1 = np.random.uniform(0, 1, size=Np)
        x2 = np.random.uniform(0, 1, size=Np) 
        x3 = np.random.uniform(0, 1, size=Np) 
        
        xt = np.sqrt(2)*erfinv(2*x1-1)
        yt = np.sqrt(2)*erfinv(2*x2-1)
        tt = np.sqrt(2)*erfinv(2*x3-1)
        
        # cut
        mask = (xt/c_sig_x)**2 +(yt/c_sig_y)**2 < 1
        x = xt[mask]
        y = yt[mask]
        mask = tt < c_sig_t
        t = tt[mask]
        
        x = sigx*x
        y = sigy*y
        t = sigt*t+dt
        
        # make sure they have the same size
        if len(t) > len(x):
            indices = np.random.choice(len(t),len(x),replace=False)
            t = t[indices]
        else:
            indices = np.random.choice(len(x),len(t),replace=False)
            x = x[indices]
            y = y[indices]
            
        xyt = np.column_stack((x,y,t))    
        return xyt
    
    def com2bunch(self):
        mass = 0.511001e6
        gam = (self.Ekin+mass)/mass
        gambet = np.sqrt(gam**2-1)
        self.bet0 = gambet/gam
 
        pha1 = self.gen_gauss3d(**self.param1)
        pha2 = self.gen_gauss3d(**self.param2)

        x = np.concatenate((pha1[:,0],pha2[:,0]))
        y = np.concatenate((pha1[:,1],pha2[:,1]))
        t = np.concatenate((pha1[:,2],pha2[:,2]))

        #shift to behind the wall
        t -=np.max(t)
        self.length_t = np.max(t)-np.min(t)

        nptol =len(x)
        px=py = np.zeros(nptol)
        pz = np.ones(nptol)*gambet
        z = self.bet0*clight*t
        self.sigz_ct = np.std(clight*t)
        
        #x(m), gambetx, y(m), gambety, z(m), gambetz
        phase = np.column_stack((x,px,y,py,z,pz))
        
        return nptol, phase
    
    def plotphase(self):
        phase = self.phase
        
        plt.figure(figsize=(23,6))
        
        plt.subplot(1,3,1)
        plt.hist2d(phase[:,0]*1e3,phase[:,2]*1e3,bins=100,cmap="jet")
        # plt.axis("equal")
        plt.xlabel("x (mm)")
        plt.ylabel("y (mm)")
        
        plt.subplot(1,3,2)
        plt.hist2d(phase[:,4]/self.bet0/clight*1e12,phase[:,0]*1e3,bins=100,cmap="jet")
        # plt.axis("equal")
        plt.xlabel("z (ps)")
        plt.ylabel("x (mm)")
        
        # plt.figure()
        plt.subplot(1,3,3)
        plt.hist(phase[:,4]/self.bet0/clight*1e12,100)
        plt.xlabel("z (ps)")
        plt.savefig('inidis.png', bbox_inches='tight')
        
        plt.close()

if __name__ =="__main__":

    Q1= 1e-9
    Q2= 100e-12
    
    Np1=int(1e5)
    Np2=int(Q2/Q1*Np1)
    
    sigt1 = 2.9724e-12
    sigt2 = sigt1/5
    dt1 = 0
    dt2 = -sigt1*1.5
    
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
    
    tmp = twobeam(param1, param2, Ekin)
    
    tmp.plotphase()
    
    tmp.dump_partcl_impt()





















