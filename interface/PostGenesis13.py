# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 14:57:04 2020

@author: lixiangk
"""

from .Misc import *
from scipy.interpolate import interp1d

#%%
try:
    import h5py
except Exception as err:
    print(err)
    print('Please install h5py!')
    
    
from numpy.fft import fftshift,fft

def calcSpectrum(amp, phase = None, lambda0 = 100e-6, sample = 1, freq0 = None):
    '''
    Calculate the spectrum from samples

    Parameters
    ----------
    amp : 1D or 2D array
        Amplitudes of samples of the signal. In the case of 2D, the first dimension 
        is along the slices and the second dimension is along the undulator. 
        Can also be the complex fields, then `phase` is not needed.
    phase : 1D or 2D array
        Phases of samples of the signal. The default is None, if the amp is given as complex fields.
    lambda0 : double, optional
        Seperation of sampling (usually the wavelength) in meter. The default is 100e-6.
    freq0 : double, optional
        Sampling frequency. If defined, it dominates `lambda0`. The default is None.
        
    Returns
    -------
    wavelength : 1D array
        Wavelength of the signal transformed in frequency domain.
    spectrum : 1D or 2D array
        Spectra intensity of the signal transformed in frequency domain.

    '''
    
    if phase is not None:
        signal = np.sqrt(amp)*(np.cos(phase)+np.sin(phase)*1j)
    else:
        if np.iscomplexobj(amp):
            signal = amp
        else:
            print('The amp should be the complex fields if no phase is given!')
    
    nsample = len(signal) # number of samples
    
    axis = 0
    spectrum = np.abs(fftshift(fft(signal, nsample, axis), axis))
    spectrum = spectrum*spectrum
    
    if freq0 is None and lambda0 is not None:
        freq0 = 1./(lambda0/sample)*g_c # sampling frequency
    
    F = 1.0*np.arange(-nsample/2, nsample/2,)/nsample*freq0+freq0 # Frequency
    wavelength = g_c/F # Frequency to wavelength
    
    return wavelength, spectrum
    
class PostGenesis13:
    version = 4 # default
    def __init__(self, fname = None, version = 3, **kwargs):
        '''
        
        Parameters
        ----------
        fname : str
            Genesis V4 main output file.
        **kwargs : TYPE
            DESCRIPTION.
            
        Returns
        -------
        None.

        '''
        
        debug = 0
        self.fieldname = 'Field'
        if len(kwargs)>0:
            if 'debug' in kwargs.keys():
                debug = kwargs['debug']    
            if 'harmonic' in kwargs.keys():
                harmonic = kwargs['harmonic']
                self.fieldname = 'Field%d' % harmonic
                
        if fname is None:
            fname = get_file('.h5'); print(fname)
        
        _, _, ext = fileparts(fname); print(ext)
        
        if ext.upper() in ['.H5']:
            if version == 3:
                self.version = 3
                self.load3(fname)
            else:
                self.version = 4
                self.load4(fname)
                
        elif ext.upper() in ['.OUT']:
            self.version = 2
            self.load2(fname)    
        else:
            print('Unknown file extension!')
            return
        
        self.Lsat = self.get_satpos()
        if debug:  
            if 'fig_ext' in kwargs.keys():
                fig_ext = kwargs['fig_ext']
            else:
                fig_ext = '.png'

            plt.figure(figsize=(18,10))
            plt.subplot(2,3,1)
            self.plot_current(fig_ext = fig_ext, fig=False)
            plt.subplot(2,3,2)
            self.plot_power(fig_ext = fig_ext, fig=False)
            plt.subplot(2,3,3)
            self.plot_energy(fig_ext = fig_ext, fig=False)
            plt.subplot(2,3,4)
            self.plot_tpower(fig_ext = fig_ext, fig=False, at=self.Lsat)
            plt.subplot(2,3,5)
            self.plot_spectrum(fig_ext = fig_ext, fig=False, at=self.Lsat)
       
    def print_all(self, prop = None, order = 2):
        if self.version == 2:
            print(self.outputs)
            return
        
        print('./')
        for key in self.file.keys():
            if prop is None:
                print('  - %s' % key)
                if order == 1:
                    continue
                for subkey in self.file.get(key).keys():
                    print('    - %s' % subkey)
            elif prop.upper() == key.upper():
                print('  - %s' % key)
                for subkey in self.file.get(key).keys():
                    print('    - %s' % subkey)
                break
            
    def load4(self, fname):
        file = h5py.File(fname, 'r')
        
        tmp = {}
        for k, v in file.get('Beam').items():
            tmp.update({k:v[:]})
        self.beam = tmp
        
        tmp = {}
        for k, v in file.get(self.fieldname).items():
            tmp.update({k:v[:]})
        self.field = tmp
        
        nstep, nslice = file.get(self.fieldname).get('power').shape
        
        meta = file.get('Meta').get('InputFile')[0].decode().split()
        
        kv = {}
        for _, a in enumerate(meta):
            tmp = a.split('=')
            if len(tmp)>1:
                kv.update({tmp[0]:tmp[1]})
                
        self.file = file
        self.nslice = nslice
        self.nstep = nstep
        self.kv = kv
        
        self.outputs = [temp for temp in self.file.get(self.fieldname).keys()] # Output properties for radiation fields
        
        self.sample = file.get('Global/sample')[0]
        self.lambdaref = file.get('Global/lambdaref')[0]
        self.gamma0 = file.get('Global/gamma0')[0]
        self.one4one = file.get('Global/one4one')[0]
        self.scan = file.get('Global/scan')[0]
        self.slen = file.get('Global/slen')[0]
        self.time = file.get('Global/time')[0]
        
        lslice = self.lambdaref*self.sample
        
        self.lslice = lslice # length of a slice
        #self.zbunch = np.linspace(lslice/2, self.slen-lslice/2, nslice)
        self.zbunch = np.arange(nslice)*self.lslice+self.lslice/2
        self.current = file.get('Beam/current')[:].flatten()
        
        #self.bunching = file.get('Beam/bunching')
        
        self.zplot = file.get('Lattice/zplot')[:] # z coordinates along the undulator
        self.zaw = file.get('Lattice/aw')[:]
        
        self.zstep = self.zplot[1]-self.zplot[0]  
        
        self.zpower = np.sum(file.get(self.fieldname).get('power')[:], axis = 1) # power along the undulator
        self.zenergy = self.zpower*self.lslice/g_c # energy along the undulator
        
        self.power = np.sum(self.zpower)    # Total power
        self.energy = np.sum(self.zenergy)  # Total energy
        
        # spectrum at the end of the undulator
        self.wavelength, self.spectrum = self.get_spectrum(self.zplot[-1])

        # get cum_dE_lsc
        self.get_cum_dE_lsc(delz=0.015)    
            
        return
    
    load = load4
    
    def float_no_error(self, ss):
        try:
            r = float(ss)
        except Exception as err:
            print(err)
            r = 0
        return r
    
    def load2(self, fname):
        
        '''
        Read the standard output (in ascii format and stored in `fname`) 
        into memories for further analysis
        '''
        
        file = open(fname, 'r')
        line = file.readline()

        # First, read the input namelist from the header of the file 
        kv = {}
        while True:
            if not line:
                break

            match = re.search(r'^[^=]*=[^=]*$' , line)
            if match:
                # print line
                key, value = line.split('=')
                match = re.search(r'file', key)
                if match:
                    kv.update({key.strip():value})
                    line = file.readline()
                    continue

                if len(value.split()) == 1:
                    value = float(value.replace("D", "E"))
                else:
                    value = [float(v.replace("D", "E")) for v in value.split()]
                kv.update({key.strip():value})          

            match = re.search(r'^[ ]*z\[m\][ ]*aw[ ]*qfld[ ]*$', line)
            if match:
                # print line
                nstep = int(kv['zstop']/kv['delz']/kv['xlamd'])+1
                #nstep = 17
                field = np.zeros((nstep, 3))
                for i in np.arange(nstep):
                    line = file.readline()
                    field[i] = np.array([float(v.replace("D", "E")) for v in line.split()])
                break
            line = file.readline()
        
        nslice = int(kv['nslice'])
        if nslice == 0:
            nslice = 3110
        nc = int(np.sum(kv['lout'])); #print nc # number of output items

        # Initialize the arrays for storing the data blocks
        current = np.zeros((nslice, 1))
        data = np.zeros((nslice, nstep, nc)); #print data.shape
        
        # Then read the blocks of data from the file
        islice = 0
        while True:
            if not line:
                break
            match = re.search(r'[ E]* current', line)
            if match:
                icurrent, tmp = line.split()
                icurrent = float(icurrent)
                current[islice, 0] = icurrent
        
                line = file.readline()
                line = file.readline()
                line = file.readline()

                outputs = line.split()

                for j in np.arange(nstep):
                    line = file.readline(); #print line
                    line = line.replace('NaN', '  0')
                    data[islice, j] = np.array([self.float_no_error(v.replace("D", "E")) for v in line.split()])
                #break
                islice += 1
            line = file.readline()
            
        file.close()
        
        self.kv = kv # dict with key-value pairs
        self.nslice = nslice; self.nstep = nstep; self.nc = nc
        self.field = field
        self.lattice = field
        
        self.outputs = outputs # Output properties for radiation fields
        
        #self.data = data
        # to make the first two dimensions consistent with version 4
        # now nstep, nslice, nc
        self.data = np.transpose(data, [1, 0, 2]) 
        
        self.gamma0 = self.kv['gamma0']
        self.time = self.kv['itdp']
        self.lambdaref = self.kv['xlamds']
        self.sample = 1.0/self.kv['zsep']
        
        self.lslice = self.lambdaref/self.sample # length of a slice
        
        self.zbunch = np.arange(nslice)*self.lslice+self.lslice/2
        self.current = current
        
        self.zplot = field[:,0]; self.zstep = self.zplot[1]-self.zplot[0] # z coordinates along the undulator
        self.zaw = field[:,1]
        
        self.zpower = self.data[:,:,0].sum(axis = 1) # power along the undulator
        self.zenergy = self.zpower*self.lslice/g_c   # energy along the undulator
        
        self.power = np.sum(self.zpower)    # Total power
        self.energy = np.sum(self.zenergy)  # Total energy
        
        self.wavelength, self.spectrum = self.get_spectrum(at = -999) # spectrum at the end of the undulator
        return
    
    def load3(self, fname):
        
        '''
        Read the standard output (in ascii format and stored in `fname`) 
        into memories for further analysis
        '''
        
        file = h5py.File(fname, 'r')
    
        # First, read the input namelist from the header of the file 
        kv = {}
        
        for key in file.get('input').keys():
            value = file.get('input').get(key)[0]
            kv.update({key:value})
        
        self.kv = kv # dict with key-value pairs
        
        aw = file.get('lattice').get('aw')[:]
        zz = file.get('lattice').get('z')[:]
        qfld = file.get('lattice').get('qfld')[:]
        field = np.zeros((len(aw), 3))
        field[:,0] = zz
        field[:,1] = aw
        field[:,2] = qfld
        
        self.field = field
        nslice = self.kv.get('nslice')
        self.nslice = nslice
        
        shape = file.get('power')[:].shape
        data = np.zeros((shape[0], shape[1], 10))
        
        keys = ['power', 'increment', 'signalamp', 'signalphase', 'radsize',
                'energy', 'bunching', 'xbeamsize', 'ybeamsize', 'error']
        for i, key in enumerate(keys):
            data[:,:,i] = file.get(key)[:]
        
        #self.outputs = outputs # Output properties for radiation fields
        
        self.outputs = ['power', 'increment', 'p_mid', 'phi_mid', 'r_size', 
                        'energy', 'bunching', 'xrms', 'yrms', 'error']
        
        self.data = data
        # to make the first two dimensions consistent with version 4
        # now nstep, nslice, nc
        #self.data = np.transpose(data, [1, 0, 2]) 
        self.nstep = self.data.shape[0]
        
        self.gamma0 = self.kv['gamma0']
        self.time = self.kv['itdp']
        self.lambdaref = self.kv['xlamds']
        self.sample = 1.0/self.kv['zsep']
        
        self.lslice = self.lambdaref/self.sample # length of a slice
        
        self.zbunch = np.arange(nslice)*self.lslice+self.lslice/2
        current = file.get('current')[:]
        self.current = current
        
        self.zplot = field[:,0]; self.zstep = self.zplot[1]-self.zplot[0] # z coordinates along the undulator
        
        self.zpower = self.data[:,:,0].sum(axis = 1) # power along the undulator
        self.zenergy = self.zpower*self.lslice/g_c   # energy along the undulator
        
        self.power = np.sum(self.zpower)    # Total power
        self.energy = np.sum(self.zenergy)  # Total energy

        #self.get_cum_dE_lsc(delz=0.015)    
        
        self.wavelength, self.spectrum = self.get_spectrum(at = -999) # spectrum at the end of the undulator
        return

    def get_fielddata(self, name, at = None):
        '''
        Get the field data/property defined by `name` along the slices at position z = at
        Parameters
          name: string
              Name of a property such as 'power', 'increment', 'p_mid', and so on for version 2
              and 'power', 'intensity-farfield', 'phase-farfield' and so on for version 4
          at: float or None
              Longitudinal position along the undulator, used to calculate the nth step of 
              interest; if None, return all
        Returns
        -------
        data : 1D array
            The requested data along the slices
        '''
        if at is not None:
            if at < 0:
                col = -1
            else:
                #col = int(at/self.zstep)
                col = np.where(np.abs(self.zplot-at)<self.zstep/2)[0][0]
        else:
            col = slice(0, self.nstep)
        #print('col = ', col)
        
        props = [temp.upper() for temp in self.outputs]
        if name.upper() not in props: # check if `name` is included in the outputs
            print('Available properties are: ', self.outputs)
            return None
            
        if self.version<4:
            
            third = self.outputs.index(name)
            return self.data[col,:,third]
        else:
            return self.get_data(self.fieldname, name)[:][col]

    def get_beamdata(self, name, at = None):
        '''
        Get the beam data/property defined by `name` along the slices at position z = at
        Parameters
          name: string
              Name of a property such as 'power', 'increment', 'p_mid', and so on for version 2
              and 'power', 'intensity-farfield', 'phase-farfield' and so on for version 4
          at: float or None
              Longitudinal position along the undulator, used to calculate the nth step of 
              interest; if None, return all
        Returns
        -------
        data : 1D array
            The requested data along the slices
        '''
        if at is not None:
            if at < 0:
                col = -at
            else:
                #col = int(at/self.zstep)
                col = np.where(np.abs(self.zplot-at)<self.zstep/2)[0][0]
        else:
            col = slice(0, self.nstep)
        #print(col)
        
        if self.version < 4:
            props = [temp.upper() for temp in self.outputs]
        else:
            props = [temp.upper() for temp in self.file.get('Beam').keys()]
            
        if name.upper() not in props: # check if `name` is included in the outputs
            print('Available properties are: ', self.outputs)
            return None
            
        if self.version<4:
            third = self.outputs.index(name)
            return self.data[col,:,third]
        else:
            return self.get_data('Beam', name)[:][col]
    
    def get_spectrum(self, at = None, nearfield = False):
        '''
        Get the radiation spectrum at position z = at
        Parameters
          at: float or None
              Longitudinal position along the undulator, used to calculate the
              nth step of interest; if None, return all
        Returns
        -------
        wavelength : 1D array
            Wavelength of the signal transformed in frequency domain.
        spectrum : 1D or 2D array
            Intensity of the signal transformed in frequency domain. if at == None, 
            then return a 2D array, with the first dimension along the slices and
            the second dimension along the undulator.
        '''
        
        if self.version<4:
            n1, n2 = 'p_mid', 'phi_mid'
        else:
            if nearfield:
                n1, n2 = 'intensity-nearfield', 'phase-nearfield'
            else:
                n1, n2 = 'intensity-farfield', 'phase-farfield'
                
        amp = self.get_fielddata(n1, at)
        phase = self.get_fielddata(n2, at)
        print(n1, n2, at)
        
        return calcSpectrum(amp, phase, self.lslice, self.sample)
    
    def get_data(self, path, *args):
        '''
        Get data by the `path` from the hdf5 file, only valid for version 4

        Parameters
        ----------
        path : string
            Path of the variable of interest, e.g., 'Field/power'.
        *args: string
            subpath, e.g., get_data('Field', 'power') is equiverlent to get_data('Field/power')
            
        Returns
        -------
        data : array
            The dimension of the array depends on the variable type.

        '''
        if self.version == 2 or self.version == 3:
            print('Only work for version 4!')
            return
        
        if self.version == 4:
            if len(args)>0:
                for temp in args:
                    path = path+'/'+temp
            try:
                data = self.file[path][:]
            except Exception as err:
                print(err)
                data = None
        return data
        
    def plot_current(self, x = 'time', fig_ext = '.png', fig=True):
        if fig==True: 
            plt.figure(figsize=(8,6))
        
        if x.upper() == 'TIME':
            plt.plot(self.zbunch/g_c*1e12, self.current, '-')
            plt.xlabel(r'Time (ps)')
        elif x.upper() == 'LENGTH':
            plt.plot(self.zbunch, self.current, '-')
            plt.xlabel(r'Position (mm)')
        elif x.upper() == 'SLICE':
            plt.plot(self.current, '-')
            plt.xlabel(r'slice #')
            
        plt.ylabel(r'Current (A)')
        plt.grid()
        if fig_ext != None:
            plt.savefig('current'+fig_ext)
        
    def plot_spectrum(self, fig_ext = '.png', at=-1, fig=True, xrange=[80,120]):
        self.wavelength, self.spectrum = self.get_spectrum(at=at)
        if fig==True: 
            plt.figure(figsize=(8,6))
 
        plt.plot(self.wavelength*1e6, self.spectrum, '-',label=f"@{at:.2f} m")
        plt.xlabel(r'Wavelength ($\mu$m)')
        plt.ylabel(r'Intensity (arb. unit)')
        plt.legend()
        plt.xlim(xrange)
        plt.grid()
        #plt.title(f"@{at:.2f} m")
        if fig_ext != None:
            plt.savefig('spectrum'+fig_ext)

    def get_tpower(self, at=-1):
        self.tpower = self.get_fielddata("power",at)
        self.t = self.zbunch/2.998e8
        
        return self.tpower

    def plot_tpower(self,x="time", fig_ext='.png', at=-1, fig=True, xrange=None, linelabel=None, linestyle=None):
        if linelabel == None:
            linelabel = f"s={at:.2f} m"        
        #else:
        #    linelabel += f", s={at:.2f} m"  

        if linestyle==None:
            linestyle="-"

        tpower = self.get_fielddata("power",at)
        t = self.zbunch/2.998e8
        if fig==True: 
            plt.figure(figsize=(8,6))
 
        if x=="SLICE":
            plt.plot(tpower,linestyle,label=linelabel)
            plt.xlabel(r'slice #')
        else:
            plt.plot(t*1e12, tpower/1e6,linestyle,label=linelabel)
            plt.xlabel(r't (ps)')
        plt.ylabel(r'power (MW)')
        plt.legend()
        if xrange != None:
            plt.xlim(xrange)
        plt.grid()
        #plt.title(f"@{at:.2f} m")
        if fig_ext != None:
            plt.savefig('t_power'+fig_ext)
       
    def plot_power(self, fig_ext = '.png', fig=True):
        if fig==True: 
            plt.figure(figsize=(8,6))
          
        plt.plot(self.zplot, self.zpower/1e6, '-')
        plt.xlabel(r'$z$ (m)')
        plt.ylabel(r'Power (MW)')
        plt.yscale('log')
        plt.grid()
        if fig_ext != None:
            plt.savefig('power-z'+fig_ext)
        
    def plot_energy(self, fig_ext = '.png', fig=True):
        if fig==True: 
            plt.figure(figsize=(8,6))
       
        plt.plot(self.zplot, self.zenergy*1e6, '-')
        plt.xlabel(r'$z$ (m)')
        plt.ylabel(r'Energy ($\mu$J)')
        plt.yscale('log')
        plt.grid()
        if fig_ext != None:
            plt.savefig('energy-z'+fig_ext)

    def get_satpos(self):
        interp_func_Lsat= interp1d(self.zenergy, self.zplot, kind='linear', fill_value="extrapolate")
        Lsat = interp_func_Lsat(0.8*np.max(self.zenergy))

        return Lsat

    def get_cum_dE_lsc(self, delz=0.015):
        #dz, step in gen4.in, i.e. delz
        lsc = self.file["Beam"]["LSCfield"]
        current = self.file["Beam"]["current"][:].flatten()
        
        nstep = lsc.shape[0]
        
        dE_lsc =[]
        for j in np.arange(nstep):
        
            lsc_s = lsc[j].flatten()
            #dE_lsc.append( lsc_s*1.6e-19 *current*self.lslice/2.998e8/1.6e-19 *delz *1e6) #[uJ/slice]
            dE_lsc.append( lsc_s*delz *1.6e-19*1e6) #[uJ/electron]
        
        dE_lsc = np.array(dE_lsc)
        self.cum_dE_lsc = np.cumsum(dE_lsc, axis=0)

PostGenesis = PostGenesis13
