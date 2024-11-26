# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 16:56:20 2020

@author: lixiangk
"""

from interface import *

#%% Define the analyzing function
def Analyze(start, end, z0 = 3.6, ext = '.txt', rootname = 'pithz', version = 4):
    
    remove_files(['power-z-ipseed'+ext, 'energy-z-ipseed'+ext,
                  'spectrum-lamds-ipseed'+ext, 'power-t-ipseed'+ext, 'peak-power-z-ipseed'+ext])
    count = 0
    for i in np.arange(start, end):
        try:
            ipseed = i+0; print(ipseed, end = ' ')
            fname = '%s.%d.out.h5' % (rootname, ipseed)
            # fname = 'pithz.%dA.%.1fps.%03d.out' % (Ipeak, FWHM, ipseed)
            
            pg = PostGenesis(fname, version = version)

            zplot = pg.zplot
            zpower = pg.zpower
            zenergy = pg.zenergy
            
            ppower = np.max(pg.get_fielddata('power'), axis = 1)
            
            tt = pg.zbunch/g_c
            #z0 = 1.2; 
            ss = '_%.2fm' % z0
            tspec = pg.get_fielddata('power', at = z0)
            ww, fspec = pg.get_spectrum(at = z0)
            
            count += 1
            if count == 1:
                with open('./power-z-ipseed'+ext, 'w') as f_handle:
                    np.savetxt(f_handle,np.atleast_2d(zplot),fmt='%15.6E')
                with open('./energy-z-ipseed'+ext, 'w') as f_handle:
                    np.savetxt(f_handle,np.atleast_2d(zplot),fmt='%15.6E')
                with open('./peak-power-z-ipseed'+ext, 'w') as f_handle:
                    np.savetxt(f_handle,np.atleast_2d(zplot),fmt='%15.6E')
                    
                with open('./spectrum-lamds-ipseed'+ss+ext, 'w') as f_handle:
                    np.savetxt(f_handle,np.atleast_2d(ww),fmt='%15.6E')
                with open('./power-t-ipseed'+ss+ext, 'w') as f_handle:
                    np.savetxt(f_handle,np.atleast_2d(tt),fmt='%15.6E')
                
            with open('./power-z-ipseed'+ext, 'a') as f_handle:
                np.savetxt(f_handle,np.atleast_2d(zpower),fmt='%15.6E')
            with open('./energy-z-ipseed'+ext, 'a') as f_handle:
                np.savetxt(f_handle,np.atleast_2d(zenergy),fmt='%15.6E')
            with open('./peak-power-z-ipseed'+ext, 'a') as f_handle:
                np.savetxt(f_handle,np.atleast_2d(ppower),fmt='%15.6E')
                
            with open('./spectrum-lamds-ipseed'+ss+ext, 'a') as f_handle:
                np.savetxt(f_handle,np.atleast_2d(fspec),fmt='%15.6E')
            with open('./power-t-ipseed'+ss+ext, 'a') as f_handle:
                np.savetxt(f_handle,np.atleast_2d(tspec),fmt='%15.6E')
                
        except Exception as err:
            print(err)
            print('Error got at: ', i, )
            pass
    print('\n')
Analyze(1, 101, 3.6, '.txt', rootname = 'g4', version = 4)
#for z0 in np.linspace(0.3, 3.6, 12):
    #Analyze(1, 21, z0, '.txt', rootname = 'pithz', version = 3)
#exit()

#%% Astra2hdf5
astra2hdf5()
#%% Load an output file
workdir = r'\\afs\ifh.de\group\pitz\data\lixiangk\sim\2021\Genesis-demo\Minimal'
#workdir = r'\\afs\ifh.de\group\pitz\data\lixiangk\sim3\2024\Biaobin'
#os.chdir(workdir)

fname = 'g4.1.out.h5'
version = 4
pg = PostGenesis(fname, debug = 1, version = version, fig_ext = '.png')

#%%% Plot electron beam size along the undulator
z = pg.zplot
current = pg.current[1:]

if version < 4: # version 2 or 3
    xx = pg.get_fielddata('xrms')
    yy = pg.get_fielddata('yrms')
    xrms = (xx[:,1:] @ current)/np.sum(current); #print(xrms)
    yrms = (pg.get_fielddata('yrms')[:,1:] @ current)/np.sum(current)
    
    #xrms = np.mean(xx[:,1:], axis = 1); #print(xrms)
    #yrms = np.mean(yy[:,1:], axis = 1); #print(xrms)

else: # version 4
    nn = np.sum(current>0)
    xrms = (pg.get_data('Beam', 'xsize')[:,:nn] @ current[:nn])/np.sum(current)
    yrms = (pg.get_data('Beam', 'ysize')[:,:nn] @ current[:nn])/np.sum(current)

fig, ax = plt.subplots(figsize = (5, 3))

ax.plot(z, xrms*1e3, 'r-')
ax.plot(z, yrms*1e3, 'b-')
ax.plot(z, np.sqrt(xrms*yrms)*1e3, 'g-')

ax.grid()

ax.set_xlabel(r'$z$ (m)')
ax.set_ylabel(r'RMS size (mm)')
ax.legend([r'$x$', r'$y$'])

#%%% Plot the radiation size along the undulator
zplot = pg.zplot
current = pg.current

#name = 'r_size' # version 2
name = 'xsize' # or 'ysize' for version 3
r_size = (pg.get_fielddata(name) @ current)/np.sum(current)

fig, ax = plt.subplots(figsize = (5, 3))

ax.plot(zplot, r_size*1e3, 'r-')
ax.grid()

ax.set_xlabel(r'$z$ (m)')
ax.set_ylabel(r'Radiation trans. size (mm)')
#%% Comapre near field spectrum and far field spectrum, version 3?
#pg = PostGenesis('test_input_slices.out.h5', debug = 0, version = 4, fig_ext = '.png')

fig, [ax1, ax2] = plt.subplots(ncols = 2, figsize = (6, 4))
temp = pg.get_fielddata('intensity-nearfield')[-1]
ax1.plot(temp/np.max(temp))

temp = pg.get_fielddata('intensity-farfield')[-1]
ax1.plot(temp/np.max(temp))

temp = pg.get_fielddata('power')[-1]
ax1.plot(temp/np.max(temp))

amp = pg.get_fielddata('intensity-nearfield')[-1]
phase = pg.get_fielddata('phase-nearfield')[-1]

ww, ss = calcSpectrum(amp, phase)

ax2.plot(ww*1e6, ss/np.max(ss))

amp = pg.get_fielddata('intensity-farfield')[-1]
phase = pg.get_fielddata('phase-farfield')[-1]
ww, ss = calcSpectrum(amp[67:], phase[67:])

ax2.plot(ww*1e6, ss/np.max(ss))

fig.tight_layout()
#%% Read the 3D fields, version 3
import h5py

#f1 = h5py.File('pithz.211.out.dfl.h5', 'r')
f0 = h5py.File('pithz.2111.out.fld.h5', 'r')
f1 = f0['step000240']

nx, ny, nz = 251, 251, 187
shape = (nx, ny, 2)
fields = np.zeros((nx, ny, nz, 2))
for islice in np.arange(0, nz):
    slicename = 'slice%06d/field' % (islice+1)
    temp = f1[slicename][:]
    temp = np.reshape(temp, shape)
    
    fields[:,:,islice,:] = temp[...]
complex_field = fields[:,:,:,0]+1j*fields[:,:,:,1]

### Plot 2D projection of the fields
profile2d = np.sum(np.abs(complex_field)**2, axis=2)
plt.figure()
plt.imshow(profile2d)

#%% Load fields saved in version 4

i = 240
fname = f'g4.21.{i:d}.fld.h5'
f1 = h5py.File(fname, 'r')

#fields = np.zeros((nx, ny, nz, 2))
#intensity_xy = np.zeros((nx, ny))

# shape = (nx, ny)
# for islice in np.arange(0, nz):
#     slicename = 'slice%06d/field-real' % (islice+1)
#     temp1 = f1[slicename][:]
#     temp1 = np.reshape(temp1, shape)
    
#     #fields[:,:,islice,0] = temp[...]
    
#     slicename = 'slice%06d/field-imag' % (islice+1)
#     temp2 = f1[slicename][:]
#     temp2 = np.reshape(temp2, shape)
    
#     #fields[:,:,islice,1] = temp[...]
    
#     intensity_xy += np.abs(temp1+1j*temp2)**2
#     print(islice+1, end = ', ')
# #complex_field = fields[:,:,:,0]+1j*fields[:,:,:,1]


### Plot 2D projection of the fields
#profile2d = np.sum(np.abs(complex_field)**2, axis=2)
#plt.figure()
#plt.imshow(profile2d)

nx, nz = f1['gridpoints'][0], f1['slicecount'][0]
ny = nx

dx = f1['gridsize'][0]
dy = dx

xx = np.arange(nx)*dx*1e3 # mm
xx -= xx.mean()

yy = xx
aX, aY = np.meshgrid(xx, yy)

int_xy = f1['int_xy'][:].reshape((nx, ny))
#%
fig, ax = plt.subplots(figsize = (5,4))
#cax = ax.pcolormesh(aX, aY, int_xy.T, shading = 'auto')
int_max = int_xy.max()
select = int_xy<int_max*0.001

select = (np.abs(aX)<=5.5)*(np.abs(aY)<=2.5)
tmp = np.copy(int_xy).T

#tmp[~select] = np.nan
tmp = tmp/int_max
#tmp[tmp<0.005] = np.nan

cax = ax.imshow(tmp, extent = [xx.min(), xx.max(), yy.min(), yy.max()])

cbar = fig.colorbar(cax, fraction = 0.046, pad = 0.04)
#cbar.set_ticks(np.linspace(0, 900, 10))
cbar.ax.tick_params(labelsize = 10, pad = 2, rotation = 60)
cbar.ax.set_ylabel('Intensity (arb. unit)', labelpad = 3)

ax.set_aspect('equal')
#ax.set_rasterization_zorder(-1)

#ax.grid()
#ax.set_title(title2, fontsize = 11)
ax.set_xlabel(r'$x$ (mm)', fontsize = 12)
ax.set_ylabel(r'$y$ (mm)', fontsize = 12)
mpl.rc('xtick', labelsize = 11, direction = 'in', top   = True)
mpl.rc('ytick', labelsize = 11, direction = 'in', right = True)

xmax = xx.max();
#xmax = 2
ax.set_xlim(-xmax, xmax)
ax.set_ylim(-xmax, xmax)
ax.set_yticks(np.linspace(-xmax, xmax, 5))

fig.savefig(f'int_xy_{i:d}.png')

#%% Load particles and make 2D plot of LPS, version 4
#f2 = h5py.File('g4.21.240.par.h5', 'r')
f2 = h5py.File('g4.101.0.par.h5', 'r')
#%
slicecount = f2.get('slicecount')[0]
slicelength = f2.get('slicelength')[0]
slicespacing = f2.get('slicespacing')[0]
refposition = f2.get('refposition')[0]

fig, ax = plt.subplots(ncols = 1, figsize = (5, 4))

nps = 8192*4
nz = 88
zall = np.zeros(nps*nz)
gall = np.zeros(nps*nz)
wall = np.zeros(nps*nz)

res = []
for i in np.arange(nz):
    
    s1 = f2.get('slice%06d' % (i+1))
    theta1 = s1.get('theta')[:]
    zz = (theta1/np.pi/2+0.5)*slicelength+refposition
    gamma = s1.get('gamma')[:]
    current = s1.get('current')[0]
    
    bf = np.mean(np.exp(-1j*theta1))
    res.append([refposition, np.abs(bf), current]) 
    
    range = [[refposition*1e3, 1e3*(refposition+slicespacing)], [32, 35]]
    ax.hist2d(zz*1e3, gamma, bins = (4, 15), range = range, cmin = 1e-9, 
               weights = current*np.ones(zz.shape))
    
    refposition += slicespacing
    
    zall[i*nps:i*nps+nps] = zz[::1]
    gall[i*nps:i*nps+nps] = gamma[::1]
    wall[i*nps:i*nps+nps] = current
    
ax.set_xlim(0, refposition*1e3)

#%
xx = (zall-zall.mean())*1e3
yy = gall*g_mec2

weights = wall

flag = 'zpz'
ymin, ymax, xmin, xmax = -5, 5, 16, 18
xmin, xmax, ymin, ymax = -5, 5, 16, 18
xbin, ybin = 720, 40

bins = int(xmax*10); bins = 400
bins = [xbin, ybin]

range = [[xmin, xmax], [ymin, ymax]]
extent = [xmin, xmax, ymin, ymax]
cmin = 1

#range = None

fig = plt.figure(figsize=(6, 4))

l, b, w, h = 0.15, 0.15, 0.75, 0.75
l, b, w, h = 0.135, 0.15, 0.70, 0.70
rectC = [l, b, w, h]
rectL = [l, b, 0.20, h]
rectB = [l, b, w, 0.20]
rectCB= [l+w, b, 0.03, h]

ax  = fig.add_axes(rectC); ax.patch.set_alpha(0)
axl = fig.add_axes(rectL); axl.patch.set_alpha(0)
axb = fig.add_axes(rectB); axb.patch.set_alpha(0)
axc = fig.add_axes(rectCB);axc.patch.set_alpha(0)

hh, xedges, yedges = np.histogram2d(xx, yy, bins = bins, range = range, 
                                    weights = weights); 
print(np.sum(hh))
hmax = hh.max(); hh /= hmax; hh[hh==0] = np.nan
im = ax.imshow(hh.T[::-1], extent = extent, aspect = 'auto', interpolation='none'); 
hh[hh==np.nan] = 0

ax.set_xlim(xmax, xmin)
ax.set_ylim(ymin, ymax)


### colorbar
#ax4.axis('off')
cbar = fig.colorbar(im, cax = axc)#, fraction = 0.046, pad = 0.04)
#cbar.set_lim(0, 1000)
#cbar.set_ticks(np.linspace(0, 1000, 6))
cbar.ax.tick_params(labelsize = 10, pad = 2, rotation = 60)
#cbar.ax.set_ylim(0, 1000)
cbar.ax.set_ylabel('Intensity (arb. unit)', fontsize = 12)
###

#%
ticks = np.linspace(xmin, xmax, 5)
#ticks = [-0.5, 0, 0.5, 1]
#ax.set_xticks(ticks)
ticks = np.linspace(ymin, ymax, 5)
#ticks = [-0.5, 0, 0.5]
#ax.set_yticks(ticks)
if xmax==ymax:
    print('Apply equal aspect ratio...')
    #ax.set_aspect('equal')
    
if flag == 'xy':
    ax.set_xlabel(r'$x$ (mm)')
    ax.set_ylabel(r'$y$ (mm)', labelpad = -10)

if flag == 'zpz':
    ax.set_xlabel(r'$z$ (mm)')
    #ax.set_ylabel(r'$\Delta P_z$ (keV/c)')
    ax.set_ylabel(r'$P$ (MeV/c)')
    
if flag == 'rpz':
    ax.set_xlabel(r'$r$ (mm)')
    ax.set_ylabel(r'$\Delta P_z$ (keV/c)')

# xx 1d hist
hx1 = axb.hist(xx, bins = xbin, histtype = r'step', density = True, weights = weights,
          range = [xmin, xmax])
#axb.set_xlim(xmax, xmin)
#axb.set_xticks(np.linspace(15.5, 17.5, 5)[::-1])

# yy 1d hist
axl.hist(yy, bins = ybin, histtype = r'step', density = True, weights = weights,
          range = [ymin, ymax], orientation = 'horizontal')
axl.set_ylim(ymin, ymax)
axl.set_xlim(0, 3)

axl.axis('off')
axb.axis('off')

###
ax.set_xlim(xmax, xmin)
ax.set_ylim(ymin, ymax)
ax.set_xticks(np.linspace(xmin, xmax, 5)[::-1])
ax.set_yticks(np.linspace(ymin, ymax, 5))
axb.set_xlim(xmax, xmin)
axl.set_ylim(ymin, ymax)
###

plt.tight_layout()

fig.savefig('zpz-240_.png')

#%% Comapre the on-axis field amp and phase to that saved in the general output file
fig, [ax1, ax2] = plt.subplots(nrows = 2, figsize = (4, 6), sharex = True)

temp = pg.get_fielddata('p_mid')[-1]
ax1.plot(temp/np.max(temp))

temp = np.abs(complex_field[nx//2,nx//2,:])**2
ax1.plot(temp/np.max(temp))

ax2.plot(pg.get_fielddata('phi_mid')[-1])
ax2.plot(np.angle(complex_field[nx//2,nx//2,:]))


#%% Check on-axis and off-axis spectra
wavelength, spectrum = calcSpectrum(complex_field[125,125,:])
plt.figure();plt.plot(wavelength*1e6, spectrum)

wavelength, spectrum = calcSpectrum(complex_field[136,125,:])
plt.plot(wavelength*1e6, spectrum)

plt.plot(wavelength*1e6, np.mean(np.mean(spectrums,axis = 0), axis = 0))

plt.figure()
plt.plot(np.abs(complex_field[125,125,:])**2)
plt.plot(np.abs(complex_field[136,125,:])**2)

intensity = np.abs(complex_field)**2
intensity_z = np.sum(np.sum(intensity, axis = 0), axis = 0)/251/251*4

plt.plot(intensity_z)

#%% Apply 3D FFT to 3D fields, just for test

from numpy.fft import fftshift,fft,fftn

axis = 2
spectrums = np.abs(fftshift(fft(complex_field, nz, axis), axis))
#spectrums = np.abs(fftshift(fftn(complex_field)))

spectrums = spectrums*spectrums

sum_z = np.sum(spectrums, axis=2)

from scipy.signal import find_peaks

# Find peaks in the 2D spectrum
peaks_x, _ = find_peaks(sum_z.max(axis=1))
peaks_y, _ = find_peaks(sum_z.max(axis=0))

# Main frequency components (assuming the strongest peak is the main component)
main_freq_x = peaks_x[np.argmax(sum_z[peaks_x, :].max(axis=1))]
main_freq_y = peaks_y[np.argmax(sum_z[:, peaks_y].max(axis=0))]

### Make plots
plt.figure();
wavelength, temp = calcSpectrum(complex_field[nx//2,ny//2,:])
plt.plot(wavelength*1e6, temp/np.max(temp))

temp = np.mean(np.mean(spectrums,axis = 0), axis = 0)
plt.plot(wavelength*1e6, temp/np.max(temp))

temp = spectrums[nx//2, ny//2, :]
plt.plot(wavelength*1e6, temp/np.max(temp))

temp = spectrums[nx//2, 114, :]
plt.plot(wavelength*1e6, temp/np.max(temp))

temp = spectrums[nx//2, 136, :]
plt.plot(wavelength*1e6, temp/np.max(temp))

# wavelength, temp = calcSpectrum(complex_field[main_freq_x,main_freq_y,:])
# plt.plot(wavelength*1e6, temp/np.max(temp))

#%% Load fields dumped after the simulation, version 3
import h5py

f1 = h5py.File('pithz.2111.out.dfl.h5', 'r')

#%% Load an output from version 2, for comparison
workdir = r'\\afs\ifh.de\group\pitz\data\lixiangk\sim3\2024\DLW-XYZ\beam_with_BC'
workdir = r'\\afs\ifh.de\group\pitz\data\lixiangk\sim\2021\G4-SC-test'
os.chdir(workdir)

fname = 'pithz.211.out.h5'
fname = 'test_sc_on_v2.out'
pg = PostGenesis(fname, debug = 1, version = 2, fig_ext = '.png')

#%% Animation of field amp and phase along the undulator
from IPython import display
fig, [ax1, ax2] = plt.subplots(nrows = 2, figsize = (4, 6), sharex = True)
plt.show()

if version == 3:
    p_mid = pg.get_fielddata('p_mid')
    phi_mid = pg.get_fielddata('phi_mid')

if version == 4:
    name1 = 'intensity-farfield'
    name2 = 'phase-farfield'
    
    ppower = pg.get_fielddata('power')
    p_mid = pg.get_fielddata(name1)
    phi_mid = pg.get_fielddata(name2)
    gamma = pg.get_beamdata('energy')
    bunching = pg.get_beamdata('bunching')


temp = p_mid
step = 10
for i in np.arange(10, temp.shape[0], step):
    
    ax1.plot(p_mid[i]/np.sum(p_mid[i]), '-')

    ax1.set_title('z = : %.3f m' % pg.zplot[i])
    # ax1.set_xlabel('# of slice')
    # ax1.set_ylabel('Power (MW)')
    # #ax.set_ylim(0, ymax)
    # ax1.grid()
    #ax1.set_xlim(0, 120)
    
    ax1.set_ylabel(r'Far field amp. (arb. unit)')
    
    ax2.plot(phi_mid[i], '-')
    # ax2.set_xlabel(r'$z$ (m)')
    # ax2.set_ylabel(r'Peak power (MW)')
    # ax2.grid()
    ax2.set_ylabel(r'Far field phase')
    
    arg_max = np.argmax(p_mid[i])
    #ax2.set_xlim(0, 250)
    if arg_max<60 or i//step <= 3:
        arg_max = 60
    ax2.set_xlim(arg_max-60, arg_max+60)
    
    display.display(plt.gcf())
    display.clear_output(wait=True)
    
    fig.tight_layout()
    fig.savefig('phase@z_%.2fm' % pg.zplot[i] + '.png')    
    
    plt.pause(1)
    if i<temp.shape[0]-step:
        ax1.cla()
        ax2.cla()
        
#%% Animation of field amp, phase, energy loss and spectrum
from IPython import display
fig, [ax1, ax2, ax3, ax4, ax5] = plt.subplots(nrows = 5, figsize = (4, 10))
ax12 = ax1.twinx()
ax32 = ax3.twinx()

plt.show()

name1 = 'intensity-farfield'
name2 = 'phase-farfield'

ppower = pg.get_fielddata('power')
p_mid = pg.get_fielddata(name1)
phi_mid = pg.get_fielddata(name2)
gamma = pg.get_beamdata('energy')
bunching = pg.get_beamdata('bunching')

step = 10
for i in np.arange(0, ppower.shape[0], step):
    
    ax1.plot(p_mid[i]/p_mid[i].max(), 'r-')
    ax12.plot(ppower[i]/ppower[i].max(), 'b-')
    
    title = 'z = %.3f m' % pg.zplot[i]
    ax1.set_title(title)
    ax1.set_xlabel(r'# of slice')
    ax1.set_ylabel(r'On-axis power (au)')
    #ax12.set_ylabel(r'Power (au)')
    #ax.set_ylim(0, 2500)
    # ax1.grid()
    ax1.set_xlim(0, 150)
    ax1.set_ylim(0., 1.1)
    #ax1.set_yscale('log')
    
    ax12.set_ylim(0.0, 1.1)
    #ax12.set_yscale('log')
    
    ax2.plot(gamma[0], 'r-')
    ax2.plot(gamma[i], 'b-')
    ax2.set_xlabel('# of slice')
    ax2.set_ylabel(r'$\gamma$')
    # ax2.grid()
    ax2.set_ylim(31, 35)
    ax2.set_xlim(0, 150)
    
    DE = (gamma[i]-gamma[0])*pg.current
    #ax3.plot(DE/np.max(np.abs(DE)), 'r-')
    ax3.plot(DE, 'r-')
    ax3.set_xlim(0, 150)
    ax3.set_ylim(-2, 2)
    ax3.set_xlabel('# of slice')
    ax3.set_ylabel(r'Energy loss')
    
    ax32.plot(pg.current, 'k-')
    #ax32.set_ylabel(r'Current (A)')
    
    ax4.plot(bunching[i], 'r-')
    ax4.set_xlabel('# of slice')
    ax4.set_ylabel('Bunching')
    ax4.set_ylim(1e-5, 1)
    ax4.set_yscale('log')
    ax4.set_xlim(0, 150)
    
    w, s = pg.get_spectrum(pg.zplot[i])
    #ax5.plot(w*1e6, s/np.max(s), '-')
    ax5.plot(w*1e6, s, '-')
    
    # w, s = calcSpectrum(p_mid[i], phi_mid[i])
    # ax5.plot(w*1e6, s/np.max(s), '-')
    
    ax5.set_xlim(80, 120)
    #ax5.set_ylim(0, 1.1)
    ax5.set_ylim(0, 2e13)
    
    ax5.set_xlabel(r'Wavelength ($\mu$m)')
    ax5.set_ylabel(r'Intensity (arb. unit)')
    
    display.display(plt.gcf())
    display.clear_output(wait=True)

    fig.tight_layout()
    fig.savefig('z_%.2fm' % pg.zplot[i] + '.png')    
    
    plt.pause(1)
    if i<ppower.shape[0]-step:
        ax1.cla()
        ax12.cla()
        ax2.cla()
        ax3.cla()
        ax32.cla()
        ax4.cla()
        ax5.cla()
        
#%% Check the energy loss along the slices and so on
gamma = pg.get_beamdata('energy')
ppower = pg.get_fielddata('power')

current = pg.current
zenergy = pg.zenergy

slices = [40, 41, 42]

fig, ax = plt.subplots(nrows = 3, figsize = (5, 10), sharex = True)
ax[0].plot(gamma[0], '-')
ax[0].plot(gamma[-1], '-')
ax[0].set_ylabel(r'$\gamma$')
ax[0].legend(['Before und.', 'After und.'])

steps = [10, 120, -1]
for istep in steps:
    ax[1].plot(ppower[istep]/np.max(ppower[istep]), '-')
ax[1].legend(['Step = %d' % i for i in steps])

ax[2].plot([], [], 'r-')
ax[2].plot([], [], 'b-')

ax[2].plot(current, 'r-')
ax2 = ax[2].twinx()
ax2.plot((gamma[-1]-gamma[0])*current, 'b-')
ax[2].legend(['Current', 'Energy change'])

for axis in ax:
    axis.grid()

ax[2].set_xlabel(r'# of slice')

fig.savefig('Along-slices.png')

NN = pg.zplot/0.03 # number of undulator period

fig, ax = plt.subplots(nrows = 4, figsize = (5, 10), sharex = True)

# Slice energy along the undulator
for islice in slices:
    ax[0].plot(NN, gamma[:,islice], '-')
ax[0].legend(['%0d' % islice for islice in slices])

# Fit the slope of slice energy
islice = 21
select = (NN>20)*(NN<44) # 41
linear = lambda x, a, b:a+b*x
popt, pcov = curve_fit(linear, NN[select], gamma[select,islice]); print(popt)

#ax[0].plot(NN, linear(NN, *popt), '-')


ax[1].plot(NN, zenergy/1e6, '-')
ax[1].set_ylabel(r'$E$ ($\mu$J)')

try:
    aa = pg.zaw; aa = np.hstack(([aa[0]], aa))
except Exception as err:
    aa = pg.field[:,-1]
    
ccc1 = np.sqrt(2)
#ccc1 = 1
# Resonant wavelength change due to energy change
for islice in slices:
    l1 = resonant_wavelength(aa*ccc1, 3e-2, gamma[:,islice])*1e6
    ax[2].plot(NN, l1, '-') 
ax[2].set_ylim(90, 110)
ax[2].legend(['%0d' % islice for islice in slices])

islice = slices[0]
a2 = resonant_undulator_parameter(100e-6, 4e-2, gamma[:,islice]+pg.gamma0)
ax[3].plot(NN, aa/aa.max(), '-')
ax[3].plot(NN, a2/a2.max(), '-')

# Print the "tapered" undulator strength
for i in np.arange(len(a2)//2):
    #print(i, a2[i*2]/a2.max()*2.4678)
    ss = str.format('%-6s%12.4E%12.4E%6d' % ('AW', a2[i*2]/a2.max()*2.4678, 3e-2, 1))
    #print(ss)

for axis in ax:
    axis.grid()
ax[3].set_xlabel(r'# of undulator period')

fig.savefig('Along-undulator.png')

gamma_ = np.mean(gamma[:,:] @ current)/np.sum(current)

#%% Check the energy loss and beam size along the slices and so on
gamma = pg.get_beamdata('energy')
ppower = pg.get_fielddata('power')
ysize = pg.get_beamdata('ysize')

current = pg.current
zenergy = pg.zenergy

slices = [30, 41, 50, 61]

fig, ax = plt.subplots(nrows = 3, figsize = (5, 10), sharex = True)
ax[0].plot(gamma[0], '-')
ax[0].plot(gamma[-1], '-')
ax[0].set_ylabel(r'$\gamma$')
ax[0].legend(['Before und.', 'After und.'])

steps = [10, 120, -1]
for istep in steps:
    ax[1].plot(ppower[istep]/np.max(ppower[istep]), '-')
ax[1].legend(['Step = %d' % i for i in steps])

ax[2].plot([], [], 'r-')
ax[2].plot([], [], 'b-')

ax[2].plot(current, 'r-')
ax2 = ax[2].twinx()
ax2.plot((gamma[-1]-gamma[0])*current, 'b-')
ax[2].legend(['Current', 'Energy change'])

for axis in ax:
    axis.grid()

ax[2].set_xlabel(r'# of slice')

fig.savefig('Along-slices.png')

NN = pg.zplot/0.03 # number of undulator period

fig, ax = plt.subplots(nrows = 4, figsize = (5, 10), sharex = True)

# Slice energy along the undulator
for islice in slices:
    ax[0].plot(NN, gamma[:,islice], '-')
    ax[3].plot(NN, ysize[:,islice], '-')
ax[0].legend(['%0d' % islice for islice in slices])

# Fit the slope of slice energy
islice = 21
select = (NN>20)*(NN<44) # 41
linear = lambda x, a, b:a+b*x
popt, pcov = curve_fit(linear, NN[select], gamma[select,islice]); print(popt)

#ax[0].plot(NN, linear(NN, *popt), '-')


ax[1].plot(NN, zenergy/1e6, '-')
ax[1].set_ylabel(r'$E$ ($\mu$J)')

try:
    aa = pg.zaw; aa = np.hstack(([aa[0]], aa))
except Exception as err:
    aa = pg.field[:,-1]
    
ccc1 = np.sqrt(2)
#ccc1 = 1
# Resonant wavelength change due to energy change
for islice in slices:
    l1 = resonant_wavelength(aa*ccc1, 3e-2, gamma[:,islice])*1e6
    ax[2].plot(NN, l1, '-') 
ax[2].set_ylim(90, 110)
ax[2].legend(['%0d' % islice for islice in slices])

islice = slices[0]
a2 = resonant_undulator_parameter(100e-6, 4e-2, gamma[:,islice]+pg.gamma0)
#ax[3].plot(NN, aa/aa.max(), '-')
#ax[3].plot(NN, a2/a2.max(), '-')

# Print the "tapered" undulator strength
for i in np.arange(len(a2)//2):
    #print(i, a2[i*2]/a2.max()*2.4678)
    ss = str.format('%-6s%12.4E%12.4E%6d' % ('AW', a2[i*2]/a2.max()*2.4678, 3e-2, 1))
    #print(ss)

for axis in ax:
    axis.grid()
ax[3].set_xlabel(r'# of undulator period')

fig.savefig('Along-undulator.png')

#%% 3D visualization
#%%% Plot the spectrum along a third axis of z

from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure(figsize = (5, 4))
ax = fig.add_subplot(111, projection='3d')

zplot = pg.zplot
zenergy = pg.zenergy
ax.plot(zplot, zenergy/zenergy.max(), zs = 0.9, zdir = 'x', linestyle = '-', 
        marker = 'none')
# Plot each 1D plot at a different time (z-axis)
for i, zi in enumerate(zplot[::40][::-1]):
    w, s = pg.get_spectrum(zi)
    i = int(zi/0.015)
    
    select = (100/(w*1e6)>0.9)*(100/(w*1e6)<1.1)
    smax = 60742192099236.7
    ax.plot(100/(w[select]*1e6), s[select]/s.max()*i/len(zplot), zs=zi, 
            zdir='y', linestyle = '-', marker = 'none', label=f'$z$={zi:.1f}m')

ax.set_xlim(0.9, 1.1)
#ax.set_xticks([0.8, 0.9, 1, 1.1, 1.2])
ax.set_xlabel('$f/f_0$', labelpad = 5)
ax.set_ylabel(r'$z$ (m)', labelpad = 3)
ax.set_zlabel('Norm. intensity', fontsize = 11, labelpad = 3)
ax.set_yticks([0, 1, 2, 3, 4, 5])

fig.tight_layout()
fig.savefig('Spectrum_along_undulator.png')
#%%% Plot the energy losses along a third axis of z
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure(figsize = (5, 4))
ax = fig.add_subplot(111, projection='3d')

current = pg.current
zplot = pg.zplot
zenergy = pg.zenergy
gamma = pg.get_beamdata('energy')
slices = np.arange(len(current))

ax.plot(zplot, zenergy/zenergy.max(), zs = 0, zdir = 'x', linestyle = '-', 
        marker = 'none')
# Plot each 1D plot at a different time (z-axis)
for i, zi in enumerate(zplot[::-40]):
    i = int(zi/0.015)
    ax.plot(slices, (gamma[i]-gamma[0])*current, zs=zi, zdir='y', 
            linestyle = '-', marker = 'none', label=f'$z$={zi:.1f}m')

#ax.set_xlim(80, 120)
ax.set_xlabel('# of slices')
ax.set_ylabel(r'$z$ (m)')
ax.set_zlabel('Energy loss')

fig.tight_layout()
#%%% Plot the THz temeporal profiles along a third axis of z
name = 'power'
name = 'intensity-farfield'
ppower = pg.get_fielddata(name)

fig = plt.figure(figsize = (5, 4))
ax = fig.add_subplot(111, projection='3d')

current = pg.current
zplot = pg.zplot
zenergy = pg.zenergy
gamma = pg.get_beamdata('energy')
slices = np.arange(len(current))

xdata = slices
xlabel = '# of slices'
use_time = 1
if use_time:
    xdata = (slices*100e-6/g_c*1e12)
    xlabel = r'Time (ps)'

select = (xdata<=40)
xdata = xdata[select]
ppower = ppower[:,select]
current = current[select]
    
pmax = ppower.max()
ax.plot(xdata, current/current.max(), zs = zplot[-1], zdir = 'y', 
        linestyle = '--', color = 'k', marker = 'none')
ax.plot(zplot, zenergy/zenergy.max(), zs = 0, zdir = 'x', linestyle = '-', 
        marker = 'none')

# Plot each 1D plot at a different time (z-axis)
for i, zi in enumerate(list(zplot[:1:-40])+[0.15]):
#for i, zi in enumerate([0.15, 0.6, 1.2, 1.8, 2.4, 3.0, 3.6]):
#for i, zi in enumerate([0.15, 0.9, 1.8, 2.7, 3.6, 4.5, 5.4]):
    i = int(zi/0.015); print(zi)
    ax.plot(xdata, ppower[i]/ppower[i].max()*(i/len(zplot))**0.3, zs=zi, zdir='y',
            linestyle = '-', marker = 'none', label=f'$z$={zi:.1f}m')

#ax.set_zscale('log')

ax.set_xlim(0, 40)#xdata[-1])
ax.set_xlabel(xlabel, labelpad = 3)
ax.set_ylabel(r'$z$ (m)', labelpad = 3)
ax.set_zlabel('Norm. intensity', fontsize = 11, labelpad = 3)
ax.set_yticks([0, 1, 2, 3, 4, 5])
ax.set_ylim(0, zplot[-1])

fig.tight_layout()
fig.savefig(name+'_along_undulator.png')

#%%% Plot the bunching along a third axis of z
data = pg.get_beamdata('bunching')

fig = plt.figure(figsize = (5, 4))
ax = fig.add_subplot(111, projection='3d')

current = pg.current
zplot = pg.zplot
zenergy = pg.zenergy
gamma = pg.get_beamdata('energy')
slices = np.arange(len(current))

select = (slices<=40)
slices = slices[select]
data = data[:,select]
current = current[select]

ax.plot(zplot, zenergy/zenergy.max(), zs = 0, zdir = 'x', linestyle = '-', 
        marker = 'none')
# Plot each 1D plot at a different time (z-axis)
for i, zi in enumerate(zplot[::-40]):
    i = int(zi/0.015)
    ax.plot(slices, data[i]/data[i].max()*i/len(zplot), zs=zi, zdir='y', 
            linestyle = '-', marker = 'none', label=f'$z$={zi:.1f}m')

#ax.set_zscale('log')

#ax.set_xlim(80, 120)
ax.set_xlabel('# of slices')
ax.set_ylabel(r'$z$ (m)')
ax.set_zlabel('Norm. intensity')

fig.tight_layout()
#%% Plot bunching amp and phase of slices along z
current = pg.current
zplot = pg.zplot; print(zplot.shape)

# print(pg.outputs) # check the name of the property of interest
temp = pg.get_beamdata('bunching')
#temp1 = pg.get_beamdata('phi_mid')
temp1 = pg.get_beamdata('bunchingphase')

fig, [ax1, ax2, ax3] = plt.subplots(nrows = 3, figsize = (4, 6))

ax1.plot(current, 'k-')
ax1.set_xlabel(r'# of slices')
ax1.set_ylabel(r'Current (A)')
ax1.grid()

slices = [40, 41, 42]
#@slices = [45, 47, 49]
for i in slices:
    ax2.plot(zplot, temp[:,i], '--')
    ax3.plot(zplot, temp1[:,i], '--')
    print(temp[0,i])

ax2.set_xlabel(r'$z$ (m)')
ax2.set_ylabel(r'Bunching factor')
ax3.set_ylabel(r'phase')
ax2.legend(['# %d' % i for i in slices])
ax3.legend(['# %d' % i for i in slices])

ax2.grid()
ax3.grid()
ax2.set_yscale('log')
ax2.set_xlim(0, 3.6)
ax3.set_xlabel(r'$z$ (m)')

fig.tight_layout()
#fig.savefig('bunching-factor-vs-z.png')

#%% ANALYZE results from parameter scan
#workdir = os.path.join(simdir, r'LCLS-I-3', 'currentScan9')
#os.chdir(workdir)

### Scan 
currents = np.linspace(50, 200, 4)
currents = [50, 100, 150, 200]
currents = [112]
charges = np.linspace(1.5, 5, 8)
charges = [2]
seeds = np.arange(20)+1
#seeds = [1]

combi = np.array([[v1, v2, v3] for v1 in currents for v2 in charges for v3 in seeds])
combi = [[78, 1], [100, 1.5], [112, 2.0], [180, 4.0], [179, 4.0]]

PE = []
for v in combi:
    
    Ipeak, Qtot = v[:] # A, nC
    direc = 'beam_%.0fA_%.1fnC' % (Ipeak, Qtot)
    os.chdir(direc); print(direc)
    
    try:
        ext = '.txt'
        #ext = '-%dA-%.1fps.txt' % (Ipeak, FWHM)
        pp = np.loadtxt('power-z-ipseed'+ext); pp = pp.T; #print pp.shape
        EE = np.loadtxt('energy-z-ipseed'+ext); EE = EE.T; #print EE.shape

        end = -1
        PE.append([Ipeak, Qtot, np.mean(pp[-1,1:end])/1e6, np.std(pp[-1,1:end])/1e6, \
                  np.mean(EE[-1,1:end])*1e6, np.std(EE[-1,1:end])*1e6]) # 'MW' and 'uJ'
    except Exception as err:
        print(err)
        pass
    os.chdir('..')
        
PE = np.array(PE)

fname = 'power-and-energy-vs-Ipeak-Qtot_20240208.dat'
header = '%14s%16s%16s%16s%16s%16s' % ('current (A)', 'charge (nC)', 
                                       '<Power> (MW)', 'rms_P (MW)',
                                       '<Energy> (uJ)', 'rms_E (uJ)')
np.savetxt(fname, PE, fmt = '%15.6E', header = header)

tr = 0.58
PE_exp = np.array([[1, 6.118, 0.784],
                   [1.5, 15.153, 1.019],
                   [2, 21.436, 2.129]])

PE_exp1 = np.array([[2, 33.517, 3.574]])

fig, ax = plt.subplots()
ax.errorbar(PE[:3,1], PE[:3,4], yerr = PE[:3,5], 
            capsize = 4, fmt = 'r-*')
ax.errorbar(PE_exp[:,0], PE_exp[:,1]/tr, yerr = PE_exp[:,2]/tr, 
            capsize = 4, fmt = 'b-*')
ax.errorbar(PE_exp1[:,0], PE_exp1[:,1]/tr, yerr = PE_exp1[:,2]/tr, 
            capsize = 4, fmt = 'bD')

ax.legend(['Simulation, SASE', 'Measured, SASE', 
           'Simulation, Seeded', 'Measured, Seeded'],
          labelspacing = 1, fontsize = 12, loc = 'upper left')

ax.set_xlabel(r'$Q$ (nC)')
ax.set_ylabel(r'Pulse energy ($\mu$J)')
ax.set_ylim(0, 500)

fig.savefig('energy-vs-charge_SASE_vs_Seeded.png')

#%% 2D scan with Apple II undulator
var1 = [17, 22]
#var1 = [17]
var2 = [20e-6, 30e-6, 60e-6, 100e-6]
#var2 = [100e-6]
seeds = np.arange(20)+1
#seeds = [1]

combi = np.array([[v1, v2] for v1 in var1 for v2 in var2])
#combi = [[78, 1], [100, 1.5]]
#combi = [[112, 2]]

PE = []
for v in combi:
    P0, lambda0 = v[:]
    direc = 'P0_%.0fMeV_c_lam_%.0fum_opt' % (P0, lambda0*1e6)
    
    os.chdir(direc); print(direc)
        
    try:
        ext = '.txt'
        #ext = '-%dA-%.1fps.txt' % (Ipeak, FWHM)
        pp = np.loadtxt('power-z-ipseed'+ext); pp = pp.T; #print pp.shape
        EE = np.loadtxt('energy-z-ipseed'+ext); EE = EE.T; #print EE.shape

        end = -1
        PE.append([P0, lambda0, np.mean(pp[-1,1:end])/1e6, np.std(pp[-1,1:end])/1e6, \
                  np.mean(EE[-1,1:end])*1e6, np.std(EE[-1,1:end])*1e6]) # 'MW' and 'uJ'
    except Exception as err:
        print(err)
        pass
    os.chdir('..')
    
PE = np.array(PE)

fname = 'power-and-energy-vs-Ipeak-Qtot_APPLE2_20240212.csv'
header = '%14s%16s%16s%16s%16s%16s' % ('current (A)', 'charge (nC)', 
                                       '<Power> (MW)', 'rms_P (MW)',
                                       '<Energy> (uJ)', 'rms_E (uJ)')
np.savetxt(fname, PE, fmt = '%15.6E', header = header, delimiter = ',')


fig, ax = plt.subplots()
ax.errorbar(PE[:4,1]*1e6, PE[:4,4], yerr = PE[:4,5], 
            capsize = 4, fmt = 'r-*')
ax.errorbar(PE[4:,1]*1e6, PE[4:,4], yerr = PE[4:,5], 
            capsize = 4, fmt = 'b-*')

ax.legend(['17 MeV/c', '22 MeV/c'],
          labelspacing = 1, fontsize = 12, loc = 'upper left')

ax.set_xlabel(r'$Q$ (nC)')
ax.set_ylabel(r'Pulse energy ($\mu$J)')
#ax.set_ylim(0, 500)

fig.savefig('energy-vs-wavelength_APPLE2.png')

#%% 2D scan of peak current and bunch charge
currents = np.linspace(50, 200, 4)
charges = [0.5, 1, 1.5, 2, 3, 4]
charges = [1, 1.5, 2]

fname = 'power-and-energy-vs-Ipeak-Qtot.dat'
PE = np.loadtxt(fname)

fig, ax = plt.subplots()

vars1 = charges
for xi in vars1:
    select = (PE[:,1] == xi)
    ax.errorbar(PE[select,0], PE[select,2], yerr = PE[select,3], capsize = 4)

ax.grid()
ax.set_xlabel(r'$I$ (A)')
ax.set_ylabel(r'$P$ (MW)')
#ax.legend(['cor_Ek = 0.5 %', 'cor_Ek = 1.0%'])
ax.legend(['%.1f nC' % xi for xi in vars1], labelspacing = 0.3)
fig.savefig('power-vs-current_.png')

fig, ax = plt.subplots(figsize = (4, 4))

for xi in vars1:
    select = (PE[:,1] == xi)
    ax.errorbar(PE[select,0], PE[select,4], yerr = PE[select,5], 
                linestyle = '-', capsize = 4)

ax.grid()
ax.set_xlabel(r'$I$ (A)')
ax.set_ylabel(r'$E$ ($\mu$J)')
ax.set_ylim(0, 500)
#ax.legend(['cor_Ek = 0.5 %', 'cor_Ek = 1.0%'])
ax.legend(['%.1f nC' % xi for xi in vars1], labelspacing = 0.3)
fig.savefig('energy-vs-current_.png')

#%% 1D plot: vs charge

fig, ax = plt.subplots()

vars1 = [0]
for xi in vars1:
    select = (PE[:,0] == PE[:,0])
    ax.errorbar(PE[select,1], PE[select,2], yerr = PE[select,3], capsize = 4)

ax.grid()
ax.set_xlabel(r'$Q$ (nC)')
ax.set_ylabel(r'$P$ (MW)')
#ax.legend(['cor_Ek = 0.5 %', 'cor_Ek = 1.0%'])
#ax.legend(['%.0f A' % xi for xi in vars1], labelspacing = 0.3)
fig.savefig('power-vs-charge_.png')

fig, ax = plt.subplots(figsize = (4, 4))

for xi in vars1:
    select = (PE[:,0] == PE[:,0])
    
    ax.errorbar(PE[select,1], PE[select,4], yerr = PE[select,5], capsize = 4, fmt = '-')

ax.grid()
ax.set_xlabel(r'$Q$ (nC)')
ax.set_ylabel(r'$E$ ($\mu$J)')
ax.set_ylim(0, 500)
#ax.legend(['cor_Ek = 0.5 %', 'cor_Ek = 1.0%'])
#ax.legend(['%.1f nC' % xi for xi in vars1], labelspacing = 0.3)
fig.savefig('energy-vs-charge_.png')

fig, ax = plt.subplots()

vars1 = [0]
for xi in vars1:
    select = (PE[:,0] == PE[:,0])
    ax.plot(PE[select,1], PE[select,0])

ax.grid()
ax.set_xlabel(r'$Q$ (nC)')
ax.set_ylabel(r'$I_{\rm peak}$ (A)')
#ax.legend(['cor_Ek = 0.5 %', 'cor_Ek = 1.0%'])
#ax.legend(['%.0f nC' % xi for xi in vars1], labelspacing = 0.3)
fig.savefig('peak-current-vs-charge_.png')

#%% Make plots for batch analysis
#%%% power along the undulator for various seed numbers
ext = '.txt'
pp = np.loadtxt('power-z-ipseed'+ext); pp = pp.T; print (pp.shape)
start = 1
end = pp.shape[1]

output = f'P =: {np.mean(pp[-1,1:end])/1e6:.2f} +/- {np.std(pp[-1,1:end])/1e6:.2f} MW\n'
print(output)

with open('stat.dat', 'w') as f_handle:
    f_handle.write(output)

fig, ax1 = plt.subplots(figsize = (4, 4))
for i in np.arange(start, end):
    ax1.plot(pp[:,0], pp[:,i]/1e6, '-', color='grey') # energy: MW
ax1.plot(pp[:,0], np.mean(pp[:,1:], 1)/1e6, 'k-')

ax1.set_xlabel(r'$z$ (m)')
ax1.set_ylabel(r'power (MW)')
ax1.set_yscale('log')
#ax1.grid()
fig.savefig('ppower-vs-z-ipseed'+ext+'.png', dpi = 300)


#%%% power and energy along the undulator for various seed numbers
ext = '.txt'

pp = np.loadtxt('power-z-ipseed'+ext); pp = pp.T; print (pp.shape)
end = pp.shape[1]

output = f'P =: {np.mean(pp[-1,1:end])/1e6:.2f} +/- {np.std(pp[-1,1:end])/1e6:.2f} MW\n'
print(output)

fig, ax1 = plt.subplots(figsize = (4, 4))
for i in np.arange(1, end):
    ax1.plot(pp[:,0], pp[:,i]/1e6, '-', color='grey') # energy: MW
ax1.plot(pp[:,0], np.mean(pp[:,1:], 1)/1e6, 'k-')

ax1.set_xlabel(r'$z$ (m)')
ax1.set_ylabel(r'power (MW)')
ax1.set_yscale('log')
#ax1.grid()
fig.savefig('power-vs-z-ipseed.png')

#
ext = '.txt'
EE = np.loadtxt('energy-z-ipseed'+ext); EE = EE.T; print (EE.shape)
end = EE.shape[1]

output = f'E =: {np.mean(EE[-1,1:end])*1e6:.2f} +/- {np.std(EE[-1,1:end])*1e6:.2f} uJ\n'
print(output)

fig, ax2 = plt.subplots(figsize = (4, 4))
# ax2.plot([], [], 'r--')
# ax2.plot([], [], 'b--')
# ax2.plot([], [], 'k--')
for i in np.arange(1, end):
    print(i, end = ' ')
    ax2.plot(EE[:,0]-0., EE[:,i]*1e6, '-', color='grey') # energy: uJ
ax2.plot(EE[:,0], np.mean(EE[:,1:], 1)*1e6, 'k-')
#res.append(EE[:,0])
#res.append(np.mean(EE[:,1:], 1))

ax2.set_xlabel(r'$z$ (m)', fontsize = 12)
ax2.set_ylabel(r'Pulse energy ($\mu$J)', fontsize = 12)
ax2.set_yscale('log')
#ax2.set_xlim(0.5, 3.5)
#ax2.set_yticks([1e-7, 1e-6, 1e-6, 1e])
#ax2.grid()

#locmin = mpl.ticker.LogLocator(base = 10.0, subs = (0.2, 0.4, 0.6, 0.8), numticks = 12)
#ax2.yaxis.set_minor_locator(locmin)
#ax2.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())

#ax2.text(0.925, 0.95, '($b$)', horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes)

fig.savefig('energy-vs-z-ipseed.png')

#%%% spectrum for various seed numbers
ext = '_5.40m.txt'
#ext = '.txt'
ss = np.loadtxt('spectrum-lamds-ipseed'+ext); ss = ss.T; print (ss.shape)

end = ss.shape[1]


xlamds = []; width = []
fig, ax = plt.subplots(figsize = (4, 4))
for i in np.arange(1, end-1):
    if i<20:
        ax.plot(ss[:,0]*1e6, ss[:,i]/1e13, '-', color='grey') # energy: MW
    
    xlamds.append([weighted_mean(ss[:,0], ss[:,i])*1e6])
    width.append([weighted_std(ss[:,0], ss[:,i])*1e6])
ax.plot(ss[:,0]*1e6, np.mean(ss[:,1:], 1)/1e13, 'k--')

output =  f'Centre wavelength: {np.mean(xlamds):.2f} +/- {np.std(xlamds):.2f} um\n'
output += f'Spectrum width: {np.mean(width):.2f} +/- {np.std(width):.2f} um\n'
print(output)
ax.set_xlabel(r'$\lambda_s$ ($\mu$m)', fontsize = 12)
ax.set_ylabel(r'Intensity (arb. unit)', fontsize = 12, labelpad = 1)
#ax.set_yscale('log')
ax.set_xlim(90, 120)
#ax.set_xlim(50, 70)
#ax.set_ylim(-0.05, 1.1)
#ax.grid()

#ax.text(0.925, 0.95, '($b$)', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
fig.savefig('spectrum-vs-lamds-ipseed.png')
#%%% temporal profile
def fmt_two_digits(x, pos):
    return f'{x:.2f}'

flag = '_5.40m'
ext = flag+'.txt'
ss = np.loadtxt('power-t-ipseed'+ext); ss = ss.T; print (ss.shape)
end = ss.shape[1]

duration = []
arrival = []
fig, ax = plt.subplots(figsize = (4, 4))
for i in np.arange(1, end-1):
    if i<20:
        ax.plot(ss[:,0]*1e12, ss[:,i]/1e8, '-', color='grey') # energy: MW
    arrival.append([weighted_mean(ss[:,0], ss[:,i])*1e12])
    duration.append([weighted_std(ss[:,0], ss[:,i])*1e12])
ax.plot(ss[:,0]*1e12, np.mean(ss[:,1:], 1)/1e8, 'k--')

output = f'Arrival time: {np.mean(arrival):.2f} +/- {np.std(arrival):.2f} ps\n'
output += f'Pulse duration: {np.mean(duration):.2f} +/- {np.std(duration):.2f} ps\n'

print(output)

with open('stat.dat', 'a') as f_handle:
    f_handle.write(output)

ax.set_xlabel(r'$t$ (ps)')
ax.set_ylabel(r'Intensity (arb. unit)', labelpad = 1)

ax.yaxis.set_major_formatter(lambda x, pos: f'{x:.2f}')
#ax.set_yscale('log')
ax.set_xlim(0, 60)
#ax.set_ylim(-0.005, 0.3)
#ax.grid()
fig.savefig('power-vs-t-ipseed'+flag+'_.png')

#%%% spectrum in THz

flag = '_5.40m'
ext = flag+'.txt'
ss = np.loadtxt('spectrum-lamds-ipseed'+ext); ss = ss.T; print (ss.shape)

end = ss.shape[1]
#end = 20
def lambda2mu(l, v = g_c, f0 = 3e12):
    return v/l/f0
def lambda2mu(l, v = g_c, f0 = 1):
    return v/l/f0

xlamds = []; width = []
fig, ax = plt.subplots(figsize = (4, 4))
for i in np.arange(1, end-1):
    if i<20:
        ax.plot(lambda2mu(ss[:,0])/3e12, ss[:,i]/np.sum(ss[:,i]), '-', color='grey') # energy: MW
    xlamds.append([weighted_mean(lambda2mu(ss[:,0]), ss[:,i])])
    width.append([weighted_std(lambda2mu(ss[:,0]), ss[:,i])])
pmean = np.mean(ss[:,1:], 1)
ax.plot(lambda2mu(ss[:,0])/3e12, pmean/np.sum(pmean), 'k--')

output =  f'Centre wavelength: {np.mean(xlamds)/1e12:.3f} +/- {np.std(xlamds)/1e12:.3f} THz\n'
output += f'Spectrum width: {np.mean(width)/1e12:.3f} +/- {np.std(width)/1e12:.3f} THz\n'
print(output)

ax.set_xlabel(r'$f/f_0$', fontsize = 12)
ax.set_ylabel(r'Intensity (arb. unit)', fontsize = 12, labelpad = 1)
#ax.set_yscale('log')
ax.set_xlim(0.9, 1.1)
def fmt_two_digits(x, pos):
    return f'{x:.2f}'
ax.yaxis.set_major_formatter(fmt_two_digits)
ax.set_yticks([0, 0.05, 0.1, 0.15, 0.2])
#ax.set_xlim(50, 70)
#ax.set_ylim(-0.05, 1.1)
#ax.grid()

#ax.text(0.925, 0.95, '($b$)', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
fig.savefig('freq-vs-lamds-ipseed'+flag+'_.png')

#%%% Call a cell in a loop, batch processing
for z0 in np.linspace(0.3, 3.6, 12):
    flag = '_%.2fm' % z0
    runcell('Make plots for batch analysis: spectrum in THz', '//afs/ifh.de/group/pitz/data/lixiangk/sim3/2023/THzComm2/post.py')
    #runcell('Make plots for batch analysis: temporal profile', '//afs/ifh.de/group/pitz/data/lixiangk/sim3/2023/THzComm2/post.py')

#%% Compare bunching vs seed from input and from Genesis output

fig, ax1 = plt.subplots(nrows = 1, figsize = (5, 4))

for seed in np.arange(2, 21):
    fname = 'bunching_32768@scan.%d.out.par.h5' % seed
    #fname = 'bunching@scan.%d.out.par' % seed
    stat = np.loadtxt(fname)
    
    ax1.plot(stat[:,0], stat[:,-1], 'r-')
    
    # fname = 'pithz.%d.out.h5' % seed
    # f1 = h5py.File(fname, 'r')
    # bunching = f1.get('bunching')
    # ax1.plot(bunching[0]**2)
    
    fname = '../beam_112A_2.0nC_gaussian/bunching@scan.%d.out.par' % seed
    stat = np.loadtxt(fname)
    
    ax1.plot(stat[:,0], stat[:,-1], 'b-')
    
ax1.set_yscale('log')
ax1.grid()
ax1.set_xlabel(r'Number of slices')
ax1.set_ylabel(r'$|b_f|^2$')
ax1.set_xlim(0, 100)

fig.tight_layout()
#fig.savefig('bunching@'+outputName+'.png')

#%% Compare bunching vs seed from input and from Genesis output

fig, ax = plt.subplots(nrows = 3, figsize = (5, 6), sharex = True)

seed = 5
fname = 'bunching_32768@scan.%d.out.par.h5' % seed
stat = np.loadtxt(fname)

ilast = 100
ax1, ax3, ax4 = ax.flatten()
ax1.plot(stat[:ilast,0], stat[:ilast,3], '-')
ax1.plot(stat[:ilast,0], 1/stat[:ilast,5], 'k--')
ax1.set_yscale('log')
ax1.grid()
#ax1.set_xlabel(r'Number of slices')
ax1.set_ylabel(r'$|b_s|^2$')
ax1.set_ylim(1e-9, 1)
ax1.legend([r'$|b_s|^2$', r'$1/N_{\rm e,slice}$'])

# ax2.plot(stat[:ilast,0], stat[:ilast,-3], '-'); 
# print('Sum of bunching: ', np.sum(stat[:ilast,-3]))
# ax2.set_yscale('log')
# ax2.set_ylabel(r'$N^2|b_f|^2$')

#ax2.plot(stat[:ilast,0], stat[:ilast,8]/stat[:ilast,7], '-')
# ax2.plot(stat[1:ilast,0], stat[1:ilast,6], '-')
# ax2.set_ylabel('Phase (rad)')

ax3.plot([], [], 'r-')
ax3.plot([], [], 'k--')
ax3.grid()
ax3.plot(stat[:ilast,0], stat[:ilast,3]*stat[:ilast,2]**1, '-')
ax32 = ax3.twinx()
ax32.plot(stat[:ilast,0], stat[:ilast,2], 'k--')

ax32.set_ylabel(r'$I_{\rm peak}$ (A)')
ax3.set_ylabel(r'$I\cdot|b_s|^2$')

ax3.set_ylim(1e-9, 1)
ax3.set_yscale('log')

ax3.legend([r'$I|b_s|^2$', r'$I_{\rm peak}$'])

pg = PostGenesis13(fname = '../beam_112A_2.0nC_shotnoise/pithz.5.out.h5', version = 3)
gamma = pg.get_beamdata('energy')
current = pg.current
ax4.plot((gamma[-1]-gamma[0])*current, '-')
ax4.set_ylabel(r'$I\cdot\Delta\gamma_s$')
ax4.grid()

ax4.set_xlim(0, ilast)
ax4.set_xlabel(r'Number of slices')

fig.tight_layout()
fig.savefig(fname+'_.png')

#%% Check the distribution of the slices generated by Astra2GenesisSlices
# Load Genesis slices
f1 = h5py.File('scan.1.out.par.h5', 'r')

slicecount = f1.get('slicecount')[0]
slicelength = f1.get('slicelength')[0]
slicespacing = f1.get('slicespacing')[0]
refposition = f1.get('refposition')[0]

fig, [ax1, ax2, ax3] = plt.subplots(ncols = 3, figsize = (12, 4))
bf_all = 0
current_all = 0
for i in np.arange(86):
    
    s1 = f1.get('slice%06d' % (i+1))
    theta1 = s1.get('theta')[:]
    zz = (theta1/np.pi/2+0.5)*slicelength+refposition
    gamma = s1.get('gamma')[:]
    current = s1.get('current')[0]
    
    bf_all += np.mean(np.exp(-1j*theta1))*current
    current_all += current
    
    range = [[refposition*1e3, 1e3*(refposition+slicespacing)], [32.5, 34]]
    ax1.hist2d(zz*1e3, gamma, bins = (4, 150), range = range, cmin = 1e-9, 
               weights = current*np.ones(zz.shape))
    ax3.hist(zz*1e3, bins = 8, range = range[0], 
             weights = current*np.ones(zz.shape), 
             histtype = r'step')
    refposition += slicespacing
    
ax1.set_xlim(0, refposition*1e3)

# # Load Astra file
# data = pd_loadtxt('../368A.2809.002.1000')
# data[1:,2] -= data[1:,2].min()
# ax2.hist2d(data[1:,2]*1e3, data[1:,5]/510000, bins = (320, 150), cmin = 1e-9)
# ax2.set_xlim(0, refposition*1e3)
# ax2.set_ylim(-1, 0.5)
# ax3.hist(data[1:,2]*1e3, bins = 80*4)

# Load Genesis slices
f1 = h5py.File('scan.1.out.par.h5', 'r')

slicecount = f1.get('slicecount')[0]
slicelength = f1.get('slicelength')[0]
slicespacing = f1.get('slicespacing')[0]
refposition = f1.get('refposition')[0]

fig, [ax1, ax2, ax3] = plt.subplots(ncols = 3, figsize = (12, 4))
bf_all = 0
current_all = 0
for i in np.arange(40, 43):
    
    s1 = f1.get('slice%06d' % (i+1))
    theta1 = s1.get('theta')[:]
    zz = ((theta1/np.pi/2+0.5)*slicelength+refposition)
    #zz /= 100e-3
    gamma = s1.get('gamma')[:]
    current = s1.get('current')[0]
    
    bf_all += np.mean(np.exp(-1j*theta1))*current
    current_all += current
    
    range = [[refposition*1e3, 1e3*(refposition+slicespacing)], [32.5, 34]]
    ax1.hist2d(zz*1e3, gamma, bins = (4, 150), range = range, cmin = 1e-9, 
               weights = current*np.ones(zz.shape))
    ax3.hist(zz*1e3, bins = 8, range = range[0], 
             weights = current*np.ones(zz.shape), 
             histtype = r'step')
    refposition += slicespacing
    
ax1.set_xlim(0, refposition*1e3)

#%% Check the evolution of bunching factors along the undulator
f1 = h5py.File('pithz.1.out.par.h5', 'r')

steps = f1.keys()
steps = [step for step in steps]

bf_ = []
for step in steps:
    if 'step' in step:
        beam = f1.get(step)
        islice = 45
        sname = 'slice%06d' % islice
        s1 = beam.get(sname)
        
        theta1 = s1.get('theta')[:]
        
        current = s1.get('current')[0]
        
        bf = np.mean(np.exp(-1j*theta1))
        bf_.append([bf])
    
amp = np.abs(bf_)

fig, ax = plt.subplots()
ax.plot(amp)
ax.set_yscale('log')
ax.legend(['#15', '#30', '#45'])
ax.set_ylabel('Bunching')
ax.set_xlabel(r'# of steps')

#%%  Load slices and calculate the bunching factors
#f1 = h5py.File('scan.1.out.par.h5', 'r')
f1 = h5py.File('g4.101.0.par.h5', 'r')
res = []; 
for i in np.arange(1, 100):
    sss = str.format('slice%06d/current' % i); print(sss)
    c = f1.get(sss)[0]; 
    
    sss = str.format('slice%06d/theta' % i); print(sss)
    theta = f1.get(sss)[:]
    
    bf = np.abs(np.mean(np.exp(1j*theta)))
    
    theta_c = np.mean(theta)
    
    res.append([i, c, bf, theta_c, bf**2*(c*100e-6/g_c/g_qe)**2])
    
res = np.array(res)
#%
fig, ax = plt.subplots(nrows = 3, figsize = (6, 6), sharex = True)
ax1, ax2, ax3 = ax.flatten()[:]

res = res[:100,:]
ax1.plot(res[:,0], res[:,1], '-')

ax2.plot(res[:,0], res[:,2]**2, 'b-'); ax2.set_yscale('log')
ax21 = ax2.twinx()
ax21.plot(res[:,0], res[:,4]/res[:,4].max(), '-'); ax21.set_yscale('log')

ax3.plot(res[:,0], res[:,3], '-')

ax1.set_ylabel(r'Current (A)')
ax2.set_ylabel(r'$|b_f|^2$')
ax21.set_ylabel(r'$|b_fN_e|^2$')

ax3.set_ylabel(r'$\langle \theta \rangle$ (rad)')
ax3.set_xlabel(r'# of slice')

#ax3.set_ylim([-0.1, 0.1])

ax2.legend([r'$|b_f|^2$'])
ax21.legend([r'$|b_fN_e|^2$'], loc = 'lower right')

fig.tight_layout()
#fig.savefig('Bunching_.png')
#%% Manually compress the bunch with R56
dist = pd_loadtxt('ast.1826.002')

delta = dist[1:,5]/dist[0,5]
z0 = dist[1:,2]

R56 = 0.2
z1 = z0+delta*R56

plt.figure()
plt.hist(z0, bins = 100, histtype = r'step')
plt.hist(z1, bins = 100, histtype = r'step')


#%% Bunching paper gain curve SASE vs seeding

workdir = r'\\afs\ifh.de\group\pitz\data\lixiangk\sim3\2023\THzComm2\G2-quite_loading'
os.chdir(workdir)

fig, ax = plt.subplots(figsize = (4, 4))
ax.plot([], [], 'r--*')
ax.plot([], [], 'b-')
ax.plot([], [], 'g-')
#ax.plot([], [], 'b-')

###
measure = r'\\afs\ifh.de\group\pitz\doocs\measure\THz\2023\PulseEnergy\2023\20230415M'
basename = measure+os.sep+'GainCurve_1429'

fname = basename+'_.mat'
f1 = h5py.File(fname, 'r')

#NoP = f1.get('sNOP')[0,0]
NoP = np.array([1, 1, 1, 1, 1, 1, 2, 3])

zz = f1.get('UeffLenghtArray')[:,0]*1e-3
PE = f1.get('avgPE')[:,0]/NoP
errPE = f1.get('errorPE')[:,0]/NoP
stdPE = errPE/PE*100

#ax.errorbar(zz, PE*1e6, yerr = errPE*1e6, linestyle = '-', marker = '<')
ax.errorbar(zz, PE*1e6*2, yerr = errPE*1e6*2, linestyle = '--', marker = '*')
###

ext = '_.txt'
EE = np.loadtxt('energy-z-ipseed'+ext); EE = EE.T; print (EE.shape)

end = EE.shape[1]

output = f'E =: {np.mean(EE[-1,1:end])*1e6:.2f} +/- {np.std(EE[-1,1:end])*1e6:.2f} uJ\n'
print(output)

for i in np.arange(1, end):
    print(i, end = ' ')
    ax.plot(EE[8:-8,0]-0.1, EE[8:-8,i]*1e6, '-', color='grey') # energy: uJ
ax.plot(EE[8:-8,0]-0.1, np.mean(EE[8:-8,1:], 1)*1e6, 'b-')
#ax2.plot(EE[8:-8,0]-0.1, np.std(EE[8:-8,1:], 1)/np.mean(EE[8:-8,1:], 1)*100, 'b--')


ax.set_xlabel(r'$z$ (m)', fontsize = 12)
ax.set_ylabel(r'Pulse energy ($\mu$J)', fontsize = 12)
ax.set_yscale('log')
ax.set_xlim(0., 3.5)
ax.set_ylim(1e-8, 1000)

ax.legend(['Measured', 'Quiet loading'], fontsize = 11)
fig.savefig('energy-vs-z-ipseed_SASE_vs_1.png')

ext = '.txt'
EE = np.loadtxt('..'+os.sep+'beam_112A_2.0nC'+os.sep+'energy-z-ipseed'+ext); EE = EE.T; print (EE.shape)

end = EE.shape[1]

output = f'E =: {np.mean(EE[-1,1:end])*1e6:.2f} +/- {np.std(EE[-1,1:end])*1e6:.2f} uJ\n'
print(output)

Dz = 0
for i in np.arange(1, end):
    print(i, end = ' ')
    ax.plot(EE[8:-8,0]-0.1+Dz, EE[8:-8,i]*1e6, '-', color='grey') # energy: uJ
ax.plot(EE[8:-8,0]-0.1+Dz, np.mean(EE[8:-8,1:], 1)*1e6, 'g-')

ax.legend(['Measured', 'Quiet loading', 'Local bunching'], fontsize = 11, loc = 'upper left')

fig.savefig('energy-vs-z-ipseed_SASE_vs_2.png')

#%% Genesis1.3 V3 vs V4 comparison
#%%% energy vs z

fname = 'energy-z-ipseed.txt'
files = ['beam_112A_2.0nC_shotnoise_ffspec0'+os.sep+fname,
         'beam_112A_2.0nC_shotnoise_v4'+os.sep+fname]

fig, ax = plt.subplots(figsize = (4, 4))
ax.plot([], [], 'r-')
ax.plot([], [], 'b--')
ax.plot([], [], 'g-')
styles = ['r-', 'b--']

for ic, file in enumerate(files):
    EE = np.loadtxt(file); EE = EE.T; print (EE.shape)
    end = EE.shape[1]
    
    output = f'E =: {np.mean(EE[-1,1:end])*1e6:.2f} +/- {np.std(EE[-1,1:end])*1e6:.2f} uJ\n'
    print(output)
    
    for i in np.arange(1, end):
        print(i, end = ' ')
        ax.plot(EE[8:-8,0]-0.1, EE[8:-8,i]*1e6, '-', color='grey') # energy: uJ
    ax.plot(EE[8:-8,0]-0.1, np.mean(EE[8:-8,1:], 1)*1e6, styles[ic])

ax.set_xlabel(r'$z$ (m)', fontsize = 12)
ax.set_ylabel(r'Pulse energy ($\mu$J)', fontsize = 12)
ax.set_yscale('log')
ax.set_xlim(0., 3.5)
ax.set_ylim(1e-3, 1e3)

ax.legend(['V3', 'V4'], fontsize = 11)
fig.savefig('energy-vs-z-ipseed_V3_vs_V4.png')

#%%% power vs z
fname = 'peak-power-z-ipseed.txt'
files = ['beam_112A_2.0nC_shotnoise_ffspec0'+os.sep+fname,
         'beam_112A_2.0nC_shotnoise_v4'+os.sep+fname]

fig, ax = plt.subplots(figsize = (4, 4))
ax.plot([], [], 'r-')
ax.plot([], [], 'b--')
ax.plot([], [], 'g-')
styles = ['r-', 'b--']

for ic, file in enumerate(files):
    pp = np.loadtxt(file); pp = pp.T; print (pp.shape)
    end = pp.shape[1]
    
    output = f'P =: {np.mean(pp[-1,1:end])/1e6:.2f} +/- {np.std(pp[-1,1:end])/1e6:.2f} MW\n'
    print(output)
    
    for i in np.arange(1, end):
        ax.plot(pp[:,0], pp[:,i]/1e6, '-', color='grey') # energy: MW
    ax.plot(pp[:,0], np.mean(pp[:,1:], 1)/1e6, styles[ic])

ax.set_xlabel(r'$z$ (m)')
ax.set_ylabel(r'peak power (MW)')
ax.set_yscale('log')
ax.grid()

ax.set_xlim(0., 3.5)
ax.set_ylim(1e-3, 1e3)

ax.legend(['V3', 'V4'], fontsize = 11)
fig.savefig('peak-power-vs-z-ipseed_V3_vs_V4.png')

#%%% temporal profiles
fname = 'power-t-ipseed_3.60m.txt'
files = ['beam_112A_2.0nC_shotnoise_ffspec0'+os.sep+fname,
         'beam_112A_2.0nC_shotnoise_v4'+os.sep+fname]

fig, ax = plt.subplots(figsize = (4, 4))
ax.plot([], [], 'r-')
ax.plot([], [], 'b--')
ax.plot([], [], 'g-')
styles = ['r-', 'b--']

for ic, file in enumerate(files):
    pp = np.loadtxt(file); pp = pp.T; print (pp.shape)
    end = pp.shape[1]
    
    jitter = []
    for i in np.arange(1, end):
        ax.plot(pp[:,0]*1e12-0, pp[:,i]/1e6, '-', color='grey') # energy: MW
        jitter.append([weighted_mean(pp[:,0], pp[:,i])*1e12-0])
    ax.plot(pp[:,0]*1e12-0, np.mean(pp[:,1:], 1)/1e6, styles[ic])

    output += f'Arrival time jitter: {np.mean(jitter):.2f} +/- {np.std(jitter):.2f} ps\n'
    print(output)

ax.set_xlabel(r'$t$ (ps)')
ax.set_ylabel(r'intensity (arb. unit)')
#ax.set_yscale('log')
#ax.set_xlim(-20, 20)
#ax.set_ylim(-0.005, 0.6)
ax.grid()

ax.legend(['V3', 'V4'], fontsize = 11)
fig.savefig('power-vs-t-ipseed_V3_vs_V4.png')

#%% Bunching factors vs quiet loading
workdir = r'\\Afs\ifh.de\group\pitz\data\lixiangk\sim3\2023\THzComm2\beam_112A_2.0nC_shotnoise_v4'
os.chdir(workdir)

fname = 'g4.2.out.h5'
version = 4
pg1 = PostGenesis(fname, debug = 0, version = version, fig_ext = '.png')

workdir = r'\\Afs\ifh.de\group\pitz\data\lixiangk\sim3\2023\THzComm2\beam_112A_2.0nC_shotnoise_v4_quiet2'
os.chdir(workdir)

fname = 'g4.24.out.h5'
version = 4
pg2 = PostGenesis(fname, debug = 0, version = version, fig_ext = '.png')

from IPython import display
fig, [ax1, ax2] = plt.subplots(nrows = 2, figsize = (4, 6), sharex = True)
plt.show()

if version == 4:
    name1 = 'intensity-farfield'
    name2 = 'phase-farfield'
    
    ppower1 = pg1.get_fielddata('power')
    p_mid1 = pg1.get_fielddata(name1)
    phi_mid1 = pg1.get_fielddata(name2)
    gamma1 = pg1.get_beamdata('energy')
    bunching1 = pg1.get_beamdata('bunching')
    
    ppower2 = pg2.get_fielddata('power')
    p_mid2 = pg2.get_fielddata(name1)
    phi_mid2 = pg2.get_fielddata(name2)
    gamma2 = pg2.get_beamdata('energy')
    bunching2 = pg2.get_beamdata('bunching')

temp = p_mid1
step = 1
#for i in np.arange(10, temp.shape[0], step):
for i in np.arange(0, temp.shape[0], step):   
    
    ax1.plot(p_mid1[i]/np.max(p_mid1[i])*3, '-')
    ax1.plot(p_mid2[i]/np.max(p_mid2[i])*3, '-')
    
    
    ax1.set_title('z = : %.3f m' % pg.zplot[i])
    # ax1.set_xlabel('# of slice')
    # ax1.set_ylabel('Power (MW)')
    # #ax.set_ylim(0, ymax)
    # ax1.grid()
    ax1.set_xlim(0, 120)
    
    ax2.plot(phi_mid1[i], '-')
    ax2.plot(phi_mid2[i], '-')
    
    # ax2.set_xlabel(r'$z$ (m)')
    # ax2.set_ylabel(r'Peak power (MW)')
    # ax2.grid()
    ax2.set_xlim(0, 120)
    
    display.display(plt.gcf())
    display.clear_output(wait=True)
    
    plt.pause(1)
    if i<temp.shape[0]-step:
        ax1.cla()
        ax2.cla()
        
#%% Compare simulations with local bunching factors and with shot noises

#%%% Energy
flag = '_gaussian'
flag = ''
fname = 'energy-z-ipseed.txt'
files = ['beam_112A_2.0nC_shotnoise_v4'+flag+os.sep+fname,
         'beam_112A_2.0nC_shotnoise_v4'+flag+'_quiet2e'+os.sep+fname]

         #'beam_112A_2.0nC_shotnoise_v4_gaussian_quiet2_'+os.sep+fname,
         
fig, ax = plt.subplots(figsize = (4, 4))
ax.plot([], [], 'r-')
ax.plot([], [], 'b--')
#ax.plot([], [], 'g-')
styles = ['r-', 'b--', 'g:']

for ic, file in enumerate(files):
    EE = np.loadtxt(file); EE = EE.T; print (EE.shape)
    end = EE.shape[1]
    
    output = f'E =: {np.mean(EE[-1,1:end])*1e6:.2f} +/- {np.std(EE[-1,1:end])*1e6:.2f} uJ\n'
    print(output)
    
    for i in np.arange(1, end):
        print(i, end = ' ')
        ax.plot(EE[8:-8,0]-0.1, EE[8:-8,i]*1e6, '-', color='grey') # energy: uJ
    ax.plot(EE[8:-8,0]-0.1, np.mean(EE[8:-8,1:], 1)*1e6, styles[ic])

ax.set_xlabel(r'$z$ (m)', fontsize = 12)
ax.set_ylabel(r'Pulse energy ($\mu$J)', fontsize = 12)
ax.set_yscale('log')
#ax.set_xlim(0., 3.5)
ax.set_ylim(1e-6, 1e3)

#ax.legend(['V3', 'V4'], fontsize = 11)
fig.savefig('energy-vs-z-ipseed'+flag+'_bunching_vs_shot_nosie.png')

#%%% Spectrum
flag = '_gaussian'
flag = ''

fname = 'energy-z-ipseed.txt'
files = ['beam_112A_2.0nC_shotnoise_v4'+flag+os.sep+'spectrum-lamds-ipseed_3.60m.txt',
         'beam_112A_2.0nC_shotnoise_v4'+flag+'_quiet2e'+os.sep+'spectrum-lamds-ipseed_5.40m.txt']

fig, ax = plt.subplots(figsize = (4, 4))
ax.plot([], [], 'r-')
ax.plot([], [], 'b--')
#ax.plot([], [], 'g-')
styles = ['r-', 'b--', 'g:']

sampling = [2, 10]
for ic, file in enumerate(files):
    ss = np.loadtxt(file); ss = ss.T; print (ss.shape)
    end = ss.shape[1]
    
    xlamds = []; width = []
    for i in np.arange(1, end-1):
        if i%sampling[ic] == 0:
            ax.plot(ss[:,0]*1e6, ss[:,i]/np.mean(ss[:,i]), styles[ic], alpha = 0.3) # energy: MW
        
    xlamds.append([weighted_mean(ss[:,0], ss[:,i])*1e6])
    width.append([weighted_std(ss[:,0], ss[:,i])*1e6])
    ax.plot(ss[:,0]*1e6, np.mean(ss[:,1:], 1)/np.mean(np.mean(ss[:,1:], 1)), styles[ic], linewidth = 2)

    output =  f'Centre wavelength: {np.mean(xlamds):.2f} +/- {np.std(xlamds):.2f} um\n'
    output += f'Spectrum width: {np.mean(width):.2f} +/- {np.std(width):.2f} um\n'
    print(output)
    
ax.set_xlabel(r'$\lambda_s$ ($\mu$m)', fontsize = 12)
ax.set_ylabel(r'Intensity (arb. unit)', fontsize = 12, labelpad = 1)
#ax.set_yscale('log')
ax.set_xlim(90, 115)

#ax.grid()

#ax.text(0.925, 0.95, '($b$)', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
fig.savefig('spectrum-vs-lamds-ipseed'+flag+'_bunching_vs_shot_nosie.png')

#%%% Temporal profile
flag = '_gaussian'
flag = ''

files = ['beam_112A_2.0nC_shotnoise_v4'+flag+os.sep+'power-t-ipseed_3.60m.txt',
         'beam_112A_2.0nC_shotnoise_v4'+flag+'_quiet2e'+os.sep+'power-t-ipseed_5.40m.txt']
 

def fmt_two_digits(x, pos):
    return f'{x:.1f}'

fig, ax = plt.subplots(figsize = (4, 4))
ax.plot([], [], 'r-')
ax.plot([], [], 'b--')
#ax.plot([], [], 'g-')
styles = ['r-', 'b--', 'g:']

sampling = [2, 10]
for ic, file in enumerate(files):
    ss = np.loadtxt(file); ss = ss.T; print (ss.shape)
    end = ss.shape[1]
    
    duration = []
    arrival = []
    for i in np.arange(1, end-1):
        if i%sampling[ic] == 0:
            ax.plot(ss[:,0]*1e12, ss[:,i]/np.mean(ss[:,i]), styles[ic], alpha = 0.3) # energy: MW
        arrival.append([weighted_mean(ss[:,0], ss[:,i])*1e12])
        duration.append([weighted_std(ss[:,0], ss[:,i])*1e12])
    ax.plot(ss[:,0]*1e12, np.mean(ss[:,1:], 1)/np.mean(np.mean(ss[:,1:], 1)), styles[ic])

output = f'Arrival time: {np.mean(arrival):.2f} +/- {np.std(arrival):.2f} ps\n'
output += f'Pulse duration: {np.mean(duration):.2f} +/- {np.std(duration):.2f} ps\n'

print(output)

with open('stat.dat', 'a') as f_handle:
    f_handle.write(output)

ax.set_xlabel(r'$t$ (ps)')
ax.set_ylabel(r'Intensity (arb. unit)', labelpad = 1)

ax.yaxis.set_major_formatter(fmt_two_digits)
#ax.set_yscale('log')
ax.set_xlim(0, 60)
#ax.set_ylim(-0.005, 0.3)
#ax.grid()
fig.savefig('power-vs-t-ipseed'+flag+'_bunching_vs_shot_nosie.png')

#%%% Initial bunching amp.
fig, ax1 = plt.subplots(figsize = (4, 4), nrows = 1, sharex = True)
color_cycle = ['r', 'b', 'g', 'c', 'm', 'y', 'k']

cwd = os.getcwd()

os.chdir('beam_112A_2.0nC_shotnoise_v4')

for i in np.arange(1, 6):
    
    seed = i
    fname = f'g4.{seed}.out.h5'
    version = 4
    pg = PostGenesis(fname, debug = 0, version = version, fig_ext = '.png')
    
    data = pg.get_beamdata('bunching')
    iz = int(zi/0.015)
    
    ax1.plot(data[iz,1:86]**2, '-', color = color_cycle[i]) # energy: MW

os.chdir(cwd)
os.chdir('beam_112A_2.0nC_shotnoise_v4_quiet2e')

ax2=ax1
for i in np.arange(1, 6):
    
    seed = i
    fname = f'g4.{seed}.out.h5'
    version = 4
    pg = PostGenesis(fname, debug = 0, version = version, fig_ext = '.png')
    
    data = pg.get_beamdata('bunching')
    iz = int(zi/0.015)
    
    ax2.plot(data[iz,1:86]**2, '--', color = color_cycle[i]) # energy: MW

os.chdir(cwd)

for axis in [ax1]:
    axis.set_ylabel(r'$|B_f|^2$', labelpad = -7)
    axis.yaxis.set_major_formatter(fmt_two_digits)
    axis.set_yscale('log')
    
ax2.set_xlabel(r'# of slices')

#fig.tight_layout()
fig.savefig('Bunching'+flag+'_bunching_vs_shot_nosie.png')

#%%% Initial bunching phase
fig, [ax1, ax2] = plt.subplots(figsize = (4, 4), nrows = 2, sharex = True, sharey = True)
color_cycle = ['r', 'b', 'g', 'c', 'm', 'y', 'k']

cwd = os.getcwd()

os.chdir('beam_112A_2.0nC_shotnoise_v4')

for i in np.arange(1, 6):
    
    seed = i
    fname = f'g4.{seed}.out.h5'
    version = 4
    pg = PostGenesis(fname, debug = 0, version = version, fig_ext = '.png')
    
    data = pg.get_beamdata('bunchingphase')
    iz = int(zi/0.015)
    iz = 0
    
    ax1.plot(data[iz,1:86], '-', color = color_cycle[i]) # energy: MW

os.chdir(cwd)
os.chdir('beam_112A_2.0nC_shotnoise_v4_quiet2e')

for i in np.arange(1, 6):
    
    seed = i
    fname = f'g4.{seed}.out.h5'
    version = 4
    pg = PostGenesis(fname, debug = 0, version = version, fig_ext = '.png')
    
    data = pg.get_beamdata('bunchingphase')
    iz = int(zi/0.015)
    iz = 0
    
    ax2.plot(data[iz,1:86], '--', color = color_cycle[i]) # energy: MW

os.chdir(cwd)

for axis in [ax1, ax2]:
    #axis.set_ylabel(r'Bunch. phase (rad)', labelpad = 0)
    axis.yaxis.set_major_formatter(lambda x, pos:f'{x:.1f}')
    #axis.set_yscale('log')
    
ax2.set_xlabel(r'# of slices')
fig.text(0.04, 0.5, 'Bunching phase (rad)', va='center', ha='center', rotation='vertical', 
         fontsize = 12)

#fig.tight_layout()
fig.savefig('Bunching-phase'+flag+'_bunching_vs_shot_nosie.png')