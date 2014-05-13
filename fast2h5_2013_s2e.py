# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

#$%pylab inline
#Importing necessary modules:
import sys
import os
import errno

#@
sys.path.insert(0,'/data/S2E/packages/WPG')
#$sys.path.insert(0,'..')

import shutil
import numpy
import h5py

#Import base wavefront class
from wpg import Wavefront

# <codecell>

def mkdir_p(path):
    """
    Create directory tree, if not exists (mkdir -p)
    """
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

# <codecell>

# set name of original FEL input file
# i.e. input data file for FORTRAN <fast2xy>.exe, for details see README
# the <fast2xy>.exe takes its input parameters from config file <FAST2XY>.dat
def set_iname_2013(namg,ifb,nz):
    # suppose that undulator position identifier, nz, has 2 digits
    name = namg+'T'+str(ifb)+'0'+str(nz) 
    return name

# <codecell>

# build the names of <fast2xy>.exe output files, namt={'FXY', 'E', 'PXY'}
# for details see README
def set_oname_2013(namg,namt,ifb,nz):
    # suppose that undulator position identifier, nz, has 2 digits
    name = namg+namt+str(ifb)+'0'+str(nz) 
    return name

# <codecell>

# set name of FEL output h5 file, i.e. add a suffix to a base name 
# bname base name   
# tmin is starting point of pulse chunk
# return bname+'_'+str(tmin)+'fs'
def set_FELout_name(bname,tmin):
    strd1=str(int(numpy.floor(tmin)))
    if (trd1+0.5) >(int(trd1) +1.):
        strd1=str(int(numpy.ceil(trd1)))
        fname = bname +'_'+strd1+'fs'
    if int(trd1) < 10:
        fname = bname+'_00'+strd1+'fs'
    if  (int(trd1) >= 10) and (int(trd1) < 100):
        fname = bname+'_0'+strd1+'fs'
    return fname

# <codecell>

def parse_toe_file(f_name):
    """ Parse *.res file to list of strings """
    #TODO: now all rows keeped in memory. We should store only slices we need.
    rows = []
    params = []
    with open(f_name,'r') as f:
        for il,l in enumerate(f.readlines()):    
            try:
                t=numpy.array(l.split(),dtype='float32')
                if il<2:
                    params.append(t)
                else:
                    rows.append(t)
            except ValueError:
                print 'Error pasing line #%s in file %s' %(il,f_name)
    return params, rows


#Example of usage:
#params, rows=parse_toe_file(toe_file_name)

# <codecell>

def create_numpy_array_from_rows(rows,slices=None):
    # slice size (Re, Im)
    N=len(rows[0])/2
    
    if slices is None:
        slices=range(len(rows)/N)
    slice_count=len(slices)

    y = numpy.array(rows,dtype='float32').reshape((slice_count,N,N,2))
    
    print y.shape
    y=numpy.swapaxes(y,0,2)
    print y.shape
    return y

# <codecell>

#save 3d array to hdf5
def store_wavefront_hdf5(wf_struct,file_name):
    """ Store wavefront structure in hdf5 file"""
    def store_group(group,pearent_goup):
        for k,v in group.items():
            if isinstance(v,dict):
                tmp_group=pearent_goup.create_group(k)
                store_group(v,tmp_group)
            else:
                store_value(k,v,pearent_goup)
    
    def store_value(name,value,group):
        filed_data=value[0]
        filed_type=value[1]
        if filed_type=='s':
            group.create_dataset(name,data=filed_data)
        else:
            try:
                group.create_dataset(name,data=filed_data,dtype=filed_type,chunks=True,
                             compression='gzip', compression_opts=9)
            except TypeError:
                group.create_dataset(name,data=filed_data,dtype=filed_type)
                
    with h5py.File(file_name, 'w') as res_file:
        store_group(wf_struct,res_file)

# <codecell>

def load_wavefront_hdf5(file_name):
    """Return dictionary with fields of wavefront"""
    def load_group(group,d):
        for k,v in dict(group).items():
            if isinstance(v,h5py.Group):
                d[k.encode()]={}
                load_group(v,d[k.encode()])
            elif isinstance(v,h5py.Dataset):
                d[k.encode()]=v.value
            else:
                raise TypeError
    
    wf={}
    with h5py.File(file_name, 'r') as h5_file:
        load_group(h5_file,wf)
        
    return wf

# <codecell>

# Working directory for FEL data convertion:
#@
thepath = '/data/S2E/data/FELsource/' 
#$thepath = '/diskmnt/a/exflwgs03/lsamoylv/code/WPG-develop/samples/' 
os.chdir(thepath)
#Path to original FEL data:  
fel_data_path = thepath
# input parameters for updating FAST2XY_2013.dat
trd1 = 1.9                    # Start time for conversion (fs)
trd2 = 10.9                 # End time for conversion (fs)
nxy = 100                     # Mesh parameter, number of nodes = 2*NXY=+1
nskip = 8                     # skip nskip slices	
ifb = 1001                    # File No. (FAST FEL code run number)
#@
nzc = 33                      # Output point No., defines active undulator length (see below)
#$nzc = 15                      # Output point No., defines active undulator length (see below)
nharm = 1                     #
namg = "SASE1_5keV_14GeV_"	# Generic file name prefix
namt = 'FXY'+str(nharm)+'_'   # part of FEL output ASCII file name 

fast2xydat='FAST2XY_2013.DAT'
fast2xyexe='fast2xy_2013_wo_name.exe'
tmp_dir = os.path.join(thepath,set_FELout_name(set_oname_2013(namg,namt,ifb,nzc),trd1))
mkdir_p(tmp_dir)
print 'All temporary files are saved in \n'+ tmp_dir+'/'
shutil.copy(os.path.join(thepath,fast2xydat), tmp_dir+'/')
shutil.copy(os.path.join(thepath,fast2xyexe), tmp_dir+'/')

# <codecell>

os.getcwd()

# <codecell>

os.chdir(tmp_dir)
f_in = open(fast2xydat,'r')
a = f_in.readlines()
strInputPar = numpy.empty(len(a), dtype=object)
strComment  = numpy.empty(len(a), dtype=object)
print '==Beforehand=='
for idx in range(len(a)): 
    strInputPar[idx] =  a[idx].split(' ',1)[0]
    strComment[idx]  =  a[idx].split(' ',1)[1]
    print strInputPar[idx]+strComment[idx].rstrip()
f_in.close()
strInputPar[0] = str(trd1)
strInputPar[1] = str(trd2)
strInputPar[2] = str(nxy)
strInputPar[3] = str(nskip)
strInputPar[4] = str(ifb)
strInputPar[5] = str(nzc)
strInputPar[6] = namg
print '==Afterwards=='
for idx in range(len(a)): 
    a[idx] = strInputPar[idx]+' '+strComment[idx].rstrip()
    print a[idx]
numpy.savetxt( fast2xydat,a,fmt='%s')

# <codecell>

os.getcwd()

# <codecell>

#def set_iname_2013(namg,ifb,nz):
ifname=set_iname_2013(namg,ifb,nzc)+'.RES'
fulltname=fel_data_path+ifname
if os.system("ln -s "+fulltname+" fort.4") == 0:
    print 'Creating symbolic link (fort.4): OK'
#shutil.move('fort.4', tmp_dir+'/')
outp=os.popen(fast2xyexe).read()

fxy_data_file=set_oname_2013(namg,namt,ifb,nzc)+'.RES'
p_data_file = set_oname_2013(namg,'PXY'+str(nharm)+'_',ifb,nzc)+'.RES'
e_data_file = set_oname_2013(namg,'E'+str(nharm),ifb,nzc)+'.RES'
pz_data_file = 'PZ1001.RES'
print fulltname
print fast2xyexe

# <codecell>

print (outp)
print fulltname
print tmp_dir

# <codecell>

#loading data from FXY*.RES file
print 'Wave field: \tReading text data from %s ...' %fxy_data_file
params, rows=parse_toe_file(os.path.join(tmp_dir,fxy_data_file))
print '\t...done'

print 'Time structure: \tReading text data from %s ...' %e_data_file
e_data = numpy.loadtxt(os.path.join(tmp_dir,e_data_file))
print '\t...done'

print 'Gain curve: \tReading text data from %s ...' %pz_data_file
pz_data = numpy.loadtxt(os.path.join(fel_data_path,pz_data_file))
print '\t...done'

print 'P( x, y ): \tReading text data from %s ...' %p_data_file
prows = numpy.loadtxt(os.path.join(tmp_dir,p_data_file))
N = len(prows)
p_data = numpy.array(prows,dtype='float32').reshape((N,N))
print p_data.shape
del prows
print '\t...done'

# <codecell>

wl = params[0][0]*1e-2
photonEnergy = 12.4e3/(wl*1e10)
slMin0 = 0.
slStep = params[0][1]
xStep = params[0][2]*1e-2
nx = int(params[0][3])
nStart = 1
nEnd = numpy.size(rows)/(2*nx*nx)
#ratio = 6.4550128060534480 
xMin = -xStep*(nx-1)/2
xMax =  xStep*(nx-1)/2
slMin = slMin0 + (nStart-1)*slStep
slMax = slMin0 + (nEnd - 1)*slStep
print 'wl [nm], Eph [keV]',wl*1e9,photonEnergy*1e-3
print 'nStart,nEnd: ',nStart,nEnd
print 'slMin,slMax [fs]: ',slMin*1e15,slMax*1e15
print 'nx,xMin,xMax [um]',nx,xMin*1e6,xMax*1e6

# <codecell>

RK1 = params[0][4]
print RK1

# <codecell>

#build numpy arrays from list of rows
wf_data=create_numpy_array_from_rows(rows,slices=range(nStart-1,nEnd))
del rows

# <codecell>

#wavefront structure based on glossary
wf_struct={'version':(0.1,'f')}
wf_struct['params']={
    'photonEnergy':(photonEnergy,'f'),
    'wDomain':('time','s'),
    'wSpace':('R-space','s'),
    'wEFieldUnit':('sqrt(W/mm^2)','s'),
    'Rx':(2.,'f'),
    'Ry':(2.,'f'),
    'dRx':(0.125,'f'),
    'dRy':(0.125,'f'),
    'xCentre':(0,'f'),
    'yCentre':(0,'f'),
    'nval':(2,'i')
    }
  
wf_struct['params']['Mesh']={
    'nx':(wf_data.shape[1],'i'),
    'ny':(wf_data.shape[0],'i'),
    'xMin':(xMin,'f'),
    'xMax':(xMax,'f'),
    'yMin':(xMin,'f'),
    'yMax':(xMax,'f'),    
    'nSlices':(wf_data.shape[2],'i'),
    'sliceMin':(slMin,'f'),
    'sliceMax':(slMax,'f'),
    'zCoord':(0.0,'d')
   }

wf_struct['data']={
    'arrEver':(wf_data*numpy.sqrt(RK1),'f'),
    'arrEhor':(numpy.zeros(shape=wf_data.shape,dtype='float32'),'f')
    }
wf_struct['misc']={
    }

# <codecell>

#wf_struct_fel={'version':(0.1,'f')}
wf_struct['history/parent/detail/info']={
    'contact':(r'''Mikhail Yurkov
Evgeny Schneidmiller

Deutsches Elektronen Synchrotron (DESY)
Notkestrasse 85
22607 Hamburg
Germany

E-mail: mikhail.yurkov@desy.de
E-mail: evgeny.schneidmiller@desy.de
Tel. +49 40 8998 2676''','s'),
   'data_description':(r'''SASE12, 5 keV, 14 GeV. 
Current pulse lengthis 2.5 ps, or about 10^4 coherence time (0.25 fs) - 
should be sufficient for any statistical analysis. 
Test files are as usually temporal profile, intensity distributions in the near and far zone.
Output data points for s2e simulations:
(see gain_curve):
No.    z(m)
15     35.5 (beginning of saturation regime)
18     42.7 (saturation regime)
25     59.5 (deep nonlinear)
33     78.7  -/-
42    100.3  -/- 
====
Generic file name for simulation run is
SASE1_5keV_14GeV_

Parameters of the run:
Frequency harmonics: w1, w3
Azimuthal harmonics: -2...2
File names:
<generic_filename>??XXXXZZZ
XXXX    -    No. of statistical run
ZZZ    -    No. of output point along undulator 
??     stands for:
P1/P3    -    temporal structure (fs,W) for the 1st/3rd harmonic
S1/S3    -    temporal structure (dw/w [%],spectral power) for the 1st/3rd harmonic
PZ1/PZ3    -    energy in the radiation pulse for the 1st/3rd harmonic harmonic
column 2: z(cm)
column3: total power W
column 4: power in relevant frequency harmonic
columns 5-9: power in azimuthal modes -2...2
N1/N3    -    Intensity distribution in the near zone (cm,intensity)
F1/F3    -    Intensity distribution in the near zone (rad,intensity)
PZ n files:
1 :  No. poit of output
2 : z-coordinate of output (cm)
3 : total radiation power (all freq, harm+all azimuthal modes)
4 : control information (for internal purpose)
5 :  power in n-th frequency harm.
6-10: partial contribution of azimuthal modes -2:2 to n-th freq. harm.
Pn files:
1 : time [fs]
2: radiation power (all freq, harm+all azimuthal modes) [W]
3: control information (for internal purpose)
4:  power in n-th frequency harm.
5-9: partial contribution of azimuthal modes -2:2 to n-th freq. harm.
SFn
1: dw/w [%]
2: spectral power oh n-th frequency harmonic (far zone, zero angle)
Spectrum is normalized to 1: \int P(w)dw = 1 
all other columns as well as SNn files are internal control information.''','s'),
   'method_description':('Code: FAST (Schneidmiller Yurkov)','s'),
   'source_data_path': ('/data1/yurkov/2013/SASE12_5keV_14GeV_long_pulse/','s')}
wf_struct['history/parent/detail/misc']={
    'FAST2XY_DAT':(a,'s'),                                     
    'temporal_struct':(e_data,'f'),
    'gain_curve':(pz_data,'f'),
    'power_xy':(p_data,'f')}

# <codecell>

#store wavefront in hdf5 file
fname0 = set_FELout_name(set_oname_2013(namg,namt,ifb,nzc),trd1)
#print "Store wavefront in hdf5 file: "+fname0+'.h5'
store_wavefront_hdf5(wf_struct,os.path.join(tmp_dir,fname0+'.h5'))

# <codecell>

# use srwlib glossary to add attributes to wavefront datasets 
# ATTENTION: external groups are dissapeared !
in_fname  = os.path.join(tmp_dir,fname0+'.h5')
bare_fname = os.path.join(tmp_dir,fname0+'_bare.h5')
print('Loading wavefront data from the file:     '+in_fname)
wf_struct=Wavefront()
wf_struct.load_hdf5(in_fname)
wfr = wf_struct._srwl_wf
wf_struct = Wavefront(wfr)
print('Saving the wavefront data with attributes:'+bare_fname)
wf_struct.store_hdf5(bare_fname)
print('Replacing data with attributes:  '+bare_fname)
with h5py.File(bare_fname) as h2:
    with h5py.File(in_fname) as h1:
        try:
            del h1['params']  # delete group
        except KeyError:
            pass
        h2.copy('params',h1) #copy h2['params'] to h1

# <codecell>

out_fname = os.path.join(thepath,fname0+'_out.h5')
shutil.copy(in_fname,out_fname)
print 'out_fname:',out_fname
os.chdir(thepath)

# <codecell>

os.getcwd()

# <codecell>

shutil.copy(out_fname, os.path.join('/data/S2E/data/','prop_in.h5'))

