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

def set_pzname_2014(namg,ifb,nz):
    # suppose that undulator position identifier, nz, has 2 digits
    name = namg+'PZ1_000'+str(ifb)+'000'
    return name

# <codecell>

# set name of FEL input file
# i.e. output of pproc-fast2xy-2013-v2-06.exe, for details see
# its input (parameter) file PPROC_FAST2XY_2013.dat
def set_iname_2014(namg,ifb,nz):
    # suppose that undulator position identifier, nz, has 2 digits
    name = namg+'T'+str(ifb)+'0'+str(nz) 
    return name

# <codecell>

# set name of FEL output ASCII file
# i.e. output of pproc-fast2xy-2013-v2-06.exe, for details see
# its input (parameter) file PPROC_FAST2XY_2013.dat
def set_oname_2014(namg,namt,ifb,nz):
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
    
    print 'y=rows.reshape((nSlices,nx,nx,2)',y.shape
    y=numpy.swapaxes(y,0,2)
    print 'return swapaxes(y,0,2)',y.shape
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

def _resample(wf, axis, data, x0, x1):
    if axis.lower()=='x':
        y = data[data.shape[0]/2,:]
        x = numpy.linspace(wf.params.Mesh.xMin, wf.params.Mesh.xMax, y.shape[0])
    elif axis.lower()=='y':
        y = data[:,data.shape[1]/2]
        x = numpy.linspace(wf.params.Mesh.yMin, wf.params.Mesh.yMax, y.shape[0])
    else:
        raise ValueError(
            'Wrong axis {}, should be "x" or "y"'.format(axis))
    
    if not x0 is None:
        xmin = x0
    else:
        xmin = x[0]
    
    if not x1 is None:
        xmax = x1
    else:
        xmax = x[-1] 
    
    x1 = numpy.linspace(xmin,xmax,len(y))
    y1 = numpy.interp(x1, x,y)
    return x1, y1
    
def intensity_cut(wf, axis, polarization, x0=None, x1=None, M=0):
    
    if polarization.lower()  == 'v' or polarization.lower() == 'vertical':
        pol = 'vertical'
    elif polarization.lower() == 'h' or polarization.lower() == 'horizontal':
        pol = 'horizontal'
    elif polarization.lower() == 't' or polarization.lower() == 'total':
        pol = 'total'
    else:
        raise ValueError(
            'Wrong polarization {}, should be "v" or "vertical"'+
            ' or "h" or "horizontal" or "t" or "total"'.format(polarization))
    
    dx=(wf.params.Mesh.xMax - wf.params.Mesh.xMin)/(wf.params.Mesh.nx-1)
    dy=(wf.params.Mesh.yMax - wf.params.Mesh.yMin)/(wf.params.Mesh.ny-1)
    dt=(wf.params.Mesh.sliceMax - wf.params.Mesh.sliceMin)/(wf.params.Mesh.nSlices-1)

    data = wf.get_intensity(slice_number=M, polarization=pol)
    return _resample(wf, axis, data*dx*dy*dt*1e-6, x0, x1)

def phase_cut(wf, axis, polarization, x0=None, x1=None, M=0):
    
    if polarization.lower()  == 'v' or polarization.lower() == 'vertical':
        pol = 'vertical'
    elif polarization.lower() == 'h' or polarization.lower() == 'horizontal':
        pol = 'horizontal'
    else:
        raise ValueError(
            'Wrong polarization {}, should be "v" or "vertical" or "h" or "horizontal"'.format(polarization))
    
    data = wf.get_phase(slice_number=M, polarization=pol)
    return _resample(wf, axis, data, x0, x1)

# <codecell>

def plot_int_xy(wf,xscale=1.,yscale=1.,xtitle='time [s]',ytitle='[a.u.]',capture=' '):
    dx=(wf.params.Mesh.xMax - wf.params.Mesh.xMin)/(wf.params.Mesh.nx-1)
    dy=(wf.params.Mesh.yMax - wf.params.Mesh.yMin)/(wf.params.Mesh.ny-1)
    xaxis=linspace(wf.params.Mesh.sliceMin,wf.params.Mesh.sliceMax,wf.params.Mesh.nSlices)
    # plot the pulse time structure /spectrum 
    pylab.figure()    
    pylab.plot(xaxis*xscale,wf.get_intensity().sum(axis=0).sum(axis=0)*dx*dy*1e-6*yscale,'o-')
    pylab.xlabel(xtitle)
    pylab.ylabel(ytitle)
    pylab.title(capture)
    pylab.grid(True)
    
def plot_pulse_intensity(wf,subtitle):
    dx=(wf.params.Mesh.xMax - wf.params.Mesh.xMin)/(wf.params.Mesh.nx-1)
    dy=(wf.params.Mesh.yMax - wf.params.Mesh.yMin)/(wf.params.Mesh.ny-1)
    dt=(wf.params.Mesh.sliceMax-wf.params.Mesh.sliceMin)/(wf.params.Mesh.nSlices-1)
    pylab.figure()
    pylab.imshow(wf.get_intensity().sum(axis=-1)*dt*dx*dy*1e-6*1e6,extent=wf.get_limits())
    pylab.colorbar()
    pylab.title('pulse intensity [$\mu J$] '+subtitle)

# <codecell>

def update(in_fast2xydat='PPROC-FAST2XY_2013_LP.DAT',trd1=0.,trd2=None,
           nxy=None,nskip=None,ifb=None,nzc=None,namg=None):
    f_in = open(in_fast2xydat,'r')
    a = f_in.readlines()    
    strInputPar = numpy.empty(len(a), dtype=object)
    strComment  = numpy.empty(len(a), dtype=object)
    print '==Beforehand=='
    for idx in range(0,7): 
        strInputPar[idx] =  a[idx].split('#',1)[0]
        strComment [idx]  = a[idx].split('#',1)[1]
        print strInputPar[idx].rstrip()+strComment[idx].rstrip()
    f_in.close()
    strInputPar[0] = str(trd1)
    if trd2 is not None:
        strInputPar[1] = str(trd2)
    if nxy <> None:
        strInputPar[2] = str(nxy)
    if nskip <> None:
        strInputPar[3] = str(nskip)

    if ifb <> None:
        strInputPar[4] = str(ifb)
    else:
        ifb = int(float(strInputPar[4]))

    if nzc <> None:
        strInputPar[5] = str(nzc)
    else:
        nzc = int(float(strInputPar[5]))
    if namg <> None:
        strInputPar[6] = namg
    else:
        namg = strInputPar[6]
    print '==Afterwards=='
    for idx in range(0,7): 
        a[idx] = strInputPar[idx]+' #'+strComment[idx].rstrip()
#    for idx in range(0,len(a)):
    for idx in range(0,7):
        print a[idx]
    aout=a[0:31]
    numpy.savetxt( in_fast2xydat,aout,fmt='%s')  
    return namg,ifb,nzc

# <codecell>

def fill_wf_history_detail(wf_struct):
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
                                             
                                             File naming: ABCD.RES
                                             A == Name of group
                                             B == File Type
                                             C == Sequential number of statistical run
                                             D == Sequential number of output point along undulator
                                             h == Frequency harmonic, h = 1, 3, 5, ...
                                             
                                             File Type B:
                                             T == raw data from FAST
                                             FXYh == raw data from FAST, cartesian coordinate
                                             Ph == Temporal structure of the radiation pulse for frequency harmonic
                                             PZh == Energy in the radiation pulse versus undulator length
                                             SFh == Power spectrum of the radiation pulse
                                             Fh == Intensity distribution in the far zone
                                             Nh == Intensity distribution in the near zone
                                             FWHMh == Evolution along the undulator of the spot size and 
                                             angular divirgence of the radiation
                                             ***
                                             Structure of data files
                                             T == raw data from FAST
                                             
                                             FXYh == raw data from FAST, cartesian coordinate
                                             NXY x NXY x MZ complex array produced by routine fast2xy from FAST file T
                                             
                                             Ph == Temporal structure of the radiation pulse for frequency harmonic h
                                             col[0] Position along bunch  [cm]
                                             col[1] Total radiation power [W]
                                             col[2] Power loss by electrons [W]
                                             col[3] Radiation power in harmonic h [W]
                                             -n
                                             ...  radiation power for frequency harmonic h and azimuthal mode m [W], m = -n .. n
                                             n
                                             
                                             PZh == Energy in the radiation pulse versus undulator length (gain curve)
                                             col[0] Sequential number of output point along undulator
                                             col[1] Position along undulator [cm]
                                             col[2] Total energy in the radiation pulse for frequency harmonic h [J]
                                             col[3] Energy loss by electrons [J]
                                             col[4:4+2*n] 
                                             -n
                                             ...   Energy in the radiation pulse for frequency harmonic h and 
                                             ...   azimuthal harmonic m [J], m = -n .. n
                                             n
                                             
                                             SFh == Power spectrum of the radiation pulse in the far zone at zero angle
                                             col[0] frequency dw/w [%]
                                             col[1] Spectral power [normalized]
                                             
                                             Fh == Intensity distribution in the far zone
                                             col[0] Angle [rad]
                                             col[1] Radiation intensity
                                             
                                             Nh == Intensity distribution in the near zone
                                             col[0] Radial coordinate [cm]
                                             col[1] Radiation intensity
                                             
                                             FWHMh == Evolution along the undulator of the spot size and 
                                             col[0] angular divirgence of the radiation for frequency harmonic h
                                             col[1] position along the undulator [cm]
                                             col[2] FWHM spot size of the radiation [cm]
                                             col[3] FWHM angular divergence of the radiation [rad]
                                             col[4] Radiation pulse energy [J]
                                             
                                             all other columns as well as SNn files are internal control information.''','s'),
                                             'method_description':('Code: FAST (Schneidmiller Yurkov)','s'),
                                             'source_data_path': ('http://dcache-door-photon03:2980/2013-EXFEL-S1-5keV-14GeV-LongPulse/','s')}

# <codecell>

def fill_wf_params(wf_struct,params,nrows):
    wl = params[0][0]*1e-2
    photonEnergy = 12.4e3/(wl*1e10)
    slMin0 = 0.
    slStep = params[0][1]
    xStep = params[0][2]*1e-2
    nx = int(params[0][3])
    nStart=1
    nEnd=nrows/(2*nx*nx)
    xMin = -xStep*(nx-1)/2
    xMax =  xStep*(nx-1)/2
    slMin = slMin0 + (nStart-1)*slStep
    slMax = slMin0 + (nEnd - 1)*slStep
    print 'wl [nm], Eph [keV]',wl*1e9,photonEnergy*1e-3
    print 'nStart,nEnd: ',nStart,nEnd
    print 'nSlices,slMin,slMax [fs]: ',nEnd,slMin*1e15,slMax*1e15
    print 'nx,xMin,xMax [um]',nx,xMin*1e6,xMax*1e6  
    RK1 = params[0][4]
    wf_struct['params']={
                         'photonEnergy':(photonEnergy,'f'),
                         'wDomain':('time','s'),
                         'wSpace':('R-space','s'),
                         'wEFieldUnit':('sqrt(W/mm^2)','s'),
                         'Rx':(5.,'f'),
                         'Ry':(5.,'f'),
                         'dRx':(0.125,'f'),
                         'dRy':(0.125,'f'),
                         'xCentre':(0,'f'),
                         'yCentre':(0,'f'),
                         'nval':(2,'i')
                         }
  
    wf_struct['params']['Mesh']={
                                 'nx':(nx,'i'),
                                 'ny':(nx,'i'),
                                 'xMin':(xMin,'f'),
                                 'xMax':(xMax,'f'),
                                 'yMin':(xMin,'f'),
                                 'yMax':(xMax,'f'),
                                 'nSlices':(nEnd,'i'),
                                 'sliceMin':(slMin,'f'),
                                 'sliceMax':(slMax,'f'),
                                 'zCoord':(0.0,'d')
                                 }
    #xStep*1e-3: |E|^2 in W/mm^2 but not W/m^2
    return RK1/(xStep*1e-3)**2,nStart,nEnd

# <codecell>

def add_wf_attributes(fname0):
    # use srwlib glossary to add attributes to wavefront datasets 
    in_fname   = fname0+'.h5'
    bare_fname = fname0+'_bare.h5'
    print('Loading wavefront data from the file:     '+in_fname)
    wf_struct=Wavefront()
    wf_struct.load_hdf5(in_fname)
    wfr = wf_struct._srwl_wf
    wf_struct = Wavefront(wfr)
    print('Saving the wavefront data with attributes:'+bare_fname)
    wf_struct.store_hdf5(bare_fname)
    print('Replacing data with attributes from  '+bare_fname)
    with h5py.File(bare_fname) as h2:
        with h5py.File(in_fname) as h1:
            try:
                del h1['params']  # delete group
            except KeyError:
                pass
            h2.copy('params',h1) #copy h2['params'] to h1

# <codecell>

def convert_fast2h5(fel_data_path,fast2xyexe,fast2xydat):
    ifname=set_iname_2014(namg,ifb,nzc)+'.RES'
    fulltname=fel_data_path+ifname
    print 'fulltname:',fulltname
    print 'Create sym link:',os.system("ln -s "+fulltname+" fort.4")
    outp=os.popen(fast2xyexe).read()
    print 'Remove sym link:',os.system("rm fort.4")
    print 'fast2xy.exe:\n',outp

    #loading data from text toe*.res file
    nharm=1 # for fundamental
    namt = 'FXY'+str(nharm)+'_'   # part of FEL output ASCII file name 
    print 'Reading text data from: '+ set_oname_2014(namg,namt,ifb,nzc)+'.RES ...'
    params, rows=parse_toe_file(set_oname_2014(namg,namt,ifb,nzc)+'.RES')
    print '...done'

    #wavefront structure based on glossary
    wf_struct={'version':(0.1,'f')}
    RK1,nStart,nEnd = fill_wf_params(wf_struct, params,size(rows))
    #build numpy arrays from list of rows
    wf_data=create_numpy_array_from_rows(rows,slices=range(nStart-1,nEnd))
    del rows
    wf_struct['data']={
                       'arrEver':(wf_data*numpy.sqrt(RK1),'f'),
                       'arrEhor':(numpy.zeros(shape=wf_data.shape,dtype='float32'),'f')
                       }
    wf_struct['misc']={
                       }
    fill_wf_history_detail(wf_struct)
    f_in = open(fast2xydat,'r')
    wf_struct['history/parent/detail/misc']={
                                             'FAST2XY_DAT':(f_in.readlines(),'s')
                                             }
    # < read complimentary FEL data >
    #wf_struct['history/parent/detail/misc']={
    #                                         'temporal_struct':(e_data,'f'),
    #                                         'gain_curve':(pz_data,'f'),
    #                                         'fwhm_curve':(fwhm_data,'f'),
    #                                         'power_xy':(p_data,'f')
    #                                         }
    
    #store wavefront in hdf5 file
    fname0 = set_oname_2014(namg,namt,ifb,nzc)
    print "Store wavefront in hdf5 file: "+fname0+'.h5'
    store_wavefront_hdf5(wf_struct,fname0+'.h5')
    add_wf_attributes(fname0)
    return fname0

# <codecell>

def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-i", "--input-parameter-file",dest="in_fast2xydat", help="Input fast2xy parameter file", )
    parser.add_option("-t", "--time-start",          dest="trd1", help="Start time value for reading the pulse, fs")
    parser.add_option("-q", "--time-end",            dest="trd2", help="End time value for reading the pulse, fs")
    parser.add_option("-n", "--nxy",                 dest="nxy",  help="Number of xy nodes")
    parser.add_option("-s", "--skip-nslices",        dest="nskip",help="skip slices (1 - no skip)")
    parser.add_option("-r", "--run-number",          dest="ifb",  help="XXXX Run number (Input file name TXXXXYYY.RES)")
    parser.add_option("-z", "--zc-point-num",        dest="nzc",  help="YYY	Number of output point along undulator, zc")
    parser.add_option("-p", "--prefix",              dest="namg", help="Prefix of FEL data files")

    (options, args) = parser.parse_args()
    
    if not options.in_fast2xydat:   # if filename is not given
        parser.error('Input parameter filename not given')
    else:
        in_fast2xydat=options.in_fast2xydat
   
    if not options.trd1:   # if time value not given
        parser.error('Start time value not given')
    else:
        trd1=float(options.trd1)

    if not options.nzc:   # if nzc value not given
        nzc=33 #default value
    else:
        nzc=int(float(options.nzc))


    #copy-paste from else block of the cell below
    out_fast2xydat='PPROC-FAST2XY_2013.DAT'
    #$thepath = '/diskmnt/a/exflwgs03/lsamoylv/code/WPG-develop/samples/' 
    #@
    thepath = '/data/S2E/modules/FELsource/' 
    fel_data_path='/pnfs/desy.de/exfel/disk/XFEL/2013/SIM/FAST/2013-EXFEL-S1-5keV-14GeV-LongPulse/'
    os.chdir(thepath)

    namg,ifb,nzc = update(in_fast2xydat,trd1,trd2=trd1+9,nzc=nzc)
    #print namg,ifb,nz

    fast2xyexe='pproc-fast2xy-2013-v2-06-wo-fname.exe'
    work_dir = os.path.join(thepath,namg+str(ifb)+'0'+str(nzc))
    tmp_dir = work_dir+'/tmp'

    mkdir_p(work_dir)
    print 'The result hdf5 file will be saved in \n'+ work_dir+'/'
    mkdir_p(tmp_dir)
    print 'All temporary files will be saved in \n'+ tmp_dir+'/'
    shutil.copy(os.path.join(thepath,in_fast2xydat), tmp_dir+'/'+out_fast2xydat)
    shutil.copy(os.path.join(thepath,fast2xyexe), tmp_dir+'/')

    os.chdir(tmp_dir)
    in_fname=convert_fast2h5(fel_data_path,fast2xyexe,out_fast2xydat)
    out_fname=set_FELout_name(in_fname,trd1)
    print 'in_fname,out_fname:',in_fname,out_fname
    os.system('chmod a+rw '+ tmp_dir+'/*.*')    
    
    shutil.copy(in_fname+'.h5', os.path.join(work_dir,out_fname+'.h5'))
    os.system('rm '+tmp_dir+'/*.*')
    os.chdir(thepath)

# <codecell>

if __name__ == '__main__':
        main()
else:
    # typical command line parameters:
    # fast2xy_new.py -i'PPROC-FAST2XY_2013_LP.DAT' --time-start=3. --skip-nslices=8 --zc-point-num=33
    in_fast2xydat='PPROC-FAST2XY_2013_LP.DAT';trd1=12.;nskip=8;nzc=33
    
    out_fast2xydat='PPROC-FAST2XY_2013.DAT'
    #$thepath = '/diskmnt/a/exflwgs03/lsamoylv/code/WPG-develop/samples/' 
    #@
    thepath = '/data/S2E/modules/FELsource/' 
    fel_data_path='/pnfs/desy.de/exfel/disk/XFEL/2013/SIM/FAST/2013-EXFEL-S1-5keV-14GeV-LongPulse/'
    os.chdir(thepath)

    namg,ifb,nzc = update(in_fast2xydat,trd1,trd2=trd1+9,nzc=nzc,namg='SASE1_5keV_14GeV_LP_')
    #print namg,ifb,nz

    fast2xyexe='pproc-fast2xy-2013-v2-06-wo-fname.exe'
    work_dir = os.path.join(thepath,namg+str(ifb)+'0'+str(nzc))
    tmp_dir = work_dir+'/tmp'

    mkdir_p(work_dir)
    print 'The result hdf5 file will be saved in \n'+ work_dir+'/'
    mkdir_p(tmp_dir)
    print 'All temporary files will be saved in \n'+ tmp_dir+'/'
    shutil.copy(os.path.join(thepath,in_fast2xydat), tmp_dir+'/'+out_fast2xydat)
    shutil.copy(os.path.join(thepath,fast2xyexe), tmp_dir+'/')

    os.chdir(tmp_dir)
    in_fname=convert_fast2h5(fel_data_path,fast2xyexe,out_fast2xydat)
    out_fname=set_FELout_name(in_fname,trd1)
    print 'in_fname,out_fname:',in_fname,out_fname
    os.system('chmod a+rw '+ tmp_dir+'/*.*')    
    
    shutil.copy(in_fname+'.h5', os.path.join(work_dir,out_fname+'.h5'))
    os.system('rm '+tmp_dir+'/*.*')
    os.chdir(thepath)

# <codecell>


# <codecell>


