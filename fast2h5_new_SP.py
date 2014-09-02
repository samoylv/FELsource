
# coding: utf-8

# In[ ]:

# Contact L.Samoylova <liubov.samoylova@xfel.eu>, A.Buzmakov <buzmakov@gmail.com>
# SPB S2E simulation project, European XFEL Hamburg <www.xfel.eu>
# May 2014
# Wave optics software is based on 
# SRW core library <https://github.com/ochubar/SRW>, and 
# WPG framework <https://github.com/samoylv/WPG>


# In[ ]:

isS2E = True        # True if working at S2E server, make it True before downloading as .py!
isIpynb = False       # True if working with iPython notebook
if isS2E:
    isIpynb = False  # at S2E we only run the python script 
     
isLP = False
NHARM = 1


# In[ ]:

#Importing necessary modules
import sys
import os
import errno

if isS2E:
    sys.path.insert(0,'/data/S2E/packages/WPG')
    doPrint = False
else:
    #sys.path.insert(0,'../..')
    sys.path.insert(0,'/diskmnt/a/lsamoylv/WPG')
    doPrint = True

import shutil
import uuid
import numpy
import h5py

#Import base wavefront class
from wpg import Wavefront


# In[ ]:

def dequote(s):
    """
    If a string has single or double quotes around it, remove them.
    Make sure the pair of quotes match.
    If a matching pair of quotes is not found, return the string unchanged.
    """
    if (s[0] == s[-1]) and s.startswith(("'", '"')):
        return s[1:-1]
    return s


# In[ ]:

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


# In[ ]:

def set_ABCDname(prefix,ftype,idx,nz):
    """
    Set the ABCD file name 
    of prefix, file type, run number, and number of z output point
    """
    if ftype[0]=='F' or ftype[0]=='N' or ftype[1:2]=='XY':
        name="%s%s_%s%s.RES" % (prefix,ftype,str(idx).zfill(4),str(nz).zfill(3))
    else:
        name="%s%s%s%s.RES" %(prefix,ftype,str(idx).zfill(4),str(nz).zfill(3))
    return name


# In[ ]:

def set_out_name(prefix,idx):
    """
    Set the output h5 file name 
    of prefix and 7-digit suffix
    """
    name="%s_%s.h5" %(prefix,str(idx).zfill(7))
    return name


# In[ ]:

#print set_ABCDname('S1_0.25NM_12GEV_20pC_N1_','PXY1',2,35)
print set_out_name('S1_0.25NM_12GEV_20pC_N1_',2)


# In[ ]:

def set_pzname(namg,ifb):
    """
    Set the file name with gain curve data as it is defined in FAST data set
    """
    # suppose that undulator position identifier, nz, has 2 digits
    name = namg+'PZ1_'+str(ifb).zfill(4)+'000.RES'
    return name


# In[ ]:

def set_in_felname(namg,ifb,nz):
    """
    Set input file name with FEL wavefront in cartezian coordinates
    FAST: this is output of post-processing code <pproc-fast2xy>.exe
    """
    name = namg+'T'+str(ifb).zfill(4)+str(nz).zfill(3) 
    return name


# In[ ]:

def parse_toe_file(f_name):
    """ Parse <prefix>T*.RES file to list of strings """
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


# In[ ]:

def create_numpy_array_from_rows(rows,slices=None):
    # slice size (Re, Im)
    N=len(rows[0])/2
    
    if slices is None:
        slices=range(len(rows)/N)
    slice_count=len(slices)

    y = numpy.array(rows,dtype='float32').reshape((slice_count,N,N,2))
    
    if doPrint: print 'y=rows.reshape((nSlices,nx,nx,2)',y.shape
    y=numpy.swapaxes(y,0,2)
    if doPrint: print 'return swapaxes(y,0,2)',y.shape
    return y


# In[ ]:

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


# In[ ]:

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


# In[ ]:

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


# In[ ]:

def update(in_fast2xydat,fast2xydat='PPROC-FAST2XY_2013.DAT',trd1=0.,trd2=None,
           nxy=None,nskip=None,ifb=None,nzc=None,namg=None):
    """
    Read a parameter file in_fast2xydat from FAST data set, 
    update the file with new values specified in the argument list, and save it as fast2xydat.
    This parameter file is used for extracting wavefront data output in cartesian coordinates.
    """
    f_in = open(in_fast2xydat,'r')
    a = f_in.readlines()  
    #if doPrint: print a
    strInputPar = numpy.empty(len(a), dtype=object)
    strComment  = numpy.empty(len(a), dtype=object)
    if doPrint: print '==Beforehand=='
    for idx in range(0,7): 
        strInputPar[idx] =  a[idx].split('#',1)[0]
        strComment [idx]  = a[idx].split('#',1)[1]
        if doPrint: print strInputPar[idx].rstrip()+strComment[idx].rstrip()
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
        print 'WARNING: File prefix is changed from %s to %s!' %(strInputPar[6].strip(),namg)
        strInputPar[6] = namg
    else:
        namg = strInputPar[6].strip()
    if doPrint: print '==Afterwards=='
    for idx in range(0,7): 
        a[idx] = strInputPar[idx].strip()+' #'+strComment[idx].rstrip()
#    for idx in range(0,len(a)):
    if doPrint: 
        for idx in range(0,7):
            print a[idx]
    aout=a[0:31]
    numpy.savetxt(fast2xydat,aout,fmt='%s')  
    os.system('chmod a+rw '+fast2xydat)
    return namg,ifb,nzc


# In[ ]:

def fill_wf_history(wf_struct,fast_readme,fast_internal):
    """
    fill history info/ data sets:
    contact/ info, data_description from readme file
    """
    f_readme = open(fast_readme,'r')
    f_internal = open(fast_internal,'r')
    wf_struct['history/parent/info']={
        'contact':(r'''Mikhail Yurkov
        Evgeny Schneidmiller

        Deutsches Elektronen Synchrotron (DESY)
        Notkestrasse 85
        22607 Hamburg
        Germany
                                             
        E-mail: mikhail.yurkov@desy.de
        E-mail: evgeny.schneidmiller@desy.de
        Tel. +49 40 8998 2676''','s'),
        'data_description':(f_readme.readlines(),'s'),
        'method_description':(r'''
        FEL simulation code FAST (E.A. Schneidmiller and M. V. Yurkov)

        FAST is generic name for a set of codes for analysis of the FEL amplification 
        process in the framework of 1-D and 3-D models using different techniques
        described in [1-5]. Analytical techniques implemented in these codes allow to
        analyze beam radiation modes (eigenvalue equation), and amplification process
        in the linear stage of amplification (initial-value problem). Numerical
        simulation codes allow to simulate FEL process using both, steady-state and
        time-dependent models. Algorithm of three-dimensional simulation code FAST [6]
        takes into account all important physical effects: diffraction of radiation,
        slippage of radiation (time-dependent effects), space charge, emittance
        (betatron oscillations), and energy spread in the electron beam. An
        approximation a uniform electron beam focusing is used in the simulation code.
        Field solver uses expansion of the radiation in the azimuthal modes.
        Calculation of the radiation fields is performed using retarded potentials.
        The code allows to simulate amplification process in helical and planar
        undulators. Simulation of higher odd harmonics in the planar undulator
        geometry is possible [6]. The code has been thoroughly tested in the high gain
        exponential regime using analytical results for the beam radiation modes
        (complex eigenvalues and eigenfunctions) [1-3]. Start-up from the shot noise
        and linear stage of amplification in the code FAST has been also tested using
        three-dimensional analytical results for the initial-value problem [4].

        Start-up from the shot noise in the electron beam in the code FAST can be
        simulated with an artificial ensemble [7,8], and with tracing actual number of
        electrons in the beam [9]. The latter option is straightforward and
        transparent: simulation procedure corresponds to real electrons randomly
        distributed in full 6D phase space. This allows us to avoid any artificial
        effects arising from standard procedures of macroparticle loading as it
        described in [9].

        [1] E.L. Saldin, E.A. Schneidmiller, M.V. Yurkov, ``The Physics of Free Electron Lasers'' (Springer-Verlag, Berlin, 1999).
        [2] E.L. Saldin, E.A. Schneidmiller and M.V. Yurkov, Nucl. Instrum. and Methods A 475 (2001) 86.
        [3] E.A. Schneidmiller and M.V. Yurkov, Phys. Rev. ST Accel. Beams 15 (2012) 080702.
        [4] E.L. Saldin, E.A. Schneidmiller, and M.V. Yurkov, Opt. Commun. 186(2000)185.
        [5] E.A. Schneidmiller, and M.V. Yurkov, Proc. FEL 2012 Conference, http://accelconf.web.cern.ch/AccelConf/FEL2012/papers/ mopd08.pdf.
        [6] E.L. Saldin, E.A. Schneidmiller, and M.V. Yurkov, Nucl. Instrum. and Methods A 429(1999)233.
        [7] C. Penman, B.W.J. McNeil, Optics Comm. 90 (1992) 82.
        [8] W.M. Fawley, Phys. Rev. STAB 5(2002)070701.
        [9] E.L. Saldin, E.A. Schneidmiller, and M.V. Yurkov,Opt. Commun. 281(2008)1179.
          ''','s'),
        'package_version':(r'''FAST v2.06''','s')}
    wf_struct['history/parent/detail']={
                                        'params':(f_internal.readlines(),'s')}


# In[ ]:

def fill_wf_params(wf_struct,params,nrows):
    """
    create and fill h5 wavefront params/ data sets in accordance with glossary
    """
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
    #slMin = slMin0 + (nStart-1)*slStep
    #slMax = slMin0 + (nEnd - 1)*slStep
    slMin = -slStep*(nEnd - 1)/2
    slMax =  slStep*(nEnd - 1)/2
    if doPrint: print 'wl [nm], Eph [keV]',wl*1e9,photonEnergy*1e-3
    if doPrint: print 'nStart,nEnd: ',nStart,nEnd
    if doPrint: print 'nSlices,slMin,slMax [fs]: ',nEnd,slMin*1e15,slMax*1e15
    if doPrint: print 'nx,xMin,xMax [um]',nx,xMin*1e6,xMax*1e6  
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
    return RK1/(xStep*1e3)**2,nStart,nEnd


# In[ ]:

def add_wf_attributes(fname0):
    """
    add attributes to wavefront datasets
    using WPG/SRWlib interface   
    """
    in_fname   = fname0
    bare_fname = fname0+'bare.h5'
    if doPrint: print('Loading wavefront data from the file:     '+in_fname)
    wf_struct=Wavefront()
    wf_struct.load_hdf5(in_fname)
    wfr = wf_struct._srwl_wf
    wf_struct = Wavefront(wfr)
    if doPrint: print('Saving the wavefront data with attributes:'+bare_fname)
    wf_struct.store_hdf5(bare_fname)
    if doPrint: print('Replacing data with attributes from  '+bare_fname)
    with h5py.File(bare_fname) as h2:
        with h5py.File(in_fname) as h1:
            try:
                del h1['params']  # delete group
            except KeyError:
                pass
            h2.copy('params',h1) #copy h2['params'] to h1


# In[ ]:

def convert_fast2h5(fel_data_path,fast2xyexe,fast2xydat,fast_readme,fast_internal,namg,ifb,nzc):
    """
    convert FAST native output data to hdf5 file for wavefront propagation: 
    1) convert output to binary file with wavefront in cartesian coordinate,
    2) put the wavefront data in h5 data/ dataset
    3) add params/ and history/ fields in accordance with the glossary
    4) save it into h5 file with name <prefix><nzc><ifb>.h5
    """
    ifname = set_ABCDname(namg,'T',ifb,nzc)
    fulltname=os.path.join(fel_data_path,ifname)
    if doPrint: print 'fulltname:',fulltname
    if doPrint: print 'Create sym link:'
    os.system("ln -s "+fulltname+" fort.4")
    outp=os.popen(fast2xyexe).read()
    if doPrint: print 'Remove sym link:'
    os.system("rm fort.4")
    if doPrint: print 'fast2xy.exe:\n',outp

    #loading data from text toe*.res file
    nharm=NHARM # for fundamental
    name_FXY = set_ABCDname(namg,'FXY'+str(nharm),ifb,nzc)
    print 'Reading FEL data from:'+ name_FXY+' ...'
    params, rows=parse_toe_file(name_FXY)
    if doPrint: print '...done'

    #wavefront structure based on glossary
    wf_struct={'version':(0.1,'f')}
    RK1,nStart,nEnd = fill_wf_params(wf_struct, params,numpy.size(rows))
    #build numpy arrays from list of rows
    wf_data=create_numpy_array_from_rows(rows,slices=range(nStart-1,nEnd))
    del rows
    wf_struct['data']={
                       'arrEhor':(wf_data*numpy.sqrt(RK1),'f'),
                       'arrEver':(numpy.zeros(shape=wf_data.shape,dtype='float32'),'f')
                       }
    wf_struct['misc']={
                       }
    fill_wf_history(wf_struct,fast_readme,fast_internal)
    f_in = open(fast2xydat,'r')
    
    e_data =  numpy.loadtxt(set_ABCDname(namg,'E'+str(nharm),ifb,nzc))
    ff_data = numpy.loadtxt(os.path.join(fel_data_path,set_ABCDname(namg,'F'+str(nharm),ifb,nzc)))
    nf_data = numpy.loadtxt(os.path.join(fel_data_path,set_ABCDname(namg,'N'+str(nharm),ifb,nzc)))
    pz_data_file = set_pzname(namg,ifb)
    pz_data = numpy.loadtxt(os.path.join(fel_data_path,pz_data_file))
    wf_struct['history/parent/misc']={
                                             'nzc':(nzc,'f'),
                                             'FAST2XY_DAT':(f_in.readlines(),'s'),
                                             'temporal_struct':(e_data,'f'),
                                             'spot_size':(nf_data,'f'),
                                             'angular_distribution':(ff_data,'f'),
                                             'gain_curve':(pz_data,'f'),
                                             }
                                             #'fwhm_curve':(fwhm_data,'f') #<- not defined for long pulse data
    
    #store wavefront in hdf5 file
    fname0 = set_out_name(namg+str(nzc).zfill(3),ifb)
    if doPrint: print "Store wavefront in hdf5 file: "+fname0
    store_wavefront_hdf5(wf_struct,fname0)
    add_wf_attributes(fname0)
    return fname0


# In[ ]:

def main():
    # typical command line parameters:
    # fast2xy_new.py -i'PPROC-FAST2XY_2013_LP.DAT' --time-start=3. --skip-nslices=8 --zc-point-num=33 --jmax=2
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-o", "--output-dir",          dest="out_dir", help="Output directory", )
    parser.add_option("-d", "--data-path",           dest="datapath",help="FEL data path", )
    parser.add_option("-f", "--file-id",             dest="f_id",    help="ID of the first of output files: FELsource_out_<ID+{0..JMAX-1}>", )
    parser.add_option("-j", "--jmax",                dest="jmax",    help="how many output pulses should be provided (1 by default)", )
    parser.add_option("-t", "--time-start",          dest="trd1",    help="Start time value for reading the pulse, fs")
    parser.add_option("-q", "--time-end",            dest="trd2",    help="End time value for reading the pulse, fs")
    parser.add_option("-n", "--nxy",                 dest="nxy",     help="Number of xy nodes")
    parser.add_option("-s", "--skip-nslices",        dest="nskip",   help="skip slices (1 - no skip)")
    parser.add_option("-r", "--run-number",          dest="ifb",     help="XXXX Run number (Input file name TXXXXYYY.RES)")
    parser.add_option("-z", "--zc-point-num",        dest="nzc",     help="YYY	Number of output point along undulator, zc")
    parser.add_option("-p", "--prefix",              dest="namg",    help="Prefix of FEL data files")
    parser.add_option("-C", "--e-charge",            dest="e_charge",help="Suffix for e-bunch charge")

    (options, args) = parser.parse_args()
       
    if not options.out_dir:   # if output directory not given
        parser.error('Output directory is not specified')
    else:
        out_dir=options.out_dir

    if not options.datapath:   # if output directory not given
        fel_data_path='/pnfs/desy.de/exfel/disk/XFEL/2013/SIM/FAST/'
    else:
        fel_data_path=options.datapath
        
    if not options.f_id:   # if time value not given
        f_id = 1 #default value
    else:
        f_id=int(float(options.f_id))

    if not options.trd1:   # if time value not given
        trd1 = 0.
    else:
        trd1=float(options.trd1)

    if not options.trd2:   # if time value not given
        trd2 =trd1 + 9.*2
    else:
        trd2=float(options.trd2)

    if not options.nzc:   # if nzc value not given
        nzc=35 #default value
    else:
        nzc=int(float(options.nzc))

    if not options.nskip:   # if nskip value not given
        nskip=1 # no skip
    else:
        nskip=int(float(options.nskip))
        
    if not options.jmax:   # if jmax value not given
        jmax=1 #default value
    else:
        jmax=int(float(options.jmax))

    if not options.namg: # if dir prefix is not given
        dir_prefix = '2014-05_XFEL_5keV_12GeV_'
    else:
        dir_prefix=options.namg

    if not options.e_charge: # if prefix is not given
        e_charge = '100pC'
    else:
        e_charge=options.e_charge

    fel_data_dir=dir_prefix+e_charge+'_N1'
    in_fast2xydat='PPROC-FAST2XY_2013.DAT'
    nharm = NHARM
    
    fel_data_path=fel_data_path+fel_data_dir
    in_fast2xydat = os.path.join(fel_data_path,in_fast2xydat) 
    fast_readme = os.path.join(fel_data_path,dir_prefix+e_charge+'_readme.txt')
    fast_internal = os.path.join(fel_data_path,'FAST_2013.DAT')
    
    #doPrint = True # switch on/off debug printing    
    doCopyRes = True # switch on/off copying results into FELsource/prop_in_XXX.h5    
    if isS2E: 
        work_dir ='/data/S2E/data/FELsource/';doCopyRes=True
    else:
        out_dir = '/diskmnt/a/lsamoylv/sim_data/FELsource'
        work_root_dir = '/diskmnt/a/lsamoylv/FELsource/' #working directory to process the data 
        work_dir = work_root_dir+fel_data_dir;
        
    tmp_dir = work_dir+'/tmp'
    mkdir_p(out_dir); mkdir_p(work_dir); mkdir_p(tmp_dir)    
    if doPrint: 
        print 'The data will be processed in \n'+ work_dir+'/'
        print 'All temporary files will be saved in \n'+ tmp_dir+'/'
    os.chdir(work_dir)
    fast2xyexe='pproc-fast2xy-2013-v2-06-wo-fname.exe';fast2xydat='PPROC-FAST2XY_2013.DAT'
    for idx in range(0,jmax):
        ifb = idx+f_id
        namg,ifb,nzc = update(in_fast2xydat,fast2xydat,trd1=trd1,trd2=trd2,ifb=ifb,nskip=nskip,nzc=nzc)
        if doPrint: print 'namg,ifb,nzc:',namg,ifb,nzc
        shutil.copy(os.path.join(work_dir,fast2xydat), tmp_dir)
        shutil.copy(os.path.join(work_dir,fast2xyexe), tmp_dir)  
        os.chdir(tmp_dir)
        in_fname=convert_fast2h5(fel_data_path,fast2xyexe,fast2xydat,fast_readme,fast_internal,namg,ifb,nzc)
        out_fname     = set_out_name(in_fname,       (f_id+idx))
        prop_in_fname = set_out_name('FELsource_out',(f_id+idx))
        
        if doPrint: print 'in_fname,out_fname,prop_in_fname:',in_fname,out_fname,prop_in_fname
        os.system('chmod a+rw '+ tmp_dir+'/*.*')    
        
        shutil.copy(in_fname, os.path.join(work_dir,out_fname))
        shutil.rmtree(tmp_dir)

        print 'The result hdf5 file  '+out_fname+' will be moved to '
        print out_dir+'/'+ prop_in_fname            
        shutil.move(os.path.join(work_dir,out_fname), 
                            os.path.join(out_dir,prop_in_fname))
        os.chdir(work_dir)


# In[ ]:

#$
if not isIpynb:
    if __name__ == '__main__':
        main()
else:
    # (2) the FAST data directory name pns/desy.de/... is a parameter, 
    # (3) before, outside the for-loop, copy PPROC-FAST2XY_2013.DAT from DCACHE to work_dir 
    
    # typical command line parameters:
    # fast2xy_new_SP_mod.py -i'PPROC-FAST2XY_2013.DAT' --prefix='2014-05_XFEL_5keV_12GeV_'\
    #                --time-start=1.5 --time-end=6.5 --skip-nslices=7 \
    #                --e-charge='20pC' --zc-point-num=35 --jmax=1

    fel_data_path='/pnfs/desy.de/exfel/disk/XFEL/2013/SIM/FAST/'
    dir_prefix='2014-05_XFEL_5keV_12GeV_';
    in_fast2xydat='PPROC-FAST2XY_2013.DAT'
    e_charge='20pC';f_id=1;trd1=1.5;trd2=6.5;nskip=7;nzc=35;jmax=1
    fel_data_dir=dir_prefix+e_charge+'_N1'
    nharm = NHARM
    
    fel_data_path=fel_data_path+fel_data_dir
    in_fast2xydat = os.path.join(fel_data_path,in_fast2xydat) 
    fast_readme = os.path.join(fel_data_path,dir_prefix+e_charge+'_readme.txt')
    fast_internal = os.path.join(fel_data_path,'FAST_2013.DAT')
    doPrint = True # switch on/off debug printing    
    doCopyRes = True # switch on/off copying results into FELsource/prop_in_XXX.h5    

    if isIpynb: 
        out_dir = '/diskmnt/a/lsamoylv/sim_data/FELsource'
        work_root_dir = '/diskmnt/a/lsamoylv/FELsource/' #working directory to process the data 
        work_dir = work_root_dir+fel_data_dir
    else:
        work_dir ='/data/S2E/data/FELsource/'
    tmp_dir = work_dir+'/tmp'
    mkdir_p(out_dir); mkdir_p(work_dir); mkdir_p(tmp_dir)    
    if doPrint: 
        print 'The data will be processed in \n'+ work_dir+'/'
        print 'All temporary files will be saved in \n'+ tmp_dir+'/'
    os.chdir(work_dir)
    fast2xyexe='pproc-fast2xy-2013-v2-06-wo-fname.exe';fast2xydat='PPROC-FAST2XY_2013.DAT'
    for idx in range(0,jmax):
        ifb = idx+f_id
        namg,ifb,nzc = update(in_fast2xydat,fast2xydat,trd1=trd1,trd2=trd2,ifb=ifb,nskip=nskip,nzc=nzc)
        if doPrint: print 'namg,ifb,nzc:',namg,ifb,nzc
        shutil.copy(os.path.join(work_dir,fast2xydat), tmp_dir)
        shutil.copy(os.path.join(work_dir,fast2xyexe), tmp_dir)  
        
        os.chdir(tmp_dir)
        in_fname=convert_fast2h5(fel_data_path,fast2xyexe,fast2xydat,fast_readme,fast_internal,namg,ifb,nzc)
        out_fname     = set_out_name(in_fname,       (f_id+idx))
        prop_in_fname = set_out_name('FELsource_out',(f_id+idx))
        
        if doPrint: print 'in_fname,out_fname,prop_in_fname:',in_fname,out_fname,prop_in_fname
        os.system('chmod a+rw '+ tmp_dir+'/*.*')    
        
        shutil.copy(in_fname, os.path.join(work_dir,out_fname))
        shutil.rmtree(tmp_dir)
    
        if doCopyRes:
            print 'The result hdf5 file  '+out_fname+' will be copied/moved to '
            print out_dir+'/'+ prop_in_fname
            if jmax == 1:
                shutil.move(os.path.join(work_dir,   out_fname), os.path.join(out_dir,prop_in_fname))
            else:
                shutil.copy(os.path.join(work_dir,   out_fname), os.path.join(out_dir,prop_in_fname))
        else:
            print out_fname, prop_in_fname
            print '... done'
        os.chdir(work_dir)


# In[ ]:

#%tb


# In[ ]:

#ls 


# In[ ]:

NHARM


# In[ ]:



