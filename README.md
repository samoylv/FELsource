Hello, this is the Start-to-End simulation workflow!

Sources hosted at git repository https://github.com/samoylv/prop.git

***** DESCRIPTION OF SCRIPTS *****
fast2h5_2013_s2e.py - python script that converts the FEL data
                /data/S2E/data/FELsource/SASE1_5keV_14GeV_T1001033.RES to
                HDF5 file SASE1_5keV_14GeV_FXY1_1001033_003fs_out.h5
the python script uses FORTRAN code to access the original FEL data
The script creates a temporary directory for intermediate files. 

[2014-05-09] Now the directory name is the same as base name for 
output FEL hdf5. It leads to collision: if the temporary directory and files 
exist, another user is not able to rewrite the files unless the output
FEL name is changed. The latter can be done by changing 
the ‘trd1’ (and trd2=trd1+9fs) value in fast2h5_2013_s2e.py 
==SHOULD BE FIXED ASAP==

fast2h5_2013_s2e.ipynb - IPython notebook equals to fast2h5_2013.py

fast2xy_2013_wo_name.for - FORTRAN file to take a chunk of long pulse data from
a *T*.RES file, to convert the chunk into cartesian coordinates, and
to save it in a binary file *FXY*.RES

The original *T*.RES file is specified as a symbolic link in fort.4 file,
the link is created in the python script fast2h5*.py, after building the file
name with parameters specified in config file FAST2XY_2013.DAT

Fortran file compiles with:

# print available modules installed
module avail
# add corresponding Intel FORTRAN compiler
module add intel/2013
# compiling and linking
ifort -c -132 pproc-fast2xy-2013-v2-06-wo-fname.for

ifort -o pproc-fast2xy-2013-v2-06-wo-fname.exe pproc-fast2xy-2013-v2-06-wo-fname.o
# remove the Intel module
module rm  intel/2013
