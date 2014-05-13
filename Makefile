.PHONY: all
	
all:
	# print available modules installed
	module avail
	# add corresponding Intel FORTRAN compiler
	module add intel/2013
	# compiling and linking
	ifort -c -132 fast2xy_2013_wo_name.for
	ifort -o fast2xy_2013_wo_name.exe fast2xy_2013_wo_name.o
	# remove the Intel module
	module avail intel/2013