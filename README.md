# DensityTool
DensityTool v0.2 is a FORTRAN program, designed to compute the local density of states and local spin density of states from the output of the VASP package. 

Created by Lucas Lodeiro and Tomáš Rauch.

lucas.lodeiro@ug.uchile.cl

tomas.rauch@uni-jena.de

If you use the code, please cite:
L. Lodeiro and T. Rauch: DensityTool: A post-processing tool for space- and spin-resolved
density of states from VASP, Computer Physics Communications 277, 108384 (2022),
https://doi.org/10.1016/j.cpc.2022.108384



### Description of files / folders:

DENSITYTOOL.F90 - program source code

DENSITYTOOL.IN - sample input

LICENSE - MIT License

Plot_LDOS.py and Plot_LSDOS.py - sample script for plotting the data generated with DensityTool

User_Manual.pdf - User Manual

Examples - folder with examples, including input and output data, explained in contained README file



### Installation

The program can be compiled using a FORTRAN compiler, e.g.:
ifort -O2 -o DENSITYTOOL.X DENSITYTOOL.F90
or
gfortran -O3 -o DENSITYTOOL.X DENSITYTOOL.F90



### Execution

The executable DENSITYTOOL.X has to be executed in a folder with the VASP output and the input file, DENSITYTOOL.IN, if used, as
/PROGRAMFOLDER/DENSITYTOOL.X < DENSITYTOOL.IN > DENSITYTOOL.OUT
In DENSITYTOOL.IN the parameters of the calculation can be specified. If DENSITYTOOL.IN is not used, the parameters for the calculation are asked by the program and can be entered manually.
