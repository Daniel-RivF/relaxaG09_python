#!/usr/bin/python

import tools_py_relx as tools
import os.path

filename =  raw_input('Name of the .log file (including .log): ')
step = float(raw_input('Displacement norm:  '))
nsteps     = int(raw_input('Number of steps:  '))
filename2   = raw_input('name of the .com file (generated at output):  ')
charg  = raw_input('charge of the molecule:  ')
mult  = raw_input('spin multiplicity:  ')
route = raw_input('route section (including #p):  ')

Zs =   tools.parseZxyz(filename)[0]
xyzs = tools.scan(filename, step,nsteps)
tools.writer(Zs,xyzs, filename2, route, charg, mult)



if os.path.isfile(filename2):
    print 'File was generated'
else:
    print 'Error: the .com file was not generated'

