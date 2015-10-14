#!/usr/bin/python

import tools_py_relx as tools
import os.path

name =  raw_input('Name of the .log file (including .log): ')
normd = float(raw_input('Displacement norm:  '))
n     = int(raw_input('Number of steps:  '))
inp   = raw_input('name of the .com file (generated at output):  ')
carg  = raw_input('charge of the molecule:  ')
mult  = raw_input('spin multiplicity:  ')
route = raw_input('route section (including #p):  ')


tools.scancoIn2(name,normd,n,inp,route,carg,mult)

if os.path.isfile(inp):
    print 'File was generated'
else:
    print 'Error: the .com file was not generated'

