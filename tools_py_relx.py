#!/usr/bin/python
import numpy as np
from numpy.linalg import norm
import os
from math import sin, cos


def chunks_optim(filename):
    with open(filename) as f:
        lines = f.read().splitlines()
    input_or = [i for i,k in enumerate(lines) if 'Input orientation:' in k]
    chunks = []
    a = 0 
    for ind in input_or:
        chunks.append(lines[a:ind])
        a = ind 
        if input_or.index(ind) == len(input_or)-1:
            chunks.append(lines[ind:])
    
    b = []
    for geom in chunks:
        for s in geom:
            if 'Optimization completed' in s:
             b.append(geom)

    return (b)

def parser_lists(data,first,second):
    whole_data = []
    parse = False
    for line in data:
        if str(first) in line:
            parse = True
        elif str(second) in line:
            parse = False
        if parse:
            whole_data.append(line)
    return whole_data

#####################################################
#####################################################

def parseZxyz(filename):
    lst = chunks_optim(filename)
    for i in lst:
        a = parser_lists(i,'Input orientation:','Distance matrix')[5:-1]
        b = [ i.split()[3:] for i in a]
        xyz = np.array([map(float,i) for i in b])
        Z = [ i.split()[1] for i in a]
    return (Z,xyz)


def parseDC(filename):
    dcugd = chunks_optim(filename)
    for i in dcugd:
        a_dc = parser_lists(i," Derivative Coupling"," Unscaled Gradient Difference")
        b_dc = [ i.split() for i in a_dc[1:] ]
        dc = np.array([map(float,i) for i in b_dc])
        return dc

def parseUGD(filename):
    dcugd = chunks_optim(filename)
    for i in dcugd:
        a_ugd = parser_lists(i," Unscaled Gradient Difference", " Gradient of iOther State")
        b_ugd = [ j.split() for j in a_ugd[1:] ]
        ugd = np.array([map(float,j) for j in b_ugd])
        return ugd

#### Reescribir desde aqui, 23-02-2016


def scan(filename, step,nsteps):
    step_au = step * 1.889725989
    dc_a = parseDC(filename)
    ugd_a = parseUGD(filename)
    dcn = dc_a / norm(dc_a)
    ugdn = ugd_a / norm(ugd_a)
    esc = np.vdot(ugdn,dcn)
    ugd_t = ugdn - esc * dcn
    ugd_tn = ugd_t / norm(ugd_t)
    xyz_CoIn_au = parseZxyz(filename)[1] * 1.889725989
    alpha = (2 * 3.141592654) / nsteps
    geoms_scan = []
    for i in range(nsteps):
        a1 = cos(alpha*(i+1))
        a2 = sin(alpha*(i+1))
        geom_i = ( xyz_CoIn_au + (dcn * step_au * a1) + (ugd_tn * step_au * a2) ) * 0.529177249
        geoms_scan.append(geom_i)
    return geoms_scan


def writer(Zs,xyzs, filename2, route, charg, mult):
    info = [zip(Zs,i) for i in xyzs] 
    p = '%'
    chkname = os.path.splitext(filename2)[0]
    with open(filename2, 'w') as f:
        for geom in info:
           f.write('%schk=%s \n' %(p,chkname) )  
           f.write('%snproc=4\n' % p)
           f.write('%smem=4000Mb \n' % p)
           f.write('%s \n '  % route)
           f.write('  \n')
           f.write(' comment\n')
           f.write('  \n')
           f.write ( '%s %s \n' %(charg,mult))
           for line in geom:
               f.write('%s %s \n' %( line[0], ' '.join(map(str,line[1]))))
           f.write('  \n')
           f.write('----Link1--\n')
    return

#def workF(listofgeoms,filename,modF,atomF1,atomF2):
##    listofgeoms = scancoIn2(filename,normd,n,gauss_inp,route,carg,mult)
#    geomCoIn_au = 1.889725989 * parsexyz(filename)[1]
#    m = len(listofgeoms[0])
##    alldirections = []
#    workl = []
#    for i in listofgeoms:
#        i = np.array(i)
#        # Force vector is z
#        z = np.zeros((m,3),dtype=float)
#        newrow1 = ( (i[atomF1-1]-i[atomF2-1]) / norm(i[atomF1-1]-i[atomF2-1]) )* modF/82.387
#        newrow2 = -newrow1
#        z[atomF1-1] = newrow1
#        z[atomF2-1] = newrow2
#        Fmat = z 
# #       alldirections.extend(z)
#        # Calculate the difference of the scan points with the
#        displacement = (geomCoIn_au -  i*1.889725989  )
#        work_au = np.vdot( Fmat,displacement)
#        workl.append(work_au)
#        # Modificar para que escriba  workl en un fichero.
#    f = open('works.out','w')
#    for i in workl:
#        f.write('%s \n' % i )
#    f.close()
#    return
