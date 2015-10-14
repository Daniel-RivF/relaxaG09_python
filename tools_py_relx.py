#!/usr/bin/python
import numpy as np
from numpy.linalg import *
import os
from math import sin, cos

def parseALL(filename,first,second):
    parsing = False
    name = []
    f = open(filename)
    for line in f.readlines():
        if line.startswith(first):
            parsing = True
            continue
        elif line.startswith(second):
            parsing = False
        if parsing:
            a = line.split()
            name.extend(a)
#    name = map(float,name)
    return name

def parseZ(filename):
    return parseALL(filename,"Atomic numbers","Nuclear charges")

def parsexyz(filename):
    lst = parseALL(filename," Number     Number       Type","         ")
    lines = map(float,lst[1:-1])
    xyz = [lines[x:x+3] for x in xrange(0, len(lines), 3)]
    b = xyz[1::2]
    Z = map(int,lines[1::6])
    return (Z,np.array(b)) 
    

def parseDCUGD(filename):
    aa = map(float,parseALL(filename," Derivative Coupling"," Unscaled Gradient Difference"))
    bb = map(float,parseALL(filename," Unscaled Gradient Difference", " Gradient of iOther State"))
    aaa = [aa[x:x+3] for x in xrange(0, len(aa), 3)] 
    bbb = [bb[x:x+3] for x in xrange(0, len(bb), 3)]
    return ( np.array(aaa), np.array(bbb) )


   
def scancoIn2(filename,normd,n,gauss_inp,route,carg,mult):
    'filename: name of the file to read (log), normd : displacement norm \
      n = number of points , gauss_inp: name of the generated input, charge and multiplicity'
    normd = normd*1.889725989
    DC = parseDCUGD(filename)[0]
    UGD = parseDCUGD(filename)[1]
    DCn = DC / norm(DC)
    UGDn = UGD / norm(UGD)
    geom_CoIn_au = 1.889725989 * parsexyz(filename)[1]
    Zlist = parsexyz(filename)[0]
    esc = np.vdot(UGDn,DCn) 
    UGDT = UGDn - esc*DCn
    UGDTn = UGDT / norm(UGDT)
    #    UGDT = (UGD - DC) / norm( UGD - DC )
    alpha = (2 * 3.141592654) / n
    inpfile = open(gauss_inp,'w')
    p = '%'
    all_geoms=[]
    for i in range(n):
        chk = os.path.splitext(gauss_inp)[0]
        inpfile.write('%smem=4000Mb  \n' % p)
        inpfile.write('%snproc=2 \n' % p)
        inpfile.write('%schk=%s \n' % (p,chk))
        inpfile.write('%s \n' % route)
        inpfile.write('\n')
        inpfile.write('comment \n')
        inpfile.write('\n')
        inpfile.write('%s %s \n' %(carg,mult))
#        geom_p = ( geom_CoIn_au + ( cos(alpha*(i+1))*normd*DCn + sin(alpha* (i+1) )*normd*UGDTn) ) * 0.529177249   
        geom_p = ( geom_CoIn_au +  DCn*normd*cos(alpha*(i+1)) + UGDTn*normd*sin(alpha*(i+1)) ) * 0.529177249
        geoml = geom_p.tolist()
        aa = zip(Zlist,geoml)
        aaa = []
        for j in aa:
            aaa.append((map(str,[j[0]]+j[1])))
        for k in aaa:
            inpfile.write('%s %s %s %s\n' % (k[0],k[1],k[2],k[3]))
        inpfile.write('\n')
        inpfile.write(' 0.5       0.5\n')
        inpfile.write('--Link1-- \n')
        all_geoms.extend(geoml)
    inpfile.close()

    length=len(all_geoms)/n
    new_list = []
    m=0
    while m < len(all_geoms):
        new_list.append(all_geoms[m:m+length])
        m += length
    return new_list 



# listofgeoms is a nested list which contains all the geometries of the scan (returned by the previous function)

def workF(listofgeoms,filename,modF,atomF1,atomF2):
#    listofgeoms = scancoIn2(filename,normd,n,gauss_inp,route,carg,mult)
    geomCoIn_au = 1.889725989 * parsexyz(filename)[1]
    m = len(listofgeoms[0])
#    alldirections = []
    workl = []
    for i in listofgeoms:
        i = np.array(i)
        # Force vector is z
        z = np.zeros((m,3),dtype=float)
        newrow1 = ( (i[atomF1-1]-i[atomF2-1]) / norm(i[atomF1-1]-i[atomF2-1]) )* modF/82.387
        newrow2 = -newrow1
        z[atomF1-1] = newrow1
        z[atomF2-1] = newrow2
        Fmat = z 
 #       alldirections.extend(z)
        # Calculate the difference of the scan points with the
        displacement = (geomCoIn_au -  i*1.889725989  )
        work_au = np.vdot( Fmat,displacement)
        workl.append(work_au)
        # Modificar para que escriba  workl en un fichero.
    f = open('works.out','w')
    for i in workl:
        f.write('%s \n' % i )
    f.close()
    return

   #[0.00010110834137117427, 9.1813386449662184e-05, 7.3511330586352206e-05, 4.7987851734991272e-05, 1.7742639990230262e-05, -1.4255794359727751e-05, -4.4863735586282939e-05, -7.1074383239433966e-05, -9.0316244937549132e-05, -0.00010070705081810215, -0.0001012375120993741, -9.1867132589642573e-05, -7.3523827684323186e-05, -4.8009030993233425e-05, -1.7819116159512612e-05, 1.4098528095873198e-05, 4.4631044971167415e-05, 7.0800443439781368e-05, 9.0050986977326439e-05, 0.00010049708970987196]

