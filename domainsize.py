import numpy as np
import math
import fileinput

nconf = 1
filename = "3_ht.lammpstrj"
zcount =np.zeros((1000))
def MSD(nconf):
    conf = 1
    rs = open('result1.data', 'w')
    step = 0
    mod=1
    k=0
    zmax=zmin=zc=count=0
    fp=fileinput.input(filename)
    while True:
        Line=fp.readline()
        if Line=='':
            break
        k=k+1       
        if mod==1:              
            if 'ITEM: NUMBER OF ATOMS' in Line:
                Line=fp.readline()
                Natom=int(Line)                         
            if 'ITEM: BOX BOUNDS' in Line:
                Line=fp.readline()
                [xmin, xmax] = map(float,Line.split())
                Line=fp.readline()                   
                [ymin, ymax] = map(float,Line.split())
                Line=fp.readline()
                [zmin, zmax] = map(float,Line.split())                    
            if "ITEM: ATOMS" in Line:                   
                step = step + 1
                mod=2   
                Line=fp.readline()  
        if "ITEM: TIMESTEP" in Line:
            mod=3                           
        if mod==3:
            z=zmax-zmin
            if step >=conf:
                zc+=zmax-zmin
                count+=1
                zcount[step]=(zmax-zmin)*0.35/2
            mod=1
            if step >= nconf:
                rs.write('%4i %7.5f\n' % (step, z) )
                print("doing step: %3i" %step)
    print("The average domain size after %1ith frame is: %4.3f nm" %  (nconf, np.average(zcount[nconf:step]))   )
    print("The standardard dev is: %4.3f nm" %np.std(zcount[nconf:step]))


MSD(nconf)
  


    
