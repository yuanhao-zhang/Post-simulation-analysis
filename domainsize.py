# -------------------------------------------------------

from re import I, X
import sys,string,fileinput
#wfrom tkinter import Menubutton
import numpy as np
from math import *

# INPUT PARAMETERs
nconf =    50
nskip =    50 #nskip =0 means starting from the first frame
M=np.array([0,1.0,1.0,1.0,1.0])
M_ps =9287
M_poem=9287

nbin =     26   # number of bins in a layer
nlayers =  2
mg = nbin*nlayers
half =     True
infiles = ['3_ht.lammpstrj']
if half:
    file = 'klong.data'
else:
    file = 'denprof.csv'


ddnum1 = np.zeros(nbin+1,np.float32)
ddnum2 = np.zeros(nbin+1,np.float32)
ddnum3 = np.zeros(nbin+1,np.float32)
ddnum4 = np.zeros(nbin+1,np.float32)

def process_bar(num, total):
    rate = float(num)/total
    ratenum = int(100*rate)
    r = '\r[{}{}]{}%'.format('*'*ratenum,' '*(100-ratenum), ratenum)
    sys.stdout.write(r)
    sys.stdout.flush()

def dennum():
    for ii in range(1,natoms+1):
        xr = xc[ii]
        if xr >= xcm1:
            ig = int((xr - xcm1) / drg) + 1
            if typea[ii] == 1:
                dnum1[ig] += 1
            elif typea[ii] == 2:
                dnum2[ig] += 1
            elif typea[ii] == 3:
                dnum3[ig] += 1
            elif typea[ii] == 4:
                dnum4[ig] += 1

        elif xr < xcm1:
            ig = int((xr+xbox - xcm1) / drg) + 1
            if typea[ii] == 1:
                dnum1[ig] += 1
            elif typea[ii] == 2:
                dnum2[ig] += 1
            elif typea[ii] == 3:
                dnum3[ig] += 1
            elif typea[ii] == 4:
                dnum4[ig] += 1

IN = fileinput.input(infiles)
for loopnum in range(0,nskip):
    IN.readline()
    IN.readline()
    IN.readline()
    line = IN.readline()      # number of atoms
    fields = str.split(line)
    natoms = int(fields[0])
    IN.readline()
    IN.readline()      
    IN.readline()      
    IN.readline()      
    IN.readline()
    for j in range(natoms):
        IN.readline()

avgLz = 0 
xcom0=ycom0=zcom0=xcom=ycom=zcom=Mtotal=0
st=np.zeros((nconf,1))

for kconf in range(0,nconf):
    process_bar(kconf+1,nconf)
    IN.readline()
    line = IN.readline()
    st[kconf] = int(line.split()[0]) 
  
    IN.readline()
    line = IN.readline()      # number of atoms
    fields = str.split(line)
    if kconf == 0:
        natoms = int(fields[0])
        dim=natoms+1
        zc = [0]*dim
        zz = [0]*dim
        yc = [0]*dim
        xc = [0]*dim
        xb = [0]*dim
        yb = [0]*dim
        zb = [0]*dim
        x0 = [0]*dim
        y0 = [0]*dim
        z0 = [0]*dim
        delt = [0]*dim
        nc=[0]*8  
        msd=np.zeros((nconf, 3), dtype=float)      
    cirzcm = 0
    cirycm = 0
    xcom=0
    ycom=0
    zcom=0
    typea = [0]*dim
    IN.readline()
    line = IN.readline()      
    [xm,xp] = map(float,line.split())
    line = IN.readline()      
    [ym,yp] = map(float,line.split())
    line = IN.readline()
    [zm,zp] = map(float,line.split())
    line = IN.readline()
    xbox = xp - xm
    ybox = yp - ym
    zbox = zp - zm
    drg = xbox/mg

    vol = xbox*ybox
    
    xboxi = xbox/nlayers
    for j in range(1,dim):
        line = IN.readline()
        [ii,molj,typej,q,x1,x2,x3,n1,n2,n3] = str.split(line)
        k = int(ii) # atom index
        typea[k] = int(typej)
        xc[k] = xbox*float(x1)
        xc[k] = xbox*float(x1)+xbox*float(n1)
        yc[k] = ybox*float(x2)+ybox*float(n2)
        zz[k] = zbox*float(x3)+zbox*float(n3)
        if kconf == 0:
            x0[k]=xbox*float(x1)+xbox*float(n1)
            y0[k]=ybox*float(x2)+ybox*float(n2)
            z0[k]=zbox*float(x3)+zbox*float(n3)
            if int(typej) <=2:
                xcom0 += xbox*float(x1)*M[int(typej)]
                ycom0 += ybox*float(x2)*M[int(typej)]
                zcom0 += zbox*float(x3)*M[int(typej)]
                Mtotal += M[int(typej)]
            
        if (kconf>=1 and int(typej) <=2):
            xcom +=xbox*float(x1)*M[int(typej)]/Mtotal
            ycom +=ybox*float(x2)*M[int(typej)]/Mtotal
            zcom +=zbox*float(x3)*M[int(typej)]/Mtotal

        if xc[k] < 0:
            xc[k] += xbox
        elif xc[k] > xbox:
            xc[k] -= xbox
        for ii in range(nlayers):      #combine the four layers of lamellae
            if xboxi*ii < xc[k] <= xboxi*(ii+1):
                xc1 = xc[k] - xboxi*ii
        if typea[k] == 2:
            theta = (2*pi/xboxi)*xc1 #circling z-axis for one layer (combined with others)
            cirz = (xboxi/(2*pi))*sin(theta) #x coordinate in cylindrical coordinates(xytheta), R=zbox/(2*pi)
            ciry = (xboxi/(2*pi))*cos(theta) #y coordinate in cylindrical coordinates(xytheta), R=zbox/(2*pi)
        if typea[k] == 2:
            cirzcm += cirz/(M_ps)
            cirycm += ciry/(M_ps)
    if kconf == 0:        
        xcom0/=Mtotal
        ycom0/=Mtotal
        zcom0/=Mtotal
        for i in range(1,4):
            nc[i] = typea.count(i)

    thetacm = atan2(-cirzcm,-cirycm) + pi
    xcm1 = (xboxi/(2*pi))*thetacm
    avgLz += xboxi/nconf
    dnum1 = [0]*(mg+1)
    dnum2 = [0]*(mg+1)
    dnum3 = [0]*(mg+1)
    dnum4 = [0]*(mg+1)         
    dennum()
    zbreal= max(zz)-min(zz)
    
    binsize = (vol*zbreal)/mg
    print(zbreal,vol*zbreal)


    if kconf == 2:
        print(dnum1)
    for ig in range(1,mg+1):
        bin_id = ((ig-1)%nbin)+1
        ddnum1[bin_id] += float(dnum1[ig])/binsize
        ddnum2[bin_id] += float(dnum2[ig])/binsize
        ddnum3[bin_id] += float(dnum3[ig])/binsize
        ddnum4[bin_id] += float(dnum4[ig])/binsize

print ('\r')
print ('Average spacing:', avgLz)

ddnum1 /= nconf*nlayers
ddnum2 /= nconf*nlayers
ddnum3 /= nconf*nlayers
ddnum4 /= nconf*nlayers

halfbin = int(nbin/2)+1

if half:
    # Save original half data for variance calculation
    ddnum1_left = ddnum1[halfbin:]
    ddnum1_right = ddnum1[1:halfbin][::-1]
    ddnum2_left = ddnum2[halfbin:]
    ddnum2_right = ddnum2[1:halfbin][::-1]
    ddnum3_left = ddnum3[halfbin:]
    ddnum3_right = ddnum3[1:halfbin][::-1]
    ddnum4_left = ddnum4[halfbin:]
    ddnum4_right = ddnum4[1:halfbin][::-1]
  
    # Calculate symmetrized averages
    ddnum1 = (ddnum1_left + ddnum1_right)/2
    ddnum2 = (ddnum2_left + ddnum2_right)/2
    ddnum3 = (ddnum3_left + ddnum3_right)/2
    ddnum4 = (ddnum4_left + ddnum4_right)/2

    
    # Calculate variances
    var_num1 = np.std([ddnum1_left, ddnum1_right], axis=0)
    var_num2 = np.std([ddnum2_left, ddnum2_right], axis=0)
    var_num3 = np.std([ddnum3_left, ddnum3_right], axis=0)
    var_num4 = np.std([ddnum4_left, ddnum4_right], axis=0)

    
    nbin = halfbin-1
else:
    ddnum1 = np.append(ddnum1[halfbin:],ddnum1[1:halfbin])
    ddnum2 = np.append(ddnum2[halfbin:],ddnum2[1:halfbin])
    ddnum3 = np.append(ddnum3[halfbin:],ddnum3[1:halfbin])

OUT = open(file, 'w')
if half:
#    OUT.write("z, den(EO), den(PS), den(Li), D(Li), D(Tf1), D(Tf2)\n")
    for i in range(nbin):
        z = (i+0.5)/(2*nbin)
#        OUT.write("%1.4f %8.9f %8.9f \n" % 
#                 (z, ddnum1[i], ddnum2[i]))
    slop = np.gradient(ddnum2,z)    
    print(slop[int(nbin/4)])
    print(slop)
#    OUT.write('\n')
#    OUT.write("var(EO), var(PS), var(Li), var_D(Li), var_D(Tf1), var_D(Tf2)\n") 
#    for i in range(nbin):
#        z = (i+0.5)/(2*nbin)
#        OUT.write("%1.4f %8.9f %8.9f %8.9f %8.9f %8.9f %8.9f\n" % 
#                 (z,var_num1[i], var_num2[i], var_num4[i], var_li[i], var_tf1[i], var_tf2[i]))    
else:
    OUT.write("z, den(EO), den(PS), den(Li), D(Li), D(Tf1), D(Tf2)\n")
    for i in range(nbin):
        z = (i+0.5)/nbin-0.5
        OUT.write("%1.4f, %8.9f, %8.9f, %8.9f, %8.9f, %8.9f, %8.9f\n" % 
                 (z, ddnum1[i], ddnum2[i], ddnum4[i], ddli[i], ddtf1[i], ddtf2[i]))

OUT.write('\n')
OUT.write('\n')
OUT.close()

