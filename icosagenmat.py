#!/usr/bin/python
import os
from sys import *
from math import *
from string import atof


def howtouse():
      print("./icosagenmat.py \n")
      print("./icosagenmat.py 90 0 0 \n")
      print("./icosagenmat.py 20 30 10 115.2 123.1 134.1\n")
      print("./icosagenmat.py 20 30 10 115.2 123.1 134.1 \"FIXR ON FIXB ON\" \n")
      print("matrices numbering is the same as viperdb.scripps.edu one \n")



Rad=pi/180.
Deg=1./Rad

xcen,ycen,zcen = 0.0,0.0,0.0
cenlist=[xcen,ycen,zcen]

rot = [0.]*10
rot2 = [0.]*10
rot3 = [0.]*10
rotr = [0.]*10

axis5= [0.,0.,0.]*13

icosamat=[[0., 0.,0.,0.,0.,   0.,0.,0.,0.,   0.,0.,0.,0. ] for ico in range(61)]

#icosamat[1][1][1]=1.
#icosamat[1][2][2]=1.
#icosamat[1][3][3]=1.

def icosamat_axis3permu_zxyz(ifrom,ito):
  icosamat[ito][1: 4]=icosamat[ifrom][9:12]  
  icosamat[ito][5: 8]=icosamat[ifrom][1: 4]  
  icosamat[ito][9:12]=icosamat[ifrom][5: 8]  

def icosamat_axis2_Zrot(ifrom,ito):
  icosamat[ito][1: 4]=map(lambda x: -x, icosamat[ifrom][1: 4] ) 
  icosamat[ito][5: 8]=map(lambda x: -x, icosamat[ifrom][5: 8] ) 
  icosamat[ito][9:12]=icosamat[ifrom][9:12]  

def icosamat_axis2_Xrot(ifrom,ito):
  icosamat[ito][1: 4]=icosamat[ifrom][1: 4]  
  icosamat[ito][5: 8]=map(lambda x: -x, icosamat[ifrom][5: 8] ) 
  icosamat[ito][9:12]=map(lambda x: -x, icosamat[ifrom][9:12] )




if __name__ == '__main__':
  a=(-1.+sqrt(5.))/2.

#  5-fold axises w/o normilize
  axis5[ 1]=[ 0., a , 1.] #1-2-3-4-5     
  axis5[ 2]=[ 1., 0., a ] #27-10-9-8-26
  axis5[ 3]=[ a , 1., 0.] #28-12-11-59-41

  axis5[ 4]=[ 0.,-a , 1.] #6-18-22-23-7
  axis5[ 5]=[ 1., 0 ,-a ] #55-54-47-57-60 
  axis5[ 6]=[-a , 1., 0.] #13-29-45-39-14

  axis5[ 7]=[ 0.,-a ,-1.] #48-53-33-37-49 -axis5[ 1]
  axis5[ 8]=[-1., 0 ,-a ] #36-38-40-44-50 -axis5[ 2]           
  axis5[ 9]=[-a ,-1., 0.] #31-21-20-35-34 -axis5[ 3]

  axis5[10]=[ 0., a ,-1.] #51-43-42-58-46 -axis5[ 4]
  axis5[11]=[-1., 0., a ] #15-30-19-17-16 -axis5[ 5]
  axis5[12]=[ a ,-1., 0.] #24-32-52-56-25 -axis5[ 6]


  ang2to5=atan2(a,1)
  lll=sin(ang2to5)*cos(pi/2.)
  mmm=sin(ang2to5)*sin(pi/2.)
  nnn=cos(ang2to5)

# generate 1 to 5 by 5-fold axis5[1]
  for i5 in range(1,6):
    angle=2.*pi*(i5-1)/5.

    icosamat[i5] = [0.,  lll**2+(mmm**2+nnn**2)*cos(angle),      lll*mmm*(1.-cos(angle))-nnn*sin(angle), nnn*lll*(1.-cos(angle))+mmm*sin(angle), 0.,
                         lll*mmm*(1.-cos(angle))+nnn*sin(angle), mmm**2+(lll**2+nnn**2)*cos(angle),      mmm*nnn*(1.-cos(angle))-lll*sin(angle), 0.,
                         nnn*lll*(1.-cos(angle))-mmm*sin(angle), mmm*nnn*(1.-cos(angle))+lll*sin(angle), nnn**2+(lll**2+mmm**2)*cos(angle),      0.]
                    
#  print icosamat[5][1:4]
#  print icosamat[5][5:8]
#  print icosamat[5][9:11]

#  Generating 3 fold symetries by z->x->y->z
# 1 to 27
  icosamat_axis3permu_zxyz(1,27)
#27 to 28
  icosamat_axis3permu_zxyz(27,28)

# 2 to 10
  icosamat_axis3permu_zxyz(2,10)
#10 to 12
  icosamat_axis3permu_zxyz(10,12)

# 3 to 9
  icosamat_axis3permu_zxyz(3,9)
# 9 to 11
  icosamat_axis3permu_zxyz(9,11)

# 4 to 8
  icosamat_axis3permu_zxyz(4,8)
# 8 to 59
  icosamat_axis3permu_zxyz(8,59)

# 5 to 26
  icosamat_axis3permu_zxyz(5,26)
# 26 to 41
  icosamat_axis3permu_zxyz(26,41)

# Generating Zaxis-2fold  x -> -x, y -> -y, z -> z

  icosamat_axis2_Zrot( 1,6)
  icosamat_axis2_Zrot(27,30)
  icosamat_axis2_Zrot(28,31)

  icosamat_axis2_Zrot( 2,18)
  icosamat_axis2_Zrot(10,19)
  icosamat_axis2_Zrot(12,21)

  icosamat_axis2_Zrot( 3,22)
  icosamat_axis2_Zrot( 9,17)
  icosamat_axis2_Zrot(11,20)

  icosamat_axis2_Zrot( 4,23)
  icosamat_axis2_Zrot( 8,16)
  icosamat_axis2_Zrot(59,35)

  icosamat_axis2_Zrot( 5, 7)
  icosamat_axis2_Zrot(26,15)
  icosamat_axis2_Zrot(41,34)

# Generating Xaxis-2fold  x -> x, y -> -y, z -> -z

  icosamat_axis2_Xrot( 1,48)
  icosamat_axis2_Xrot(27,55)
  icosamat_axis2_Xrot(28,32)

  icosamat_axis2_Xrot( 6,51)
  icosamat_axis2_Xrot(31,29)
  icosamat_axis2_Xrot(30,40)


  icosamat_axis2_Xrot( 2,53)
  icosamat_axis2_Xrot(10,54)
  icosamat_axis2_Xrot(12,52)

  icosamat_axis2_Xrot(18,43)
  icosamat_axis2_Xrot(19,44)
  icosamat_axis2_Xrot(21,45)


  icosamat_axis2_Xrot( 3,33)
  icosamat_axis2_Xrot( 9,47)
  icosamat_axis2_Xrot(11,56)

  icosamat_axis2_Xrot(22,42)
  icosamat_axis2_Xrot(17,50)
  icosamat_axis2_Xrot(20,39)


  icosamat_axis2_Xrot( 4,37)
  icosamat_axis2_Xrot( 8,57)
  icosamat_axis2_Xrot(59,25)

  icosamat_axis2_Xrot(23,58)
  icosamat_axis2_Xrot(16,36)
  icosamat_axis2_Xrot(35,14)


  icosamat_axis2_Xrot( 5,49)
  icosamat_axis2_Xrot(26,60)
  icosamat_axis2_Xrot(41,24)

  icosamat_axis2_Xrot( 7,46)
  icosamat_axis2_Xrot(15,38)
  icosamat_axis2_Xrot(34,13)

#  print icosamat[13][1:4]
#  print icosamat[13][5:8]
#  print icosamat[13][9:12]
#  for ico in range(61):
#    sum=0.
#    for index in range(1,12):
#     sum+=abs(icosamat[ico][index])
#    print ico, sum

  alpha,beta,gamma,xcen,ycen,zcen=0.,0.,0.,0.,0.,0.
  phaser_str=" "

  argc=len(argv)
  if argc > 1 :
    try:
      alpha,beta,gamma=atof(argv[1])*Rad,atof(argv[2])*Rad,atof(argv[3])*Rad
    except:
      howtouse()
    if argc > 4 :
      try:
         xcen,ycen,zcen=atof(argv[4]),atof(argv[5]),atof(argv[6])
      except:
         howtouse()
      if argc > 6:
         phaser_str=argv[7]
              
    
# Euler rotation from input

  rotn=[0,cos(alpha)*cos(beta)*cos(gamma) - sin(alpha)* sin(gamma), -cos(alpha)*cos(beta)*sin(gamma) - sin(alpha)*cos(gamma),cos(alpha)*sin(beta),
           sin(alpha)*cos(beta)*cos(gamma) + cos(alpha)* sin(gamma), -sin(alpha)*cos(beta)*sin(gamma) + cos(alpha)*cos(gamma),sin(alpha)*sin(beta),
                     -sin(beta)*cos(gamma)                         ,             sin(beta)*sin(gamma)                        ,cos(beta)            ]


#rotr is inverse matrix  of rotn
  rotr[0]                 = rotn[1]*rotn[5]*rotn[9] + rotn[4]*rotn[8]*rotn[3] + rotn[7]*rotn[2]*rotn[6] - rotn[1]*rotn[8]*rotn[6] - rotn[7]*rotn[5]*rotn[3] - rotn[4]*rotn[2]*rotn[9] 
  rotr[1],rotr[2],rotr[3] = rotn[5]*rotn[9] - rotn[6]*rotn[8], rotn[3]*rotn[8] - rotn[2]*rotn[9], rotn[2]*rotn[6] - rotn[3]*rotn[5]
  rotr[4],rotr[5],rotr[6] = rotn[6]*rotn[7] - rotn[4]*rotn[9], rotn[1]*rotn[9] - rotn[3]*rotn[7], rotn[3]*rotn[4] - rotn[1]*rotn[6]
  rotr[7],rotr[8],rotr[9] = rotn[4]*rotn[8] - rotn[5]*rotn[7], rotn[2]*rotn[7] - rotn[1]*rotn[8], rotn[1]*rotn[5] - rotn[2]*rotn[4]

#
  f_refmac_m=open("refmac_ncscon_mat",'w')
  f_refmac_e=open("refmac_ncscon_eul",'w')
  f_phaser=open("phaser_sol_eul",'w')  
  f_phaser.write("SOLU SET \n")
#  
  for ico in range(1,61):

         rot2[1] = icosamat[ico][1]*rotr[1] +  icosamat[ico][2]*rotr[4] + icosamat[ico][3]*rotr[7] 
         rot2[2] = icosamat[ico][1]*rotr[2] +  icosamat[ico][2]*rotr[5] + icosamat[ico][3]*rotr[8] 
         rot2[3] = icosamat[ico][1]*rotr[3] +  icosamat[ico][2]*rotr[6] + icosamat[ico][3]*rotr[9] 
         rot2[4] = icosamat[ico][5]*rotr[1] +  icosamat[ico][6]*rotr[4] + icosamat[ico][7]*rotr[7] 
         rot2[5] = icosamat[ico][5]*rotr[2] +  icosamat[ico][6]*rotr[5] + icosamat[ico][7]*rotr[8] 
         rot2[6] = icosamat[ico][5]*rotr[3] +  icosamat[ico][6]*rotr[6] + icosamat[ico][7]*rotr[9] 
         rot2[7] = icosamat[ico][9]*rotr[1] +  icosamat[ico][10]*rotr[4] + icosamat[ico][11]*rotr[7] 
         rot2[8] = icosamat[ico][9]*rotr[2] +  icosamat[ico][10]*rotr[5] + icosamat[ico][11]*rotr[8] 
         rot2[9] = icosamat[ico][9]*rotr[3] +  icosamat[ico][10]*rotr[6] + icosamat[ico][11]*rotr[9] 

         rot3[1] = rotn[1]*rot2[1] +  rotn[2]*rot2[4] + rotn[3]*rot2[7] 
         rot3[2] = rotn[1]*rot2[2] +  rotn[2]*rot2[5] + rotn[3]*rot2[8] 
         rot3[3] = rotn[1]*rot2[3] +  rotn[2]*rot2[6] + rotn[3]*rot2[9] 
         rot3[4] = rotn[4]*rot2[1] +  rotn[5]*rot2[4] + rotn[6]*rot2[7] 
         rot3[5] = rotn[4]*rot2[2] +  rotn[5]*rot2[5] + rotn[6]*rot2[8] 
         rot3[6] = rotn[4]*rot2[3] +  rotn[5]*rot2[6] + rotn[6]*rot2[9] 
         rot3[7] = rotn[7]*rot2[1] +  rotn[8]*rot2[4] + rotn[9]*rot2[7] 
         rot3[8] = rotn[7]*rot2[2] +  rotn[8]*rot2[5] + rotn[9]*rot2[8] 
         rot3[9] = rotn[7]*rot2[3] +  rotn[8]*rot2[6] + rotn[9]*rot2[9] 

         ncenx=xcen-(xcen*rot3[1]+ycen*rot3[2]+zcen*rot3[3])
         nceny=ycen-(xcen*rot3[4]+ycen*rot3[5]+zcen*rot3[6])
         ncenz=zcen-(xcen*rot3[7]+ycen*rot3[8]+zcen*rot3[9])

         for ii in range(9):
           if -0.0000001<rot3[ii]<0.000000001:
                 rot3[ii]=0.000000

         f_refmac_m.write("ncscons matrix - \n")
         f_refmac_m.write(" %13.5f  %13.5f  %13.5f -\n"%(rot3[1],rot3[2],rot3[3]))
         f_refmac_m.write(" %13.5f  %13.5f  %13.5f -\n"%(rot3[4],rot3[5],rot3[6]))
         f_refmac_m.write(" %13.5f  %13.5f  %13.5f -\n"%(rot3[7],rot3[8],rot3[9]))

         f_refmac_m.write(" %15.3f  %15.3f  %15.3f \n"%(ncenx,nceny,ncenz))

# matrix to Euler

         if -0.0001 < rot3[9] < 0.0001: #rot33==0 beta=90.00
             nbeta=90.000
             nalpha=atan2(rot3[6],rot3[3])/Rad 
             ngamma=atan2(rot3[8],-rot3[7])/Rad      
         elif rot3[9] < -0.9999: #rot33==-1 beta=180.0
             nbeta=180.000
             nalpha=atan2(-rot3[2],rot3[5])/Rad
             ngamma=0.0000
         elif 0.9999 < rot3[9]: #rot33 == 1 beta=0.000
             nbeta=0.00000
             nalpha=atan2(rot3[4],rot3[1])/Rad
             ngamma=0.0000
         else:
             betarad=acos(rot3[9])
             sb=sin(betarad)
             nbeta=betarad/Rad
             nalpha=atan2(rot3[6]/sb,rot3[3]/sb)/Rad
             ngamma=atan2(rot3[8]/sb,-rot3[7]/sb)/Rad

         f_refmac_e.write("ncscons euler %13.4f %13.4f %13.4f %15.3f %15.3f %15.3f \n"%(nalpha,nbeta,ngamma,ncenx,nceny,ncenz))         
         f_phaser.write("SOLU 6DIM ENSE mol1 EULER %12.3f %12.3f %12.3f   ORTH %12.3f %12.3f %12.3f %s\n"%(nalpha,nbeta,ngamma,ncenx,nceny,ncenz,phaser_str))         
 
  f_refmac_m.close()
  f_refmac_e.close()
  f_phaser.close()

           
