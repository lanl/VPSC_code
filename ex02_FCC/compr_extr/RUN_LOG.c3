
 **** INPUT FILE VPSC.IN FOR THIS RUN ****
1                          number of elements (nelem)                                                                             
1                          number of phases (nph)                                                                                 
1.0  0.0                   relative vol. fract. of phases (wph(i))                                                                
*INFORMATION ABOUT PHASE #1                                                                                                       
0   0   25                      grain shape contrl, fragmentn, crit aspect ratio                                                  
1.0  1.0  1.0                 initial ellipsoid ratios (dummy if ishape=4)                                                        
0.0  0.0  0.0                 init Eul ang ellips axes (dummy if ishape=3,4)                                                      
* name and path of texture file (filetext)                                                                                        
ex02_FCC\rand500.tex                                                                                                              
* name and path of single crystal file (filecrys)                                                                                 
ex02_FCC\fcc.sx                                                                                                                   
* name and path of grain shape file (dummy if ishape=0) (fileaxes)                                                                
dummy                                                                                                                             
* name and path of diffraction file (dummy if idiff=0)                                                                            
0                                                                                                                                 
dummy                                                                                                                             
*PRECISION SETTINGS FOR CONVERGENCE PROCEDURES (default values)                                                                   
0.001 0.001 0.001 0.001    errs,errd,errm,errso                                                                                   
100 100 25     itmax:   max # of iter, external, internal and SO loops                                                            
0  2  10  2    irsvar & jrsini,jrsfin,jrstep (dummy if irsvar=0)                                                                  
*INPUT/OUTPUT SETTINGS FOR THE RUN (default is zero)                                                                              
0              irecover:read grain states from POSTMORT.IN (1) or not (0)?                                                        
0              isave:   write grain states in POSTMORT.OUT at step 'isave'?                                                       
0              icubcomp:calculate fcc rolling components?                                                                         
0              nwrite (frequency of texture downloads)                                                                            
*MODELING CONDITIONS FOR THE RUN                                                                                                  
3               interaction (0:FC,1:affine,2:secant,3:neff=10,4:tangent,5:SO)                                                     
1  1  1        iupdate: update orient, grain shape, hardening                                                                     
0              nneigh (0 for no neighbors, 1 for pairs, etc.)                                                                     
0              iflu (0: don't calc, 1: calc fluctuations)                                                                         
*NUMBER OF PROCESSES (Lij const; Lij variable; PCYS ;LANKFORD; rigid rotatn)                                                      
1                                                                                                                                 
*IVGVAR AND PATH\NAME OF FILE FOR EACH PROCESS                                                                                    
0                                                                                                                                 
ex02_FCC\compres_extrusion.3                                                                                                      

  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  0
  EULER ANGLES   180.00    90.00   180.00
  EIGNVAL   1.0000                     1.0000                     1.0000
  EIGNVEC   1.0000   0.0000   0.0000  -0.0000  -0.0000  -1.0000   0.0000   1.0000   0.0000

 **** CRYSTAL DATA FILE ****
*Material: AUSTENITIC STEEL                                                                                                       
cubic           crysym                                                                                                            
   1.0   1.0   1.0   90.   90.   90.   unit cell axes and angles                                                                  
Elastic stiffness (single crystal [GPa]; scaled=0.85xINTERPOLATED)                                                                
 205.0   138.0   138.0   000.0   000.0   000.0                                                                                    
 138.0   205.0   138.0   000.0   000.0   000.0                                                                                    
 138.0   138.0   205.0   000.0   000.0   000.0                                                                                    
 000.0   000.0   000.0   126.0   000.0   000.0                                                                                    
 000.0   000.0   000.0   000.0   126.0   000.0                                                                                    
 000.0   000.0   000.0   000.0   000.0   126.0                                                                                    
*Thermal expansion coefficients (single crystal in crystal axis):                                                                 
  1.0e-6   1.0e-6  1.0e-6   0.0e0   0.0e0   0.0e0                                                                                 
*Info about slip & twinning modes in this file:                                                                                   
  2          nmodesx    (total # of modes listed in file)                                                                         
  1          nmodes     (# of modes to be used in the calculation)                                                                
  1          mode(i)    (label of the modes to be used)                                                                           
  <111>{110} SLIP                                                                                                                 
 1  12  1   0                       modex,nsmx,isensex,itwtypex                                                                   
   1  1  1        0  1 -1           slip (n & b)                                                                                  
   1  1  1        1  0 -1                                                                                                         
   1  1  1        1 -1  0                                                                                                         
  -1  1  1        0  1 -1                                                                                                         
  -1  1  1        1  0  1                                                                                                         
  -1  1  1        1  1  0                                                                                                         
  -1 -1  1        0  1  1                                                                                                         
  -1 -1  1        1  0  1                                                                                                         
  -1 -1  1        1 -1  0                                                                                                         
   1 -1  1        0  1  1                                                                                                         
   1 -1  1        1  0 -1                                                                                                         
   1 -1  1        1  1  0                                                                                                         
  <111>{112} TWIN                                                                                                                 
 2  12  0   2                        modex,nsmx,isensex,itwtypex                                                                  
 0.707                               twshx                                                                                        
   1  1  1       -2  1  1                                                                                                         
   1  1  1        1 -2  1                                                                                                         
   1  1  1        1  1 -2                                                                                                         
  -1  1  1        2  1  1                                                                                                         
  -1  1  1       -1 -2  1                                                                                                         
  -1  1  1       -1  1 -2                                                                                                         
  -1 -1  1        2 -1  1                                                                                                         
  -1 -1  1       -1  2  1                                                                                                         
  -1 -1  1       -1 -1 -2                                                                                                         
   1 -1  1       -2 -1  1                                                                                                         
   1 -1  1        1  2  1                                                                                                         
   1 -1  1        1 -1 -2                                                                                                         
*Constitutive law                                                                                                                 
   0      Voce=0, MTS=1                                                                                                           
   1      iratesens (0:rate insensitive, 1:rate sensitive)                                                                        
   50     grsze --> grain size only matters if HPfactor is non-zero                                                               
  <111>{110} SLIP -------------------------------------------                                                                     
 20                               nrsx                                                                                            
 1.0   0.0   0.0   0.0  0.        tau0x,tau1x,thet0,thet1, hpfac                                                                  
 1.0    1.0                          hlatex(1,im),im=1,nmodes                                                                     
  <111>{112} TWIN --------------------------------------                                                                          
 20                               nrsx                                                                                            
 1.0   0.0    0.   0.   0.        tau0x,tau1x,thet0,thet1, hpfac                                                                  
 1.0    1.0                       hlatex(1,im),im=1,nmodes                                                                        
 0     0.05   0.50                isectw, thres1,thres2                                                                           
**** END OF CRYSTAL DATA FILE ****


 CHECKING THAT CELCC*SELCC-ID6=0   0.2134069E-31
 *********** PHASE   1
 RANDOM PX BULK & POISSON MODULI     160.333       0.266
 RANDOM PX ELASTIC CTES C11, C12, C44     279.000     101.000      89.000

 INSIDE SUBROUTINE DATA_CRYSTAL
 SHEAR MODULUS FOR MODE  1 IN PHASE  1 IS      64.333
 N & B FOR MODE  1 IN PHASE  1
     0.577     0.577     0.577        0.000     0.707    -0.707
     0.577     0.577     0.577        0.707     0.000    -0.707
     0.577     0.577     0.577        0.707    -0.707     0.000
    -0.577     0.577     0.577        0.000     0.707    -0.707
    -0.577     0.577     0.577        0.707     0.000     0.707
    -0.577     0.577     0.577        0.707     0.707     0.000
    -0.577    -0.577     0.577        0.000     0.707     0.707
    -0.577    -0.577     0.577        0.707     0.000     0.707
    -0.577    -0.577     0.577        0.707    -0.707     0.000
     0.577    -0.577     0.577        0.000     0.707     0.707
     0.577    -0.577     0.577        0.707     0.000    -0.707
     0.577    -0.577     0.577        0.707     0.707     0.000

 **** CRYST TEXTURE (FIRST FEW LINES) ****
dummy                                                                                                                             
dummy                                                                                                                             
random texture generated by RANDTEXT.FOR (23/01/97)                                                                               
B   500                                                                                                                           
    102.74    119.56     33.65  1.0                                                                                               
    219.06     36.21     70.51  1.0                                                                                               
    166.66     28.59     45.80  1.0                                                                                               
    149.74     86.13     38.68  1.0                                                                                               
    143.88     72.50      9.75  1.0                                                                                               
     27.71      7.88    225.69  1.0                                                                                               
     10.26     81.59    324.82  1.0                                                                                               
    213.21    156.91    352.41  1.0                                                                                               
    313.75     60.06    131.59  1.0                                                                                               
    283.74    137.64    213.65  1.0                                                                                               
    329.59    158.18    298.33  1.0                                                                                               
     61.37    171.68    104.65  1.0                                                                                               
     44.90     68.33     40.28  1.0                                                                                               
     29.40    114.77    163.14  1.0                                                                                               
     56.39    115.89    341.90  1.0                                                                                               
    271.21     76.66    171.67  1.0                                                                                               
    .........................
 **** END OF CRYST TEXTURE DATA FILE ****

  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  1
  EULER ANGLES   180.00    90.00   180.00
  EIGNVAL   1.0000                     1.0000                     1.0000
  EIGNVEC   1.0000   0.0000   0.0000  -0.0000  -0.0000  -1.0000   0.0000   1.0000   0.0000

 UB THERMAL TENSOR (VOIGT NOTATION)
  0.1000E-05  0.1000E-05  0.1000E-05  0.1482E-23 -0.1622E-23 -0.3991E-23

 UB ELASTIC STIFFNESS (VOIGT NOTATION)
  0.2795E+03  0.9994E+02  0.1016E+03  0.1676E+00  0.1479E+01  0.6049E+00
  0.9994E+02  0.2801E+03  0.1009E+03 -0.1180E+01  0.1349E+00 -0.7842E+00
  0.1016E+03  0.1009E+03  0.2785E+03  0.1012E+01 -0.1614E+01  0.1793E+00
  0.1676E+00 -0.1180E+01  0.1012E+01  0.8894E+02  0.1793E+00  0.1349E+00
  0.1479E+01  0.1349E+00 -0.1614E+01  0.1793E+00  0.8959E+02  0.1676E+00
  0.6049E+00 -0.7842E+00  0.1793E+00  0.1349E+00  0.1676E+00  0.8794E+02

 LB THERMAL TENSOR (VOIGT NOTATION)
  0.1000E-05  0.1000E-05  0.1000E-05  0.2658E-23 -0.1216E-23 -0.1258E-23

 LB ELASTIC STIFFNESS (VOIGT NOTATION)
  0.2406E+03  0.1195E+03  0.1209E+03  0.1580E+00  0.1267E+01  0.5168E+00
  0.1195E+03  0.2412E+03  0.1203E+03 -0.1012E+01  0.1034E+00 -0.6672E+00
  0.1209E+03  0.1203E+03  0.2398E+03  0.8542E+00 -0.1370E+01  0.1503E+00
  0.1580E+00 -0.1012E+01  0.8542E+00  0.5984E+02  0.1446E+00  0.1204E+00
  0.1267E+01  0.1034E+00 -0.1370E+01  0.1446E+00  0.6040E+02  0.1448E+00
  0.5168E+00 -0.6672E+00  0.1503E+00  0.1204E+00  0.1448E+00  0.5899E+02


*** LOAD CONDITIONS FOR THIS RUN
  50   3   0.02    298.   298.   nsteps  ictrl  eqincr  temp_i  temp_f                                                            
* boundary conditions                                                                                                             
    1       1       1           iudot    |    flag for vel.grad.                                                                  
    1       1       1                    |    (0:unknown-1:known)                                                                 
    1       1       1                    |                                                                                        
                                         |                                                                                        
    0.5     0.      0.          udot     |    vel.grad                                                                            
    0.      0.5     0.                   |                                                                                        
    0.      0.     -1.0                  |                                                                                        
                                         |                                                                                        
    0       0       0           iscau    |    flag for Cauchy                                                                     
            0       0                    |                                                                                        
                    0                    |                                                                                        
                                         |                                                                                        
    0.      0.      0.          scauchy  |    Cauchy stress                                                                       
            0.      0.                   |                                                                                        
                    0.                   @                                                                                        
                                                                                                                                 

*******   STEP   1      *** N_EFF=10 CALCULATION               
   10   0.927E-04   0.390E-03   0.650E-03
  INSIDE UPD_FIJ - PHASE #  0
  LIJBAR     0.5000    0.0000    0.0000
  LIJBAR     0.0000    0.5000    0.0000
  LIJBAR     0.0000    0.0000   -1.0000
  FIJ        1.0100    0.0000    0.0000
  FIJ        0.0000    1.0100    0.0000
  FIJ        0.0000    0.0000    0.9800
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  0
  EULER ANGLES     0.00     0.00     0.00
  EIGNVAL   1.0100                     1.0100                     0.9800
  EIGNVEC   1.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000   1.0000
  INSIDE UPD_FIJ - PHASE #  1
  LIJPH      0.5003    0.0000    0.0000
  LIJPH      0.0000    0.5003    0.0000
  LIJPH      0.0000    0.0000   -1.0006
  FIJ        1.0100    0.0000    0.0000
  FIJ        0.0000    1.0100    0.0000
  FIJ        0.0000    0.0000    0.9800
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  1
  EULER ANGLES   -30.55     0.00     0.00
  EIGNVAL   1.0100                     1.0100                     0.9800
  EIGNVEC   0.8612  -0.5082  -0.0000   0.5082   0.8612   0.0000  -0.0000  -0.0000   1.0000

*******   STEP   2      *** N_EFF=10 CALCULATION               
    6   0.125E-03   0.579E-03   0.176E-03
  INSIDE UPD_FIJ - PHASE #  0
  LIJBAR     0.5000    0.0000    0.0000
  LIJBAR     0.0000    0.5000    0.0000
  LIJBAR     0.0000    0.0000   -1.0000
  FIJ        1.0201    0.0000    0.0000
  FIJ        0.0000    1.0201    0.0000
  FIJ        0.0000    0.0000    0.9604
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  0
  EULER ANGLES     0.00     0.00     0.00
  EIGNVAL   1.0201                     1.0201                     0.9604
  EIGNVEC   1.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000   1.0000
  INSIDE UPD_FIJ - PHASE #  1
  LIJPH      0.5001   -0.0000    0.0000
  LIJPH     -0.0000    0.5001    0.0000
  LIJPH      0.0000    0.0000   -1.0002
  FIJ        1.0201    0.0000    0.0000
  FIJ        0.0000    1.0201    0.0000
  FIJ        0.0000    0.0000    0.9604
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  1
  EULER ANGLES    -8.24     0.00     0.00
  EIGNVAL   1.0201                     1.0201                     0.9604
  EIGNVEC   0.9897  -0.1434  -0.0000   0.1434   0.9897   0.0000  -0.0000  -0.0000   1.0000

*******   STEP   3      *** N_EFF=10 CALCULATION               
    6   0.147E-03   0.666E-03   0.265E-03
  INSIDE UPD_FIJ - PHASE #  0
  LIJBAR     0.5000    0.0000    0.0000
  LIJBAR     0.0000    0.5000    0.0000
  LIJBAR     0.0000    0.0000   -1.0000
  FIJ        1.0303    0.0000    0.0000
  FIJ        0.0000    1.0303    0.0000
  FIJ        0.0000    0.0000    0.9412
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  0
  EULER ANGLES     0.00     0.00     0.00
  EIGNVAL   1.0303                     1.0303                     0.9412
  EIGNVEC   1.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000   1.0000
  INSIDE UPD_FIJ - PHASE #  1
  LIJPH      0.5001   -0.0000    0.0000
  LIJPH     -0.0000    0.5001    0.0000
  LIJPH      0.0000    0.0000   -1.0003
  FIJ        1.0303    0.0000    0.0000
  FIJ        0.0000    1.0303    0.0000
  FIJ        0.0000    0.0000    0.9412
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  1
  EULER ANGLES    -2.42     0.00     0.00
  EIGNVAL   1.0303                     1.0303                     0.9412
  EIGNVEC   0.9991  -0.0423   0.0000   0.0423   0.9991   0.0000  -0.0000  -0.0000   1.0000

*******   STEP   4      *** N_EFF=10 CALCULATION               
    5   0.102E-03   0.982E-03   0.557E-03
  INSIDE UPD_FIJ - PHASE #  0
  LIJBAR     0.5000    0.0000    0.0000
  LIJBAR     0.0000    0.5000    0.0000
  LIJBAR     0.0000    0.0000   -1.0000
  FIJ        1.0406    0.0000    0.0000
  FIJ        0.0000    1.0406    0.0000
  FIJ        0.0000    0.0000    0.9224
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  0
  EULER ANGLES     0.00     0.00     0.00
  EIGNVAL   1.0406                     1.0406                     0.9224
  EIGNVEC   1.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000   1.0000
  INSIDE UPD_FIJ - PHASE #  1
  LIJPH      0.5003   -0.0000   -0.0000
  LIJPH     -0.0000    0.5003    0.0000
  LIJPH     -0.0000    0.0000   -1.0006
  FIJ        1.0406    0.0000    0.0000
  FIJ        0.0000    1.0406    0.0000
  FIJ        0.0000    0.0000    0.9223
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  1
  EULER ANGLES    -2.43     0.00     0.00
  EIGNVAL   1.0406                     1.0406                     0.9223
  EIGNVEC   0.9991  -0.0425   0.0000   0.0425   0.9991   0.0000  -0.0000  -0.0000   1.0000

*******   STEP   5      *** N_EFF=10 CALCULATION               
    5   0.119E-03   0.911E-03   0.543E-03
  INSIDE UPD_FIJ - PHASE #  0
  LIJBAR     0.5000    0.0000    0.0000
  LIJBAR     0.0000    0.5000    0.0000
  LIJBAR     0.0000    0.0000   -1.0000
  FIJ        1.0510    0.0000    0.0000
  FIJ        0.0000    1.0510    0.0000
  FIJ        0.0000    0.0000    0.9039
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  0
  EULER ANGLES     0.00     0.00     0.00
  EIGNVAL   1.0510                     1.0510                     0.9039
  EIGNVEC   1.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000   1.0000
  INSIDE UPD_FIJ - PHASE #  1
  LIJPH      0.5003   -0.0000   -0.0000
  LIJPH     -0.0000    0.5003    0.0000
  LIJPH     -0.0000    0.0000   -1.0005
  FIJ        1.0510    0.0000    0.0000
  FIJ        0.0000    1.0510    0.0000
  FIJ        0.0000    0.0000    0.9039
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  1
  EULER ANGLES    -1.82     0.00     0.00
  EIGNVAL   1.0510                     1.0510                     0.9039
  EIGNVEC   0.9995  -0.0318  -0.0000   0.0318   0.9995   0.0000  -0.0000  -0.0000   1.0000

*******   STEP   6      *** N_EFF=10 CALCULATION               
    5   0.911E-04   0.776E-03   0.531E-03
  INSIDE UPD_FIJ - PHASE #  0
  LIJBAR     0.5000    0.0000    0.0000
  LIJBAR     0.0000    0.5000    0.0000
  LIJBAR     0.0000    0.0000   -1.0000
  FIJ        1.0615    0.0000    0.0000
  FIJ        0.0000    1.0615    0.0000
  FIJ        0.0000    0.0000    0.8858
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  0
  EULER ANGLES     0.00     0.00     0.00
  EIGNVAL   1.0615                     1.0615                     0.8858
  EIGNVEC   1.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000   1.0000
  INSIDE UPD_FIJ - PHASE #  1
  LIJPH      0.5003   -0.0000   -0.0000
  LIJPH     -0.0000    0.5003    0.0000
  LIJPH     -0.0000    0.0000   -1.0005
  FIJ        1.0615    0.0000   -0.0000
  FIJ        0.0000    1.0615    0.0000
  FIJ       -0.0000    0.0000    0.8858
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  1
  EULER ANGLES    -1.16     0.00     0.00
  EIGNVAL   1.0615                     1.0615                     0.8858
  EIGNVEC   0.9998  -0.0203  -0.0000   0.0203   0.9998   0.0000   0.0000  -0.0000   1.0000

*******   STEP   7      *** N_EFF=10 CALCULATION               
    4   0.365E-03   0.802E-03   0.484E-03
  INSIDE UPD_FIJ - PHASE #  0
  LIJBAR     0.5000    0.0000    0.0000
  LIJBAR     0.0000    0.5000    0.0000
  LIJBAR     0.0000    0.0000   -1.0000
  FIJ        1.0721    0.0000    0.0000
  FIJ        0.0000    1.0721    0.0000
  FIJ        0.0000    0.0000    0.8681
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  0
  EULER ANGLES     0.00     0.00     0.00
  EIGNVAL   1.0721                     1.0721                     0.8681
  EIGNVEC   1.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000   1.0000
  INSIDE UPD_FIJ - PHASE #  1
  LIJPH      0.5003   -0.0000   -0.0000
  LIJPH     -0.0000    0.5002    0.0000
  LIJPH     -0.0000    0.0000   -1.0005
  FIJ        1.0722   -0.0000   -0.0000
  FIJ       -0.0000    1.0722    0.0000
  FIJ       -0.0000    0.0000    0.8681
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  1
  EULER ANGLES    11.64     0.00     0.00
  EIGNVAL   1.0722                     1.0722                     0.8681
  EIGNVEC   0.9794   0.2018   0.0000  -0.2018   0.9794   0.0000   0.0000  -0.0000   1.0000

*******   STEP   8      *** N_EFF=10 CALCULATION               
    4   0.649E-03   0.827E-03   0.632E-03
  INSIDE UPD_FIJ - PHASE #  0
  LIJBAR     0.5000    0.0000    0.0000
  LIJBAR     0.0000    0.5000    0.0000
  LIJBAR     0.0000    0.0000   -1.0000
  FIJ        1.0829    0.0000    0.0000
  FIJ        0.0000    1.0829    0.0000
  FIJ        0.0000    0.0000    0.8508
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  0
  EULER ANGLES     0.00     0.00     0.00
  EIGNVAL   1.0829                     1.0829                     0.8508
  EIGNVEC   1.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000   1.0000
  INSIDE UPD_FIJ - PHASE #  1
  LIJPH      0.5003    0.0000   -0.0000
  LIJPH      0.0000    0.5003    0.0000
  LIJPH     -0.0000    0.0000   -1.0006
  FIJ        1.0829   -0.0000   -0.0000
  FIJ       -0.0000    1.0829    0.0000
  FIJ       -0.0000    0.0000    0.8507
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  1
  EULER ANGLES     2.08     0.00     0.00
  EIGNVAL   1.0829                     1.0829                     0.8507
  EIGNVEC   0.9993   0.0363  -0.0000  -0.0363   0.9993   0.0000   0.0000  -0.0000   1.0000

*******   STEP   9      *** N_EFF=10 CALCULATION               
    4   0.604E-03   0.736E-03   0.597E-03
  INSIDE UPD_FIJ - PHASE #  0
  LIJBAR     0.5000    0.0000    0.0000
  LIJBAR     0.0000    0.5000    0.0000
  LIJBAR     0.0000    0.0000   -1.0000
  FIJ        1.0937    0.0000    0.0000
  FIJ        0.0000    1.0937    0.0000
  FIJ        0.0000    0.0000    0.8337
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  0
  EULER ANGLES     0.00     0.00     0.00
  EIGNVAL   1.0937                     1.0937                     0.8337
  EIGNVEC   1.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000   1.0000
  INSIDE UPD_FIJ - PHASE #  1
  LIJPH      0.5003    0.0000   -0.0000
  LIJPH      0.0000    0.5003    0.0000
  LIJPH     -0.0000    0.0000   -1.0006
  FIJ        1.0937   -0.0000   -0.0000
  FIJ       -0.0000    1.0937    0.0000
  FIJ       -0.0000    0.0000    0.8337
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  1
  EULER ANGLES     0.10     0.00     0.00
  EIGNVAL   1.0937                     1.0937                     0.8337
  EIGNVEC   1.0000   0.0017  -0.0000  -0.0017   1.0000   0.0000   0.0000  -0.0000   1.0000

*******   STEP  10      *** N_EFF=10 CALCULATION               
    4   0.954E-03   0.751E-03   0.933E-03
  INSIDE UPD_FIJ - PHASE #  0
  LIJBAR     0.5000    0.0000    0.0000
  LIJBAR     0.0000    0.5000    0.0000
  LIJBAR     0.0000    0.0000   -1.0000
  FIJ        1.1046    0.0000    0.0000
  FIJ        0.0000    1.1046    0.0000
  FIJ        0.0000    0.0000    0.8171
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  0
  EULER ANGLES     0.00     0.00     0.00
  EIGNVAL   1.1046                     1.1046                     0.8171
  EIGNVEC   1.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000   1.0000
  INSIDE UPD_FIJ - PHASE #  1
  LIJPH      0.5004    0.0000   -0.0000
  LIJPH      0.0000    0.5005   -0.0000
  LIJPH     -0.0000   -0.0000   -1.0009
  FIJ        1.1047    0.0000   -0.0000
  FIJ        0.0000    1.1047    0.0000
  FIJ       -0.0000    0.0000    0.8170
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  1
  EULER ANGLES    -0.72     0.00     0.00
  EIGNVAL   1.1047                     1.1047                     0.8170
  EIGNVEC   0.9999  -0.0125  -0.0000   0.0125   0.9999   0.0000   0.0000  -0.0000   1.0000

*******   STEP  11      *** N_EFF=10 CALCULATION               
    5   0.199E-03   0.875E-03   0.487E-03
  INSIDE UPD_FIJ - PHASE #  0
  LIJBAR     0.5000    0.0000    0.0000
  LIJBAR     0.0000    0.5000    0.0000
  LIJBAR     0.0000    0.0000   -1.0000
  FIJ        1.1157    0.0000    0.0000
  FIJ        0.0000    1.1157    0.0000
  FIJ        0.0000    0.0000    0.8007
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  0
  EULER ANGLES     0.00     0.00     0.00
  EIGNVAL   1.1157                     1.1157                     0.8007
  EIGNVEC   1.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000   1.0000
  INSIDE UPD_FIJ - PHASE #  1
  LIJPH      0.5002   -0.0000   -0.0000
  LIJPH     -0.0000    0.5002   -0.0000
  LIJPH     -0.0000   -0.0000   -1.0005
  FIJ        1.1157   -0.0000   -0.0000
  FIJ       -0.0000    1.1157    0.0000
  FIJ       -0.0000    0.0000    0.8006
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  1
  EULER ANGLES     0.72     0.00     0.00
  EIGNVAL   1.1157                     1.1157                     0.8006
  EIGNVEC   0.9999   0.0125  -0.0000  -0.0125   0.9999   0.0000   0.0000  -0.0000   1.0000

*******   STEP  12      *** N_EFF=10 CALCULATION               
    4   0.924E-03   0.924E-03   0.982E-03
  INSIDE UPD_FIJ - PHASE #  0
  LIJBAR     0.5000    0.0000    0.0000
  LIJBAR     0.0000    0.5000    0.0000
  LIJBAR     0.0000    0.0000   -1.0000
  FIJ        1.1268    0.0000    0.0000
  FIJ        0.0000    1.1268    0.0000
  FIJ        0.0000    0.0000    0.7847
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  0
  EULER ANGLES     0.00     0.00     0.00
  EIGNVAL   1.1268                     1.1268                     0.7847
  EIGNVEC   1.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000   1.0000
  INSIDE UPD_FIJ - PHASE #  1
  LIJPH      0.5005    0.0000   -0.0001
  LIJPH      0.0000    0.5005   -0.0000
  LIJPH     -0.0001   -0.0000   -1.0010
  FIJ        1.1269   -0.0000   -0.0000
  FIJ       -0.0000    1.1269    0.0000
  FIJ       -0.0000    0.0000    0.7846
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  1
  EULER ANGLES     0.41     0.00     0.00
  EIGNVAL   1.1269                     1.1269                     0.7846
  EIGNVEC   1.0000   0.0072  -0.0000  -0.0072   1.0000   0.0000   0.0000  -0.0000   1.0000

*******   STEP  13      *** N_EFF=10 CALCULATION               
    6   0.278E-03   0.709E-03   0.219E-03
  INSIDE UPD_FIJ - PHASE #  0
  LIJBAR     0.5000    0.0000    0.0000
  LIJBAR     0.0000    0.5000    0.0000
  LIJBAR     0.0000    0.0000   -1.0000
  FIJ        1.1381    0.0000    0.0000
  FIJ        0.0000    1.1381    0.0000
  FIJ        0.0000    0.0000    0.7690
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  0
  EULER ANGLES     0.00     0.00     0.00
  EIGNVAL   1.1381                     1.1381                     0.7690
  EIGNVEC   1.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000   1.0000
  INSIDE UPD_FIJ - PHASE #  1
  LIJPH      0.5001    0.0000   -0.0000
  LIJPH      0.0000    0.5001   -0.0000
  LIJPH     -0.0000   -0.0000   -1.0002
  FIJ        1.1382    0.0000   -0.0000
  FIJ        0.0000    1.1382    0.0000
  FIJ       -0.0000    0.0000    0.7689
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  1
  EULER ANGLES    -1.51     0.00     0.00
  EIGNVAL   1.1382                     1.1382                     0.7689
  EIGNVEC   0.9997  -0.0263  -0.0000   0.0263   0.9997   0.0000   0.0000  -0.0000   1.0000

*******   STEP  14      *** N_EFF=10 CALCULATION               
    6   0.206E-03   0.488E-03   0.533E-04
  INSIDE UPD_FIJ - PHASE #  0
  LIJBAR     0.5000    0.0000    0.0000
  LIJBAR     0.0000    0.5000    0.0000
  LIJBAR     0.0000    0.0000   -1.0000
  FIJ        1.1495    0.0000    0.0000
  FIJ        0.0000    1.1495    0.0000
  FIJ        0.0000    0.0000    0.7536
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  0
  EULER ANGLES     0.00     0.00     0.00
  EIGNVAL   1.1495                     1.1495                     0.7536
  EIGNVEC   1.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000   1.0000
  INSIDE UPD_FIJ - PHASE #  1
  LIJPH      0.5000    0.0000   -0.0000
  LIJPH      0.0000    0.5000   -0.0000
  LIJPH     -0.0000   -0.0000   -1.0000
  FIJ        1.1496    0.0000   -0.0000
  FIJ        0.0000    1.1496    0.0000
  FIJ       -0.0000    0.0000    0.7535
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  1
  EULER ANGLES    -3.19     0.00     0.00
  EIGNVAL   1.1496                     1.1496                     0.7535
  EIGNVEC   0.9984  -0.0557  -0.0000   0.0557   0.9984   0.0000   0.0000  -0.0000   1.0000

*******   STEP  15      *** N_EFF=10 CALCULATION               
    6   0.189E-03   0.437E-03   0.208E-03
  INSIDE UPD_FIJ - PHASE #  0
  LIJBAR     0.5000    0.0000    0.0000
  LIJBAR     0.0000    0.5000    0.0000
  LIJBAR     0.0000    0.0000   -1.0000
  FIJ        1.1610    0.0000    0.0000
  FIJ        0.0000    1.1610    0.0000
  FIJ        0.0000    0.0000    0.7386
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  0
  EULER ANGLES     0.00     0.00     0.00
  EIGNVAL   1.1610                     1.1610                     0.7386
  EIGNVEC   1.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000   1.0000
  INSIDE UPD_FIJ - PHASE #  1
  LIJPH      0.4999    0.0000   -0.0000
  LIJPH      0.0000    0.4999    0.0000
  LIJPH     -0.0000    0.0000   -0.9998
  FIJ        1.1610    0.0000   -0.0000
  FIJ        0.0000    1.1611    0.0000
  FIJ       -0.0000    0.0000    0.7385
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  1
  EULER ANGLES    -4.55     0.00     0.00
  EIGNVAL   1.1610                     1.1611                     0.7385
  EIGNVEC   0.9968  -0.0794  -0.0000   0.0794   0.9968   0.0000   0.0000  -0.0000   1.0000

*******   STEP  16      *** N_EFF=10 CALCULATION               
    6   0.186E-03   0.434E-03   0.331E-03
  INSIDE UPD_FIJ - PHASE #  0
  LIJBAR     0.5000    0.0000    0.0000
  LIJBAR     0.0000    0.5000    0.0000
  LIJBAR     0.0000    0.0000   -1.0000
  FIJ        1.1726    0.0000    0.0000
  FIJ        0.0000    1.1726    0.0000
  FIJ        0.0000    0.0000    0.7238
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  0
  EULER ANGLES     0.00     0.00     0.00
  EIGNVAL   1.1726                     1.1726                     0.7238
  EIGNVEC   1.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000   1.0000
  INSIDE UPD_FIJ - PHASE #  1
  LIJPH      0.4998    0.0000   -0.0000
  LIJPH      0.0000    0.4998    0.0000
  LIJPH     -0.0000    0.0000   -0.9997
  FIJ        1.1726    0.0000   -0.0000
  FIJ        0.0000    1.1727    0.0000
  FIJ       -0.0000    0.0000    0.7237
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  1
  EULER ANGLES    -5.64     0.00     0.00
  EIGNVAL   1.1726                     1.1727                     0.7237
  EIGNVEC   0.9952  -0.0982  -0.0000   0.0982   0.9952   0.0000   0.0000  -0.0000   1.0000

*******   STEP  17      *** N_EFF=10 CALCULATION               
    6   0.192E-03   0.445E-03   0.424E-03
  INSIDE UPD_FIJ - PHASE #  0
  LIJBAR     0.5000    0.0000    0.0000
  LIJBAR     0.0000    0.5000    0.0000
  LIJBAR     0.0000    0.0000   -1.0000
  FIJ        1.1843    0.0000    0.0000
  FIJ        0.0000    1.1843    0.0000
  FIJ        0.0000    0.0000    0.7093
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  0
  EULER ANGLES     0.00     0.00     0.00
  EIGNVAL   1.1843                     1.1843                     0.7093
  EIGNVEC   1.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000   1.0000
  INSIDE UPD_FIJ - PHASE #  1
  LIJPH      0.4998    0.0000   -0.0000
  LIJPH      0.0000    0.4998    0.0000
  LIJPH     -0.0000    0.0000   -0.9996
  FIJ        1.1844    0.0000   -0.0000
  FIJ        0.0000    1.1844    0.0000
  FIJ       -0.0000    0.0000    0.7092
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  1
  EULER ANGLES    -6.53     0.00     0.00
  EIGNVAL   1.1844                     1.1844                     0.7092
  EIGNVEC   0.9935  -0.1138  -0.0000   0.1138   0.9935   0.0000   0.0000  -0.0000   1.0000

*******   STEP  18      *** N_EFF=10 CALCULATION               
    6   0.207E-03   0.445E-03   0.491E-03
  INSIDE UPD_FIJ - PHASE #  0
  LIJBAR     0.5000    0.0000    0.0000
  LIJBAR     0.0000    0.5000    0.0000
  LIJBAR     0.0000    0.0000   -1.0000
  FIJ        1.1961    0.0000    0.0000
  FIJ        0.0000    1.1961    0.0000
  FIJ        0.0000    0.0000    0.6951
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  0
  EULER ANGLES     0.00     0.00     0.00
  EIGNVAL   1.1961                     1.1961                     0.6951
  EIGNVEC   1.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000   1.0000
  INSIDE UPD_FIJ - PHASE #  1
  LIJPH      0.4997    0.0000   -0.0000
  LIJPH      0.0000    0.4998    0.0000
  LIJPH     -0.0000    0.0000   -0.9995
  FIJ        1.1962    0.0000   -0.0000
  FIJ        0.0000    1.1962    0.0000
  FIJ       -0.0000    0.0000    0.6951
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  1
  EULER ANGLES    -7.29     0.00     0.00
  EIGNVAL   1.1962                     1.1962                     0.6951
  EIGNVEC   0.9919  -0.1269  -0.0000   0.1269   0.9919   0.0000   0.0000  -0.0000   1.0000

*******   STEP  19      *** N_EFF=10 CALCULATION               
    6   0.228E-03   0.434E-03   0.534E-03
  INSIDE UPD_FIJ - PHASE #  0
  LIJBAR     0.5000    0.0000    0.0000
  LIJBAR     0.0000    0.5000    0.0000
  LIJBAR     0.0000    0.0000   -1.0000
  FIJ        1.2081    0.0000    0.0000
  FIJ        0.0000    1.2081    0.0000
  FIJ        0.0000    0.0000    0.6812
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  0
  EULER ANGLES     0.00     0.00     0.00
  EIGNVAL   1.2081                     1.2081                     0.6812
  EIGNVEC   1.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000   1.0000
  INSIDE UPD_FIJ - PHASE #  1
  LIJPH      0.4997    0.0000    0.0000
  LIJPH      0.0000    0.4997    0.0000
  LIJPH      0.0000    0.0000   -0.9995
  FIJ        1.2082    0.0000   -0.0000
  FIJ        0.0000    1.2082    0.0000
  FIJ       -0.0000    0.0000    0.6812
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  1
  EULER ANGLES    -7.95     0.00     0.00
  EIGNVAL   1.2082                     1.2082                     0.6812
  EIGNVEC   0.9904  -0.1383  -0.0000   0.1383   0.9904   0.0000   0.0000  -0.0000   1.0000

*******   STEP  20      *** N_EFF=10 CALCULATION               
    6   0.167E-03   0.419E-03   0.335E-03
  INSIDE UPD_FIJ - PHASE #  0
  LIJBAR     0.5000    0.0000    0.0000
  LIJBAR     0.0000    0.5000    0.0000
  LIJBAR     0.0000    0.0000   -1.0000
  FIJ        1.2202    0.0000    0.0000
  FIJ        0.0000    1.2202    0.0000
  FIJ        0.0000    0.0000    0.6676
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  0
  EULER ANGLES     0.00     0.00     0.00
  EIGNVAL   1.2202                     1.2202                     0.6676
  EIGNVEC   1.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000   1.0000
  INSIDE UPD_FIJ - PHASE #  1
  LIJPH      0.4998    0.0000   -0.0000
  LIJPH      0.0000    0.4998    0.0000
  LIJPH     -0.0000    0.0000   -0.9997
  FIJ        1.2202    0.0000   -0.0000
  FIJ        0.0000    1.2203    0.0000
  FIJ       -0.0000    0.0000    0.6675
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  1
  EULER ANGLES    -8.93     0.00     0.00
  EIGNVAL   1.2202                     1.2203                     0.6675
  EIGNVEC   0.9879  -0.1553  -0.0000   0.1553   0.9879   0.0000   0.0000  -0.0000   1.0000

*******   STEP  21      *** N_EFF=10 CALCULATION               
    6   0.185E-03   0.407E-03   0.312E-03
  INSIDE UPD_FIJ - PHASE #  0
  LIJBAR     0.5000    0.0000    0.0000
  LIJBAR     0.0000    0.5000    0.0000
  LIJBAR     0.0000    0.0000   -1.0000
  FIJ        1.2324    0.0000    0.0000
  FIJ        0.0000    1.2324    0.0000
  FIJ        0.0000    0.0000    0.6543
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  0
  EULER ANGLES     0.00     0.00     0.00
  EIGNVAL   1.2324                     1.2324                     0.6543
  EIGNVEC   1.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000   1.0000
  INSIDE UPD_FIJ - PHASE #  1
  LIJPH      0.4998    0.0000   -0.0000
  LIJPH      0.0000    0.4998    0.0000
  LIJPH     -0.0000    0.0000   -0.9997
  FIJ        1.2324    0.0000   -0.0000
  FIJ        0.0000    1.2325    0.0000
  FIJ       -0.0000    0.0000    0.6542
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  1
  EULER ANGLES    -9.90     0.00     0.00
  EIGNVAL   1.2324                     1.2325                     0.6542
  EIGNVEC   0.9851  -0.1720  -0.0000   0.1720   0.9851   0.0000   0.0000  -0.0000   1.0000

*******   STEP  22      *** N_EFF=10 CALCULATION               
    6   0.206E-03   0.395E-03   0.305E-03
  INSIDE UPD_FIJ - PHASE #  0
  LIJBAR     0.5000    0.0000    0.0000
  LIJBAR     0.0000    0.5000    0.0000
  LIJBAR     0.0000    0.0000   -1.0000
  FIJ        1.2447    0.0000    0.0000
  FIJ        0.0000    1.2447    0.0000
  FIJ        0.0000    0.0000    0.6412
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  0
  EULER ANGLES     0.00     0.00     0.00
  EIGNVAL   1.2447                     1.2447                     0.6412
  EIGNVEC   1.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000   1.0000
  INSIDE UPD_FIJ - PHASE #  1
  LIJPH      0.4999    0.0000    0.0000
  LIJPH      0.0000    0.4998    0.0000
  LIJPH      0.0000    0.0000   -0.9997
  FIJ        1.2448    0.0000   -0.0000
  FIJ        0.0000    1.2448    0.0000
  FIJ       -0.0000    0.0000    0.6411
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  1
  EULER ANGLES   -10.89     0.00     0.00
  EIGNVAL   1.2448                     1.2448                     0.6411
  EIGNVEC   0.9820  -0.1889  -0.0000   0.1889   0.9820   0.0000   0.0000  -0.0000   1.0000

*******   STEP  23      *** N_EFF=10 CALCULATION               
    6   0.229E-03   0.389E-03   0.294E-03
  INSIDE UPD_FIJ - PHASE #  0
  LIJBAR     0.5000    0.0000    0.0000
  LIJBAR     0.0000    0.5000    0.0000
  LIJBAR     0.0000    0.0000   -1.0000
  FIJ        1.2572    0.0000    0.0000
  FIJ        0.0000    1.2572    0.0000
  FIJ        0.0000    0.0000    0.6283
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  0
  EULER ANGLES     0.00     0.00     0.00
  EIGNVAL   1.2572                     1.2572                     0.6283
  EIGNVEC   1.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000   1.0000
  INSIDE UPD_FIJ - PHASE #  1
  LIJPH      0.4999    0.0000    0.0000
  LIJPH      0.0000    0.4998    0.0000
  LIJPH      0.0000    0.0000   -0.9997
  FIJ        1.2572    0.0000   -0.0000
  FIJ        0.0000    1.2572    0.0000
  FIJ       -0.0000    0.0000    0.6283
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  1
  EULER ANGLES   -11.89     0.00     0.00
  EIGNVAL   1.2572                     1.2572                     0.6283
  EIGNVEC   0.9786  -0.2060  -0.0000   0.2060   0.9786   0.0000   0.0000  -0.0000   1.0000

*******   STEP  24      *** N_EFF=10 CALCULATION               
    6   0.252E-03   0.392E-03   0.282E-03
  INSIDE UPD_FIJ - PHASE #  0
  LIJBAR     0.5000    0.0000    0.0000
  LIJBAR     0.0000    0.5000    0.0000
  LIJBAR     0.0000    0.0000   -1.0000
  FIJ        1.2697    0.0000    0.0000
  FIJ        0.0000    1.2697    0.0000
  FIJ        0.0000    0.0000    0.6158
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  0
  EULER ANGLES     0.00     0.00     0.00
  EIGNVAL   1.2697                     1.2697                     0.6158
  EIGNVEC   1.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000   1.0000
  INSIDE UPD_FIJ - PHASE #  1
  LIJPH      0.4999    0.0000    0.0000
  LIJPH      0.0000    0.4998    0.0000
  LIJPH      0.0000    0.0000   -0.9997
  FIJ        1.2698    0.0000   -0.0000
  FIJ        0.0000    1.2698    0.0000
  FIJ       -0.0000    0.0000    0.6157
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  1
  EULER ANGLES   -12.83     0.00     0.00
  EIGNVAL   1.2698                     1.2698                     0.6157
  EIGNVEC   0.9750  -0.2221  -0.0000   0.2221   0.9750   0.0000   0.0000  -0.0000   1.0000

*******   STEP  25      *** N_EFF=10 CALCULATION               
    6   0.272E-03   0.398E-03   0.276E-03
  INSIDE UPD_FIJ - PHASE #  0
  LIJBAR     0.5000    0.0000    0.0000
  LIJBAR     0.0000    0.5000    0.0000
  LIJBAR     0.0000    0.0000   -1.0000
  FIJ        1.2824    0.0000    0.0000
  FIJ        0.0000    1.2824    0.0000
  FIJ        0.0000    0.0000    0.6035
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  0
  EULER ANGLES     0.00     0.00     0.00
  EIGNVAL   1.2824                     1.2824                     0.6035
  EIGNVEC   1.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000   1.0000
  INSIDE UPD_FIJ - PHASE #  1
  LIJPH      0.4999   -0.0000    0.0000
  LIJPH     -0.0000    0.4998   -0.0000
  LIJPH      0.0000   -0.0000   -0.9997
  FIJ        1.2825    0.0000   -0.0000
  FIJ        0.0000    1.2825    0.0000
  FIJ       -0.0000    0.0000    0.6034
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  1
  EULER ANGLES   -13.67     0.00     0.00
  EIGNVAL   1.2825                     1.2825                     0.6034
  EIGNVEC   0.9717  -0.2363  -0.0000   0.2363   0.9717   0.0000   0.0000  -0.0000   1.0000

*******   STEP  26      *** N_EFF=10 CALCULATION               
    6   0.286E-03   0.403E-03   0.280E-03
  INSIDE UPD_FIJ - PHASE #  0
  LIJBAR     0.5000    0.0000    0.0000
  LIJBAR     0.0000    0.5000    0.0000
  LIJBAR     0.0000    0.0000   -1.0000
  FIJ        1.2953    0.0000    0.0000
  FIJ        0.0000    1.2953    0.0000
  FIJ        0.0000    0.0000    0.5914
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  0
  EULER ANGLES     0.00     0.00     0.00
  EIGNVAL   1.2953                     1.2953                     0.5914
  EIGNVEC   1.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000   1.0000
  INSIDE UPD_FIJ - PHASE #  1
  LIJPH      0.4999   -0.0000    0.0000
  LIJPH     -0.0000    0.4998   -0.0000
  LIJPH      0.0000   -0.0000   -0.9997
  FIJ        1.2953    0.0000   -0.0000
  FIJ        0.0000    1.2953    0.0000
  FIJ       -0.0000    0.0000    0.5914
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  1
  EULER ANGLES   -14.39     0.00     0.00
  EIGNVAL   1.2953                     1.2953                     0.5914
  EIGNVEC   0.9686  -0.2485  -0.0000   0.2485   0.9686   0.0000   0.0000  -0.0000   1.0000

*******   STEP  27      *** N_EFF=10 CALCULATION               
    6   0.294E-03   0.401E-03   0.297E-03
  INSIDE UPD_FIJ - PHASE #  0
  LIJBAR     0.5000    0.0000    0.0000
  LIJBAR     0.0000    0.5000    0.0000
  LIJBAR     0.0000    0.0000   -1.0000
  FIJ        1.3082    0.0000    0.0000
  FIJ        0.0000    1.3082    0.0000
  FIJ        0.0000    0.0000    0.5796
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  0
  EULER ANGLES     0.00     0.00     0.00
  EIGNVAL   1.3082                     1.3082                     0.5796
  EIGNVEC   1.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000   1.0000
  INSIDE UPD_FIJ - PHASE #  1
  LIJPH      0.4999   -0.0000    0.0000
  LIJPH     -0.0000    0.4998   -0.0000
  LIJPH      0.0000   -0.0000   -0.9997
  FIJ        1.3082    0.0000   -0.0000
  FIJ        0.0000    1.3082    0.0000
  FIJ       -0.0000    0.0000    0.5795
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  1
  EULER ANGLES   -15.04     0.00     0.00
  EIGNVAL   1.3082                     1.3082                     0.5795
  EIGNVEC   0.9657  -0.2595  -0.0000   0.2595   0.9657   0.0000   0.0000  -0.0000   1.0000

*******   STEP  28      *** N_EFF=10 CALCULATION               
    6   0.292E-03   0.383E-03   0.332E-03
  INSIDE UPD_FIJ - PHASE #  0
  LIJBAR     0.5000    0.0000    0.0000
  LIJBAR     0.0000    0.5000    0.0000
  LIJBAR     0.0000    0.0000   -1.0000
  FIJ        1.3213    0.0000    0.0000
  FIJ        0.0000    1.3213    0.0000
  FIJ        0.0000    0.0000    0.5680
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  0
  EULER ANGLES     0.00     0.00     0.00
  EIGNVAL   1.3213                     1.3213                     0.5680
  EIGNVEC   1.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000   1.0000
  INSIDE UPD_FIJ - PHASE #  1
  LIJPH      0.4999   -0.0000    0.0000
  LIJPH     -0.0000    0.4998   -0.0000
  LIJPH      0.0000   -0.0000   -0.9997
  FIJ        1.3213    0.0000   -0.0000
  FIJ        0.0000    1.3213    0.0000
  FIJ       -0.0000    0.0000    0.5679
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  1
  EULER ANGLES   -15.72     0.00     0.00
  EIGNVAL   1.3213                     1.3213                     0.5679
  EIGNVEC   0.9626  -0.2709  -0.0000   0.2709   0.9626   0.0000   0.0000  -0.0000   1.0000

*******   STEP  29      *** N_EFF=10 CALCULATION               
    6   0.276E-03   0.341E-03   0.397E-03
  INSIDE UPD_FIJ - PHASE #  0
  LIJBAR     0.5000    0.0000    0.0000
  LIJBAR     0.0000    0.5000    0.0000
  LIJBAR     0.0000    0.0000   -1.0000
  FIJ        1.3345    0.0000    0.0000
  FIJ        0.0000    1.3345    0.0000
  FIJ        0.0000    0.0000    0.5566
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  0
  EULER ANGLES     0.00     0.00     0.00
  EIGNVAL   1.3345                     1.3345                     0.5566
  EIGNVEC   1.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000   1.0000
  INSIDE UPD_FIJ - PHASE #  1
  LIJPH      0.4998   -0.0000    0.0000
  LIJPH     -0.0000    0.4998   -0.0000
  LIJPH      0.0000   -0.0000   -0.9996
  FIJ        1.3345    0.0000   -0.0000
  FIJ        0.0000    1.3345    0.0000
  FIJ       -0.0000    0.0000    0.5566
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  1
  EULER ANGLES   -16.51     0.00     0.00
  EIGNVAL   1.3345                     1.3345                     0.5566
  EIGNVEC   0.9588  -0.2842  -0.0000   0.2842   0.9588   0.0000   0.0000  -0.0000   1.0000

*******   STEP  30      *** N_EFF=10 CALCULATION               
    6   0.239E-03   0.267E-03   0.502E-03
  INSIDE UPD_FIJ - PHASE #  0
  LIJBAR     0.5000    0.0000    0.0000
  LIJBAR     0.0000    0.5000    0.0000
  LIJBAR     0.0000    0.0000   -1.0000
  FIJ        1.3478    0.0000    0.0000
  FIJ        0.0000    1.3478    0.0000
  FIJ        0.0000    0.0000    0.5455
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  0
  EULER ANGLES     0.00     0.00     0.00
  EIGNVAL   1.3478                     1.3478                     0.5455
  EIGNVEC   1.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000   1.0000
  INSIDE UPD_FIJ - PHASE #  1
  LIJPH      0.4998   -0.0000    0.0000
  LIJPH     -0.0000    0.4997   -0.0000
  LIJPH      0.0000   -0.0000   -0.9995
  FIJ        1.3479    0.0000   -0.0000
  FIJ        0.0000    1.3479    0.0000
  FIJ       -0.0000    0.0000    0.5455
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  1
  EULER ANGLES   -17.48     0.00     0.00
  EIGNVAL   1.3479                     1.3479                     0.5455
  EIGNVEC   0.9538  -0.3004  -0.0000   0.3004   0.9538   0.0000   0.0000  -0.0000   1.0000

*******   STEP  31      *** N_EFF=10 CALCULATION               
    6   0.178E-03   0.175E-03   0.657E-03
  INSIDE UPD_FIJ - PHASE #  0
  LIJBAR     0.5000    0.0000    0.0000
  LIJBAR     0.0000    0.5000    0.0000
  LIJBAR     0.0000    0.0000   -1.0000
  FIJ        1.3613    0.0000    0.0000
  FIJ        0.0000    1.3613    0.0000
  FIJ        0.0000    0.0000    0.5346
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  0
  EULER ANGLES     0.00     0.00     0.00
  EIGNVAL   1.3613                     1.3613                     0.5346
  EIGNVEC   1.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000   1.0000
  INSIDE UPD_FIJ - PHASE #  1
  LIJPH      0.4997   -0.0000    0.0000
  LIJPH     -0.0000    0.4997   -0.0000
  LIJPH      0.0000   -0.0000   -0.9993
  FIJ        1.3613    0.0000   -0.0000
  FIJ        0.0000    1.3613    0.0000
  FIJ       -0.0000    0.0000    0.5346
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  1
  EULER ANGLES   -18.70     0.00     0.00
  EIGNVAL   1.3613                     1.3613                     0.5346
  EIGNVEC   0.9472  -0.3205  -0.0000   0.3205   0.9472   0.0000   0.0000  -0.0000   1.0000

*******   STEP  32      *** N_EFF=10 CALCULATION               
    5   0.545E-03   0.950E-03   0.173E-03
  INSIDE UPD_FIJ - PHASE #  0
  LIJBAR     0.5000    0.0000    0.0000
  LIJBAR     0.0000    0.5000    0.0000
  LIJBAR     0.0000    0.0000   -1.0000
  FIJ        1.3749    0.0000    0.0000
  FIJ        0.0000    1.3749    0.0000
  FIJ        0.0000    0.0000    0.5239
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  0
  EULER ANGLES     0.00     0.00     0.00
  EIGNVAL   1.3749                     1.3749                     0.5239
  EIGNVEC   1.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000   1.0000
  INSIDE UPD_FIJ - PHASE #  1
  LIJPH      0.5001   -0.0000    0.0000
  LIJPH     -0.0000    0.5001   -0.0000
  LIJPH      0.0000   -0.0000   -1.0002
  FIJ        1.3750    0.0000   -0.0000
  FIJ        0.0000    1.3750    0.0000
  FIJ       -0.0000    0.0000    0.5239
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  1
  EULER ANGLES   -22.51     0.00     0.00
  EIGNVAL   1.3750                     1.3750                     0.5239
  EIGNVEC   0.9238  -0.3828  -0.0000   0.3828   0.9238   0.0000   0.0000  -0.0000   1.0000

*******   STEP  33      *** N_EFF=10 CALCULATION               
    5   0.524E-03   0.865E-03   0.860E-04
  INSIDE UPD_FIJ - PHASE #  0
  LIJBAR     0.5000    0.0000    0.0000
  LIJBAR     0.0000    0.5000    0.0000
  LIJBAR     0.0000    0.0000   -1.0000
  FIJ        1.3887    0.0000    0.0000
  FIJ        0.0000    1.3887    0.0000
  FIJ        0.0000    0.0000    0.5134
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  0
  EULER ANGLES     0.00     0.00     0.00
  EIGNVAL   1.3887                     1.3887                     0.5134
  EIGNVEC   1.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000   1.0000
  INSIDE UPD_FIJ - PHASE #  1
  LIJPH      0.5001   -0.0000    0.0000
  LIJPH     -0.0000    0.5000   -0.0000
  LIJPH      0.0000   -0.0000   -1.0001
  FIJ        1.3887    0.0000   -0.0000
  FIJ        0.0000    1.3887    0.0000
  FIJ       -0.0000    0.0000    0.5134
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  1
  EULER ANGLES   -37.67     0.00     0.00
  EIGNVAL   1.3887                     1.3887                     0.5134
  EIGNVEC   0.7915  -0.6112  -0.0000   0.6112   0.7915   0.0000   0.0000  -0.0000   1.0000

*******   STEP  34      *** N_EFF=10 CALCULATION               
    5   0.522E-03   0.673E-03   0.217E-03
  INSIDE UPD_FIJ - PHASE #  0
  LIJBAR     0.5000    0.0000    0.0000
  LIJBAR     0.0000    0.5000    0.0000
  LIJBAR     0.0000    0.0000   -1.0000
  FIJ        1.4026    0.0000    0.0000
  FIJ        0.0000    1.4026    0.0000
  FIJ        0.0000    0.0000    0.5031
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  0
  EULER ANGLES     0.00     0.00     0.00
  EIGNVAL   1.4026                     1.4026                     0.5031
  EIGNVEC   1.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000   1.0000
  INSIDE UPD_FIJ - PHASE #  1
  LIJPH      0.4999   -0.0000    0.0000
  LIJPH     -0.0000    0.4999   -0.0000
  LIJPH      0.0000   -0.0000   -0.9998
  FIJ        1.4026    0.0000   -0.0000
  FIJ        0.0000    1.4026    0.0000
  FIJ       -0.0000    0.0000    0.5031
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  1
  EULER ANGLES   101.21     0.00     0.00
  EIGNVAL   1.4026                     1.4026                     0.5031
  EIGNVEC  -0.1945   0.9809   0.0000  -0.9809  -0.1945   0.0000   0.0000  -0.0000   1.0000

*******   STEP  35      *** N_EFF=10 CALCULATION               
    5   0.177E-03   0.342E-03   0.474E-03
  INSIDE UPD_FIJ - PHASE #  0
  LIJBAR     0.5000    0.0000    0.0000
  LIJBAR     0.0000    0.5000    0.0000
  LIJBAR     0.0000    0.0000   -1.0000
  FIJ        1.4166    0.0000    0.0000
  FIJ        0.0000    1.4166    0.0000
  FIJ        0.0000    0.0000    0.4931
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  0
  EULER ANGLES     0.00     0.00     0.00
  EIGNVAL   1.4166                     1.4166                     0.4931
  EIGNVEC   1.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000   1.0000
  INSIDE UPD_FIJ - PHASE #  1
  LIJPH      0.4998   -0.0000   -0.0000
  LIJPH     -0.0000    0.4998   -0.0000
  LIJPH     -0.0000   -0.0000   -0.9995
  FIJ        1.4166    0.0000   -0.0000
  FIJ        0.0000    1.4166    0.0000
  FIJ       -0.0000    0.0000    0.4931
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  1
  EULER ANGLES    99.55     0.00     0.00
  EIGNVAL   1.4166                     1.4166                     0.4931
  EIGNVEC  -0.1659   0.9862   0.0000  -0.9862  -0.1659   0.0000   0.0000  -0.0000   1.0000

*******   STEP  36      *** N_EFF=10 CALCULATION               
    5   0.358E-03   0.217E-03   0.679E-03
  INSIDE UPD_FIJ - PHASE #  0
  LIJBAR     0.5000    0.0000    0.0000
  LIJBAR     0.0000    0.5000    0.0000
  LIJBAR     0.0000    0.0000   -1.0000
  FIJ        1.4308    0.0000    0.0000
  FIJ        0.0000    1.4308    0.0000
  FIJ        0.0000    0.0000    0.4832
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  0
  EULER ANGLES     0.00     0.00     0.00
  EIGNVAL   1.4308                     1.4308                     0.4832
  EIGNVEC   1.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000   1.0000
  INSIDE UPD_FIJ - PHASE #  1
  LIJPH      0.4997    0.0000   -0.0000
  LIJPH      0.0000    0.4997   -0.0000
  LIJPH     -0.0000   -0.0000   -0.9993
  FIJ        1.4308    0.0000   -0.0000
  FIJ        0.0000    1.4308    0.0000
  FIJ       -0.0000    0.0000    0.4832
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  1
  EULER ANGLES   105.24     0.00     0.00
  EIGNVAL   1.4308                     1.4308                     0.4832
  EIGNVEC  -0.2629   0.9648   0.0000  -0.9648  -0.2629   0.0000   0.0000  -0.0000   1.0000

*******   STEP  37      *** N_EFF=10 CALCULATION               
    5   0.292E-03   0.243E-03   0.546E-03
  INSIDE UPD_FIJ - PHASE #  0
  LIJBAR     0.5000    0.0000    0.0000
  LIJBAR     0.0000    0.5000    0.0000
  LIJBAR     0.0000    0.0000   -1.0000
  FIJ        1.4451    0.0000    0.0000
  FIJ        0.0000    1.4451    0.0000
  FIJ        0.0000    0.0000    0.4735
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  0
  EULER ANGLES     0.00     0.00     0.00
  EIGNVAL   1.4451                     1.4451                     0.4735
  EIGNVEC   1.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000   1.0000
  INSIDE UPD_FIJ - PHASE #  1
  LIJPH      0.4997    0.0000   -0.0000
  LIJPH      0.0000    0.4997   -0.0000
  LIJPH     -0.0000   -0.0000   -0.9995
  FIJ        1.4451    0.0000   -0.0000
  FIJ        0.0000    1.4451    0.0000
  FIJ       -0.0000    0.0000    0.4736
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  1
  EULER ANGLES   115.89     0.00     0.00
  EIGNVAL   1.4451                     1.4451                     0.4736
  EIGNVEC  -0.4366   0.8997   0.0000  -0.8997  -0.4366   0.0000   0.0000  -0.0000   1.0000

*******   STEP  38      *** N_EFF=10 CALCULATION               
    5   0.343E-03   0.380E-03   0.582E-03
  INSIDE UPD_FIJ - PHASE #  0
  LIJBAR     0.5000    0.0000    0.0000
  LIJBAR     0.0000    0.5000    0.0000
  LIJBAR     0.0000    0.0000   -1.0000
  FIJ        1.4595    0.0000    0.0000
  FIJ        0.0000    1.4595    0.0000
  FIJ        0.0000    0.0000    0.4641
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  0
  EULER ANGLES     0.00     0.00     0.00
  EIGNVAL   1.4595                     1.4595                     0.4641
  EIGNVEC   1.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000   1.0000
  INSIDE UPD_FIJ - PHASE #  1
  LIJPH      0.4997    0.0000   -0.0000
  LIJPH      0.0000    0.4997   -0.0000
  LIJPH     -0.0000   -0.0000   -0.9994
  FIJ        1.4595    0.0000   -0.0000
  FIJ        0.0000    1.4595    0.0000
  FIJ       -0.0000    0.0000    0.4641
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  1
  EULER ANGLES   -40.94     0.00     0.00
  EIGNVAL   1.4595                     1.4595                     0.4641
  EIGNVEC   0.7554  -0.6552  -0.0000   0.6552   0.7554  -0.0000   0.0000  -0.0000   1.0000

*******   STEP  39      *** N_EFF=10 CALCULATION               
    4   0.926E-04   0.922E-03   0.259E-03
  INSIDE UPD_FIJ - PHASE #  0
  LIJBAR     0.5000    0.0000    0.0000
  LIJBAR     0.0000    0.5000    0.0000
  LIJBAR     0.0000    0.0000   -1.0000
  FIJ        1.4741    0.0000    0.0000
  FIJ        0.0000    1.4741    0.0000
  FIJ        0.0000    0.0000    0.4548
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  0
  EULER ANGLES     0.00     0.00     0.00
  EIGNVAL   1.4741                     1.4741                     0.4548
  EIGNVEC   1.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000   1.0000
  INSIDE UPD_FIJ - PHASE #  1
  LIJPH      0.4999   -0.0000    0.0000
  LIJPH     -0.0000    0.4999   -0.0000
  LIJPH      0.0000   -0.0000   -0.9997
  FIJ        1.4741    0.0000   -0.0000
  FIJ        0.0000    1.4741    0.0000
  FIJ       -0.0000   -0.0000    0.4548
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  1
  EULER ANGLES   120.02     0.00     0.00
  EIGNVAL   1.4741                     1.4741                     0.4548
  EIGNVEC  -0.5004   0.8658   0.0000  -0.8658  -0.5004   0.0000   0.0000   0.0000   1.0000

*******   STEP  40      *** N_EFF=10 CALCULATION               
    4   0.599E-04   0.842E-03   0.375E-03
  INSIDE UPD_FIJ - PHASE #  0
  LIJBAR     0.5000    0.0000    0.0000
  LIJBAR     0.0000    0.5000    0.0000
  LIJBAR     0.0000    0.0000   -1.0000
  FIJ        1.4889    0.0000    0.0000
  FIJ        0.0000    1.4889    0.0000
  FIJ        0.0000    0.0000    0.4457
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  0
  EULER ANGLES     0.00     0.00     0.00
  EIGNVAL   1.4889                     1.4889                     0.4457
  EIGNVEC   1.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000   1.0000
  INSIDE UPD_FIJ - PHASE #  1
  LIJPH      0.4998   -0.0000    0.0000
  LIJPH     -0.0000    0.4998   -0.0000
  LIJPH      0.0000   -0.0000   -0.9996
  FIJ        1.4888    0.0000   -0.0000
  FIJ        0.0000    1.4888   -0.0000
  FIJ       -0.0000   -0.0000    0.4457
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  1
  EULER ANGLES   107.19     0.00     0.00
  EIGNVAL   1.4888                     1.4888                     0.4457
  EIGNVEC  -0.2956   0.9553  -0.0000  -0.9553  -0.2956   0.0000   0.0000   0.0000   1.0000

*******   STEP  41      *** N_EFF=10 CALCULATION               
    4   0.857E-04   0.681E-03   0.483E-03
  INSIDE UPD_FIJ - PHASE #  0
  LIJBAR     0.5000    0.0000    0.0000
  LIJBAR     0.0000    0.5000    0.0000
  LIJBAR     0.0000    0.0000   -1.0000
  FIJ        1.5038    0.0000    0.0000
  FIJ        0.0000    1.5038    0.0000
  FIJ        0.0000    0.0000    0.4368
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  0
  EULER ANGLES     0.00     0.00     0.00
  EIGNVAL   1.5038                     1.5038                     0.4368
  EIGNVEC   1.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000   1.0000
  INSIDE UPD_FIJ - PHASE #  1
  LIJPH      0.4998   -0.0000    0.0000
  LIJPH     -0.0000    0.4998   -0.0000
  LIJPH      0.0000   -0.0000   -0.9995
  FIJ        1.5037    0.0000   -0.0000
  FIJ        0.0000    1.5037   -0.0000
  FIJ       -0.0000   -0.0000    0.4368
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  1
  EULER ANGLES   105.39     0.00     0.00
  EIGNVAL   1.5037                     1.5037                     0.4368
  EIGNVEC  -0.2653   0.9642  -0.0000  -0.9642  -0.2653   0.0000   0.0000   0.0000   1.0000

*******   STEP  42      *** N_EFF=10 CALCULATION               
    5   0.316E-03   0.357E-03   0.554E-03
  INSIDE UPD_FIJ - PHASE #  0
  LIJBAR     0.5000    0.0000    0.0000
  LIJBAR     0.0000    0.5000    0.0000
  LIJBAR     0.0000    0.0000   -1.0000
  FIJ        1.5188    0.0000    0.0000
  FIJ        0.0000    1.5188    0.0000
  FIJ        0.0000    0.0000    0.4281
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  0
  EULER ANGLES     0.00     0.00     0.00
  EIGNVAL   1.5188                     1.5188                     0.4281
  EIGNVEC   1.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000   1.0000
  INSIDE UPD_FIJ - PHASE #  1
  LIJPH      0.4997    0.0000    0.0000
  LIJPH      0.0000    0.4997   -0.0000
  LIJPH      0.0000   -0.0000   -0.9994
  FIJ        1.5187    0.0000   -0.0000
  FIJ        0.0000    1.5187   -0.0000
  FIJ        0.0000   -0.0000    0.4281
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  1
  EULER ANGLES   101.98     0.00     0.00
  EIGNVAL   1.5187                     1.5187                     0.4281
  EIGNVEC  -0.2076   0.9782  -0.0000  -0.9782  -0.2076   0.0000  -0.0000   0.0000   1.0000

*******   STEP  43      *** N_EFF=10 CALCULATION               
    4   0.311E-03   0.986E-03   0.182E-03
  INSIDE UPD_FIJ - PHASE #  0
  LIJBAR     0.5000    0.0000    0.0000
  LIJBAR     0.0000    0.5000    0.0000
  LIJBAR     0.0000    0.0000   -1.0000
  FIJ        1.5340    0.0000    0.0000
  FIJ        0.0000    1.5340    0.0000
  FIJ        0.0000    0.0000    0.4195
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  0
  EULER ANGLES     0.00     0.00     0.00
  EIGNVAL   1.5340                     1.5340                     0.4195
  EIGNVEC   1.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000   1.0000
  INSIDE UPD_FIJ - PHASE #  1
  LIJPH      0.5001   -0.0000    0.0000
  LIJPH     -0.0000    0.5000   -0.0000
  LIJPH      0.0000   -0.0000   -1.0002
  FIJ        1.5339   -0.0000    0.0000
  FIJ       -0.0000    1.5339   -0.0000
  FIJ        0.0000   -0.0000    0.4195
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  1
  EULER ANGLES    88.01     0.00     0.00
  EIGNVAL   1.5339                     1.5339                     0.4195
  EIGNVEC   0.0347   0.9994  -0.0000  -0.9994   0.0347  -0.0000  -0.0000   0.0000   1.0000

*******   STEP  44      *** N_EFF=10 CALCULATION               
    4   0.960E-04   0.407E-03   0.420E-03
  INSIDE UPD_FIJ - PHASE #  0
  LIJBAR     0.5000    0.0000    0.0000
  LIJBAR     0.0000    0.5000    0.0000
  LIJBAR     0.0000    0.0000   -1.0000
  FIJ        1.5493    0.0000    0.0000
  FIJ        0.0000    1.5493    0.0000
  FIJ        0.0000    0.0000    0.4111
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  0
  EULER ANGLES     0.00     0.00     0.00
  EIGNVAL   1.5493                     1.5493                     0.4111
  EIGNVEC   1.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000   1.0000
  INSIDE UPD_FIJ - PHASE #  1
  LIJPH      0.4998    0.0000    0.0000
  LIJPH      0.0000    0.4998   -0.0000
  LIJPH      0.0000   -0.0000   -0.9996
  FIJ        1.5493   -0.0000    0.0000
  FIJ       -0.0000    1.5493   -0.0000
  FIJ        0.0000   -0.0000    0.4111
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  1
  EULER ANGLES    88.32     0.00     0.00
  EIGNVAL   1.5493                     1.5493                     0.4111
  EIGNVEC   0.0294   0.9996  -0.0000  -0.9996   0.0294  -0.0000  -0.0000   0.0000   1.0000

*******   STEP  45      *** N_EFF=10 CALCULATION               
    4   0.247E-03   0.607E-03   0.314E-03
  INSIDE UPD_FIJ - PHASE #  0
  LIJBAR     0.5000    0.0000    0.0000
  LIJBAR     0.0000    0.5000    0.0000
  LIJBAR     0.0000    0.0000   -1.0000
  FIJ        1.5648    0.0000    0.0000
  FIJ        0.0000    1.5648    0.0000
  FIJ        0.0000    0.0000    0.4029
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  0
  EULER ANGLES     0.00     0.00     0.00
  EIGNVAL   1.5648                     1.5648                     0.4029
  EIGNVEC   1.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000   1.0000
  INSIDE UPD_FIJ - PHASE #  1
  LIJPH      0.4999   -0.0000    0.0000
  LIJPH     -0.0000    0.4998   -0.0000
  LIJPH      0.0000   -0.0000   -0.9997
  FIJ        1.5648   -0.0000    0.0000
  FIJ       -0.0000    1.5648   -0.0000
  FIJ        0.0000   -0.0000    0.4029
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  1
  EULER ANGLES    87.98     0.00     0.00
  EIGNVAL   1.5648                     1.5648                     0.4029
  EIGNVEC   0.0352   0.9994  -0.0000  -0.9994   0.0352  -0.0000  -0.0000   0.0000   1.0000

*******   STEP  46      *** N_EFF=10 CALCULATION               
    4   0.232E-03   0.543E-03   0.431E-03
  INSIDE UPD_FIJ - PHASE #  0
  LIJBAR     0.5000    0.0000    0.0000
  LIJBAR     0.0000    0.5000    0.0000
  LIJBAR     0.0000    0.0000   -1.0000
  FIJ        1.5805    0.0000    0.0000
  FIJ        0.0000    1.5805    0.0000
  FIJ        0.0000    0.0000    0.3948
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  0
  EULER ANGLES     0.00     0.00     0.00
  EIGNVAL   1.5805                     1.5805                     0.3948
  EIGNVEC   1.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000   1.0000
  INSIDE UPD_FIJ - PHASE #  1
  LIJPH      0.4998   -0.0000    0.0000
  LIJPH     -0.0000    0.4998   -0.0000
  LIJPH      0.0000   -0.0000   -0.9996
  FIJ        1.5804   -0.0000    0.0000
  FIJ       -0.0000    1.5804   -0.0000
  FIJ        0.0000   -0.0000    0.3949
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  1
  EULER ANGLES    88.10     0.00     0.00
  EIGNVAL   1.5804                     1.5804                     0.3949
  EIGNVEC   0.0331   0.9995  -0.0000  -0.9995   0.0331  -0.0000  -0.0000   0.0000   1.0000

*******   STEP  47      *** N_EFF=10 CALCULATION               
    4   0.119E-03   0.417E-03   0.411E-03
  INSIDE UPD_FIJ - PHASE #  0
  LIJBAR     0.5000    0.0000    0.0000
  LIJBAR     0.0000    0.5000    0.0000
  LIJBAR     0.0000    0.0000   -1.0000
  FIJ        1.5963    0.0000    0.0000
  FIJ        0.0000    1.5963    0.0000
  FIJ        0.0000    0.0000    0.3869
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  0
  EULER ANGLES     0.00     0.00     0.00
  EIGNVAL   1.5963                     1.5963                     0.3869
  EIGNVEC   1.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000   1.0000
  INSIDE UPD_FIJ - PHASE #  1
  LIJPH      0.4998    0.0000    0.0000
  LIJPH      0.0000    0.4998   -0.0000
  LIJPH      0.0000   -0.0000   -0.9996
  FIJ        1.5962   -0.0000    0.0000
  FIJ       -0.0000    1.5962   -0.0000
  FIJ        0.0000   -0.0000    0.3870
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  1
  EULER ANGLES    88.49     0.00     0.00
  EIGNVAL   1.5962                     1.5962                     0.3870
  EIGNVEC   0.0264   0.9997  -0.0000  -0.9997   0.0264  -0.0000  -0.0000   0.0000   1.0000

*******   STEP  48      *** N_EFF=10 CALCULATION               
    3   0.178E-03   0.671E-03   0.401E-03
  INSIDE UPD_FIJ - PHASE #  0
  LIJBAR     0.5000    0.0000    0.0000
  LIJBAR     0.0000    0.5000    0.0000
  LIJBAR     0.0000    0.0000   -1.0000
  FIJ        1.6122    0.0000    0.0000
  FIJ        0.0000    1.6122    0.0000
  FIJ        0.0000    0.0000    0.3792
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  0
  EULER ANGLES     0.00     0.00     0.00
  EIGNVAL   1.6122                     1.6122                     0.3792
  EIGNVEC   1.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000   1.0000
  INSIDE UPD_FIJ - PHASE #  1
  LIJPH      0.5002   -0.0000    0.0001
  LIJPH     -0.0000    0.5002   -0.0000
  LIJPH      0.0001   -0.0000   -1.0004
  FIJ        1.6122   -0.0000    0.0000
  FIJ       -0.0000    1.6122   -0.0000
  FIJ        0.0000   -0.0000    0.3792
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  1
  EULER ANGLES    87.53     0.00     0.00
  EIGNVAL   1.6122                     1.6122                     0.3792
  EIGNVEC   0.0430   0.9991  -0.0000  -0.9991   0.0430  -0.0000  -0.0000   0.0000   1.0000

*******   STEP  49      *** N_EFF=10 CALCULATION               
    4   0.157E-03   0.570E-03   0.657E-03
  INSIDE UPD_FIJ - PHASE #  0
  LIJBAR     0.5000    0.0000    0.0000
  LIJBAR     0.0000    0.5000    0.0000
  LIJBAR     0.0000    0.0000   -1.0000
  FIJ        1.6283    0.0000    0.0000
  FIJ        0.0000    1.6283    0.0000
  FIJ        0.0000    0.0000    0.3716
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  0
  EULER ANGLES     0.00     0.00     0.00
  EIGNVAL   1.6283                     1.6283                     0.3716
  EIGNVEC   1.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000   1.0000
  INSIDE UPD_FIJ - PHASE #  1
  LIJPH      0.4997    0.0000    0.0001
  LIJPH      0.0000    0.4997   -0.0000
  LIJPH      0.0001   -0.0000   -0.9993
  FIJ        1.6283   -0.0000    0.0000
  FIJ       -0.0000    1.6283   -0.0000
  FIJ        0.0000   -0.0000    0.3716
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  1
  EULER ANGLES    87.91     0.00     0.00
  EIGNVAL   1.6283                     1.6283                     0.3716
  EIGNVEC   0.0364   0.9993  -0.0000  -0.9993   0.0364  -0.0000  -0.0000   0.0000   1.0000

*******   STEP  50      *** N_EFF=10 CALCULATION               
    4   0.174E-03   0.677E-03   0.452E-03
  INSIDE UPD_FIJ - PHASE #  0
  LIJBAR     0.5000    0.0000    0.0000
  LIJBAR     0.0000    0.5000    0.0000
  LIJBAR     0.0000    0.0000   -1.0000
  FIJ        1.6446    0.0000    0.0000
  FIJ        0.0000    1.6446    0.0000
  FIJ        0.0000    0.0000    0.3642
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  0
  EULER ANGLES     0.00     0.00     0.00
  EIGNVAL   1.6446                     1.6446                     0.3642
  EIGNVEC   1.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000   1.0000
  INSIDE UPD_FIJ - PHASE #  1
  LIJPH      0.4998    0.0000    0.0001
  LIJPH      0.0000    0.4998   -0.0000
  LIJPH      0.0001   -0.0000   -0.9996
  FIJ        1.6446   -0.0000    0.0000
  FIJ       -0.0000    1.6445   -0.0000
  FIJ        0.0000   -0.0000    0.3642
  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #  1
  EULER ANGLES    88.32     0.00     0.00
  EIGNVAL   1.6445                     1.6446                     0.3642
  EIGNVEC   0.0293   0.9996  -0.0000  -0.9996   0.0293  -0.0000  -0.0000   0.0000   1.0000

*******   STEP  51      *** N_EFF=10 CALCULATION               
    4   0.198E-03   0.836E-03   0.439E-03

 TIME ELAPSED    6.36secs
