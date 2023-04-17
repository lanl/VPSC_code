
 **** INPUT FILE VPSC.IN FOR THIS RUN ****
1                          VP (1) or ELAST (-1) regime                                                                            
1                          number of phases (nph)                                                                                 
1.0  0.0                   relative vol. fract. of phases (wph(i))                                                                
*INFORMATION ABOUT PHASE #1                                                                                                       
0   0   25                      grain shape contrl, fragmentn, crit aspect ratio                                                  
1.0  1.0  1.0                 initial ellipsoid ratios (dummy if ishape=4)                                                        
0.0  0.0  0.0                 init Eul ang ellips axes (dummy if ishape=3,4)                                                      
* name and path of texture file (filetext)                                                                                        
ex02_FCC\rand500.tex                                                                                                              
* name and path of single crystal file (filecrys)                                                                                 
ex02_FCC\FCC.sx                                                                                                                   
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
3  10          interaction (0:FC,1:aff,2:sec,3:neff or RDC,4:tan,5:SO), neffx                                                     
1  1  1        iupdate: update orient, grain shape, hardening                                                                     
0              nneigh (0 for no neighbors, 1 for pairs, etc.)                                                                     
0              iflu (0: don't calc, 1: calc fluctuations)                                                                         
*NUMBER OF PROCESSES (Lij const; Lij variable; PCYS ;LANKFORD; rigid rotatn)                                                      
1                                                                                                                                 
*IVGVAR AND PATH\NAME OF FILE FOR EACH PROCESS                                                                                    
0                                                                                                                                 
ex02_FCC\COMPRES_extrusion.3                                                                                                      


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


*** RUNNING INTERACTION NEFF=   10.00


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

*******   STEP   1      *** N_EFF=cte CALCULATION              
   10   0.927E-04   0.390E-03   0.650E-03

*******   STEP   2      *** N_EFF=cte CALCULATION              
    6   0.125E-03   0.579E-03   0.176E-03

*******   STEP   3      *** N_EFF=cte CALCULATION              
    6   0.147E-03   0.666E-03   0.265E-03

*******   STEP   4      *** N_EFF=cte CALCULATION              
    5   0.102E-03   0.982E-03   0.557E-03

*******   STEP   5      *** N_EFF=cte CALCULATION              
    5   0.119E-03   0.911E-03   0.543E-03

*******   STEP   6      *** N_EFF=cte CALCULATION              
    5   0.911E-04   0.776E-03   0.531E-03

*******   STEP   7      *** N_EFF=cte CALCULATION              
    4   0.365E-03   0.802E-03   0.484E-03

*******   STEP   8      *** N_EFF=cte CALCULATION              
    4   0.649E-03   0.827E-03   0.632E-03

*******   STEP   9      *** N_EFF=cte CALCULATION              
    4   0.604E-03   0.736E-03   0.597E-03

*******   STEP  10      *** N_EFF=cte CALCULATION              
    4   0.954E-03   0.751E-03   0.933E-03

*******   STEP  11      *** N_EFF=cte CALCULATION              
    5   0.199E-03   0.875E-03   0.487E-03

*******   STEP  12      *** N_EFF=cte CALCULATION              
    4   0.924E-03   0.924E-03   0.982E-03

*******   STEP  13      *** N_EFF=cte CALCULATION              
    6   0.278E-03   0.709E-03   0.219E-03

*******   STEP  14      *** N_EFF=cte CALCULATION              
    6   0.206E-03   0.488E-03   0.533E-04

*******   STEP  15      *** N_EFF=cte CALCULATION              
    6   0.189E-03   0.437E-03   0.208E-03

*******   STEP  16      *** N_EFF=cte CALCULATION              
    6   0.186E-03   0.434E-03   0.331E-03

*******   STEP  17      *** N_EFF=cte CALCULATION              
    6   0.192E-03   0.445E-03   0.424E-03

*******   STEP  18      *** N_EFF=cte CALCULATION              
    6   0.207E-03   0.445E-03   0.491E-03

*******   STEP  19      *** N_EFF=cte CALCULATION              
    6   0.228E-03   0.434E-03   0.534E-03

*******   STEP  20      *** N_EFF=cte CALCULATION              
    6   0.167E-03   0.419E-03   0.335E-03

*******   STEP  21      *** N_EFF=cte CALCULATION              
    6   0.185E-03   0.407E-03   0.312E-03

*******   STEP  22      *** N_EFF=cte CALCULATION              
    6   0.206E-03   0.395E-03   0.305E-03

*******   STEP  23      *** N_EFF=cte CALCULATION              
    6   0.229E-03   0.389E-03   0.294E-03

*******   STEP  24      *** N_EFF=cte CALCULATION              
    6   0.252E-03   0.392E-03   0.282E-03

*******   STEP  25      *** N_EFF=cte CALCULATION              
    6   0.272E-03   0.398E-03   0.276E-03

*******   STEP  26      *** N_EFF=cte CALCULATION              
    6   0.286E-03   0.403E-03   0.280E-03

*******   STEP  27      *** N_EFF=cte CALCULATION              
    6   0.294E-03   0.401E-03   0.297E-03

*******   STEP  28      *** N_EFF=cte CALCULATION              
    6   0.292E-03   0.383E-03   0.332E-03

*******   STEP  29      *** N_EFF=cte CALCULATION              
    6   0.276E-03   0.341E-03   0.397E-03

*******   STEP  30      *** N_EFF=cte CALCULATION              
    6   0.239E-03   0.267E-03   0.502E-03

*******   STEP  31      *** N_EFF=cte CALCULATION              
    6   0.178E-03   0.175E-03   0.657E-03

*******   STEP  32      *** N_EFF=cte CALCULATION              
    5   0.545E-03   0.950E-03   0.173E-03

*******   STEP  33      *** N_EFF=cte CALCULATION              
    5   0.524E-03   0.865E-03   0.860E-04

*******   STEP  34      *** N_EFF=cte CALCULATION              
    5   0.522E-03   0.673E-03   0.217E-03

*******   STEP  35      *** N_EFF=cte CALCULATION              
    5   0.177E-03   0.342E-03   0.474E-03

*******   STEP  36      *** N_EFF=cte CALCULATION              
    5   0.358E-03   0.217E-03   0.679E-03

*******   STEP  37      *** N_EFF=cte CALCULATION              
    5   0.292E-03   0.243E-03   0.546E-03

*******   STEP  38      *** N_EFF=cte CALCULATION              
    5   0.343E-03   0.380E-03   0.582E-03

*******   STEP  39      *** N_EFF=cte CALCULATION              
    4   0.926E-04   0.922E-03   0.259E-03

*******   STEP  40      *** N_EFF=cte CALCULATION              
    4   0.599E-04   0.842E-03   0.375E-03

*******   STEP  41      *** N_EFF=cte CALCULATION              
    4   0.857E-04   0.681E-03   0.483E-03

*******   STEP  42      *** N_EFF=cte CALCULATION              
    5   0.316E-03   0.357E-03   0.554E-03

*******   STEP  43      *** N_EFF=cte CALCULATION              
    4   0.311E-03   0.986E-03   0.182E-03

*******   STEP  44      *** N_EFF=cte CALCULATION              
    4   0.960E-04   0.407E-03   0.420E-03

*******   STEP  45      *** N_EFF=cte CALCULATION              
    4   0.247E-03   0.607E-03   0.314E-03

*******   STEP  46      *** N_EFF=cte CALCULATION              
    4   0.232E-03   0.543E-03   0.431E-03

*******   STEP  47      *** N_EFF=cte CALCULATION              
    4   0.119E-03   0.417E-03   0.411E-03

*******   STEP  48      *** N_EFF=cte CALCULATION              
    3   0.178E-03   0.671E-03   0.401E-03

*******   STEP  49      *** N_EFF=cte CALCULATION              
    4   0.157E-03   0.570E-03   0.657E-03

*******   STEP  50      *** N_EFF=cte CALCULATION              
    4   0.174E-03   0.677E-03   0.452E-03

*******   STEP  51      *** N_EFF=cte CALCULATION              
    4   0.198E-03   0.836E-03   0.439E-03

 TIME ELAPSED    5.36secs
