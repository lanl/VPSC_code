
 **** INPUT FILE VPSC.IN FOR THIS RUN ****
1                          number of elements (nelem)                           
1                          number of phases (nph)                               
1.0  0.0                   relative vol. fract. of phases (wph(i))              
*INFORMATION ABOUT PHASE #1                                                     
0   0   25                    grain shape contrl, fragmentn, crit aspect ratio  
1.0  1.0  1.0                 initial ellipsoid ratios (dummy if ishape=4)      
0.0  0.0  0.0                 init Eul ang ellips axes (dummy if ishape=3,4)    
* name and path of texture file (filetext)                                      
ex00_elast\Zr_bar.tex                                                           
* name and path of single crystal file (filecrys)                               
ex00_elast\zr_76K.sx                                                            
* name and path of grain shape file (dummy if ishape=0) (fileaxes)              
shape1.100                                                                      
* name and path of diffraction file (dummy if idiff=0)                          
0                                                                               
ex00_elast\hcp.dif                                                              
*PRECISION SETTINGS FOR CONVERGENCE PROCEDURES (default values)                 
1.e-3  1.e-3  1.e-3   0.01     errs,errd,errm,errso                             
100 100 25     itmax:   max # of iter, external, internal and SO loops          
0  2  10  2    irsvar & jrsini,jrsfin,jrstep (dummy if irsvar=0)                
*INPUT/OUTPUT SETTINGS FOR THE RUN (default is zero)                            
0              irecover:read grain states from POSTMORT.IN (1) or not (0)?      
0              isave:   write grain states in POSTMORT.OUT at step 'isave'?     
0              icubcomp:calculate fcc rolling components?                       
0              nwrite (frequency of texture downloads)                          
*MODELING CONDITIONS FOR THE RUN                                                
-1             interaction (-1:elast,0:FC,1:aff,2:sec,3:neff=10,4:tang,5:SO)    
1  1  1        iupdate: update orient, grain shape, hardening                   
0              nneigh (0 for no neighbors, 1 for pairs, etc.)                   
0              iflu (0: don't calc, 1: calc fluctuations)                       
*NUMBER OF PROCESSES (Lij const; Lij variable; PCYS ;LANKFORD; rigid rotatn)    
1                                                                               
*IVGVAR AND PATH\NAME OF FILE FOR EACH PROCESS (dummy if ivgvar=2,3)            
0                  -->   cooling from 900 K to 300 K with small E33             
ex00_elast\TENSION_PLUS_COOLING                                                 


 **** CRYSTAL DATA FILE ****
*Material: zirconium                                                            
HEXAGONAL                                crysym                                 
1.  1.  1.594     90.  90.  120.         cdim(i),cang(i)                        
*Elastic stiffness of single crystal [MPa]                                      
  143.5e3    72.5e3    65.4e3     0.0     0.0     0.0                           
   72.5e3   143.5e3    65.4e3     0.0     0.0     0.0                           
   65.4e3    65.4e3   164.9e3     0.0     0.0     0.0                           
    0.0     0.0     0.0    32.1e3     0.0e3     0.0e3                           
    0.0     0.0     0.0     0.0e3    32.1e3     0.0e3                           
    0.0     0.0     0.0     0.0e3     0.0e3    35.5e3                           
*Thermal expansion coefficients of single crystal[K^(-1)]                       
  5.7e-6   5.7e-6  10.3e-6   0.0e0   0.0e0   0.0e0                              
SLIP AND TWINNING MODES                                                         
6                               nmodesx                                         
4                               nmodes                                          
1  4  5  6                      mode(i)                                         
PRISMATIC <a>                                                                   
  1    3   1   0                       modex,nsmx,isensex,itwtypex              
 1  0 -1  0    -1  2 -1  0                                                      
 0 -1  1  0     2 -1 -1  0                                                      
 -1  1  0  0    -1 -1  2  0                                                     
BASAL <a>                                                                       
  2    3   1   0                       modex,nsmx,isensex,itwtypex              
 0  0  0  1     2 -1 -1  0                                                      
 0  0  0  1    -1  2 -1  0                                                      
 0  0  0  1    -1 -1  2  0                                                      
PYRAMIDAL<a>                                                                    
  3    6   1   0                       modex,nsmx,isensex,itwtypex              
 1  0 -1  1    -1  2 -1  0                                                      
 0 -1  1  1     2 -1 -1  0                                                      
 -1  1  0  1    -1 -1  2  0                                                     
 -1  0  1  1    -1  2 -1  0                                                     
 0  1 -1  1     2 -1 -1  0                                                      
 1 -1  0  1     1  1 -2  0                                                      
PYRAMIDAL<c+a>                                                                  
  4   12   1   0                       modex,nsmx,isensex,itwtypex              
 1  0 -1  1    -1 -1  2  3                                                      
 1  0 -1  1    -2  1  1  3                                                      
 0 -1  1  1     1  1 -2  3                                                      
 0 -1  1  1    -1  2 -1  3                                                      
 -1  1  0  1     2 -1 -1  3                                                     
 -1  1  0  1     1 -2  1  3                                                     
 -1  0  1  1     2 -1 -1  3                                                     
 -1  0  1  1     1  1 -2  3                                                     
 0  1 -1  1    -1 -1  2  3                                                      
 0  1 -1  1     1 -2  1  3                                                      
 1 -1  0  1    -2  1  1  3                                                      
 1 -1  0  1    -1  2 -1  3                                                      
TENSILE TWIN {10-12}                                                            
  5    6   0   2                       modex,nsmx,isensex,itwtypex              
 0.167                                 twshx                                    
 1  0 -1  2    -1  0  1  1                                                      
 0  1 -1  2     0 -1  1  1                                                      
 -1  1  0  2     1 -1  0  1                                                     
 -1  0  1  2     1  0 -1  1                                                     
 0 -1  1  2     0  1 -1  1                                                      
 1 -1  0  2    -1  1  0  1                                                      
COMPRESSIVE TWIN {11-22}                                                        
  6    6   0   2                       modex,nsmx,isensex,itwtypex              
 0.225                                 twshx                                    
 2 -1 -1  2     2 -1 -1 -3                                                      
 1  1 -2  2     1  1 -2 -3                                                      
 -1  2 -1  2    -1  2 -1 -3                                                     
 -2  1  1  2    -2  1  1 -3                                                     
 -1 -1  2  2    -1 -1  2 -3                                                     
 1 -2  1  2     1 -2  1 -3                                                      
*Constitutive law                                                               
   0      Voce=0, MTS=1                                                         
   1      iratesens (0:rate insensitive, 1:rate sensitive)                      
   25     grsze --> grain size only matters if HPfactor is non-zero             
 PRISMATIC <a> SLIP -------------------------------------------                 
 20                               nrsx                                          
   45.   42.  1290. 25.   0.      tau0x,tau1x,thet0,thet1, hpfac                
 1.0    1.0   10.0    2.0         hlatex(1,im),im=1,nmodes                      
 PYRAMIDAL <c+a> SLIP -------------------------------------------               
 20                               nrsx                                          
  495. 100. 1000. 15.   0.        tau0x,tau1x,thet0,thet1, hpfac                
 1.0    1.0    2.0   2.0          hlatex(1,im),im=1,nmodes                      
 {10-12} TENSILE TWIN --------------------------------------                    
 20                               nrsx                                          
  102.  17.  100. 30.   0.        tau0x,tau1x,thet0,thet1, hpfac                
 1.0    1.0    10.0   16.0        hlatex(1,im),im=1,nmodes                      
  1   1   0.10   0.50             isectw,itwinlaw,thres1,thres2                 
  {11-22} COMPRESSIVE TWIN --------------------------------------               
 20                               nrsx                                          
  270.  30.  1000. 178.   0.      tau0x,tau1x,thet0,thet1, hpfac                
 1.0    1.0    10.0   5.0         hlatex(1,im),im=1,nmodes                      
  1   1   0.10  0.50              isectw,itwinlaw,thres1,thres2                 
                                                                               
**** END OF CRYSTAL DATA FILE ****


 CHECKING THAT CELCC*SELCC-ID6=0   0.2319549E-15
 *********** PHASE   1
 RANDOM PX BULK & POISSON MODULI   95388.889       0.330
 RANDOM PX ELASTIC CTES C11, C12, C44  144073.333   71046.667   36513.333

 SHEAR MODULUS FOR MODE  1 IN PHASE  1 IS   35500.000
 N & B FOR MODE  1 IN PHASE  1
     0.866     0.500     0.000       -0.500     0.866     0.000
     0.000    -1.000     0.000        1.000     0.000     0.000
    -0.866     0.500     0.000       -0.500    -0.866     0.000
 SHEAR MODULUS FOR MODE  2 IN PHASE  1 IS   40331.548
 N & B FOR MODE  2 IN PHASE  1
     0.761     0.439     0.477       -0.266    -0.460     0.847
     0.761     0.439     0.477       -0.531     0.000     0.847
     0.000    -0.879     0.477        0.266     0.460     0.847
     0.000    -0.879     0.477       -0.266     0.460     0.847
    -0.761     0.439     0.477        0.531     0.000     0.847
    -0.761     0.439     0.477        0.266    -0.460     0.847
    -0.761    -0.439     0.477        0.531     0.000     0.847
    -0.761    -0.439     0.477        0.266     0.460     0.847
     0.000     0.879     0.477       -0.266    -0.460     0.847
     0.000     0.879     0.477        0.266    -0.460     0.847
     0.761    -0.439     0.477       -0.531     0.000     0.847
     0.761    -0.439     0.477       -0.266     0.460     0.847
 SHEAR MODULUS FOR MODE  3 IN PHASE  1 IS   44315.532
 N & B FOR MODE  3 IN PHASE  1
     0.586     0.339     0.736       -0.637    -0.368     0.677
     0.000     0.677     0.736        0.000    -0.736     0.677
    -0.586     0.339     0.736        0.637    -0.368     0.677
    -0.586    -0.339     0.736        0.637     0.368     0.677
     0.000    -0.677     0.736        0.000     0.736     0.677
     0.586    -0.339     0.736       -0.637     0.368     0.677
 SHEAR MODULUS FOR MODE  4 IN PHASE  1 IS   42070.802
 N & B FOR MODE  4 IN PHASE  1
     0.847     0.000     0.531        0.531     0.000    -0.847
     0.424     0.734     0.531        0.266     0.460    -0.847
    -0.424     0.734     0.531       -0.266     0.460    -0.847
    -0.847     0.000     0.531       -0.531     0.000    -0.847
    -0.424    -0.734     0.531       -0.266    -0.460    -0.847
     0.424    -0.734     0.531        0.266    -0.460    -0.847

 **** CRYST TEXTURE (FIRST FEW LINES) ****
AXES OF THE REPRESENTATIVE ELLIPSOID                                            
   1.0   1.0   1.0                                                              
DISCRETE TEXTURE FROM ODF FILE  8101.sod                                        
B    1944                                                                       
      5.00      7.07      5.00     0.00000513                                   
      5.00      7.07     15.00     0.00000626                                   
      5.00      7.07     25.00     0.00000840                                   
      5.00      7.07     35.00     0.00000793                                   
      5.00      7.07     45.00     0.00000569                                   
      5.00      7.07     55.00     0.00000493                                   
      5.00     15.79      5.00     0.00000522                                   
      5.00     15.79     15.00     0.00000567                                   
      5.00     15.79     25.00     0.00000628                                   
      5.00     15.79     35.00     0.00000630                                   
      5.00     15.79     45.00     0.00000542                                   
      5.00     15.79     55.00     0.00000511                                   
      5.00     25.46      5.00     0.00000656                                   
      5.00     25.46     15.00     0.00000468                                   
      5.00     25.46     25.00     0.00000521                                   
      5.00     25.46     35.00     0.00000705                                   
    .........................
 **** END OF CRYST TEXTURE DATA FILE ****


 UB THERMAL TENSOR (VOIGT NOTATION)
  0.8087E-05  0.8020E-05  0.5715E-05 -0.9959E-09 -0.3171E-07  0.1811E-07

 UB ELASTIC STIFFNESS (VOIGT NOTATION)
  0.1469E+06  0.7138E+05  0.6993E+05  0.4825E+01 -0.1429E+03  0.1062E+03
  0.7138E+05  0.1465E+06  0.7006E+05  0.9160E+01  0.2597E+02 -0.9630E+01
  0.6993E+05  0.7006E+05  0.1424E+06 -0.1507E+02  0.3894E+02 -0.4662E+02
  0.4825E+01  0.9160E+01 -0.1507E+02  0.3501E+05 -0.3368E+02  0.5794E+01
 -0.1429E+03  0.2597E+02  0.3894E+02 -0.3368E+02  0.3493E+05  0.4544E+01
  0.1062E+03 -0.9630E+01 -0.4662E+02  0.5794E+01  0.4544E+01  0.3783E+05

 LB THERMAL TENSOR (VOIGT NOTATION)
  0.7874E-05  0.7817E-05  0.6009E-05 -0.3502E-09 -0.2508E-07  0.1608E-07

 LB ELASTIC STIFFNESS (VOIGT NOTATION)
  0.1451E+06  0.7189E+05  0.7033E+05  0.4446E+01 -0.1162E+03  0.9253E+02
  0.7189E+05  0.1448E+06  0.7043E+05  0.8077E+01  0.1725E+02 -0.1176E+02
  0.7033E+05  0.7043E+05  0.1420E+06 -0.1251E+02  0.3844E+02 -0.3878E+02
  0.4446E+01  0.8077E+01 -0.1251E+02  0.3464E+05 -0.3000E+02  0.9816E+01
 -0.1162E+03  0.1725E+02  0.3844E+02 -0.3000E+02  0.3456E+05  0.3762E+01
  0.9253E+02 -0.1176E+02 -0.3878E+02  0.9816E+01  0.3762E+01  0.3668E+05


  ELASTIC ESHELBY TENSOR (VOIGT NOTATION)
     0.53628     0.06387     0.06285     0.00001    -0.00033     0.00025
     0.06393     0.53579     0.06338     0.00003     0.00011    -0.00007
     0.06491     0.06539     0.53642    -0.00004     0.00018    -0.00016
     0.00001     0.00003    -0.00004     0.22975    -0.00010     0.00003
    -0.00030     0.00014     0.00021    -0.00010     0.22946     0.00001
     0.00024    -0.00009    -0.00017     0.00003     0.00001     0.23654

SC ELASTIC STIFFNESS (VOIGT NOTATION)
  0.1459E+06  0.7163E+05  0.7013E+05  0.4636E+01 -0.1293E+03  0.9980E+02
  0.7163E+05  0.1457E+06  0.7025E+05  0.8617E+01  0.2169E+02 -0.1080E+02
  0.7013E+05  0.7025E+05  0.1422E+06 -0.1371E+02  0.3870E+02 -0.4275E+02
  0.4636E+01  0.8617E+01 -0.1371E+02  0.3482E+05 -0.3174E+02  0.7851E+01
 -0.1293E+03  0.2169E+02  0.3870E+02 -0.3174E+02  0.3474E+05  0.4152E+01
  0.9980E+02 -0.1080E+02 -0.4275E+02  0.7851E+01  0.4152E+01  0.3724E+05

SC ELASTIC COMPLIANCE (VOIGT NOTATION)
  0.1010E-04 -0.3366E-05 -0.3317E-05 -0.8856E-09  0.2169E-07 -0.1592E-07
 -0.3366E-05  0.1013E-04 -0.3347E-05 -0.1697E-08 -0.7567E-08  0.4060E-08
 -0.3317E-05 -0.3347E-05  0.1032E-04  0.2655E-08 -0.1088E-07  0.9883E-08
 -0.8856E-09 -0.1697E-08  0.2655E-08  0.7179E-05  0.6556E-08 -0.1512E-08
  0.2169E-07 -0.7567E-08 -0.1088E-07  0.6556E-08  0.7196E-05 -0.8400E-09
 -0.1592E-07  0.4060E-08  0.9883E-08 -0.1512E-08 -0.8400E-09  0.6713E-05

SC THERMAL EXPANSION TENSOR (VOIGT NOTATION)
  0.7977E-05  0.7915E-05  0.5867E-05 -0.6344E-09 -0.2830E-07  0.1724E-07

*** LOAD CONDITIONS FOR THIS RUN
  10  3   0.0001   900  300      nsteps   ictrl   eqincr   temp_i   temp_f      
* boundary conditions                                                           
    0       0       0           iudot     |    flag for vel.grad.               
    1       0       0                     |    (0:unknown-1:known)              
    1       1       1                     |                                     
                                          |                                     
   -0.001   0.      0.          udot      |    vel.grad                         
    0.     -0.001   0.                    |                                     
    0.      0.     0.001                 |                                      
                                          |                                     
    1       1        1           iscau    |    flag for Cauchy                  
            1        1                    |                                     
                     0                    |                                     
                                          |                                     
    0.      0.       0.          scauchy  |    Cauchy stress                    
            0.       0.                   |                                     
                     0.                   @                                     

*******   STEP   1      *** ELASTO-THERMAL CALCULATION         

*******   STEP   2      *** ELASTO-THERMAL CALCULATION         

*******   STEP   3      *** ELASTO-THERMAL CALCULATION         

*******   STEP   4      *** ELASTO-THERMAL CALCULATION         

*******   STEP   5      *** ELASTO-THERMAL CALCULATION         

*******   STEP   6      *** ELASTO-THERMAL CALCULATION         

*******   STEP   7      *** ELASTO-THERMAL CALCULATION         

*******   STEP   8      *** ELASTO-THERMAL CALCULATION         

*******   STEP   9      *** ELASTO-THERMAL CALCULATION         

*******   STEP  10      *** ELASTO-THERMAL CALCULATION         
