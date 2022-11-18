
 **** INPUT FILE VPSC.IN FOR THIS RUN ****
1                          number of elements (nelem)                                                                             
1                          number of phases (nph)                                                                                 
1.0  0.0                   relative vol. fract. of phases (wph(i))                                                                
*INFORMATION ABOUT PHASE #1                                                                                                       
0   0   25                    grain shape contrl, fragmentn, crit aspect ratio                                                    
1.0  1.0  1.0                 initial ellipsoid ratios (dummy if ishape=4)                                                        
0.0  0.0  0.0                 init Eul ang ellips axes (dummy if ishape=3,4)                                                      
* name and path of texture file (filetext)                                                                                        
ex08_Zr\zr_5678.tex                                                                                                               
* name and path of single crystal file (filecrys)                                                                                 
ex08_Zr\zr_76K.sx                                                                                                                 
* name and path of grain shape file (dummy if ishape=0) (fileaxes)                                                                
shape1.100                                                                                                                        
* name and path of diffraction file (dummy if idiff=0) (filediff)                                                                 
0                                                                                                                                 
dummy                                                                                                                             
*PRECISION SETTINGS FOR CONVERGENCE PROCEDURES (default values)                                                                   
1.e-3 1.e-3 1.e-3 0.01     errs,errd,errm,errso                                                                                   
100 100 25     itmax:   max # of iter, external, internal and SO loops                                                            
0  2  10  2   irsvar & jrsini,jrsfin,jrstep (dummy if irsvar=0)                                                                   
*INPUT/OUTPUT SETTINGS FOR THE RUN (default is zero)                                                                              
0              irecover:read grain states from POSTMORT.IN (1) or not (0)?                                                        
0              isave:   write grain states in POSTMORT.OUT at step 'isave'?                                                       
0              icubcomp:calculate fcc rolling components?                                                                         
0              nwrite (frequency of texture downloads)                                                                            
*MODELING CONDITIONS FOR THE RUN                                                                                                  
3              interaction (0:FC,1:affine,2:secant,3:neff=10,4:tangent,5:SO)                                                      
1  1  1        iupdate: update orient, grain shape, hardening                                                                     
0              nneigh (0 for no neighbors, 1 for pairs, etc.)                                                                     
0              iflu (0: don't calc, 1: calc fluctuations)                                                                         
*NUMBER OF PROCESSES (Lij const; Lij variable; PCYS ;LANKFORD; rigid rotatn)                                                      
2                                                                                                                                 
*IVGVAR AND PATH\NAME OF FILE FOR EACH PROCESS (dummy if ivgvar=2,3)                                                              
0                  -->   20% in-plane compression                                                                                 
ex08_Zr\IPT20.1                                                                                                                   
2             ivgvar=2 will calculate PCYS                                                                                        
1 2                -->   (S1,S2) projection of stress space                                                                       


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
  1   2   0.1   0.5               isectw,itwlaw,thres1,thres2                                                                     
  {11-22} COMPRESSIVE TWIN --------------------------------------                                                                 
 20                               nrsx                                                                                            
  270.  30.  1000. 178.   0.      tau0x,tau1x,thet0,thet1, hpfac                                                                  
 1.0    1.0    10.0   5.0         hlatex(1,im),im=1,nmodes                                                                        
  1   2   0.1   0.5               isectw,itwlaw,thres1,thres2                                                                     
                                                                                                                                 
**** END OF CRYSTAL DATA FILE ****


 CHECKING THAT CELCC*SELCC-ID6=0   0.2297933E-15
 *********** PHASE   1
 RANDOM PX BULK & POISSON MODULI   95388.889       0.330
 RANDOM PX ELASTIC CTES C11, C12, C44  144073.333   71046.667   36513.333

 INSIDE SUBROUTINE DATA_CRYSTAL
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
DISCRETE TEXTURE FROM ODF FILE  5678.sod                                                                                          
B    1944                                                                                                                         
      5.00      7.07      5.00     0.00028969                                                                                     
      5.00      7.07     15.00     0.00029481                                                                                     
      5.00      7.07     25.00     0.00029881                                                                                     
      5.00      7.07     35.00     0.00030268                                                                                     
      5.00      7.07     45.00     0.00029263                                                                                     
      5.00      7.07     55.00     0.00028680                                                                                     
      5.00     15.79      5.00     0.00081995                                                                                     
      5.00     15.79     15.00     0.00080314                                                                                     
      5.00     15.79     25.00     0.00080906                                                                                     
      5.00     15.79     35.00     0.00081521                                                                                     
      5.00     15.79     45.00     0.00078358                                                                                     
      5.00     15.79     55.00     0.00079817                                                                                     
      5.00     25.46      5.00     0.00102771                                                                                     
      5.00     25.46     15.00     0.00110657                                                                                     
      5.00     25.46     25.00     0.00108154                                                                                     
      5.00     25.46     35.00     0.00108865                                                                                     
    .........................
 **** END OF CRYST TEXTURE DATA FILE ****


*** RUNNING INTERACTION NEFF=   10.00


 UB THERMAL TENSOR (VOIGT NOTATION)
  0.6193E-05  0.6284E-05  0.9318E-05  0.1405E-06  0.2533E-06  0.3489E-07

 UB ELASTIC STIFFNESS (VOIGT NOTATION)
  0.1419E+06  0.7123E+05  0.7033E+05 -0.3970E+02 -0.4075E+02  0.6534E+01
  0.7123E+05  0.1421E+06  0.7030E+05  0.1977E+02 -0.1135E+03 -0.2461E+02
  0.7033E+05  0.7030E+05  0.1508E+06  0.3858E+03  0.8229E+03  0.1084E+03
 -0.3970E+02  0.1977E+02  0.3858E+03  0.3647E+05  0.1318E+03  0.5950E+02
 -0.4075E+02 -0.1135E+03  0.8229E+03  0.1318E+03  0.3645E+05  0.5495E+02
  0.6534E+01 -0.2461E+02  0.1084E+03  0.5950E+02  0.5495E+02  0.3534E+05

 LB THERMAL TENSOR (VOIGT NOTATION)
  0.6354E-05  0.6419E-05  0.8926E-05  0.1177E-06  0.2151E-06  0.2906E-07

 LB ELASTIC STIFFNESS (VOIGT NOTATION)
  0.1415E+06  0.7128E+05  0.7090E+05 -0.1584E+02 -0.8449E+02 -0.1421E+01
  0.7128E+05  0.1416E+06  0.7088E+05 -0.1117E+02 -0.6395E+02 -0.2727E+02
  0.7090E+05  0.7088E+05  0.1482E+06  0.3219E+03  0.6940E+03  0.1025E+03
 -0.1584E+02 -0.1117E+02  0.3219E+03  0.3556E+05  0.1086E+03  0.9214E+01
 -0.8449E+02 -0.6395E+02  0.6940E+03  0.1086E+03  0.3556E+05  0.2441E+02
 -0.1421E+01 -0.2727E+02  0.1025E+03  0.9214E+01  0.2441E+02  0.3510E+05


*** LOAD CONDITIONS FOR THIS RUN
   20   1   0.01    300   300   nsteps  ictrl  ctrlincr  temp_i  temp_f                                                           
* boundary conditions                                                                                                             
    1       1       1           iudot     |    flag for vel.grad.                                                                 
    1       0       1                     |    (0:unknown-1:known)                                                                
    1       1       0                     |                                                                                       
                                          |                                                                                       
    0.001   0.      0.          udot      |    vel.grad                                                                           
    0.     -0.0008  0.                    |                                                                                       
    0.      0.     -0.0002                |                                                                                       
                                          |                                                                                       
    0       0        0           iscau    |    flag for Cauchy                                                                    
            1        0                    |                                                                                       
                     1                    |                                                                                       
                                          |                                                                                       
    0.      0.       0.          scauchy  |    Cauchy stress                                                                      
            0.       0.                   |                                                                                       
                     0.                   @                                                                                       

*******   STEP   1      *** N_EFF=cte CALCULATION              
   63   0.654E-04   0.950E-03   0.426E-04

*******   STEP   2      *** N_EFF=cte CALCULATION              
   71   0.677E-04   0.952E-03   0.466E-04

*******   STEP   3      *** N_EFF=cte CALCULATION              
   60   0.701E-04   0.993E-03   0.474E-04

*******   STEP   4      *** N_EFF=cte CALCULATION              
   51   0.681E-04   0.968E-03   0.461E-04

*******   STEP   5      *** N_EFF=cte CALCULATION              
   43   0.670E-04   0.949E-03   0.441E-04

*******   STEP   6      *** N_EFF=cte CALCULATION              
   35   0.681E-04   0.983E-03   0.610E-04

*******   STEP   7      *** N_EFF=cte CALCULATION              
   32   0.767E-04   0.992E-03   0.623E-04

*******   STEP   8      *** N_EFF=cte CALCULATION              
   28   0.103E-03   0.999E-03   0.107E-03

*******   STEP   9      *** N_EFF=cte CALCULATION              
   26   0.705E-04   0.987E-03   0.130E-03

*******   STEP  10      *** N_EFF=cte CALCULATION              
   25   0.676E-04   0.954E-03   0.809E-04

*******   STEP  11      *** N_EFF=cte CALCULATION              
   22   0.128E-03   0.958E-03   0.132E-03

*******   STEP  12      *** N_EFF=cte CALCULATION              
   23   0.102E-03   0.957E-03   0.537E-04

*******   STEP  13      *** N_EFF=cte CALCULATION              
   24   0.686E-04   0.922E-03   0.688E-04

*******   STEP  14      *** N_EFF=cte CALCULATION              
   23   0.906E-04   0.973E-03   0.463E-04

*******   STEP  15      *** N_EFF=cte CALCULATION              
   21   0.131E-03   0.992E-03   0.175E-03

*******   STEP  16      *** N_EFF=cte CALCULATION              
   17   0.166E-03   0.949E-03   0.173E-03

*******   STEP  17      *** N_EFF=cte CALCULATION              
   18   0.133E-03   0.959E-03   0.126E-03

*******   STEP  18      *** N_EFF=cte CALCULATION              
   21   0.986E-04   0.982E-03   0.159E-03

*******   STEP  19      *** N_EFF=cte CALCULATION              
   18   0.130E-03   0.974E-03   0.134E-03

*******   STEP  20      *** N_EFF=cte CALCULATION              
   17   0.136E-03   0.950E-03   0.222E-03

*******   STEP  21      *** N_EFF=cte CALCULATION              
   18   0.136E-03   0.989E-03   0.172E-03

 ******* STARTS A PCYS CALCULATION *******

 NON CENTRO-SYMMETRIC YIELD SURFACE

NORM OF STRAIN-RATE PROBES IS HARDWIRED TO  0.100E+01

WILL CALCULATE A   1   2   PCYS PROJECTION
 INDICES OF THE COMPONENTS DEFINING STRESS SPACE
    1 = (S22-S11)/SQRT(2)
    2 = (2*S33-S22-S11)/SQRT(6)
    3 = S23*SQRT(2)
    4 = S13*SQRT(2)
    5 = S12*SQRT(2)


*******   STEP   1      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS    1   -0.0000   -1.0000    0.0000    0.0000    0.0000
   17   0.110E-03   0.870E-03   0.217E-03

*******   STEP   2      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS    2   -0.0872   -0.9962    0.0000    0.0000    0.0000
   10   0.124E-03   0.725E-03   0.174E-03

*******   STEP   3      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS    3   -0.1736   -0.9848    0.0000    0.0000    0.0000
   10   0.131E-03   0.681E-03   0.129E-03

*******   STEP   4      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS    4   -0.2588   -0.9659    0.0000    0.0000    0.0000
   10   0.833E-04   0.693E-03   0.200E-03

*******   STEP   5      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS    5   -0.3420   -0.9397    0.0000    0.0000    0.0000
   10   0.952E-04   0.804E-03   0.264E-03

*******   STEP   6      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS    6   -0.4226   -0.9063    0.0000    0.0000    0.0000
   11   0.103E-03   0.680E-03   0.173E-03

*******   STEP   7      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS    7   -0.5000   -0.8660    0.0000    0.0000    0.0000
   12   0.902E-04   0.482E-03   0.123E-03

*******   STEP   8      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS    8   -0.5736   -0.8192    0.0000    0.0000    0.0000
   12   0.106E-03   0.584E-03   0.157E-03

*******   STEP   9      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS    9   -0.6428   -0.7660    0.0000    0.0000    0.0000
   12   0.109E-03   0.681E-03   0.203E-03

*******   STEP  10      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   10   -0.7071   -0.7071    0.0000    0.0000    0.0000
   12   0.111E-03   0.771E-03   0.219E-03

*******   STEP  11      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   11   -0.7660   -0.6428    0.0000    0.0000    0.0000
   12   0.978E-04   0.869E-03   0.214E-03

*******   STEP  12      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   12   -0.8192   -0.5736    0.0000    0.0000    0.0000
   12   0.113E-03   0.980E-03   0.234E-03

*******   STEP  13      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   13   -0.8660   -0.5000    0.0000    0.0000    0.0000
   13   0.938E-04   0.726E-03   0.230E-03

*******   STEP  14      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   14   -0.9063   -0.4226    0.0000    0.0000    0.0000
   13   0.129E-03   0.969E-03   0.256E-03

*******   STEP  15      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   15   -0.9397   -0.3420    0.0000    0.0000    0.0000
   15   0.106E-03   0.808E-03   0.343E-03

*******   STEP  16      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   16   -0.9659   -0.2588    0.0000    0.0000    0.0000
   18   0.928E-04   0.712E-03   0.160E-03

*******   STEP  17      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   17   -0.9848   -0.1736    0.0000    0.0000    0.0000
   28   0.650E-04   0.785E-03   0.108E-03

*******   STEP  18      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   18   -0.9962   -0.0872    0.0000    0.0000    0.0000
   53   0.576E-04   0.873E-03   0.415E-04

*******   STEP  19      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   19   -1.0000    0.0000    0.0000    0.0000    0.0000
   39   0.117E-03   0.753E-03   0.538E-04

*******   STEP  20      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   20   -0.9962    0.0872    0.0000    0.0000    0.0000
   22   0.944E-04   0.845E-03   0.892E-04

*******   STEP  21      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   21   -0.9848    0.1736    0.0000    0.0000    0.0000
   19   0.134E-03   0.737E-03   0.835E-04

*******   STEP  22      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   22   -0.9659    0.2588    0.0000    0.0000    0.0000
   17   0.210E-03   0.889E-03   0.589E-03

*******   STEP  23      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   23   -0.9397    0.3420    0.0000    0.0000    0.0000
   15   0.974E-04   0.791E-03   0.191E-03

*******   STEP  24      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   24   -0.9063    0.4226    0.0000    0.0000    0.0000
   13   0.156E-03   0.849E-03   0.391E-03

*******   STEP  25      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   25   -0.8660    0.5000    0.0000    0.0000    0.0000
   12   0.101E-03   0.509E-03   0.167E-03

*******   STEP  26      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   26   -0.8192    0.5736    0.0000    0.0000    0.0000
   12   0.171E-03   0.698E-03   0.222E-03

*******   STEP  27      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   27   -0.7660    0.6428    0.0000    0.0000    0.0000
   12   0.254E-03   0.933E-03   0.583E-03

*******   STEP  28      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   28   -0.7071    0.7071    0.0000    0.0000    0.0000
   13   0.185E-03   0.866E-03   0.202E-03

*******   STEP  29      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   29   -0.6428    0.7660    0.0000    0.0000    0.0000
   13   0.243E-03   0.971E-03   0.417E-03

*******   STEP  30      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   30   -0.5736    0.8192    0.0000    0.0000    0.0000
   14   0.169E-03   0.905E-03   0.145E-03

*******   STEP  31      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   31   -0.5000    0.8660    0.0000    0.0000    0.0000
   16   0.118E-03   0.753E-03   0.541E-03

*******   STEP  32      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   32   -0.4226    0.9063    0.0000    0.0000    0.0000
   17   0.315E-03   0.969E-03   0.799E-03

*******   STEP  33      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   33   -0.3420    0.9397    0.0000    0.0000    0.0000
   19   0.250E-03   0.568E-03   0.151E-03

*******   STEP  34      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   34   -0.2588    0.9659    0.0000    0.0000    0.0000
   16   0.252E-03   0.686E-03   0.144E-03

*******   STEP  35      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   35   -0.1736    0.9848    0.0000    0.0000    0.0000
   15   0.173E-03   0.679E-03   0.235E-03

*******   STEP  36      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   36   -0.0872    0.9962    0.0000    0.0000    0.0000
   15   0.209E-03   0.737E-03   0.140E-03

*******   STEP  37      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   37    0.0000    1.0000    0.0000    0.0000    0.0000
   11   0.144E-03   0.981E-03   0.984E-03

*******   STEP  38      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   38    0.0872    0.9962    0.0000    0.0000    0.0000
   14   0.238E-03   0.985E-03   0.246E-03

*******   STEP  39      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   39    0.1736    0.9848    0.0000    0.0000    0.0000
   12   0.191E-03   0.837E-03   0.213E-03

*******   STEP  40      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   40    0.2588    0.9659    0.0000    0.0000    0.0000
   13   0.192E-03   0.661E-03   0.780E-03

*******   STEP  41      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   41    0.3420    0.9397    0.0000    0.0000    0.0000
   13   0.276E-03   0.760E-03   0.908E-03

*******   STEP  42      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   42    0.4226    0.9063    0.0000    0.0000    0.0000
   17   0.126E-03   0.758E-03   0.342E-03

*******   STEP  43      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   43    0.5000    0.8660    0.0000    0.0000    0.0000
   20   0.264E-03   0.570E-03   0.685E-03

*******   STEP  44      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   44    0.5736    0.8192    0.0000    0.0000    0.0000
   17   0.314E-03   0.508E-03   0.455E-03

*******   STEP  45      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   45    0.6428    0.7660    0.0000    0.0000    0.0000
   28   0.870E-04   0.676E-03   0.760E-04

*******   STEP  46      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   46    0.7071    0.7071    0.0000    0.0000    0.0000
   20   0.233E-03   0.952E-03   0.192E-03

*******   STEP  47      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   47    0.7660    0.6428    0.0000    0.0000    0.0000
   15   0.163E-03   0.683E-03   0.307E-03

*******   STEP  48      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   48    0.8192    0.5736    0.0000    0.0000    0.0000
   13   0.113E-03   0.815E-03   0.124E-03

*******   STEP  49      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   49    0.8660    0.5000    0.0000    0.0000    0.0000
   12   0.154E-03   0.789E-03   0.113E-03

*******   STEP  50      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   50    0.9063    0.4226    0.0000    0.0000    0.0000
   11   0.188E-03   0.695E-03   0.221E-03

*******   STEP  51      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   51    0.9397    0.3420    0.0000    0.0000    0.0000
   11   0.172E-03   0.524E-03   0.160E-03

*******   STEP  52      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   52    0.9659    0.2588    0.0000    0.0000    0.0000
   12   0.178E-03   0.690E-03   0.101E-03

*******   STEP  53      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   53    0.9848    0.1736    0.0000    0.0000    0.0000
   14   0.168E-03   0.892E-03   0.317E-03

*******   STEP  54      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   54    0.9962    0.0872    0.0000    0.0000    0.0000
   22   0.195E-03   0.711E-03   0.268E-03

*******   STEP  55      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   55    1.0000    0.0000    0.0000    0.0000    0.0000
   54   0.528E-04   0.950E-03   0.388E-04

*******   STEP  56      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   56    0.9962   -0.0872    0.0000    0.0000    0.0000
   27   0.135E-03   0.766E-03   0.751E-04

*******   STEP  57      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   57    0.9848   -0.1736    0.0000    0.0000    0.0000
   16   0.114E-03   0.973E-03   0.200E-03

*******   STEP  58      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   58    0.9659   -0.2588    0.0000    0.0000    0.0000
   15   0.683E-04   0.746E-03   0.361E-03

*******   STEP  59      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   59    0.9397   -0.3420    0.0000    0.0000    0.0000
   13   0.154E-03   0.982E-03   0.220E-03

*******   STEP  60      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   60    0.9063   -0.4226    0.0000    0.0000    0.0000
   12   0.119E-03   0.850E-03   0.183E-03

*******   STEP  61      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   61    0.8660   -0.5000    0.0000    0.0000    0.0000
   11   0.106E-03   0.849E-03   0.151E-03

*******   STEP  62      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   62    0.8192   -0.5736    0.0000    0.0000    0.0000
   10   0.133E-03   0.910E-03   0.247E-03

*******   STEP  63      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   63    0.7660   -0.6428    0.0000    0.0000    0.0000
   10   0.116E-03   0.740E-03   0.217E-03

*******   STEP  64      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   64    0.7071   -0.7071    0.0000    0.0000    0.0000
   10   0.123E-03   0.755E-03   0.225E-03

*******   STEP  65      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   65    0.6428   -0.7660    0.0000    0.0000    0.0000
   10   0.134E-03   0.820E-03   0.256E-03

*******   STEP  66      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   66    0.5736   -0.8192    0.0000    0.0000    0.0000
   11   0.108E-03   0.572E-03   0.159E-03

*******   STEP  67      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   67    0.5000   -0.8660    0.0000    0.0000    0.0000
   11   0.110E-03   0.753E-03   0.236E-03

*******   STEP  68      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   68    0.4226   -0.9063    0.0000    0.0000    0.0000
   11   0.893E-04   0.731E-03   0.166E-03

*******   STEP  69      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   69    0.3420   -0.9397    0.0000    0.0000    0.0000
   10   0.105E-03   0.995E-03   0.233E-03

*******   STEP  70      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   70    0.2588   -0.9659    0.0000    0.0000    0.0000
   10   0.115E-03   0.793E-03   0.185E-03

*******   STEP  71      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   71    0.1736   -0.9848    0.0000    0.0000    0.0000
   10   0.894E-04   0.712E-03   0.234E-03

*******   STEP  72      *** N_EFF=cte CALCULATION              
 step & DBAR in PCYS   72    0.0872   -0.9962    0.0000    0.0000    0.0000
   10   0.892E-04   0.708E-03   0.269E-03

 TIME ELAPSED   74.38secs
