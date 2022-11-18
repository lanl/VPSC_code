
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
ex08_Zr\zr_76K_mc.sx                                                                                                              
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
ex08_Zr\IPC20.1                                                                                                                   
0                  -->   20% in-plane tension                                                                                     
ex08_Zr\IPT20.1                                                                                                                   


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
   20   1   0.01    76.   76      nsteps  ictrl  ctrlincr  temp_i  temp_f                                                         
* boundary conditions                                                                                                             
    1       1       1           iudot     |    flag for vel.grad.                                                                 
    1       0       1                     |    (0:unknown-1:known)                                                                
    1       1       0                     |                                                                                       
                                          |                                                                                       
   -0.001   0.      0.          udot      |    vel.grad                                                                           
    0.      0.0008  0.                    |                                                                                       
    0.      0.      0.0002                |                                                                                       
                                          |                                                                                       
    0       0        0           iscau    |    flag for Cauchy                                                                    
            1        0                    |                                                                                       
                     1                    |                                                                                       
                                          |                                                                                       
    0.      0.       0.          scauchy  |    Cauchy stress                                                                      
            0.       0.                   |                                                                                       
                     0.                   @                                                                                       

*******   STEP   1      *** N_EFF=cte CALCULATION              
   76   0.689E-04   0.935E-03   0.465E-04

*******   STEP   2      *** N_EFF=cte CALCULATION              
   54   0.120E-03   0.956E-03   0.112E-03

*******   STEP   3      *** N_EFF=cte CALCULATION              
   42   0.124E-03   0.896E-03   0.112E-03

*******   STEP   4      *** N_EFF=cte CALCULATION              
   34   0.118E-03   0.943E-03   0.908E-04

*******   STEP   5      *** N_EFF=cte CALCULATION              
   26   0.117E-03   0.970E-03   0.138E-03

*******   STEP   6      *** N_EFF=cte CALCULATION              
   19   0.226E-03   0.951E-03   0.373E-03

*******   STEP   7      *** N_EFF=cte CALCULATION              
   15   0.469E-03   0.826E-03   0.536E-03

*******   STEP   8      *** N_EFF=cte CALCULATION              
   16   0.476E-03   0.818E-03   0.485E-03

*******   STEP   9      *** N_EFF=cte CALCULATION              
   16   0.456E-03   0.841E-03   0.471E-03

*******   STEP  10      *** N_EFF=cte CALCULATION              
   21   0.132E-03   0.997E-03   0.183E-03

*******   STEP  11      *** N_EFF=cte CALCULATION              
   24   0.106E-03   0.845E-03   0.199E-03

*******   STEP  12      *** N_EFF=cte CALCULATION              
   26   0.137E-03   0.923E-03   0.330E-03

*******   STEP  13      *** N_EFF=cte CALCULATION              
   24   0.110E-03   0.889E-03   0.141E-03

*******   STEP  14      *** N_EFF=cte CALCULATION              
   25   0.114E-03   0.913E-03   0.278E-03

*******   STEP  15      *** N_EFF=cte CALCULATION              
   26   0.161E-03   0.999E-03   0.458E-03

*******   STEP  16      *** N_EFF=cte CALCULATION              
   26   0.119E-03   0.889E-03   0.283E-03

*******   STEP  17      *** N_EFF=cte CALCULATION              
   25   0.988E-04   0.847E-03   0.156E-03

*******   STEP  18      *** N_EFF=cte CALCULATION              
   25   0.115E-03   0.959E-03   0.233E-03

*******   STEP  19      *** N_EFF=cte CALCULATION              
   26   0.156E-03   0.927E-03   0.234E-03

*******   STEP  20      *** N_EFF=cte CALCULATION              
   24   0.125E-03   0.943E-03   0.293E-03

*******   STEP  21      *** N_EFF=cte CALCULATION              
   24   0.127E-03   0.912E-03   0.282E-03

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
   63   0.831E-04   0.926E-03   0.129E-03

*******   STEP   2      *** N_EFF=cte CALCULATION              
   48   0.853E-04   0.929E-03   0.142E-03

*******   STEP   3      *** N_EFF=cte CALCULATION              
   41   0.862E-04   0.919E-03   0.159E-03

*******   STEP   4      *** N_EFF=cte CALCULATION              
   45   0.993E-04   0.973E-03   0.169E-03

*******   STEP   5      *** N_EFF=cte CALCULATION              
   35   0.117E-03   0.930E-03   0.183E-03

*******   STEP   6      *** N_EFF=cte CALCULATION              
   27   0.143E-03   0.969E-03   0.355E-03

*******   STEP   7      *** N_EFF=cte CALCULATION              
   17   0.287E-03   0.987E-03   0.381E-03

*******   STEP   8      *** N_EFF=cte CALCULATION              
   30   0.107E-03   0.930E-03   0.145E-03

*******   STEP   9      *** N_EFF=cte CALCULATION              
   30   0.102E-03   0.961E-03   0.139E-03

*******   STEP  10      *** N_EFF=cte CALCULATION              
   23   0.159E-03   0.963E-03   0.176E-03

*******   STEP  11      *** N_EFF=cte CALCULATION              
   22   0.156E-03   0.985E-03   0.226E-03

*******   STEP  12      *** N_EFF=cte CALCULATION              
   27   0.103E-03   0.999E-03   0.148E-03

*******   STEP  13      *** N_EFF=cte CALCULATION              
   25   0.923E-04   0.100E-02   0.192E-03

*******   STEP  14      *** N_EFF=cte CALCULATION              
   14   0.132E-03   0.808E-03   0.127E-03

*******   STEP  15      *** N_EFF=cte CALCULATION              
   13   0.394E-03   0.942E-03   0.485E-03

*******   STEP  16      *** N_EFF=cte CALCULATION              
   18   0.224E-03   0.901E-03   0.234E-03

*******   STEP  17      *** N_EFF=cte CALCULATION              
   13   0.348E-03   0.836E-03   0.382E-03

*******   STEP  18      *** N_EFF=cte CALCULATION              
   15   0.129E-03   0.885E-03   0.161E-03

*******   STEP  19      *** N_EFF=cte CALCULATION              
   24   0.109E-03   0.992E-03   0.904E-04

*******   STEP  20      *** N_EFF=cte CALCULATION              
   13   0.426E-03   0.938E-03   0.478E-03

*******   STEP  21      *** N_EFF=cte CALCULATION              
   16   0.782E-04   0.867E-03   0.199E-03

 TIME ELAPSED   83.69secs
