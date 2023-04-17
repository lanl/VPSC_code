*Material: Al5182 with MTS parameters at the end
CUBIC                 crysym
   1.0   1.0   1.0   90.   90.   90.   unit cell axes and angles
Elastic stiffness (single crystal [GPa]) (aluminum 300C, Simmons & Wang)
 106.8    60.7    60.7   000.0   000.0   000.0
  60.7   106.8    60.7   000.0   000.0   000.0
  60.7    60.7   106.8   000.0   000.0   000.0
 000.0   000.0   000.0    28.2   000.0   000.0
 000.0   000.0   000.0   000.0    28.2   000.0
 000.0   000.0   000.0   000.0   000.0    28.2
*Thermal expansion coefficients:
  0.0e-6   0.0e-6   0.0e-6   0.0e0   0.0e0   0.0e0            "alfacc"
*Info about slip & twinning modes in this file:
   1              nmodesx (total # of modes listed in the file)
   1              nmodes  (# of modes to be used in the calculation)
   1              mode(i) (label of the modes to be used)
   <110>(111) SLIP
 1  12  1   0                       modex,nsmx,isensex,itwtypex
   1    1   -1      0    1    1     slip (n & b)
   1    1   -1      1    0    1
   1    1   -1      1   -1    0
   1   -1   -1      0    1   -1
   1   -1   -1      1    0    1
   1   -1   -1      1    1    0
   1   -1    1      0    1    1
   1   -1    1      1    0   -1
   1   -1    1      1    1    0
   1    1    1      0    1   -1
   1    1    1      1    0   -1
   1    1    1      1   -1    0
*Constitutive law
   1      Voce=0, MTS=1
   0      iratesens (0:rate insensitive, 1:rate sensitive)
MTS PARAMETERS FOR Al5182 (from KOK et al.)
  0.5899  28815    0.1194   215        KOVB3 [MPa/K] / MU0 [MPa] / D0 / T0 [K]
  3.333   99.075   0.0      800   2.0  TAUa,TAUi,TAUeini,TH0 [MPa] / KAP
  1.196   1.e7     1.5      0.5        G0i / ED0i [1/s] / QQi / PPi
  1.6     1.e7     1.0      0.666666   G0e / ED0e [1/s] / QQe / PPe
  0.1058  1.e7     665.3               G0esat / EDesat0 [1/s] / TAUesat0 [MPa]
  0.00    16.6                         PSI    /RHO [Mg/m^3]        --> Ta!
  0.1455  9544.e-9 -68.9               Cp1,Cp2,Cp3 [MPa m^3/Mg/K]  --> Ta!
