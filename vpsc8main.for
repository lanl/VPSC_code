

C ****************************************************************************
C                       code VPSC version VPSC8                              *
C ****************************************************************************
C  AUTHORS
C ****************************************************************************
C  RICARDO LEBENSOHN                 |  CARLOS TOME                          *
C  T-3 - LANL - MS B216              |  MST-8 - LANL - MS G755               *
C  LOS ALAMOS - NM 87545 - USA       |  LOS ALAMOS - NM 87545 - USA          *
C  e-mail: lebenso@lanl.gov          |  e-mail: tome@lanl.gov                *
C ****************************************************************************
C  GENERAL REFERENCES                                                        *
C ****************************************************************************
C  R.A.Lebensohn & C.N.Tomé, Acta metall mater 41, 2611 (1993)               *
C  C.N.Tomé & R.A.Lebensohn, in "Continuum Scale Simulation of Engineering   *
C    Materials", edited by D. Raabe et al., Wiley-VCH Verlag GmbH, Weinheim, *
C    473-499(2004)                                                           *

C ****************************************************************************
C U.S. GOVERNMENT RIGHTS & COPYRIGHT NOTICE                                  *
C ****************************************************************************
C 2022. Triad National Security, LLC. All rights reserved.                   *

C This program was produced under U.S. Government contract 89233218CNA000001 *
C for Los Alamos National Laboratory (LANL), which is operated by Triad      *
C National Security, LLC for the U.S. Department of Energy/National Nuclear  *
C Security Administration. All rights in the program are reserved by Triad   *
C National Security, LLC, and the U.S. Department of Energy/National Nuclear *
C Security Administration. The Government is granted for itself and others   *
C acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license*
C in this material to reproduce, prepare derivative works, distribute copies *
C to the public, perform publicly and display publicly, and to permit others *
C to do so.                                                                  *
 
C This program is open source under the BSD-3 License. Redistribution and    *
C use in source and binary forms, with or without modification, are          *
C permitted provided that the following conditions are met:                  *
C  1. Redistributions of source code must retain the above copyright notice, *
C     this list of conditions and the following disclaimer.                  *
C  2. Redistributions in binary form must reproduce the above copyright      *
C     notice, this list of conditions and the following disclaimer in the    *
C     documentation and/or other materials provided with the distribution.   *
C  3. neither the name of the copyright holder nor the names of its          *
C     contributors may be used to endorse or promote products derived from   *
C     this software without specific prior written permission.               *

C This software is provided by the copyright holders and contributors "as is"*
C and any express or implied warranties, including, but not limited to, the  *
C implied warranties of merchantability and fitness for a particular purpose *
C are disclaimed. In no event shall the copyright holder or contributors be  *
C liable for any direct, indirect, incidental, special, exemplary, or        *
C consequential damages (including, but not limited to, procurement of       *
C substitute goods or services; loss of use, data, or profits; or business   *
C interruption) however caused and on any theory of liability, whether in    *
C contract, strict liability, or tort (including negligence or otherwise)    *
C arising in any way out of the use of this software, even if advised of the *
C possibility of such damage.                                                *

C ****************************************************************************
C                              # FEATURES #                                  *
C ****************************************************************************
C  - SINGLE PHASE & MULTIPHASE AGGREGATES.                                   *
C  - SELF CONSISTENT TREATMENT WITH MIXED RATE-SENSITIVITY EXPONENTS         *
C  - ESHELBY AND PRESSURE TENSORS FOR FULLY INCOMPRESSIBLE INCLUSION & HEM.  *
C  - MIXED BOUNDARY CONDITIONS ON VELOCITY GRADIENT AND STRESS COMPONENTS.   *
C  - OPTION FOR ANALYZING ROLLING FCC COMPONENTS (icubcom=1)                 *
C  - 'VOCE' HARDENING FOR EACH SYSTEM (ihardlaw=0) OR                        *
C  - 'MTS'  HARDENING LAW FOR EACH SYSTEM (ihardlaw=1) OR                    *
C  - DISLOCATION DENSITY EVOLUTION HARDENING (ihardlaw=20-29) OR             * 
C  - IRRADIATION INDUCED CREEP & GROWTH (ihardlaw=30-32)                     * 
C  - PREDOMINANT TWIN REORIENTATION AND VOLUME FRACTION TRANSFER SCHEMES.    *
C  - GRAIN SHAPE EVOLUTION: OVERALL, PER PHASE, PER INDIVIDUAL GRAIN.        *
C  - INCREMENTAL CONTROL: VON MISES STRAIN, STRAIN COMPONENT, STRESS         *
C    COMPONENT, TIME, TEMPERATURE.                                           *
C  - THERMO-ELASTIC SIMULATION AND SIMULATION OF DIFFRACTION STRAIN.         *
C  - POSTMORTEM CAPABILITY: RECOVERS ARRAYS FROM PREVIOUS RUN                *
C  - POLYCRYSTAL YIELD SURFACE AND LANKFORD COEFFICIENTS CALCULATION.        *
C  - VARIABLE DEFORMATION HISTORY SIMULATION.                                *
C  - ROTATION COUPLING BETWEEN NEIGHBOR ORIENTATIONS                         *
C  - USES Bunge CONVENTION (phi1,Phi,phi2 -> phi,theta,omega) FOR            *
C    EULER ANGLES. BUT READS ALSO Kocks and Roe ANGLES AS INPUT.             *
C                                                                            *
C ****************************************************************************
C                         # INPUT FILES #                                    *
C ****************************************************************************
C ALWAYS NEEDED:                                                             *
C                                                                            *
C       VPSC8.IN     - general input and parameters                          *
C       "FILECRYS"   - elastic moduli,slip & twin modes,CRSS (one per phase) *
C       "FILETEXT"   - initial texture file (one per phase)                  *
C       "FILEPROC"   - deformation process to be simulated                   *
C                                                                            *
C SOMETIMES NEEDED:                                                          *
C                                                                            *
C  if ishape>1     "FILEAXES"  - initial orientation of ellipsoids           *
C                               (morphologic texture, one per phase)         *
C  if irecover=1   POSTMORT.IN - initial state from previous run             *
C  if icubcom=1    CUBCOMP.IN  - ideal orientations - fcc rolling            *
C  if ivgvar=1     "FILEHIST"  - strain history                              *
C  if ihardlaw>=0  hardening law parameters at the end of "FILECRYS"         *
C  if idiff=1      "FILEDIFF"  - diffracting planes and directions           *
C                                                                            *
C ****************************************************************************
C                         # OUTPUT FILES #                                   *
C ****************************************************************************
C DEPENDING ON CASE RUN:                                                     *
C        RERR.OUT     - convergence history                                  *
C        STR_STR_STATS.OUT    - stress components statistics                 *
C        ACT_PHn.OUT  - activity of deformation modes                        *
C        TEX_PHn.OUT  - cryst. texture for each phase/elem at dif steps      *
C        MOR_PHn.OUT  - morph. texture for each phase/elem at dif steps      *
C        STR_STR.OUT  - macroscopic stress-strain components                 *
C        RUN_LOG.OUT  - copy of input files and some control output          *
C        STAT_AXES.OUT - statistics on grain axes                            *
C        DIS_DEN.OUT   - dislocation density evolution in slip modes         *
C     if idiff=1    DIF_PHn.OUT  - internal strains in diffracting planes    *
C     if isave=1    POSTMORT.OUT - final state in grains and PX              *
C     if icubcom=1  CUBCOMPn.OUT - rolling components for fcc                *
C                                                                            *
C   COMMON PARAMETERS, VARIABLES, ARRAYS DEFINED IN VPSC8DIM.FOR             *
C                                                                            *
C ****************************************************************************

      PROGRAM VPSC_MAIN

      USE VPSC8DIM

      CHARACTER FLABEL*2,ITERLBL*40, TLABEL*2
      DIMENSION ROTMAT(3,3)

C ***************************************************************************
C     ASSIGNS # TO UNITS AND OPENS I/O FILES.

      UR0= 9     ! VPSC8.IN      (OPEN/CLOSE IN MAIN)
      UR1= 1     ! FILECRYS      (OPEN/READ/CLOSE IN VPSC_INPUT)
      UR2= 2     ! FILETEXT      (OPEN/READ/CLOSE IN VPSC_INPUT)
      UR3= 3     ! FILEAXES      (OPEN/READ/CLOSE IN VPSC_INPUT)
      UR6= 99    ! FILEHIST      (OPEN/CLOSE IN MAIN)
C     UR7= 7     ! FILEDIFF      (OPEN IN DIFF_PLANES)
      UR4= 4     ! POSTMORT.IN   (OPEN/CLOSE IN MAIN)
      UR5= 98    ! CUBCOMP.IN    (OPEN/CLOSE IN MAIN)

C         50+    ! ACT_PHn.OUT      (OPEN IN MAIN)
C         60+    ! CUBCOMPn.OUT     (OPEN IN MAIN)
C         20-22  ! DIF_STR_PHn.OUT  (OPEN IN DIFF_PLANES IF IDIFF=1)
C         23-25  ! DIF_WGT_PHn.OUT  (OPEN IN DIFF_PLANES IF IDIFF=1)
C         59     ! DIS_DEN.OUT      (OPEN IN UPDATE_CRSS_DD)
C         777    ! EL-TH-MODULI.OUT (OPEN IN ELSC)
C         83     ! FLUCT.OUT        (OPEN IN VPSC_INPUT IF IFLU=1)
C         84     ! FLCUB.OUT        (OPEN IN VPSC_INPUT IF INTERACT=5)
C         15     ! LANKFORD.OUT     (OPEN IN MAIN)
C         40+    ! MOR_PHn.OUT      (OPEN IN MAIN)
C         45+    ! STAT_AXES_PHn.OUT  (OPEN IN MAIN)
C         14     ! PCYS.OUT         (OPEN IN MAIN)
C         141    ! PCYSwRATES.OUT   (OPEN IN PCYS)
C         142    ! PCYSwNORMALS.OUT (OPEN IN PCYS)
      UW2=19     ! POSTMORT.OUT     (OPEN/CLOSE IN MAIN IF ISAVE=n)
      UW1=12     ! RERR.OUT
C         10     ! RUN_LOG.OUT      (OPEN IN MAIN)
C         97     ! SO.OUT           (OPEN IN VPSC_INPUT IF INTERACT=5)
C         13     ! STR_STR.OUT      (OPEN IN MAIN)
C         11     ! STR_STR_STATS.OUT  (OPEN IN MAIN)
C         30+    ! TEX_PHn.OUT      (OPEN IN MAIN)

C ***************************************************************************
      CALL CPU_TIME(START_TIME)
C ***************************************************************************

      OPEN(10, FILE='RUN_LOG.OUT',STATUS='UNKNOWN')
      OPEN(11, FILE='STR_STR_STATS.OUT'  ,STATUS='UNKNOWN')
      OPEN(13, FILE='STR_STR.OUT',STATUS='UNKNOWN')
      OPEN(UW1,FILE='RERR.OUT'   ,STATUS='UNKNOWN')

C ***************************************************************************
C     CALL SUBROUTINE VPSC_INPUT FOR READING CRYSTAL, GRAIN AND TEXTURE DATA.
C     INITIALIZES ARRAYS.

      OPEN(UR0,FILE='VPSC8.in',STATUS='OLD')

      CALL VPSC_INPUT

      CALL ELSC (0)     ! CALC UB ELASTIC & THERMAL MODULI
      CALL ELSC (1)     ! CALC SC ELASTIC & THERMAL MODULI & LOCALIZATION TENSORS
     
C ***************************************************************************
C     READS SX & PX STATE FROM PREVIOUS RUN (sav,xmsec,sg,ag,crss,..)
C     WHEN IRECOVER=1.

      IF(IRECOVER.EQ.1) THEN
        OPEN(UR4,file='POSTMORT.IN',form='UNFORMATTED',
     #           access='SEQUENTIAL',status='OLD')
          CALL POSTMORTEM (1)
        CLOSE(UNIT=UR4)
      ENDIF

      IF(ICUBCOM.EQ.1) THEN
        OPEN(UR5,file='CUBCOMP.IN',status='OLD')
          CALL CUBCOMP(0,0)          ! reads ideal rolling components
        CLOSE(UNIT=UR5)
      ENDIF

      
C ***********************************************************************
C     OPEN OUTPUT FILES
C ***********************************************************************

      IF(NELEM.EQ.1) THEN
        NFILES=NPH
        FLABEL='PH'
      ELSE IF(NELEM.GT.1) THEN
        NFILES=NELEM
        FLABEL='EL'
      ENDIF

      DO I=1,NFILES
        IF(IVGVAR.LE.1) THEN     ! may be redundant
          IUNIT=30+I
          OPEN(IUNIT,FILE='TEX_'//FLABEL//CHAR(48+I)//'.OUT',
     #        STATUS='UNKNOWN')
        ENDIF

        IF(ISHAPE(I).NE.0) THEN
          IUNIT=40+I
          OPEN(IUNIT,FILE='MOR_'//FLABEL//CHAR(48+I)//'.OUT',
     #        STATUS='UNKNOWN')
          IUNIT=45+I
          OPEN(IUNIT,FILE='STAT_AXES_'//FLABEL//CHAR(48+I)//'.OUT',
     #        STATUS='UNKNOWN')
        ENDIF

        IF(ICUBCOM.EQ.1) THEN
          IUNIT=60+I
          OPEN(IUNIT,FILE='CUBCOMP'//CHAR(48+I)//'.OUT',
     #        STATUS='UNKNOWN')
        ENDIF

        IUNIT=50+I
        OPEN(IUNIT,FILE='ACT_'//FLABEL//CHAR(48+I)//'.OUT',
     #        STATUS='UNKNOWN')

       ENDDO         ! end of DO I=1,NFILES

c     cubcomp analysis for initial texture
      IF(ICUBCOM.EQ.1) CALL CUBCOMP (0,1)

C *******************************************************************
C     DO LOOP OVER PROCESSES
C *******************************************************************

      TIME=0.
      CONVERGED=.TRUE.   ! indicates if convergence tolerance is not reached
	  
C *** SET TO ZERO ACCUMULATED STRAIN COMPONENTS
      EPSACU=0.
      EPSVM =0.
      EPSTOTc(:,:)=0.
      EPSTOTv(:)  =0.

      EELGR(:,:)=0.
      EELAV(:)  =0.
      ETHAV(:)  =0.
	  
C *** READS NUMBER OF PROCESSES TO RUN SEQUENTIALLY. IT CAN BE ANY
C     COMBINATION OF UNIFORM LOAD (IVGVAR=0), VARIABLE LOAD (IVGVAR=1),
C     PCYS (IVGVAR=2), LANKFORD (IVGVAR=3), RIGID ROTATION (IVGVAR=4)

      READ(UR0,'(A)') PROSA
      READ(UR0,*)     NPROC
      READ(UR0,'(A)') PROSA

      DO 3000 IPROC=1,NPROC

C *** READS LOAD CONDITIONS: STRESS, VELOCITY GRADIENT, TEMP, INCR, NSTEPS.
C *** CALCULATES SYMMETRIC STRAIN RATE 'DBARc'

      READ(UR0,*) IVGVAR

      IPCYSOPEN=0
      IPCYSSKIP=0
      ILANKOPEN=0

      IF(IVGVAR.EQ.0) THEN
        READ(UR0,'(A)') FILEHIST
        OPEN(UR6,file=FILEHIST,status='OLD')
        CALL LOAD_CONDITIONS (UR6)
        CLOSE(UNIT=UR6)
C --> EXTRA STEP REQUIRED FOR UPDATING MACRO STRESS AT FINAL STRAIN
        NSTEPSX=NSTEPS+1      
        IF(INTERACTION.EQ.-1) THEN
C --> FOR ELASTIC LOADING MACRO STRESS CORRESPONDS TO END OF THE STEP 
          NSTEPSX=NSTEPS 
        ENDIF
      ELSE IF(IVGVAR.EQ.1) THEN
        READ(UR0,'(A)') FILEHIST
        OPEN(UR6,file=FILEHIST,status='UNKNOWN')
        CALL VAR_VEL_GRAD(0)
        NSTEPSX=NSTEPS
      ELSE IF(IVGVAR.EQ.2) THEN
        IF(IPCYSOPEN.EQ.0) THEN
          OPEN(14 ,file='PCYS.OUT',status='UNKNOWN')
          IPCYSOPEN=1
        ENDIF
        READ(UR0,*) INDX,INDY              ! components of 2D subspace
        CALL PCYS(IDUM,INDX,INDY,0)        ! generates 2D-PCYS probes
        NSTEPSX=NSTEPS
      ELSE IF(IVGVAR.EQ.-2) THEN
        IF(IPCYSOPEN.EQ.0) THEN
          OPEN(14,file='PCYS.OUT',status='UNKNOWN')
          IPCYSOPEN=1
        ENDIF
        CALL PCYS_IT(IDUM,IPCYSSKIP,0)     ! generates 5D-PCYS probes
        NSTEPSX=NSTEPS
      ELSE IF(IVGVAR.EQ.3) THEN
        IF(ILANKOPEN.EQ.0) then
          OPEN(15,file='LANKFORD.OUT',status='UNKNOWN')
          ILANKOPEN=1
        ENDIF
        READ(UR0,*) DELTALANK              ! angular increment for probes
        CALL LANKFORD(ISTEP,DELTALANK,0)   ! initializes arrays
        NSTEPSX=NSTEPS
      ELSE IF(IVGVAR.EQ.4) THEN
        READ(UR0,'(A)') FILEHIST
        OPEN(UR6,file=FILEHIST,status='OLD')
          READ(UR6,'(A)') PROSA
          READ(UR6,*) ((ROTMAT(I,J),J=1,3),I=1,3)
        CLOSE(UNIT=UR6)
        CALL TEXTURE_ROTATION(ROTMAT)
        CALL WRITE_TEXTURE
        GO TO 3000
      ENDIF

      IF(IVGVAR.LE.1 .AND. NSTEPS.LT.NWRITE) THEN
        WRITE(*,'(/,'' WARNING *** NWRITE='',I3,''  > NSTEPS='',I3,
     #        '' --> WILL NOT WRITE TEXTURE !!!'',/)') NWRITE,NSTEPS
        print *, 'enter c to continue'
        read  *
      ENDIF

C *******************************************************************
C     DO LOOP OVER DEFORMATION STEPS
C *******************************************************************

      IF(INTERACTION.EQ.-1) ITERLBL=' *** THERMO-ELASTIC CALCULATION'
      IF(INTERACTION.EQ. 0) ITERLBL=' *** FULL CONSTRAINT CALCULATION'
      IF(INTERACTION.EQ. 1) ITERLBL=' *** AFFINE CALCULATION'
      IF(INTERACTION.EQ. 2) ITERLBL=' *** SECANT CALCULATION'
      IF(INTERACTION.EQ. 3) THEN
        IF(.NOT.RDCrun) ITERLBL=' *** N_EFF=cte CALCULATION'
        IF     (RDCrun) ITERLBL=' *** N_EFF FROM RDC CALCULATION'
      ENDIF
      IF(INTERACTION.EQ. 4) ITERLBL=' *** TANGENT CALCULATION'
      IF(INTERACTION.EQ. 5) ITERLBL=' *** SECOND ORDER CALCULATION'

      DO 2500 ISTEP=1,NSTEPSX

        NWRITEX=ISTEP      ! used inside SO_SDPX

        if(iflu.eq.1) then
          write(83,'(a,i7)') ' STEP =',istep
          if(icubcom.eq.1) write(84,'(a,i7)') ' STEP =',istep
        endif

        if(interaction.eq.5) then
          write(97,'(a,i7)') ' STEP =',istep
          write(97,'(a)') '    IRS  ITSO  ERRASO      ERRESO'
        endif

        WRITE(UW1,'(/,''*******   STEP'',I4,5X,A)') ISTEP,ITERLBL
        WRITE( 10,'(/,''*******   STEP'',I4,5X,A)') ISTEP,ITERLBL
        WRITE(  *,'(/,''*******   STEP'',I4,5X,A)') ISTEP,ITERLBL
        WRITE(  *,'(''   ITSGR    SIGGR      SIGAV        DAV'')')

CC        WRITE(  *,'(''   ITSGR    SIGGR      SIGAV        DAV'',
CC     #             ''    ITTAN     MTAN    ITSEC     MSEC'')')

C *** IMPOSE VELOCITY GRADIENT AT EACH STEP:
C       IF IMPOSING THERMO-ELASTIC MONOTONIC LOADING    (IVGVAR=-1)
C       IF IMPOSING VISCO-PLASTIC MONOTONIC LOADING     (IVGVAR= 0)
C       IF IMPOSING A VISCO-PLASTIC HISTORY             (IVGVAR= 1)
C       IF PROBING THE POLYCRYSTAL YIELD SURFACE        (IVGVAR= 2)
C       IF CALCULATING A 5D POLYCRYSTAL YIELD SURFACE   (IVGVAR=-2)
C       IF PROBING FOR LANKFORD COEFFICIENTS            (IVGVAR= 3)

C ***********************************************************************
C *** FOR THERMO-ELASTIC LOADING SKIP PLASTICITY-RELATED CALCULATIONS

        IF(INTERACTION.EQ.-1) GO TO 1800   
C ***********************************************************************

        IF(IVGVAR.EQ.1)  CALL VAR_VEL_GRAD(1)
        IF(IVGVAR.EQ.2)  CALL PCYS(ISTEP,INDX,INDY,1)      ! 2D PCYS
        IF(IVGVAR.EQ.-2) THEN                              ! 5D PCYS
          CALL PCYS_IT(ISTEP,IPCYSSKIP,1)
          IF(IPCYSSKIP.EQ.1) GO TO 2500       ! SKIP ALREADY CALCULATED STATE
        ENDIF
        IF(IVGVAR.EQ.3) CALL LANKFORD(ISTEP,DELTALANK,1)

C ***********************************************************************
C *** RECALCULATES SCHMID TENSORS OF EACH SYSTEM EXPRESSED IN SAMPLE AXES
 
        CALL UPDATE_SCHMID
C ***********************************************************************
C *** IF IRATESENS=0 SCALES GAMD0 TO ELIMINATE RATE SENSITIVITY INDUCED BY 'n'

        IF(IHARDLAW.EQ.0) THEN
          IF(IRATESENS.EQ.0) CALL SCALE_GAMD0S (ISTEP)   ! rate-insensitive response
        ENDIF
C ***********************************************************************
C     FOR MTS SCALE GAMD0 TO ELIMINATE RATE SENSITIVITY INDUCED BY 'n'.
C     CALCULATE CRRS's ASSOCIATED WITH TEMPERATURE AND RATE FOR THIS STEP

        IF(IHARDLAW.EQ.1) THEN      ! compatibilize with ihardlaw=20,23
          CALL SCALE_GAMD0S (ISTEP)
          EDOT=SQRT(2./3.)*VNORM(DBAR,5)      ! MTS uses von Mises norm
          CALL UPDATE_CRSS_MTS (EDOT,2)
        ENDIF
C **********************************************************************
C *** FOR DISLOC DENSITY LAW MAKES GAMD0=NORM OF THE RATE TENSOR IMPOSED
C     IN THE STEP TO ELIMINATE RATE SENSITIVITY ASSOCIATED WITH 'n'.
C     CALCULATES CRRS's ASSOCIATED WITH PROCESS TEMPERATURE FOR 1st STEP.
C *** IOPTION=1 RESETS DISLOC DENSITIES & CRSS, WHICH DOES NOT CORRESPOND 
C     WHEN SIMULATING PATH CHANGES --> USE ONLY FOR STEP #1 OF PROCESS #1.

        IF(IHARDLAW.EQ.20) THEN    ! what about MTS & DD_REV?
          CALL SCALE_GAMD0S (ISTEP)                 ! makes Gamd0s=norm(DBAR)
          IF(IPROC.EQ.1 .and. ISTEP.EQ.1) THEN
C --> quick&dirty fix when 1st process is a PCYS or LANKFORD calculation and
C --> the process temperature has not been read from process file.
C --> Besides, PCYS equidisipation procedure relies on n-power law and is not
C --> applicable to DD or MTS constitutive laws.
            IF(IVGVAR.EQ.2 .OR. IVGVAR.EQ.3) TEMP=300.
		    CALL UPDATE_CRSS_DD (1)
          ENDIF
        ENDIF

        IF(IHARDLAW.EQ.23) THEN    
          CALL SCALE_GAMD0S (ISTEP)                 ! makes Gamd0s=norm(DBAR)
        ENDIF

C **********************************************************************        
C *** UPDATES TRANSFORMATION TENSOR DG_TRANS(I,KKK)
C *** IN THE CASE OF IRRADIATION GROWTH IT IS THE GROWTH TENSOR,  
C     FUNCTION OF THE EVOLVING LOOP DENSITY 

        IF(IHARDLAW.EQ.30) CALL UPDATE_GROWTH_RATE (2)          ! CNT
        IF(IHARDLAW.EQ.31) CALL UPDATE_GROWTH_ZIRC2(2,IDUMMY)   ! AP:
        IF(IHARDLAW.EQ.32) CALL UPDATE_GROWTH_ZR   (2,IDUMMY)   ! AP:
        
C **********************************************************************
C     IF IMPOSING STRAIN MAKES A TAYLOR GUESS FOR THE FIRST STEP WHEN
C     IRECOVER=0 OR STARTS FROM RECOVERED STRESS STATE 'SG' WHEN IRECOVER=1.
C     IF IMPOSING STRESS MAKES A SACHS GUESS FOR EVERY STEP.

        IF(ISTEP.EQ.1.AND.IRECOVER.EQ.0 .OR. STRAIN_CONTROL.EQ.0) THEN
          CALL INITIAL_STATE_GUESS
        ENDIF
        IF(ISTEP.EQ.1.AND.IRECOVER.EQ.1) THEN
          KGX=1
          DO IPH=IPHBOT,IPHTOP
            IPHEL=IPH-IPHBOT+1
          DO KKK=NGR(IPH-1)+1,NGR(IPH)
            CALL GRAIN_RATE_AND_MODULI (1,KGX,KKK,IPHEL,IPH)
            KGX=KGX+1
          ENDDO
          ENDDO
        ENDIF
		
C **********************************************************************
 1800 CONTINUE      ! 1800 
C ***************************************************************************
C     THIS IS THE CORE OF THE CODE !!
C     * VISCO PLASTIC FULL CONSTRAINT CALCULATION (TAYLOR) IF INTERACTION=0
C     * VISCO PLASTIC SELF CONSISTENT CALCULATION IF INTERACTION > 0
C     * ELASTIC-THERMAL SELF CONSISTENT CALCULATION IF INTERACTION = -1

      IF(INTERACTION.EQ.-1) CALL ELSC (2)     ! UPDATES ELASTIC STRAIN/STRESS IN GRAINS

      IF(INTERACTION.EQ.0) CALL VPFC (ISTEP)

      IF(INTERACTION.GT.0) CALL VPSC (ISTEP)

C ***************************************************************************
C *** MACROSCOPIC RIGID ROTATION AND VELOCITY GRADIENT UPDATE
C        L_ij = D_ij + W_ij    -->    LBAR = DBAR+WBAR   
C     SUBROUTINES VPFC & VPSC PROVIDE THE FULL MACROSCOPIC STRAIN RATE 'D_ij'  
C     ACCOUNTING FOR BC's ENFORCED ON MACROSCOPIC STRESS COMPONENTS.

C     CALCULATES MACROSCOPIC ROTATION RATE 'WBAR' USING  W_ij=L_ij-D_ij
C     --> straightforward when all Lij components are enforced
C     --> need to use boundary conditions 'ILBAR' on 'L_ij' otherwise
C     --> off-diagonal ILBAR_ij & ILBAR_ji cannot be zero simultaneously 
C         because macroscopic rigid rotation has to be enforced externally

      DO I=1,3
       WBARc(I,I)=0.      ! TAKES CARE OF DIAGONAL ELEMENTS
        DO J=1,3
          IF(ILBAR(I,J).EQ.1 .AND. ILBAR(J,I).EQ.1) THEN
            WBARc(I,J)=(LIJBARc(I,J)-LIJBARc(J,I))/2.
          ENDIF
          IF(ILBAR(I,J).EQ.1 .AND. ILBAR(J,I).EQ.0) THEN
            WBARc(I,J)=LIJBARc(I,J)-DBARc(I,J)
            WBARc(J,I)=-WBARc(I,J)
          ENDIF
          LIJBARc(I,J)=DBARc(I,J)+WBARc(I,J)
        ENDDO
      ENDDO
C *************************************************************************

C     IF INTERACTION=0 : DBAR IS IMPOSED AND SBAR=SAV.
C     IF INTERACTION>0 : DBAR & SBAR FOLLOW FROM THE CALL TO SUBR STATE6x6
C        INSIDE SUBR VPSC.
C     IF INTERACTION=-1: DBAR & SBAR FOLLOW FROM THE CALL TO SUBR STATE6x6
C        INSIDE SUBR ELSC.
C     VON MISES STRAIN-RATE & STRESS: FUNCTION OF DEVIATORIC COMPONENTS ONLY

      CALL VOIGT(DBARv,DBARc,AUX66,AUX3333,2)
      CALL VOIGT(SBARv,SBARc,AUX66,AUX3333,2)
      SVM=0.
      DVM=0.
      DO I=1,5
        SVM=SVM+SBAR(I)*SBAR(I)
        DVM=DVM+SBAR(I)*DBAR(I)
      ENDDO
      SVM=SQRT(SVM*3./2.)
      DVM=ABS(DVM)             ! TO AVOID <0 WHEN GROWTH DOMINATES CREEP
      IF(SVM.NE.0.) DVM=DVM/SVM       ! WORK CONJUGATE DEFINITION OF DVM
      IF(SVM.EQ.0.) DVM=SQRT(2./3.)*VNORM(DBAR,5)

C *******************************************************************
C     WRITES STRESS-STRAIN AND STATISTICS FOR THE STEP
C *******************************************************************

      IF(INTERACTION.EQ.-1) THEN      ! RELEVANT ONLY TO ELASTICITY
        DO I=1,3
        DO J=1,3
          EPSTOTc(I,J)=EPSTOTc(I,J)+DBARc(I,J)*TIME_INCR
        ENDDO
        ENDDO
        CALL VOIGT(EPSTOTv,EPSTOTc,AUX66,AUX3333,2)
        IF(ICTRL.LE.6) EPSACU=EPSTOTv(ICTRL)
        IF(ICTRL.EQ.8) EPSACU=VNORM(EPSTOTv,6)
        EPSVM=EPSVM+DVM*TIME_INCR      ! necessary? used?
        TEMP=TEMP+TEMP_INCR
        CALL STAT_STRESS_STRAIN
          CALL WRITE_STRESS_STRAIN (ISTEP)

C *** IF IDIFF=1 CALCULATES INTERNAL STRAIN ON DIFFRACTING PLANES USING 
C     'SG(1:6,KGR)' CALCULATED INSIDE SUBROUT ELSC.
        DO IPH=IPHBOT,IPHTOP             
          IF(IDIFF(IPH).EQ.1) CALL DIFF_PLANES(IPH,ICRYSYMPH(IPH),1)
        ENDDO

        GO TO 2500                ! GOES TO THE END OF INCREMENTAL STEP
      ENDIF
	  
      IF(INTERACTION.NE.-1) THEN      ! RELEVANT ONLY TO PLASTICITY

        IF(IVGVAR.EQ.-2) THEN                              ! 5D PCYS
          CALL PCYS_IT (ISTEP,IPCYSSKIP,2)
          GO TO 2500              ! GOES TO THE END OF INCREMENTAL STEP
        ENDIF
		
        IF(IVGVAR.LE.1) THEN
          ISHAPESUM=0
		  DO IPH=1,NPH
		    ISHAPESUM=ISHAPESUM+ISHAPE(IPH)
		  ENDDO
          IF(ISHAPESUM.NE.0) CALL STAT_GRAIN_SHAPE (ISTEP)
          IF(ICUBCOM.EQ.1.AND.ISTEP.GT.1) CALL CUBCOMP (ISTEP,1)
		ENDIF
        IF(IVGVAR.EQ.2) CALL PCYS (ISTEP,INDX,INDY,2)      ! 2D PCYS
        IF(IVGVAR.EQ.3) CALL LANKFORD (ISTEP,DELTALANK,2)

        CALL STAT_SHEAR_ACTIVITY
          CALL WRITE_SHEAR_ACTIVITY(ISTEP)
        CALL GRAIN_INFO      ! CALCULATES GRAIN & PX TAYLOR FACTOR & PLASTIC WORK
        CALL STAT_STRESS_STRAIN
		  CALL WRITE_STRESS_STRAIN (ISTEP)

C *** CALCULATES INTERNAL STRAIN ON DIFFRACTING PLANES IF IDIFF.NE.0 
C *** IN THE VISCO-PLASTIC CASE USES 'SBAR' TO CALCULATE 'SG(6,KGR)'
C     USING THE ELASTIC LOCALIZATION EQUATION (APPROXIMATE !).

        DO IPH=IPHBOT,IPHTOP             
          IF(IDIFF(IPH).EQ.1) CALL DIFF_PLANES(IPH,ICRYSYMPH(IPH),1)
        ENDDO

        IF(IVGVAR.GE.2) GO TO 2500       ! GOES TO THE END OF INCREMENTAL STEP

        IF(ISTEP.EQ.NSTEPSX) GO TO 2000  ! IN LAST VP STEP ONLY UPDATES STRESS
		
      ENDIF      ! END OF IF(INTERACTION.NE.-1)

C *********************************************************************
C     GIVEN THE OVERALL STRAIN RATE AND STRESS THEN:
C     *  IF ICTRL=0 CALCULATES TIME INCREMENT NECESSARY TO ACHIEVE
C        THE IMPOSED VON MISES STRAIN INCREMENT.
C     *  IF ICTRL=1,6 AND STRAIN_CONTROL=1 IT CORRESPONDS TO AN
C        IMPOSED STRAIN COMPONENT AND STRAIN INCREMENT.  FROM KNOWN
C        STRAIN RATE CALCULATES TIME INCREMENT REQUIRED TO ACHIEVE IT.
C     *  IF ICTRL=1,6 AND STRAIN_CONTROL=0 IT CORRESPONDS TO A CREEP
C        TEST WITH IMPOSED STRESS COMPONENT AND TIME INCREMENT.  FROM KNOWN
C        STRAIN RATE CALCULATES THE STRAIN INCREMENT THAT TAKES PLACE.
C     *  IF ICTRL=8 IT CORRESPONDS TO A THERMO-ELASTIC CASE WHERE THE
C        TEMPERATURE INCREMENT IS IMPOSED. THE PROCESS IS CONTROLLED 
C        INSIDE SUBROUTINE ELSC.
C ***********************************************************************

      IF(ICTRL.EQ.0) THEN
        if(IVGVAR.NE.1) TIME_INCR =EVMINCR/DVM
        if(IVGVAR.EQ.1) EVMINCR   =DVM*TIME_INCR      ! VARIABLE VEL GRAD CASE
        EPSACU =EPSACU+EVMINCR
      ELSE IF(ICTRL.GE.1.AND.ICTRL.LE.6 .AND. STRAIN_CONTROL.EQ.1) THEN
        TIME_INCR =EIJINCR/ABS(DBARCTRL)
        EPSACU =EPSACU+EIJINCR
      ELSE IF(ICTRL.GE.1.AND.ICTRL.LE.6 .AND. STRAIN_CONTROL.EQ.0) THEN
        EVMINCR=DVM*TIME_INCR
        EPSACU =EPSACU+EVMINCR
      ENDIF

      EPSVM=EPSVM+DVM*TIME_INCR
      TEMP=TEMP+TEMP_INCR
      TIME=TIME+TIME_INCR
      DO I=1,3
      DO J=1,3
        EPSTOTc(I,J)=EPSTOTc(I,J)+DBARc(I,J)*TIME_INCR
      ENDDO
      ENDDO
      CALL VOIGT(EPSTOTv,EPSTOTc,AUX66,AUX3333,2)

c --> consider defining here Green strain tensor: E=(F_t:F-I)/2
 
C ***********************************************************************
C     UPDATES ORIENTATION, HARDENING AND SHAPE OF EVERY GRAIN USING A
C     FORWARD EXTRAPOLATION: 'CALCULATED RATE' x 'TIME_INCR'.
C     REORIENTS GRAINS BY TWINNING AND MAKES TWINNING STATISTICS.
C     SKIPS UPDATE FOR PCYS OR LANKFORD RUN.
C ***********************************************************************

        DO IPH=IPHBOT,IPHTOP
          IPHEL=IPH-IPHBOT+1
          IF(NTWMOD(IPHEL).NE.0) THEN
             IF(ITWINLAW.LE.2) THEN
			   CALL UPDATE_TWINNING_PTR (IPH)  
               CALL WRITE_TWIN_STATS (ISTEP,IPH)            ! WRITES TWINNING STATS
             ELSE IF (ITWINLAW.EQ.3 .AND. IPH.EQ.1)  THEN
			   CALL UPDATE_TWINNING_VFT           ! INTERNALLY LOOPS OVER PH1 & PH2
               CALL WRITE_TWIN_STATS (ISTEP,IPH)            ! WRITES TWINNING STATS
			 ENDIF
           ENDIF
        ENDDO
              
        IF(IUPDORI.EQ.1) CALL UPDATE_ORIENTATION     ! UPDATES CRYSTALLOGR ORIENT

        CALL UPDATE_FIJ(0)                   ! UPDATES AVERAGE DEFORMATION TENSOR
        CALL UPDATE_SHAPE(0)   ! UPDATES AVERAGE ELLIPSOID SHAPE
        DO IPH=IPHBOT,IPHTOP
          CALL UPDATE_FIJ(IPH)          ! UPDATES DEFORM TENSOR OF PHASE & GRAINS
          IF(IUPDSHP.EQ.1) CALL UPDATE_SHAPE(IPH)    ! UPDATES SHAPE OF PH & GRNS
        ENDDO

        IF(IUPDHAR.EQ.1) THEN

C *** VOCE HARDENING PLUS PREDOMINANT TWIN REORIENTATION SCHEME
          IF(IHARDLAW.EQ.0)  CALL UPDATE_CRSS_VOCE (2,IDUMMY)
C *** MECHANICAL THRESHOLD STRESS HARDENING (NO TWINNING, ONLY ONE PHASE)
          IF(IHARDLAW.EQ.1) THEN
            EDOT=SQRT(2./3.)*VNORM(DBAR,5)
            CALL UPDATE_CRSS_MTS (EDOT,3)
          ENDIF
C *** DISLOCATION DENSITY STRESS HARDENING (ONLY ONE PHASE)
          IF(IHARDLAW.EQ.20) CALL UPDATE_CRSS_DD (2)
          IF(IHARDLAW.EQ.23) CALL UPDATE_CRSS_DD_REV (2)

        ENDIF

 2000 CONTINUE      ! 2000

C *******************************************************************
C     END OF FORWARD UPDATE OF STRAIN, SHAPE, ORIENTATION & CRSS
C *******************************************************************

C *******************************************************************
C     WRITES TEXTURE FILES FOR EACH PHASE

      IWRITE=0
      IF(NWRITE.EQ.0 .AND. ISTEP.EQ.NSTEPS) IWRITE=1
      IF(NWRITE.NE.0) THEN
        IF(MOD(ISTEP,NWRITE).EQ.0) IWRITE=1
      ENDIF
 
      IF(IWRITE.EQ.1) THEN
        IF(IVGVAR.EQ.0) CALL WRITE_TEXTURE
        ICOROT=0      ! hardwires the non-corotational option (default)
        ICOROT=1      ! hardwires the corotational option
c *** non co-rotational algorithm
        IF(IVGVAR.EQ.1 .AND. ICOROT.EQ.0) CALL WRITE_TEXTURE 
c *** co-rotational algorithm
        IF(IVGVAR.EQ.1 .AND. ICOROT.EQ.1) THEN         
          write (*,'('' XXXXXXXXXXX   calling texture_rotation'')')
          CALL TEXTURE_ROTATION(ROT_PD)      ! PUT TEXTURE IN SAMPLE AXES
          CALL WRITE_TEXTURE
          CALL TEXTURE_ROTATION(ROT_PD_t)    ! PUT TEXTURE BACK IN CO-VARIANT AXES
c *** RESET TEXTURE TO SAMPLE AXES AT THE END OF THE VAR_VEL_GRAD PROCESS
          IF(ISTEP.EQ.NSTEPS) CALL TEXTURE_ROTATION(ROT_PD)
        ENDIF
      ENDIF

C *******************************************************************
C     WRITES GRAIN'S & PX STATES INTO BINARY FILE 'POSTMORT.OUT'

      IF(ISAVE.EQ.ISTEP) THEN
        OPEN(UW2,file='POSTMORT.OUT',form='UNFORMATTED',
     #           access='SEQUENTIAL',status='UNKNOWN')
        CALL POSTMORTEM (2)
        CLOSE(UNIT=UW2)
        ISAVE=-1
      ENDIF

C ************************************************************************
C *** AFTER CRYSTAL REORIENTATION UPDATES PX ELASTIC MODULI
      IF(ABS(IVGVAR).LE.1) CALL ELSC (1) 

 2500 CONTINUE      ! END OF DO 2500 ISTEP=1,NSTEPSX OVER DEFORMATION INCREMENTS
 
 3000 CONTINUE      ! END OF DO 3000 OVER PROCESSES

C ************************************************************************
      IF(.NOT.CONVERGED) THEN
	    WRITE( *,'(/,'' CONVERGENCE FAILED IN ONE OR MORE STEPS'',/,
     #    '' CHECK RELSGR,RELS,RELD= IN FILE RUN_LOG.OUT'')') 
	    WRITE(10,'(/,'' CONVERGENCE FAILED IN ONE OR MORE STEPS'',/,
     #    '' CHECK RELSGR,RELS,RELD= IN FILE RUN_LOG.OUT'')') 
      ENDIF

      CALL CPU_TIME(END_TIME)
        TIME_ELAPSED= END_TIME-START_TIME
      WRITE( *,'(/,'' TIME ELAPSED'',F8.2,''secs'')') TIME_ELAPSED
      WRITE(10,'(/,'' TIME ELAPSED'',F8.2,''secs'')') TIME_ELAPSED
C ************************************************************************

      STOP
      END

C ************************************************************************
