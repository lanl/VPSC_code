C
C *******************************************************************
C ******    file VPSC8SUB.FOR with subroutines used by VPSC8   ******
C *******************************************************************
C
C     SUBROUTINE CHECK_VOCE    
C     SUBROUTINE CHG_BASIS    
C     SUBROUTINE CRYSTAL_SYMMETRY   
C     SUBROUTINE CUBCOMP     
C     SUBROUTINE DATA_CRYSTAL
C     SUBROUTINE DATA_GRAIN
C     SUBROUTINE DIFF_PLANES 
C     SUBROUTINE EIGEN_SORT
C     SUBROUTINE EIGEN_VAL
C     SUBROUTINE ELSC        
C
C  ESHELBY TENSOR RELATED SUBROUTINES:
C     SUBROUTINE ESHELBY_TENSOR  
C     SUBROUTINE ESHELBY         
C     SUBROUTINE ESH_GAUSS_LEGENDRE
C     SUBROUTINE ESH_INV3_VOIGT  
C     SUBROUTINE ESH_INV4_VOIGT  
C     SUBROUTINE ESH_MULT_VOIGT

C     SUBROUTINE EULER
C     SUBROUTINE GRAIN_INFO    
C     SUBROUTINE GRAIN_RATE_AND_MODULI
C     SUBROUTINE GRAIN_STRESS         
C     SUBROUTINE INITIAL_STATE_GUESS  
C     SUBROUTINE LANKFORD             
C     SUBROUTINE LOAD_CONDITIONS      

C  MATRIX ALGEBRA SUBROUTINES
C     SUBROUTINE LU_BACKSUBS
C     SUBROUTINE LU_INVERSE(A,N)
C     SUBROUTINE LU_EQSYSTEM
C     SUBROUTINE LU_DECOMP

C     SUBROUTINE N_EFFECTIVE    
C     SUBROUTINE NEIGHBOURS     
C     SUBROUTINE NEWTON_RAPHSON 
C     SUBROUTINE PCYS   
C     SUBROUTINE PCYS_IT        
C     SUBROUTINE POSTMORTEM     
C     SUBROUTINE RODRIGUES      
C     SUBROUTINE SCALE_GAMD0S   
C
C  SECOND-ORDER & FLUCTUATIONS SUBROUTINES:
C     SUBROUTINE SO_FLUCTUATIONS
C     SUBROUTINE SO_LINSOLVER25
C     SUBROUTINE SO_GET_THEFLU
C     SUBROUTINE SO_GET_GAMDOT
C     SUBROUTINE SO_PROCEDURE
C     SUBROUTINE SO_MOD
C     SUBROUTINE SO_SDPX
C     SUBROUTINE SO_GRAIN_STRESS_ALT
C     SUBROUTINE SO_EXTRAPOL
C     SUBROUTINE SO_VOIGT10
C
C     SUBROUTINE STATE_NxN      
C     SUBROUTINE STATE_6x6      
C     SUBROUTINE STAT_GRAIN_SHAPE
C     SUBROUTINE STAT_SHEAR_ACTIVITY      
C     SUBROUTINE STAT_STRESS_STRAIN   
C     SUBROUTINE STAT_TWINNING        
C     SUBROUTINE TEXTURE_ROTATION     
C     SUBROUTINE TWIN_ORIENTATION
C     
C     SUBROUTINE UPDATE_CRSS_DD       
C     SUBROUTINE UPDATE_CRSS_DD_REV   
C     SUBROUTINE UPDATE_CRSS_MTS      
C     SUBROUTINE UPDATE_CRSS_VOCE     
C     SUBROUTINE UPDATE_FIJ  
C     SUBROUTINE UPDATE_GROWTH_RATE     
C     SUBROUTINE UPDATE_GROWTH_ZR     
C     SUBROUTINE UPDATE_GROWTH_ZIRC2  
C     SUBROUTINE UPDATE_ORIENTATION   
C     SUBROUTINE UPDATE_SCHMID
C     SUBROUTINE UPDATE_SHAPE  
C     SUBROUTINE UPDATE_TWINNING_PTR 
C     SUBROUTINE UPDATE_TWINNING_VFT 
C     SUBROUTINE VAR_VEL_GRAD        
C     SUBROUTINE VOIGT    
C     SUBROUTINE VPFC     
C     SUBROUTINE VPSC     
C     SUBROUTINE VPSC_INPUT

C  SUBROUTINES FOR OUTPUT WRITING
C     SUBROUTINE WRITE_SHEAR_ACTIVITY  
C     SUBROUTINE WRITE_STRESS_STRAIN  
C     SUBROUTINE WRITE_TEXTURE        
C     SUBROUTINE WRITE_TWIN_STATS     
C
C  FUNCTIONS:
C     FUNCTION DET
c     FUNCTION RANDOM2
C     FUNCTION VNORM       
C     FUNCTION TNORM5x5    
C     FUNCTION TNORM6x6    
C     FUNCTION VMISMATCH   
C     FUNCTION TMISMATCH   
C *******************************************************************

C *******************************************************************
C     SUBROUTINE CHECK_VOCE      --->      VERSION 12/JUN/2001
C
C     CHECKS WHETHER VOCE PARAMETERS ARE KOSHER:
C        TAU0>0 , TAU1 >= 0 , THET0 >= THET1 >= 0
C        TAU1=0   CORRESPONDS TO LINEAR HARDENING.
C        THETA0=0 FORCES NO-HARDENING.
C     IF VOCE PARAMETERS ARE NON-KOSHER CHECKS FOR ILL-POSED HARDENING.
C *********************************************************************

      SUBROUTINE CHECK_VOCE (IMODE,IPH,TAU0X,TAU1X,THET0X,THET1X)

   70 FORMAT(' *** MODE',I3,'   IN PHASE',I3)
   71 FORMAT('     TAU0.LE.0 --> ILL-POSED HARDENING  !!!')
   72 FORMAT('     THETA1<0 --> NON-KOSHER HARDENING MAY GIVE TAU<0')
   73 FORMAT('     TAU1<0  --> NON-KOSHER HARDENING MAY GIVE TAU<0')
   75 FORMAT('     THETA0=0 --> WILL RESET TAU1=THETA1=0')
   76 FORMAT('     |THETA1|.GE.|THETA0| --> NON-KOSHER HARDENING')
   77 FORMAT('     THETA0<0 --> HARDENING DEPENDS ON |THETA0| ANYWAY')

      TINY=1.E-4*TAU0X
      IF(TAU0X.LE.0.0) THEN
        WRITE(*,70) IMODE,IPH
        WRITE(*,71)
        STOP
      ENDIF
      IF(ABS(THET1X).LE.TINY) THEN
        THET1X=0.
      ELSE IF(THET1X.LT.0.0) THEN
        WRITE(*,70) IMODE,IPH
        WRITE(*,72)
        print *, 'enter c to continue'
        read  *
      ENDIF
      IF(TAU1X.LT.0) THEN
        WRITE(*,70) IMODE,IPH
        WRITE(*,73)
        print *, 'enter c to continue'
        read  *
      ENDIF
      IF(ABS(THET0X).LE.TINY) THEN
        IF(ABS(TAU1X).LE.TINY) THEN
          TAU1X =0.
          THET0X=THET1X
        ENDIF
        IF(ABS(TAU1X).GT.TINY) THEN
          WRITE(*,70) IMODE,IPH
          WRITE(*,75)
          print *, 'enter c to continue'
          read  *
          TAU1X =0.
          THET0X=0.
          THET1X=0.
        ENDIF
      ENDIF
      IF(ABS(TAU1X).LE.TINY) THEN
        TAU1X =0.0
        THET0X=THET1X
      ENDIF
      IF(THET0X.LT.0.0) THEN
        WRITE(*,70) IMODE,IPH
        WRITE(*,77)
        print *, 'enter c to continue'
        read  *
        THET0X=ABS(THET0X)
      ENDIF
      IF(TAU1X.NE.0.) THEN
        IF(ABS(THET1X).GE.ABS(THET0X)) THEN
          WRITE(*,70) IMODE,IPH
          WRITE(*,76)
          print *, 'enter c to continue'
          read  *
        ENDIF
      ENDIF

      RETURN
      END
C
C ************************************************************************
C     SUBROUTINE CHG_BASIS    --->   VERSION 02/APR/2023
C
C     KDIM=5 or 6, FOR DEVIATORIC or DEV+HYDROST TENSORS, RESPECTIVELY.
C     IOPT=0: INITIALIZES A b-BASIS OF 6 SECOND ORDER TENSORS B(3,3,N).
C     IOPT=1: CALCULATES SECOND ORDER TENSOR C2(3,3) AS AN EXPANSION IN TERMS
C             OF VECTOR COMPONENTS CE2(KDIM) AND THE b-BASIS TENSORS.
C     IOPT=2: CALCULATES COMPONENTS OF C2 AS A VECTOR CE2(KDIM).
C     IOPT=3: CALCULATES FOURTH ORDER TENSOR C4(3,3,3,3) AS AN EXPANSION IN 
C             TERMS OF MATRIX COMPONENTS CE4(KDIM,KDIM) AND THE b-BASIS TENSORS.
C     IOPT=4: CALCULATES MATRIX COMPONENTS CE4(KDIM,KDIM) OF TENSOR 'C4'.
C **************************************************************************

      SUBROUTINE CHG_BASIS(CE2,C2,CE4,C4,IOPT,KDIM)

      USE CHANGE_BASIS

      DIMENSION CE2(KDIM),C2(3,3),CE4(KDIM,KDIM),C4(3,3,3,3)

C *** INITIALIZES b-BASIS TENSORS B(N)
      IF(IOPT.EQ.0) THEN

        B(1:3,1:3,1:6)=0.0

        B(1,1,2)=-1./SQRT(6.)
        B(2,2,2)=-1./SQRT(6.)
        B(3,3,2)= 2./SQRT(6.)

        B(1,1,1)=-1./SQRT(2.)
        B(2,2,1)= 1./SQRT(2.)

        B(2,3,3)= 1./SQRT(2.)
        B(3,2,3)= 1./SQRT(2.)

        B(1,3,4)= 1./SQRT(2.)
        B(3,1,4)= 1./SQRT(2.)

        B(1,2,5)= 1./SQRT(2.)
        B(2,1,5)= 1./SQRT(2.)

        B(1,1,6)= 1./SQRT(3.)
        B(2,2,6)= 1./SQRT(3.)
        B(3,3,6)= 1./SQRT(3.)

      ENDIF

C *** CALCULATES CARTESIAN SECOND ORDER TENSOR FROM b-COMPONENTS VECTOR.
      IF(IOPT.EQ.1) THEN
        DO I=1,3
        DO J=1,3
          C2(I,J)=0.0
          DO N=1,KDIM
            C2(I,J)=C2(I,J)+CE2(N)*B(I,J,N)
          ENDDO
        ENDDO
        ENDDO
      ENDIF

C *** CALCULATES KDIMx1 b-COMPONENTS VECTOR FROM SECOND ORDER TENSOR.
      IF(IOPT.EQ.2) THEN
        DO N=1,KDIM
          CE2(N)=0.0
          DO I=1,3
          DO J=1,3
            CE2(N)=CE2(N)+C2(I,J)*B(I,J,N)
          ENDDO
          ENDDO
        ENDDO
      ENDIF

C *** CALCULATES FOURTH ORDER TENSOR FROM KDIMxKDIM b-COMPONENTS MATRIX.
      IF(IOPT.EQ.3) THEN
        DO I=1,3
        DO J=1,3
        DO K=1,3
        DO L=1,3
          C4(I,J,K,L)=0.0
          DO N=1,KDIM
          DO M=1,KDIM
            C4(I,J,K,L)=C4(I,J,K,L)+CE4(N,M)*B(I,J,N)*B(K,L,M)
          ENDDO
          ENDDO
        ENDDO
        ENDDO
        ENDDO
        ENDDO
      ENDIF

C *** CALCULATES KDIMxKDIM b-COMPONENTS MATRIX FROM FOURTH ORDER TENSOR.
      IF(IOPT.EQ.4) THEN
        DO N=1,KDIM
        DO M=1,KDIM
          CE4(N,M)=0.0
          DO I=1,3
          DO J=1,3
          DO K=1,3
          DO L=1,3
            CE4(N,M)=CE4(N,M)+C4(I,J,K,L)*B(I,J,N)*B(K,L,M)
          ENDDO
          ENDDO
          ENDDO
          ENDDO
        ENDDO
        ENDDO
      ENDIF

      RETURN
      END

c ***********************************************************************
C     SUBROUTINE CRYSTAL_SYMMETRY   --->   version 17/DEC/2020
c
c *** If IOPTION=1:
c     Reads crystal symmetry 'icrysym' and unit cell parameters.
c     Generates vectors 'cvec(i,n)' of the unit cell.
c     Generates symmetry operators 'h(i,j,nsymop)' for all crystal symmetries.
c *** If IOPTION=2:
c     Reads Miller indices of systems in 3 or 4-index notation 'isn(i)'
c     & 'isb(i)'. Calculates normal & burgers vectors 'sn(i)' & 'sb(i)'
c *** If IOPTION=3:
c     Generates 'nequiv' crystallographically equivalent orientations sneq(i,n)
c     of normal vector sn(i) by applying all the symmetry operations to it.
c     Discards repeated orientations and defines 'nequiv'.
c
c *** Simmetry parameter ICRYSYM:
c        1: CUBIC
c        2: HEXAGONAL
c        3: TRIGONAL
c        4: TETRAGONAL
c        5: ORTHORHOMBIC
c        6: MONOCLINIC
c        7: TRICLINIC
c ***********************************************************************

      SUBROUTINE CRYSTAL_SYMMETRY (ioption,ur1,icrysym,
     #                             isn,sn,sneq,isb,sb,nequiv)

      USE CRYSTAL_SYM

      dimension isn(4),sn(3),sneq(3,24),isb(4),sb(3)
      dimension hx(3,3,6),itag(24)
      dimension cdim(3),cang(3)

      integer ur1
      character crysym*5

      pi=4.*atan(1.)

c ****************************************************************************
      if(ioption.eq.1) then
c ****************************************************************************

        read(ur1,*)
        read(ur1,'(a)') crysym
        icrysym=0
        if(crysym.eq.'cubic' .or. crysym.eq.'CUBIC') icrysym=1
        if(crysym.eq.'hexag' .or. crysym.eq.'HEXAG') icrysym=2
        if(crysym.eq.'trigo' .or. crysym.eq.'TRIGO') icrysym=3
        if(crysym.eq.'tetra' .or. crysym.eq.'TETRA') icrysym=4
        if(crysym.eq.'ortho' .or. crysym.eq.'ORTHO') icrysym=5
        if(crysym.eq.'monoc' .or. crysym.eq.'MONOC') icrysym=6
        if(crysym.eq.'tricl' .or. crysym.eq.'TRICL') icrysym=7
        if(icrysym.eq.0) then
          write(*,*) ' *** CANNOT RECOGNIZE THE CRYSTAL SYMMETRY'
          stop
        endif

        READ(UR1,*) (CDIM(i),i=1,3),(CANG(i),i=1,3)
        DO I=1,3
          CANG(I)=CANG(I)*PI/180.
        ENDDO

c *** assumes 'c' coincident with 'z' and 'a' in the plane 'xz'
        CVEC(1,1)=SIN(CANG(2))
        CVEC(2,1)=0.
        CVEC(3,1)=COS(CANG(2))
        CVEC(1,2)=(COS(CANG(3))-COS(CANG(1))*COS(CANG(2)))/SIN(CANG(2))
        CVEC(3,2)=COS(CANG(1))
        CVEC(2,2)=SQRT(1.-CVEC(1,2)**2-CVEC(3,2)**2)
        CVEC(1,3)=0.
        CVEC(2,3)=0.
        CVEC(3,3)=1.

        DO J=1,3
        DO I=1,3
          CVEC(I,J)=CDIM(J)*CVEC(I,J)
        ENDDO
        ENDDO

        HX(:,:,:)=0.d0
        H (:,:,:)=0.d0

c *** identity operation ---> triclinic & all symmetries
      do i=1,3
        h(i,i,1)=1.d0
      enddo
      nsymop=1

c *** 180 deg rotation around (001) ---> orthorhombic, monoclinic
      if(icrysym.eq.5 .or. icrysym.eq.6) then
        h(1,1,2)= cos(pi)
        h(2,2,2)= cos(pi)
        h(3,3,2)= 1.d0
        h(1,2,2)=-sin(pi)
        h(2,1,2)= sin(pi)
        nsymop=2
      endif

c *** x-mirror & y-mirror ---> orthorhombic
      if(icrysym.eq.5) then
        h(1,1,3)=-1.d0
        h(2,2,3)= 1.d0
        h(3,3,3)=-1.d0      ! invert axis 3 to preserve right-handed system

        h(1,1,4)= 1.d0
        h(2,2,4)=-1.d0
        h(3,3,4)=-1.d0      ! invert axis 3 to preserve right-handed system
        nsymop=4
      endif

c *** cubic symmetry
      if(icrysym.eq.1) then

c *** rotations of 120 and 240 deg around <111>
        hx(1,3,1)= 1.d0
        hx(2,1,1)= 1.d0
        hx(3,2,1)= 1.d0

        hx(1,2,2)= 1.d0
        hx(2,3,2)= 1.d0
        hx(3,1,2)= 1.d0

        do m=1,2
          do n=1,nsymop
            mn=m*nsymop+n
            do i=1,3
            do j=1,3
            do k=1,3
              h(i,j,mn)=h(i,j,mn)+hx(i,k,m)*h(k,j,n)
            enddo
            enddo
            enddo
          enddo
        enddo
        nsymop=mn

c *** mirror across the plane (110)
        hx(1,2,3)= 1.d0
        hx(2,1,3)= 1.d0
        hx(3,3,3)=-1.d0      ! invert axis 3 to preserve right-handed system

        do n=1,nsymop
          mn=nsymop+n
            do i=1,3
            do j=1,3
            do k=1,3
              h(i,j,mn)=h(i,j,mn)+hx(i,k,3)*h(k,j,n)
            enddo
            enddo
            enddo
        enddo
        nsymop=mn

c *** rotations of 90, 180, 270 around x3

        do m=1,3
          ang=pi/2.*float(m)
          hx(1,1,m)= cos(ang)
          hx(2,2,m)= cos(ang)
          hx(3,3,m)= 1.0
          hx(1,2,m)=-sin(ang)
          hx(2,1,m)= sin(ang)
          hx(1,3,m)= 0.0
          hx(3,1,m)= 0.0
          hx(2,3,m)= 0.0
          hx(3,2,m)= 0.0
        enddo

        do m=1,3
          do n=1,nsymop
            mn=m*nsymop+n
              do i=1,3
              do j=1,3
              do k=1,3
                h(i,j,mn)=h(i,j,mn)+hx(i,k,m)*h(k,j,n)
              enddo
              enddo
              enddo
          enddo
        enddo
        nsymop=mn

      endif                    !end of condition for icrysym=1

c *** hexagonal, trigonal and tetragonal symmetry

      if(icrysym.ge.2 .and. icrysym.le.4) then
        if(icrysym.eq.2) nrot=6
        if(icrysym.eq.3) nrot=3        
        if(icrysym.eq.4) nrot=4

c *** mirror plane at 30 deg or 60 deg or 45 deg with respect to x1
        ang=pi/float(nrot)
        h(1,1,2)= cos(ang)**2-sin(ang)**2
        h(2,2,2)=-h(1,1,2)
        h(3,3,2)=-1.d0      ! invert axis 3 to preserve right-handed system
        h(1,2,2)= 2.*cos(ang)*sin(ang)
        h(2,1,2)= h(1,2,2)
        nsymop=2

c *** rotations of 2*pi/6 around axis <001> for hexagonals.
c *** rotations of 2*pi/3 around axis <001> for trigonals.
c *** rotations of 2*pi/4 around axis <001> for tetragonals.
        do nr=1,nrot-1
          ang=nr*2.*pi/nrot
          hx(1,1,nr)= cos(ang)
          hx(2,2,nr)= cos(ang)
          hx(3,3,nr)= 1.d0
          hx(1,2,nr)=-sin(ang)
          hx(2,1,nr)= sin(ang)
        enddo

        do m=1,nrot-1
          do n=1,nsymop
            mn=m*nsymop+n
            do i=1,3
            do j=1,3
            do k=1,3
              h(i,j,mn)=h(i,j,mn)+hx(i,k,m)*h(k,j,n)
            enddo
            enddo
            enddo
          enddo
        enddo
        nsymop=mn

      endif               !end of condition for icrysym= 2,3,4

c      write(10,*)
c      write(10,'(''  # of symmetry operations='',i4)') nsymop
c      do n=1,nsymop
c        write(10,'(''  symmetry matrix #'',i5)') n
c        write(10,'(3f7.3)') ((h(i,j,n),j=1,3),i=1,3)
c     enddo

      endif               !end of condition for ioption=1

c **************************************************************************
c   Converts Miller-Bravais indices of plane normal and slip direction
c   into normalized vectors sn(i) and sb(i), respectively.
c   Indices for cubic (1), tetragonal (4), orthorhombic (5), monoclinic (6)
c   & triclinic (7) systems are in 3-index notation.
c   For hexagonal (2) & trigonal (3) systems uses 4-index notation.
c **************************************************************************
      if (ioption.eq.2) then
c **************************************************************************

        if(icrysym.eq.2 .or. icrysym.eq.3) then
          isn(3)=isn(4)
          isb(1)=isb(1)-isb(3)
          isb(2)=isb(2)-isb(3)
          isb(3)=isb(4)
        endif

c *** assumes 'c' coincident with 'z' and 'a' in the plane 'xz'
        sn(3)= isn(3)/cvec(3,3)
        sn(1)=(isn(1)-cvec(3,1)*sn(3))/cvec(1,1)
        sn(2)=(isn(2)-cvec(1,2)*sn(1)-cvec(3,2)*sn(3))/cvec(2,2)

        snnor=sqrt(sn(1)**2+sn(2)**2+sn(3)**2)
        do j=1,3
          sn(j)=sn(j)/snnor
          if(abs(sn(j)).lt.1.e-03) sn(j)=0.
        enddo

c *** this block specific for EPSC & VPSC

        do i=1,3
          sb(i)=isb(1)*cvec(i,1)+isb(2)*cvec(i,2)+isb(3)*cvec(i,3)
        enddo
        sbnor=sqrt(sb(1)**2+sb(2)**2+sb(3)**2)
        do j=1,3
          sb(j)=sb(j)/sbnor
          if(abs(sb(j)).lt.1.e-03) sb(j)=0.
        enddo

      endif      ! end of if(ioption.eq.2)

c **************************************************************************
      IF(IOPTION.EQ.3) THEN                                                    
c **************************************************************************

        NIND=3                                                                 
        IF (ICRYSYM.EQ.2 .OR. ICRYSYM.EQ.3) NIND=4 
                                                                               
c *** This block required when CRYSTAL_SYMMETRY called from DIFF_PLANES.
c *** Reads diffraction directions defined by polar angles referred to sample axes.
c *** Calculate cartesian components in sample axes.
c *** Has to be commented out when called by POLE8.
                                      
        READ(UR1,*) (ISN(I),I=1,NIND),CHI,ETA                                
        ETA=ETA*PI/180.0                                                     
        CHI=CHI*PI/180.0                                                     
        SB(1)=COS(ETA)*SIN(CHI)                                            
        SB(2)=SIN(ETA)*SIN(CHI)                                            
        SB(3)=         COS(CHI)                                              

c *** This block required when called by POLE8 and by DIFF_PLANES.
c *** Generates all symmetry related vectors sneq(i,n) with z>0.
c *** Eliminates redundant poles: coincidents and opposites.
                            
        IF(NIND.EQ.4) ISN(3)=ISN(4)                                            
        SN(1)= ISN(1)/CVEC(1,1)                                                
        SN(2)=(ISN(2)-CVEC(1,2)*SN(1))/CVEC(2,2)                               
        SN(3)=(ISN(3)-CVEC(1,3)*SN(1)-CVEC(2,3)*SN(2))/CVEC(3,3)               
        SNNOR=SQRT(SN(1)**2+SN(2)**2+SN(3)**2)                             
        DO J=1,3                                                               
          SN(J)=SN(J)/SNNOR                                                    
          IF(ABS(SN(J)).LT.1.D-03) SN(J)=0.D0                              
        ENDDO                                                                  

C *** calculates all crystallographically equivalent plane normals                                                                               
        DO N=1,NSYMOP                                                          
          ITAG(N)=0                                                            
          DO I=1,3                                                             
          SNEQ(I,N)=0.D0                                                       
            DO J=1,3                                                           
              SNEQ(I,N)=SNEQ(I,N)+H(I,J,N)*SN(J)                               
            ENDDO                                                              
          ENDDO                                                                
        ENDDO                                                                  

C *** keeps only non-repeated plane normals.                                                                               
        IF(ICRYSYM.NE.7) THEN      ! NSYMOP=1 FOR TRIGONAL                     
          DO M=1,NSYMOP-1                                                      
            IF(ITAG(M).EQ.0) THEN                                              
              DO N=M+1,NSYMOP                                                  
                SNDIF=ABS(SNEQ(1,M)-SNEQ(1,N))+ABS(SNEQ(2,M)-SNEQ(2,N))       
     #              +ABS(SNEQ(3,M)-SNEQ(3,N))                                 
                IF(SNDIF .LE. 0.0001) ITAG(N)=1                                
                SNDIF=ABS(SNEQ(1,M)+SNEQ(1,N))+       
     #               ABS(SNEQ(2,M)+SNEQ(2,N))+ABS(SNEQ(3,M)+SNEQ(3,N))                               
                IF(SNDIF .LE. 0.0001) ITAG(N)=1                                
              ENDDO                                                            
            ENDIF                                                              
          ENDDO                                                                
        ENDIF                                                                  

c *** takes the opposite if n(3)<0.                                                                                
        NEQUIV=0                                                               
        DO N=1,NSYMOP                                                          
          IF(ITAG(N).EQ.0) THEN                                                
            NEQUIV=NEQUIV+1                                                    
            ISIGN=1                                                            
            IF(SNEQ(3,N).LT.0.) ISIGN=-1                                       
            SNEQ(1,NEQUIV)=ISIGN*SNEQ(1,N)                                     
            SNEQ(2,NEQUIV)=ISIGN*SNEQ(2,N)                                     
            SNEQ(3,NEQUIV)=ISIGN*SNEQ(3,N)                                     
          ENDIF                                                                
        ENDDO                                                                  
                                                                               
      ENDIF            !END OF IF(IOPTION=3)                                   
C **************************************************************************   
 
      RETURN
      END

C ***************************************************************************
C     SUBROUTINE CUBCOMP      --->      VERSION 19/SET/00
C
C     ASSIGNS ORIENTATIONS TO ONE OF THE IDEAL FCC ROLLING COMPONENTS:
C     (CUBE/rotated CUBE/GOSS/BRASS/COPPER/S/other).
C     WORKS FOR MULTIPHASE BUT ALL PHASES HAVE TO BE FCC FOR THE RESULTS
C     TO MAKE SENSE.
C ***************************************************************************

      SUBROUTINE CUBCOMP (ISTEP,IOPTION)

      USE VPSC8DIM
      USE CUB_COMP

      DIMENSION dmin(NGRMX),pidmod(0:NIDMODMX)
      DIMENSION aa(3,3),amis(3,3)

C **********************************************************************
C     READS ORIENTATION MATRICES OF IDEAL COMPONENTS ALREADY TRANSPOSED
C **********************************************************************

      IF(IOPTION.EQ.0) THEN

        read(UR5,*) nidmod
        idlabel(0)='OTH'
        nor=0
        do im=1,nidmod
          read(UR5,*) normod(im)
          read(UR5,'(a3)') idlabel(im)
          do j=1,normod(im)
            nor=nor+1
            do j1=1,3
              read(UR5,*) (aidort(i1,j1,nor),i1=1,3)
            enddo
            itype(nor)=im
          enddo
        enddo

        RETURN
      ENDIF

C **********************************************************************
C     IDENTIFIES ORIENTATIONS IN EACH PHASE ALIGNED WITH ONE OF THE IDEAL
C     COMPONENTS AND ACCUMULATES
C **********************************************************************

      DO IPH=IPHBOT,IPHTOP

      totwgt=0.
      dminav=0.
      do im=0,nidmod
        widmod(im)=0.
      enddo

      DO KKK=ngr(iph-1)+1,ngr(iph)
      totwgt=totwgt+wgt(kkk)

       do i=1,3
       do j=1,3
         aa(i,j)=ag(i,j,kkk)
       enddo
       enddo

      if(aa(3,3).lt.0.) then
        do i=2,3
        do j=1,3
          aa(i,j)=-aa(i,j)
        enddo
        enddo
      endif

      dmin(kkk)=500.
      do ior=1,nor

        do i=1,3
        do j=1,3
          amis(i,j)=0.
          do k=1,3
            amis(i,j)=amis(i,j)+aidort(i,k,ior)*aa(k,j)
          enddo
        enddo
        enddo

        trace=(amis(1,1)+amis(2,2)+amis(3,3))
        arg=(trace-1)/2.
        if(arg.gt.1) arg=1.
        if(arg.lt.-1) arg=-1.
        angmis=acos(arg)

        if(abs(angmis).lt.dmin(kkk)) then
          dmin(kkk)=abs(angmis)
          igrtype(kkk)=itype(ior)
        endif

      enddo

      if(dmin(kkk).gt.(15.*pi/180.)) igrtype(kkk)=0
      dminav=dminav+dmin(kkk)*wgt(kkk)
      widmod(igrtype(kkk))=widmod(igrtype(kkk))+wgt(kkk)

      enddo    ! end of loop over grains in the phase

      dminav=(dminav/totwgt)*(180./pi)

      do im=0,nidmod
        pidmod(im)=widmod(im)*100./totwgt
      enddo

      IUNIT=60+IPH
      IF(ISTEP.EQ.0) WRITE(IUNIT,'(13a7)')
     #               '    EPS',(idlabel(im),im=0,nidmod),' AVMISO'
      WRITE(IUNIT,'(f7.4,12(2x,f5.1))')
     #               EPSACU,(PIDMOD(IM),IM=0,NIDMOD),DMINAV

      enddo      ! end of loop over all the phases

      end

C *****************************************************************************
C     SUBROUTINE DATA_CRYSTAL       --->     VERSION 02/JUN/2018

C *** READS ELASTIC & THERMAL MODULI OF CRYSTAL. READS CRYSTALLOGRAPHIC SYSTEMS
C     TO BE USED. CALCULATES NORMAL AND SHEAR VECTORS OF SLIP & TWIN SYSTEMS.
C     CALCULATES SCHMID TENSORS. (ALL TENSORS REFERRED TO CRYSTAL AXES)
C *****************************************************************************

      SUBROUTINE DATA_CRYSTAL (IPH)

      USE VPSC8DIM

      DIMENSION MODE(50)
      DIMENSION ISN(4),SN(3),SNEQ(3,24),ISB(4),SB(3)

C -----------------------------------------------------------------
C     WRITES CRYSTAL DATA FILE INTO 'RUN_LOG.OUT' FILE
      WRITE(10,*)
      WRITE(10,'('' **** CRYSTAL DATA FILE ****'')')
      DO IDUM=1,200
        READ(UNIT=UR1,END=99,FMT='(A)') PROSA
        WRITE(10,'(A)') PROSA
      ENDDO
   99 REWIND UR1
      WRITE(10,'(''**** END OF CRYSTAL DATA FILE ****'')')
      WRITE(10,*)
C -----------------------------------------------------------------

C *** READS CRYSTAL SYMMETRY & UNIT CELL PARAMETERS. CALCULATES CELL VECTORS.
C *** GENERATES ALL SYMMETRY OPERATIONS ASSOCIATED WITH 'CRYSYM'.

      CALL CRYSTAL_SYMMETRY (1,UR1,ICRYSYM,ISN,SN,SNEQ,ISB,SB,NPOLES)
      ICRYSYMPH(IPH)=ICRYSYM
      NIND=3
      IF(ICRYSYM.EQ.2 .OR. ICRYSYM.EQ.3) NIND=4

C *** READS SINGLE CRYSTAL ELASTIC STIFFNESS (in Voigt notation)
      READ(UR1,'(A)') PROSA
      READ(UR1,*)     ((CELCCv(I,J,IPH),J=1,6),I=1,6)
C *** READS SINGLE CRYSTAL THERMAL EXPANSION COEFFICIENTS (in Voigt notation)
      READ(UR1,'(A)') PROSA
      READ(UR1,*)     (ATHCCv(I,IPH),I=1,6)

C *** CALCULATES & STORES SX MODULI (b-BASIS & VOIGT NOTATION)
      CALL VOIGT(AUX6,AUX33,CELCCv(:,:,IPH),AUX3333,3)
      CALL CHG_BASIS(AUX6,AUX33,CELCC(:,:,IPH),AUX3333,4,6)   ! STIFFNESS

      SELCC(:,:,IPH)=CELCC(:,:,IPH)
        CALL LU_INVERSE(SELCC(:,:,IPH),6)
      CALL CHG_BASIS(AUX6,AUX33,SELCC(:,:,IPH),AUX3333,3,6)
      CALL VOIGT(AUX6,AUX33,SELCCv(:,:,IPH)   ,AUX3333,4)     ! COMPLIANCE

      CALL VOIGT(ATHCCv(:,IPH)   ,AUX33,AUX66,AUX3333,1)
      CALL CHG_BASIS(ATHCC(:,IPH),AUX33,AUX66,AUX3333,2,6)    ! THERMAL
	  
C --> VERIFICATION OF CORRECT STIFFNESS INVERSION
      RESIDUAL=0.
      DO I=1,6
      DO J=1,6
        AUX66(I,J)=0.
        DO K=1,6
          AUX66(I,J)=AUX66(I,J)+CELCC(I,K,IPH)*SELCC(K,J,IPH)
        ENDDO
        RESIDUAL=RESIDUAL+ABS(AUX66(I,J)-(I/J)*(J/I))
      ENDDO
      ENDDO
      WRITE(10,*)
      WRITE(10,'('' CHECKING THAT CELCC*SELCC-ID6=0 '',
     #              E15.7)')  RESIDUAL 

C *** UPPER BOUND ISOTROPIC ELASTIC CONSTANTS FOR A RANDOM PX
C *** UPPER BOUND BULK MODULUS (Delta_P/Delta_V/V) FOR A TEXTURED POLYCRYSTAL

      ALF=(CELCCv(1,1,IPH)+CELCCv(2,2,IPH)+CELCCv(3,3,IPH))/3.
      BET=(CELCCv(1,2,IPH)+CELCCv(1,3,IPH)+CELCCv(2,3,IPH))/3.
      GAM=(CELCCv(4,4,IPH)+CELCCv(5,5,IPH)+CELCCv(6,6,IPH))/3.
      BULK=(3.*ALF+6.*BET)/9.
      C11ISO=(3.*ALF+2.*BET+4.*GAM)/5.
      C12ISO=(   ALF+4.*BET-2.*GAM)/5.
      C44ISO=(C11ISO-C12ISO)/2.
      XNU=(C11ISO-2.*C44ISO)/(2.*C11ISO-2.*C44ISO)

      WRITE(10,'('' *********** PHASE'',I4)') IPH
      WRITE(10,'('' RANDOM PX BULK & POISSON MODULI'',2F12.3)')
     #              BULK,XNU
      WRITE(10,'('' RANDOM PX ELASTIC CTES C11, C12, C44'',3F12.3)')
     #              C11ISO,C12ISO,C44ISO

C ********************************************************************
C *** READS INFORMATION ABOUT SLIP AND TWINNING SYSTEMS.  
C *** THIS INFORMATION IS NOT NECESSARY or USED FOR ELASTIC SIMULATION

      IF(INTERACTION.EQ.-1) RETURN
C ********************************************************************

      READ(UR1,'(A)') PROSA
      READ(UR1,*)     NMODESX
      READ(UR1,*)     NMODES(IPH)
      READ(UR1,*)     (MODE(I),I=1,NMODES(IPH))

      IF(NMODES(IPH).GT.NMODMX) THEN
        WRITE(*,'('' NMODES IN PHASE'',I3,'' IS'',I3)') IPH,NMODES(IPH)
        WRITE(*,'('' CHANGE PARAMETER NMODMX IN VPSC.DIM'')')
        STOP
      ENDIF

      ICS=1       ! ICS=1 --> centro_symmetric SCYS
      NSLMOD(iph)=0
      NTWMOD(iph)=0
      NSYST(iph) =0
      NSLSYS(iph)=0
      NTWSYS(iph)=0
      MCOUNT=1    ! counter for the number of modes

C *** READS DEFORMATION MODES AND ASSOCIATED PARAMETERS FROM FILECRYS

      DO 100 MLOOP=1,NMODESX

        READ(UR1,'(a)') PROSA
        READ(UR1,*)     MODEX,NSMX,ISENSEX,ITWTYPEX
        IF(MODEX.NE.MLOOP) THEN
          WRITE(*,*) ' WARNING !!!'
          WRITE(*,*) ' MODE NUMBERS MUST BE SEQUENTIAL IN CRYSTAL FILE'
          STOP
        ENDIF
        IF(ITWTYPEX.NE.0) THEN
          IF(ISENSEX.EQ.1) THEN
            WRITE(*,*) ' WARNING: ISENSEX HAS TO BE 0 FOR TWIN SYSTEMS'
            STOP
          ENDIF
          READ(UR1,*) TWSHX
        ENDIF

c       write(10,'(''mloop,mcount,mode(mcount)'',3i5)')
c    #                     mloop,mcount,mode(mcount)

C *** SKIPS MODE IF IT IS NOT IN THE LIST
        IF(MLOOP.NE.MODE(MCOUNT)) THEN
          DO IS=1,NSMX
            READ(UR1,*)
          ENDDO
          GO TO 100
        ENDIF

        IF(ISENSEX.EQ.0) ICS=0      ! ICS=0 --> non-centro-sym SCYS

        NSM(MCOUNT,IPH)=NSMX

C ***********************************************************************
C *** VERIFIES THAT SLIP SYSTEMS PRECEDE TWINNING SYSTEMS IN THE SEQUENCE

        IF(ITWTYPEX.EQ.0. .AND. NTWMOD(IPH).NE.0) THEN
          WRITE(*,'(''  SLIP MODE'',I4,''  IN PHASE'',I4)') MLOOP,IPH
          WRITE(*,*) '  SLIP MODES MUST PRECEDE TWIN MODES'
          WRITE(*,*) '  -->   REORDER CRYSTAL FILE'
          STOP
        ENDIF

        NSYST(IPH)=NSYST(IPH)+NSMX
        IF(ITWTYPEX.EQ.0) THEN
          NSLMOD(IPH)=NSLMOD(IPH)+1
          NSLSYS(IPH)=NSLSYS(IPH)+NSMX
        ELSE IF(ITWTYPEX.NE.0) THEN
          NTWMOD(IPH)=NTWMOD(IPH)+1
          NTWSYS(IPH)=NTWSYS(IPH)+NSMX
          TWSH(NTWMOD(IPH),IPH) =TWSHX
        ENDIF

c     write(10,'(''nsyst(iph),nslsys(iph),ntwsys(iph)'',3i5)')
c    #             nsyst(iph),nslsys(iph),ntwsys(iph)

        IF(NSYST(IPH).GT.NSYSMX) THEN
          WRITE(*,'('' NSYST IN PHASE'',I3,'' IS'',I3)') IPH,NSYST(IPH)
          WRITE(*,'('' --> CHANGE PARAMETER NSYSMX IN VPSC.DIM'')')
          STOP
        ENDIF

        IF(NTWMOD(IPH).GT.NTWMMX) THEN
          WRITE(*,'('' NTWMOD IN PHASE'',I3,'' IS'',I3)')
     #                 IPH,NTWMOD(IPH)
          WRITE(*,'('' --> CHANGE PARAMETER NTWMMX IN VPSC.DIM'')')
          STOP
        ENDIF

        IF(NTWSYS(IPH).GT.NTWSMX) THEN
          WRITE(*,'('' NTWSYS IN PHASE'',I3,'' IS'',I3)') IPH,
     #                 NTWSYS(IPH)
          WRITE(*,'('' --> CHANGE PARAMETER NTWSMX IN VPSC.DIM'')')
          STOP
        ENDIF
C ********************************************************************

      ISBOT=NSYST(IPH)-NSMX+1
      ISTOP=NSYST(IPH)
      DO 200 ISYS=ISBOT,ISTOP

        ISENSE (ISYS,IPH)=ISENSEX
        ITWTYPE(ISYS,IPH)=ITWTYPEX

C *** READS MILLER INDICES AND CALCULATES CARTESIAN COMPONENTS OF NORMAL
C *** (SN) AND SHEAR (SB) VECTORS OF SLIP OR TWIN SYSTEMS.
C *** CALCULATES SCHMID TENSORS IN CRYSTAL AXES 'SCHCA' FOR EACH SYSTEM.

        READ(UR1,*) (ISN(I),I=1,NIND),(ISB(I),I=1,NIND)

        CALL CRYSTAL_SYMMETRY (2,UR1,ICRYSYM,ISN,SN,SNEQ,ISB,SB,NPOLES)

        PROD=SN(1)*SB(1)+SN(2)*SB(2)+SN(3)*SB(3)
        IF(ABS(PROD) .GE. 1.E-4) THEN
          WRITE(*,'('' SYSTEM IS NOT ORTHOGONAL !!'')')
          WRITE(*,'('' ISN='',4I7)') (ISN(J),J=1,NIND)
          WRITE(*,'('' ISB='',4I7)') (ISB(J),J=1,NIND)
          WRITE(*,'(''   N='',3F7.3)') (SN(J),J=1,3)
          WRITE(*,'(''   B='',3F7.3)') (SB(J),J=1,3)
          STOP
        ENDIF

        DO I=1,3
          DNCA(I,ISYS,IPH)=SN(I)
          DBCA(I,ISYS,IPH)=SB(I)
          DO J=1,3
            AUX33(I,J)=0.5*(SB(I)*SN(J)+SB(J)*SN(I))
          ENDDO
        ENDDO
        CALL CHG_BASIS(AUX5,AUX33,AUX55,AUX3333,2,5)
        DO I=1,5
          SCHCA(I,ISYS,IPH)=AUX5(I)
        ENDDO

  200 CONTINUE    ! END OF LOOP OVER DEFORMATION MODES

C *** USE INFO FROM LAST SYSTEM TO CALCULATE SHEAR MODULUS OF THE SLIP MODE
      CALL VOIGT (AUX6,AUX33,CELCCv(:,:,IPH),AUX3333,3)
      MU_MODE(MCOUNT,IPH)=0.
      DO I=1,3
      DO J=1,3
      DO K=1,3
      DO L=1,3
        MU_MODE(MCOUNT,IPH)=MU_MODE(MCOUNT,IPH)+
     #                      AUX33(I,J)*AUX3333(I,J,K,L)*AUX33(K,L)
      ENDDO 
      ENDDO
      ENDDO
      ENDDO

      MCOUNT=MCOUNT+1

  100 CONTINUE    ! END OF LOOP OVER ALL MODES IN PHASE 'IPH'

C *** CHECKS WHETHER THE SINGLE CRYSTAL YIELD SURFACE IS OPEN

      DO ICOMP=1,5
        ICLOSEPOS=0
        ICLOSENEG=0
        DO NS=1,NSLSYS(IPH)
          IF(ABS(SCHCA(ICOMP,NS,IPH)).GT.1.E-3) THEN
            ICLOSEPOS=1
            ICLOSENEG=1
          ENDIF
        ENDDO
        IF(NTWSYS(IPH).NE.0) THEN
          DO NS=NSLSYS(IPH)+1,NSYST(IPH)
            IF(SCHCA(ICOMP,NS,IPH).GT. 1.E-3) ICLOSEPOS=1
            IF(SCHCA(ICOMP,NS,IPH).LT.-1.E-3) ICLOSENEG=1
          ENDDO
        ENDIF
        IF(ICLOSEPOS.NE.1 .OR. ICLOSENEG.NE.1) THEN
          WRITE(*,'('' WARNING ! THE SCYS IS OPEN FOR PHASE'',I5,
     #              '' ALONG DIRECTION'',I5)') IPH,ICOMP
          print *, 'enter c to continue'
          read  *
        ENDIF
      ENDDO

      WRITE(10,*)
      WRITE(10,'('' INSIDE SUBROUTINE DATA_CRYSTAL'')')
      I=0
      DO IM=1,NMODES(IPH)
        WRITE(10,'('' SHEAR MODULUS FOR MODE'',I3,'' IN PHASE'',I3,
     #             '' IS'',F12.3)') IM,IPH,MU_MODE(IM,IPH)
        WRITE(10,'('' N & B FOR MODE'',I3,'' IN PHASE'',I3)') IM,IPH
        DO IS=1,NSM(IM,IPH)
          I=I+1
          WRITE(10,'(3F10.3,3X,3F10.3)') (DNCA(J,I,IPH),J=1,3),
     #                                   (DBCA(J,I,IPH),J=1,3)
        ENDDO
      ENDDO

C     WRITE(10,*)
C     IS=0
C     DO IM=1,NMODES(IPH)
C       AUX55(:,:)=0.D0
C       DO ISX=1,NSM(IM,IPH)
C         IS=IS+1
C         DO I=1,5
C         DO J=1,5
C           AUX55(I,J)=AUX55(I,J)+SCHCA(I,IS,IPH)*SCHCA(J,IS,IPH)
C         ENDDO
C         ENDDO
C       ENDDO
C       WRITE(10,'(''  MATRIX sum(Mi*Mj) FOR MODE'',I3)') IM
C       WRITE(10,'(5F10.5)') ((AUX55(I,J),J=1,5),I=1,5)
C     ENDDO

      RETURN
      END

C
C *****************************************************************************
C     SUBROUTINE DATA_GRAIN        --->      VERSION 18/SEP/2020

C     READS GRAIN EULER ANGLES AND CALCULATES ROTATION MATRIX 'AG' WHICH
C     TRANSFORMS FROM CRYSTAL AXES TO SAMPLE AXES.
C     READS GRAIN ELLIPSOID EULER ANGLES AND SHAPE OF EACH PHASE (OR, IF
C     ISHAPE=3, OF EACH INDIVIDUAL GRAIN)
C     CALCULATES INITIAL DEFORMATION GRADIENT OF EACH PHASE ('FIJPH')
C     OR OF EACH GRAIN IF ISHAPE=3 ('FIJGR')
C *****************************************************************************

      SUBROUTINE DATA_GRAIN (IPH)

      USE VPSC8DIM

      DIMENSION AA(3,3),B(3,3),FIJEA(3,3),EULANG(3),AX(3)
      CHARACTER NOMEN*1

C -----------------------------------------------------------------
C     WRITES TEXTURE FILE INTO 'RUN_LOG.OUT' FILE
      WRITE(10,*)
      WRITE(10,'('' **** CRYST TEXTURE (FIRST FEW LINES) ****'')')
      DO IDUM=1,20
        READ(UNIT=UR2,END=98,FMT='(A)') PROSA
        WRITE(10,'(A)') PROSA
      ENDDO
   98 REWIND UR2
      WRITE(10,'(''    .........................'')')
      WRITE(10,'('' **** END OF CRYST TEXTURE DATA FILE ****'')')
      WRITE(10,*)
C -----------------------------------------------------------------

      READ(UR2,'(a)') PROSA
      READ(UR2,'(a)') PROSA
      READ(UR2,'(a)') PROSA
      READ(UR2,  *  ) NOMEN,NGRAIN

      NGR(IPH)=NGR(IPH-1)+NGRAIN
      if(ngr(iph).gt.ngrmx) then
        write(*,'('' number of grains exceeds dimension !!'')')
        write(*,'('' --> increase parameter NGRMX to'',i7)') ngr(iph)
        stop
      endif

C ***************************************************************************
C     READS EULER ANGLES, CONVERTS TO BUNGE NOTATION, CALCULATES ROT MATRIX

      TOTWGT=0.
      DO KKK=NGR(IPH-1)+1,NGR(IPH)
        READ(UR2,*) (EULANG(I),I=1,3),WGT(KKK)
        TOTWGT=TOTWGT+WGT(KKK)
        if(nomen.eq.'B' .or. nomen.eq.'b') then
          eul1=  eulang(1)
          eul2=  eulang(2)
          eul3=  eulang(3)
        else if(nomen.eq.'K' .or. nomen.eq.'k') then
          eul1= (eulang(1)-90.)
          eul2= -eulang(2)
          eul3=(-eulang(3)-90.)
        else if(nomen.eq.'R' .or. nomen.eq.'r') then
          eul1= (eulang(1)+90.)
          eul2=  eulang(2)           ! fixed 29/apr/02
          eul3= (eulang(3)-90.)
        else
          write(*,'(/,'' CANNOT IDENTIFY EULER ANGLE CONVENTION  !!'')')
          stop
        endif

C ***********************************************************************
C *** OPTIONAL: RANDOM PERTURBATION OF EULER ANGLES BY +/- 'randrange'
c       if(KKK.EQ.NGR(IPH-1)+1)
c     #          print *, ' RANDOM SHIFT OF EULER ANGLES ACTIVATED  !!'
c        randrange= 2.
c       eul1= eul1+ (2.*random2(jran)-1.)*randrange
c       eul2= eul2+ (2.*random2(jran)-1.)*randrange
c       eul3= eul3+ (2.*random2(jran)-1.)*randrange
C ***********************************************************************
C     EULER CALCULATES THE TRANSFORMATION MATRIX AA WHICH TRANSFORMS FROM
C     SAMPLE TO CRYSTAL. STORES AG, WHICH TRANSFORMS FROM CRYSTAL TO SAMPLE.

        CALL EULER(2,EUL1,EUL2,EUL3,AA)
        DO J=1,3
        DO K=1,3
          AG(J,K,KKK)=AA(K,J)
        ENDDO
        ENDDO

      ENDDO     ! END OF LOOP OVER GRAINS IN EACH PHASE

C     DOUBLE RENORMALIZATION OF WEIGHTS: FIRST NORMALIZE THE WEIGHTS WITHIN
C     EACH PHASE, NEXT RENORMALIZE TO THE VOLUME FRACTION OF THE PHASE.

      DO KKK=NGR(IPH-1)+1,NGR(IPH)
        WGT(KKK)=WGT(KKK)/TOTWGT*WPH(IPH)
      ENDDO

C ***********************************************************************
C       INITIAL F TENSOR AND EIGENVECTORS:

C     * IF ISHAPE=0 ASSUMES SAME INITIAL SHAPE (axisph) AND ORIENTATION
C       (eulerph) FOR ALL THE GRAIN ELLPSOIDS IN THE PHASE.
C       CALCULATES ESHELBY TENSOR FOR THE AVERAGE PHASE SHAPE.

C     * IF ISHAPE=1 ASSUMES SAME INITIAL SHAPE (axisph) AND ORIENTATION
C       (eulerph) FOR ALL THE GRAIN ELLIPSOIDS IN THE PHASE.
C       CALCULATES ESHELBY TENSOR FOR THE AVERAGE GRAIN PHASE.
C       CALCULATES INDIVIDUAL GRAIN ELLIPSOID SHAPE EVOLUTION.

C     * IF ISHAPE=2 ASSUMES SAME INITIAL SHAPE (axisph) AND ORIENTATION
C       (eulerph) FOR ALL THE GRAIN ELLIPSOIDS IN THE PHASE.
C       CALCULATES INDIVIDUAL GRAIN ELLIPSOID SHAPE EVOLUTION.
C       CALCULATES ESHELBY TENSOR FOR INDIVIDUAL GRAINS.

C     * IF ISHAPE=3 READS INDIVIDUAL INITIAL GRAIN AXES ORIENTATION
C       AND SHAPE FROM fileaxes (da,db,dc,ax(1),ax(2),ax(3)).
C       CALCULATES ESHELBY TENSOR FOR INDIVIDUAL GRAIN SHAPE AND
C       KEEPS TRACK OF INDIVIDUAL GRAIN SHAPE.

C       'eulerph' ANGLES OF (EA) WRT (SA).
C       'aa'    TRANSFORMS FROM (SA) TO (EA)
C       'fijea' COLUMNS ARE GRAIN AXES EXPRESSED IN ELLIPSOID SYSTEM
C       'fijph' COLUMNS ARE GRAIN AXES EXPRESSED IN SAMPLE SYSTEM
C ***********************************************************************

      da=EULERPH(1,iph)
      db=EULERPH(2,iph)
      dc=EULERPH(3,iph)

      call euler(2,da,db,dc,aa)

      do i=1,3
      do j=1,3
        FIJEA(I,J)=XID3(I,J)*AXISPH(0,I,IPH)
        B(I,J)         =AA(J,I)
		AXISPH(I,J,IPH)=AA(I,J)
      enddo
      enddo

      do j=1,3
      do i=1,3
        fijph(i,j,iph)=0.
        do k=1,3
        do l=1,3
          fijph(i,j,iph)=fijph(i,j,iph)+B(I,K)*B(J,L)*fijea(k,l)
        enddo
        enddo
      enddo
      enddo

C *** ISHAPE=1,2: ASSIGNS DEF GRAD OF PHASE (ELLIPSOID SHAPE & ORIENT) TO EVERY GRAIN

      IF(ISHAPE(IPH).EQ.1 .OR. ISHAPE(IPH).EQ.2) THEN
        DO KKK=NGR(IPH-1)+1,NGR(IPH)
          FIJGR(:,:,KKK)=FIJPH(:,:,IPH)
        ENDDO
      ENDIF

C *** ISHAPE=3: READ EULER ANGLES & ASPECT RATIOS OF EACH ELLIPSOID
C *** eul1,eul2,eul3: EULER ANGLES OF ELLIPSOID AXES WITH RESPECT TO SAMPLE
C *** AA TRANSFORMS FROM SAMPLE TO ELLIPSOID AXES. AAt FROM ELLIPSOID TO SAMPLE

      IF(ISHAPE(IPH).EQ.3) THEN

C -----------------------------------------------------------------
C     WRITES MORPHOLOGIC TEXTURE FILE INTO 'RUN_LOG.OUT' FILE
        WRITE(10,*)
        WRITE(10,'('' **** MORPH TEXTURE (FIRST 20 LINES) ****'')')
        DO IDUM=1,20
          READ(UNIT=UR3,END=99,FMT='(A)') PROSA
          WRITE(10,'(A)') PROSA
        ENDDO
   99   REWIND UR3
        WRITE(10,'(''    .........................'')')
        WRITE(10,'('' **** END OF MORPH TEXTURE DATA FILE ****'')')
        WRITE(10,*)
C -----------------------------------------------------------------

        DO IDUMMY=1,6
		  READ(UR3,'(A)') PROSA
        ENDDO

        DO KKK=ngr(iph-1)+1,ngr(iph)

          READ(UR3,*) eul1,eul2,eul3, xwgt, ax(1),ax(2),ax(3)

          CALL EULER(2,eul1,eul2,eul3,AA)
		  
          do i=1,3
            AXISGR(0,I,KKK)=AX(I)
            do j=1,3
              FIJEA(i,j)=XID3(I,J)*AX(I)
              B(I,J)    =AA(J,I)
		      AXISGR(I,J,IPH)=AA(I,J)
            enddo
          enddo

          do j=1,3
          do i=1,3
            FIJGR(i,j,KKK)=0.
            do k=1,3
            do l=1,3
              fijgr(i,j,KKK)=fijgr(i,j,KKK)+B(I,K)*B(J,L)*FIJEA(K,L)
            enddo
            enddo
          enddo
          enddo

        ENDDO      ! END OF LOOP OVER GRAINS
      ENDIF      ! END OF IF ISHAPE=3

      RETURN
      END

C **********************************************************************
C     SUBROUTINE DIFF_PLANES    -->   version 04/MAY/2018
C **********************************************************************
C *** IOPTION=0: FOR PHASE 'IPH' OPEN DIFFRACTION FILE & READ DIFF 
C                DIRECTIONS & Miller INDICES OF DIFFRACTING PLANES
C
C *** IOPTION=1: FOR PHASE 'IPH' DETERMINE GRAINS IN EACH HKL SET AND  
C                CALCULATE AVERAGE LATTICE STRAIN & INTENSITY FOR THE SET
C ***********************************************************************
C   NDIFF(IPH): # OF CRYSTALLOGRAPHIC DIFFRACTION PLANES hkl TO PROBE FOR 
C               IN PHASE 'IPH'
C   ND=1,NDIFF(IPH): INDEX OF INDIVIDUAL DIFFRACTION DIRECTION AND
C                    CORRESPONDING PLANE hkl
C   NFAMILY(ND,IPH): # OF EQUIV XTALL PLANES hkl TO SCAN FOR DIRECTION ND
C   DIFF_VS(3,ND,IPH): DIRECTION IN SAMPLE AXES TO PROBE FOR PLANE hkl 
C   DIFF_VC(3,NF,ND,IPH): DIRECTION IN CRYSTAL AXES OF EQUIV PLANES hkl 
C   EPS_W(ND,IPH): AVERAGE NORMAL STRAIN FOR ALL PLANES hkl ALONG DIR ND 
C ***********************************************************************

      SUBROUTINE DIFF_PLANES (IPH,ICRYSYM,IOPTION)

      USE VPSC8DIM
      USE DIFFRACT

      DIMENSION PS(3),ETELCS(6),RG(3,3),AA(3,3),N4(4)
      DIMENSION ISN(4),SN(3),SNEQ(3,24),ISB(4),SB(4)                   
      CHARACTER*12 PLANE_INDICES(100)
	  
      PI=4.D0*DATAN(1.D0)


C **********************************************************************
      IF (IOPTION.EQ.0) THEN 
C **********************************************************************
C *** READS INDICES OF DIFFRACTING PLANE AND ANGLES OF DIFFRACTION DIRECTION.
C *** GENERATES CRYSTALLOGRAPHICALLY EQUIVALENT PLANE NORMALS 'SNEQ' IN CAxes.
C *** CALCULATES DIFFRACTION DIRECTION 'SB' IN SAxes

        UR7=7
        OPEN(UNIT=UR7,FILE=FILEDIFF,STATUS='OLD')
        IUNIT=19+IPH
        OPEN(IUNIT,FILE='DIF_STR_PH'//CHAR(48+IPH)//'.OUT',
     #                   STATUS='UNKNOWN')
        IUNIT=22+IPH
        OPEN(IUNIT,FILE='DIF_WGT_PH'//CHAR(48+IPH)//'.OUT',
     #                   STATUS='UNKNOWN')
C -----------------------------------------------------------------
C     WRITES DIFFRACTION FILE INTO 'RUN_LOG.OUT' FILE
        WRITE(10,*)
        WRITE(10,'('' **** DIFFRACTION FILE ****'')')
        DO IDUM=1,100
          READ(UNIT=UR7,END=98,FMT='(A)') PROSA
          WRITE(10,'(A)') PROSA
        ENDDO

   98   REWIND UR7
        WRITE(10,'('' **** END OF DIFFRACTION FILE ****'')')
        WRITE(10,*)

C -----------------------------------------------------------------
C     WRITES MILLER INDICES OF DIFFRACTING PLANES IN OUTPUT FILES
        READ(UR7,'(A)') PROSA
        READ(UR7,'(A)') PROSA
        READ(UR7,*) NDIFF(IPH)
        READ(UR7,'(A)') PROSA
        READ(UR7,'(A)') PROSA
        DO NDIFFX=1,NDIFF(IPH)
          READ(UR7,'(A12)') PLANE_INDICES(NDIFFX)
        ENDDO
        IUNIT=19+IPH
        WRITE(IUNIT,'(''    EPSTOT11    EPSTOT22    EPSTOT33'',
C    #                ''    EPSTOT23    EPSTOT13    EPSTOT12'',
     #                ''      SBAR11      SBAR22      SBAR33   TEMP'',
C    #                ''      SBAR23      SBAR13      SBAR22'',
     #        6X,100A12)') (PLANE_INDICES(I),I=1,NDIFF(IPH))
        IUNIT=22+IPH
        WRITE(IUNIT,'(''    EPSTOT11    EPSTOT22    EPSTOT33'',
C    #                ''    EPSTOT23    EPSTOT13    EPSTOT12'',
     #                ''      SBAR11      SBAR22      SBAR33   TEMP'',
C    #                ''      SBAR23      SBAR13      SBAR22'',
     #        6X,100A12)') (PLANE_INDICES(I),I=1,NDIFF(IPH))

        REWIND UR7

C -----------------------------------------------------------------

        READ(UR7,'(A)') PROSA
        READ(UR7,'(A)') PROSA
        READ(UR7,*) NDIFF(IPH),ANGDETECTOR
        IF (NDIFF(IPH).GT.NDIFMX) THEN
          WRITE(*,'('' NUMBER OF DIFFRACTION PLANES IS'',I4,/,
     #    '' GREATER THAN CODE DIMENSION'',I4)') NDIFF(IPH),NDIFMX
          WRITE(*,*)
          WRITE(*,'(1H ,''STOP IN ROUTINE *** DIF_PLANES ***'')')
          STOP
        ENDIF
        READ(UR7,'(A)') PROSA
        READ(UR7,'(A)') PROSA

        DO ND=1,NDIFF(IPH)

          CALL CRYSTAL_SYMMETRY(3,UR7,ICRYSYM,ISN,SN,SNEQ,ISB,SB,NEQUIV)  
          NFAMILY(ND,IPH)=NEQUIV

          DO I=1,3
            DIFF_VS(I,ND,IPH)=SB(I)           ! DIFFACTION DIRECTION IN SA
            DO NF=1,NFAMILY(ND,IPH)
              DIFF_VC(I,NF,ND,IPH)=SNEQ(I,NF)       ! PLANE NORMAL IN CA
            ENDDO
          ENDDO
        ENDDO                              

        CLOSE(UNIT=UR7)

      ENDIF   ! END OF INITIALIZATION RUN FOR IOPTION=0

C **********************************************************************
      IF (IOPTION.EQ.1) THEN 
C **********************************************************************
C *** IDENTIFIES AND TAGS THOSE GRAINS THAT CONTRIBUTE TO EACH DIFFRACTION
C     DIRECTION (ND) IN THE PHASE (IPH) --> IGRSET(ND,NGRSETPH,IPH)

        COSANGDET=ABS(COS(ANGDETECTOR*PI/180.0)) 

        DO ND=1,NDIFF(IPH)

          NGRSETPH(ND,IPH)=0
          WGTSETPH(ND,IPH)=0.

          DO NG=NGR(IPH-1)+1,NGR(IPH)

            DO IPL=1,NFAMILY(ND,IPH)

              DO I=1,3
                PS(I)=0.0
                DO J=1,3
                  PS(I)=PS(I)+AG(I,J,NG)*DIFF_VC(J,IPL,ND,IPH)
                ENDDO
              ENDDO

              COSPSVS=0.
              DO I=1,3
	        COSPSVS=COSPSVS+DIFF_VS(I,ND,IPH)*PS(I)
	      ENDDO
              COSPSVS=ABS(COSPSVS)

              IF(COSPSVS.GE.COSANGDET) THEN
                NGRSETPH(ND,IPH)=NGRSETPH(ND,IPH)+1
                NGRX=NGRSETPH(ND,IPH)
                IGRSET(ND,NGRX,IPH)=NG
                WGTSETPH(ND,IPH)=WGTSETPH(ND,IPH)+WGT(NG)
              ENDIF

            ENDDO    ! END OF DO IPL=1,NFAMILY(ND)
          ENDDO      ! END OF DO NG =NGR(IPH-1)+1,NGR(IPH)
        ENDDO        ! END OF DO ND =1,NDIFF

      IPRINT=0
	  IF(IPRINT.EQ.1) THEN
        WRITE(19+iph,'(''* # OF GRAINS IN EACH SET: '',20i10)')
     #                  (ngrsetph(ndx,iph),ndx=1,ndiff(iph))
        WRITE(19+iph,'(''* VOL FRACTION OF EACH SET:'',20f10.5)')
     #                  (wgtsetph(ndx,iph),ndx=1,ndiff(iph))
      ENDIF

C **********************************************************************
C *** AVERAGE CRYSTALLOGRAPHIC STRAIN FOR GRAINS IN EACH SET

        DO ND=1,NDIFF(IPH)
          PARA_W(ND)=0.D0

          DO NG=1,NGRSETPH(ND,IPH)
            IGRSETX=IGRSET(ND,NG,IPH)

C *** WHAT FOLLOWS ARE DIFFERENT WAYS OF CALCULATING ELASTIC STRAIN IN EACH
C     GRAIN DEPENDING ON THE TYPE OF PROCESS BEING SIMULATED
C
C     1- FOR AN ELASTIC PROCES THE STRAIN 'EELGR' HAS BEEN CALCULATED IN THE
C        SUBROUTINE 'ELSC' AND CAN BE USED WITHOUT REPEATING THE CALCULATION
C     2- OPTIONALLY, ONE CAN CALCULATE 'EELGR' FROM THE KNOWN GRAIN STRESS 'SG'
C        USING THE GRAIN ELASTIC COMPLIANCE 'SELCS'.
C     3- TO ESTIMATE GRAIN ELASTIC STRAINS WHEN A VISCO-PLASTIC PROCESS IS
C        SIMULATED ONE HAS THE DEVIATORIC STRESS SG(1:5) FROM VP, BUT NOT THE
C        HYDROSTATIC COMPONENT SG(6). AN APPROXIMATION IS TO GET IT FROM
C        ELASTICITY, ASSUMING THAT THE HYDROSTATIC COMPONENT IN EACH GRAIN IS
C        NOT AFFECTED BY VISCO-PLASTICITY: FROM THE KNOWN AVERAGE STRESS 'SAV'
C        & THE STRESS CONCENTRATION TENSOR --> 'SG(6,NG)=BG_EL(6,J,NG)*SAV(J)'

C *** CASE 3 ***
          IF(INTERACTION.NE.-1) THEN
            SG(6,IGRSETX)=0.
            DO J=1,6
              SG(6,IGRSETX)=SG(6,IGRSETX)+BG_EL(6,J,IGRSETX)*SBAR(J)
            ENDDO
          ENDIF

C getting full stress tensor from elastic concentration equation
C           SG(:,IGRSETX)=0.
C           DO I=1,6
C           DO J=1,6
C              SG(I,IGRSETX)=SG(I,IGRSETX)+BG_EL(I,J,IGRSETX)*SBAR(J)
C           ENDDO
C           ENDDO

C *** CASE 2 & 3 ***
            EELGR(:,IGRSETX)=0.
            DO I=1,6
            DO J=1,6
              EELGR(I,IGRSETX)=EELGR(I,IGRSETX)+SELCS(I,J,IGRSETX)*
     #                                               SG(J,IGRSETX)
            ENDDO
            ENDDO

C *** CASE 1 & 2 & 3 *** CALCULATION OF AVERAGE NORMAL STRESS FOR DIFFRACTING SET 
          ISKIP=0
          IF(ISKIP.EQ.0) THEN
            CALL CHG_BASIS(EELGR(:,IGRSETX),AUX33,AUX66,AUX3333,1,6)
            EPS=0.D0
            DO I=1,3
            DO J=1,3
              EPS=EPS+DIFF_VS(I,ND,IPH)*DIFF_VS(J,ND,IPH)*AUX33(I,J) 
            ENDDO
            ENDDO
            PARA_W(ND)=PARA_W(ND)+EPS*WGT(IGRSETX)       ! weighted by grain fraction
          ENDIF

C     ALTERNATIVE CALCULATION OF AVERAGE NORMAL STRESS FOR DIFFRACTING SET 
          IF(ISKIP.EQ.1) THEN
            CALL CHG_BASIS(SG(:,IGRSETX),AUX33,AUX66,AUX3333,1,6)
            SIG=0.D0
            DO I=1,3
            DO J=1,3
              SIG=SIG+DIFF_VS(I,ND,IPH)*DIFF_VS(J,ND,IPH)*AUX33(I,J) 
            ENDDO
            ENDDO
            PARA_W(ND)=PARA_W(ND)+SIG*WGT(IGRSETX)       ! weighted by grain fraction
          ENDIF
		  
          ENDDO      ! END OF 'DO NG' OVER GRAINS IN SUBSET

          IF(WGTSETPH(ND,IPH).NE.0.D0) 
     #      EPS_W(ND,IPH)=PARA_W(ND)/WGTSETPH(ND,IPH)    ! renormalized using set fraction

        ENDDO      ! END OF 'DO ND' OVER DIFFRACTING PLANES
		
        IUNIT=19+IPH
        WRITE(IUNIT,'(6e12.4,F7.0,4x,100E12.4)')  
     #    EPSTOTc(1,1),EPSTOTc(2,2),EPSTOTc(3,3),
C    #    EPSTOTc(2,3),EPSTOTc(1,3),EPSTOTc(1,2),
     #    SBARc(1,1) ,SBARc(2,2) ,SBARc(3,3) , TEMP,
C    #    SBARc(2,3) ,SBARc(1,3) ,SBARc(1,2) ,
     #    (EPS_W(ND,IPH),ND=1,NDIFF(IPH))
       
        IUNIT=22+IPH
        WRITE(IUNIT,'(6e12.4,F7.0,4x,100E12.4)') 
     #    EPSTOTc(1,1),EPSTOTc(2,2),EPSTOTc(3,3),
C    #    EPSTOTc(2,3),EPSTOTc(1,3),EPSTOTc(1,2),
     #    SBARc(1,1) ,SBARc(2,2) ,SBARc(3,3) , TEMP,
C    #    SBARc(2,3) ,SBARc(1,3) ,SBARc(1,2) ,
     #    (WGTSETPH(ND,IPH),ND=1,NDIFF(IPH))

      ENDIF            !END OF IOPTION=1
C **********************************************************************

      RETURN
      END

c *****************************************************************************
      SUBROUTINE EIGEN_SORT(d,v,n,np)

C *** SORTS EIGENVALUES 'd' FROM LARGER TO SMALLER AND REARRANGES EIGENVECTORS 'v'

      INTEGER n,np
      REAL d(np),v(np,np)
      INTEGER i,j,k
      REAL p

      do i=1,n-1
        k=i
        p=d(i)
        do j=i+1,n
          if(d(j).ge.p)then
            k=j
            p=d(j)
          endif
        enddo
        if(k.ne.i)then
          d(k)=d(i)
          d(i)=p
          do j=1,n
            p=v(j,i)
            v(j,i)=v(j,k)
            v(j,k)=p
          enddo
        endif
      enddo
      return
      END

c *****************************************************************************
      SUBROUTINE EIGEN_VAL (a,n,np,d,v,nrot,ier)

C *** CALCULATES EIGENVALUES AND EIGENVECTORS OF SYMMETRIC MATRIX A=FIJ*FIJ_t

      INTEGER n,np,nrot,NMAX
      REAL a(np,np),d(np),v(np,np)
      PARAMETER (NMAX=500)
      INTEGER i,ip,iq,j
      REAL c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)

      do ip=1,n
        do iq=1,n
          v(ip,iq)=0.
        enddo
        v(ip,ip)=1.
      enddo
      do ip=1,n
        b(ip)=a(ip,ip)
        d(ip)=b(ip)
        z(ip)=0.
      enddo
      nrot=0
      do i=1,50
        sm=0.
        do ip=1,n-1
          do iq=ip+1,n
            sm=sm+abs(a(ip,iq))

          enddo
        enddo

        if(sm.eq.0.)then
        ier=0
        return
        endif

        if(i.lt.4)then
          tresh=0.2*sm/n**2
        else
          tresh=0.
        endif
        do ip=1,n-1
          do iq=ip+1,n
            g=100.*abs(a(ip,iq))
            if((i.gt.4).and.(abs(d(ip))+
     *  g.eq.abs(d(ip))).and.(abs(d(iq))+g.eq.abs(d(iq))))then
              a(ip,iq)=0.
            else if(abs(a(ip,iq)).gt.tresh)then
              h=d(iq)-d(ip)
              if(abs(h)+g.eq.abs(h))then
                t=a(ip,iq)/h
              else
                theta=0.5*h/a(ip,iq)
                t=1./(abs(theta)+sqrt(1.+theta**2))
                if(theta.lt.0.)t=-t
              endif
              c=1./sqrt(1+t**2)
              s=t*c
              tau=s/(1.+c)
              h=t*a(ip,iq)
              z(ip)=z(ip)-h
              z(iq)=z(iq)+h
              d(ip)=d(ip)-h
              d(iq)=d(iq)+h
              a(ip,iq)=0.
              do j=1,ip-1
                g=a(j,ip)
                h=a(j,iq)

                a(j,ip)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
              enddo
              do j=ip+1,iq-1
                g=a(ip,j)
                h=a(j,iq)
                a(ip,j)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
              enddo
              do j=iq+1,n
                g=a(ip,j)
                h=a(iq,j)
                a(ip,j)=g-s*(h+g*tau)
                a(iq,j)=h+s*(g-h*tau)
              enddo
              do j=1,n
                g=v(j,ip)
                h=v(j,iq)

                v(j,ip)=g-s*(h+g*tau)
                v(j,iq)=h+s*(g-h*tau)
              enddo
              nrot=nrot+1
            endif
          enddo
        enddo
        do ip=1,n
          b(ip)=b(ip)+z(ip)
          d(ip)=b(ip)
          z(ip)=0.
        enddo
      enddo
c      pause 'too many iterations in eigen_val'
      ier=1

      return
      END

C ************************************************************************
C     SUBROUTINE ELSC      -->      VERSION 08/JUL/2020
C
C     CALCULATES ITERATIVELY SELF-CONSISTENT ELASTIC MODULI OF POLYCRYSTAL
C         Csc=<Cg:Ag>:<Ag>^(-1)
C         with Ag =<(Cg+Cg~)^(-1)>:(Csc+Cg~)
C         with Cg~=Csc:(Eshg^(-1)-I)
C     CALCULATES SELF-CONSISTENT THERMAL MODULI OF POLYCRYSTAL
C         ATHsc=Csc^(-1)::<Ag^(-1)>^(-1):=<Ag^(-1):Cg:ATHg>
C
C     ITERATIVE PROCEDURE. USES b-BASIS REPRESENTATION OF TENSORS.
C       IOPTION=0 : CALCULATES UPPER BOUND ELASTIC AND THERMAL TENSORS
C                  'CELUB' AND 'ATHUB'
C       IOPTION=1 : CALCULATES SELF CONSISTENT ELASTIC AND THERMAL TENSORS
C                  'CELAV' AND 'ATHAV'
C       IOPTION=2 : CALCULATES STRESS IN GRAINS INDUCED BY ELASTIC LOAD
C                   AND/OR TEMPERATURE CHANGE
C ************************************************************************

      SUBROUTINE ELSC (IOPTION)

      USE VPSC8DIM

      DIMENSION C4XC(3,3,3,3),C4XS(3,3,3,3),S4XC(3,3,3,3),
     #          S4XS(3,3,3,3),A2XC(3,3),A2XS(3,3)
      DIMENSION CELUB(6,6),CELUBv(6,6),CELAVx(6,6),CTILDE(6,6)
	  DIMENSION CELLB(6,6),CELLBv(6,6),SELLB(6,6),SELUB(6,6),
     #          ATHUB(6),ATHUBv(6),ATHLB(6),ATHLBv(6)
      DIMENSION C4SA(3,3,3,3),E4SA(3,3,3,3),R4SA(3,3,3,3),
     #          AXB(3),EIGB(3,3),ESHELv(6,6)
      DIMENSION AGX(6,6),AGTX(6,6),DGX(6,6),TGX(6,6),
     #          TG_ELx(6,6,NGRPEL),DG_ELx(6,6,NGRPEL),ATHAUX(6),
     #          AGAV(6,6),AGAVINV(6,6),AGTAV(6,6),AGTAVINV(6,6),
     #          CGAG(6,6),CGAGAV(6,6)

      DIMENSION dSBAR(3,3),dEBARel(3,3)
      DIMENSION DELGR(6,NGRPEL),DELAV(6),DELAVc(3,3),DTHAV(6)
      DIMENSION EELcINCR(3,3),SAVcINCR(3,3),SELAVc(3,3,3,3)

      IPRINT=0     ! controls detailed printing

      IF(IOPTION.EQ.0) THEN
        OPEN (777,file='EL-TH-MODULI.OUT',status='unknown')
        WRITE(777,'(5x,''Evm'',9x,''C11  C22  C33  C23  C13  C12   '',
     #                 ''C44  C55  C66'',67x,''A1  A2  A3'')')
c       open (777,file='pressure_test.out',status='unknown')
c       write(777,'(''PHASE'',8X,''SAVPH'',72X,''SAVPHDEV'')')  
      ENDIF

C ************************************************************************
      IF(IOPTION.EQ.0 .OR. IOPTION.EQ.1) THEN 
C ************************************************************************
C *** ROTATES SINGLE CRYSTAL STIFFNESS, COMPLIANCE AND THERMAL EXPANSION 
C     TENSORS (CELCC, SELCC, ATHCC) TO SAMPLE AXES.
C *** CALCULATES UPPER BOUND (VOIGT) STIFFNESS 'CELUB' FOR THE POLYCRYSTAL.
C *** CALCULATES UPPER BOUND THERMAL EXPANSION TENSOR 'ATHUB' FOR THE PX.
C *** BOUNDS USED AS INITIAL GUESS FOR SELF-CONSISTENT CALCULATION
C ************************************************************************

      CELUB(:,:)=0.
      SELLB(:,:)=0.
      ATHUB(:)  =0.
      ATHLB(:)  =0.

      DO IPH=IPHBOT,IPHTOP

C *** TRANSFORMS (ATHCC,CELCC,SELCC) FROM b-BASIS TO CARTESIAN (A2XC,C4XC,S4XC)
C     AND ROTATES INTO SAMPLE AXES (A2XS,C4XS,S4XS) FOR EVERY GRAIN

        CALL CHG_BASIS(ATHCC(:,IPH),A2XC,AUX66,AUX3333,1,6) 
        CALL CHG_BASIS(AUX6,AUX33,CELCC(:,:,IPH),C4XC,3,6)
        CALL CHG_BASIS(AUX6,AUX33,SELCC(:,:,IPH),S4XC,3,6)  

        DO KKK=NGR(IPH-1)+1,NGR(IPH)

          DO I=1,3
          DO J=1,3
            A2XS(I,J)=0.
            DO I1=1,3
            DO J1=1,3
              A2XS(I,J)=A2XS(I,J)+AG(I,I1,KKK)*AG(J,J1,KKK)*A2XC(I1,J1)                
            ENDDO 
            ENDDO
          ENDDO
          ENDDO

          DO 85 I=1,3
          DO 85 J=1,3
          DO 85 K=1,3
          DO 85 L=1,3
            C4XS(I,J,K,L)=0.
            S4XS(I,J,K,L)=0.
            DO 80 I1=1,3
            DO 80 J1=1,3           
            DO 80 K1=1,3
            DO 80 L1=1,3
              C4XS(I,J,K,L)=C4XS(I,J,K,L)+AG(I,I1,KKK)*AG(J,J1,KKK)*
     #                      AG(K,K1,KKK)*AG(L,L1,KKK)*C4XC(I1,J1,K1,L1)
              S4XS(I,J,K,L)=S4XS(I,J,K,L)+AG(I,I1,KKK)*AG(J,J1,KKK)*
     #                      AG(K,K1,KKK)*AG(L,L1,KKK)*S4XC(I1,J1,K1,L1)
   80       ENDDO 
   85     ENDDO

C *** FROM CARTESIAN TO B-BASIS
          CALL CHG_BASIS(ATHCS(:,KKK),A2XS,AUX66,AUX3333,2,6)
          CALL CHG_BASIS(AUX6,AUX33,CELCS(:,:,KKK),C4XS,4,6)
          CALL CHG_BASIS(AUX6,AUX33,SELCS(:,:,KKK),S4XS,4,6)

          DO I=1,6
          DO J=1,6
            ATHUB(I)=ATHUB(I)+CELCS(I,J,KKK)*ATHCS(J,KKK)*WGT(KKK)
		  ENDDO
		  ENDDO

		  ATHLB(:)  =ATHLB(:)  +ATHCS(:,KKK)  *WGT(KKK)		  
          CELUB(:,:)=CELUB(:,:)+CELCS(:,:,KKK)*WGT(KKK)
          SELLB(:,:)=SELLB(:,:)+SELCS(:,:,KKK)*WGT(KKK)
		  
        ENDDO        ! END OF DO KKK
      ENDDO          ! END OF DO IPH


      SELUB(:,:)=CELUB(:,:)	  
      CALL LU_INVERSE(SELUB,6)
      AUX6(:) =ATHUB(:)
      ATHUB(:)=0.
      DO I=1,6
      DO J=1,6
        ATHUB(I)=ATHUB(I)+SELUB(I,J)*AUX6(J)
	  ENDDO
	  ENDDO

      IF(IOPTION.EQ.0) THEN

C *** SAVES UPPER BOUND IN B-BASIS FOR USING AS NEXT GUESS
      CELAV(:,:)=CELUB(:,:) 

      CELLB(:,:)=SELLB(:,:)	  
      CALL LU_INVERSE(CELLB,6)

C *** TRANSFORM FROM B-BASIS TO CARTESIAN TO VOIGT FOR OUTPUT
      CALL CHG_BASIS(ATHUB ,AUX33,AUX66,AUX3333,1,6)
      CALL VOIGT    (ATHUBv,AUX33,AUX66,AUX3333,2)
      CALL CHG_BASIS(AUX6,AUX33,CELUB ,AUX3333,3,6)
      CALL VOIGT    (AUX6,AUX33,CELUBv,AUX3333,4)	
	  
      CALL CHG_BASIS(ATHLB ,AUX33,AUX66,AUX3333,1,6)
      CALL VOIGT    (ATHLBv,AUX33,AUX66,AUX3333,2)
      CALL CHG_BASIS(AUX6,AUX33,CELLB ,AUX3333,3,6)
      CALL VOIGT    (AUX6,AUX33,CELLBv,AUX3333,4)	

      WRITE(10,*)
      WRITE(10,'('' UB THERMAL TENSOR (VOIGT NOTATION)'')')
      WRITE(10,'(6E12.4)') (ATHUBv(I),I=1,6)
      WRITE(10,*)
      WRITE(10,'('' UB ELASTIC STIFFNESS (VOIGT NOTATION)'')')
      WRITE(10,'(6E12.4)') ((CELUBv(I,J),J=1,6),I=1,6)

      WRITE(10,*)
      WRITE(10,'('' LB THERMAL TENSOR (VOIGT NOTATION)'')')
      WRITE(10,'(6E12.4)') (ATHLBv(I),I=1,6)
      WRITE(10,*)
      WRITE(10,'('' LB ELASTIC STIFFNESS (VOIGT NOTATION)'')')
      WRITE(10,'(6E12.4)') ((CELLBv(I,J),J=1,6),I=1,6)
      WRITE(10,*)

      ENDIF

      ENDIF      ! END OF IF IOPTION = 0 or 1

C *************************************************************************
C    IF IOPTION=1 : USES PREVIOUS CONVERGED MODULI 'CELAV' AS INITIAL 
C    GUESS FOR A SELF-CONSISTENT CALCULATION OF THE ELASTIC STIFFNESS AND  
C    THE THERMAL MODULI
C *************************************************************************
      IF(IOPTION.EQ.1) THEN 
C *************************************************************************

      WRITE(*,*)
      IT=0
      RER=2*1.E-05      ! hardwired for elasticity
      ITMAXINT=100      ! hardwired for elasticity
      DO WHILE(RER.GT.1.E-05 .AND. IT.LT.ITMAXINT)

        IT=IT+1
        CGAGAV(:,:)  =0.
        AGAV(:,:)    =0.
        AGTAV(:,:)   =0.
        ATHAUX(:)    =0.

        CALL CHG_BASIS(AUX6,AUX33,CELAV,C4SA,3,6)

        KGX=1
        DO IPH=IPHBOT,IPHTOP
		
C *** SAME SHAPE FOR ALL GRAINS IN THE PHASE 
        IF(ISHAPE(IPH).LE.1) THEN
          AXB(:)     =AXISPH(0,:,IPH)
          EIGB(1:3,:)=AXISPH(1:3,:,IPH)
        ENDIF

          DO KKK=NGR(IPH-1)+1,NGR(IPH)

C *** DIFFERENT SHAPE FOR EACH GRAIN IN THE PHASE 
          IF(ISHAPE(IPH).GE.2) THEN 
            AXB(:)     =AXISGR(0,:,KGX)
            EIGB(1:3,:)=AXISGR(1:3,:,KGX)
          ENDIF

C *** ISHAPE=0,1 --> SAME SHAPE FOR ALL GRAINS IN THE PHASE 
C *** ISHAPE=2,3 -->DIFFERENT SHAPE FOR EACH GRAIN IN THE PHASE 

          IF(ISHAPE(IPH).GE.2 .OR. KKK.EQ.NGR(IPH-1)+1) THEN

            CALL ESHELBY_TENSOR (AXB,EIGB,C4SA,E4SA,R4SA,1)      ! R4SA is not used

            CALL CHG_BASIS(AUX6,AUX33,AUX66,E4SA,4,6)
            CALL LU_INVERSE(AUX66,6)                   ! ESHg^(-1)

            DO I=1,6
            DO J=1,6
              CTILDE(I,J)=0.
              DO K=1,6
               CTILDE(I,J)=CTILDE(I,J)+CELAV(I,K)*(AUX66(K,J)-XID6(K,J))
              ENDDO
            ENDDO
            ENDDO

          ENDIF     

            DO I=1,6
            DO J=1,6
              DGX(I,J)=CELCS(I,J,KGX)+CTILDE(I,J)      ! (Cg+C~)
            ENDDO
            ENDDO
            CALL LU_INVERSE(DGX,6)
            DO I=1,6
            DO J=1,6
              AGX(I,J) =0. 
              AGTX(I,J)=0. 
              TGX(I,J) =0.
              DO K=1,6
                AGX(I,J) =AGX(I,J) +DGX(I,K)*(CELAV(K,J)+CTILDE(K,J))    ! Ag =(Cg+C~)^(-1):(Cav+C~)
				AGTX(I,J)=AGTX(I,J)+(CELAV(I,K)+CTILDE(I,K))*DGX(K,J)    ! AgT=(Cav+C~):(Cg+C~)^(-1)
                TGX(I,J)=TGX(I,J)+DGX(I,K)*CTILDE(K,J)                   ! TgX=(Cg+C~)^(-1): C~
              ENDDO
            ENDDO
            ENDDO

C      IF(kgx.eq.77) THEN         ! testing if AG=AGT (symmetric)
C        write(10,'('' DGX  '',6e12.4)') DGX
C        write(10,'('' C~   '',6e12.4)') CTILDE
C        write(10,'('' AG   '',6e12.4)') AGX
C        write(10,'('' AGT  '',6e12.4)') AGTX
C      ENDIF

            DO I=1,6
            DO J=1,6
              CGAG(I,J)=0.
              DO K=1,6
                CGAG(I,J)=CGAG(I,J)+CELCS(I,K,KGX)*AGX(K,J)
              ENDDO
              CGAGAV(I,J )=CGAGAV(I,J)+CGAG(I,J) *WGT(KGX)      ! <Cg:Ag>
              AGAV(I,J)   =AGAV(I,J)  +AGX(I,J)  *WGT(KGX)      ! <Ag>
              AGTAV(I,J)  =AGTAV(I,J) +AGTX(I,J) *WGT(KGX)      ! <AgT>
            ENDDO
            ENDDO

C *** CALCULATES CONCENTRATION TENSORS RELATED TO THE AVERAGE THERMAL MODULI
            DO I=1,6
            DO J=1,6
              TG_ELx(I,J,KGX)=0.      ! STRESS CONCENTRATION TENSOR 'TG_EL'
			  DO K=1,6
                TG_ELx(I,J,KGX)=TG_ELx(I,J,KGX)+CELCS(I,K,KGX)*TGX(K,J)
              ENDDO
              DO K=1,6
                ATHAUX(I)=ATHAUX(I)+AGTX(I,J)*CELCS(J,K,KGX)* 
     #                                ATHCS(K,KGX) *WGT(KGX)
C --> AG & AGT (TRANSPOSE) EMPIRICALLY PROVED TO BE EQUAL  !!! 
C                ATHAUX(I)=ATHAUX(I)+AGX(J,I)*CELCS(J,K,KGX)*      
C     #                                ATHCS(K,KGX) *WGT(KGX)
              ENDDO

              AG_EL(I,J,KGX) =AGX(I,J)      ! STRAIN CONCENTRATION TENSOR 'AG_EL'
              DG_ELx(I,J,KGX)=DGX(I,J)      ! AUX TO STRAIN CONCTN TENSOR 'DG_EL'
            ENDDO
            ENDDO

            KGX=KGX+1
          ENDDO      ! END OF DO KKK
        ENDDO      ! END OF DO IPH

        AGAVINV(:,:)=AGAV(:,:)
        CALL LU_INVERSE(AGAVINV,6)
        AGTAVINV(:,:)=AGTAV(:,:)
        CALL LU_INVERSE(AGTAVINV,6)
		
        DO I=1,6
        DO J=1,6
          CELAVx(I,J)=0.
          DO K=1,6
            CELAVx(I,J)=CELAVx(I,J)+CGAGAV(I,K)*AGAVINV(K,J)     ! <Cg:Ag>:<Ag>^(-1)
          ENDDO
        ENDDO
        ENDDO

        RER=TMISMATCH(CELAV,CELAVx,6)

C      WRITE(10,*)
C      WRITE(10,'(''CELAVx ELASTIC STIFFNESS (b-basis NOTATION)'')')
C      WRITE(10,'(6E12.4)') ((CELAVx(I,J),J=1,6),I=1,6)
C      WRITE(10,'(''ITER'',I4,''   REL ERROR'',E11.3)') IT,RER
C      WRITE( *,'(''ITER'',I4,''   REL ERROR'',E11.3)') IT,RER

C *** SAVES THE SYMMETRIZED TENSOR AS THE NEW GUESS
        DO I=1,6
          DO J=I,6
            CELAV(I,J)=0.5*(CELAVx(I,J)+CELAVx(J,I))
            CELAV(J,I)=CELAV(I,J)
          ENDDO
        ENDDO

        IF (IT. EQ. ITMAXINT) THEN
          WRITE( *,'(''  NO ELASTIC SELF-CONSISTENT CONVERGENCE'',
     #    '' --> ITERATION'',I4,''  RELATIVE ERR>1.e-5'',E11.3)') IT,RER
          STOP
        ENDIF

      ENDDO      ! END OF (DO..WHILE)

      SELAV(:,:) = CELAV(:,:)
      CALL LU_INVERSE(SELAV,6)

C *** CALCULATES THE PX THERMAL EXPANSION TENSOR 'ATHAV'
      DO I=1,6
        ATHAV(I)=0.
        DO J=1,6
        DO K=1,6
          ATHAV(I)=ATHAV(I)+SELAV(I,J)*AGTAVINV(J,K)*ATHAUX(K) 
C --> AG TRANSPOSE IS EMPIRICALLY SHOWN TO BE EQUAL TO AGT, AND GIVES THE SAME ANSWER
C         ATHAV(I)=ATHAV(I)+SELAV(I,J)*AGAVINV(K,J)*ATHAUX(K)
        ENDDO
        ENDDO


      ENDDO      
      CALL CHG_BASIS(ATHAV,ATHAVc,AUX66,AUX3333,1,6)
	  
C *** CALCULATES THE STRESS CONCENTRATION TENSORS 'TG_EL' AND 'DG_EL'
C     (REQUIRES PREVIOUS CALCULATION OF 'ATHAV') AND 'BG_EL'
        KGX=1
        DO IPH=IPHBOT,IPHTOP
        DO KKK=NGR(IPH-1)+1,NGR(IPH)

          TG_EL(:,KGX)=0.
          DO I=1,6
            DO K=1,6
			  TG_EL(I,KGX)=TG_EL(I,KGX)+TG_ELx(I,K,KGX)*
     #                                 (ATHAV(K)-ATHCS(K,KGX))
            ENDDO
          ENDDO
		  
          DG_EL(:,KGX)=0.
          DO I=1,6
		    DO K=1,6
            DO L=1,6
			  DG_EL(I,KGX)=DG_EL(I,KGX)+DG_ELx(I,K,KGX)*
     #             (CELCS(K,L,KGX)*ATHCS(L,KGX)-CELAV(K,L)*ATHAV(L))
            ENDDO
            ENDDO
          ENDDO
		  
          BG_EL(:,:,KGX)=0.		  
          DO I=1,6
          DO J=1,6
            DO K=1,6
            DO L=1,6
              BG_EL(I,J,KGX)=BG_EL(I,J,KGX)+
     #              CELCS(I,K,KGX)*AG_EL(K,L,KGX)*SELAV(L,J)           
            ENDDO
            ENDDO
          ENDDO
          ENDDO
		  
          KGX=KGX+1
        ENDDO      ! END OF DO KKK
        ENDDO      ! END OF DO IPH

C *** CALCULATES VOIGT REPRESENTATION OF TENSORS
      CALL VOIGT    (AUX6,AUX33,ESHELv,E4SA,4)
      CALL CHG_BASIS(AUX6,AUX33,CELAV ,AUX3333,3,6)
      CALL VOIGT    (AUX6,AUX33,CELAVv,AUX3333,4)
      CALL CHG_BASIS(AUX6,AUX33,SELAV ,AUX3333,3,6)
      CALL VOIGT    (AUX6,AUX33,SELAVv,AUX3333,4)
      CALL CHG_BASIS(ATHAV, ATHAVc,AUX66,AUX3333,1,6)
      CALL VOIGT    (ATHAVv,ATHAVc,AUX66,AUX3333,2)

	  IF(IPRINT.EQ.1) THEN
	  
C *** WRITE VOIGT MATRIX REPRESENTATION OF ESHELBY TENSOR
      WRITE(10,*)
      WRITE(10,*) ' ELASTIC ESHELBY TENSOR (VOIGT NOTATION)'
      WRITE(10,'(6F12.5)') ((AUX66(I,J),J=1,6),I=1,6)

C *** WRITE VOIGT MATRIX REPRESENTATION OF STIFFNESS TENSOR
      WRITE(10,*)
      WRITE(10,'(''SC ELASTIC STIFFNESS (VOIGT NOTATION)'')')
      WRITE(10,'(6E12.4)') ((CELAVv(I,J),J=1,6),I=1,6)

C --> WRITE VOIGT MATRIX REPRESENTATION OF COMPLIANCE TENSOR
      WRITE(10,*)
      WRITE(10,'(''SC ELASTIC COMPLIANCE (VOIGT NOTATION)'')')
      WRITE(10,'(6E12.4)') ((SELAVv(I,J),J=1,6),I=1,6)

C --> WRITE VOIGT MATRIX REPRESENTATION OF THERMAL TENSOR
      WRITE(10,*)
      WRITE(10,'(''SC THERMAL EXPANSION TENSOR (VOIGT NOTATION)'')')
      WRITE(10,'(6E12.4)') (ATHAVv(I),I=1,6)

	  ENDIF      ! END OF IF(IPRINT.eq.1)
	  
      WRITE(777,'(F12.5,3X,9E12.4,3X,3E12.4)') EPSVM,
     #     CELAVv(1,1),CELAVv(2,2),CELAVv(3,3),CELAVv(2,3),CELAVv(1,3),
     #     CELAVv(1,2),CELAVv(4,4),CELAVv(5,5),CELAVv(6,6),
     #     ATHAVv(1)  ,ATHAVv(2)  ,ATHAVv(3)
      
      ENDIF      ! END OF IF(IOPTION.EQ.1)

C ************************************************************************
C *** IOPTION=2: CALCULATES STRESS, ELASTIC STRAIN, THERMAL STRAIN IN 
C                POLYCRYSTAL AND GRAINS BASED ON LOAD CONDITIONS
C ************************************************************************
      IF(IOPTION.EQ.2) THEN 
C ************************************************************************

C *** THE FOLLOWING BC's ARE READ or SET INSIDE SUBROUTINE LOAD_CONDITIONS
C     -->  DBARc(:,:) , SBARc(:,:) , TIME_INCR , TEMP_INCR 

C *** IN SUBROUTINE STATE_6x6 TENSORS 'DBAR & SBAR & SELAV' ENTER IN 
C     CARTESIAN, GET INTERNALLY CONVERTED TO VOIGT, AND EXIT IN CARTESIAN.
C *** IMPOSED 'DBAR & SBAR' COMPONENTS MAY BE MIXED.
C *** THE INCREMENTAL EQUATION TO BE SOLVED IS:
C     -->  delta_E_el = S_el*delta_Sig = delta_E_tot-delta_E_th

C *** THE SOLUTION IS TRIVIAL IF ALL STRESS COMPONENTS (ZERO OR NON ZERO)
C     ARE APPLIED AS BOUNDARY CONDITION (IN WHICH CASE delta_Sig=0) 		

        ISBARSUM=ISBARv(1)+ISBARv(2)+ISBARv(3)+ISBARv(4)+
     #           ISBARv(5)+ISBARv(6)
        IF (ISBARSUM.EQ.6) THEN
		  DBARc(:,:)=ATHAVc(:,:)*TEMP_INCR/TIME_INCR
		  dSBAR(:,:)=0.
          dEBARel(:,:)=0.
		ELSE
          dEBARel(:,:)=DBARc(:,:)*TIME_INCR-ATHAVc(:,:)*TEMP_INCR     ! known comps enforced
		  CALL VOIGT(AUX6,AUX33,SELAVv,SELAVc,3)
          CALL STATE_6x6 (IDBARv,dEBARel,ISBARv,dSBAR,SELAVc)
		  DBARc(:,:)=(dEBARel(:,:)+ATHAVc(:,:)*TEMP_INCR) /TIME_INCR
		ENDIF
        SBARc(:,:)=SBARc(:,:)+dSBAR(:,:)

		CALL CHG_BASIS(DBAR,DBARc,AUX66,AUX3333,2,6)
        CALL CHG_BASIS(SBAR,SBARc,AUX66,AUX3333,2,6)

        KGX=1
        DO IPH=IPHBOT,IPHTOP
        DO KKK=NGR(IPH-1)+1,NGR(IPH)

          DO I=1,6
          SG(I,KGX)=TG_EL(I,KGX)*(TEMP+TEMP_INCR-TEMP_INI)      ! TEMP is updated later in main
          DO J=1,6
            SG(I,KGX)=SG(I,KGX)+BG_EL(I,J,KGX)*SBAR(J)
          ENDDO
          ENDDO

		  EELGR(:,KGX)=0.
		  DO I=1,6
          DO J=1,6
		   EELGR(I,KGX)=EELGR(I,KGX)+SELCS(I,J,KGX)*SG(J,KGX)
          ENDDO
          ENDDO

          KGX=KGX+1
        ENDDO      ! END OF DO KKK
        ENDDO      ! END OF DO IPH

C *** THE FOLLOWING USED FOR ANALYZING PRESSURE IN 2-phase AGGREGATE

        SAVPH(:,:)=0.
        SAVPHDEV(:,:)=0.
        KGX=1
        DO IPH=IPHBOT,IPHTOP
        DO KKK=NGR(IPH-1)+1,NGR(IPH)
          WGTXX= WGT(KGX)/WPH(IPH) 
          SAVPH(:,IPH)   =SAVPH(:,IPH)   +sg(:,KGX)*wgtxx
          SAVPHDEV(:,IPH)=SAVPHDEV(:,IPH)+sg(:,KGX)**2*wgtxx
          KGX=KGX+1
        ENDDO      ! END OF DO KKK
        SAVPHDEV(:,IPH)=sqrt(ABS(SAVPHDEV(:,IPH)-SAVPH(:,IPH)**2))

	    IPRINT=0
        IF (IPRINT.EQ.1) THEN
          write(777,'(''SAVPH for IPH   ='',i2,5x,6e12.4)')
     #                iph,savph(:,iph)
          write(777,'(''SAVPHDEV for IPH='',i2,5x,6e12.4)')
     #                iph,savphdev(:,iph)
C --> case of one grain per phase
		  CALL CHG_BASIS(SG(:,IPH),AUX33,AUX66,AUX3333,1,6)
          write(777,'(''SGRc for IPH   ='',i2,/,(5x,3e12.4) )')
     #                iph,AUX33
		  CALL CHG_BASIS(EELGR(:,IPH),AUX33,AUX66,AUX3333,1,6)
          write(777,'(''EELGRc for IPH   ='',i2,/,(5x,3e12.4) )')
     #                iph,AUX33
        ENDIF

        ENDDO      ! END OF DO IPH

        IPRINT=0
        NGBOT=1
        NGTOP=1
        IF (IPRINT.EQ.1) THEN
          DO KGX=NGBOT,NGTOP
		  ELENERGY=0.
		  DO I=1,6
            ELENERGY=ELENERGY+0.5*SG(I,kgx)*EELGR(I,kgx)
          ENDDO
          CALL CHG_BASIS(EELGR(:,kgx),AUX33,AUX66,AUX3333,1,6)
          CALL VOIGT    (AUX6,        AUX33,AUX66,AUX3333,2)
          WRITE(777,'(''EL ENERGY & EL TENSOR for grain #'',
     #                I3,E12.4,3x,6F11.6)') KGX,ELENERGY,AUX6
          ENDDO
        ENDIF

      ENDIF      ! END OF IF(IOPTION.EQ.2)

      RETURN
      END

C **********************************************************************
C     SUBROUTINE ESHELBY_TENSOR      --->      VERSION 02/APR/2023

C *** ESHELBY CALCULATION FOR EVERY PHASE OR FOR EVERY GRAIN
C *** ROTATES STIFNESS TO ELLIPSOID AXES 'EIGB'
C *** CALLS ESHELBY TO CALCULATE DISTORTION AND ROTATION ESHELBY TENSORS
C     'E4GA' & 'R4GA' IN ELLIPSOID AXES
C *** ROTATES THEM BACK TO SAMPLE AXES 'E4SA' & 'R4SA' 
C **********************************************************************

      SUBROUTINE ESHELBY_TENSOR (AXB,EIGB,C4SA,E4SA,R4SA,IOPTION)

      USE C4GAFLU
	  
      DIMENSION AUX6(6),AUX33(3,3),AUX66(6,6),AUX3333(3,3,3,3)
      DIMENSION C4SA(3,3,3,3),C4GA(3,3,3,3)
      DIMENSION E4SA(3,3,3,3),E4GA(3,3,3,3),R4SA(3,3,3,3),R4GA(3,3,3,3)
      DIMENSION AXB(3),EIGB(3,3)

C *** ROTATION OF STIFFNESS 'C4' FROM SAMPLE TO ELLIPSOID AXES
            DO I=1,3
            DO J=1,3
            DO M=1,3
            DO N=1,3
              DUMMY=0.
              DO I1=1,3
              DO J1=1,3
              DO M1=1,3
              DO N1=1,3
                DUMMY=DUMMY+EIGB(I1,I)*EIGB(J1,J)*EIGB(M1,M)
     #               *EIGB(N1,N)*C4SA(I1,J1,M1,N1)
              ENDDO
              ENDDO
              ENDDO
              ENDDO
              C4GA(I,J,M,N)=DUMMY
              C4GA_FLU(I,J,M,N)=DUMMY     ! to be used later by GET_THEFLU
            ENDDO
            ENDDO
            ENDDO
            ENDDO

            CALL ESHELBY (AXB,C4GA,0.,E4GA,R4GA,AUX33,AUX33,
     #                   PDIL,AUX3333,AUX3333,IOPTION)

C *** ROTATES THE ESHELBY TENSOR FOR PHASE OR GRAIN BACK INTO SAMPLE AXES.
            DO I=1,3
            DO J=1,3
            DO M=1,3
            DO N=1,3
              DUMMY1=0.
              DUMMY2=0.
              DO I1=1,3
              DO J1=1,3
              DO M1=1,3
              DO N1=1,3
                DUMMY1=DUMMY1+EIGB(I,I1)*EIGB(J,J1)*EIGB(M,M1)
     #                *EIGB(N,N1)*E4GA(I1,J1,M1,N1)
                DUMMY2=DUMMY2+EIGB(I,I1)*EIGB(J,J1)*EIGB(M,M1)
     #                *EIGB(N,N1)*R4GA(I1,J1,M1,N1)
              ENDDO
              ENDDO
              ENDDO
              ENDDO
              E4SA(I,J,M,N)=DUMMY1
              R4SA(I,J,M,N)=DUMMY2
            ENDDO
            ENDDO
            ENDDO
            ENDDO

c      call voigt (aux6,aux33,aux66,r4sa,4)
c      write(10,'(''   R4SA in Voigt notation'')')
c      write(10,'(6f12.4)') aux66
	    
	  RETURN 
	  END
C
C ***********************************************************************
C     SUBROUTINE ESHELBY      --->      VERSION 15/NOV/07
C
C     IOPTION=0: Initialize arrays assoc. with Gauss integration points.
C     IOPTION=1: Calculate elastic Eshelby tensor for elastic inclusion.
C     IOPTION=2: Calculate incompressible Eshelby tensors ESIM (strain-
C                rate) & ESCR (spin-rate) associated with the visco-
C                plastic inclusion.
C     IOPTION=3: Calculate incompressible and hydrostatic Eshelby tensors
C                PESH (deviatoric pressure), PDIL (spherical pressure) &
C                DESH (dilatation) for visco-plastic inclusion.
C     IOPTION=4: Calculates 1st term in d(S)/d(M)=d(T)/d(M)*L+T*d(L)/d(M)
C     IOPTION=5: Calculates 2nd term in d(S)/d(M)=d(T)/d(M)*L+T*d(L)/d(M)
C
C     Options 2-3-4-5 admit a non-zero visco-plastic bulk modulus KEFF.
C
C     Algorithms are based in Lebensohn et al, MSMSE 6 (1998) p.447.
C
C     Uses explicit matrix inversion and explicit Voigt notation (when
C     possible) to optimize computing time.
C
C     Modified oct/2005 to adapt number of integration points to the shape
C     of the ellipsoid in order to keep Eshelby tensor within a certain
C     tolerance (based on analysis done by Gwenaelle Proust).
C     Aspect ratio criterion was adopted for the case AXIS(2)>AXIS(1)>AXIS(3)
C
C     Modified oct/2007 to use a Gauss-Lobatto integration with fix (ten)
C     integration points and weights, chosen to optimize the integration
C     within (eleven) different domains of the ellipsoid aspect ratios.
C     (based on analysis made by Laurent Capolungo).
C        if IGAUSSLEG=0 uses Gauss Lobatto  (cases 1 to 11) (hardwired)
C        if IGAUSSLEG=1 uses Gauss Legendre (case 12) (kept as an option)
C ***********************************************************************

      SUBROUTINE ESHELBY (axis,c4,keff,esim,escr,
     #                    desh,pesh,pdil,dldm,dsddm,ioption)
      USE ESH123

      DIMENSION c4(3,3,3,3),esim(3,3,3,3),escr(3,3,3,3)
      DIMENSION p(3,3),pesh(3,3),desh(3,3)
      DIMENSION c2(6,6),gamma2(6,6),gamma4(3,3,3,3)
      DIMENSION axis(3),x1(10),a1(10),a1inv(10)
      DIMENSION aa1x(6),aa2x(3,3),aaww1x(6),aaww2x(3,3)

      dimension dldm(3,3,3,3),dsddm(3,3,3,3),dldm2(6,6)
      dimension da1(6),da(4,4),dainv(4,4),ainv(4,4)
      dimension aux33(3,3),aux66(6,6),aux3333(3,3,3,3),aux44(4,4)

      DIMENSION xph(ngaumx),xth(ngaumx),wph(ngaumx),wth(ngaumx)

      REAL      keff
      INTEGER   case

      pi=4.*atan(1.0)

      IGAUSSLEG=0     ! hardwires the Gauss-Lobatto option
C     IGAUSSLEG=1     ! hardwires the Gauss-Legendre option
      NGAU_c12 =40    ! make sure that ngau_c12 for Gauss-Leg is .le. ngaumx

C ***********************************************************************
C     INITIALIZATION RUN
C     Calculates Gauss-Legendre integration points and weights in the
c     interval [0,pi]. Initializes Gauss-Lobatto points and weights.
C     Initializes arrays associated with each point to avoid repeating
C     its calculation at every call. All values of the points and
C     weights were calculated to optimize the error.
C ***********************************************************************

      if(ioption.eq.0) then

        do i=1,11
          ngaussph(i)=10        ! Gauss-Lobatto
          ngaussth(i)=10
        end do
        ngaussph(12)=ngau_c12   ! Gauss-Legendre
        ngaussth(12)=ngau_c12

c****************************
c  CASE 1
c****************************
        punti(1,1)=4.71236594E-02
        punti(2,1)=0.241774723
        punti(3,1)=0.565131843
        punti(4,1)=0.968887568
        punti(5,1)=1.37937832
        punti(6,1)=1.76221442
        punti(7,1)=2.17270517
        punti(8,1)=2.57646084
        punti(9,1)=2.89981818
        punti(10,1)=3.09446883

        pesi(1,1)=0.120191820
        pesi(2,1)=0.264987558
        pesi(3,1)=0.373805553
        pesi(4,1)=0.420841277
        pesi(5,1)=0.390970200
        pesi(6,1)=0.390970260
        pesi(7,1)=0.420841366
        pesi(8,1)=0.373805553
        pesi(9,1)=0.264987499
        pesi(10,1)=0.120192111

c****************************
c  CASE 2
c****************************
        punti(1,2)=1.57080423E-02
        punti(2,2)=0.144995824
        punti(3,2)=0.425559640
        punti(4,2)=0.829968274
        punti(5,2)=1.31460333
        punti(6,2)=1.82698941
        punti(7,2)=2.31162453
        punti(8,2)=2.71603298
        punti(9,2)=2.99659705
        punti(10,2)=3.12588477

        pesi(1,2)=5.41692823E-02
        pesi(2,2)=0.207461149
        pesi(3,2)=0.348739326
        pesi(4,2)=0.452716887
        pesi(5,2)=0.507709801
        pesi(6,2)=0.507709682
        pesi(7,2)=0.452716798
        pesi(8,2)=0.348738998
        pesi(9,2)=0.207461327
        pesi(10,2)=5.41692935E-02

c****************************
c  CASE 3
c****************************
        punti(1,3)=3.76990959E-02
        punti(2,3)=0.198626831
        punti(3,3)=0.483041346
        punti(4,3)=0.871647120
        punti(5,3)=1.32964790
        punti(6,3)=1.81194484
        punti(7,3)=2.26994562
        punti(8,3)=2.65855122
        punti(9,3)=2.94296598
        punti(10,3)=3.10389376

        pesi(1,3)=9.68142375E-02
        pesi(2,3)=0.224478707
        pesi(3,3)=0.341134071
        pesi(4,3)=0.430180043
        pesi(5,3)=0.478189558
        pesi(6,3)=0.478189170
        pesi(7,3)=0.430180043
        pesi(8,3)=0.341134191
        pesi(9,3)=0.224478647
        pesi(10,3)=9.68143344E-02

c****************************
c  CASE 4
c****************************
        punti(1,4)=3.45576368E-02
        punti(2,4)=0.187556863
        punti(3,4)=0.468425453
        punti(4,4)=0.859980166
        punti(5,4)=1.32527423
        punti(6,4)=1.81631863
        punti(7,4)=2.28161263
        punti(8,4)=2.67316723
        punti(9,4)=2.95403576
        punti(10,4)=3.10703516

        pesi(1,4)=8.95763785E-02
        pesi(2,4)=0.217725381
        pesi(3,4)=0.341026783
        pesi(4,4)=0.435772508
        pesi(5,4)=0.486694932
        pesi(6,4)=0.486695170
        pesi(7,4)=0.435772508
        pesi(8,4)=0.341026902
        pesi(9,4)=0.217725128
        pesi(10,4)=8.95764604E-02

c****************************
c  CASE 5
c****************************
        punti(1,5)= 3.14158052E-02
        punti(2,5)=0.177928671
        punti(3,5)= 0.457155794
        punti(4,5)= 0.851592362
        punti(5,5)= 1.32222414
        punti(6,5)= 1.81936860
        punti(7,5)=2.29000044
        punti(8,5)=2.68443704
        punti(9,5)=2.96366405
        punti(10,5)=3.11017680

        pesi(1,5)=8.26927349E-02
        pesi(2,5)=0.213228315
        pesi(3,5)=0.342008322
        pesi(4,5)=0.440196186
        pesi(5,5)=0.492670894
        pesi(6,5)=0.492670983
        pesi(7,5)=0.440195888
        pesi(8,5)=0.342008322
        pesi(9,5)=0.213227972
        pesi(10,5)=8.26930404E-02

c****************************
c  CASE 6
c****************************
        punti(1,6)= 2.98452154E-02
        punti(2,6)=0.173592165
        punti(3,6)=0.452448040
        punti(4,6)=0.848216832
        punti(5,6)=1.32101476
        punti(6,6)=1.82057810
        punti(7,6)= 2.29337597
        punti(8,6)=2.68914461
        punti(9,6)=2.96800065
        punti(10,6)=3.11174774

        pesi(1,6)=7.93928578E-02
        pesi(2,6)=0.211627841
        pesi(3,6)=0.342669785
        pesi(4,6)=0.442057431
        pesi(5,6)=0.495048553
        pesi(6,6)=0.495048642
        pesi(7,6)=0.442057490
        pesi(8,6)=0.342670023
        pesi(9,6)=0.211627468
        pesi(10,6)=7.93929026E-02

c****************************
c  CASE 7
c****************************
        punti(1,7)=2.67036632E-02
        punti(2,7)=0.165752888
        punti(3,7)=0.444431901
        punti(4,7)=0.842614472
        punti(5,7)=1.31902647
        punti(6,7)= 1.82256627
        punti(7,7)=2.29897833
        punti(8,7)=2.69716072
        punti(9,7)=2.97583985
        punti(10,7)=3.11488938

        pesi(1,7)=7.30879456E-02
        pesi(2,7)=0.209402516
        pesi(3,7)=0.344104946
        pesi(4,7)=0.445234656
        pesi(5,7)=0.498966068
        pesi(6,7)= 0.498966306
        pesi(7,7)=0.445234746
        pesi(8,7)= 0.344104946
        pesi(9,7)=0.209402665
        pesi(10,7)=7.30878562E-02

c****************************
c  CASE 8
c****************************
        punti(1,8)=2.67036632E-02
        punti(2,8)=0.165752888
        punti(3,8)=0.444431901
        punti(4,8)=0.842614472
        punti(5,8)=1.31902647
        punti(6,8)=1.82256627
        punti(7,8)=2.29897833
        punti(8,8)=2.69716072
        punti(9,8)=2.97583985
        punti(10,8)=3.11488938

        pesi(1,8)=7.30879456E-02
        pesi(2,8)=0.209402516
        pesi(3,8)=0.344104946
        pesi(4,8)=0.445234656
        pesi(5,8)=0.498966068
        pesi(6,8)=0.498966306
        pesi(7,8)=0.445234746
        pesi(8,8)=0.344104946
        pesi(9,8)=0.209402665
        pesi(10,8)= 7.30878562E-02

c****************************
c  CASE 9
c****************************
        punti(1,9)=2.43473575E-02
        punti(2,9)=0.160516247
        punti(3,9)=0.439386278
        punti(4,9)=0.839168847
        punti(5,9)=1.31781363
        punti(6,9)=1.82377899
        punti(7,9)=2.30242372
        punti(8,9)=2.70220637
        punti(9,9)=2.98107672
        punti(10,9)=3.11724544

        pesi(1,9)=6.86219111E-02
        pesi(2,9)=0.208388865
        pesi(3,9)=0.345189095
        pesi(4,9)=0.447236270
        pesi(5,9)=0.501360059
        pesi(6,9)=0.501359940
        pesi(7,9)=0.447236151
        pesi(8,9)=0.345189214
        pesi(9,9)=0.208388969
        pesi(10,9)=6.86219335E-02

c****************************
c  CASE 10
c****************************
        punti(1,10)=2.19910536E-02
        punti(2,10)=0.155757755
        punti(3,10)=0.434985727
        punti(4,10)=0.836206555
        punti(5,10)=1.31677616
        punti(6,10)= 1.82481658
        punti(7,10)=2.30538607
        punti(8,10)=2.70660710
        punti(9,10)=2.98583508
        punti(10,10)=3.11960149

        pesi(1,10)=6.43825606E-02
        pesi(2,10)=0.207786217
        pesi(3,10)=0.346235514
        pesi(4,10)=0.448981822
        pesi(5,10)=0.503410578
        pesi(6,10)= 0.503410578
        pesi(7,10)=0.448981792
        pesi(8,10)=0.346235693
        pesi(9,10)=0.207785636
        pesi(10,10)= 6.43827692E-02

c****************************
c  CASE 11
c****************************
        punti(1,11)=2.04204638E-02
        punti(2,11)=0.152822554
        punti(3,11)=0.432348520
        punti(4,11)=0.834448099
        punti(5,11)=1.31616223
        punti(6,11)=1.82543063
        punti(7,11)=2.30714464
        punti(8,11)=2.70924401
        punti(9,11)=2.98877001
        punti(10,11)=3.12117243

        pesi(1,11)=6.16818815E-02
        pesi(2,11)=0.207559645
        pesi(3,11)=0.346902698
        pesi(4,11)=0.450027168
        pesi(5,11)=0.504624724
        pesi(6,11)= 0.504624426
        pesi(7,11)=0.450027317
        pesi(8,11)=0.346902847
        pesi(9,11)=0.207559645
        pesi(10,11)=6.16819337E-02

C ******************************************************************
C  CASE 12: GAULEG generates the integration points
C ******************************************************************

        CALL ESH_GAUSS_LEGENDRE (0.0,pi,puntigl,pesigl,ngaussph(12))

C ******************************************************************
C *** Calculates and saves arrays that depend on integration points

       do case=1,12

         if (case.eq.12) then
           do i=1,ngau_c12
             xph(i)=puntigl(i)
             xth(i)=puntigl(i)
             wph(i)= pesigl(i)
             wth(i)= pesigl(i)
           end do
         else
           do i=1,10
             xph(i)=punti(i,case)
             xth(i)=punti(i,case)
             wph(i)= pesi(i,case)
             wth(i)= pesi(i,case)
           end do
         end if
c--------------------------------------------------------------
c *** integration [0,pi][0,pi] adds a factor 2 in Eqs. B11 & B14.

         do ith=1,ngaussth(case)
           sinth=sin(xth(ith))
           costh=cos(xth(ith))
           simbtet=wth(ith)*sinth/(2.0*pi)

           do iph=1,ngaussph(case)
             ny=iph+(ith-1)*ngaussph(case)
             ww(case,ny)=simbtet*wph(iph)
             alpha(case,1,ny)=sinth*cos(xph(iph))
             alpha(case,2,ny)=sinth*sin(xph(iph))
             alpha(case,3,ny)=costh

             do i=1,3
             do j=1,3
               aa2x(i,j)  =alpha(case,i,ny)*alpha(case,j,ny)
               aaww2x(i,j)=aa2x(i,j)*ww(case,ny)
             enddo
             enddo
             call voigt(aa1x  ,aa2x  ,c2,c4,2)
             call voigt(aaww1x,aaww2x,c2,c4,2)
             do i=1,6
               aa1(case,i,ny)  =aa1x(i)
               aaww1(case,i,ny)=aaww1x(i)
             enddo

c *** Array AWW is used only if ICAUCHY=1.
c            do i=1,3
c              aww(case,i,ny)=alpha(case,i,ny)*ww(case,ny)
c            enddo

           enddo
         enddo

        enddo      ! end of do case=1,12

      endif      ! ENDIF FOR IOPTION=0

C************************************************************************
C     End of initialization
C************************************************************************

C ***********************************************************************
C     CALCULATION OF ESHELBY TENSORS FOR STIFFNESS 'C4' AND ELLIPSOID
C     AXES 'AXIS'
C     ASSUMES: AXIS2 > AXIS1 > AXIS3  --> RATIO1 > RATIO2 > 1
C ***********************************************************************

      if(ioption.ge.1) then

        abc=axis(1)*axis(2)*axis(3)
        ratio1=axis(2)/axis(3)
        ratio2=axis(1)/axis(3)

        if (igaussleg.eq.1) then
          case=12
        else
          dte(0) = 0.0
          dte(1) =-0.7*ratio1+7
          dte(2) =-ratio1+17
          dte(3) =-ratio1+23
          dte(4) =-ratio1+26
          dte(5) =-ratio1+29.3
          dte(6) =-ratio1+32
          dte(7) =-ratio1+34.85
          dte(8) =-ratio1+37
          dte(9) =-ratio1+41.9
          dte(10)=-ratio1+44.5
          case=11
          do i=1,10
            if(ratio2.ge.dte(i-1) .and. ratio2.lt.dte(i)) case=i
          enddo
        endif

c--------------------------------------------------------

        npoints=ngaussph(case)*ngaussth(case)

        pdil=0.
        do j=1,3
        do i=1,3
          p(i,j)=0.
        enddo
        enddo
        do j=1,6
        do i=1,6
          gamma2(i,j)=0.
        enddo
        enddo

        call voigt(aa1x,aa2x,c2,c4,4)
        IF(IOPTION.EQ.5) CALL VOIGT(AA1X,AA2X,DLDM2,DLDM,4)

        do ny=1,npoints

c   Compute Voigt components A1(1)-A(6) of tensor A(3,3) defined by Eq.B3:
c   --->  A(i,j)=L(i,j,k,l)*a(j)*a(l)

          do i=1,6
            aa1x(i)=aa1(case,i,ny)
          enddo
          call esh_mult_voigt(c2,aa1x,a1(1:6))   ! c2*aa1x=a1

      IF(IOPTION.EQ.1) THEN

c   If solving an elastic inclusion invert the system
c   --> A(3,3) x X(3,3) = C(3,3)
c   Inverts A(3,3) using explicit Voigt notation.
c   Uses explicit form of C(3,3) to calculate solution in Voigt notation.

            call esh_inv3_voigt(a1,a1inv)
            do i=1,6
              x1(i)=a1inv(i)
            enddo

      ELSE IF(IOPTION.GE.2) THEN

c   If solving a visco-plastic inclusion defines components A1(7) to A1(10).
c   Solves the system given by Eq.B4 --> A(4,4) x X(4,4) = C(4,4)
c   Inverts A(4,4) using explicit Voigt notation.
c   Uses explicit form of C(4,4) to calculate solution in Voigt notation.
c   The solution is symmetric. Numerical deviation from symmetry is averaged.

            a1(7) = alpha(case,1,ny)
            a1(8) = alpha(case,2,ny)
            a1(9) = alpha(case,3,ny)
            a1(10)= 0.
            if(keff.gt.0.) a1(10)=-1./keff

            call esh_inv4_voigt(a1,a1inv)

            do i=1,10
              x1(i)=a1inv(i)
            enddo

      ENDIF

      IF(IOPTION.EQ.5) THEN

        CALL SO_VOIGT10(A1INV,AINV,1)

        call esh_mult_voigt(dldm2,aa1x,da1)        ! dldm2*aa1x=da1
        call voigt(da1,aux33,aux66,aux3333,1)

        do i=1,3
        do j=1,3
          da(i,j)=aux33(i,j)
        enddo
        enddo
c
        do i=1,3
        da(i,4)=0.
        da(4,i)=0.
        enddo
        da(4,4)=0.
c
c           dAinv/dp = - Ainv  : dA/dp : Ainv
c                x1  = - a1inv :  da1  : a1inv
c
        do i=1,4
        do j=1,4
        dummy=0.
        do k=1,4
        dummy=dummy+da(i,k)*ainv(k,j)
        enddo
        aux44(i,j)=dummy
        enddo
        enddo

        do i=1,4
        do j=1,4
        dummy=0.
        do k=1,4
          dummy=dummy+ainv(i,k)*aux44(k,j)
        enddo
        dainv(i,j)=-dummy
        enddo
        enddo

        CALL SO_VOIGT10(X1,DAINV,2)

      ENDIF

          ro3=((alpha(case,1,ny)*axis(1))**2+
     #         (alpha(case,2,ny)*axis(2))**2+
     #         (alpha(case,3,ny)*axis(3))**2)**1.5
          abcoro3=abc/ro3

c   Compute the Eshelby integral Eq.B11 defining:
c         Gamma(m,j,n,i)=T(m,n,i,j)=a(m)*a(j)*G(n,i)
c   with the property:
c         Gamma(m,j,n,i)=Gamma(j,m,n,i)=Gamma(m,j,i,n)

          do i=1,6
          do j=1,6
            gamma2(i,j)=gamma2(i,j)+aaww1(case,i,ny)*x1(j)*abcoro3
          enddo
          enddo

c   Compute the pressure related Eshelby integral Eq.B14
          if(ioption.eq.3) then
            do j=1,3
            do i=1,3
              p(i,j)=p(i,j)+aww(case,j,ny)*x1(i+6)*abcoro3
            enddo
            enddo
            pdil=pdil+ww(case,ny)*x1(10)*abcoro3
          endif

        end do   ! end of loop over double integration

c ********************************************************************
c   Go back to the 3*3*3*3 notation
        call voigt(aa1x,aa2x,gamma2,gamma4,3)

c   Compute symmetric (distortion) Eshelby tensor from Eq.B9.
c       esim(n,m,k,l)=0.5*(gamma(m,j,n,i)+gamma(n,j,m,i))*c4(i,j,k,l)
c   Compute anti-symmetric (rotation) Eshelby tensor from Eq.B9.
c       escr(n,m,k,l)=0.5*(gamma(m,j,n,i)-gamma(n,j,m,i))*c4(i,j,k,l)

        do l=1,3
        do k=1,3
        do m=1,3
        do n=1,3
          dumsim=0.
          dumscr=0.
        do j=1,3
        do i=1,3

        IF(IOPTION.NE.4) THEN
          dumsim=dumsim+(gamma4(m,j,n,i)+gamma4(n,j,m,i))*c4(i,j,k,l)
          dumscr=dumscr+(gamma4(m,j,n,i)-gamma4(n,j,m,i))*c4(i,j,k,l)
        ELSE
          dumsim=dumsim+(gamma4(m,j,n,i)+gamma4(n,j,m,i))*dldm(i,j,k,l)
          dumscr=dumscr+(gamma4(m,j,n,i)-gamma4(n,j,m,i))*dldm(i,j,k,l)
        ENDIF

        enddo
        enddo

        IF(IOPTION.LT.4) THEN
          esim(n,m,k,l)=0.5*dumsim
          escr(n,m,k,l)=0.5*dumscr
        ELSE
          dsddm(n,m,k,l)=0.5*dumsim
        ENDIF

        enddo
        enddo
        enddo
        enddo

c   Compute pressure & dilatation related Eshelby tensors (Eq.B13)

        IF(IOPTION.EQ.3) THEN
          do l=1,3
          do k=1,3
            pesh(k,l)=0.
            do j=1,3
            do i=1,3
              pesh(k,l)=pesh(k,l)+p(i,j)*c4(i,j,k,l)
            end do
            end do
          end do
          end do
          do j=1,3
          do i=1,3
            desh(i,j)=(p(i,j)+p(j,i))/2.
          end do
          end do
        ENDIF

      endif      !  endif for IOPTION.GE.1

      RETURN
      END

C ********************************************************************
      SUBROUTINE ESH_GAUSS_LEGENDRE (x1,x2,x,w,n)

      dimension x(n),w(n)
      parameter(eps=1.e-07)
      pi=4.0*atan(1.0)
	  
      m=(n+1)/2
      xm=0.5*(x1+x2)
      xl=0.5*(x2-x1)
      xn=n
      do i=1,m
      xi=i
      z=cos(pi*(xi-.25)/(xn+0.5))
c
      iter=0
1     continue
      iter=iter+1
c
c     R.L. 8/2/97
c
      if(iter.gt.10000) then
      write(*,*)'GAULEG WARNING: TOL 1.e-07 NEVER REACHED - ERR = ',
     #           abs(z-z1)
      return
      endif
c
      p1=1.0
      p2=0.0
      do j=1,n
      xj=j
      p3=p2
      p2=p1
      p1=((2.0*j-1.0)*z*p2-(xj-1.0)*p3)/xj
      enddo
      pp=n*(z*p1-p2)/(z*z-1.0)
      z1=z
      z=z1-p1/pp

      if(abs(z-z1).gt.eps) go to 1
      x(i)=xm-xl*z
      x(n+1-i)=xm+xl*z
      w(i)=2.0*xl/((1.0-z*z)*pp*pp)
      w(n+1-i)=w(i)
      enddo
      return
      end
C
C *************************************************************************
C     SUBROUTINE ESH_INV3_VOIGT   --->   version 23/jul/01
C
C     Inverts the 3x3 symmetric matrix 'A' using explicit Voigt notation:
C     11->1, 22->2, 33->3, 23=32->4, 31=13->5, 12=21->6
C *************************************************************************

      SUBROUTINE ESH_INV3_VOIGT (A,AINV)

      DIMENSION A(10),AINV(10)

      DET = A(1)*A(2)*A(3) + 2*A(4)*A(5)*A(6) - A(1)*A(4)*A(4)
     #     - A(2)*A(5)*A(5) - A(3)*A(6)*A(6)

      AINV(1) = ( A(2)*A(3) - A(4)*A(4))/DET
      AINV(2) = ( A(1)*A(3) - A(5)*A(5))/DET
      AINV(3) = ( A(1)*A(2) - A(6)*A(6))/DET
      AINV(4) = (-A(1)*A(4) + A(5)*A(6))/DET
      AINV(5) = ( A(4)*A(6) - A(2)*A(5))/DET
      AINV(6) = (-A(3)*A(6) + A(4)*A(5))/DET

      RETURN
      END
C
C **********************************************************************
C     SUBROUTINE ESH_INV4_VOIGT   --->   VERSION 20/JUL/01

C     Inverts the 4*4 symmetric matrix 'A' using explicit Voigt notation:
C     11-->1, 22-->2, 33-->3, 23=32-->4, 31=13-->5, 12=21-->6
C     14-->7, 24-->8, 34-->9, 44-->10.
C **********************************************************************

      SUBROUTINE ESH_INV4_VOIGT (A,AINV)

      DIMENSION A(10),AINV(10)


      ainv(1) = a(2)*a(3)*a(10)+2*a(4)*a(8)*a(9) -
     #          a(2)*a(9)*a(9)-a(3)*a(8)*a(8)-a(4)*a(4)*a(10)

      ainv(2) = a(1)*a(3)*a(10)+2*a(5)*a(7)*a(9) -
     #          a(1)*a(9)*a(9)-a(3)*a(7)*a(7)-a(5)*a(5)*a(10)

      ainv(3) = a(1)*a(2)*a(10)+2*a(6)*a(7)*a(8) -
     #          a(1)*a(8)*a(8)-a(2)*a(7)*a(7)-a(6)*a(6)*a(10)

      ainv(4) = a(1)*a(4)*a(10)+a(5)*a(7)*a(8)+a(6)*a(7)*a(9) -
     #          a(1)*a(8)*a(9)-a(4)*a(7)*a(7)-a(5)*a(6)*a(10)
      ainv(4) =-ainv(4)

      ainv(5) = a(4)*a(6)*a(10)+a(2)*a(7)*a(9)+a(5)*a(8)*a(8) -
     #          a(4)*a(7)*a(8)-a(6)*a(8)*a(9)-a(2)*a(5)*a(10)

      ainv(6) = a(3)*a(6)*a(10)+a(5)*a(8)*a(9)+a(4)*a(7)*a(9) -
     #          a(3)*a(7)*a(8)-a(6)*a(9)*a(9)-a(4)*a(5)*a(10)
      ainv(6) =-ainv(6)

      ainv(7) = a(4)*a(6)*a(9)+a(4)*a(5)*a(8)+a(2)*a(3)*a(7) -
     #          a(4)*a(4)*a(7)-a(2)*a(5)*a(9)-a(3)*a(6)*a(8)
      ainv(7) =-ainv(7)

      ainv(8) = a(1)*a(4)*a(9)+a(5)*a(5)*a(8)+a(3)*a(6)*a(7) -
     #          a(4)*a(5)*a(7)-a(5)*a(6)*a(9)-a(1)*a(3)*a(8)

      ainv(9) = a(1)*a(2)*a(9)+a(5)*a(6)*a(8)+a(4)*a(6)*a(7) -
     #          a(2)*a(5)*a(7)-a(6)*a(6)*a(9)-a(1)*a(4)*a(8)
      ainv(9) =-ainv(9)

      ainv(10)=a(1)*a(2)*a(3)+2*a(4)*a(5)*a(6) -
     #         a(1)*a(4)*a(4)-a(2)*a(5)*a(5)-a(3)*a(6)*a(6)

      det=   a(1)*ainv(1)+   a(2)*ainv(2)+   a(3)*ainv(3)+
     #    2.*a(4)*ainv(4)+2.*a(5)*ainv(5)+2.*a(6)*ainv(6)+
     #    2.*a(7)*ainv(7)+2.*a(8)*ainv(8)+2.*a(9)*ainv(9)+
     #       a(10)*ainv(10)
      det=   det/4.

      do i=1,10
        ainv(i)=ainv(i)/det
      enddo

      return
      end
C
C ***********************************************************
C     SUBROUTINE ESH_MULT_VOIGT   -->   version 16/MAY/2012
C ***********************************************************
      SUBROUTINE ESH_MULT_VOIGT(B,C,A)

C     Performs the multiplication:
C        A(i,k)=B(i,j,k,l)*C(j,l) using Voigt's notation
C        B is a 6*6 symmetric matrix
C        C is a 3*3 symmetric tensor
C        A will be a 3*3 symmetric tensor

      DIMENSION B(6,6),C(6),A(6)

      A(1)=B(1,1)*C(1)+B(6,6)*C(2)+B(5,5)*C(3)
     #    +2*(B(5,6)*C(4)+B(1,5)*C(5)+B(1,6)*C(6))

      A(2)=B(6,6)*C(1)+B(2,2)*C(2)+B(4,4)*C(3)
     #    +2*(B(2,4)*C(4)+B(4,6)*C(5)+B(2,6)*C(6))

      A(3)=B(5,5)*C(1)+B(4,4)*C(2)+B(3,3)*C(3)
     #    +2*(B(3,4)*C(4)+B(3,5)*C(5)+B(4,5)*C(6))

      A(4)=B(5,6)*C(1)+B(2,4)*C(2)+B(3,4)*C(3)
     #      +(B(2,3)+B(4,4))*C(4)
     #      +(B(3,6)+B(4,5))*C(5)
     #      +(B(4,6)+B(2,5))*C(6)

      A(5)=B(1,5)*C(1)+B(4,6)*C(2)+B(3,5)*C(3)
     #      +(B(3,6)+B(4,5))*C(4)
     #      +(B(1,3)+B(5,5))*C(5)
     #      +(B(1,4)+B(5,6))*C(6)

      A(6)=B(1,6)*C(1)+B(2,6)*C(2)+B(4,5)*C(3)
     #      +(B(4,6)+B(2,5))*C(4)
     #      +(B(1,4)+B(5,6))*C(5)
     #      +(B(1,2)+B(6,6))*C(6)

      RETURN
      END
C
C *****************************************************************************
C     subroutine euler      -->   version 2022 02 09   ************************
C *****************************************************************************
C     CALCULATES THE EULER ANGLES ASSOCIATED WITH THE TRANSFORMATION
C     MATRIX A(I,J) IF IOPT=1 AND VICEVERSA IF IOPT=2
C *** A(i,j) TRANSFORMS FROM SAMPLE AXES TO CRYSTAL AXES BUT
C     AG(i,j) INSIDE VPSC TRANSFORMS FROM CRYSTAL TO SAMPLE AXES !!
C *** ph,th,om ARE THE BUNGE EULER ANGLES (in degrees) OF ca REFERRED TO sa:
C       * ph : FIRST ROTATION AROUND sa X3
C       * th : SECOND ROTATION AROUND ca X1
C       * om : THIRD ROTATION AROUND ca X3
C *****************************************************************************
      subroutine euler (iopt,ph,th,om,a)

      dimension a(3,3)
      pi=4.*atan(1.d0)

      if(iopt.eq.1) then
        if(abs(a(3,3)).gt.1.0001) then
          write(*,'('' --> a(3,3)>1.0001 in euler'',e12.4)') a(3,3)
          stop
        endif
        if(a(3,3) .ge. 0.9999) then
		  a(3,3)=1.0
          th=acos(a(3,3))
          om=0. 
          ph=atan2(a(1,2),a(1,1))
        else if(a(3,3) .le. -0.9999) then
		  a(3,3)=-1.0
          th=acos(a(3,3))
          om=0. 
          ph=atan2(a(1,2),a(1,1)) 
        else
          th=acos(a(3,3))
          sth=sin(th)
          om=atan2(a(1,3)/sth,a(2,3)/sth)
          ph=atan2(a(3,1)/sth,-a(3,2)/sth)
        endif
        th=th*180./pi
        ph=ph*180./pi
        om=om*180./pi

      else if(iopt.eq.2) then
        sph=sin(ph*pi/180.)
        cph=cos(ph*pi/180.)
        sth=sin(th*pi/180.)
        cth=cos(th*pi/180.)
        som=sin(om*pi/180.)
        com=cos(om*pi/180.)
        a(1,1)=com*cph-sph*som*cth
        a(2,1)=-som*cph-sph*com*cth
        a(3,1)=sph*sth
        a(1,2)=com*sph+cph*som*cth
        a(2,2)=-sph*som+cph*com*cth
        a(3,2)=-sth*cph
        a(1,3)=sth*som
        a(2,3)=com*sth
        a(3,3)=cth
      endif

      return
      end

C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     FUNCTIONS OF VPSC
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c **********************************************************************
      function det(a)
	  
      dimension a(3,3)

      det=a(1,1)*a(2,2)*a(3,3)
     #   +a(1,2)*a(2,3)*a(3,1)
     #   +a(1,3)*a(2,1)*a(3,2)
     #   -a(1,3)*a(2,2)*a(3,1)
     #   -a(2,3)*a(3,2)*a(1,1)
     #   -a(1,2)*a(2,1)*a(3,3)
      return
      end

C ***********************************************************************
      FUNCTION random2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL random2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *  IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *  NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy

      idum2=123456789
      iv=0
      iy=0
      
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
        enddo
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      random2=min(AM*iy,RNMX)
      return
      END

C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     FUNCTION TMISMATCH   ---->   VERSION OF 13/SEP/2019
C
C     CALCULATES RELATIVE DIFFERENCE BETWEEN TWO 5x5 MATRICES.
C     THE DIFFERENCE IS RELATIVE TO THE NORM OF THE ARITHMETIC AVERAGE.
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      FUNCTION TMISMATCH (T1,T2,N)

      DIMENSION T1(N*N),T2(N*N)
      DIMENSION T_dif(N*N),T_ave(N*N)

      DO I=1,N*N
        T_dif(I) =T1(I)-T2(I)
        T_ave(I)=0.5*(T1(I)+T2(i))
      ENDDO
      TMISMATCH=VNORM(T_dif,N*N)/VNORM(T_ave,N*N)

      RETURN
      END

C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     FUNCTION TNORM    ---->   VERSION OF 05/sep/2019
C
C     CALCULATES THE NORM OF A NxN MATRIX
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      FUNCTION TNORM (T,N)

      DIMENSION T(N*N)

      TNORM=0.
      DO I=1,N*N
        TNORM=TNORM+T(I)*T(I)
      ENDDO
      TNORM=SQRT(TNORM)

      RETURN
      END

C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     FUNCTION VMISMATCH   ---->   VERSION OF 23/NOV/2011
C
C     CALCULATES RELATIVE DIFFERENCE BETWEEN TWO VECTORS WITH N COMPONENTS.
C     THE DIFFERENCE IS RELATIVE TO THE NORM OF THE ARITHMETIC AVERAGE.
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      FUNCTION VMISMATCH (V1,V2,N)

      DIMENSION V1(N),V2(N)
      DIMENSION V_dif(N),V_ave(N)

      DO I=1,N
        V_dif(I)=V1(I)-V2(I)
        V_ave(I)=0.5*(V1(I)+V2(i))
      ENDDO
      VMISMATCH =VNORM(V_dif,N)/VNORM(V_ave,N)

      RETURN
      END

C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     FUNCTION VNORM    ---->   VERSION OF 23/NOV/2011
C
C     CALCULATES THE NORM OF A VECTOR V WITH N COMPONENTS
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      FUNCTION VNORM (V,N)

      DIMENSION V(N)

      VNORM=0.
      DO I=1,N
        VNORM=VNORM+V(I)*V(I)
      ENDDO
      VNORM=SQRT(VNORM)

      RETURN
      END

C ***************************************************************************
C     SUBROUTINE GRAIN_INFO     --->   VERSION 20/APR/2022
C
C     FOR EACH GRAIN CALCULATES ACCUMULATED VON MISES STRAIN, 
C     ACCUMULATED PLASTIC WORK, AND 'TAYLOR' FACTOR DEFINED AS:
C          Mg=(Sg(i).DBAR(i))/TAU(1,g)/|DBAR|
C     CALCULATES PLASTIC WORK RATE USING MACROSCOPIC (wrate1) AND AVERAGE OF
C     GRAIN STRESS & RATE TENSORS (wrate2)
C ***************************************************************************

      SUBROUTINE GRAIN_INFO

      USE VPSC8DIM

      DBARNORM=VNORM(DBAR,5)    ! only deviatoric components
      TAYLAV=0.
      WORKAV=0.

      WRATE1=0.
      DO I=1,5                  ! only deviatoric components
        WRATE1=WRATE1+SAV(I)*DAV(I)
      ENDDO

      WRATE2=0.
      KGX=1
      DO IPH=IPHBOT,IPHTOP
      DO KKK=NGR(IPH-1)+1,NGR(IPH)
        SGVMX=0.
        WRATE=0.
        TAYLX=0.
        DO I=1,5                ! only deviatoric components
          SGVMX=SGVMX+SG(I,KKK)**2
          WRATE=WRATE+SG(I,KKK)*DG(I,KGX)
          TAYLX=TAYLX+SG(I,KKK)*DBAR(I)
        ENDDO
        SGVM(KGX)=SQRT(SGVMX*3./2.)              ! deviatoric Von Mises
        IF(SGVMX.NE.0.) DGVMX=WRATE/SGVM(KGX)    ! defined as the conjugate VM
        IF(SGVMX.EQ.0.) DGVMX=0.

        EPSVMGR(KKK) =EPSVMGR(KKK)+DGVMX *TIME_INCR     ! deviatoric Von Mises
        WORKGR(KKK)  =WORKGR(KKK) +WRATE *TIME_INCR     ! deviatoric non Von Mises
        TAYLGR(KKK)  =TAYLX/CRSS(1,KKK)/DBARNORM        ! deviatoric non Von Mises

C *** THE FOLLOWING ARE ALL non VM. GET WRITTEN IN SUBR WRITE_STRESS_STRAIN		
        TAYLAV=TAYLAV+TAYLGR(KKK)*WGT(KKK)
        WORKAV=WORKAV+WORKGR(KKK)*WGT(KKK)
        WRATE2=WRATE2+WRATE      *WGT(KKK)
		
        KGX=KGX+1
      ENDDO
      ENDDO

      RETURN
      END

C ************************************************************************************
C     SUBROUTINE GRAIN_RATE_AND_MODULI  -->  VERSION 19/OCT/2020
C
C     GIVEN THE STRESS 'X' IN GRAIN 'KKK', CALCULATES STRAIN-RATE AND
C     VISCO-PLASTIC MODULI USING THE RATE SENSITIVITY KINEMATIC LAW
C
C     CALLED FROM SUBR INITIAL_STATE_GUESS 
C     CALLED FROM SUBROUTINE VPSC AFTER INNER LOOP CONVERGENCE
C *************************************************************************************

      SUBROUTINE GRAIN_RATE_AND_MODULI (JSC,KGX,KKK,IPHEL,IPH)

      USE VPSC8DIM
      USE NLNR

      DIMENSION SGX(5),RSS(NSYSMX),RSSX(NSYSMX)

C     COPY MAIN ARRAYS INTO AUXILIAR ARRAYS FOR COMPUTATIONAL EFFICIENCY

      NSYSTX=NSYST(IPHEL)
      DO IS=1,NSYSTX
        ISENSEX(IS)=ISENSE(IS,IPHEL)
        NRSX(IS)   =NRS(IS,IPHEL)
        if(jsc.eq.1.and.irsvar.eq.1) NRSX(IS)=JXRS     ! used when NRS is evolved for convergence
        TAUX(IS)  =CRSS(IS,KKK)
        GAMD0X(IS)=GAMD0S(IS,KGX)
        DO J=1,5
          SCHX(J,IS)=SCH(J,IS,KGX)
        ENDDO
      ENDDO

      DO I=1,5
        SGX(I)=SG(I,KKK)
      ENDDO

C     GETS RESOLVED SHEAR STRESSES 'rss' AND SHEAR RATES 'gamdot'.
C       SIGN(GAMDOT)=SIGN(RSS).
C       NRS CAN BE EVEN OR ODD.
C       RSS IS THE POWER (n-1) AND IS USED TO CALCULATE VISCOUS COMPLIANCE.

      DO IS=1,NSYSTX
        RSSX(IS)=SGX(1)*SCHX(1,IS)+SGX(2)*SCHX(2,IS)+SGX(3)*SCHX(3,IS)+
     #       SGX(4)*SCHX(4,IS)+SGX(5)*SCHX(5,IS)
        IF(.NOT.(RSSX(IS).GT.0 .OR. ISENSEX(IS).EQ.1)) RSSX(IS)=0.
        RSSX(IS)=RSSX(IS)/TAUX(IS)
		
        RSS(IS)  =GAMD0X(IS)* ABS (RSSX(IS)**(NRSX(IS)-1)) /TAUX(IS)  
        GAMDOT(IS,KGX)=GAMD0X(IS)*ABS(RSSX(IS)**NRSX(IS))
     #                                 *SIGN(1.,RSSX(IS))
      ENDDO

C     CALCULATE STRAIN-RATE IN GRAIN 'dg' (INCLUDING TRANSFORMATION RATE)

      DO I=1,5
        DG(I,KGX)=DG_TRANS(I,KKK)   
        DO IS=1,NSYSTX
          DG(I,KGX)=DG(I,KGX)+SCHX(I,IS)*GAMDOT(IS,KGX)
        ENDDO
      ENDDO

C *** CALCULATE CRYSTAL COMPLIANCE Mg_tg AND EPSg_DOT_0
C *** ALGORITMS ARE WRITTEN IN TERMS OF MGTG, WHICH EXCEPT FOR AFFINE & 2nd ORDER,
C *** IS REALLY MGSEC. THE 'n' VALUE 1,n_eff,n IS ACCOUNTED FOR IN THE INTERATION EQ

        DO I=1,5
        DO J=1,5
          MGSEC(I,J,KGX)=0.
          MGTG (I,J,KGX)=0.
          DO IS=1,NSYSTX
            MGSEC(I,J,KGX)=MGSEC(I,J,KGX)+               ! for SX secant equation
     #            SCHX(I,IS)* SCHX(J,IS)*RSS(IS)     
            MGTG (I,J,KGX)=MGTG(I,J,KGX)+ NRSX(IS)*      ! for SX tangent equation
     #            SCHX(I,IS)* SCHX(J,IS)*RSS(IS)
          ENDDO
        ENDDO
        ENDDO

      IF(INTERACTION.EQ.1 .OR. INTERACTION.EQ.5) THEN       ! affine,SO       

        DO I=1,5
          DG0(I,KGX)=0.
        DO J=1,5
          DG0(I,KGX)=DG0(I,KGX)+
     #       (MGSEC(I,J,KGX)-MGTG(I,J,KGX))*SGX(J)          ! (Mg_sec-Mg_tg)*sg
        ENDDO
        ENDDO

      ELSEIF(INTERACTION.EQ.0 .OR. INTERACTION.EQ.2 .OR.    ! fc,sec,neff,tg
     #       INTERACTION.EQ.3 .OR. INTERACTION.EQ.4) THEN 

        DO I=1,5
          DG0(I,KGX)=0.
        DO J=1,5
          MGTG(I,J,KGX)=MGSEC(I,J,KGX)      
        ENDDO
        ENDDO
		
      ENDIF        ! INTERACTION ENDIF

      if (kgx.le.0) then
      write(10,'(i3,''  MGSEC(1,i)'',5e12.4)') 
     #                   kgx,(mgsec(1,i,kgx),i=1,5)
      write(10,'(i3,''  MGTG(1,i) '',5e12.4)') 
     #                   kgx,(mgtg(1,i,kgx),i=1,5)
      endif

      IF (INTERACTION.EQ.5) THEN      ! second order model
        DO IS=1,NSYSTX
          aso(is,kgx)=nrsx(is)* rss(is)       
          eso(is,kgx)=(1-nrsx(is))*gamdot(is,kgx)
        ENDDO
      ENDIF
 
      RETURN
      END
C
C **********************************************************************************
C     SUBROUTINE GRAIN_STRESS  --->   version of 25/AUG/2019
C
C     GIVEN A GUESS STRESS 'SGX(I)' AND GRAIN INDEX 'KKK', SOLVES INTERACTION EQ
C     TO FIND THE STRESS 'SGX' COMPATIBLE WITH MACROSCOPIC COMPLIANCE.
C
C     IF INTERACTION= 0 SOLVES NON-LINEAR POWER LAW FOR A TAYLOR CASE.
C        Dg-Dbar = 0   or, written explicitly: 
C        FNR=sum(mij*(m*sg/crss)^n) +dg_trans -Dbar =0
C     IF INTERACTION> 0 SOLVES INTERACTION EQUATION FOR A SELF-CONSISTENT CASE.
C        Dg-Dbar = -M~ * (Sg-Sbar)   or, written explicitly: 
C        FNR=sum(mij*(m*sg/crss)^n) +dg_trans -Mbar*Sbar-D0bar + M~*Sg - M~*Sbar =0
C ***********************************************************************************

      SUBROUTINE GRAIN_STRESS (KGX,KKK,IPHEL,IPH)

      USE VPSC8DIM
      USE NLNR

      DIMENSION SGX(6),SGXOLD(6)

C     EMPIRIC ALGORITHM TO GET TAULIM FOR NR SUBROUTINE (RL: 01/FEB/00)

      taulim=2.*(Vnorm(dbar,5)/gamd0s(1,1))**0.05
      if(taulim.lt.2.) taulim=2.

      iprint=0   ! controls diagnostic print-out
      ngg=699
      iprx=iprint*(ngg/kkk)*(kkk/ngg)+
     #     iprint*((ngg+ngr(1))/kkk)*(kkk/(ngg+ngr(1)))

      if(iprx.eq.1) then
        write(10,'('' INSIDE GRAIN_STRESS: GRAIN'',i5 )') KKK
      endif

C     COPY MAIN ARRAYS INTO AUXILIAR ARRAYS FOR COMPUTATIONAL EFFICIENCY
C     AND TO MAKE 'NEWTON_RAPHSON' A 'STAND-ALONE' SUBROUTINE

      NSYSTX=NSYST(IPHEL)
      DO IS=1,NSYSTX
        ISENSEX(IS)=ISENSE(IS,IPHEL)
        NRSX(IS)   =NRS(IS,IPHEL)
        if(irsvar.eq.1) NRSX(IS)=JXRS        ! used when NRS is evolved for convergence
        TAUX(IS)  =CRSS(IS,KKK)
        GAMD0X(IS)=GAMD0S(IS,KKK)
        DO J=1,5
          SCHX(J,IS)=SCH(J,IS,KGX)
        ENDDO
      ENDDO

      DO I=1,5
        SGX(I)=SGTRY(I,KGX)
      ENDDO

C *** CORRECTS STRESS 'X' IF IT EXCEEDS THE YIELD SURFACE.

      TAUMAX=0.
      DO IS=1,NSYSTX
        RSSX=SGX(1)*SCHX(1,IS)+SGX(2)*SCHX(2,IS)+SGX(3)*SCHX(3,IS)+
     #            SGX(4)*SCHX(4,IS)+SGX(5)*SCHX(5,IS)
        IF(.NOT.(RSSX.GT.0 .OR. ISENSEX(IS).EQ.1)) RSSX=0.D0
        RSSX=RSSX/TAUX(IS)
        IF(ABS(RSSX).GT.TAUMAX) TAUMAX=ABS(RSSX)
      ENDDO

      if(iprx.eq.1) then
        write(10,'('' GRAIN'',i5,''  TAUMAX='',F10.3)') KKK,TAUMAX
      endif

      IF(TAUMAX. LT. 1.E-10) THEN
        WRITE(*,'('' TAUMAX<1e-10 inside subroutine GRAIN_STRESS'')')
        WRITE(*,'('' GRAIN #'',i5)') KKK
        WRITE(*,'('' SGTRY  ='',5e12.4)') X
        WRITE(*,'('' CRSS  ='',6E12.4)') (TAUX(IS),IS=1,NSYST(IPHEL))
        STOP
      ENDIF

      IF(TAUMAX.GT.TAULIM .OR. INTERACTION.EQ.0) THEN
        DO I=1,5
          SGX(I)=SGX(I)/TAUMAX
        ENDDO
      ENDIF

      IF(INTERACTION.EQ.0) THEN     ! interaction=-1 was used for RC in the past!
        DO I=1,5
          F0NR(I)=DG_TRANS(I,KKK)-DBAR(I)
          DO J=1,5
            XMASTX(I,J)=0.      ! allows to use the same algorithm in Newt-Raph for FC & SC
          ENDDO
        ENDDO
      ELSE IF(INTERACTION.GT.0) THEN

	  IF(ISHAPE(IPH).LE.1) THEN
          DO I=1,5
          DO J=1,5
            XMASTX(I,J)=MASTPH(I,J,IPH)      ! this is M~ calc in VPSC & includes neff factor 
          ENDDO
          ENDDO
        ELSE IF(ISHAPE(IPH).GT.1) THEN
          DO I=1,5
          DO J=1,5
            XMASTX(I,J)=MASTGR(I,J,KKK)      ! this is M~ calc in VPSC & includes neff factor 
          ENDDO
          ENDDO
        ENDIF

        DO I=1,5
          F0NR(I)=DG_TRANS(I,KKK)-(DASTBAR(I)+DBAR0(I))        ! DASTBAR=MBARTG*SASTBAR, calc inside VPSC
          DO J=1,5
            F0NR(I)=F0NR(I)-XMASTX(I,J)*SASTAV(J)
          ENDDO
        ENDDO

      ENDIF

C *** CALLS Newton-Raphson SUBROUTINE TO CALCULATE GRAIN STRESS

      DO I=1,5
        SGXOLD(I)=SGX(I)
      ENDDO
      ITMX=1000
      EPS=5.e-04      ! maximum relative error for convergence

c     iprx=0
c     if(kgx.eq.195 .or. kgx.eq.1288) then
c       write(10,'('' stress in #'',i5,'' before entering N-R'',
c    #             5e12.3)') kgx,x
c       write( *,'('' stress in #'',i5,'' before entering N-R'',
c    #             5e12.3)') kgx,x
c       iprx=1
c     endif

      CALL NEWTON_RAPHSON (SGX,ITMX,EPS,TAULIM,IERROR,iprx)

      IF(IERROR.GT.0) THEN
        IF(IERROR.EQ.1) WRITE(*,'('' SINGULAR SYSTEM IN NEWTRAPH -->'',
     #       '' CANNOT SOLVE GRAIN'',I6,'' IN PHASE'',I6)') KKK,IPH
        IF(IERROR.EQ.2) WRITE(*,'('' ITMAX WAS REACHED IN NEWTRAPH -->'',
     #       '' CANNOT SOLVE GRAIN'',I6,'' IN PHASE'',I6)') KKK,IPH
        WRITE(*,'('' THE INPUT STRESS IS RETAINED'',5e12.3)') SGX(1:5)
        DO I=1,5
          SGX(I)=SGXOLD(I)
        ENDDO
      ENDIF

      DO I=1,5
        SG(I,KKK)=SGX(I)
      ENDDO

      RETURN
      END
C
C ********************************************************************************
C     SUBROUTINE INITIAL_STATE_GUESS     --->     VERSION 15/JAN/2020

C     IF STRAIN IS IMPOSED, USES A STRESS GUESS COLINEAR WITH THE
C           STRAIN RATE AND CALCULATES A GRAIN STRESS GIVEN BY TAYLOR.
C     IF STRESS IS IMPOSED, SETS THE GRAIN STRESS EQUAL TO THE MACROSCOPIC.
C     CALCULATES GRAIN STRAIN RATES 'DG' CALLING EITHER SUBR
C            GRAIN_STRESS OR GRAIN _RATE_AND_MODULI.
C     CALCULATES AVERAGE STRESS 'SAV' AND AVERAGE STRAIN-RATE 'DAV'.
C     CALCULATES VISCO-PLASTIC MODULI 'MGTG' FOR EVERY GRAIN.
C     CALCULATES INITIAL GUESS FOR THE MACROSCOPIC VISCO-PLASTIC MODULUS 'MBARTG'
C ********************************************************************************

      SUBROUTINE INITIAL_STATE_GUESS

      USE VPSC8DIM


      DBARNORM=VNORM(DBAR,6)
	  
      KGX=1
      DO IPH=IPHBOT,IPHTOP
        IPHEL=IPH-IPHBOT+1
        DO KKK=NGR(IPH-1)+1,NGR(IPH)

          IF(STRAIN_CONTROL.EQ.1) THEN

C *** NEED TO ADD A BACK-STRESS ASSOCIATED WITH AN EIGEN-RATE ??
C             AUX6(1:6)=DBAR(1:6)-DG_TRANS(1:6,KGX)
C             DBARNORM=VNORM(AUX6,6)
			
            DO J=1,6
              SGTRY(J,KGX)=DBAR(J)/DBARNORM
            ENDDO
            INTERACTION_SAVE=INTERACTION
            INTERACTION=0
            CALL GRAIN_STRESS (KGX,KKK,IPHEL,IPH)   ! Taylor guess
            INTERACTION=INTERACTION_SAVE

          ELSE IF(STRAIN_CONTROL.EQ.0) THEN

            DO J=1,6
              SG(J,KKK)=SBAR(J)                     ! Sachs guess
            ENDDO
          ENDIF
	  
c *** for each grain calculates DG (using power law), MGTG and DG0  
          CALL GRAIN_RATE_AND_MODULI (0,KGX,KKK,IPHEL,IPH)

C        if (kgx.lt.5) then
C		  write(10,*) ' inside INITIAL_STATE_GUESS   -->   grain', kgx
C		  write(10,'(''  sg='',6f10.3)')  sg(:,kgx)
C		  write(10,'(''  dg='',6f10.3)')  dg(:,kgx)
C		  write(10,'(''  MGTG='',5e12.3)')  MGTG(:,:,kgx)
C		endif	
          KGX=KGX+1
        ENDDO
      ENDDO

      DO I=1,MDIM
        DAV(I)=0.
        SAV(I)=0.
        DO IPH=IPHBOT,IPHTOP
        DO KKK=NGR(IPH-1)+1,NGR(IPH)
          DAV(I)=DAV(I)+DG(I,KKK)*WGT(KKK)
          SAV(I)=SAV(I)+SG(I,KKK)*WGT(KKK)
        ENDDO
        ENDDO
      ENDDO

C     CALCULATE INITIAL GUESS FOR MACROSCOPIC MODULI 'Mtg' AS THE INVERSE
C     OF THE AVERAGE OF THE GRAIN'S STIFFNESSES --> MBARTG=<LBARTGc>^-1
C     ACTUALLY, WHEN INTERACTION=2,3,4 IT AVERAGES THE SECANT MODULI

      DO I=1,MDIM
        DBAR0(I)=0.
        DO J=1,MDIM
          LBARTG(I,J)=0.
          MBARTG(I,J)=0.
        ENDDO
      ENDDO

      KGX=1
      DO IPH=IPHBOT,IPHTOP
      DO KKK=NGR(IPH-1)+1,NGR(IPH)
        DO I=1,MDIM
        DO J=1,MDIM
          AUX66(I,J)=MGTG(I,J,KGX)        ! MGTG=MGSEC if INTERACTION=0,2,3,4
        ENDDO
        ENDDO

        CALL LU_INVERSE(AUX66(1:MDIM,1:MDIM), MDIM)

        DO I=1,MDIM
          DBAR0(I)=DBAR0(I)+(DG0(I,KGX)+DG_TRANS(I,KGX))*WGT(KKK)      
          DO J=1,MDIM
            LBARTG(I,J)=LBARTG(I,J)+AUX66(I,J)*WGT(KKK)
          ENDDO
        ENDDO

        KGX=KGX+1

      ENDDO
      ENDDO

C *** INVERTS LBARTG TO OBTAIN A GUESS FOR THE TANGENT COMPLIANCE MBARTG
C *** MBARTG IS REALLY M_sec FOR THE SECANT, NEFF and TANGENT INTERACTIONS

      MBARTG(1:MDIM,1:MDIM)=LBARTG(1:MDIM,1:MDIM)

      CALL LU_INVERSE(MBARTG(1:MDIM,1:MDIM), MDIM)

C      write(10,*)
C      write(10,'('' leaving initial_state_guess'')')
C      write(10,'('' sav='',6e12.3)') sav
C      write(10,'('' dav='',6e12.3)') dav
C      write(10,'('' MBARTG ='',5e12.3)') MBARTG
C      write(10,*)

      RETURN
      END

C ***************************************************************************
C     SUBROUTINE LANKFORD     --->      VERSION OF 27/JUL/2022
C
C     ROTATES CRYSTALLOGRAPHIC AND MORPHOLOGIC TEXTURE AND REFERS THEM TO THE 
C     TESTING SYSTEM USED FOR LANKFORD PROBES. ASSUMES X1=TENSILE DIRECTION.
C     CALCULATES LANKFORD COEFFICIENT DBAR(2,2)/DBAR(3,3) VERSUS THE ANGLE 
C     BETWEEN THE ROLLING DIRECTION AND THE TENSILE DIRECTION.
C     AT THE END OF LAST STEP ROTATES TEXTURE BACK TO ORIGINAL SAMPLE SYSTEM.
C     CALCULATES DIRECTIONAL ELASTIC YOUNG MODULUS ALONG TENSILE DIRECTION.
C ***************************************************************************

      SUBROUTINE LANKFORD (ISTEP,DELTALANK,IOPTION)

      USE VPSC8DIM

      DIMENSION ROTALPHA(3,3)

      IU=15      ! UNIT FOR WRITING OUTPUT

C ************************************************************************
C    IF IOPTION=0 INITIALIZES RUN AND BOUNDARY CONDITIONS
C ************************************************************************

      IF(IOPTION.EQ.0) THEN

        IF(INTERACTION.EQ.0) THEN
          WRITE(*,*)
          WRITE(*,*) 'CANNOT CALCULATE LANKFORD FOR FC CASE'
          WRITE(*,*) ' --> WILL RESET INTERACTION TO SECANT CASE'
          INTERACTION=2
          print *, 'enter c to continue'
          read  *
        ENDIF

        WRITE(10,*)
        WRITE(10,'(''*** STARTS LANKFORD CALCULATION PROBING EVERY'',
     #               F10.2,'' DEGREES'')') DELTALANK

        NSTEPS=90./DELTALANK
        TINY=ABS(NSTEPS*DELTALANK-90.)
        IF(TINY .GT. 0.1) THEN
          DELTALANK=90./NSTEPS
          WRITE(*,*)
          WRITE(*,*) 'LANKFORD INTERVALS MUST ADD TO 90 DEGREES'
          WRITE(*,*) ' --> WILL RESET DELTALANK TO', DELTALANK
          print *, 'enter c to continue'
          read  *
        ENDIF
        NSTEPS=NSTEPS+1

        ICTRL=1                ! DBARc(1,1) IS THE CONTROL COMPONENT
        STRAIN_CONTROL=1

        ILBAR(:,:)=1           ! ENFORCES ALL SHEARS & D11
        ILBAR(2,2)=0
        ILBAR(3,3)=0

        LIJBARc(:,:)= 0.       ! OFF-DIAGONAL COMPONENTS ENFORCED TO ZERO
        LIJBARc(1,1)= 1.
        LIJBARc(2,2)=-0.5
        LIJBARc(3,3)=-0.5

        IDBARv(:)=1
        IDBARv(2)=0
        IDBARv(3)=0

        ISBARv(:)=0
        ISBARv(2)=1
        ISBARv(3)=1

C --> OPTION: ACTIVATE NEXT 3 LINES FOR ALLOWING IN-PLANE SHEAR (DEC 2021)
c        ILBAR(1,2)=0
c        IDBARv(6) =0
c        ISBARv(6) =1
		
        DO I=1,3
        DO J=1,3
          DBARc(I,J)=(LIJBARc(I,J)+LIJBARc(J,I))/2.
        ENDDO
        ENDDO
        CALL CHG_BASIS (DBAR,DBARc,AUX55,AUX3333,2,6)

        SBARc(:,:)  = 0.
        CALL CHG_BASIS (SBAR,SBARc,AUX55,AUX3333,2,6)
		
        DBARNORM=VNORM(DBAR,6)
        DBAR(:)=DBAR(:)/DBARNORM
        CALL CHG_BASIS (DBAR,DBARc,AUX55,AUX3333,1,6)

      ENDIF

C ************************************************************************
C    IF IOPTION=1 EXPRESSES TEXTURE AND GRAIN SHAPE IN TENSILE AXES
C    AT AN ANGLE 'ALPHA' WITH RESPECT TO ROLLING DIRECTION.
C    MATRIX 'ROTALPHA' ROTATES TEXTURE BY 'ALPHA'.
C    ASSUMES THAT RD=axis1, TD=axis2, ND=axis3
C ************************************************************************

      IF(IOPTION.EQ.1) THEN

        IF(ISTEP.EQ.1) ALPHA= 0.
        IF(ISTEP.GT.1) ALPHA= DELTALANK*PI/180.

        ROTALPHA(1,1)= COS(ALPHA)
        ROTALPHA(1,2)= SIN(ALPHA)
        ROTALPHA(1,3)= 0.
        ROTALPHA(2,1)=-SIN(ALPHA)
        ROTALPHA(2,2)= COS(ALPHA)
        ROTALPHA(2,3)= 0.
        ROTALPHA(3,1)= 0.
        ROTALPHA(3,2)= 0.
        ROTALPHA(3,3)= 1.

        CALL TEXTURE_ROTATION (ROTALPHA)

      ENDIF

C ************************************************************************
C *** IF IOPTION=2 WRITES LANKFORD COEFFICIENT, DIAGONAL COMPONENTS OF THE
C     STRAIN RATE, AND TENSILE AND SHEAR STRESS (IN TENSILE TEST SYSTEM).
C *** CALCULATES DIRECTIONAL ELASTIC YOUNG MODULUS (OCT/9/03). THIS REQUIRES 
C     TO RECALCULATE THE PX SC ELASTIC MODULI FOR EACH ROTATED TEXTURE
C *** FOR THE LAST STEP RESETS TEXTURE AND GRAIN SHAPE BACK TO THE ORIGINAL
C     ROLLING AXES
C ************************************************************************

      IF(IOPTION.EQ.2) THEN

        CALL ELSC(1)      ! SC ELASTIC COMPLIANCE AND STIFFNESS TENSORS

        IF(ISTEP.EQ.1)
     #    WRITE(IU,'('' ANGLE     YOUNG     LANKF       D(1,1)'',
     #      ''    D(2,2)    D(3,3)    D(1,2)       S(1,1)    S(1,2)'')')

        ANGLE=DELTALANK*(ISTEP-1)
        IF(ABS(DBARc(3,3)/DBARc(1,1)).LE.1.E-06) THEN
          RLANK=999999.
        ELSE
          RLANK =DBARc(2,2)/DBARc(3,3)
        ENDIF

C *** CALCULATION OF TAYLOR FACTOR (2022 09 05) **********

      CRSSAVE=0.
      DO IPH=1,NPH
      DO KGX=NGR(IPH-1)+1,NGR(IPH)
      DO IS=1,NSYST(IPH)
        CRSSAVE=CRSSAVE+CRSS(IS,KGX)*WGT(KGX)/NSYST(IPH)
      ENDDO
      ENDDO
      ENDDO
      TAYLORFAC=0.
      DO I=1,3
      DO J=1,3
        TAYLORFAC=TAYLORFAC+SBARc(I,J)*DBARc(I,J)
      ENDDO
      ENDDO
      DBARNORM=VNORM(DBAR,6)
      TAYLORFAC=TAYLORFAC/CRSSAVE/DBARNORM*SQRT(3./2.)

C ********************************************************

        WRITE(IU,'(F6.1,F10.1,F10.2,3X,4F10.4,3x,3F10.3,3X,7F10.3)')
     #      ANGLE, 1./SELAVv(1,1), RLANK, DBARc(1,1),DBARc(2,2),
     #      DBARc(3,3),DBARc(1,2),SBARc(1,1),SBARc(1,2)   ,TAYLORFAC
C     #      ,CRSS(1,1),CRSSAVE,STRNORM,TAYLORFAC

        IF(ISTEP .EQ. NSTEPS) THEN
          ALPHA= -PI/2.
          ROTALPHA(1,1)= COS(ALPHA)
          ROTALPHA(1,2)= SIN(ALPHA)
          ROTALPHA(1,3)= 0.
          ROTALPHA(2,1)=-SIN(ALPHA)
          ROTALPHA(2,2)= COS(ALPHA)
          ROTALPHA(2,3)= 0.
          ROTALPHA(3,1)= 0.
          ROTALPHA(3,2)= 0.
          ROTALPHA(3,3)= 1.
          CALL TEXTURE_ROTATION (ROTALPHA)
        ENDIF

      ENDIF
C ************************************************************************

      RETURN
      END

C ************************************************************************
C     SUBROUTINE LOAD_CONDITIONS      -->     VERSION SEP/10/2021
C
C     READS BOUNDARY CONDITIONS ON STRAIN RATE AND STRESS:
C        * NUMBER OF DEF STEPS, CONTROL VARIABLE, INCREMENT, TEMPERATURE
C        * IMPOSED AND RELAXED COMPONENTS OF VELOCITY GRADIENT AND STRESS
C     CALCULATES SYMMETRIC STRAIN RATE 'DBARc'
C     CALCULATES 5-DIM (TENTATIVE) VECTORS 'DBAR' AND 'SBAR'.
C     CHECKS WHETHER THE BOUNDARY CONDITIONS ARE CONSISTENT.
C *** SLIGHTLY MODIFIED FOR ELASTIC BOUNDARY CONDITIONS (02/SEP/2017)
C ************************************************************************

      SUBROUTINE LOAD_CONDITIONS (UNIT)

      USE VPSC8DIM

      INTEGER UNIT
      DIMENSION PROFAC(6)

      profac(1)=1.
      profac(2)=1.
      profac(3)=1.
      profac(4)=2.
      profac(5)=2.
      profac(6)=2.

C *** WRITES LOAD CONDITIONS FILE INTO 'RUN_LOG.OUT' FILE
      WRITE(10,*)
      WRITE(10,'(''*** LOAD CONDITIONS FOR THIS RUN'')')
      DO IDUM=1,100
        READ(UNIT=UNIT,END=100,FMT='(A)') PROSA
        WRITE(10,'(A)') PROSA
      ENDDO
  100 REWIND UNIT

C *** READS BOUNDARY CONDITIONS ON DEFORMATION RATE, STRESS COMPONENTS
C     AND TEMPERATURE. A TEMPERATURE INCREMENT CAN BE IMPOSED IN EACH
C     STEP, BUT TEMPERATURE CONTROL (ICTROL=8) CAN ONLY BE USED FOR THE
C     THERMO-ELASTIC CASE (IVGVAR=-1)

      READ(UNIT,*) NSTEPS,ICTRL,CTRLINCR,TEMP_INI,TEMP_FIN
	  TEMP=TEMP_INI
      IF(ICTRL.LT.0 .OR. ICTRL.EQ.7 .OR. ICTRL.GT.8) THEN
        WRITE(*,'('' THE VALUE OF ICTRL IS NOT VALID'')')
        STOP
      ENDIF
	  IF(ICTRL.NE.8) TEMP_INCR=(TEMP_FIN-TEMP_INI)/NSTEPS
	  IF(ICTRL.EQ.8) TEMP_INCR= CTRLINCR
      IF(CTRLINCR.LE.0 .AND. ICTRL.NE.8) THEN
        WRITE(*,'('' STRAIN OR TIME INCREMENT (CTRLINCR) HAS TO BE'')')
        WRITE(*,'('' POSITIVE  --> THE SIGN OF THE DEFORMATION IS'')')
        WRITE(*,'('' CONTROLLED BY THE IMPOSED RATE COMPONENT'')')
        STOP
      ENDIF

      READ(UNIT,*)
      DO I=1,3
        READ(UNIT,*) (ILBAR(I,J),J=1,3)
      ENDDO
      READ(UNIT,*)
      DO I=1,3
        READ(UNIT,*) (LIJBARc(I,J),J=1,3)
      ENDDO
      READ(UNIT,*)
      READ(UNIT,*) ISBARv(1),ISBARv(6),ISBARv(5)
      READ(UNIT,*)           ISBARv(2),ISBARv(4)
      READ(UNIT,*)                     ISBARv(3)
      READ(UNIT,*)
      READ(UNIT,*) SBARc(1,1),SBARc(1,2),SBARc(1,3)
      READ(UNIT,*)            SBARc(2,2),SBARc(2,3)
      READ(UNIT,*)                       SBARc(3,3)
        SBARc(3,2)=SBARc(2,3)
        SBARc(3,1)=SBARc(1,3)
        SBARc(2,1)=SBARc(1,2)

C **********************************************************
C *** CHECKS WHETHER THE BOUNDARY CONDITIONS ARE CONSISTENT
C **********************************************************

      ISBARSUM=ISBARv(1)+ISBARv(2)+ISBARv(3)+
     #         ISBARv(4)+ISBARv(5)+ISBARv(6)
      IF(INTERACTION.EQ.0 .AND. ISBARSUM.NE.0) THEN
        WRITE(*,*) ' CANNOT IMPOSE STRESS COMPONENTS FOR FC CASE'
        STOP
      ENDIF

      IF(INTERACTION.NE.-1) THEN      ! elastic deform is compressible 

      IF(ILBAR(1,1)+ILBAR(2,2)+ILBAR(3,3).eq.2) then
        write(*,*) ' CHECK DIAGONAL BOUNDARY CONDITIONS ILBAR:'
        write(*,*) ' ENFORCING TWO DIAGONAL COMPONENTS FIXES THE'
        write(*,*) ' REMAINING COMPONENT BECAUSE OF INCOMPRESSIBILITY'
        STOP
      ENDIF

      DILAT =LIJBARc(1,1)+LIJBARc(2,2)+LIJBARc(3,3)
      IDILAT=ILBAR(1,1)*ILBAR(2,2)*ILBAR(3,3)
      IF(IDILAT*DILAT.GT.1.E-6) THEN
        WRITE(*,*) 'CHECK DIAGONAL STRAIN RATE COMPONENTS LIJBARc'
        WRITE(*,*) '--> THE IMPOSED RATE IS NOT INCOMPRESSIBLE'
        STOP
      ENDIF

      ENDIF      ! end of IF(INTERACTION.NE.-1)

      DO I=1,2
      DO J=I+1,3
        if(ILBAR(i,j)+ILBAR(j,i).eq.0) then
          write(*,*) ' CHECK OFF-DIAGONAL BOUNDARY CONDITIONS ILBAR'
          write(*,*) ' CANNOT RELAX BOTH OFF-DIAGONAL VEL GRAD COMPS'
          stop
        endif
      ENDDO
      ENDDO

C *** CALCULATES SYMMETRIC STRAIN-RATE COMPONENTS FROM VELOCITY GRADIENT.
C *** TAGS IMPOSED & UNKNOWN CARTESIAN COMPONENTS

      DO I=1,3
        DO J=1,3
          DBARc(I,J)=(LIJBARc(I,J)+LIJBARc(J,I))/2.
        ENDDO
      ENDDO

      IDBARv(1)=ILBAR(1,1)
      IDBARv(2)=ILBAR(2,2)
      IDBARv(3)=ILBAR(3,3)
      IDBARv(4)=ILBAR(2,3)*ILBAR(3,2)
      IDBARv(5)=ILBAR(1,3)*ILBAR(3,1)
      IDBARv(6)=ILBAR(1,2)*ILBAR(2,1)

      IDBARSUM=IDBARv(1)+IDBARv(2)+IDBARv(3)+
     #         IDBARv(4)+IDBARv(5)+IDBARv(6)

      DO I=1,6
        IF((ISBARv(I)*IDBARv(I)).NE.0 .OR. 
     #     (ISBARv(I)+IDBARv(I)).NE.1) THEN
          WRITE(*,*) ' BOUNDARY CONDITIONS ON STRAIN-RATE AND STRESS'
          WRITE(*,*) ' COMPONENTS MUST BE COMPLEMENTARY'
          WRITE(*,'(''   IDBARv = '',6I3)') IDBARv
          WRITE(*,'(''   ISBARv = '',6I3)') ISBARv
          STOP
        ENDIF
      ENDDO

C ***********************************************************************	  
C *** CALCULATES b-BASIS & Voigt REPRESENTATION OF MACRO STRAIN RATE AND
C     STRESS CARTESIAN TENSORS 'DBARc' AND 'SBARc'.
C *** while ONLY ONE COMPLEMENTARY COMPONENT IN EACH TENSOR IS ENFORCED,
C     INPUT IS USED TO MAKE 1st STEP GUESS FOR THE SC PROCEDURE.
C *** FOR THE ELASTIC CASE A GUESS IS NOT REQUIRED AND THE NON-IMPOSED  
C     COMPONENTS ARE SET TO ZERO TO AVOID 1st STEP INCONSISTENCIES.

      CALL VOIGT (SBARv,SBARc,AUX66,AUX3333,2)
      CALL VOIGT (DBARv,DBARc,AUX66,AUX3333,2)

      KDIM=5
      IF(INTERACTION.EQ.-1) THEN        ! elasticity option
        KDIM=6
        SBARv(:) = SBARv(:)*ISBARv(:)	
        DBARv(:) = DBARv(:)*IDBARv(:)	
        CALL VOIGT (SBARv,SBARc,AUX66,AUX3333,1)
        CALL VOIGT (DBARv,DBARc,AUX66,AUX3333,1)
      ENDIF

      CALL CHG_BASIS(DBAR,DBARc,AUX66(1:KDIM,1:KDIM),AUX3333,2,KDIM)
      CALL CHG_BASIS(SBAR,SBARc,AUX66(1:KDIM,1:KDIM),AUX3333,2,KDIM)
	  
C ***********************************************************************
C     IN CONNECTION WITH THE OVERALL STRAIN RATE AND STRESS:
C     *  IF ICTRL=0 WILL IMPOSE A VON MISES STRAIN INCREMENT. THIS OPTION
C        VALID ONLY FOR FULLY IMPOSED STRAIN TENSOR AND PLASTICITY. 
C        WILL CALCULATE THE TIME INCREMENT NECESSARY TO ACHIEVE IT
C     *  IF ICTRL=1,6 WILL EITHER
C        --> IMPOSE A STRAIN RATE COMPONENT AND A STRAIN INCREMENT & 
C            CALCULATE THE TIME INCREMENT REQUIRED TO ACHIEVE IT OR,
C        --> IMPOSE A STRESS COMPONENT (CREEP) AND A TIME INCREMENT. 
C            THE STRAIN INCREMENT THAT RESULTS FOLLOWS FROM THE VPSC
C            CALCULATED STRAIN RATE TIMES THE TIME INCREMENT.
C     *  THE OPTION ICTRL=7 MEANT TO READ AN IMPOSED TIME INCREMENT WAS
C        ELIMINATED SINCE IMPOSING A STRAIN RATE AND A STRAIN INCREMENT
C        ALLOWS TO CALCULATE THE REQUIRED TIME INCREMENT.
C     *  IF ICTRL=8 WILL IMPOSE A TEMPERATURE INCREMENT (=CTRLINCR) AT
C           EACH STEP. BC's STRESS AND STRAIN COMPONENTS MAY STILL BE
C           IMPOSED. THIS OPTION IS VALID ONLY FOR THERMO-ELASTIC CASE.
C ***********************************************************************

      IF(ICTRL.EQ.0) THEN
        STRAIN_CONTROL=1
        IF(IDBARSUM.NE.6. OR. INTERACTION.EQ.-1) THEN
          WRITE(*,*) 'CAN CONTROL VON MISES ONLY FOR PLASTIC SIMULATION'
          WRITE(*,*) 'AND WHEN ALL STRAIN COMPONENTS ARE IMPOSED !'
          STOP
        ENDIF
        EVMINCR  =CTRLINCR
        TIME_INCR=EVMINCR/VNORM(DBAR,5)
      ELSE IF(ICTRL.GE.1 .AND. ICTRL.LE.6) THEN
        IF(IDBARv(ICTRL).EQ.1) STRAIN_CONTROL=1
        IF(ISBARv(ICTRL).EQ.1) STRAIN_CONTROL=0
        IF (STRAIN_CONTROL.EQ.1)THEN
		  DBARCTRL =DBARv(ictrl)
          IF(DBARCTRL.EQ.0.) THEN
            WRITE(*,*) 'STRAIN CONTROL COMPONENT CANNOT BE ZERO !'
            STOP
          ENDIF
          EIJINCR=CTRLINCR
          TIME_INCR=EIJINCR/ABS(DBARCTRL)
        ELSE IF (STRAIN_CONTROL.EQ.0)THEN
		  SBARCTRL =SBARv(ictrl)
          IF(SBARCTRL.EQ.0.) THEN
            WRITE(*,*) 'STRESS CONTROL COMPONENT CANNOT BE ZERO !'
            STOP
          ENDIF
		  TIME_INCR=CTRLINCR
		ENDIF
      ELSE IF(ICTRL.EQ.8) THEN       ! TENTATIVE / VALID FOR THERMO-ELASTIC ONLY
        TEMP_INCR=CTRLINCR
        TIME_INCR=1.
      ENDIF
  
C *** IN ELASTICITY THE INITIAL STRESS-STRAIN CONDITIONS AT THE BEGINNING 
C     OF STEP #1 ARE KNOWN, AND USED TO DEFINE STRESS-STRAIN AT TIME=0.
C *** ASSUMES THAT EPSTOTv(i) IS ONLY ELASTIC (NO THERMAL COMP) AT TIME=0.

      IF(INTERACTION.EQ.-1) THEN
        EPSTOTv(:)=0.
        DO I=1,6
        DO J=1,6
          EPSTOTv(I)=EPSTOTv(I)+SELAVv(I,J)*SBARv(J)*PROFAC(J)
        ENDDO
        ENDDO
        CALL VOIGT(EPSTOTv,EPSTOTc,AUX66,AUX3333,1)

	    CALL WRITE_STRESS_STRAIN (0)
        DO IPH=IPHBOT,IPHTOP             
          IF(IDIFF(IPH).EQ.1) CALL DIFF_PLANES(IPH,ICRYSYMPH(IPH),1)
        ENDDO
      ENDIF

      RETURN
      END

C *****************************************************************************
      SUBROUTINE LU_BACKSUBS(a,n,np,indx,b)

      INTEGER n,np,indx(n)
      REAL a(np,np),b(n)
      INTEGER i,ii,j,ll
      REAL sum

      ii=0
      do i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do j=ii,i-1
            sum=sum-a(i,j)*b(j)
          enddo
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
      enddo
      do i=n,1,-1
        sum=b(i)
        do j=i+1,n
          sum=sum-a(i,j)*b(j)
        enddo
        b(i)=sum/a(i,i)
      enddo
      return
      END

C *****************************************************************************
      SUBROUTINE LU_DECOMP (a,n,np,indx,d,isingular)

      INTEGER n,np,indx(n),NMAX
      REAL d,a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1.0e-20)
      INTEGER i,imax,j,k,isingular
      REAL aamax,dum,sum,vv(NMAX)

      d=1.
      do i=1,n
        aamax=0.
        do j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
        enddo

c        if (aamax.eq.0.) pause 'singular matrix in lu_decomp'

        if(aamax.eq.0.) then
        isingular=1
        return
        endif

        vv(i)=1./aamax
      enddo
      do j=1,n
        do i=1,j-1
          sum=a(i,j)
          do k=1,i-1
            sum=sum-a(i,k)*a(k,j)
          enddo
          a(i,j)=sum
         enddo
        aamax=0.

        do i=j,n
          sum=a(i,j)
          do k=1,j-1
            sum=sum-a(i,k)*a(k,j)
          enddo
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
        enddo
        if (j.ne.imax)then
          do k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
          enddo
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
c
c        if(a(j,j).eq.0.) a(j,j)=TINY
c
        if(a(j,j).eq.0.) then
        isingular=1
        return
        endif
c
        if(j.ne.n)then
          dum=1./a(j,j)
          do i=j+1,n
            a(i,j)=a(i,j)*dum
         enddo
        endif
      enddo
c
      isingular=0
c
      return
      END

C *************************************************************************
      SUBROUTINE LU_EQSYSTEM (A,B,N,ISINGULAR)

C *** SOLVES A*X=B USING LU DECOMPOSITION

      DIMENSION A(N,N),B(N),INDX(N)     

      CALL LU_DECOMP(A,N,N,INDX,D,ISINGULAR)

      IF(ISINGULAR.EQ.1) RETURN

      CALL LU_BACKSUBS(A,N,N,INDX,B)

      RETURN
      END

C *************************************************************************
      SUBROUTINE LU_INVERSE (A,N)

C *** INVERTS A MATRIX USING LU DECOMPOSITION

      DIMENSION A(N,N),Y(N,N),INDX(N),ax(n,n)     

c     write(*,*) 'A(i,j) matrix inside lu_inverse'
c     write(*,'(5e12.3)') ((a(i,j),j=1,n),i=1,n)
c     pause

C **************************************************************
C   NORMALIZATION ADDED 03/DEC/05 TO AVOID NUMERICALLY SINGULAR MATRIX
      AMAX=0.
      DO I=1,N
      DO J=1,N
        DUM=ABS(A(I,J))
        IF(DUM .GT. AMAX) AMAX=DUM
      ENDDO
      ENDDO
      DO I=1,N
      DO J=1,N
        A(I,J)=A(I,J)/AMAX      ! normalize the matrix
      ENDDO
      ENDDO
	  AX(:,:)=A(:,:)
C **************************************************************

      DO I=1,N
        DO J=1,N
          Y(I,J)=0.
        ENDDO
        Y(I,I)=1.
      ENDDO

      CALL lu_decomp(A,N,N,INDX,D,ISINGULAR)
	  
      IF(ISINGULAR.EQ.1) THEN
        WRITE(*,*) ' *** SINGULAR MATRIX IN LU_INVERSE !!'
        WRITE(*,'(5E12.3)') ((AX(I,J),J=1,N),I=1,N)
        WRITE(10,*) ' *** SINGULAR MATRIX IN LU_INVERSE !!'
        WRITE(10,'(5E12.3)') ((AX(I,J),J=1,N),I=1,N)
        STOP
      ENDIF

      DO J=1,N
        CALL lu_backsubs(A,N,N,INDX,Y(1,J))
      ENDDO

      DO I=1,N
      DO J=1,N
        A(I,J)=Y(I,J) /AMAX      ! renormalize the inverse
      ENDDO
      ENDDO

      RETURN
      END

C ***************************************************************************
C     SUBROUTINE NEIGHBOURS      --->      VERSION 25/JUL/2002
C
C *** IDENTIFIES 'NNEIGH' NEIGHBORS FOR EVERY GRAIN AND ASSIGNS WEIGHTS
C     'WNEIGH' WHICH ARE GOING TO BE USED EXCLUSIVELY FOR CO-ROTATION
C     PURPOSES INSIDE SUBROUTINE UPDATE_ORIENTATION.
C *** IF IOPTION=0 DEFINES NEIGHBORS FOR EACH GRAIN PICKED AT RANDOM FROM
C     ORIENTATION FILE (REORIENTATION WORKS ALSO FOR MULTIELEMENT RUN)
C *** IF IOPTION=1 DEFINES NEIGHBORS FOR A TWO PHASE MATERIAL WHERE THE
C     SECOND PHASE IS CRYSTALLOGRAPHICALLY ASSOCIATED TO THE FIRST PHASE
C     (i.e.: MATRIX/TWIN or PARENT/CHILD IN RECRYSTALLIZATION)
C ***************************************************************************

      SUBROUTINE NEIGHBOURS (IOPTION)

      USE VPSC8DIM

      DIMENSION WAUX(0:NNEIMX)


      DO I=0,NNEIGH
        WAUX(I)=1./(NNEIGH+1)       ! DEFINES EQUAL CO-ROTATION WGT's
      ENDDO

      DO KGX=1,NGR(NPH)             ! SCANS ALL GRAINS IN ALL PHASES

        NEIGH(0,KGX)=KGX            ! NEIGHBOR ZERO IS THE GRAIN ITSELF
        WTOT=WGT(KGX)*WAUX(0)

        IF(NNEIGH.GT.0) THEN

          IF(IOPTION.EQ.0) THEN             ! PICK NEIGHBOURS AT RANDOM
            DO I=1,NNEIGH
              RAND=random2(JRAN)
              NEIGH(I,KGX)=INT(RAND*NGR(NPH)+1)
              WTOT=WTOT+WGT(NEIGH(I,KGX))*WAUX(I)
            ENDDO
          ELSE IF(IOPTION.EQ.1) THEN        ! PICK NEIGHBOURS IN ORDER
            IF (KGX.LE.NGR(1)) NEIGH(1,KGX)=KGX+NGR(1)
            IF (KGX.GT.NGR(1)) NEIGH(1,KGX)=KGX-NGR(1)
            WTOT=WTOT+WGT(NEIGH(1,KGX))*WAUX(1)
          ENDIF

        ENDIF

C *** NORMALIZATION OF CO-ROTATION WGT's

        IF(WTOT.LT.1.E-10) THEN
          DO I=0,NNEIGH
            WNEIGH(I,KGX)=0.
          ENDDO
        ELSE
          DO I=0,NNEIGH
            WNEIGH(I,KGX)=WAUX(I)*WGT(NEIGH(I,KGX))/WTOT
          ENDDO
        ENDIF

      ENDDO       !  END LOOP OVER ORIENTATIONS

      END

C ********************************************************************
C     SUBROUTINE N_EFFECTIVE  --->   VERSION 19/OCT/2021
C
C     CALCULATES THE COMPONENT OF THE MACROSCOPIC SECANT COMPLIANCE
C     THAT INTERACTS WITH THE MACROSCOPIC STRAIN RATE TENSOR. 
C     CALCULATES THE COMPONENT OF EACH GRAIN SECANT COMPLIANCE
C     THAT INTERACTS WITH THE MACROSCOPIC STRAIN RATE TENSOR.
C     USES THE RATIO OF THE TWO (RELATIVE DIRECTIONAL COMPLIANCE) TO
C     ASSIGN AN EFFECTIVE VALUE OF 'n' FOR THE INTERACTION EQUATION.
C     SUCH VALUE IS IN THE INTERVAL: 1 < NEFFGR < nrs(1,1)
C ********************************************************************

      SUBROUTINE N_EFFECTIVE (ISTEP)

      USE VPSC8DIM
      USE NEFFSAVE

      DIMENSION MSEC_PROJ(5,5),MG_PROJ(5,5),PROJ(5,5)
      DIMENSION RDCg(NGRPEL),RR(NGRPEL)
	  REAL*8    MSEC_PROJ,MG_PROJ
	  
C *** CALCULATES PROJECTION TENSOR OVER STRAIN-RATE DIRECTION
      DAVNORM=VNORM(DAV,5)
      DO I=1,5
      DO J=1,5
        PROJ(I,J)=DAV(I)*DAV(J)/DAVNORM**2
      ENDDO
      ENDDO

C *** CALCULATES PROJECTION OF OVERALL COMPLIANCE AGAINST STRAIN-RATE.
C --> CAVEAT! FOR SEC, TAN & AFFINE 'MBARTG' IS ACTUALLY THE MACROSCOPIC SECANT
      DO I=1,5
      DO J=1,5
        MSEC_PROJ(I,J)=0.
        DO K=1,5
          MSEC_PROJ(I,J)=MSEC_PROJ(I,J)+PROJ(I,K)*MBARTG(K,J)
        ENDDO
      ENDDO
      ENDDO
      XMSNORM=TNORM(MSEC_PROJ,5)

C *** CALCULATES PROJECTION OF GRAIN'S COMPLIANCE AGAINST STRAIN-RATE.
C *** FOR EACH GRAIN CALCULATES RELATIVE DIR COMPLIANCE 'RC' & RELATIVE RATE 
C     'RR' WITH RESPECT TO MACRO COMPLIANCE AND MACRO STRAIN RATE
C *** CALCULATES AVERAGE AND STANDARD DEVIATION OF REL DIR COMPLIANCE OVER
C     ALL THE GRAINS.

      RDCave=0.
      RDCdev=0.
      KGX=1
      DO IPH=1,NPH
      DO KKK=NGR(IPH-1)+1,NGR(IPH)
        RR(KGX)=0.
        DO I=1,5
        RR(KGX)=RR(KGX)+DG(I,KGX)*DAV(I)/DAVNORM**2
        DO J=1,5
          MG_PROJ(I,J)=0.
          DO K=1,5
            MG_PROJ(I,J)=MG_PROJ(I,J)+PROJ(I,K)*MGSEC(K,J,KGX)
          ENDDO
        ENDDO
        ENDDO
        RDCg(KGX)=TNORM(MG_PROJ,5)/XMSNORM
        RDCave=RDCave+RDCg(KGX)    *WGT(KKK)
        RDCdev=RDCdev+RDCg(KGX)**2 *WGT(KKK)
        KGX=KGX+1
      ENDDO
      ENDDO
      RDCdev=SQRT(RDCdev-RDCave**2)

C *** ONLY RDCmax & RDCmin ARE REQUIRED BY THE REL DIR COMPLIANCE CRITERION.
C *** THE CHOICE MADE BELOW IS TO USE THE COMPLIANCE VALUES ONTAINED AFTER 
C     THE FIRST STEP AND TO USE THEM THROUGH THE REST OF THE SIMULATION.
C     A QUICK & DIRTY CORRECTION (below) PREVENTS GOING OUT OF THOSE BOUNDS.

      IF(ISTEP.EQ.1) THEN
        RDCmax=RDCg(1)
        RDCmin=RDCg(1)
        KGX=1
        DO IPH=1,NPH
        DO KKK=NGR(IPH-1)+1,NGR(IPH)
          IF(RDCg(KGX).LE.RDCmin) THEN
		    RDCmin=RDCg(KGX)
			KGMIN=KGX
          ELSE IF(RDCg(KGX).GE.RDCmax) THEN 
		    RDCmax=RDCg(KGX)
            KGMAX=KGX
          ENDIF
          KGX=KGX+1
        ENDDO
        ENDDO
        write(10,'(''kg,NEFFGR,RC_min,RelRate'',i7,3f10.3)')
     #         KGMIN,NEFFGR(KGMIN),RDCg(KGMIN),rr(KGMIN)
        write(10,'(''dg(KGMIN)'',5f10.5)') (dg(i,KGMIN),i=1,5)
        write(10,'(''kg,NEFFGR,RC_max,RelRate'',i7,3f10.3)')
     #         KGMAX,NEFFGR(KGMAX),RDCg(KGMAX),rr(KGMAX)
        write(10,'(''dg(KGMAX)'',5f10.5)') (dg(i,KGMAX),i=1,5)
      ENDIF

C *** DOES A POWER INTERPOLATION BETWEEN 1 AND NRS(1,1)
C   * WHEN RC=RDCmax --> MAX GRAIN COMPLIANCE --> NEFFGR=1        --> SECANT
C   * WHEN RC=RDCmin --> MIN GRAIN COMPLIANCE --> NEFFGR=NRS(1,1) --> TANGENT
C
      IPOWER=1
      AVNEFF=0.
      KGX=1
      DO IPH=1,NPH
      DO KKK=NGR(IPH-1)+1,NGR(IPH)
C        IF (RDCg(KKK).GE.RDCmax) RDCg(KKK)=RDCmax      ! should be equivalent
C        IF (RDCg(KKK).LE.RDCmin) RDCg(KKK)=RDCmin      ! should be equivalent
        X=(RDCmax-RDCg(KKK))/(RDCmax-RDCmin)
        IF (X.LE.0.) X=0.
        IF (X.GE.1.) X=1.
        NEFFGR(KGX)=1+(NRS(1,1)-1)*X**IPOWER
        AVNEFF=AVNEFF+NEFFGR(KGX)*WGT(KKK)
        KGX=KGX+1
      ENDDO
      ENDDO

        write(10,'(''nrs,avneff,RDC_min,RDCave,RDCmax,RDCdev'',
     #         i4,f7.3,4e12.4)')
     #         nrs(1,1),avneff,RDC_min,RDCave,RDCmax,RDCdev
      iskip=1
      if(iskip.eq.0) then
        write(10,'(''kg,NEFFGR,RC_min,RelRate'',i7,3f10.3)')
     #         KGMIN,NEFFGR(KGMIN),RDCg(KGMIN),rr(KGMIN)
        write(10,'(''dg(KGMIN)'',5f10.5)') (dg(i,KGMIN),i=1,5)
        write(10,'(''kg,NEFFGR,RC_max,RelRate'',i7,3f10.3)')
     #         KGMAX,NEFFGR(KGMAX),RDCg(KGMAX),rr(KGMAX)
        write(10,'(''dg(KGMAX)'',5f10.5)') (dg(i,KGMAX),i=1,5)
      endif

      RETURN
      END
	  
C *********************************************************************************
C     SUBROUTINE NEWTON_RAPHSON    --->   VERSION 25/AUG/2019
C
C     CALLED BY SUBR GRAIN_STRESS. ALL RELEVANT CRYSTAL ARRAYS ARE PASSED THROUGH
C     MODULE 'NLNR'.
C
C     IF INTERACTION= 0 SOLVES NON-LINEAR POWER LAW FOR A TAYLOR CASE.
C        Dg-Dbar = 0   or, written explicitly: 
C        FNR=sum(mij*(m*sg/crss)^n) +dg_trans -Dbar =0
C     IF INTERACTION> 0 SOLVES INTERACTION EQ FOR A SELF-CONSISTENT CASE.
C        Dg-Dbar = -M~ * (Sg-Sbar)   or, written explicitly: 
C        FNR=sum(mij*(m*sg/crss)^n) +dg_trans -Mbar*Sbar-D0bar + M~*Sg - M~*Sbar =0
C
C     GIVEN AN INPUT GRAIN STRESS 'SGX' SOLVES THE VISCOPLASTIC EQUATION
C     USING NEWTON-RAPHSON LINEARIZATION AND ITERATING UNTIL CONVERGENCE.
C     - Fi(X)   : ARE THE FUNCTIONS TO MINIMIZE.
C     - FGRADij : ARE THE DERIVATIVES OF Fi WITH RESPECT TO Xj.
C     CORRECTIONS TO 'X' AND RELAXATION OF TOLERANCE ELIMINATED (FEB/2000).
C     RELAXED CONSTRAINTS (IRC=1) ADDED BUT NOT A STANDARD OPTION) (APR/2003).
C
C ***************************************************************************************************

      SUBROUTINE NEWTON_RAPHSON (SGX,KMAX,EPS,TAULIM,IERROR ,iprx)

      USE VPSC8DIM
      USE NLNR
	  

      DIMENSION FGRAD(5,5),F(5),FGRADX(3,3),FX(3),SGX(5),SGXOLD(5)
      DIMENSION RSS(NSYSMX),GD(NSYSMX)

      COEF=0.2
      IERROR=0

      SGXOLD(1:5)=SGX(1:5)
	  
      DO 1000 K=1,KMAX

C *** RESOLVED SHEARS CALCULATION & OUTSIDE Y.S. ERROR MANAGING BLOCK.
C *** NRS MAY BE EVEN OR ODD:
C     GD ALWAYS > 0 --> GD IS USED TO GET DERIVATIVES TO BUILD THE
C     COEFFICIENT MATRIX FOR N-R METHOD IN F(I) CALCULATION (INDEPENDENT
C     TERM FOR N-R METHOD)   RSS*GD=GAMDOT -> SIGN(RSS*GD)=SIGN(RSS)

        DO IS=1,NSYSTX
          RSS(IS)=SCHX(1,IS)*SGX(1)+SCHX(2,IS)*SGX(2)+SCHX(3,IS)*SGX(3)+
     #                 SCHX(4,IS)*SGX(4)+SCHX(5,IS)*SGX(5)
          IF(.NOT.(RSS(IS).GT.0. .OR. ISENSEX(IS).EQ.1)) RSS(IS)=0.
          RSS(IS)=RSS(IS)/TAUX(IS)
          IF(ABS(RSS(IS)).LT.1.E-10) RSS(IS)=0.

          GD(IS)=GAMD0X(IS)
          IF(NRSX(IS).NE.1) GD(IS)=GD(IS)*ABS(RSS(IS)**(NRSX(IS)-1))

          if(ihardlaw.ne.31 .and. ihardlaw.ne.32) then     ! AP:
            IF(ABS(RSS(IS)).GT.TAULIM) THEN
              DO I=1,5
                SGX(I)=SGXOLD(I)+COEF*(SGX(I)-SGXOLD(I))
              ENDDO
              GO TO 1000
            ENDIF
          endif
        ENDDO

c       if(iprx.eq.1) then
c         write( *,'('' db & rss='',5f10.4)') db
c         write( *,'(6f10.4)') (rss(is),is=1,nsystx)
c         write(10,'('' db & rss='',5f10.4)') db
c         write(10,'(6f10.4)') (rss(is),is=1,nsystx)
c         write(*,'(''nrsx & xmast='',30i3)') (nrsx(is),is=1,nsystx)
c         write(*,'(5f13.7)') ((xmastx(i,j),j=1,5),i=1,5)
c       endif

        DO I=1,5
          F(I)=F0NR(I)            ! F0NR is calculated inside GRAIN_STRESS
          DO J=1,5
            F(I)=F(I)+XMASTX(I,J)*SGX(J)
          ENDDO
        ENDDO
        DO I=1,5
          DO IS=1,NSYSTX
            F(I)=F(I)+SCHX(I,IS)*RSS(IS)*GD(IS)
          ENDDO
        ENDDO
        DO I=1,5
          DO J=1,5
            FGRAD(I,J)=XMASTX(I,J)
          ENDDO
        ENDDO

        DO IS=1,NSYSTX
          SFACTOR=NRSX(IS)*GD(IS)/TAUX(IS)
          DO I=1,5
          DO J=1,5
            FGRAD(I,J)=FGRAD(I,J)+SFACTOR*SCHX(I,IS)*SCHX(J,IS)
          ENDDO
          ENDDO
        ENDDO

C *** SOLVES 5X5 LINEAR SYSTEM TO GET DEVIATORIC STRESS COMPONENTS
C          F(SG)=F(SG0)+FGRAD(SG0)*(SG-SG0) = 0
C          (SG-SG0) = (-FGRAD(SG0))^(-1)  F(SG0)

        if(iprx.eq.1) then
          write(10,'('' f    ='',5e12.3)') f
          write(10,'('' fgrad='',5e12.3,/,(6x,5e12.3))') fgrad
        endif

		SGXOLD(1:5)=SGX(1:5)
        FGRAD(1:5,1:5) = -FGRAD(1:5,1:5)
        CALL LU_EQSYSTEM (FGRAD,F,5,ISINGULAR)      ! input F(sg0) is overriten with (sg-sg0) 
        
		IF(ISINGULAR.EQ.1) THEN
          if(iprx.eq.1) then
            write( *,'('' f    ='',5e12.3)') f
            write( *,'('' fgrad='',5e12.3,/,(6x,5e12.3))') fgrad
            write(10,'('' f    ='',5e12.3)') f
            write(10,'('' fgrad='',5e12.3,/,(6x,5e12.3))') fgrad
          endif
          IERROR=1
          RETURN
        ENDIF
			
        DO I=1,5
          SGX(I)=SGXOLD(I)+F(I)
        ENDDO

C *** BOUNDS THE STRESS CORRECTION TO AVOID LARGE OSCILATIONS IN CONVERGENCE
c         if(iprx.eq.1) then
c         rcorr=VNORM(F,5)/VNORM(SGXOLD,5)
c         write(10,'(''*** NR correction'',5f9.3,f12.5)')
c    #                 (f(i),i=1,5),rcorr
c         endif

        RERROR=VMISMATCH(SGX,SGXOLD,5)
        IF(RERROR.LT.EPS) RETURN

 1000   CONTINUE      ! end of DO 1000 K=1,KMAX
C *******************************************************************

      IERROR=2

      RETURN
      END

C ************************************************************************
C     SUBROUTINE PCYS      --->      VERSION 03/APR/2023
C
C  IF IOPTION=0 GENERATES AND STORES EQUISPACED STRAIN-RATE VECTORS TO
C     PROBE THE POLYCRYSTAL YIELD SURFACE IN A SUBSPACE
C  IF IOPTION=1 FEEDS 'DBAR' FOR EACH PROBING STEP
C  IF IOPTION=2 NORMALIZES STRESS 'SPCYS' AND RATE 'DPCYS' IN ORDER FOR
C     EACH POINT OF THE PCYS TO REPRESENT THE SAME DISSIPATION RATE,
C     ARBITRARILY CHOSEN HERE TO BE WRATEREF=DAV(i)*SAV(i)=1.
C
C   STRAIN RATE PROBING IS HARDWIRED --> STRAINR_PROBE=.TRUE.
C   STRESS PROBING ADDED BY CNT ON APRIL 30 2020 BUT NOT THOROUGHLY TESTED
C      --> SWITCH BELOW TO STRAINR_PROBE=.FALSE.
C   ADDED ALTERNATIVE DEFINITION OF 'WRATEAV' ON APRIL 20 2022
C ************************************************************************

      SUBROUTINE PCYS (ISTEP,INDX,INDY,IOPTION)

      USE VPSC8DIM
      USE PCYSSAVE

      DIMENSION ITHETA(2,5),UNITDIR(5)      
      DIMENSION SUMTH(0:4),XCOEF(3,3),XVECT(3),XSOL(3),
     #          YSTANG(2),YSNORM(2),TICK(2,NPROBE),ANGN(NPROBE),
     #          RADS(-2:NPROBE+3),ANGS(-2:NPROBE+3),
     #          RADD(-2:NPROBE+3),ANGD(-2:NPROBE+3) 
 	  LOGICAL   STRAINR_PROBE

      IU=14      ! UNIT FOR WRITING OUTPUT INTO 'PCYS.OUT'

      STRAINR_PROBE=.TRUE.      ! strain rate probing --> hardwired
c      STRAINR_PROBE=.FALSE.    ! activate for stress probing 
	  
C *********************************************************************
      IF(IOPTION.EQ.0) THEN

        OPEN(141,file='PCYSwRATES.OUT'  ,status='UNKNOWN')
        OPEN(142,file='PCYSwNORMALS.OUT',status='UNKNOWN')

        WRITE(10,*)
        WRITE(10,*) '******* STARTS A PCYS CALCULATION *******'
        WRITE(10,*)
        IF(ICS.EQ.1) WRITE(10,*) 'CENTRO-SYMMETRIC YIELD SURFACE'
        IF(ICS.EQ.0) WRITE(10,*) 'NON CENTRO-SYMMETRIC YIELD SURFACE'

        IF(STRAINR_PROBE) THEN
          DMOD=1.
          WRITE(10,*)
          WRITE(10,'(''NORM OF STRAIN-RATE PROBES IS HARDWIRED TO'',
     #               E11.3)') DMOD
        ENDIF
        IF(.NOT.STRAINR_PROBE) THEN
          SMOD=13.              ! quick & dirty fixed scaling
          SMOD=2.5*CRSS(1,1)    ! system #1 in grain #1 --> quick & dirty scaling
          WRITE(10,'(''NORM OF STRESS PROBES HARDWIRED TO'',
     #               E11.3)') SMOD
        ENDIF

        WRITE(10,*)
        WRITE(10,'(''WILL CALCULATE A'',2I4,''   PCYS PROJECTION'')')
     #             INDX,INDY
        WRITE(10,*) 'INDICES OF THE COMPONENTS DEFINING STRESS SPACE'
        WRITE(10,*) '   1 = (S22-S11)/SQRT(2)'
        WRITE(10,*) '   2 = (2*S33-S22-S11)/SQRT(6)'
        WRITE(10,*) '   3 = S23*SQRT(2)'
        WRITE(10,*) '   4 = S13*SQRT(2)'
        WRITE(10,*) '   5 = S12*SQRT(2)'
        WRITE(10,*)

        IF(INDX.EQ.INDY) THEN
          WRITE(*,*) 'PCYS INDICES HAVE TO BE DIFFERENT'
          STOP
        ELSE IF(INDX.GT.5 .OR. INDY.GT.5) THEN
          WRITE(*,*) 'PCYS INDICES CANNOT BE LARGER THAN 5'
          STOP
        ELSE IF(INDX.GT.INDY) THEN
          IAUX=INDX
          INDX=INDY
          INDY=IAUX
        ENDIF

        NPART=NPROBE/4
        DANG =PI/2.D0/NPART  
ccc        DANG2=DANG/50.        ! APR2020: activate for fine probing

C *********************************************************************
C     GENERATES EQUISPACED UNIT STRAIN-RATE VECTORS (b-BASIS REPRESENT).
C        NPART= # OF PARTITIONS IN EACH INTERVAL [0,PI/2]'
C        FOR NON-CENTRO SYMMETRIC PROPERTIES (ICS=0) SCANS -PI<THETA<PI
C        FOR     CENTRO SYMMETRIC PROPERTIES (ICS=1) SCANS   0<THETA<PI
C        WHEN THETA(1,i)=THETA(2,i)=PI/2  --> UNITDIR(i)=0
C *********************************************************************

      DO I=2,5
        ITHETA(1,I)=NPART
        ITHETA(2,I)=NPART
      ENDDO
      ITHETA(1,INDX)=0
      ITHETA(2,INDX)=0
      ITHETA(1,INDY)=2*NPART*(-1+ICS)
      ITHETA(2,INDY)=2*NPART-1

      NPROBEX=0

      DO ITH5=ITHETA(1,5),ITHETA(2,5)
        THETA5=ITH5*DANG
        COSTH5=COS(THETA5)
        SINTH5=SIN(THETA5)

      DO ITH4=ITHETA(1,4),ITHETA(2,4)
        THETA4=ITH4*DANG
        COSTH4=COS(THETA4)
        SINTH4=SIN(THETA4)

      DO ITH3=ITHETA(1,3),ITHETA(2,3)
        THETA3=ITH3*DANG
        COSTH3=COS(THETA3)
        SINTH3=SIN(THETA3)

      DO ITH2=ITHETA(1,2),ITHETA(2,2)
        THETA2=ITH2*DANG
ccc        THETA2=ITH2*DANG2       ! APR2020: activate for fine probing
        COSTH2=COS(THETA2)
        SINTH2=SIN(THETA2)

        IF(INDX.EQ.1) COSTH1=1.
        IF(INDX.NE.1) COSTH1=0.
	  
        NPROBEX=NPROBEX+1

C *** CALCULATES NORMALIZED VECTOR COMPONENTS IN 5-D SPACE.

        UNITDIR(1)= COSTH1*SINTH2*SINTH3*SINTH4*SINTH5
        UNITDIR(2)=        COSTH2*SINTH3*SINTH4*SINTH5
        UNITDIR(3)=               COSTH3*SINTH4*SINTH5
        UNITDIR(4)=                      COSTH4*SINTH5
        UNITDIR(5)=                             COSTH5

        DO I=1,5
          STRAINR(I,NPROBEX)=DMOD*UNITDIR(I)
          STRESS (I,NPROBEX)=SMOD*UNITDIR(I)
        ENDDO

      ENDDO      ! end of DO ITH2
      ENDDO      ! end of DO ITH3
      ENDDO      ! end of DO ITH4
      ENDDO      ! end of DO ITH5

C *** HARDWIRED BC RUN OPTIONS: FULLY IMPOSED STRAIN RATE OR STRESS PROBES

      IF(STRAINR_PROBE) THEN
        STRAIN_CONTROL=1
        ILBAR(:,:)=1
        SBARc(:,:)=0.
        IDBARv(:)=1
        ISBARv(:)=0
      ENDIF

      IF(.NOT. STRAINR_PROBE) THEN
        STRAIN_CONTROL=0
        ILBAR(:,:)=0
        DBARc(:,:)=0.
        IDBARv(:)=0
        ISBARv(:)=1
      ENDIF

      NSTEPS=NPROBEX

      ENDIF      ! END OF IOPTION=0 / DEFINES EQUISPACED RATE VECTORS & BC's

C ***************************************************************************
C     DEFINE THE STRAIN-RATE OR STRESS PROBE ('DBAR' OR 'SBAR') FOR THE STEP

      IF(IOPTION.EQ.1) THEN

      IF(STRAINR_PROBE) THEN
        DBAR(1:5)=STRAINR(1:5,ISTEP)
        DBAR(6)  =0.
        CALL CHG_BASIS(DBAR,DBARc,AUX55,AUX3333,1,5)
      ENDIF
      write(10,'('' step & DBAR in PCYS'',i5,5f10.4)') istep,dbar(1:5)
      write( *,'('' step & DBAR in PCYS'',i5,5f10.4)') istep,dbar(1:5)
      IF(.NOT. STRAINR_PROBE) THEN
        SBAR(1:5)=STRESS(1:5,ISTEP)
        SBAR(6)  =0.
        CALL CHG_BASIS(SBAR,SBARc,AUX55,AUX3333,1,5)
      ENDIF

      ENDIF      ! END OF IOPTION=1

C ***************************************************************************
C *** STORES THE RESPONSE 'SAV' TO IMPOSED STRAIN-RATE 'DBAR'.
C *** POINTS OF THE PCYS ARE NORMALIZED TO THE SAME DISSIPATION RATE
C     (WRATEREF) ARBITRARILY TAKEN AS WRATEREF=1.
C *** SUPPRESSED NORMALIZING THE STRESS TO THE FIRST POINT IN THE SECTION
C     BECAUSE IT PREVENTS COMPARING LOCI OF DIFFERENT SECTIONS.
C *** A POSSIBILITY IS TO NORMALIZE TO THE INITIAL CRSS OF THE FIRST SYSTEM
C     IN THE FIRST PHASE: STNORM=TAU(1,0,1)

      IF(IOPTION.EQ.2) THEN

        WRATE1=0.0
        WRATE2=0.0
        WRATEBAR=0.0
        DO I=1,5
          WRATE1=WRATE1+DAV(I)*SAV(I)
          WRATEBAR=WRATEBAR+DBAR(I)*SBAR(I)
        ENDDO

C *** calculating overall plastic rate as average over grains' rate 
        WRATE2=0.0
        KGX=1
        DO IPH=IPHBOT,IPHTOP
        DO KKK=NGR(IPH-1)+1,NGR(IPH)
          DO I=1,5
            WRATE2=WRATE2+SG(I,KKK)*DG(I,KGX)*WGT(KKK)
          ENDDO
          KGX=KGX+1
        ENDDO
        ENDDO

CC        WRITE(IU,'(''DBAR'',I5,6F10.3)') ISTEP,DBAR
CC        WRITE(IU,'(''SBAR'',I5,6F10.3)') ISTEP,SBAR
CC        WRITE(IU,'(''DAV '',I5,6F10.3)') ISTEP,DAV
CC        WRITE(IU,'(''SAV '',I5,6F10.3)') ISTEP,SAV
CC        WRITE(IU,'(I5,''   WRATEBAR'',E12.4,''     WRATE1'',E12.4,
CC     #          ''     WRATE2'',E12.4)') ISTEP,WRATEBAR,WRATE1,WRATE2

        WRATEREF=1.
        STNORM=1.
C        STNORM=CRSS(1,1)      ! normalize to CRSS of system #1 in grain #1

        NRSX=NRS(1,1)
        IF(NUNIQUE.EQ.1) THEN
          SCALEF=(WRATEREF/WRATE1)**(NRSX/(1.+NRSX))   ! scaling with Sav*Dav
cc          SCALEF=(WRATEREF/WRATE2)**(NRSX/(1.+NRSX))   ! scaling with <Sg*Dg>
        ELSE
          SCALEF=1.
        ENDIF

C *** CARTESIAN COMPONENTS OF RATE AND STRESS DEFINING PCYS PROJECTION
        DPCYS(1,ISTEP)= DAV(INDX) *SCALEF
        DPCYS(2,ISTEP)= DAV(INDY) *SCALEF
        SPCYS(1,ISTEP)= SAV(INDX) *SCALEF**(1./NRSX) /STNORM
        SPCYS(2,ISTEP)= SAV(INDY) *SCALEF**(1./NRSX) /STNORM

C *** STORES OPPOSITE STATE FOR CENTRO-SYMMETRIC YIELD SURFACE.
        IF(ICS.EQ.1) THEN
          DPCYS(1,ISTEP+NPROBEX)=-DPCYS(1,ISTEP)
          DPCYS(2,ISTEP+NPROBEX)=-DPCYS(2,ISTEP)
          SPCYS(1,ISTEP+NPROBEX)=-SPCYS(1,ISTEP)
          SPCYS(2,ISTEP+NPROBEX)=-SPCYS(2,ISTEP)
        ENDIF

C *** CALCULATES STRESS AND RATE VECTORS LENGTH & ANGLE: RADS,ANGS,RADD,ANGD     
        IF(ISTEP.EQ.NSTEPS) THEN

          DO J=1,NPROBE
            RADS(J)=SQRT(SPCYS(1,J)**2+SPCYS(2,J)**2)
            ANGS(J)=ATAN2(SPCYS(2,J),SPCYS(1,J))
            RADD(J)=SQRT(DPCYS(1,J)**2+DPCYS(2,J)**2)
            ANGD(J)=ATAN2(DPCYS(2,J),DPCYS(1,J))     
            TICK(:,J)=DPCYS(:,J)/RADD(J)*RADS(J)/8.

C *** REDEFINE STRAIN RATE ANGLE TO MAKE IT CONSISTENT WITH STRESS VECTOR ANGLE
            IF(SPCYS(1,J).LT.0.) THEN
              ANGX= 2.*PI-ABS(ANGD(J))
              IF(SPCYS(2,J).GT.0. .AND. ANGD(J).LT.0.) THEN
                ANGD(J)= ANGX
              ELSEIF(SPCYS(2,J).LT.0. .AND. ANGD(J).GT.0.) THEN
                ANGD(J)=-ANGX
              ENDIF
            ENDIF
          ENDDO

C *** WRITES STRESS COORDINATES OF THE PCYS AND OF THE STRAIN RATE VECTOR ('TICK')
C     FOR ADDING TO THE PLOT AND CHECK NORMALITY RULE

          DO J=1,NPROBE
            WRITE(141,'(2E12.4)') SPCYS(:,J)
            WRITE(141,'(2E12.4)') SPCYS(:,J)+TICK(:,J)
            WRITE(141,'(2E12.4)') SPCYS(:,J)
          ENDDO
            WRITE(141,'(2E12.4)') SPCYS(:,1)     ! closes the yield surface plot
            WRITE(141,'(''   -   -   '')')

        ENDIF      ! END OF IF(ISTEP.EQ.NSTEPS)

      ENDIF      ! END OF IOPTION=2

C ***************************************************************************
C *** THE FOLLOWING SECTION IS NOT REQUIRED FOR PCYS CALCULATION AND NEEDS
C     MORE RESEARCH INTO IT.
C *** FITS A QUADRATIC AT EACH POINT OF THE PCYS (see CT notes of 2021 07 16) 
C     TESTS IF CONVEXITY IS VIOLATED.
C     CALCULATES THE NORMAL AT EACH POINT FOR COMPARING AGAINST THE STRAIN  
C     RATE GIVEN BY VPSC.
C ***************************************************************************

      IF (IOPTION.EQ.2 .AND. ISTEP.EQ.NSTEPS) THEN

C *** DEFINES EXTRAPOLATED POINTS AT EACH END FOR FITTING A QUADRATIC: S(ANGS)
C *** LOOPs BETWEEN +NEIX & -NEIX NEIGHBORS OF EACH PCYS POINT AND FITS A QUADRATIC.

        NEIX=1      ! number of neighbors used. must not exceed dimension.
        DO K=0,0+NEIX
          SPCYS(:,-K)=SPCYS(:,NPROBE-K)     
          RADS(-K)   =RADS(NPROBE-K)
          ANGS(-K)   =ANGS(NPROBE-K)
          DPCYS(:,-K)=DPCYS(:,NPROBE-K)
          ANGD(-K)   =ANGD(NPROBE-K)
        ENDDO
        DO K=1,1+NEIX
          SPCYS(:,NPROBE+K)=SPCYS(:,K)     
          RADS(NPROBE+K)   =RADS(K)
          ANGS(NPROBE+K)   =ANGS(K)
          DPCYS(:,NPROBE+K)=DPCYS(:,K)
          ANGD(NPROBE+K)   =ANGD(K)
        ENDDO

        DO NPOINT=1,NPROBE      ! goes over each PCYS point

          XCOEF(:,:)=0.
          XVECT(:)  =0.
          SUMTH(:)  =0.
		  XANG=ANGS(NPOINT)
          DO NEIGX=NPOINT-NEIX,NPOINT+NEIX     ! goes over +NEIX & -NEIX neighbors of point
            XVECT(1)=XVECT(1)+RADS(NEIGX)
            XVECT(2)=XVECT(2)+RADS(NEIGX)*ANGS(NEIGX)
            XVECT(3)=XVECT(3)+RADS(NEIGX)*ANGS(NEIGX)**2
            SUMTH(0)=1+2*NEIX
            DO NPOW=1,4
              SUMTH(NPOW)=SUMTH(NPOW)+ANGS(NEIGX)**NPOW
            ENDDO
          ENDDO
          DO I=1,3
          DO J=1,3
            XCOEF(I,J)=SUMTH(I+J-2)
          ENDDO
          ENDDO
          CALL LU_INVERSE(XCOEF,3)
          DO I=1,3
            XSOL(I)=0.
          DO J=1,3
            XSOL(I)=XSOL(I)+XCOEF(I,J)*XVECT(J)
          ENDDO
          ENDDO
          XFACTOR=XSOL(2)+2.*XSOL(3)*XANG
          YSTANG(1)=-RADS(NPOINT)*SIN(XANG)+XFACTOR*COS(XANG)
          YSTANG(2)= RADS(NPOINT)*COS(XANG)+XFACTOR*SIN(XANG)
          YSNORM(1)= YSTANG(2)
          YSNORM(2)=-YSTANG(1)
		  TICK(:,NPOINT)=YSNORM(:)/SQRT(YSNORM(1)**2+YSNORM(2)**2)
		  ANGN(NPOINT)=ATAN2(YSNORM(2),YSNORM(1))

C *** REDEFINE NORMAL ANGLE TO MAKE IT CONSISTENT WITH STRESS VECTOR ANGLE
          IF(SPCYS(1,NPOINT).LT.0.) THEN
            ANGX= 2.*PI-ABS(ANGN(NPOINT))
            IF(SPCYS(2,NPOINT).GT.0. .AND. ANGN(NPOINT).LT.0.) THEN
              ANGN(NPOINT)= ANGX
            ELSEIF(SPCYS(2,NPOINT).LT.0. .AND. ANGN(NPOINT).GT.0.) THEN
              ANGN(NPOINT)=-ANGX
            ENDIF
          ENDIF

C *** writes 20 points of the quadratic associated with each PCYS point 
      iskip=1
	  if (iskip.eq.0) then 
          WRITE(142,'(''   -   -   '')')
          DO IQUAD=1,20
		   XANG=ANGS(NPOINT-2)+(ANGS(NPOINT+2)-ANGS(NPOINT-2))/20.*IQUAD
		   RFIT=XSOL(1)+XSOL(2)*XANG+XSOL(3)*XANG**2
		   XFIT=RFIT*COS(XANG)
		   YFIT=RFIT*SIN(XANG)
		   WRITE(142,'(2F12.3)') XFIT,YFIT
		  ENDDO
      endif
  
        ENDDO     ! end of DO NPOINT=1,NPROBE

C *** WRITE POINTS OF PCYS AND NORMAL VECTOR FROM QUAD FITTING
        DO J=1,NPROBE
          TICK(:,J)=TICK(:,J) *RADS(J)/8.
          WRITE(142,'(2E12.4)') SPCYS(:,J)
          WRITE(142,'(2E12.4)') SPCYS(:,J)+TICK(:,J)
          WRITE(142,'(2E12.4)') SPCYS(:,J)
        ENDDO
          WRITE(142,'(2E12.4)') SPCYS(:,1)     ! closes the yield surface plot
          WRITE(142,'(''   -   -   '')')

C *** WRITES 2D STRESS AND STRAIN RATE STATES DEFINING THE PCYS
C *** ADDS FIRST STATE AT THE END TO 'CLOSE' THE PCYS FOR GRAPHING

        WRITE(IU,'(10X,''S'',I1,10X,''S'',I1,10X''D'',I1,10X,''D'',I1,
     #    6X,''angS'',6X,''angD'',6X,''angN'')') INDX,INDY,INDX,INDY
        WRITE(IU,'(4E12.4,3F10.2)') 
     #    (SPCYS(1,J),SPCYS(2,J),DPCYS(1,J),DPCYS(2,J),
     #     ANGS(J)*180./PI,ANGD(J)*180./PI,ANGN(J)*180./PI,J=1,NPROBE)       
        WRITE(IU,'(4E12.4,3F10.2)')             ! closes the yield surface plot   
     #     SPCYS(1,1),SPCYS(2,1),DPCYS(1,1),DPCYS(2,1),
     #     ANGS(1)*180./PI,ANGD(1)*180./PI,ANGN(1)*180./PI  

      ENDIF      ! END OF IOPTION=2
C **************************************************************************

      RETURN
      END

C **************************************************************************
C     SUBROUTINE PCYS_IT      --->      VERSION 26/NOV/2022
C
C     IOPTION=0: GENERATES AND STORES EQUISPACED STRAIN-RATE OR STRESS
C                VECTORS IN 5D DEVIATORIC SPACE TO PROBE THE
C                POLYCRYSTAL YIELD SURFACE.
C     IOPTION=1: DEFINES 'DBAR' OR 'SBAR' TO BE IMPOSED AT EACH PROBING STEP.
C     IOPTION=2: STORES RATE OT STRESS (REACTION) IN INTERPOLATION TABLE.
C **************************************************************************

      SUBROUTINE PCYS_IT (ISTEP,ISKIP,IOPTION)

      USE VPSC8DIM
      USE PCYSTABLE

      DIMENSION ITHETA(2,4)
 	  LOGICAL   STRAINR_PROBE

      IU=14      ! UNIT FOR WRITING OUTPUT INTO 'PCYS.OUT'

      STRAINR_PROBE=.TRUE.      ! strain rate probing --> hardwired
C      STRAINR_PROBE=.FALSE.    ! activate for stress probing 

      IF(STRAINR_PROBE)      ITABLE=0 
      IF(.NOT.STRAINR_PROBE) ITABLE=1 

C ********************************************************************
      IF(IOPTION.EQ.0) THEN

        WRITE(10,*)
        WRITE(10,*) '******* STARTS A PCYS CALCULATION *******'
        WRITE(10,*)
        IF(ICS.EQ.1) WRITE(10,*) 'CENTRO-SYMMETRIC YIELD SURFACE'
        IF(ICS.EQ.0) WRITE(10,*) 'NON CENTRO-SYMMETRIC YIELD SURFACE'

C *** HARDWIRES THE ANGULAR INTERVAL USED FOR PROBING. ACTIVATE READ FOR CHOOSING.
        NPART=4

C        WRITE(*,*)
C        WRITE(*,*) 'ENTER PI/2 PARTITION npart'
C        WRITE(*,*) '  1,2,3,4,5,6  for  90,45,30,22.5,18,15 deg --> '
C        READ (*,*) NPART
        IF(NPART.GT.6) THEN
          WRITE(*,*) 'INCREASE DIMENSION OF ACTION & REACTION ARRAYS !'
          STOP
        ENDIF

        DANG=PI/2.D0/FLOAT(NPART)

        WRITE(IU,'('' PCYS INTERPOLATION TABLE'')')
        WRITE(IU,'('' CRYSTAL FILE USED: '',A)') FILECRYS
        WRITE(IU,'('' TEXTURE FILE USED: '',A)') FILETEXT
        WRITE(IU,'('' INTERACTION USED : '',I2)') INTERACTION
        IF (ITABLE.EQ.0)
     #    WRITE(IU,'('' INTERP TABLE BASED ON STRAIN RATE PROBES'')')
        IF (ITABLE.EQ.1)
     #    WRITE(IU,'('' INTERP TABLE BASED ON STRESS PROBES'')')
        WRITE(IU,'(I5)') ITABLE
        WRITE(IU,'('' PARTITION OF PI/2 INTERVAL'')')
        WRITE(IU,'(I5)') NPART
        WRITE(IU,'('' N POWER IN RATE SENSITIVE LAW'')')
        WRITE(IU,'(I5)') NRS(1,1)

C *********************************************************************
C     GENERATES EQUISPACED UNIT STRAIN-RATE VECTORS (b-BASIS REPRESENT).
C        NPART= # OF PARTITIONS IN EACH INTERVAL [0,PI/2]'
C        FOR NON-CENTRO SYMMETRIC PROPERTIES (ICS=0) SCANS -PI<THETA<PI
C        FOR     CENTRO SYMMETRIC PROPERTIES (ICS=1) SCANS   0<THETA<PI
C *********************************************************************

      DO I=1,4
        ITHETA(1,I)=0         ! use -NPART to avoid ortho symmetry ?
        ITHETA(2,I)=NPART
      ENDDO
      ITHETA(2,1)=2*NPART

      NPROBE=0

      DO ITH4=ITHETA(1,4),ITHETA(2,4)
        THETA4=ITH4*DANG
        COSTH4=COS(THETA4)
        SINTH4=SIN(THETA4)

      DO ITH3=ITHETA(1,3),ITHETA(2,3)
        THETA3=ITH3*DANG
        COSTH3=COS(THETA3)
        SINTH3=SIN(THETA3)

      DO ITH2=ITHETA(1,2),ITHETA(2,2)
        THETA2=ITH2*DANG
        COSTH2=COS(THETA2)
        SINTH2=SIN(THETA2)

      DO ITH1=ITHETA(1,1),ITHETA(2,1)
        THETA1=ITH1*DANG
        COSTH1=COS(THETA1)
        SINTH1=SIN(THETA1)

        NPROBE=NPROBE+1
        IF(NPROBE.GT.NPRBMX) THEN
          WRITE(*,*) 'MAXIMUM PROBING DIMENSION EXCEEDED'
          WRITE(*,*) 'NPROBE=', NPROBE
          STOP
        ENDIF

C *** CALCULATES NORMALIZED VECTOR COMPONENTS IN 5-D SPACE.

        ACTION(ITH1,ITH2,ITH3,ITH4,1)= SINTH1*SINTH2*SINTH3*SINTH4
        ACTION(ITH1,ITH2,ITH3,ITH4,2)= COSTH1*SINTH2*SINTH3*SINTH4
        ACTION(ITH1,ITH2,ITH3,ITH4,3)=        COSTH2*SINTH3*SINTH4
        ACTION(ITH1,ITH2,ITH3,ITH4,4)=               COSTH3*SINTH4
        ACTION(ITH1,ITH2,ITH3,ITH4,5)=                      COSTH4

        THGRID(1,NPROBE)=THETA1*180./PI
        THGRID(2,NPROBE)=THETA2*180./PI
        THGRID(3,NPROBE)=THETA3*180./PI
        THGRID(4,NPROBE)=THETA4*180./PI

        ITHGRID(1,NPROBE)=ITH1
        ITHGRID(2,NPROBE)=ITH2
        ITHGRID(3,NPROBE)=ITH3
        ITHGRID(4,NPROBE)=ITH4

C       WRITE(IU,'('' ACTION'',I3,3X,5F10.5)')
C    #      NPROBE,(ACTION(ITH1,ITH2,ITH3,ITH4,i),i=1,5)

      ENDDO      ! end of DO ITH1
      ENDDO      ! end of DO ITH2
      ENDDO      ! end of DO ITH3
      ENDDO      ! end of DO ITH4

C *** HARDWIRE RUN PARAMETERS

      IF(ITABLE.EQ.0) THEN      ! PROBES PCYS USING STRAIN RATES
        STRAIN_CONTROL=1
        DO I=1,3
        DO J=1,3
          ILBAR(I,J)=1
          SBARc(I,J)=0.
        ENDDO
        ENDDO
        DO I=1,6
          IDBARv(I)=1
          ISBARv(I)=0
        ENDDO
      ELSE IF(ITABLE.EQ.1) THEN      ! PROBES PCYS USING STRESS
        STRAIN_CONTROL=0
        DO I=1,3
        DO J=1,3
          ILBAR(I,J)=0
        ENDDO
        ENDDO
        ILBAR(2,1)=1
        ILBAR(3,1)=1
        ILBAR(3,2)=1
        DO I=1,6
          IDBARv(I)=0
          ISBARv(I)=1
        ENDDO
      ENDIF

      NSTEPS=NPROBE

      WRITE(IU,'('' NUMBER OF PROBE STATES USED'')')
      WRITE(IU,'(I7)') NPROBE
      WRITE(IU,*)
      IF(ITABLE.EQ.0) WRITE(IU,'(''  THE1  THE2  THE3  THE4''
     #      ,11X,''S1'',11X,''S2'',11X,''S3'',11X,''S4'',11X,''S5'')')
      IF(ITABLE.EQ.1) WRITE(IU,'(''  THE1  THE2  THE3  THE4''
     #      ,11X,''D1'',11X,''D2'',11X,''D3'',11X,''D4'',11X,''D5'')')

      ENDIF      ! END OF IOPTION=0 TO GENERATE EQUISPACED UNIT VECTORS

C ***************************************************************************
      IF(IOPTION.EQ.1) THEN

C *** DEFINE THE PROBE VECTOR 'DBAR' OR 'SBAR' FOR THE STEP
        IF(ISTEP.EQ.1) NPROBERED=0

        ITH1=ITHGRID(1,ISTEP)
        ITH2=ITHGRID(2,ISTEP)
        ITH3=ITHGRID(3,ISTEP)
        ITH4=ITHGRID(4,ISTEP)

C *** AVOID RECALCULATING FOR REDUNDANT COMBINATIONS
C     (i.e. IF THETA(4)=0 THEN ACTION1=ACTION2=ACTION3=ACTION4=0,
C     INDEPENDENT OF THE VALUE OF ITH1,ITH2,ITH3)

      ISKIP=0
      IF(ITH4.EQ.0 .AND. (ITH3+ITH2+ITH1).NE.0) THEN
        DO I=1,5
          REACTION(ITH1,ITH2,ITH3,ITH4,I)=REACTION(0,0,0,0,I)
        ENDDO
        ISKIP=1
      ELSE IF(ITH3.EQ.0 .AND. (ITH2+ITH1).NE.0) THEN
        DO I=1,5
          REACTION(ITH1,ITH2,ITH3,ITH4,I)=REACTION(0,0,0,ITH4,I)
        ENDDO
        ISKIP=1
      ELSE IF(ITH2.EQ.0 .AND. (ITH1).NE.0) THEN
        DO I=1,5
          REACTION(ITH1,ITH2,ITH3,ITH4,I)=REACTION(0,0,ITH3,ITH4,I)
        ENDDO
        ISKIP=1
      ENDIF

      IF(ISKIP.EQ.0) THEN
        NPROBERED=NPROBERED+1
        IF(ITABLE.EQ.0) THEN
          DO I=1,5
            DBAR(I)=ACTION(ITH1,ITH2,ITH3,ITH4,I)
          ENDDO
          CALL CHG_BASIS(DBAR,DBARc,AUX55,AUX3333,1,5)
        ELSE IF(ITABLE.EQ.1) THEN
          DO I=1,5
            SBAR(I)=ACTION(ITH1,ITH2,ITH3,ITH4,I)
          ENDDO
          CALL CHG_BASIS(SBAR,SBARc,AUX55,AUX3333,1,5)
        ENDIF
      ENDIF

      IF(ISKIP.EQ.1) THEN
        WRITE(IU,'(4F6.1,5E13.5)') (THGRID(I,ISTEP),I=1,4),
     #               (REACTION(ITH1,ITH2,ITH3,ITH4,I),I=1,5)

        IF(ISTEP.EQ.NSTEPS) THEN
          WRITE(IU,*)
          WRITE(IU,'('' NPROBERED ='',I7)') NPROBERED
        ENDIF
      ENDIF

      ENDIF      ! END OF IOPTION=1
C ***************************************************************************
      IF(IOPTION.EQ.2) THEN

C     STORES THE REACTION 'SAV' TO STRAIN-RATE 'DBAR' (IF ITABLE=0)
C     OR THE REACTION 'DBAR TO THE STRESS 'SBAR' (IF ITABLE=1)

        ITH1=ITHGRID(1,ISTEP)
        ITH2=ITHGRID(2,ISTEP)
        ITH3=ITHGRID(3,ISTEP)
        ITH4=ITHGRID(4,ISTEP)

        IF(ITABLE.EQ.0) THEN
          SNORM=SBAR(1)**2+SBAR(2)**2+SBAR(3)**2+SBAR(4)**2+SBAR(5)**2
          SNORM=SQRT(SNORM)
          DO I=1,5
            IF(ABS(SBAR(I)/SNORM) .LE. 1.E-5) SBAR(I)=0.
            REACTION(ITH1,ITH2,ITH3,ITH4,I)=SBAR(I)
          ENDDO
        ELSE IF(ITABLE.EQ.1) THEN
          DNORM=DBAR(1)**2+DBAR(2)**2+DBAR(3)**2+DBAR(4)**2+DBAR(5)**2
          DNORM=SQRT(DNORM)
          DO I=1,5
            IF(ABS(DBAR(I)/DNORM) .LE. 1.E-5) DBAR(I)=0.
            REACTION(ITH1,ITH2,ITH3,ITH4,I)=DBAR(I)
          ENDDO
        ENDIF

        WRITE(IU,'(4F6.1,5E13.5)') (THGRID(I,ISTEP),I=1,4),
     #               (REACTION(ITH1,ITH2,ITH3,ITH4,I),I=1,5)

        IF(ISTEP.EQ.NSTEPS) THEN
          WRITE(IU,*)
          WRITE(IU,'('' NPROBERED ='',I7)') NPROBERED
        ENDIF

      ENDIF      ! END OF IOPTION=2
C **********************************************************************

      RETURN
      END

C ***********************************************************************
C     SUBROUTINE POSTMORTEM     --->      VERSION 23/MAR/2007
C
C     SAVES ARRAYS ASSOCIATED WITH A GIVEN STATE OF DEFORMATION IN
C     GRAINS AND POLYCRYSTAL. PERMITS TO START RUN FROM SUCH STATE.
C     THIS VERSION IS FOR A SINGLE ELEMENT RUN.
C     THIS VERSION READS AND WRITES A BINARY FILE.
C
C     IF IOPTION=1 READS GRAIN AND PX STATES FROM 'POSTMORT.IN'
C     IF IOPTION=2 SAVES GRAIN AND PX STATES INTO 'POSTMORT.OUT'
C
C     IF ISAVE=isaveprevious=1: SAVES MACRO TENSORS AND INITIAL GRAIN STRESS.
C        MEANT FOR RUNNING ONLY ONE STEP IN ORDER TO GET QUICK
C        INITIAL CONVERGENCE OR RESTART WITH AN OPEN SCYS.
C     IF ISAVE=isaveprevious=ISTEP: SAVES FULL GRAIN STATES ASSOCIATED WITH
C        ISTEP AND ALLOWS TO CONTINUE DEFORMATION FROM A PREVIOUS RUN.
C ************************************************************************
C
      SUBROUTINE POSTMORTEM (ioption)

      USE VPSC8DIM

      ngr1=1
      ngr2=ngr(nph)

      IF (IOPTION.EQ.1) THEN

        read(UR4) isaveprevious

        read(UR4) (DAV(i),i=1,5)
        read(UR4) (SAV(i),i=1,5)
        read(UR4) (DBAR0(i),i=1,5)
        read(UR4) ((MBARTG(i,j),j=1,5),i=1,5)
        read(UR4) ((LBARTG(i,j),j=1,5),i=1,5)
        read(UR4) ((SG(i,KKK),i=1,5),KKK=ngr1,ngr2)

        if(isaveprevious.eq.1) return

        read(UR4) ((EPSTOTc(i,j),i=1,3),j=1,3),epsvm,epsacu
        read(UR4) ((axisph(i,j,0),i=0,3),j=1,3)
        read(UR4) (( fijph(i,j,0),i=1,3),j=1,3)

        do iph=1,nph
          ngr1=ngr(iph-1)+1
          ngr2=ngr(iph)
          read(UR4) ((axisph(i,j,iph),i=0,3),j=1,3)
          read(UR4) (( fijph(i,j,iph),i=1,3),j=1,3)

          read(UR4) (((ag(i,j,kkk),i=1,3),j=1,3),KKK=ngr1,ngr2)
          read(UR4) (wgt(kkk),KKK=ngr1,ngr2)
          read(UR4) (gtotgr(kkk),KKK=ngr1,ngr2)
          if(ihardlaw.eq.0) read(UR4) ((crss(is,kkk),is=1,nsyst(iph)),
     #                  KKK=ngr1,ngr2)
          if(ihardlaw.eq.1) read(UR4) ((taue(is,kkk),is=1,nsyst(iph)),
     #                  KKK=ngr1,ngr2)

c   read twin-related arrays for each phase if appropriate.
          if(ntwmod(iph).ne.0) then
            read(UR4)  pritw(iph),sectw(iph)
            read(UR4) (ntwevents(kkk),KKK=ngr1,ngr2)
            read(UR4) (twfrph(itwm,iph),itwm=1,ntwmod(iph))
            read(UR4) (eftwfr(itwm,iph),itwm=1,ntwmod(iph))
            read(UR4) (ktwsmx(kkk),KKK=ngr1,ngr2)
            read(UR4) ((twfrsy(is,kkk),is=1,nsyst(iph)),KKK=ngr1,ngr2)
          endif
        enddo

      ENDIF

      IF (IOPTION.EQ.2) THEN

        write(UW2) isave

        write(UW2) (DAV(i),i=1,5)
        write(UW2) (SAV(i),i=1,5)
        write(UW2) (DBAR0(i),i=1,5)
        write(UW2) ((MBARTG(i,j),j=1,5),i=1,5)
        write(UW2) ((LBARTG(i,j),j=1,5),i=1,5)
        write(UW2) ((SG(i,KKK),i=1,5),KKK=ngr1,ngr2)

        if(isave.eq.1) return

        write(UW2) ((EPSTOTc(i,j),i=1,3),j=1,3),epsvm,epsacu
        write(UW2) ((axisph(i,j,0),i=0,3),j=1,3)
        write(UW2) (( fijph(i,j,0),i=1,3),j=1,3)

        do iph=1,nph
          ngr1=ngr(iph-1)+1
          ngr2=ngr(iph)
          write(UW2) ((axisph(i,j,iph),i=0,3),j=1,3)
          write(UW2) (( fijph(i,j,iph),i=1,3),j=1,3)

          write(UW2) (((ag(i,j,kkk),i=1,3),j=1,3),KKK=ngr1,ngr2)
          write(UW2) (wgt(kkk),KKK=ngr1,ngr2)
          write(UW2) (gtotgr(kkk),KKK=ngr1,ngr2)
          if(ihardlaw.eq.0) write(UW2) ((crss(is,kkk),is=1,nsyst(iph)),
     #                  KKK=ngr1,ngr2)
          if(ihardlaw.eq.1) write(UW2) ((taue(is,kkk),is=1,nsyst(iph)),
     #                  KKK=ngr1,ngr2)

c   read twin-related arrays for each phase if appropriate.
          if(ntwmod(iph).ne.0) then
            write(UW2)  pritw(iph),sectw(iph)
            write(UW2) (ntwevents(kkk),KKK=ngr1,ngr2)
            write(UW2) (twfrph(itwm,iph),itwm=1,ntwmod(iph))
            write(UW2) (eftwfr(itwm,iph),itwm=1,ntwmod(iph))
            write(UW2) (ktwsmx(kkk),KKK=ngr1,ngr2)
            write(UW2) ((twfrsy(is,kkk),is=1,nsyst(iph)),KKK=ngr1,ngr2)
          endif
        enddo

      ENDIF

      RETURN
      END

C ******************************************************************************
C     SUBROUTINE RODRIGUES    --->   VERSION 11/AUG/2020

C     BUILDS INCREMENTAL ROTATION MATRIX 'EXP_OMEGA' BASED ON RODRIGUES FORMULA.
C     'OMEGA'    :INCREMENTAL LATTICE SPIN TENSOR (SKEW-SYMMETRIC)
C             --> OMEGA=WBAR+ROTLOC-skew(LIJ0)/2 
C     'EXP_OMEGA':INCREMENTAL TRANSFORMATION FROM INITIAL TO FINAL ORIENTATION.
C             --> EXP_OMEGA=I+sin(omega)*OMEGA+(1-cos(omega))*OMEGA^2
C ******************************************************************************

      SUBROUTINE RODRIGUES (OMEGA,EXP_OMEGA)

      DIMENSION OMEGA(3,3),OMEGA2(3,3),EXP_OMEGA(3,3),v(3),xid33(3,3)

      xid33=0.
      do i=1,3
        xid33(i,i)=1.
      enddo

c ** v(i) is the Rodrigues spin axis. The norm is the rotation angle omega.
      v(1)=OMEGA(3,2)
      v(2)=OMEGA(1,3)
      v(3)=OMEGA(2,1)
      vnorm=sqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))

      if(vnorm.lt.1.e-06) then
        do i=1,3
        do j=1,3
          EXP_OMEGA(i,j)=xid33(i,j)
        enddo
        enddo
        return
      endif

      coef1=sin(vnorm)/vnorm
      coef2=(1.-cos(vnorm))/vnorm**2

      do i=1,3
      do j=1,3
        OMEGA2(i,j)=0.
        do k=1,3
          OMEGA2(i,j)=OMEGA2(i,j)+OMEGA(i,k)*OMEGA(k,j)
        enddo
      enddo
      enddo

      do i=1,3
      do j=1,3
        EXP_OMEGA(i,j)=xid33(i,j)+coef1*OMEGA(i,j)+coef2*OMEGA2(i,j)
      enddo
      enddo

      return
      end
C
C *************************************************************************
C     SUBROUTINE SCALE_GAMD0S   ---->   VERSION 20/MAY/2015
C
C     THIS SUBROUTINE IS MEANT TO BE USED IN COMBINATION WITH MTS OR DD 
C     MODELS WHERE RATE EFFECTS ARE ACCOUNTED FOR IN THE FUNCTIONAL FORM 
c     OF THE CRSS.
C     GAMD0G (SAME FOR EVERY SYSTEM IN THE GRAIN) WAS REPLACED BY GAMD0S 
C     (POTENTIALLY DIFFERENT FOR EVERY SYSTEM IN THE GRAIN)
C
C     BY REDEFINING THE FACTOR 'GAMD0' TO BE OF THE ORDER OF THE STRAIN
C     RATE THE RATIO (RSS/CRSS) WILL BE OF ORDER ONE AND THERE ARE WEAK 
C     RATE SENSITIVITY EFFECTS ASSOCIATED WITH THE N'th POWER.
C     IN SUCH CASE THE N'th POWER IS ONLY A WAY TO AVOID AMBIGUITIES AND
C     TO GIVE SHEAR RATES IN EVERY SLIP AND TWIN SYSTEM.
C *************************************************************************

      SUBROUTINE SCALE_GAMD0S (ISTEP)

      USE VPSC8DIM

      REFRATE=VNORM(DBAR,5)

      KGX=1
      DO IPH=IPHBOT,IPHTOP
        IPHEL=IPH-IPHBOT+1
        DO KKK=NGR(IPH-1)+1,NGR(IPH)
          DO IS=1,NSYST(IPHEL)
            GAMD0S(IS,KKK)=REFRATE
          ENDDO
          KGX=KGX+1
        ENDDO
      ENDDO

      RETURN
      END
C
C *****************************************************************************
C     SUBROUTINE STATE_NxN      --->      VERSION 28/AUG/2017

C     ex SUBROUTINE STATE_5x5 EXTENDED TO SOLVE ndim x ndim SYSTEMS WITH
C     COMPONENTS EXPRESSED IN b-BASIS

C     BASED ON SECANT EQUATION: D(i)-DBAR0(i)=MS(i,j)*S(j)
C     SHIFTS KNOWN COMPONENTS TO INDEPENDENT TERM, CALCULATES COEFFICIENTS
C     FOR UNKNOWN COMPONENTS AND INVERTS THE LINEAR SYSTEM.
C *** CURRENTLY CALLED AT THE END OF EACH INCREMENTAL STEP IN SUBROUTINE VPSC
C     IF IDBARv1*IDBARv2*IDBARv3=1 (FULLY ENFORCED DIAGONAL)
C *****************************************************************************
C *****************************************************************************
C
C   SUBROUTINES FOR DOING FLUCTUATIONS AND SECOND ORDER (SO)
C   by RICARDO LEBENSOHN
C
C      SUBROUTINE SO_FLUCTUATIONS
C      SUBROUTINE SO_LINSOLVER25
C      SUBROUTINE SO_GET_GAMDOT
C      SUBROUTINE SO_GET_THEFLU
C      SUBROUTINE SO_PROCEDURE
C      SUBROUTINE SO_MOD
C      SUBROUTINE SO_SDPX
C      SUBROUTINE SO_GRAIN_STRESS_ALT
C      SUBROUTINE SO_EXTRAPOL
C      SUBROUTINE SO_VOIGT10
C
C *****************************************************************************
C *****************************************************************************
C
C     SUBROUTINE SO_FLUCTUATIONS      --->     VERSION 4/JAN/07
C
C     CALCULATES DERIVATIVES OF EFFECTIVE MODULI WRT LOCAL COMPLIANCES
C     TO OBATAIN STRESS AND STRAIN-RATE SECOND-ORDER MOMENT TENSORS
C     AND CALLS SOP (SECOND-ORDER PROCEDURE) TO OBTAIN IMPROVED
C     LINEARIZED LOCAL MODULI
C ******************************************************************

      subroutine so_fluctuations(erreso,erraso)

      USE VPSC8DIM

      DIMENSION ALFFLU(5,5,5,5),OMEFLU(5,5,5,5)
      DIMENSION THEFLU(5,5,5,5,NPHPEL)
      DIMENSION GAMFLU1A(5),GAMFLU1B(5),GAMFLU3(5)
      DIMENSION XIAUX1(5),XIAUX2(5),XIAUX3(5)
      DIMENSION PIFLU(5,5),DELAUX(5),DELFLU(5,5,5),XKAPFLU(5)
      DIMENSION DMAUX(5,5)
      DIMENSION DMDMR(5,5,5,5),DEAUX(5),DEDMR(5,5,5)
      DIMENSION DMTILDE(5,5,NPHPEL),DGDMR(5,5)

      dimension asonew(nsysmx),esonew(nsysmx)

      write(*,*)
      WRITE(*,*) 'SUBROUTINE FLUCTUATIONS'
      write(*,*)
c
c     INITIALIZE OMEFLU, DELFLU & XKAPFLU
c
      DO I=1,5
      XKAPFLU(I)=0.
      DO J=1,5
      DO K=1,5
      DELFLU(I,J,K)=0.
      DO L=1,5
c!      OMEFLU(I,J,K,L)=XID5(I,K)*XID5(J,L)
      OMEFLU(I,J,K,L)=0.5*(XID5(I,K)*XID5(J,L)+XID5(I,L)*XID5(J,K))
      ENDDO
      ENDDO
      ENDDO
      ENDDO
c
      KGX=1
      DO IPH=IPHBOT,IPHTOP

      CALL SO_GET_THEFLU(THEFLU,IPH)

      DO KKK=NGR(IPH-1)+1,NGR(IPH)
c
c     ALFFLU
c
      DO I=1,5
      DO J=1,5
      DO K=1,5
      DO L=1,5

      DUMMY=0.5*(XID5(I,K)*XID5(J,L)+XID5(I,L)*XID5(J,K))

      DO M=1,5
      DUMMY=DUMMY+THEFLU(I,M,K,L,IPH)*(XID5(M,J)-BG(M,J,KGX))
      ENDDO
c
      ALFFLU(I,J,K,L)=DUMMY
c
      ENDDO
      ENDDO
      ENDDO
      ENDDO
c
c     OMEFLU
c
      DO I=1,5
      DO J=1,5
      DO K=1,5
      DO L=1,5

      DUMMY=0.

      DO M=1,5
      DUMMY=DUMMY+BETFLU(I,M,KGX)*ALFFLU(M,J,K,L)
      ENDDO

      OMEFLU(I,J,K,L)=OMEFLU(I,J,K,L)-DUMMY*WGT(KKK)

      ENDDO
      ENDDO
      ENDDO
      ENDDO
c
c     DELFLU & XKAPFLU
c
      DO I=1,5
      DUMMY=0.
      DO J=1,5
      DUMMY=DUMMY+DG0(J,KGX)*CHIFLU(J,I,KGX)
      ENDDO
      DELAUX(I)=DUMMY
      XKAPFLU(I)=XKAPFLU(I)+DUMMY*WGT(KKK)
      ENDDO
c
      DO I=1,5
      DO J=1,5
      DO K=1,5
c
      DUMMY=0.
c
      DO M=1,5
      DUMMY=DUMMY+DELAUX(M)*ALFFLU(M,I,J,K)
      ENDDO
c
      DELFLU(I,J,K)=DELFLU(I,J,K)+DUMMY*WGT(KKK)
c
      ENDDO
      ENDDO
      ENDDO
c
      KGX=KGX+1
      ENDDO
      ENDDO
c
      erreso=0.
      erraso=0.

      utilde=0.

      KGX=1
      DO IPH=IPHBOT,IPHTOP
      IPHEL=IPH-IPHBOT+1

      DO KKK=NGR(IPH-1)+1,NGR(IPH)
c
cw      write(*,'(1H+,a,i5)') 'GRAIN ',kgx
cw      write(*,*) 'GRAIN ',kgx
c
      DO IU=1,5
      DO IV=1,5
c
c     dM/dMr(kgx)(iu,iv)
c
c     PIFLU

      DO I=1,5
      DO J=1,5
      PIFLU(I,J)=WGT(KKK)/2.*
     #           ((XID5(I,IU)-BETFLU(I,IU,KGX))*BG(IV,J,KGX)
     #           +(XID5(I,IV)-BETFLU(I,IV,KGX))*BG(IU,J,KGX))
      ENDDO
      ENDDO

      CALL SO_LINSOLVER25(OMEFLU,PIFLU,DMAUX)

      DO I=1,5
      DO J=1,5
      DMDMR(I,J,IU,IV)=DMAUX(I,J)
      ENDDO
      ENDDO
c
c     dM~/dMr TO BE USED BELOW FOR dg~/dMr
c
      do iph2=iphbot,iphtop

      DO I=1,5
      DO J=1,5
      DUMMY=0.
      DO K=1,5
      DO L=1,5
      DUMMY=DUMMY+THEFLU(I,J,K,L,IPH2)*DMAUX(K,L)
      ENDDO
      ENDDO
      DMTILDE(I,J,IPH2)=DUMMY
      ENDDO
      ENDDO

      enddo  ! (iph2)
c
c     dE~/dMr(kgx)(iu,iv)
c
c     GAMFLU1A
c
      GFACT=0.
      DO I=1,5
      GFACT=GFACT+DG0(I,KGX)*CHIFLU(I,IU,KGX)
      ENDDO

      GFACT=-WGT(KKK)/2.*GFACT

      DO J=1,5
      GAMFLU1A(J)=GFACT*BG(IV,J,KGX)
      ENDDO
c
c     GAMFLU1B
c
      GFACT=0.
      DO I=1,5
      GFACT=GFACT+DG0(I,KGX)*CHIFLU(I,IV,KGX)
      ENDDO

      GFACT=-WGT(KKK)/2.*GFACT

      DO J=1,5
      GAMFLU1B(J)=GFACT*BG(IU,J,KGX)
      ENDDO
c
c     GAMFLU3(i)= DELFLU(ijk) DMAUX(jk)
c
      DO I=1,5
      DUMMY=0.
      DO J=1,5
      DO K=1,5
      DUMMY=DUMMY+DELFLU(I,J,K)*DMAUX(J,K)
      ENDDO
      ENDDO
      GAMFLU3(I)=DUMMY
      ENDDO
c
c     dE~/dMr(kgx)(iu,iv)
c
      DO I=1,5
      DEAUX(I)=GAMFLU1A(I)+GAMFLU1B(I)+GAMFLU3(I)
      DEDMR(I,IU,IV)=DEAUX(I)
      ENDDO
c
c     dg~/dMr(kgx)(iu,iv)
c
c     XIFLU1A
c
      XFACT1=0.
      DO L=1,5
      XFACT1=XFACT1+CHIFLU(IV,L,KGX)*(DBAR0(L)-DG0(L,KGX))
      ENDDO

      XFACT2=0.
      DO I=1,5
      XFACT2=XFACT2+DG0(I,KGX)*CHIFLU(I,IU,KGX)
      ENDDO

      XIFLU1A=-WGT(KKK)/2*XFACT1*XFACT2
c!
c     XIFLU1B
c
      XFACT1=0.
      DO L=1,5
      XFACT1=XFACT1+CHIFLU(IU,L,KGX)*(DBAR0(L)-DG0(L,KGX))
      ENDDO

      XFACT2=0.
      DO I=1,5
      XFACT2=XFACT2+DG0(I,KGX)*CHIFLU(I,IV,KGX)
      ENDDO

      XIFLU1B=-WGT(KKK)/2*XFACT1*XFACT2
c
c     XIFLU2
c
      XIFLU2=0.

      KGX2=1
      DO IPH2=IPHBOT,IPHTOP
      DO KKK2=NGR(IPH2-1)+1,NGR(IPH2)

      DO I=1,5
      DUMMY=0.
      DO J=1,5
      DUMMY=DUMMY+CHIFLU(I,J,KGX2)*(DBAR0(J)-DG0(J,KGX2))
      ENDDO
      XIAUX1(I)=DUMMY
      ENDDO

      DO I=1,5
      DUMMY=0.
      DO J=1,5
      DUMMY=DUMMY+DMTILDE(I,J,IPH2)*XIAUX1(J)
      ENDDO
      XIAUX2(I)=DUMMY
      ENDDO

      DO I=1,5
      DUMMY=0.
      DO J=1,5
      DUMMY=DUMMY+CHIFLU(I,J,KGX2)*XIAUX2(J)
      ENDDO
      XIAUX3(I)=DUMMY
      ENDDO

      DUMMY=0.
      DO I=1,5
      DUMMY=DUMMY+DG0(I,KGX2)*XIAUX3(I)
      ENDDO
      XIFLU2=XIFLU2-DUMMY*WGT(KKK2)

      KGX2=KGX2+1
      ENDDO
      ENDDO
c
c     XIFLU3
c
      XIFLU3=0.
      DO I=1,5
      XIFLU3=XIFLU3+XKAPFLU(I)*DEAUX(I)
      ENDDO

      DGAUX=XIFLU1A+XIFLU1B+XIFLU2+XIFLU3
      DGDMR(IU,IV)=DGAUX

      ENDDO      ! IU ENDDO
      ENDDO      ! IV ENDDO
c
c     SECOND MOMENT STRESS TENSOR (kgx)
c
      DO K=1,5
      DO L=1,5
      DERIVM=0.
      DERIVE=0.
      DO I=1,5
      DERIVE=DERIVE+DEDMR(I,K,L)*SBAR(I)
      DO J=1,5
      DERIVM=DERIVM+DMDMR(I,J,K,L)*SBAR(I)*SBAR(J)
      ENDDO
      ENDDO
      SECMOM5(K,L,KGX)=(DERIVM+2*DERIVE+DGDMR(K,L))/WGT(KKK)
      ENDDO
      ENDDO
c
cw      write(*,*)'GRAIN',kgx
cw      write(*,*)
c
cw      write(*,*) 'SECMOM='
cw      do i=1,5
cw      write(*,'(5E12.3)')(secmom5(i,j,kgx),j=1,5)
cw      enddo
cw      pause
c
      if(interaction.eq.5) then

      call so_procedure(iphel,kgx,kkk,asonew,esonew)

      esonorm=0.
      asonorm=0.
      erreso1=0.
      erraso1=0.

      do is=1,nsyst(iphel)
cc
cc      write(94,'(2i5,4e12.3)')
cc     # kgx,is,aso(is,kgx),asonew(is),eso(is,kgx),esonew(is)
cc
cx      if(kgx.eq.1.and.is.eq.1) write(94,'(2e12.3)')
cx     # asonew(is),esonew(is)
c
      esonorm=esonorm+((eso(is,kgx)+esonew(is))/2)**2
      asonorm=asonorm+((aso(is,kgx)+asonew(is))/2)**2
      erreso1=erreso1+(eso(is,kgx)-esonew(is))**2
      erraso1=erraso1+(aso(is,kgx)-asonew(is))**2

      eso(is,kgx)=2./3.*eso(is,kgx)+esonew(is)/3.
      aso(is,kgx)=2./3.*aso(is,kgx)+asonew(is)/3.

      enddo

      if(jxrs.ne.1) then
      erreso=erreso+sqrt(erreso1/esonorm)*wgt(kkk)
      endif
      erraso=erraso+sqrt(erraso1/asonorm)*wgt(kkk)

      endif
c
c     SECOND MOMENT STRAIN-RATE TENSOR (kgx)
c
      DO I=1,5
      DO J=1,5
      DUMMY=0.
      DO K=1,5
      DO L=1,5
      DUMMY=DUMMY+MGTG(I,K,KGX)*MGTG(J,L,KGX)*SECMOM5(K,L,KGX)
      ENDDO
      ENDDO
c
      if(interaction.eq.3) then
c
      factn1=10**2
      factn2=2*(1-10)-(1-10)**2
      SECMOM5D(I,J,KGX)=factn1*DUMMY+factn2*DG(I,KGX)*DG(J,KGX)
c
      else if(interaction.eq.4) then

      factn1=jxrs**2
      factn2=2*(1-jxrs)-(1-jxrs)**2
c
      SECMOM5D(I,J,KGX)=factn1*DUMMY+factn2*DG(I,KGX)*DG(J,KGX)
c
      else
c
      SECMOM5D(I,J,KGX)=DUMMY+DG(I,KGX)*DG0(J,KGX)+
     #   DG(J,KGX)*DG0(I,KGX)-DG0(I,KGX)*DG0(J,KGX)
c
      endif
c
      ENDDO
      ENDDO
c
      KGX=KGX+1
      ENDDO
      ENDDO
c
      if(interaction.eq.5) then
      write(*,*)
      write(*,*) 'ERRESO  = ',erreso
      write(*,*) 'ERRASO  = ',erraso
c      write(*,*) 'UTILDE  = ',utilde
c      write(*,*) 'DVM     = ',dvm
c      write(*,*) 'DVMSO   = ',dvmso
c
c      write(83,*) 'ERRESO = ',erreso
c      write(83,*) 'ERRASO = ',erraso
      endif

      return

      end
C
C *****************************************************************
C     SUBROUTINE LINSOLVER25      --->     VERSION 4/JAN/07
C
C     SOLVER OF A 25X25 LINEAR SYSTEM, FOR FLUCTUATION CALCULATION
C ******************************************************************

      SUBROUTINE SO_LINSOLVER25(A,B,X)
c
      DIMENSION A(5,5,5,5),B(5,5),X(5,5)
      DIMENSION A1(25,25),B1(25)
      DIMENSION IC(5,5)

      K=0
      DO I=1,5
      DO J=1,5
      K=K+1
      IC(I,J)=K
      ENDDO
      ENDDO
c
      DO I=1,5
      DO J=1,5
c
      I1=IC(I,J)
      B1(I1)=B(I,J)
c
      DO K=1,5
      DO L=1,5
c
      J1=IC(K,L)
      A1(I1,J1)=A(I,J,K,L)

      ENDDO
      ENDDO
      ENDDO
      ENDDO
c
c      open(78,file='trace.out',status='unknown')
c      write(78,*) 'A1='
c      do i=1,25
c      write(78,'(25E10.2)') (a1(i,j),j=1,25)
c      enddo
c      write(78,*) 'B1='
c      write(78,'(25E10.2)') (b1(j),j=1,25)

       CALL LU_EQSYSTEM(A1,B1,25,ISINGULAR)

c      write(78,*) 'B1='
c      write(78,'(25E10.2)') (b1(j),j=1,25)
c
      if(isingular.eq.1) then
       write(*,*) 'SINGULAR SYSTEM IN LINSOLVER25'
       stop
      endif
c
      K=0
      DO I=1,5
      DO J=1,5
      K=K+1
      X(I,J)=B1(K)
      ENDDO
      ENDDO
c
      RETURN
      END
C
C *****************************************************************
C     SUBROUTINE GET_GAMDOT     --->     VERSION 26/DEC/2017
C
C     CALCULATES GAMDOT's FOR GIVEN VALUES OF GRAIN STRES SG.
C     THE CALL TO 'GET_GAMDOT' IS RELATED TO THE USE OF 'GRAIN_STRESS_ALT'
C     AND IT IS CARRIED OUT AT THE END OF EACH DEFORMATION STEP ONLY
C     IF INTERACTION=5 (SO CALCULATION).
C     OTHERWISE, THE GAMDOT's ARE CALCULATED IN GRAIN_STRESS
C *****************************************************************

      SUBROUTINE SO_GET_GAMDOT

      USE VPSC8DIM

      DIMENSION SX(5),GAMD0X(NSYSMX)
      DIMENSION SCHX(5,NSYSMX),TAUX(NSYSMX),NRSX(NSYSMX),ISENSEX(NSYSMX)

      KGX=1
      DO IPH=IPHBOT,IPHTOP
        IPHEL=IPH-IPHBOT+1
        DO KKK=NGR(IPH-1)+1,NGR(IPH)

          DO IS=1,NSYST(IPHEL)
cnt            GAMD0X(IS)=GAMD0G(KGX)
            GAMD0X(IS)=GAMD0S(IS,KGX)
            ISENSEX(IS)=ISENSE(IS,IPHEL)
            NRSX(IS)=   NRS(IS,IPHEL)
            TAUX(IS)=   CRSS(IS,KKK)
            DO J=1,5
              SCHX(J,IS)=SCH(J,IS,KGX)
            ENDDO
          ENDDO

          DO I=1,5
            SX(I)=SG(I,KKK)
          ENDDO

C     GETS RESOLVED SHEAR STRESSES 'rssx' AND SHEAR RATES 'gamdot'.
C       SIGN(GAMDOT)=SIGN(RSSX)
C       NRS CAN BE EVEN OR ODD
C       RSS IS ALWAYS > 0 AND IS USED TO CALCULATE VISCOUS COMPLIANCE.

          DO IS=1,NSYST(IPHEL)
            RSSX=SCHX(1,IS)*SX(1)+SCHX(2,IS)*SX(2)+SCHX(3,IS)*SX(3)+
     #           SCHX(4,IS)*SX(4)+SCHX(5,IS)*SX(5)
            IF(.NOT.(RSSX.GT.0 .OR. ISENSEX(IS).EQ.1)) RSSX=0.
            RSSX=RSSX/TAUX(IS)
            GAMDOT(IS,KGX)=GAMD0X(IS)*ABS(RSSX)**NRSX(IS)*SIGN(1.,RSSX)
          ENDDO

          KGX=KGX+1

        ENDDO      ! END OF DO KKK
      ENDDO      ! END OF DO IPH

      RETURN
      END
C
C *****************************************************************
C     SUBROUTINE GET_THEFLU     --->     VERSION 4/JAN/07
C
C     GETS THETA (i.e. THE TERM INVOLVING ESHELBY TENSOR DERIVATIVES)
C     NEEDED FOR FLUCTUATION CALCULATION
C ******************************************************************

      SUBROUTINE SO_GET_THEFLU(THEFLU,IPH)

      USE VPSC8DIM
      USE C4GAFLU
	  
      DIMENSION THEFLU(5,5,5,5,NPHPEL)
      DIMENSION DLDM5(5,5)
      DIMENSION DLDMSA(3,3,3,3),DLDMGA(3,3,3,3)

      DIMENSION EIGB(3,3),AXB(3)
      DIMENSION DS(3,3,3,3),DS5(5,5)
      DIMENSION THE1(5,5,5,5),THE2(5,5)
      DIMENSION DSDDM1(3,3,3,3),DSDDM2(3,3,3,3)
c
cw      DIMENSION SSKEW(3,3)
c!
      DO IR=1,5
      DO IS=1,5
c
      DO I=1,5
      DO J=1,5
      DLDM5(I,J)=-0.5*
     #    (LBARTG(I,IR)*LBARTG(IS,J)+LBARTG(I,IS)*LBARTG(IR,J))
      ENDDO
      ENDDO
c!
c!      CALL CHG_BASIS(AUX5,AUX33,DL5,DLSA,3,5)
      CALL CHG_BASIS(AUX5,AUX33,DLDM5,DLDMSA,3,5)
C
C     ROTATION OF DLSA TO ELLIPSOID (GRAIN) PRINCIPAL AXES
C
cc        if(ishape(iph).le.1) then
          do j=1,3
            axb(j)=axisph(0,j,iph)
            do i=1,3
              eigb(i,j)=axisph(i,j,iph)
            enddo
          enddo
cc        endif
c
      DO 95 I=1,3
      DO 95 J=1,3
      DO 95 M=1,3
      DO 95 N=1,3
        DUMMYB=0.
        DO 90 I1=1,3
        DO 90 J1=1,3
        DO 90 M1=1,3
        DO 90 N1=1,3
c!
          DUMMYB=DUMMYB+EIGB(I1,I)*EIGB(J1,J)*EIGB(M1,M)
     #           *EIGB(N1,N)*DLDMSA(I1,J1,M1,N1)
c!
   90   CONTINUE
c!
        DLDMGA(I,J,M,N)=DUMMYB
c!
   95 CONTINUE

      CALL ESHELBY(AXB,C4GA_FLU,0.,AUX3333,AUX3333,AUX33,AUX33,PDIL,
     #             DLDMGA,DSDDM1,4)
c
      CALL ESHELBY(AXB,C4GA_FLU,0.,AUX3333,AUX3333,AUX33,AUX33,PDIL,
     #             DLDMGA,DSDDM2,5)
cw
      DO 130 I=1,3
      DO 130 J=1,3
      DO 130 M=1,3
      DO 130 N=1,3
        DUMMYE=0.
        DO 120 I1=1,3
        DO 120 J1=1,3
        DO 120 M1=1,3
        DO 120 N1=1,3
          DUMMYE=DUMMYE+EIGB(I,I1)*EIGB(J,J1)*EIGB(M,M1)
     #        *EIGB(N,N1)*(DSDDM1(I1,J1,M1,N1)+DSDDM2(I1,J1,M1,N1))
  120   CONTINUE
        DS(I,J,M,N)=DUMMYE
  130 CONTINUE
c
      CALL CHG_BASIS(AUX5,AUX33,DS5,DS,4,5)
c
c     THEFLU
c
      DO I=1,5
      DO L=1,5
      DUMMY=0.
      DO K=1,5
      DUMMY=DUMMY+ImSINVPH(I,K,IPH)*DS5(K,L)
      ENDDO
      THE1(I,L,IR,IS)=DUMMY
      ENDDO
      ENDDO
c
      DO I=1,5
      DO J=1,5
      DUMMY=0.
      DO K=1,5
      DUMMY=DUMMY+(FSPH(I,K,IPH)+XID5(I,K))*MBARTG(K,J)
      ENDDO
      THE2(I,J)=DUMMY
      ENDDO
      ENDDO
c
      DO I=1,5
      DO J=1,5
      DUMMY=FSPH(I,IR,IPH)*XID5(J,IS)
      DO L=1,5
      DUMMY=DUMMY+THE1(I,L,IR,IS)*THE2(L,J)
      ENDDO
      THEFLU(I,J,IR,IS,iph)=DUMMY
      ENDDO
      ENDDO
c
      ENDDO     ! IS ENDDO
      ENDDO     ! IR ENDDO
c
      RETURN
      END
c
C *****************************************************************
C     SUBROUTINE SO_PROCEDURE
C
C     BUG RELATED TO TWINNING UNIDIRECTIONALITY
C     CORRECTED --> RL, DEC 24, 2008
C
C     USING THE STRESS SECOND-ORDER MOMENT TENSOR OF GRAIN KKK,
C     CALCULATES THE ALPHA'S AND E'S FOR EACH SYSTEM, NEEDED
C     TO DEFINE THE SECOND-ORDER LINEARIZATION, ACCORDING TO
C     THE 'CONSTITUTIVE RELATION' VERSION OF PONTE CASTANEDA'S
C     SECOND-ORDER FORMULATION
C ******************************************************************
c
      subroutine so_procedure(iph1,kgx,kkk,asonew,esonew)
c
      USE VPSC8DIM

      dimension asonew(nsysmx),esonew(nsysmx)
      dimension csigma(5,5),schc1(5)
      dimension sc(5,nsysmx),taux(nsysmx),nrsx(nsysmx)
      dimension isensex(nsysmx)
c
      DO IS=1,NSYST(IPH1)
       ISENSEX(IS)=ISENSE(IS,IPH1)

        if(irsvar.eq.1) then
        NRSX(IS)=JXRS
        else
        NRSX(IS)=NRS(IS,IPH1)
        endif
c
        TAUX(IS)=   CRSS(IS,KKK)
        DO J=1,5
          SC(J,IS)=SCH(J,IS,KGX)
        ENDDO
      ENDDO
c
      do i=1,5
      do j=1,5
      csigma(i,j)=secmom5(i,j,kgx)-sg(i,kgx)*sg(j,kgx)
      enddo
      enddo
c
      DO IS=1,NSYST(IPH1)
c
      do j=1,5
      schc1(j)=0.
      do i=1,5
cc      schc1(j)=schc1(j)+sch(i,is,kgx)*csigma(i,j)
      schc1(j)=schc1(j)+sc(i,is)*csigma(i,j)
      enddo
      enddo
c
      schc=0.
      do i=1,5
cc      schc=schc+schc1(i)*sch(i,is,kgx)
      schc=schc+schc1(i)*sc(i,is)
      enddo
c
      if(schc.lt.0) then
c      write(*,*)
      write(*,*) ' *** WARNING ***  m:C:m < 0 ',kgx,is,schc
      write(*,*)
cc
      schc=abs(schc)
      endif
cc
ccc      asonew(is)=aso(is,kgx)
ccc      esonew(is)=eso(is,kgx)
ccc      else
c
      schc=sqrt(schc)
      taubar=0.
      do i=1,5
cc      taubar=taubar+sch(i,is,kgx)*sg(i,kgx)
      taubar=taubar+sc(i,is)*sg(i,kgx)
      enddo
      tauhat=taubar+sign(1.,taubar)*schc
cc
      phipbar=sign(1.,taubar)*(abs(taubar)/taux(is))**nrsx(is)
      phiphat=sign(1.,tauhat)*(abs(tauhat)/taux(is))**nrsx(is)
c
      phihat=taux(is)/(nrsx(is)+1)*(abs(tauhat)/taux(is))**(nrsx(is)+1)
cv
      if(isensex(is).eq.0.and.taubar.lt.0) then
      asonew(is)=0.
      esonew(is)=0.
      else
      asonew(is)=(phiphat-phipbar)/(tauhat-taubar)
      esonew(is)=phipbar-asonew(is)*taubar
      endif
cv
c      write(*,'(5e12.3)') taubar,tauhat,phipbar,phiphat,asonew(is)
c
ccc      endif
cv
      if(.not.(isensex(is).eq.0.and.taubar.lt.0)) then
      utilde=utilde+wgt(kkk)*(phihat+phipbar*(taubar-tauhat))
      endif
cv
      ENDDO
c
      return
      end

C *****************************************************************
C     SUBROUTINE SO_MOD     --->     VERSION 4/JAN/07
C
C     WITH THE ALPHA'S AND E'S FOR EACH SYSTEM AND EACH GRAIN
C     CALCULATES THE LOCAL LINEARIZED MODULI, ACCORDING TO
C     THE 'CONSTITUTIVE RELATION' VERSION OF PONTE CASTANEDA'S
C     SECOND-ORDER FORMULATION
C ******************************************************************

      SUBROUTINE SO_MOD
c
      USE VPSC8DIM

      DIMENSION AUXTG(5,5)
c
      KGX=1
      DO IPH=IPHBOT,IPHTOP
      IPHEL=IPH-IPHBOT+1
c
      DO KKK=NGR(IPH-1)+1,NGR(IPH)
c
      DO I=1,5
      DO J=1,5
      AUXTG(I,J)=0.
      DO IS=1,NSYST(IPHEL)
      AUXTG(I,J)=AUXTG(I,J)+SCH(I,IS,KGX)*SCH(J,IS,KGX)*ASO(IS,KGX)
      ENDDO
      ENDDO
      ENDDO
c
      DO I=1,5
      DG0(I,KGX)=0.
      DO IS=1,NSYST(IPHEL)
      DG0(I,KGX)=DG0(I,KGX)+SCH(I,IS,KGX)*ESO(IS,KGX)
      ENDDO
      ENDDO
c
      do i=1,5
      do j=1,5
      MGTG(i,j,kgx)=auxtg(i,j)
      enddo
      enddo
c
      KGX=KGX+1
      ENDDO
c
      ENDDO
c
      return
      end

C *****************************************************************
C     SUBROUTINE SO_SDPX     --->     VERSION 4/JAN/07
C
C     CALCULATES THE STANDARD DEVIATIONS OF THE EQUIVALENT STRESS
C     AND THE STRAIN-RATE IN THE PX, PER GRAIN, AND PER ROLLING FCC
C     TEXTURE COMPONENT, IF APPROPRIATE
C ******************************************************************

      SUBROUTINE SO_SDPX

      USE VPSC8DIM
      USE CUB_COMP

      dimension seq1cub(0:6),seq2cub(0:6),deq1cub(0:6),deq2cub(0:6)

      SEQ1PX=0.
      SEQ2PX=0.

      DEQ1PX=0.
      DEQ2PX=0.

      IWRITE=0                         ! CNT added 2022 07 21
      IF(NWRITE.EQ.0 .AND. NWRITEX.EQ.NSTEPS) IWRITE=1 
      IF(NWRITE.NE.0 .AND. MOD(NWRITEX,NWRITE).EQ.0) IWRITE=1 

      if(icubcom.eq.1) then
      do i=0,6
      seq1cub(i)=0.
      seq2cub(i)=0.
      deq1cub(i)=0.
      deq2cub(i)=0.
      enddo
      endif

      IF(IWRITE.EQ.1) write(83,'(a)')
     #'   GR#       SEQ1      SEQ2      SDSEQ        DEQ1      DEQ2
     # SDDEQ'

      KGX=1
      DO IPH=IPHBOT,IPHTOP
      DO KKK=NGR(IPH-1)+1,NGR(IPH)
c
c     double bar equivalent
c
      DUMMY=0.
      DUMMYD=0.
      DO K=1,5
      DUMMY=DUMMY+SECMOM5(K,K,KGX)
      DUMMYD=DUMMYD+SECMOM5D(K,K,KGX)
      ENDDO
      SEQ2(KGX)=SQRT(3./2.*DUMMY)
      DEQ2(KGX)=SQRT(2./3.*DUMMYD)
c
      SEQ2PX=SEQ2PX+3./2.*DUMMY*WGT(KKK)
      DEQ2PX=DEQ2PX+2./3.*DUMMYD*WGT(KKK)
c
      if(icubcom.eq.1) then
      SEQ2CUB(IGRTYPE(KKK))=SEQ2CUB(IGRTYPE(KKK))+
     #        3./2.*DUMMY*WGT(KKK)/WIDMOD(IGRTYPE(KKK))
      DEQ2CUB(IGRTYPE(KKK))=DEQ2CUB(IGRTYPE(KKK))+
     #        2./3.*DUMMYD*WGT(KKK)/WIDMOD(IGRTYPE(KKK))
      endif
c
c     single bar equivalent
c
      DUMMY=0.
      DUMMYD=0.
      DO K=1,5
      DUMMY=DUMMY+SG(K,KKK)**2
      DUMMYD=DUMMYD+DG(K,KGX)**2
      ENDDO
c
      SEQ1=SQRT(3./2.*DUMMY)
      DEQ1=SQRT(2./3.*DUMMYD)
c
      SEQ1PX=SEQ1PX+3./2.*DUMMY*WGT(KKK)
      DEQ1PX=DEQ1PX+2./3.*DUMMYD*WGT(KKK)
c
      if(icubcom.eq.1) then
      SEQ1CUB(IGRTYPE(KKK))=SEQ1CUB(IGRTYPE(KKK))+
     #        3./2.*DUMMY*WGT(KKK)/WIDMOD(IGRTYPE(KKK))
      DEQ1CUB(IGRTYPE(KKK))=DEQ1CUB(IGRTYPE(KKK))+
     #        2./3.*DUMMYD*WGT(KKK)/WIDMOD(IGRTYPE(KKK))
      endif
c
      if(seq2(kgx).gt.seq1) then
      SDSEQ=SQRT(SEQ2(KGX)**2-SEQ1**2)/SVM
      else
      SDSEQ=0.
      endif
c
      if(deq2(kgx).gt.deq1) then
      SDDEQ=SQRT(DEQ2(KGX)**2-DEQ1**2)/DVM
      else
      SDDEQ=0.
      endif

c      WRITE(*,'(a,i5,a,3f10.3)')
c     # ' GRAIN = ',KGX,' - SEQ1,SEQ2,SD = ',SEQ1,SEQ2(KGX),SDSEQ

c      WRITE(*,'(a,i5,a,3f10.3)')
c    # ' GRAIN = ',KGX,' - DEQ1,DEQ2,SD = ',DEQ1,DEQ2(KGX),SDDEQ

c      WRITE(*,'(1H+,i5,3x,3f10.3,3x,3f10.3)')
c     #      KGX,SEQ1,SEQ2(KGX),SDSEQ,DEQ1,DEQ2(KGX),SDDEQ

c      WRITE(83,'(a,i5,a,3f10.3)')
c     # ' GRAIN = ',KGX,' - SEQ1,SEQ2,SD = ',SEQ1,SEQ2(KGX),SDSEQ

      IF(IWRITE.EQ.1) WRITE(83,'(i5,3x,2e10.3,f10.3,3x,2e10.3,f10.3)')
     #      KGX,SEQ1,SEQ2(KGX),SDSEQ,DEQ1,DEQ2(KGX),SDDEQ

      KGX=KGX+1
      ENDDO
      ENDDO

C --> CNT changed because occasional NaN error message (03/03/2017)
      SDSEQINTRA=MAX(SEQ2PX-SVM**2, 0.)
      SDSEQINTRA=SQRT(SDSEQINTRA)/SVM
      SDSEQINTER=MAX(SEQ1PX-SVM**2, 0.)
      SDSEQINTER=SQRT(SDSEQINTER)/SVM

      SDDEQINTRA=MAX(DEQ2PX-DVM**2, 0.)
      SDDEQINTRA=SQRT(SDDEQINTRA)/DVM
      SDDEQINTER=MAX(DEQ1PX-DVM**2, 0.)
      SDDEQINTER=SQRT(SDDEQINTER)/DVM
c
      write(*,*)
      write(*,*) 'S SDPX (inter+intra) =',SDSEQINTRA
      write(*,*) 'S SDPX (inter)       =',SDSEQINTER
      write(*,*) 'D SDPX (inter+intra) =',SDDEQINTRA
      write(*,*) 'D SDPX (inter)       =',SDDEQINTER
c
      write(83,'(a,f10.3)') 'S SDPX (inter+intra) = ',SDSEQINTRA
      write(83,'(a,f10.3)') 'S SDPX (inter)       = ',SDSEQINTER
      write(83,'(a,f10.3)') 'D SDPX (inter+intra) = ',SDDEQINTRA
      write(83,'(a,f10.3)') 'D SDPX (inter)       = ',SDDEQINTER
c
      if(icubcom.eq.1) then
        write(84,'(a,f10.3)') 
     # 'S SDPX (inter+intra) = ',SQRT(SEQ2PX-SVM**2)
        write(84,'(a,f10.3)')
     * 'D SDPX (inter+intra) = ',SQRT(DEQ2PX-DVM**2)
        write(84,'(a)')
     #'  CMP        SEQ1      SEQ2      SDSEQ        DEQ1      DEQ2
     # SDDEQ'

      do i=0,6
        SDSEQCUB=SQRT(SEQ2CUB(I)-SEQ1CUB(I))
        SDDEQCUB=SQRT(DEQ2CUB(I)-DEQ1CUB(I))
        write(84,'(2x,a3,3x,3f10.3,3x,3f10.3)')
     #      idlabel(i),sqrt(seq1cub(i)),sqrt(seq2cub(i)),sdseqcub,
     #        sqrt(deq1cub(i)),sqrt(deq2cub(i)),sddeqcub
      enddo
      endif
c
      return
      end

C *****************************************************************
C     SUBROUTINE SO_GRAIN_STRESS_ALT     --->     VERSION 4/JAN/07
C
C     IN THE SO CASE, CALCULATES THE GRAIN STRESSES
c     AND STRAIN-RATES USING THE LOCALIZATION EQUATION
C ******************************************************************

      SUBROUTINE SO_GRAIN_STRESS_ALT (KGX,KKK,IPHEL,IPH)

      USE VPSC8DIM
c
      do i=1,5
      sg(i,kgx)=SG0(i,kgx)
      do j=1,5
      sg(i,kgx)=sg(i,kgx)+BG(i,j,kgx)*sastbar(j)
      enddo
      enddo
c
      do i=1,5
      dg(i,kgx)=DG0(i,kgx)
      do j=1,5
      dg(i,kgx)=dg(i,kgx)+MGTG(i,j,kgx)*sg(j,kgx)
      enddo
      enddo

      RETURN
      END

C *****************************************************************
C     SUBROUTINE SO_EXTRAPOL     --->     VERSION 4/JAN/07
C
C     CALCULATES THE ALPHA'S AND E'S FOR A GIVEN RATE-SENSITIVITY
C     EXPONENT USING THE VALUES OF THE THREE PREVIOUS RS EXPONENTS
c     (QUADRATIC EXTRAPOLATION)
C *****************************************************************

      SUBROUTINE SO_EXTRAPOL (X,IRS,IOPTION)

      USE VPSC8DIM
      USE EXTRAPSO

      IF(IOPTION.EQ.1) THEN

        IF(IRS.LE.3) THEN
          XMRSOLD(IRS)=X
        ELSE
          XMRSOLD(1)=XMRSOLD(2)
          XMRSOLD(2)=XMRSOLD(3)
          XMRSOLD(3)=X
        ENDIF

        KGX=1
        DO IPH=IPHBOT,IPHTOP
        IPHEL=IPH-IPHBOT+1
        DO KKK=NGR(IPH-1)+1,NGR(IPH)
          DO IS=1,NSYST(IPHEL)
            IF(IRS.LE.3) THEN
              ESOOLD(IS,KGX,IRS)=ESO(IS,KGX)
              ASOOLD(IS,KGX,IRS)=ASO(IS,KGX)
            ELSE
              ESOOLD(IS,KGX,1)=ESOOLD(IS,KGX,2)
              ASOOLD(IS,KGX,1)=ASOOLD(IS,KGX,2)
              ESOOLD(IS,KGX,2)=ESOOLD(IS,KGX,3)
              ASOOLD(IS,KGX,2)=ASOOLD(IS,KGX,3)
              ESOOLD(IS,KGX,3)=ESO(IS,KGX)
              ASOOLD(IS,KGX,3)=ASO(IS,KGX)
            ENDIF
          ENDDO
          KGX=KGX+1
        ENDDO
        ENDDO

      ENDIF      ! END OF IF(IOPTION.EQ.1)

      IF(IOPTION.EQ.2) THEN

        KGX=1
        DO IPH=IPHBOT,IPHTOP
        IPHEL=IPH-IPHBOT+1
        DO KKK=NGR(IPH-1)+1,NGR(IPH)
          DO IS=1,NSYST(IPHEL)

          ASOTRY=
     # (X-XMRSOLD(2))*(X-XMRSOLD(3))/
     # (XMRSOLD(1)-XMRSOLD(2))/(XMRSOLD(1)-XMRSOLD(3))*ASOOLD(IS,KGX,1)+
     # (X-XMRSOLD(1))*(X-XMRSOLD(3))/
     # (XMRSOLD(2)-XMRSOLD(1))/(XMRSOLD(2)-XMRSOLD(3))*ASOOLD(IS,KGX,2)+
     # (X-XMRSOLD(1))*(X-XMRSOLD(2))/
     # (XMRSOLD(3)-XMRSOLD(1))/(XMRSOLD(3)-XMRSOLD(2))*ASOOLD(IS,KGX,3)

          IF(ASOTRY.GT.0) THEN
            ASO(IS,KGX)=ASOTRY
            ESO(IS,KGX)=
     # (X-XMRSOLD(2))*(X-XMRSOLD(3))/
     # (XMRSOLD(1)-XMRSOLD(2))/(XMRSOLD(1)-XMRSOLD(3))*ESOOLD(IS,KGX,1)+
     # (X-XMRSOLD(1))*(X-XMRSOLD(3))/
     # (XMRSOLD(2)-XMRSOLD(1))/(XMRSOLD(2)-XMRSOLD(3))*ESOOLD(IS,KGX,2)+
     # (X-XMRSOLD(1))*(X-XMRSOLD(2))/
     # (XMRSOLD(3)-XMRSOLD(1))/(XMRSOLD(3)-XMRSOLD(2))*ESOOLD(IS,KGX,3)
          ELSE
            ASO(IS,KGX)=ASOOLD(IS,KGX,3)
            ESO(IS,KGX)=ESOOLD(IS,KGX,3)
          ENDIF

          IF(ASO(IS,KGX).LT.0.) THEN
            WRITE(*,*) '  ASO(IS,KGX) NEGATIVE IN EXTRAPOL'
            WRITE(*,*) KGX,IS,ASO(IS,KGX)
            print *, 'enter c to continue'
            read  *
          ENDIF

          ENDDO
          KGX=KGX+1
        ENDDO
        ENDDO

      ENDIF      ! END OF IF(IOPTION.EQ.2)

      RETURN
      END

C ************************************************************
C     SUBROUTINE SO_VOIGT10     --->     VERSION 4/JAN/07
C ************************************************************

      SUBROUTINE SO_VOIGT10(T1,T2,IOPT)

      DIMENSION T1(10),T2(4,4)
      DIMENSION IJV(10,2)

      ijv(1,1)=1
      ijv(1,2)=1
      ijv(2,1)=2
      ijv(2,2)=2
      ijv(3,1)=3
      ijv(3,2)=3
      ijv(4,1)=2
      ijv(4,2)=3
      ijv(5,1)=1
      ijv(5,2)=3
      ijv(6,1)=1
      ijv(6,2)=2
      ijv(7,1)=1
      ijv(7,2)=4
      ijv(8,1)=2
      ijv(8,2)=4
      ijv(9,1)=3
      ijv(9,2)=4
      ijv(10,1)=4
      ijv(10,2)=4

      IF(IOPT.EQ.1) THEN
      DO 30 I=1,10
      I1=IJV(I,1)
      I2=IJV(I,2)
      T2(I1,I2)=T1(I)
   30 T2(I2,I1)=T1(I)
      ENDIF
C
      IF(IOPT.EQ.2) THEN
      DO 40 I=1,10
      I1=IJV(I,1)
      I2=IJV(I,2)
   40 T1(I)=T2(I1,I2)
      ENDIF
C
      RETURN
      END

C *****************************************************************************
C   END SUBROUTINES FOR DOING FLUCTUATIONS AND SECOND ORDER (SO)
C *****************************************************************************

      SUBROUTINE STATE_NxN (ibc_d,BC_D,ibc_s,BC_S,XMS,ndim)

      dimension aux11(ndim),aux21(ndim,ndim),XMS(ndim,ndim)
      dimension ibc_d(ndim),bc_d(ndim),ibc_s(ndim),bc_s(ndim)

      do i=1,ndim
        aux11(i)=-1.d0*ibc_d(i)*bc_d(i)
        do j=1,ndim
          aux11(i)=aux11(i)+xms(i,j)*ibc_s(j)*bc_s(j)
          aux21(i,j)=ibc_s(j)*(i/j)*(j/i)-ibc_d(j)*xms(i,j)
        enddo
      enddo

      CALL LU_EQSYSTEM(AUX21,AUX11,ndim,IER)

      if(ier.eq.1) then
      write(*,*) 'SINGULAR SYSTEM IN STATE_5X5'
      stop
      endif

      do i=1,ndim
        bc_d(i)=ibc_d(i)*bc_d(i)+ibc_s(i)*aux11(i)
        bc_s(i)=ibc_s(i)*bc_s(i)+ibc_d(i)*aux11(i)
      enddo

      return
      end
C
C *****************************************************************************
C     SUBROUTINE STATE_6x6      --->      VERSION 31/AUG/2017

C     SOLVES A 6x6 SYSTEMS FOR CARTESIAN STRESS-STRAIN COMPONENTS EXPRESSED
C     IN VOIGT NOTATION 

C *** BASED ON SECANT EQUATION: D(i)-D0(i)=MS(i,j)*S(j), WHERE MS IS COMPLIANCE,
C     SHIFTS KNOWN COMPONENTS TO INDEPENDENT TERM, CALCULATES COEFFICIENTS
C     FOR UNKNOWN COMPONENTS AND INVERTS THE LINEAR SYSTEM.
C *** COMPONENTS OD 'D', 'S' and 'MS' ENTER IN CARTESIAN NOTATION BUT GET 
C     TRANSFORMED INTO 6x1 OR 6x6 VOIGT MATRICES. 

C *** CALLED AT THE END OF EACH INCREMENTAL STEP IN SUBROUTINE VPSC IF
C     IDBARv1*IDBARv2*IDBARv3=0 (NON FULLY ENFORCED DIAGONAL)
C *****************************************************************************
c
      subroutine state_6x6 (ibc_d,AUX_D,ibc_s,AUX_S,AUX_MS)
c
      dimension aux6v(6),aux66v(6,6)
      dimension AUX_D(3,3),AUX_S(3,3),AUX_MS(3,3,3,3)
      dimension ibc_d(6),bc_d(6),ibc_s(6),bc_s(6)
      dimension profac(6)
      dimension aux33(3,3),aux3333(3,3,3,3),aux66(6,6)

      do i=1,6
        profac(i)=1.d0+(i/4)
      enddo

C *** CONVERTS CARTESIAN 'D','S' and 'MS' TO VOIGT.
        CALL VOIGT (BC_D,AUX_D,AUX66,AUX3333,2)
        CALL VOIGT (BC_S,AUX_S,AUX66,AUX3333,2)
        CALL VOIGT (AUX6v,AUX33,AUX66,AUX_MS,4)

      do i=1,6
        aux6v(i)=-1.d0*ibc_d(i)*bc_d(i)
        do j=1,6
          aux6v(i)=aux6v(i)+aux66(i,j)*ibc_s(j)*bc_s(j)*profac(j)
          aux66v(i,j)=ibc_s(j)*(i/j)*(j/i)-ibc_d(j)*aux66(i,j)*profac(j)
        enddo
      enddo

      CALL LU_EQSYSTEM(AUX66v,AUX6v,6,IER)

      if(ier.eq.1) then
        write(*,*) 'SINGULAR SYSTEM IN STATE_6X6'
        stop
      endif

      do i=1,6
        bc_d(i)=ibc_d(i)*bc_d(i)+ibc_s(i)*aux6v(i)
        bc_s(i)=ibc_s(i)*bc_s(i)+ibc_d(i)*aux6v(i)
      enddo
	  
C *** CONVERTS VOIGT 'D','S' TO CARTESIAN
      CALL VOIGT (BC_D,AUX_D,AUX66,AUX3333,1)
      CALL VOIGT (BC_S,AUX_S,AUX66,AUX3333,1)

      return
      end

C ***************************************************************************
C     SUBROUTINE STAT_GRAIN_SHAPE    --->    VERSION 13/SEP/2022
C     PERFORMS STATISTICS ON GRAIN AXES (AVERAGE & ST. DEV.) IF ISHAPE=1,2,3
C
C     CNT MODIFIED SET/2005:
C     --> OUTPUTTING IS MADE LOCALLY WITHIN THE SUBROUTINE.
C     --> OUTPUT FILES 'AVAXES' & 'GRAXFAC' WERE MERGED INTO 'STAT_AXES'
C ***************************************************************************

      SUBROUTINE STAT_GRAIN_SHAPE (ISTEP)

      USE VPSC8DIM


      DO IPH=IPHBOT,IPHTOP

      IUNIT=45+IPH
      IF(ISTEP.EQ.1) WRITE(IUNIT,'(''  EPS  PHASE'',3X,
     #  ''   AVAX1   AVAX2   AVAX3      SDAX1   SDAX2   SDAX3'')')  

      IF(ISHAPE(IPH).NE.0 .AND. WPH(IPH).GT.0.) THEN

        AVAX(:,IPH)=0.
        SDAX(:,IPH)=0.

        DO KKK=NGR(IPH-1)+1,NGR(IPH)
          do i=1,3
            avax(i,iph)=avax(i,iph)+AXISGR(0,I,KKK)*wgt(kkk)/wph(iph)
          enddo
        ENDDO

        DO KKK=NGR(IPH-1)+1,NGR(IPH)
          do i=1,3
            sdax(i,iph)=sdax(i,iph)+(AXISGR(0,I,KKK)-avax(i,iph))**2*
     #                wgt(kkk)/wph(iph)
          enddo
        ENDDO
        SDAX(:,IPH)=SQRT(SDAX(:,IPH))/avax(:,iph)

        WRITE(IUNIT,'(F6.2,I6,3X,3F8.3,3X,3F8.3)')
     #           EPSACU,IPH,(AVAX(i,iph),i=1,3),(SDAX(i,iph),i=1,3)

      ENDIF

      ENDDO

      RETURN
      END
C
C *****************************************************************************
C     SUBROUTINE STAT_SHEAR_ACTIVITY      --->      VERSION 2020/07/30
C
C     MAKES STATISTIC ON DEFORMATION MODE ACTIVITY.
C     CALCULATES ARRAYS:
C     - GMAX(kmo,iph) : maximum dgam for a given mode & phase in a given grain.
C     - GTOT(kmo,kgr) : sum of dgam for a given mode in a given grain.
C     - GAVGR(kgr)    : sum of all dgam for a given grain.
C     - GMAXG: maximum dgam for a given grain.
C     - ACTGR(kgr)     : number of active systems in a given grain.
C     - ACTPH(iph)     : weighted number of active systems in a given phase.
C     - GAVMOD(kmo,iph): weighted sum of dgam for a given mode in a given phase.
C     - GAVPH(iph)     : weighted sum of dgam for a given phase.
C *****************************************************************************

      SUBROUTINE STAT_SHEAR_ACTIVITY

      USE VPSC8DIM

      DIMENSION GMAX(NMODMX,NPHMX),GTOT(NMODMX,NGRMX)

      KGX=1
      DO IPH=IPHBOT,IPHTOP
        IPHEL=IPH-IPHBOT+1
      DO KGR=NGR(IPH-1)+1,NGR(IPH)

        gavgr(kgr)=0.
        ACTGR(kgr)=0.
        gmaxg=0.
        is=0
        do kmo=1,nmodes(IPHEL)
          gmax(kmo,iph)=0.
          gtot(kmo,kgr)=0.
            do ksy=1,nsm(kmo,IPHEL)
              is=is+1
              gabs=abs(gamdot(is,KGX))
              gtot(kmo,kgr)=gtot(kmo,kgr)+gabs
              if(gabs.gt.gmax(kmo,iph)) gmax(kmo,iph)=gabs
              if(gabs.gt.gmaxg) gmaxg=gabs
            enddo
          gavgr(kgr)=gavgr(kgr)+gtot(kmo,kgr)
        enddo
        gming=gmaxg*0.05      ! arbitrary threshold for labeling activity
        do isy=1,nsyst(IPHEL)
          if(abs(gamdot(isy,KGX)).ge.gming) ACTGR(kgr)=ACTGR(kgr)+1
        enddo

        KGX=KGX+1
      ENDDO      ! end of DO KGR
      ENDDO      ! end of DO IPH

      GAVPC=0.
      DO IPH=IPHBOT,IPHTOP
        IPHEL=IPH-IPHBOT+1

        ACTPH(iph)=0.
        GAVPH(IPH)=0.
        do kmo=1,nmodes(IPHEL)
          gavmod(kmo,iph)=0.
        enddo

        IF(WPH(IPH).NE.0.) THEN
          do kgr=ngr(iph-1)+1,ngr(iph)
            do kmo=1,nmodes(IPHEL)
              gavmod(kmo,iph)=gavmod(kmo,iph)+gtot(kmo,kgr)*wgt(kgr)
            enddo
            GAVPH(IPH)=GAVPH(IPH)+gavgr(kgr)*wgt(kgr)
            ACTPH(iph)=ACTPH(iph)+ACTGR(kgr)*wgt(kgr)*gavgr(kgr)
          enddo
          ACTPH(iph)=ACTPH(iph)/GAVPH(IPH)
          GAVPC=GAVPC+GAVPH(IPH)
        ENDIF
      ENDDO

      DO IPH=IPHBOT,IPHTOP
        IPHEL=IPH-IPHBOT+1
        do kmo=1,nmodes(IPHEL)
          gavmod(kmo,iph)=gavmod(kmo,iph)/GAVPC
        enddo

C     WRITE(*,'(''INSIDE STAT_SHEAR -> IPH WPH GAVMOD'',I5,
C    #   6F10.4)') IPH,WPH(IPH),
C    #   (GAVMOD(KMO,IPH),KMO=1,NMODES(IPHEL))

      ENDDO

      RETURN
      END
C
C *********************************************************************************
C     SUBROUTINE STAT_STRESS_STRAIN   --->     VERSION 24/JUL/2022
C
C     AVERAGE & STANDARD DEVIATION OF b_BASIS COMPONENTS OF STRESS, STRAIN-RATE
C     & ELASTIC STRAIN.
C     NORMALIZATION IS DONE USING THE NORM OF THE TENSOR INSTEAD OF INDIVIDUAL
C     AVERAHE COMPONENTS, WHICH MAY AVERAGE TO ZERO.
C     'D'   CALCULATIONS APPLY TO VISCOUS CASES (INTERACTION=0,1,2...). 
C     'EEL' CALCULATIONS APPLY TO THERMO-ELASTIC CASES (INTERACTION=-1).
C *********************************************************************************

      SUBROUTINE STAT_STRESS_STRAIN

      USE VPSC8DIM
	  
      sav(:)   =0.
      dav(:)   =0.
      eav(:)   =0.
      stdevs(:)=0.
      stdevd(:)=0.
      stdeve(:)=0.

      KGX=1
      DO IPH=IPHBOT,IPHTOP
      DO KGR=ngr(iph-1)+1,ngr(iph)
        do i=1,6
          sav(i)   =sav(i)+sg(i,KGR)   *wgt(KGR)
          dav(i)   =dav(i)+dg(i,KGX)   *wgt(KGR)
          eav(i)   =eav(i)+eelgr(i,KGX)*wgt(KGR)
        enddo
        KGX=KGX+1
      ENDDO
      ENDDO

      KGX=1
      DO IPH=IPHBOT,IPHTOP
      DO KGR=ngr(iph-1)+1,ngr(iph)
        do i=1,6
          difs=sg(i,KGR)-sav(i)
          difd=dg(i,KGX)-dav(i)
          dife=eelgr(i,KGX)-eav(i)
          stdevs(i)=stdevs(i)+wgt(KGR)*difs**2
          stdevd(i)=stdevd(i)+wgt(KGX)*difd**2
          stdeve(i)=stdeve(i)+wgt(KGX)*dife**2
        enddo
        KGX=KGX+1
      ENDDO
      ENDDO

C *** PROVISION MADE IN CASE SAV, or DAV, or EAV=0 (i.e.: cooling)
      STRESS_NORM =VNORM(SAV,6)
      STRAINR_NORM=VNORM(DAV,6)
      ELSTR_NORM  =VNORM(EAV,6)
c        write(11,'('' SAV'',6e12.3,''   NORM'',E12.3)') sav,stress_norm
c        write(11,'('' DAV'',6e12.3,''   NORM'',E12.3)') dav,strainr_norm
c        write(11,'('' EAV'',6e12.3,''   NORM'',E12.3)') eav,elstr_norm
      IF (STRESS_NORM .LE.1.E-10) STRESS_NORM =1.
      IF (STRAINR_NORM.LE.1.E-10) STRAINR_NORM=1.
      IF (ELSTR_NORM  .LE.1.E-10) ELSTR_NORM  =1.

      stdevs(:)=sqrt(stdevs(:))/STRESS_NORM
      stdevd(:)=sqrt(stdevd(:))/STRAINR_NORM
      stdeve(:)=sqrt(stdeve(:))/ELSTR_NORM 

      RETURN
      END
C
C **************************************************************************
C    SUBROUTINE STAT_TWINNING   -->   VERSION 24/APR/2021 
C
C --> TWINNING STATISTICS: SCHMID FACTOR & ACTIVE VARIANT DISTRIBUTION
C --> CALLED INSIDE UPDATE_TWINNING_PTR & UPDATE_TWINNING_VFT
C --> 'ITWS' : TWIN SYSTEM (ABSOLUTE INDEX) PICKED BY ONE OF THE UPDDATE_TWIN
C           ROUTINES TO EITHER SWITCH ORIENTATION OR TO CREATE A NEW GRAIN. 
C     'ITWMODE' : INDEX OF TWIN MODE TO WHICH 'ITWS' BELONGS
C **************************************************************************
C    SCH_STAT(ISLOT,ITWMODE,IPH): FREQUENCY DISTRIBUTION OF TWIN MODE
C    	             'ITWMODE' OF PHASE 'IPH' AT SCHMID FACTOR BIN 'ISLOT'
C    VAR_STAT(IVAR,ITWMODE,IPH): FREQUENCY OF VARIANT 'IVAR' IN TWIN MODE 
C                    'ITWMODE' OF PHASE 'IPH'  
C    VAR_VF(IVAR,ITWMODE,IPH): VOLUME FRACTION OF VARIANT 'IVAR' IN TWIN MODE 
C                    'ITWMODE' OF PHASE 'IPH'  
C --> THESE ARRAYS ARE WRITTEN VIA A CALL TO SUBROUTINE WRITE_TWIN_STATS.
C **************************************************************************
C     WARNING : DIMENSIONING ARRAYS TO '6' AND USING 'DO I=1,6' LOOPS  
C               IMPLICITLY ASSUMES HCP TWINS
C **************************************************************************

      SUBROUTINE STAT_TWINNING(ITWS, ITWMODE, KKK, IPH)

      USE VPSC8DIM

      DIMENSION TWB(3),TWN(3),TWBX(3),TWNX(3)
      DIMENSION TRSS(6),SGXc(3,3)
      INTEGER :: ITWMODE     

      IPHEL=IPH-IPHBOT+1      

C ********************************************************************     
C *** CALCULATES GENERALIZED SCHMID FACTOR FOR TWIN SYSTEM 'ITWS'
C *** PLACES SCHMID FACTOR OF 'ITWS' INTO 0.05 WIDE BINS FOR STATISTICS

      TWN(:) =DNCA(:,ITWS,IPHEL)
      TWB(:) =DBCA(:,ITWS,IPHEL)
      TWNX(:)=0.
      TWBX(:)=0.
      DO I=1,3
        DO J=1,3
          TWNX(I)=TWNX(I)+AG(I,J,KKK)*TWN(J)
          TWBX(I)=TWBX(I)+AG(I,J,KKK)*TWB(J)
        ENDDO
      ENDDO

C --> SDEVIAT suppressed 20190902. Test eventually whether this works.
C --> TEMPORARY FIX WHEN BC'S DO ALLOW TO CALCULATE CAUCHY STRESS SBARc.
C --> NEED TO ADD A SQRT(3/2) WHNE USING THE DEVIATORIC STRESS IN ORDER 
C     KEEP THE SCHMID FACTOR FOR TENSION IN [-0.5,0.5] INTERVAL.
C      IF(IDBARv(1)*IDBARv(2)*IDBARv(3).EQ.1) SBARc(:,:)=SDEVIAT(:,:)   

C --> OPTION FOR USING LOCAL STRESS FOR DEFINING SCHMID FACTOR.
C --> FOR COMPARISON TO EXPERIMENTS NEED TO USE THE MACROSCOPIC STRESS.
C       CALL CHG_BASIS(SG(:,KKK),SGXc,AUX66,AUX3333,1,5)
C       SGXnorm=VNORM(SGXc,9)

C --> OPTION FOR USING MACROSCOPIC STRESS FOR DEFINING SCHMID FACTOR.
C --> FOR COMPARISON TO EXPERIMENTS NEED TO USE THE MACROSCOPIC STRESS.
      SGXc(:,:)=SBARc(:,:)      ! uses macro stress
      SGXnorm=VNORM(SBARc,9)

      SCHX=0.
      DO I=1,3
        DO J=1,3
          SCHX=SCHX+TWNX(I)*TWBX(J)*SGXc(I,J) 
        ENDDO
      ENDDO
      SCHX=SCHX/SGXnorm

      ISLOT=INT(SCHX/0.05) + 1         ! islot=-5,-4, ...14,15
      IF(ISLOT.GT.15 .OR. ISLOT.LT.-5) THEN
        WRITE(*,'('' SGXc'',9F8.3)') SGXc
        WRITE(*,'('' GRAIN'',I5,'' SCHMID FACTOR'',F8.3,'' ISLOT='',I3,
     #  '' EXCEEDS DIMENSION IN SCH_STAT ARRAY'')') KKK,SCHX,ISLOT
        IF(ISLOT.GT.15) ISLOT=15      ! QUICK & DIRTY FIX FOR THE TIME BEING
        IF(ISLOT.LT.-5) ISLOT=-5      ! QUICK & DIRTY FIX FOR THE TIME BEING
      ENDIF
      SCH_STAT(ISLOT,ITWMODE,IPH)=SCH_STAT(ISLOT,ITWMODE,IPH) + 1

C ******************************************************************************
C *** CALCULATES RESOLVED SHEAR 'TRSS' FOR EVERY TWIN IN THIS MODE

      KTMODE=ITWMODE+NSLMOD(IPHEL)
	  ITWREL=0
      ITS=NSLSYS(IPHEL) + (ITWMODE-1)*6   ! beware: 6 twins per mode hardwired    
        DO NSMX=1,NSM(KTMODE,IPHEL)

	      ITWREL=ITWREL+1            ! relative counter for the twin mode
          ITS=ITS+1                  ! absolute counter over all systems
          TWN(:) =DNCA(:,ITS,IPHEL)
          TWB(:) =DBCA(:,ITS,IPHEL)
          TWNX(:)=0.
       	  TWBX(:)=0.
       	  DO I=1,3
       	  DO J=1,3
            TWNX(I)=TWNX(I)+AG(I,J,KKK)*TWN(J)
            TWBX(I)=TWBX(I)+AG(I,J,KKK)*TWB(J)
       	  ENDDO
          ENDDO

          TRSS(ITWREL)=0    
          DO I=1,3
          DO J=1,3
            TRSS(ITWREL)=TRSS(ITWREL)+TWNX(I)*TWBX(J)*SGXc(I,J)
          ENDDO
          ENDDO

        ENDDO      ! end of DO NSMX
C ******************************************************************************
C *** FIND VARIANT CORRESPONDING TO 'ITWS' AND UPDATE STATISTICS

        ITWREL=ITWS-NSLSYS(IPHEL)-(ITWMODE-1)*6   ! 6 twins per mode hardwired
        IVAR=1
        DO I=1,6
          IF (TRSS(ITWREL) .LT. TRSS(I)) IVAR=IVAR+1
        ENDDO

	    VAR_STAT(IVAR,ITWMODE,IPH)= VAR_STAT(IVAR,ITWMODE,IPH) + 1
        IF(ITWINLAW.LE.2) 
     #    VAR_VF(IVAR,ITWMODE,IPH)=VAR_VF(IVAR,ITWMODE,IPH) + WGT(KKK)
        IF(ITWINLAW.EQ.3) THEN
          ITWVFT(1,KKK)=IVAR
          ITWVFT(2,KKK)=ITWMODE
        ENDIF

      RETURN
      END	  	  
C
C ***************************************************************************
C     SUBROUTINE TEXTURE_ROTATION     -->      VERSION OF 27/JUL/2022
C
C     THE PURPOSE OF THIS SUBROUTINE IS TO ROTATE THE CRYSTALLOGRAPIC (crystal
C     axes) AND MORPHOLOGIC (ellipsoid axes) TEXTURES 'ATTACHED' TO A SAMPLE 
C     IN ORDER TO REFER THEM TO A 'TEST SYSTEM'.
C
C     THE APPLICATIONS INTENDED ARE:
C     1- THE TEXTURE IS REFERRED TO A SAMPLE SYSTEM 'SS' (S1,S2,S3) BUT THE VELOCITY
C        GRADIENT REFERS TO A TEST SYSTEM 'TS' (T1,T2,T3) WHERE THE AXES LABELS ARE 
C        DIFFERENT.e.g.: S3 IS THE AXIAL DIR OF A BAR TO BE COMPRESSED, BUT THE  
C        VELOCITY GRADIENT COMPONENTS ARE SUCH THAT T1 IS COMPRESSION.
C     2- THE LANKFORD COEFFICIENT REQUIRES TO REFER TEXTURE TO A TENSILE AXIS THAT
C        IS ON THE ROLLING PLANE AT AN ANGLE 'ALPHA' TO THE RD, WHILE THE ORIGINAL 
C        TEXTURE IS REFERRED TO SA (rd=1,td=2,nd=3)
C     3- IN ECAE THE SAMPLE IS ROTATED BEFORE BEING REINTRODUCED IN THE DIE. THE VEL
C        GRAD COMPONENTS ARE REFERRED TO DIE AXES (TS).
C
C     USER PROVIDES A ROTATION MATRIX 'ROTMAT' DEFINED AS FOLLOWS:     
C     COLUMNS OF 'ROTMAT' ARE THE ROTATED SAMPLE AXES EXPRESSED IN THE TEST AXES.
C     'ROTMAT' TRANSFORMS VECTORS AND TENSORS EXPRESED IN 'SA' INTO 'TA':
C          V_ta(i)  =ROTMAT(i,j)*V_sa(j)
C          T_ta(i,j)=ROTMAT(i,k)*T_sa(k,l)*ROTMAT(j,l)
C
C ***************************************************************************

      SUBROUTINE TEXTURE_ROTATION (ROTMAT)

      USE VPSC8DIM

      DIMENSION ROTMAT(3,3)

      DO IPH=1,NPH
        DO KGX=NGR(IPH-1)+1,NGR(IPH)

C *** 'AG' TRANSFORMS FROM CRYSTAL TO SAMPLE AXES. THE COLUMNS OF 'AG' ARE THE 
C      CRYSTAL AXES EXPRESSED IN 'SA', WHICH USUALLY COINCIDE WITH 'TA'.
C      WHEN 'SA' ARE ROTATED WITH RESPECT TO 'TA' 'ROTMAT' UPDATES THE TRANSFORMATION
C      CONNECTING THE NEW ORIENTATION OF THE SAMPLE WITH THE TESTING AXES.

          DO I=1,3
          DO J=1,3
            AUX33(I,J)=0.
            DO K=1,3
              AUX33(I,J)=AUX33(I,J)+ROTMAT(I,K)*AG(K,J,KGX)
            ENDDO
          ENDDO
          ENDDO
          DO I=1,3
          DO J=1,3
            AG(I,J,KGX)=AUX33(I,J)
          ENDDO
          ENDDO

        ENDDO      ! END OF LOOP OVER GRAINS IN THE PHASE

C *** COLUMNS OF 'FIJPH' ARE GRAIN ELLIPSOID AXES EXPRESSED IN 'TA' (which usually
C     coincide with SA).
C     SUBR 'UPDATE_SHAPE' CALCULATES EIGENVECTORS AND EIGENVALUES OF FIJPH.
C     EIGENVECTORS ARE USED INSIDE 'SUBR VPSC' TO ROTATE 'C4SA' INTO 'C4GA' AND
C     CALCULLATE ESHELBY TENSOR.

        DO I=1,3
        DO J=1,3
          AUX33(I,J)=0.
          DO K=1,3
            AUX33(I,J)=AUX33(I,J)+ROTMAT(I,K)*FIJPH(K,J,IPH)
          ENDDO
        ENDDO
        ENDDO
        DO I=1,3
        DO J=1,3
          FIJPH(I,J,IPH)=AUX33(I,J)
        ENDDO
        ENDDO

        IF(ISHAPE(IPH).GE.1) THEN
          DO KGX=NGR(IPH-1)+1,NGR(IPH)
            DO I=1,3
            DO J=1,3
              AUX33(I,J)=0.
              DO K=1,3
                AUX33(I,J)=AUX33(I,J)+ROTMAT(I,K)*FIJGR(K,J,KGX)
              ENDDO
            ENDDO
            ENDDO
            DO I=1,3
            DO J=1,3
              FIJGR(I,J,KGX)=AUX33(I,J)
            ENDDO
            ENDDO
          ENDDO
        ENDIF

        WRITE(10,*)
        WRITE(10,'(''  INSIDE SUBROUTINE TEXTURE_ROTATION'')')
        WRITE(10,*)
        WRITE(10,'(''  EULER ANGLES & ELLIPSOID AXES OF PHASE'',I3,
     #        '' before rotation'')') IPH
        WRITE(10,'(3F10.3,5x,3F10.3)') 
     #        (EULERPH(I,IPH),I=1,3),(AXISPH(0,I,IPH),I=1,3)

        CALL UPDATE_SHAPE (IPH)    ! updates also grain shape if ishape>0

        WRITE(10,*)
        WRITE(10,'(''  EULER ANGLES & ELLIPSOID AXES OF PHASE'',I3,
     #        '' after rotation'')') IPH
        WRITE(10,'(3F10.3,5x,3F10.3)') 
     #        (EULERPH(I,IPH),I=1,3),(AXISPH(0,I,IPH),I=1,3)

      ENDDO      ! END OF LOOP OVER PHASES

      RETURN
      END


C *****************************************************************************
C     SUBROUTINE TWIN_ORIENTATION    --->    VERSION OF 17/OCT/2017
C
C     GIVEN THE SHEAR DIRECTION 'TWB' AND TWIN PLANE NORMAL 'TWN', DEFINES
C     THE ORIENTATION OF THE TWIN RELATED CRYSTAL.
C
C     FOR TYPE I twins: MAKES A ROTATION OF 180 deg AROUND NORMAL TO TWIN
C       PLANE DEFINED BY 'TWN'
C     FOR TYPE II twins: MAKES A ROTATION OF 180 deg AROUND SHEAR DIRECTION
C       DEFINED BY 'TWB'
C     FOR TYPE III twins: MAKES A ROTATION OF 180 deg AROUND THE VECTOR PERPENDICULAR
C        TO BOTH, THE SHEAR DIRECTION AND THE NORMAL TO THE TWIN PLANE
C
C     USES THE RODRIGUES FORMULA FOR THE CASE phi=180 deg
C          --> R=I+sin(phi)*N+(1-cos(phi))*N^2=-I+2*V^2=-I+2*Vi*Vj
C
C     THE MATRIX 'ATWIN' TRANSFORMS FROM TWIN TO PARENT CRYSTAL AXES
C *****************************************************************************

      SUBROUTINE TWIN_ORIENTATION (ITWTYPE,TWB,TWN,ATWIN)

      DIMENSION TWB(3),TWN(3),TWBxN(3),ATWIN(3,3),XID33(3,3)

      xid33=0.
      do i=1,3
        xid33(i,i)=1.
      enddo

C *** ROTATION OF 180 deg AROUND NORMAL TO TWIN PLANE
      IF(ITWTYPE.EQ.1) THEN
        DO I=1,3
        DO J=1,3
          ATWIN(I,J)=2.*TWN(I)*TWN(J)-XID33(I,J)
        ENDDO
        ENDDO
      ENDIF

C *** ROTATION OF 180 deg AROUND SHEAR VECTOR ON TWIN PLANE
      IF(ITWTYPE.EQ.2) THEN
        DO I=1,3
        DO J=1,3
          ATWIN(I,J)=2.*TWB(I)*TWB(J)-XID33(I,J)
        ENDDO
        ENDDO
      ENDIF

C *** ROTATION OF 180 deg AROUND 3rd VECTOR: TWB x TWN
      IF(ITWTYPE.EQ.3) THEN
        TWBxN(1)=TWB(2)*TWN(3)-TWB(3)*TWN(2)
        TWBxN(2)=TWB(3)*TWN(1)-TWB(1)*TWN(3)
        TWBxN(3)=TWB(1)*TWN(2)-TWB(2)*TWN(1)
        DO I=1,3
        DO J=1,3
          ATWIN(I,J)=2.*TWBxN(I)*TWBxN(J)-XID33(I,J)
        ENDDO
        ENDDO
      ENDIF

      RETURN
      END

C *****************************************************************************
C     SUBROUTINE UPDATE_CRSS_DD      -->     VERSION 17/JUL/2022

C     Updates the CRSS taking into account dislocation density evolution.
C     
C     Version of May 2011 implements Capolungo,Beyerlein,Kaschner,Tome (MSEA 2009)  
C     which replaces 'mode densities' used in the original model of Beyerlein 
C     and Tome (IJP 2008) by 'system densities'.
C     See also application to uranium by Knezevic,McCabe,Tome et al (IJP 2013).
C     
C     Hall-Petch-effect algorithms commented out in this version.
C     Dimensioning of DD model parameters valid only for single phase PX.

C     Average shear modulus S_MOD replaced (CT: 2021 12 12) --> Now uses the 
C     anisotropic shear moduli of each mode in the shear plane & along shear
C     direction MU_MODE(KMOD,IPH) calculated from the anisotropic SX elastic 
C     ctes inside SUBR DATA_CRYSTAL. Another option could be to use the average 
C     shear modulus C44ISO corresponding to a random PX.
C *****************************************************************************

      SUBROUTINE UPDATE_CRSS_DD (IOPTION)

      USE VPSC8DIM
      USE DISDEN

      iprint=0   ! controls diagnostic print-out

      BOLTZ=1.380622E-29      ! units of MPa*m^3/K

C ***************************************************************
C *** READS PARAMETERS OF CONSTITUTIVE MODEL
C ***************************************************************  

      IF (IOPTION.EQ.0) THEN

C *** THE 'CHILD' (twin) PHASE HAS INITIAL ZERO WEIGHT, AND CRSS's 
C     HAVE TO BE REDEFINED WHEN THE TWINS ARE CREATED (NOT IMPLEMENTED)

        IPHEL=1      
        IPH  =1      ! algorithm limited to single phase case

        NSLSYSX=NSLSYS(IPHEL)
        NSLMODX=NSLMOD(IPHEL)
        NMODESX=NMODES(IPHEL)

        NRSX=20        ! power to use in 'rate sensitive' law
        DO IS=1,NSYST(IPH)
          NRS(IS,IPH)=NRSX
          DO JS=1,NSYST(IPH)
            HARD(IS,JS,IPH)=1.0
          ENDDO
        ENDDO

        OPEN(unit=59,file='DIS_DEN.OUT',status='unknown')
        WRITE(59,*) '   EPSACU   ',('Rho_for_mod ',kk=1,NSLMODX)
     #         ,'    Rho_deb     Rho_tot'

C *** READS HARDENING PARAMETERS FOR SLIP AND TWINS

          READ(UR1,*) SH_MOD,GRSZE      ! elastic shear modulus is redundant
          READ(UR1,*) CHIfor,CHIsub
          READ(UR1,*) RHO_initW,TEMPref,EDOTref
          DO KSM=1,NSLMODX
            READ(UR1,'(A)') PROSA
            READ(UR1,*) BURG(KSM),ACTENER(KSM)
            READ(UR1,*) K1gener(KSM),DRAG(KSM)
            READ(UR1,*) EDOT0(KSM)
            READ(UR1,*) RHO_INITS(KSM)
            READ(UR1,*) Atau0(KSM),Atau1(KSM),Atau2(KSM)
            READ(UR1,*) Fdeb0(KSM),Fdeb1(KSM),Fdeb2(KSM)
            READ(UR1,*) HPcoef(KSM)
          ENDDO

          IF(NTWMOD(IPHEL).NE.0) THEN
            IBEG=NSLSYSX+1
            DO KTM=NSLMODX+1,NMODESX
              KTMX=KTM-NSLMODX
              READ(UR1,'(A)') PROSA
              READ(UR1,*) BURG(KTM)
              READ(UR1,*) T0nucl(KTM),T1nucl(KTM),T2nucl(KTM)   ! nucletn threshold
              READ(UR1,*) T0prop(KTM),T1prop(KTM),T2prop(KTM)   ! propag threshold

C *** here it should read CTS1(KTMX,KSMO) & CTS2(KTMX,KSMO) giving 
C *** CTWSL=CTS1*exp((EDOT-EDOTref)/CTS2) but CTWSL=CTS1 if EDOT=EDOTref 
              READ(UR1,*) (CTS1(KTMX,KSMO),KSMO=1,NSLMODX)     ! tw-sl coupling     
C              READ(UR1,*) (CTS2(KTMX,KSMO),KSMO=1,NSLMODX)     ! add in input file   

              READ(UR1,*) ISECTWX,ITWINLAW,TWTHRES1X,TWTHRES2X 
              TWTHRES(1,KTMX,IPHEL)=TWTHRES1X
              TWTHRES(2,KTMX,IPHEL)=TWTHRES2X
              IEND=IBEG-1+NSM(KTM,IPHEL)
              ISECTW(IBEG:IEND,IPHEL)=ISECTWX
              IBEG=IEND+1
            ENDDO
          ENDIF

      ENDIF   ! END OF ioption=0

C **************************************************************************
C *** CALCULATES INITIAL CRSS (TEMPERATURE DEPENDENT) AND INITIALIZES 
C *** DISLOCATION DENSITIES IN ALL SYSTEMS OF EACH PHASE IN ALL GRAINS
C **************************************************************************

      IF(IOPTION.EQ.1 .OR. IOPTION.EQ.2) THEN

C------------------------------------------------------------------------
C   TWO RELATIONSHIPS BETWEEN SHEAR MOD AND TEMP FOR ZR (IN MPA*1000=GPA)
C       SH_MOD=(41.04298-0.02939*TEMP+9.51267E-6*TEMP^2.)*1.e+03
C       SH_MOD=(40.06-0.022*TEMP)*1.e+03
C------------------------------------------------------------------------

        IPH=1
        IPHEL=1
        NSLSYSX=NSLSYS(IPHEL)
        NSLMODX=NSLMOD(IPHEL)
        NMODESX=NMODES(IPHEL)

        EDOT=SQRT(2./3.)*VNORM(DBAR,5)
        RHO_TOT_PH(IPHEL)=0.
        RHO_DEB_PH(IPHEL)=0.
        DO KSM=1,NSLMODX
          RHO_FOR_PH(KSM,IPHEL)=0.
C --> EQ 3.17 B&T --> TAU_0(T)
          ATAU(KSM)=Atau0(KSM)+Atau1(KSM)*EXP(-TEMP/Atau2(KSM))  
        ENDDO
 
C --> EQ 3.26 B&T --> TAU_0_beta(T) FOR THE TWIN SYSTEMS
        IF(NTWMOD(IPHEL).NE.0) THEN
          DO KTM=NSLMODX+1,NMODESX
            KTMX=KTM-NSLMODX
            TAU_PROP(KTM)=T0prop(KTM)+
     #                    T1prop(KTM)*EXP(-TEMP/T2prop(KTM))
            TAU_NUCL(KTM)=T0nucl(KTM)+
     #                    T1nucl(KTM)*EXP(-TEMP/T2nucl(KTM))
            CTWSL(KTMX,1:NSLMODX)=CTS1(KTMX,1:NSLMODX)        ! tw-sl coupling
C     #            *EXP((EDOT-EDOTref)/CTS2(KTMX,1:NSLMODX))  ! not read from input file
          ENDDO
        ENDIF

      ENDIF
C **************************************************************************	  
C   ASSIGNS INITIAL RHO & CRSS VALUES FOR SLIP AND TWIN SYSTEMS
C **************************************************************************

      IF(IOPTION.EQ.1) THEN

        DO KKK=NGR(IPH-1)+1,NGR(IPH)

          ISM=0
          RHO_FOR(KKK)=0.
          RHO_DEB(KKK)=RHO_initW
          DO KMO=1,NSLMODX
            SH_MOD=MU_MODE(KMO,IPH)
            DO IS=1,NSM(KMO,IPHEL)
              ISM=ISM+1
              RHO_S(ISM,KKK)=RHO_INITS(KMO)
              RHO_FOR(KKK)  =RHO_FOR(KKK) + RHO_S(ISM,KKK)
            ENDDO
          ENDDO

          ISM=0

          DO KMO=1,NSLMODX
            DO IS=1,NSM(KMO,IPHEL)
              ISM=ISM+1
C --> Eq 3.19 & ~3.20 B&T --> TAU_0 + TAU_for + TAU_deb
              TAUE(ISM,KKK)= ATAU(KMO) + 
     #           SH_MOD*CHIfor*BURG(KMO)*SQRT(RHO_FOR(KKK))
              IF(RHO_DEB(KKK).GT. 0.0)
     #           TAUE(ISM,KKK)=TAUE(ISM,KKK)-
     #           CHIsub*SH_MOD*BURG(KMO)*SQRT(RHO_DEB(KKK))*
     #           LOG(BURG(KMO)*SQRT(RHO_DEB(KKK)))
C --> Eq 3.21 B&T; Eq 9 Capolungo et al --> TAU_0 + TAU_for + TAU_deb + TAU_hp
              CRSS(ISM,KKK)=TAUE(ISM,KKK)
     #             + HPcoef(KMO)*SH_MOD*SQRT(BURG(KMO)/GRSZE)

              RHO_FOR_PH(KMO,IPHEL)=RHO_FOR_PH(KMO,IPHEL)+
     #                              WGT(KKK)*RHO_S(ISM,KKK)
            ENDDO            ! end of IS
          ENDDO           ! end of KMO

          RHO_DEB_PH(IPHEL)=RHO_DEB_PH(IPHEL)+wgt(KKK)*RHO_DEB(KKK)

          IF(NTWMOD(IPHEL).NE.0) THEN

            DO KMO=NSLMODX+1,NMODESX
            DO IS=1,NSM(KMO,IPHEL)
              ISM=ISM+1
C --> EQ 3.26 B&T --> TAU_0_beta
              RHOnucl=0.
              TAUE(ISM,KKK)=TAU_PROP(KMO)+
     #           (TAU_NUCL(KMO)-TAU_PROP(KMO)) *EXP(-RHOnucl)
C --> EQ 3.27 B&T --> TAU_0_beta + TAU_hp
              CRSS(ISM,KKK)=TAUE(ISM,KKK)
            ENDDO
            ENDDO

          ENDIF

          ngg=10
          iprx=iprint*(ngg/kkk)*(kkk/ngg)
          if(iprx.ne.0) then
            write(10,'('' initial CRSS for grain # '', i4)') kkk
            write(10,'(12f6.0)') (crss(is,kkk),is=1,nsyst(iphel))
          endif

		ENDDO     ! end of KKK

        DO KMO=1,NSLMODX
          RHO_TOT_PH(IPHEL)=RHO_TOT_PH(IPHEL)+RHO_FOR_PH(KMO,IPHEL)
        ENDDO
          RHO_TOT_PH(IPHEL)=RHO_TOT_PH(IPHEL)+RHO_DEB_PH(IPHEL)

      WRITE(59,'(12e12.4)') EPSACU,(RHO_FOR_PH(KMO,IPHEL),KMO=1,NSLMODX)
     #                      ,RHO_DEB_PH(IPHEL),RHO_TOT_PH(IPHEL)

c      write(59,'(''#10'',e11.3,11x,12e11.3)')            ! only grain #10
c    #          rho_deb(10),(rho_s(KSY,10),KSY=1,NSLSYSX)
c      write(59,'(12e11.3)') (crss(KSY,10),KSY=1,NSLSYSX)   ! only grain #10

      ENDIF      ! END OF IOPTION=1

C **************************************************************************
C *** UPDATES DISLOCATION DENSITIES AND CRSS's OF SLIP & TWIN SYSTEMS
C **************************************************************************

      IF (IOPTION.EQ.2) THEN

        RHO_DEB_PH(IPHEL)=0.
        RHO_TOT_PH(IPHEL)=0.
        DO KMO=1,NSLMODX
          RHO_FOR_PH(KMO,IPHEL)=0.
c --> Table 1 in Capolungo et al (2009)
c         FDEB(KMO)=(Fdeb0(KMO)+Fdeb1(KMO)*log(1.+TEMP/Fdeb2(KMO)))
c    #              *Fdeb3
c --> Eq 11 Knezevic et al (2013)
          FDEB(KMO)=Fdeb0(KMO)+Fdeb1(KMO)*exp((TEMPref-TEMP)/Fdeb2(KMO))              
c --> Eq 3-12 B&T; Eq 10 Capolungo et al
          K2recom(KMO)=CHIfor*BURG(KMO)/ACTENER(KMO)
     #           *(1.-BOLTZ*TEMP/(BURG(KMO)**3)/DRAG(KMO)          ! DRAG [MPa]
     #           *LOG(EDOT/EDOT0(KMO)))                 ! BOLTZ/BURG**3 [MPa/K]
          RHOSAT(KMO)=K2recom(KMO)**(-2.0)
          K2recom(KMO)=K1gener(KMO)*K2recom(KMO)
        ENDDO

        DO KKK=NGR(IPH-1)+1,NGR(IPH)

          RHO_FOR(KKK)=0.
          RHO_DEBx=RHO_DEB(KKK)
          ISM=0
          DO KMO=1,NSLMODX
          DO ISY=1,NSM(KMO,IPHEL)
            ISM=ISM+1
c --> Eq 9 Capolungo et al
            DELRHOS(ISM)=K1gener(KMO)*SQRT(RHO_S(ISM,KKK))-
     #           K2recom(KMO)*    (RHO_S(ISM,KKK))           ! calculate with previous
c --> Eq 3-14 B&T; Eq 11 Capolungo et al
            DELRHOW(ISM)=Fdeb(KMO)*K2recom(KMO)*RHO_S(ISM,KKK)*    
     #                    BURG(KMO)*SQRT(RHO_DEBx)           ! calculate with previous
            DELGAM(ISM)=ABS(GAMDOT(ISM,KKK))*TIME_INCR
            RHO_S(ISM,KKK)=RHO_S(ISM,KKK)+DELRHOS(ISM)*DELGAM(ISM)  ! update in system
            RHO_FOR(KKK)  =RHO_FOR(KKK)  +RHO_S(ISM,KKK)
            RHO_DEB(KKK)  =RHO_DEB(KKK)  +DELRHOW(ISM)*DELGAM(ISM)  ! update in grain
            RHO_FOR_PH(KMO,IPHEL)=RHO_FOR_PH(KMO,IPHEL)+
     #                            WGT(KKK)*RHO_S(ISM,KKK)     ! update in mode & phase
		  ENDDO      ! end of ISY=1,NSM
          ENDDO      ! end of KMO=1,NSLMODX

          RHO_DEB_PH(IPHEL)=RHO_DEB_PH(IPHEL)+wgt(KKK)*RHO_DEB(KKK)  ! update in phase

C *** UPDATES crss FOR THE SLIP SYSTEMS USING RECALCULATED DISLOC DENSITIES

          ISM=0
          DO KMO=1,NSLMODX
            SH_MOD=MU_MODE(KMO,IPH)
          DO ISY=1,NSM(KMO,IPHEL)
            ISM=ISM+1
c --> Eqs 3-19 & 3-20 B&T ; Eqs. 7 & 8 Capolungo et al 
            TAUE(ISM,KKK)=ATAU(KMO) +CHIfor*BURG(KMO)*SH_MOD*   
     #          SQRT(RHO_FOR(KKK))  -CHIsub*BURG(KMO)*SH_MOD*               
     #          SQRT(RHO_DEB(KKK))*LOG(BURG(KMO)*SQRT(RHO_DEB(KKK)))
C --> Eq 3.21 B&T; Eq 9 Capolungo et al --> TAU_0 + TAU_for + TAU_deb + TAU_hp
              CRSS(ISM,KKK)=TAUE(ISM,KKK)
     #             + HPcoef(KMO)*SH_MOD*SQRT(BURG(KMO)/GRSZE) 
          ENDDO      ! end of ISY=1,NSM
          ENDDO      ! end of KMO=1,NSLMODX

C *** UPDATES crss FOR THE TWINNING SYSTEMS (COUPLED WITH DISL DENS)

          IF(NTWMOD(IPHEL).NE.0) THEN

          DO IMO=NSLMODX+1,NMODESX
            IMOX=IMO-NSLMODX
          DO IS =1,NSM(IMO,IPHEL)
            ISM=ISM+1
            TAUE(ISM,KKK) =0.0
            RHOnucl=0.
            JSM=0
            DO JMO=1,NSLMODX
              SH_MOD=MU_MODE(JMO,IPH)
            DO JS =1,NSM(JMO,IPHEL)
              JSM=JSM+1
c --> Eq 14 Capolungo et al
              TAUE(ISM,KKK)=TAUE(ISM,KKK)+SH_MOD*BURG(IMO)*
     #           BURG(JMO)*CTWSL(IMOX,JMO)*RHO_S(JSM,KKK)
              RHOnucl=RHOnucl+RHO_S(JSM,KKK)/RHOSAT(JMO)         
            ENDDO
            ENDDO

C --> Eq 3-26 B&T; Eq 15 Capolungo et al
            TAUE(ISM,KKK)=TAUE(ISM,KKK)+TAU_PROP(IMO)+
     #        (TAU_NUCL(IMO)-TAU_PROP(IMO)) *exp(-RHOnucl)       
C --> EQ 3.27 B&T --> TAU_0_beta + TAU_hp
            CRSS(ISM,KKK)=TAUE(ISM,KKK)
          ENDDO
          ENDDO

          ENDIF   ! end of NTWMOD.ne.0

      ngg=10
      iprx=iprint*(ngg/kkk)*(kkk/ngg)
      if(iprx.ne.0) then
        write(10,'('' CRSS for grain  '',i5,6x,3f6.0,/,(12f6.0))')
     #        kkk,(crss(is,kkk),is=1,nsyst(iphel))
        gamdnorm=vnorm(dbar,5)
        write(10,'('' GAMD(is,kkk)'',i5,6x,3f7.2,/,(12f7.2))')
     #        kkk,(gamdot(is,kkk)/gamdnorm,is=1,nsyst(iphel))
      endif
  
        ENDDO      ! END OF LOOP kkk OVER GRAINS

        DO KMO=1,NSLMODX
          RHO_TOT_PH(IPHEL)=RHO_TOT_PH(IPHEL)+RHO_FOR_PH(KMO,IPHEL)
        ENDDO
        RHO_TOT_PH(IPHEL)=RHO_TOT_PH(IPHEL)+RHO_DEB_PH(IPHEL)

      WRITE(59,'(12e12.4)') EPSACU,(RHO_FOR_PH(KMO,IPHEL),KMO=1,NSLMODX)
     #                        ,RHO_DEB_PH(IPHEL),RHO_TOT_PH(IPHEL)
c      write(59,'(''#10'',e11.3,11x,12e11.3)')             ! only grain #10
c     #          rho_deb(10),(rho_s(KSY,10),KSY=1,NSLSYSX)
c      write(59,'(12e11.3)') (crss(KSY,10),KSY=1,NSLSYSX)  ! only grain #10

      ENDIF      ! END OF IF(IOPTION.EQ.2)
C ************************************************************************

      RETURN
      END

C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     SUBROUTINE UPDATE_CRSS_DD_reversal   -->   VERSION 18/AUG/2022
C
C *** Hardening model based on forward and reversible dislocation
c        evolution on individual slip systems
C *** coded by k.kitayama (July 2011)  
C     --> K. Kitayama et al,Intl J of Plasticity 46 (2013) 54-69
C *** modified by CNT for revised constitutive law (July 2012)
C *** Modified by WW for Mg and cycling of LC steel (2013,2014)
C     --> W. Wen et al, Intl J of Plasticity 73 (2015) 171-183
C *** Modified by CNT for adapting to latest VPSC7d version (May 2015)
C *** Back-stress law updated by CNT: rho_total replaced by 
C        rho_s_total (01/21/2016)
C     --> W. Wen et al, shear cycling (Acta Mater, 2016)
C
C if IOPTION=0: reads parameters, calculates elastic shear moduli.
C if IOPTION=1: initializes TAU & RHO using input parameters.
C if IOPTION=2: updates disloc densities RHO and CRSS in every system

C *** Dimensioning of hardening parameters valid only for single phase PX

C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE UPDATE_CRSS_DD_REV (ioption)

      USE VPSC8DIM
      USE DISDEN_REV

      DIMENSION HSELFX(NMODMX),HLATEX(NMODMX,NMODMX),CTSX(NMODMX,NMODMX)
      REAL(8),parameter :: gamd_tiny =1.d-6
      REAL(8),parameter :: rho_tiny  =1.d7

      character(len=100) :: memo

      IPRINT=0

c-----------------------------------------------------------------------
c Reads hardening law parameters
c-----------------------------------------------------------------------

      if(ioption.eq.0) then

        open(500,file='DIS_DEN.OUT',status='replace')

!        open(501,file='srho_f.out',status='replace')
!        open(502,file='srho_r_p.out',status='replace')
!        open(503,file='srho_r_n.out',status='replace')
!        open(504,file='tauback.out',status='replace')
!        open(555,file='P_slip.out',status='replace')
!        open(556,file='D_slip.out',status='replace')
!        open(505,file='rho_r2.txt',status='replace')
!        open(506,file='TAUE.txt',status='replace')
!        open(507,file='dev_anih.txt',status='replace')
!        open(510,file='develop.txt',status='replace')

        write(500,'(4x,''EPSACU'',5x,''EPSVM'',3x,''EPSctrl''
     #             ,4x,''RHOforest'',5x,''RHOreverse'',6X,''RHOtotal''
     #             ,7x,''RHOmod1'',7x''RHOmod2''  )')
!     #             ,3x,''dREV_stored'',3x,''dREV_annihi''
!     #             ,10x,''Dave'',10x,''Pave''
!     #             ,4x,''r_for_mod1'',4x,''r_rev_mod1''
!     #             ,4x,''r_for_mod2'',4x,''r_rev_mod2'')')

!        if (labelwrite.eq.1) lun=501
!        if (labelwrite.eq.2) lun=502
!        if (labelwrite.eq.3) lun=503
!        if (labelwrite.eq.4) lun=504
!        if (labelwrite.eq.5) lun=555
!        if (labelwrite.eq.6) lun=556
!        write(lun,'(4X,''EPSACU'',8x,
!	  #          ''slip_1'',8x,''slip_2'',8x!,
!     #          ''slip_3'',8x,''slip_4'',8x!,''slip_5'',8x,
!     #          ''slip_6'',8x,''slip_7'',8x!,''slip_8'',8x,
!     #          ''slip_9'',8x,''slip_10'',8x,''slip_11''
!     #          ,8x,''slip_12'',8x,''slip_13'',8x,''slip_14''
!     #          ,8x,''slip_15'',7x,''slip_16'',7x,''slip_17''
!     #          ,7x,''slip_18'',8x,''slip_19'',8x,''slip_20''
!     #          ,8x,''slip_21'',8x,''slip_22'',8x,''slip_23''
!     #          ,8x,''slip_24'')')

c      open (511,file='p_evolution.out',status='replace')
c      write(511,'(4X,''EPSACU '',4x,''EPSVM '',10X,''P'',10X,
c     #                 ''Pdeviat'')')
c
c      write(501,'(4X,''EPSACU'',5x,''EPSVM'',5x,''EPS12'',6x,
c     #              ''FORWARD'',8x,''REVERSE'',8X,''TOTAL'',10x,
c     #              ''dREV_stored'',4X,''dREV_annihi'',4x,
c     #              ''Pave'',11x,''Pdeviat'')')


C     Dimensioning of hardening parameters valid only for single phase PX
        IPHEL=1
        IPH  =1

        NSYSTX =NSYST(IPHEL)
        NSLSYSX=NSLSYS(IPHEL)
        NTWSYSX=NTWSYS(IPHEL)

        NMODESX=NMODES(IPHEL)
        NSLMODX=NSLMOD(IPHEL)
        NTWMODX=NTWMOD(IPHEL)

        ISECTWX=0

        NRSX=20          ! power to use in 'rate sensitive' law
        DO IS=1,NSYST(IPH)
          NRS(IS,IPH)=NRSX
        ENDDO
		
        read(1,*) (fBSm(JM),JM=1,NSLMOD(IPH))   ! back stress parameter
        read(1,*) qBS,mREV   ! back stress & reversal laws power
        read(1,*) xD         ! D [m] grain size
          xlambda(:)=xD      ! same initial MFP in all grains
		  
        DO im=1,nslmodx

          READ(1,*) memo
          READ(1,*) xMu(im)     ! mu [MPa]
          READ(1,*) xB(im)      ! b [m]
          READ(1,*) tau0(im)    ! t0 at 20 deg
          READ(1,*) xK(im)      ! K at 20 deg
          READ(1,*) xf_self(im) ! f_self
          READ(1,*) xf_lat(im)  ! f_lat

          READ(1,*) xRho_min(im)
          READ(1,*) xRho_max(im)
            xrho_infl(im)=(xRho_max(im)+xRho_min(im))/2.d0

          READ(1,*)  HSELFX(im)
          READ(1,*) (HLATEX(im,jm),jm=1,nslmodx)
        ENDDO	

        IF(NTWMODX.NE.0) then
        
        DO itm=1,ntwmodx
          READ(1,*) memo
          READ(1,*) xMu(nslmodx+itm)       ! mu [MPa]
          READ(1,*) xB(nslmodx+itm)        ! b [m]
          READ(1,*) tau0(nslmodx+itm)      ! t0 at 20 deg.
          READ(1,*) xK(nslmodx+itm)        ! K at 20 deg.
          READ(1,*) ISECTWX,ITWINLAW,TWTHRES1X,TWTHRES2X
          TWTHRES(1,itm,IPH)=TWTHRES1X
          TWTHRES(2,itm,IPH)=TWTHRES2X
cnt --> ignores ISECTWX & hardwires 'non-seconday twinning allowed'         
          ISECTW(1:nsystx,iph)=0
CNT --> READS TWIN-SLIP HARDENING COUPLING COEFFICIENTS (2022 08 18)
          READ(1,*) (CTSX(itm,ism),ism=1,nslmodx)
        ENDDO

        ENDIF
 
C *** DEFINES SLIP-SLIP HARDENING MATRIX
        ISYS=1
        DO IM=1,NMODES(IPH)
        DO IS=1,NSM(IM,IPH)
          JSYS=1
          DO JM=1,NMODES(IPH)
          DO JS=1,NSM(JM,IPH)
            HARD(ISYS,JSYS,IPH)=HLATEX(IM,JM)
            JSYS=JSYS+1
          ENDDO
          ENDDO
          HARD(ISYS,ISYS,IPH)=HSELFX(IM)
          ISYS=ISYS+1
        ENDDO      ! end of DO IS
        ENDDO      ! end of DO IM

C *** ASSIGNS BACK-STRESS COEFFICIENTS TO SLIP SYSTEMS (2022 08 18)
        islips=0
        DO KMO=1,NSLMODX
          DO ISY=1,NSM(KMO,IPHEL)
            islips=islips+1
            fBSs(islips)=fBSm(kmo)
          ENDDO
        ENDDO

C *** READS TWIN-SLIP HARDENING COUPLING COEFFICIENTS (replaced 2022 08 18)
c        if(ntwsysx.ne.0) then
c          read(1,*) memo
c          do it=nslsysx+1,nsyst(iphel)
c            read(1,*) (CTS(it,is),is=1,nslsysx)
c          enddo
c        endif

C *** ASSIGNS TWIN-SLIP HARDENING COUPLING COEFFICIENTS (2022 08 18)
        if(ntwsysx.ne.0) then
          itwins=0
          do itm=1,ntwmodx
          do its=1,nsm(itm,iph)
            itwins=itwins+1
            islips=0
            do ism=1,nslmodx
            do iss=1,nsm(ism,iph)
              islips=islips+1
              CTS(itwins,islips)=CTSX(itm,ism)
            enddo
           enddo
          enddo
          enddo
        endif

      ENDIF         ! END OF IOPTION.EQ.0

c------------------------------------------------------------------------
c *** Initialization of CRSS tau for the first step
c------------------------------------------------------------------------

      if(ioption.eq.1) then

        rho_f(:,:)    =rho_tiny
        rho_r(:,:)    =rho_tiny
        rho_prev(:,:) =rho_tiny
        rho_deb(:,:)  =rho_tiny
        shear_tag(:,:)=0

        IPHEL=1
        IPH  =1

        NSYSTX =NSYST(IPHEL)
        NSLSYSX=NSLSYS(IPHEL)
        NTWSYSX=NTWSYS(IPHEL)

        NMODESX=NMODES(IPHEL)
        NSLMODX=NSLMOD(IPHEL)
        NTWMODX=NTWMOD(IPHEL)      ! CNT 2013

        XARG= -3.d0*(xRho_max(1)+xRho_min(1))/(xRho_max(1)-xRho_min(1))

        DO IS=1,NSLSYSX
        DO ISP=1,NSLSYSX

        AMSS(IS,ISP)=1-ABS(DNCA(1,IS,IPH)*DNCA(1,ISP,IPH)+DNCA(2,IS,IPH)
     #    *DNCA(2,ISP,IPH)+DNCA(3,IS,IPH)*DNCA(3,ISP,IPH))       ! 2/12/2014  

        IF (ABS(AMSS(IS,ISP)).lt.0.00001) AMSS(IS,ISP)=0.d0
       
        ENDDO
        ENDDO

        DO KKK=NGR(IPH-1)+1,NGR(IPH)
             
          ISM=0
          DO KMO=1,NMODES(IPH)
            DO IS=1,NSM(KMO,IPHEL)
              ISM=ISM+1
	          XP(ISM,KKK)=0.5d0-0.5d0 * tanh(XARG)  !2/12/2014
              TAUE(ISM,KKK)=tau0(KMO)
              tauini(ISM)=tau0(KMO)
              xBs(ism)=xB(kmo)
              CRSS(ISM,KKK)=tau0(KMO)
            ENDDO           ! end of IS
          ENDDO         ! end of KMO
        ENDDO       ! end of KKK

      ENDIF     ! END OF IOPTION=1

c-------------------------------------------------------------------------
c *** Evolves the dislocation densities
c-------------------------------------------------------------------------

      IF (IOPTION.EQ.2) THEN

        IPH=1
        IPHEL=1

        NSYSTX =NSYST(IPHEL)
        NSLSYSX=NSLSYS(IPHEL)
        NTWSYSX=NTWSYS(IPHEL)

        NMODESX=NMODES(IPHEL)
        NSLMODX=NSLMOD(IPHEL)
        NTWMODX=NTWMOD(IPHEL)      ! CNT 2013

! Dislocation densities and reversibility parameter p averaged over the grains
        rho_f_ave  =0.d0
        rho_r_ave  =0.d0
        rho_tot_ave=0.d0


        do im=1,nslmodx       ! set rhos for each mode
          rho_mod(im) =0.0d0
          for_rhos(im)=0.0d0
          rev_rhos(im)=0.0d0
        enddo

        drho_r_stored_ave=0.0
        drho_r_annihi_ave=0.0
        p_ave=0.d0
        p_std=0.d0
        debs_ave=0.0d0

        KGX=1
        DO IPH=IPHBOT,IPHTOP
          IPHEL=IPH-IPHBOT+1
          DO KKK=NGR(IPH-1)+1,NGR(IPH)

	    do is=1,nslmodx   ! set rhos for each mode
	      rhos(is)=0.d0
	      for_rhosl(is)=0.0d0
	      rev_rhosl(is)=0.0d0
        enddo

            rho_f_g  =0.d0
            rho_r_g  =0.d0
            rho_tot_g=0.d0
            drho_r_stored  =0.d0
            drho_r_annihi  =0.d0

C *** DROPS SYSTEMS WITH NEGLIGIBLE SHEAR RATES
            GAMDMAX=0.D0
            DO IS=1,NSYST(IPHEL)
              IF(ABS(GAMDOT(IS,KGX)).GT.GAMDMAX) GAMDMAX=GAMDOT(IS,KGX)
            ENDDO
            D_GAMG=0.0
            DO IS=1,NSYST(IPHEL)
              IF(ABS(GAMDOT(IS,KGX)).LT.0.001*GAMDMAX) GAMDOT(IS,KGX)=0.
              D_GAMG=D_GAMG+ABS(GAMDOT(IS,KGX))*TIME_INCR
            ENDDO

            ISM=0
            DO KMO=1,NSLMODX
              DO ISY=1,NSM(KMO,IPHEL)
                ISM=ISM+1

                bL=xB(KMO)*xlambda(KGX)

                IF(GAMDOT(ISM,KGX) > gamd_tiny) THEN

                  D_GAMS=GAMDOT(ISM,KGX)*TIME_INCR

C *** RESETS THE SHEAR TAG IF THE SENSE IS INVERTED WRT PREVIOUS STEP

                  ISTAG=SHEAR_TAG(ISM,KGX)
                  IF(ISTAG.EQ.0 .OR. ISTAG.EQ.-1) THEN
                    rho_prev(ism,kgx)=rho_f(ism,kgx)
     1                   +rho_r(2*ism-1,kgx)+rho_r(2*ism ,kgx)
                    SHEAR_TAG(ISM,KGX)=1
                  ENDIF

                  drho_f  = (1.d0-xp(ISM,kgx))*D_GAMS/bL
     1                      -xf_self(KMO)*rho_f(ISM,KGX)*D_GAMG
                  drho_r_p= xp(ISM,kgx)*D_GAMS/bL
     1                  -xf_self(KMO)*rho_r(2*ISM-1,KGX)*D_GAMG
                  drho_r_m= -D_GAMS/bL
     1                      *(rho_r(2*ISM,KGX)/rho_prev(ism,kgx))**mREV

                  drho_r_stored=drho_r_stored+drho_r_p
                  drho_r_annihi=drho_r_annihi+drho_r_m

                ELSE IF(GAMDOT(ISM,KGX) .lt. -gamd_tiny) THEN

                  D_GAMS=abs(GAMDOT(ISM,KGX))*TIME_INCR

C *** RESETS THE SHEAR TAG IF THE SENSE IS INVERTED WRT PREVIOUS STEP

                  ISTAG=SHEAR_TAG(ISM,KGX)
                  IF(ISTAG.EQ.0 .OR. ISTAG.EQ.1) THEN
                    rho_prev(ism,kgx)=rho_f(ism,kgx)
     1                   +rho_r(2*ism-1,kgx)+rho_r(2*ism ,kgx)
                    SHEAR_TAG(ISM,KGX)=-1
                  ENDIF

                  drho_f  = (1.d0-xp(ISM,kgx))*D_GAMS/bL
     1                      -xf_self(KMO)*rho_f(ISM,KGX)*D_GAMG
                  drho_r_p= -D_GAMS/bL
     1                     *(rho_r(2*ISM-1,KGX)/rho_prev(ism,kgx))**mREV
                  drho_r_m= xp(ISM,kgx)*D_GAMS/bL
     1                    -xf_self(KMO)*rho_r(2*ISM,KGX)*D_GAMG

                  drho_r_stored=drho_r_stored+drho_r_m
                  drho_r_annihi=drho_r_annihi+drho_r_p

                ELSE

C *** RESETS THE SHEAR TAG IF THE SENSE CHANGED WRT PREVIOUS STEP

                  ISTAG=SHEAR_TAG(ISM,KGX)
                  IF(ISTAG.EQ.1 .OR. ISTAG.EQ.-1) THEN
                    SHEAR_TAG(ISM,KGX)=0
                  ENDIF

                  drho_f=   -xf_lat(KMO)*rho_f(ism,kgx)    *D_GAMG
                  drho_r_p= -xf_lat(kmo)*rho_r(2*ism-1,kgx)*D_GAMG
                  drho_r_m= -xf_lat(kmo)*rho_r(2*ism,kgx)  *D_GAMG

                END IF
                
                rho_ps(ism,kgx)=rho_f(ism,kgx)
     1                   +rho_r(2*ism-1,kgx)+rho_r(2*ism ,kgx)

                rho_f(ism,kgx)     = rho_f(ism,kgx)     + drho_f
                rho_r(2*ism-1,kgx) = rho_r(2*ism-1,kgx) + drho_r_p
                rho_r(2*ism,kgx)   = rho_r(2*ism,kgx)   + drho_r_m

                if(rho_f(ism,kgx).lt.0.d0)     rho_f(ism,kgx)    =0.d0
                if(rho_r(2*ism-1,kgx).lt.0.d0) rho_r(2*ism-1,kgx)=0.d0
                if(rho_r(2*ism,kgx).lt.0.d0)   rho_r(2*ism,kgx)  =0.d0
				
                drho_tot(ISM,KGX)=rho_f(ism,kgx)+rho_r(2*ism-1,kgx)
     1                   +rho_r(2*ism ,kgx)-rho_ps(ism,kgx)        ! 2/12/2014

c        if(kgx.eq.1)then
c          write(*,*)ism,rho_ps(ism,kgx),drho_tot(ISM,KGX)
c          pause
c        endif
              
		   rhos(kmo)=rhos(kmo)+rho_f(ism,kgx)+
     #               rho_r(2*ism-1,kgx) + rho_r(2*ism,kgx)
	       for_rhosl(kmo)=for_rhosl(kmo)+rho_f(ism,kgx)
           rev_rhosl(kmo)=rev_rhosl(kmo)+
     #               rho_r(2*ism-1,kgx) + rho_r(2*ism,kgx)

         rho_f_g =rho_f_g + rho_f(ism,kgx)
         rho_r_g =rho_r_g + rho_r(2*ism-1,kgx) + rho_r(2*ism,kgx)

        end do      ! end of do isy=1,nsm
        end do        ! end of do kmo=1,nslmodx

C *** UPDATES crss of SLIP SYSTEMS USING RECALCULATED DISLOC DENSITIES

          RHO_TOT_G=RHO_F_G + RHO_R_G

          DO IS=1,NSLSYSX          ! CNT
	        RHO_EFF=0.0d0
	        DO isp=1,NSLSYSX      
                RHO_EFF=RHO_EFF+HARD(is,isp,IPH)
     #          *(rho_f(isp,kgx)+rho_r(2*isp-1,kgx)+rho_r(2*isp,kgx))
            ENDDO

            if(RHO_EFF < 0.d0) RHO_EFF=0.d0   ! avoids negative sqrt

            TAUE(IS,KGX)=tauini(is)+xMu(1)*xBs(is)*sqrt(RHO_EFF)

            if(gamdot(is,kgx).gt.gamd_tiny) then
			
c *** 2015 version of the back stress
c           CRSS(is,kgx)=TAUE(is,kgx)-fBSs(is)*TAUE(is,kgx)*
c     #                     (RHO_R(2*IS,kgx)/RHO_TOT_G)**qBS
c *** 2016 version of the back stress
          CRSS(is,kgx)=TAUE(is,kgx)-fBSs(is)*TAUE(is,kgx)*      ! ww 01/21/2016
     #    (RHO_R(2*IS,kgx)/
     #    (rho_f(is,kgx)+rho_r(2*is-1,kgx)+rho_r(2*is,kgx)))**qBS

            else if(gamdot(is,kgx).lt.-gamd_tiny)then

c *** 2015 version of the back stress			
c           CRSS(is,kgx)=TAUE(is,kgx)-fBSs(is)*TAUE(is,kgx)*  
c     #                     (RHO_R(2*IS-1,kgx)/RHO_TOT_G)**qBS
c *** 2016 version of the back stress	 
           CRSS(is,kgx)=TAUE(is,kgx)-fBSs(is)*TAUE(is,kgx)*      ! ww 01/21/2016
     #     (RHO_R(2*IS-1,kgx)/
     #     (rho_f(is,kgx)+rho_r(2*is-1,kgx)+rho_r(2*is,kgx)))**qBS

            else
              CRSS(is,kgx)=TAUE(is,kgx)
            endif
            taud(is)=CRSS(is,1)-TAUE(is,1)

c       if(abs(crss(is,kgx)/TAUE(is,kgx)).lt.0.1) then
c		write(*,*) 'is=',is, '  kgx=',kgx,
c     #    '  crss=',crss(is,kgx),'  TAUE=', Taue(is,kgx)
c       endif
			
        ENDDO

C *** UPDATES CRSS of TWINNING SYSTEMS USING Beyerlein et al equation

          IF(NTWSYSX.GT.0) THEN

C *** THIS ALGORITHM (ESPECIALLY THE rand PART) NEEDS TO BE CHECKED 
            DO IT=NSLSYSX+1,NSYSTX
              tautwint(it)=0.d0
              TAUE(IT,KGX)=tauini(IT)
              DO IS=1,NSYSTX
                tautwint(it)=tautwint(it)+
     #             xMu(1)*xBs(it)*xBs(is)*cts(it,is)*      ! assume same mu for all modes
     #             (rho_f(is,kgx)+rho_r(2*is-1,kgx)+rho_r(2*is,kgx))
              ENDDO
              rand=2*random2(JRAN)-1.0d0        ! rand=[-1,1] 
              factorx=0.	                    ! no statistical deviation
              CRSS(IT,KGX)=TAUE(IT,KGX) + rand*TAUE(IT,KGX)*factorx	 
            ENDDO

          ENDIF

C *** UPDATES MEAN FREE PATH IN THE GRAIN USING RECALCULATED DISLOC DENSITIES
          xlambda(kgx)= 1.d0/(sqrt(rho_tot_g)/xK(1)+1.d0/xD)

C *** UPDATES THE DEBRIS DISLOCATION DENSITY and THE REVERSIBILITY FACTOR P

          DO IS=1,NSLSYS(IPHEL)
            deltaD(IS,kgx)=0.d0
            DO ISP=1,NSLSYS(IPHEL)
              deltaD(IS,kgx)=deltaD(IS,kgx)+(rho_f(isp,kgx)
     #      +rho_r(2*isp-1,kgx)+rho_r(2*isp,kgx))*AMSS(IS,ISP)   ! 2/12/2014
            ENDDO
            deltaD(IS,kgx)=deltaD(IS,kgx)*drho_tot(IS,KGX)
            
C       rho_deb is no longer debris. It is D now!

            rho_deb(IS,kgx)=rho_deb(IS,kgx) + deltaD(IS,kgx)     ! 2/12/2014
            if(rho_deb(IS,kgx).lt.0.0d0)then
              rho_deb(IS,kgx)=0.01d0              !! 09/12/2014 ww 
            endif

            xp(IS,KGX)=0.5d0-0.5d0 * tanh
     #      (3.*((xRho_max(1)+xRho_min(1))/(xRho_max(1)-xRho_min(1)))
     #          *((sqrt(rho_deb(IS,kgx))-xrho_infl(1))/xrho_infl(1)))
             
          ENDDO

C *** PERFORMS SOME STATISTICS FOR OUTPUT PURPOSES
            xp_g(kgx)=0.0d0
            debs(kgx)=0.0d0

            DO IS=1, NSLSYS(IPHEL)	
              xp_g(kgx)=xp_g(kgx)+xp(IS,kgx)
              debs(kgx)=debs(kgx)+rho_deb(IS,kgx)
	        ENDDO

            xp_g(kgx)=xp_g(kgx)/NSLSYS(IPHEL)
            debs(kgx)=debs(kgx)/NSLSYS(IPHEL)

            p_ave=p_ave+xp_g(KGX)   *wgt(kgx)
            p_std=p_std+xp_g(KGX)**2*wgt(kgx) 
            debs_ave=debs_ave+debs(kgx)*wgt(kgx)   !ww 09/12/2014

            RHO_F_AVE=RHO_F_AVE+RHO_F_G      *WGT(KGX)
            RHO_R_AVE=RHO_R_AVE+RHO_R_G      *WGT(KGX)
            RHO_TOT_AVE=RHO_TOT_AVE+RHO_TOT_G*WGT(KGX)

            do is=1,NSLMODX
              rho_mod(is)=rho_mod(is)+rhos(is) *WGT(KGX)
              for_rhos(is)=for_rhos(is)+for_rhosl(is)*wgt(kgx)   ! uapt 3/12/2013
	          rev_rhos(is)=rev_rhos(is)+rev_rhosl(is)*wgt(kgx)   ! uapt 3/12/2013
            enddo

            drho_r_stored_ave=drho_r_stored_ave+drho_r_stored*wgt(kgx)
            drho_r_annihi_ave=drho_r_annihi_ave+drho_r_annihi*wgt(kgx)

            ngg=1
            iprx=iprint*(ngg/kkk)*(kkk/ngg)
            if(iprx.ne.0) then
              write(10,'('' CRSS for grain  '',i5,6x,3f6.0,/,(12f6.0))')
     #          kkk,(crss(is,kkk),is=1,nsyst(iphel))
                gamdnorm=vnorm(dbar,5)
              write(10,'('' GAMD(is,kkk)'',i5,6x,3f7.2,/,(12f7.2))')
     #          kkk,(gamdot(is,kkk)/gamdnorm,is=1,nsyst(iphel))
            endif

            KGX=KGX+1

          ENDDO      ! END OF LOOP kkk OVER GRAINS
        END DO ! END OF LOOP IPH OVER PHASES

        p_std=sqrt(abs(p_std-p_ave**2))      ! avoid negative sqrt
			
        write(500,'(3f10.4,16E14.4)') EPSACU, EPSVM, EPSTOTv(ictrl)
     #              ,RHO_F_AVE, RHO_R_AVE, RHO_TOT_AVE
     #              ,(rho_mod(i),i=1,nslmod(iphel))
!     #              ,drho_r_stored_ave,drho_r_annihi_ave
!     #              ,debs_ave,P_ave
!     #              ,(for_rhos(i),rev_rhos(i),i=1,nslmod(iphel))

!      NOG=1
!      write(501,'(1f10.4,24E15.6)') EPSACU
!     #   ,(rho_f(is,NOG),    is=1,nsyst(iphel))
!      write(502,'(1f10.4,24E15.6)') EPSACU
!     #   ,(rho_r(2*is-1,NOG),is=1,nsyst(iphel))
!      write(503,'(1f10.4,24E15.6)') EPSACU
!     #   ,(rho_r(2*is,NOG),  is=1,nsyst(iphel))
!      write(504,'(1f10.4,24E15.6)')  EPSACU
!     #   ,(taud(is),         is=1,nsyst(iphel))
!      write(555,'(1f10.4,24E15.6)') EPSACU
!     #   ,(xp(IS,NOG),      is=1,NSLSYS(IPHEL))
!      write(556,'(1f10.4,24E15.6)') EPSACU
!     #   ,(rho_deb(IS,NOG), is=1,NSLSYS(IPHEL))

!      write(505,'(25d20.10)') epsvm,rho_r(2:49:2,1)
!      write(506,'(25d20.10)') epsvm,CRSS(1:24,1)

      ENDIF      ! END OF IF(IOPTION.EQ.2)

      RETURN
      END
C
C **************************************************************************
C     SUBROUTINE UPDATE_CRSS_MTS      --->      VERSION OF 18/JAN/2022

C     IMPLEMENTS A SIMPLE 'MECHANICAL THRESHOLD STRESS' EXPRESSION AT THE
C     SLIP SYSTEM LEVEL.
C     THIS VERSION USES THE VON MISES SCALAR STRAIN RATE 'EDOT' AS THE
C     INTERNAL PARAMETER. 
C     TEMPERATURE AND RATE EFFECTS ARE ACCOUNTED FOR BY THE SCALING FACTORS
C     'Se' AND 'Si', THE SAME FOR ALL THE SYSTEMS.
C     THE RATE DEPENDENCE IS ACCOUNTED FOR BY 'CRSS', AND 'GAMMA_dot_0' IN 
C     THE POWER LAW IS RENORMALIZED TO ELIMINATE THE 'NRS' RATE SENSITIVITY

C     IOPTION=0: READS PARAMETERS OF THE MTS MODEL AND DEFINES SOME
C                HARDENING PARAMETERS
C     IOPTION=1: DEFINES INITIAL RATE INDEPENDENT CRSS 'TAUe' FOR EVERY GRAIN
C     IOPTION=2: RECALCULATES TERMS WHICH DEPEND ON RATE AND TEMPERATURE
C               (BUT WHICH ARE INDEPENDENT OF STRUCTURE EVOLUTION) AND
C                UDPDATES CRSS IN EVERY SYSTEM USING STORED VALUE OF THE
C                STRUCTURE DEPENDENT TERM 'TAUe'
C     IOPTION=3: RECALCULATES TERMS WHICH DEPEND ON RATE AND TEMPERATURE
C               (BUT WHICH ARE INDEPENDENT OF STRUCTURE EVOLUTION).
C                UPDATES STRUCTURE TERM THAT DEPENDS ON DISLOCATION HARDENING.
C                UDPATES THRESHOLD STRESS IN EVERY SYSTEM OF EVERY GRAIN.

C --> Dimensioning of hardening parameters valid only for single phase PX !
*****************************************************************************

      SUBROUTINE UPDATE_CRSS_MTS (EDOT,IOPTION)

      USE VPSC8DIM
      USE MTSSAVE

      IUNIT=UR1      ! READS MTS PARAMETERS AT THE END OF SX FILE.

      IF (IOPTION.EQ.0) THEN

        READ(IUNIT,'(A)') PROSA
        READ(IUNIT,*) KOVB3,MU0,D0,T0
        READ(IUNIT,*) TAUa,TAUi,TAUeini,TH0,KAP
        READ(IUNIT,*) G0i,ED0i,QQi,PPi
        READ(IUNIT,*) G0e,ED0e,QQe,PPe
        READ(IUNIT,*) G0esat,EDesat0,TAUesat0
        READ(IUNIT,*) PSI,RHO
        READ(IUNIT,*) Cp1,Cp2,Cp3

        IPH=IPHBOT     ! for initialization IPHBOT=IPHTOP
        NRSX=20        ! power to use in 'rate sensitive' law
        DO IS=1,NSYST(IPH)
          NRS(IS,IPH)=NRSX
          DO JS=1,NSYST(IPH)
            HARD(IS,JS,IPH)=1.0
          ENDDO
        ENDDO

        RETURN
      ENDIF

      IF(IOPTION.EQ.1) THEN
        IPH=IPHBOT     ! for initialization IPHBOT=IPHTOP
        DO KKK=NGR(IPH-1)+1,NGR(IPH)
          DO IS=1,NSYST(IPH)
            TAUe(IS,KKK)=TAUeini
          ENDDO
        ENDDO

        RETURN
      ENDIF

C*** NEXT PARAMETERS REQUIRED FOR OPTION 2 AND 3 (once TEMP & EDOT are known)

      MU=MU0*(1.-D0/(EXP(T0/TEMP)-1.))
      Cp=Cp1+Cp2*TEMP+Cp3*TEMP**(-2)

      Si=KOVB3*TEMP*LOG(ED0i/EDOT)/G0i/MU
      Si=(1.-Si**(1./QQi))**(1./PPi)
      Se=KOVB3*TEMP*LOG(ED0e/EDOT)/G0e/MU
      Se=(1.-Se**(1./QQe))**(1./PPe)

      TAUesat=TAUesat0*(EDOT/EDesat0)**(KOVB3*TEMP/G0esat/MU)
      TAUini =TAUa + MU/MU0*Si*TAUi + MU/MU0*Se*TAUeini
      HFACTOR=(KAP-1.)*TH0/(TAUesat**KAP)*(MU/MU0)

      IF (IOPTION.EQ.2) THEN
        DO IPH=IPHBOT,IPHTOP
          IPHEL=IPH-IPHBOT+1
          DO KKK=NGR(IPH-1)+1,NGR(IPH)
            DO IS=1,NSYST(IPHEL)
              CRSS(IS,KKK)=TAUini+MU/MU0*Se*TAUe(IS,KKK)
            ENDDO
          ENDDO
        ENDDO

C       WRITE(10,*)
C       WRITE(10,'(''        TEMP        EDOT'')')
C       WRITE(10,'(6E12.5)') TEMP,EDOT
C       WRITE(10,'(''          MU          Cp         Si         Se'')')
C       WRITE(10,'(6E12.5)') MU,Cp,Si,Se
C       WRITE(10,'(''     TAUesat      TAUini    HFACTOR    TFACTOR'')')
C       WRITE(10,'(6E12.5)') TAUesat,TAUini,HFACTOR,TFACTOR
C       WRITE(10,'(''     CRSS IN GRAIN '',i5)') NGR(IPHTOP)
C       WRITE(10,'(12F7.1)') (CRSS(IS,NGR(IPHTOP)),IS=1,NSYST(IPHEL))

      ENDIF       ! END OF IOPTION=2

      IF (IOPTION.EQ.3) THEN
        WORKINCR=0.
        KGX=1
        DO IPH=IPHBOT,IPHTOP
          IPHEL=IPH-IPHBOT+1
          DO KKK=NGR(IPH-1)+1,NGR(IPH)
            DO IS=1,NSYST(IPHEL)
              WORKINCR=WORKINCR+
     #            CRSS(IS,KKK)*ABS(GAMDOT(IS,KGX))*TIME_INCR*WGT(KKK)
              DTAU=0.
              DO JS=1,NSYST(IPHEL)
                DTAU=DTAU+
     #               HARD(IS,JS,IPHEL)*ABS(GAMDOT(JS,KGX))*TIME_INCR
              ENDDO
              IF(TAUesat.LT.TAUe(IS,KKK)) THEN
                TAUesat=TAUe(IS,KKK)
                write(*,'('' TAUesat<TAUe(IS,KKK) for'',2I10)') IS,KKK
                print *, 'enter c to continue'
                read  *
              ENDIF
              TAUe(IS,KKK)=TAUesat-((TAUesat-TAUe(IS,KKK))**(1.-KAP)+
     #                               HFACTOR*DTAU)**(1./(1.-KAP))
              CRSS(IS,KKK)=TAUini+MU/MU0*Se*TAUe(IS,KKK)
            ENDDO
            KGX=KGX+1
          ENDDO
        ENDDO

        TFACTOR=PSI/RHO/Cp      ! ADIABATIC HEATING
        IF(EDOT.LT.500) TFACTOR=0
        TEMP=TEMP+TFACTOR*WORKINCR

      ENDIF      ! END OF IOPTION=3

      RETURN
      END
C
C *****************************************************************************
C     SUBROUTINE UPDATE_CRSS_VOCE     --->    VERSION OF 08/MAR/2022
C
C     A VOCE HARDENING LAW, FUNCTION OF THE ACCUMULATED SHEAR IN EACH GRAIN.
C     THE UNITS AND THE STRENGTH OF THE CRSS ARE CARRIED BY THE FACTOR 'VOCE'.
C     THE SELF & LATENT COUPLING COEFFICIENTS 'HARD' ARE DIMENSIONLESS AND
C     MULTIPLY THE FACTOR 'VOCE'.
C     THE INCREMENT OF CRSS IS GIVEN BY AN ANALYTIC INTEGRAL EXPRESSION OF
C     THE VOCE FUNCTION, RATHER THAN USING FORWARD EXTRAPOLATION, AS WAS
C     DONE BEFORE 28/SET/00.
C     THE PARAMETERS IN VOCE EXPRESSION MAY ADOPT 'NON-KOSHER' VALUES (DEC/00)
C
C     IOPTION=0: READS HARDENING PARAMETERS OF SLIP AND TWIN MODES AND
C                ASSIGNS THOSE TO THE SLIP AND TWIN SYSTEMS
C     IOPTION=1: INITIALIZES THE CRSS OF EACH SYSTEM
C     IOPTION=2: UPDATES THE CRSS OF EACH SYSTEM AS A FUNCTION OF THE SHEAR
C                INCREMENTS TAKING PLACE IN ALL SYSTEMS
*******************************************************************************

      SUBROUTINE UPDATE_CRSS_VOCE (IOPTION,IPH)

      USE VPSC8DIM
      USE CRSS_VOCE
	  
      DIMENSION HSELFX(NMODMX),HLATEX(NMODMX,NMODMX)

      IPRINT=0    ! controls printing

C ***************************************************************************

      IF (IOPTION.EQ.0) THEN

        READ(UR1,*) GRSZE(IPH)
        ISYS=1
        NTWMODX=0
        DO IM=1,NMODES(IPH)

          READ(UR1,*)        ! reminder/identifier of system
          READ(UR1,*) NRSX   ! rate sensitivity parameter

C *** READ EXTENDED VOCE PARAMETERS AND LATENT HARDENING PARAMETERS H_ss'
C     THE LATTER ARE ASSUMED TO BE THE SAME FOR ALL SYSTEMS IN THE MODE
C *** SELF-HARDENING REFERS TO H_ss AND IS ARBITRARILY SET AS REFERENCE
          READ(UR1,*) TAU0X,TAU1X,THET0X,THET1X, HPFACX
          HSELFX(IM)=1.
          READ(UR1,*) (HLATEX(IM,JM),JM=1,NMODES(IPH))

C *** CHECKS WHETHER VOCE PARAMETERS ARE KOSHER:
C         TAU0>0 , TAU1 >= 0 , THET0 >= THET1 >= 0
C         TAU1=0   CORRESPONDS TO LINEAR HARDENING.
C         THETA0=0 FORCES NO-HARDENING.
C *** IF VOCE PARAMETERS ARE NON-KOSHER CHECKS FOR ILL-POSED HARDENING.

          CALL CHECK_VOCE (IM,IPH,TAU0X,TAU1X,THET0X,THET1X)

C *** READS TWIN LAW PARAMETERS FOR TWINNING SYSTEMS
          IF(ISYS.LE.NSLSYS(IPH)) ISECTWX=0
          IF(ISYS.GT.NSLSYS(IPH)) THEN
            READ(UR1,*) ISECTWX,ITWINLAW,TWTHRES1X,TWTHRES2X
            NTWMODX=NTWMODX+1
            TWTHRES(1,NTWMODX,IPH)=TWTHRES1X
            TWTHRES(2,NTWMODX,IPH)=TWTHRES2X
          ENDIF

          DO IS=1,NSM(IM,IPH)
            NRS(ISYS,IPH)   =NRSX
            TAU(ISYS,0,iph) =TAU0X
            TAU(ISYS,1,iph) =TAU1X
            THET(ISYS,0,iph)=THET0X
            THET(ISYS,1,iph)=THET1X

            ISECTW(ISYS,iph)=ISECTWX
            HPFAC(ISYS,iph) =HPFACX

            JSYS=1
            DO JM=1,NMODES(IPH)
            DO JS=1,NSM(JM,IPH)
              HARD(ISYS,JSYS,IPH)=HLATEX(IM,JM)
              JSYS=JSYS+1
            ENDDO
            ENDDO
            HARD(ISYS,ISYS,IPH)=HSELFX(IM)
            ISYS=ISYS+1

          ENDDO      ! end of DO IS
        ENDDO        ! end of DO IM

C     WRITE(10,*)
C     WRITE(10,'(''  HARDENING MATRIX FOR PHASE'',I3)') IPH
C     DO IS=1,NSYST(IPH)
C       WRITE(10,'(24F5.1)') (HARD(IS,JS,IPH),JS=1,NSYST(IPH))
C     ENDDO

      ENDIF   ! END OF IOPTION=0
C ***************************************************************************

      IF (IOPTION.EQ.1) THEN

        DO IPH=IPHBOT,IPHTOP
          IPHEL=IPH-IPHBOT+1
          DO KKK=NGR(IPH-1)+1,NGR(IPH)
            DO IS=1,NSYST(IPHEL)
              CRSS(IS,KKK)=TAU(IS,0,IPHEL)
     #                    +HPFAC(IS,IPHEL)/SQRT(GRSZE(IPH))
            ENDDO

            iprx=iprint*(1/kkk)*(kkk/1)
            if(iprx.ne.0) then
              write(10,'('' initial CRSS for grain # '', i4)') kkk
              write(10,'(9f7.1)') (crss(is,kkk),is=1,nsyst(iphel))
            endif

          ENDDO
        ENDDO

      ENDIF   ! END OF IOPTION=1
C ***************************************************************************

      IF (IOPTION.EQ.2) THEN

        KGX=1
        DO IPH=IPHBOT,IPHTOP
          IPHEL=IPH-IPHBOT+1
        DO KKK=NGR(IPH-1)+1,NGR(IPH)

          GAMTOTX=GTOTGR(KKK)
          DELTGAM=0.0
          DO IS=1,NSYST(IPHEL)
            DELTGAM=DELTGAM+ABS(GAMDOT(IS,KGX))*TIME_INCR
          ENDDO

          DO IS=1,NSYST(IPHEL)
            DTAU=0.
            DO JS=1,NSYST(IPHEL)
              DTAU=DTAU+HARD(IS,JS,IPHEL)*ABS(GAMDOT(JS,KGX))*TIME_INCR
            ENDDO
            TAU0 =TAU (IS,0,IPHEL)
            TAU1 =TAU (IS,1,IPHEL)
            THET0=THET(IS,0,IPHEL)
            THET1=THET(IS,1,IPHEL)
            TINY=1.E-4*TAU0

            VOCE=0.0
            IF(ABS(THET0).GT.TINY) THEN
              VOCE=THET1*DELTGAM
              IF(ABS(TAU1).GT.TINY) THEN
                FACT=ABS(THET0/TAU1)
                EXPINI=EXP(-GAMTOTX*FACT)
                EXPDEL=EXP(-DELTGAM*FACT)
                VOCE  =VOCE-(FACT*TAU1-THET1)/FACT*EXPINI*
     #            (EXPDEL-1.)-THET1/FACT*EXPINI*
     #            (EXPDEL*((GAMTOTX+DELTGAM)*FACT+1.)-(GAMTOTX*FACT+1.))
              ENDIF
            ENDIF
            CRSS(IS,KKK)=CRSS(IS,KKK)+DTAU*VOCE/DELTGAM
          ENDDO

          iprx=iprint*(1/kkk)*(kkk/1)
            if(iprx.ne.0) then
            write(10,'('' CRSS(is,kkk)'',i5,6x,3f6.0,/,(12f6.0))')
     #            kkk,(crss(is,kkk),is=1,nsyst(iphel))
            write(10,'('' GAMD(is,kkk)'',i5,6x,3f7.4,/,(12f7.4))')
     #            kkk,(gamdot(is,kkk),is=1,nsyst(iphel))
            endif

          GTOTGR(KKK)=GAMTOTX+DELTGAM
          KGX=KGX+1
        ENDDO
        ENDDO

      ENDIF   ! END OF IOPTION=2
C ***************************************************************************

      RETURN
      END

C *****************************************************************************
C     SUBROUTINE UPDATE_FIJ      --->      VERSION 11/JAN/2009
C
C     USES THE VELOCITY GRADIENT (AVERAGE, PHASE or GRAIN) IN THE STEP
C     TO UPDATE INCREMENTALLY THE CORRESPONDING DEFORMATION TENSOR 'FIJ'
C
C --> REPLACES PREVIOUS SUBROUTINE UPDFIJ_AVERAGE & SUBR. UPDFIJ_LOCAL  (SEP/04)
C *****************************************************************************

      SUBROUTINE UPDATE_FIJ (IPH)

      USE VPSC8DIM

      DIMENSION FNEW(3,3)


C *** UPDATES THE DEFORM GRAD IN THE ELEMENT 'FIJPH(i,j,0)' USING THE
C     MACROSCOPIC VELOCITY GRADIENT 'LIJBARc'
C *** UPDATES THE DEFORM GRAD IN EACH PHASE 'FIJPH(i,j,IPH)' USING THE
C     AVERAGE VELOCITY GRADIENT FOR THE PHASE 'LIJPH' CALCULATED INSIDE
C     SUBROUTINE UPDATE_ORIENTATION.
C *** LIJPH ACCOUNTS FOR ROTATIONS BUT NOT FOR STRETCH WHEN FREEZE_SHAPE(iph)=1.
C *** FIJPH COINCIDES WITH FIJ OF ELEMENT IF NPH=1 AND FREEZE_SHAPE(1)=0.

      DO I=1,3
      DO J=1,3
        FNEW(I,J)=0.0
        IF(IPH.EQ.0) THEN
          DO K=1,3
            FNEW(I,J)=FNEW(I,J)+(TIME_INCR*LIJBARc(I,K)+XID3(I,K))
     #                         *FIJPH(K,J,0)
          ENDDO
        ELSE IF(IPH.GT.0) THEN
          DO K=1,3
            FNEW(I,J)=FNEW(I,J)+(TIME_INCR*LIJPH(I,K,IPH)+XID3(I,K))
     #                         *FIJPH(K,J,IPH)
          ENDDO
        ENDIF
      ENDDO
      ENDDO
        
      FIJPH(:,:,IPH)=FNEW(:,:)

      iskip=1
      if(iskip.eq.0) then
       write(10,'(''  INSIDE UPD_FIJ - PHASE #'',I3)') IPH
       if(IPH.eq.0) write
     #  (10,'(''  LIJBAR '',3F10.4)') ((LIJBARc(I,J),J=1,3),I=1,3)
       if(IPH.gt.0) write      
     #  (10,'(''  LIJPH  '',3F10.4)') ((LIJPH(I,J,IPH),J=1,3),I=1,3)
       write
     #  (10,'(''  FIJ    '',3F10.4)') ((FIJPH(I,J,IPH),J=1,3),I=1,3)
      endif

C *** UPDATES THE DEFORM GRAD IN EACH GRAIN 'FIJGR(i,j,KKK)' USING THE
C     VELOCITY GRADIENT FOR THE GRAIN 'LIJGR' CALCULATED INSIDE SUBROUTINE
C     UPDATE_ORIENTATION.

      IF (ISHAPE(IPH).GT.0) THEN
        DO KKK=NGR(IPH-1)+1,NGR(IPH)

          DO I=1,3
          DO J=1,3
            FNEW(I,J)=0.0
            DO K=1,3
              FNEW(I,J)=FNEW(I,J)+(TIME_INCR*LIJGR(I,K,KKK)+XID3(I,K))
     #                           *FIJGR(K,J,KKK)
            ENDDO
          ENDDO
          ENDDO

          FIJGR(:,:,KKK)=FNEW(:,:)

        ENDDO      ! END OF DO KKK
      ENDIF

      RETURN
      END

C *****************************************************************************
C     SUBROUTINE UPDATE_GROWTH_RATE    -->   VERSION 09/JUN/2022   
C
C     USES A CONSTANT GROWTH RATE COMBINED WITH A SIMPLE CREEP LAW (NRS=1) IN 
C     A VOCE MODEL
C
C     IOPTION=0 : READS TRANSFORMATION RATE TENSOR OF SINGLE CRYSTAL
C     IOPTION=1 : INITIALIZES GAMDOS IN IRRADIATION CREEP LAW (NRS=1)
C     IOPTION=2 : ROTATES TRANSFORMATION RATE TENSOR OF EACH GRAIN FROM
C                 CRYSTAL TO SAMPLE AXES
C ******************************************************************************

      SUBROUTINE UPDATE_GROWTH_RATE (IOPTION)

      USE VPSC8DIM

      DIMENSION aux33r(3,3),SGc(3,3,10)

      IF(IOPTION.EQ.0) THEN

        READ(UR1,'(a)') PROSA
        READ(UR1,*) GAMD0
        READ(UR1,'(a)') PROSA
        READ(UR1,*) DG_TRANSCc(:,:)

        CALL CHG_BASIS(DG_TRANSC,DG_TRANSCc,AUX55,AUX3333,2,6)

        WRITE(10,'('' GAMDO IN IRRADIATION CREEP LAW'')')
        WRITE(10,'(3F12.6)') GAMD0
        WRITE(10,'('' CARTESIAN CRYSTAL GROWTH RATE'')')
        WRITE(10,'(3F12.6)') DG_TRANSCc
        WRITE(10,'('' b-BASIS CRYSTAL GROWTH RATE'')')
        WRITE(10,'(5F12.6)') DG_TRANSC

        OPEN (71,FILE='IRRAD_AUX.OUT',STATUS='UNKNOWN')
        WRITE(71,'(''  TIME     SG1ax      SG1ho      SG1ra'',     
     #             ''       SG2ax      SG2ho       SG2ra'',
     #             ''       SG3ax      SG3ho       SG3ra'')')
	 
      ENDIF

      IF(IOPTION.EQ.1) THEN
        KGX=1
        DO IPH=IPHBOT,IPHTOP
        IPHEL=IPH-IPHBOT+1
        DO KKK=NGR(IPH-1)+1,NGR(IPH)
          DO IS=1,NSYST(IPH)
            GAMD0S(IS,KGX)=GAMD0
          ENDDO
        KGX=KGX+1
        ENDDO
        ENDDO
      ENDIF

      IF(IOPTION.EQ.2) THEN

        KGX=1
        DO IPH=IPHBOT,IPHTOP
        IPHEL=IPH-IPHBOT+1
        DO KKK=NGR(IPH-1)+1,NGR(IPH)

          DO I=1,3
          DO J=1,3
            AUX33R(I,J)=0.
            DO I1=1,3
            DO J1=1,3
              AUX33R(I,J)=AUX33R(I,J)+AG(I,I1,KKK)*AG(J,J1,KKK)*
     #                                DG_TRANSCc(I1,J1)
            ENDDO
            ENDDO
          ENDDO
          ENDDO

          CALL CHG_BASIS(DG_TRANS(:,KGX),AUX33R,AUX55,AUX3333,2,5)

          KGX=KGX+1
        ENDDO
        ENDDO

C *** GET (APPROXIMATE) HYDROSTATIC COMPONENT TO CALCULATE CARTESIAN STRESS 
C     TENSOR IN A FEW GRAINS FOR PLOTTING PURPOSES
        DO KGX=1,3
          SG(6,KGX)=0.
          DO J=1,6
            SG(6,KGX)=SG(6,KGX)+BG_EL(6,J,KGX)*SBAR(J)
          ENDDO
          CALL CHG_BASIS(SG(:,KGX),SGc(:,:,KGX),AUX55,AUX3333,1,6)

C        IF(KGX.EQ.1) WRITE(71,'('' SBAR '',6F12.3)') (SBAR(J),J=1,6)
C        WRITE(71,'('' BG_EL'',6f12.3)') (BG_EL(6,J,KGX),J=1,6)
C        WRITE(71,'('' SG   '',6F12.3)') (SG(I,KGX),I=1,6)
C        WRITE(71,'('' SGcii #'',I1,3F12.3)')  KGX,(SGc(I,I,KGX),I=1,3)

        ENDDO
        WRITE(71,'(F8.3,9F11.3)')  TIME,((SGc(I,I,KGX),I=1,3),KGX=1,3)

      ENDIF

      RETURN
      END
    
C *****************************************************************************
C     SUBROUTINE UPDATE_GROWTH_ZR      -->      version 20/JUL/2017
C
C     by AP: THIS SUBROUTINE IMPLEMENTS GOLUBOV'S IRRADIATION GROWTH MODEL
C
C     Reference:
C     A.V. Barashev, S.I. Golubov, R.E. Stoller, 'Theoretical investigation of
C     microstructure evolution and deformation of zirconium under neutron 
C     irradiation', J. Nucl. Mater. 461 (2015) 85-94.
C     Refer to the paper for eq. # corresponding to various parts of the code.

C     Dimensioning of CONSTITUTIVE LAW hardening parameters valid only for 
C     single phase PX

*******************************************************************************

      SUBROUTINE UPDATE_GROWTH_ZR (IOPTION,IPH)

      USE VPSC8DIM
      USE IRR_VARS

      DIMENSION rho_reduced(4,ngrmx),rho_temp(4),EpsRate_red(4),
     & EpsRate_grain(3,3),EpsRate_global(3,3),tau_uncoupled(18),
     & grad_e(3),xk2Ge(3),xk2Gs(3),xk2Gc(3),coth(3),ctemp1(3) ! related to GB sinks
      DIMENSION rate_Nv(4),rate_Ni(4)

      iprint=0   ! controls diagnostic print-out

C ****************************************************************************
C *** Read parameters of Golubov's irradiation growth model for Zr. 

      IF (IOPTION.EQ.0) THEN    

        write(10,*) 'Start reading Golubov model parameters:'

C Hardening model parameters
C AP: The growth model in its present form does not have hardening coupled.
C     Parameters related to hardening are currently not read in.

        READ(UR1,*) rad_mean(IPH)   ! Mean grain radius
        ISYS=1
        NTWMODX=0
        DO IM=1,NMODES(IPH)
          READ(UR1,*)        ! reminder/identifier of system
          READ(UR1,*) NRSX   ! rate sensitivity parameter
          DO IS=1,NSM(IM,IPH)
            NRS(ISYS,IPH)=NRSX
            ISYS=ISYS+1
          END DO      ! end of DO IS
        END DO

        ! Read irradiation growth parameters
        READ(UR1,*) ! reminder/identifier of parameter type
        READ(UR1,*) dpa_rate
        READ(UR1,*) B_
        READ(UR1,*) rho_ref
        READ(UR1,*) f_recomb
        READ(UR1,*) f_cl
c *** specific to Zr model
        READ(UR1,*) xNmax_a
        READ(UR1,*) xNmax_c
        READ(UR1,*) dose_max_a
        READ(UR1,*) dose_crit_c
        READ(UR1,*) dose_max_c
        READ(UR1,*) A_

        READ(UR1,*) bmag_a
        READ(UR1,*) bmag_c
        READ(UR1,*) ! reminder/identifier
        READ(UR1,*) rho0(1),rho0(2),rho0(3)
        READ(UR1,*) rho0(4),rho0(5),rho0(6)
        READ(UR1,*) rho0(7:18)

!       Parameter meaning:
!       dpa_rate     ! Defect production rate NRT
!       B_           ! Irradiation creep compliance
!       rho_ref      ! Reference dislocation density for irradiation creep
!       f_recomb     ! Fraction of defects recombined in cascades
!       f_cl         ! Fraction of SIAs clustered in cascades
!       xNmax_a      ! Max. number density of sessile loops (both interst & vacancy-type) on prism planes
!       xNmax_c      ! Max. number density of sessile vacancy loops on basal planes
!       dose_max_a   ! Max. dpa dose at which nucleation of sessile loops saturates on prismatic planes
!       dose_crit_c  ! Critical dpa dose at which nucleation of sessile loops commences on basal planes
!       dose_max_c   ! Max. dpa dose at which nucleation of sessile loops saturates on basal planes
!       A_           ! Parameter for nucleation of sessile loops in eq. 25
!       bmag_a       ! Burgers vector magnitude on prismatic slip systems
!       bmag_c       ! Burgers vector magnitude on basal slip systems
!       rho0_a       ! Initial line dislocation density on prismatic planes
!       rho0_c       ! Initial line dislocation density on basal planes

        write(10,*) 'Reading Golubov model parameters complete.'

      END IF   ! END OF ioption=0

C ********************************************************************************

      IF (IOPTION.EQ.1) THEN    ! Initialize state variables
        write(10,*) 'Initializing state variables for Golubov model:'

        dpa_dose = 0.0e0
        t_count = 0.0e0
        EpsGlobal(:,:) = 0.0e0

        do IPH=IPHBOT,IPHTOP
        IPHEL=IPH-IPHBOT+1
        do igr=NGR(IPH-1)+1,NGR(IPH)
          do isl = 1,18
            rho_irr(isl,igr) = rho0(isl) ! initial line dislocation density

C AP: initialize CRSS to a dummy value to avoid NaN in subr grain_rate_and_moduli
C AP: this has to be replaced by an appropriate hardening law considering defect densities
            crss(isl,igr) = 1.0e0
          end do
          xN_i(1:4,igr) = 0.5e17  ! initial interstitial loop density
          xr_i(1:4,igr) = 1.0e-10 ! initial interstitial loop size
          xN_v(1:4,igr) = 0.5e17  ! initial vacancy loop density
          xr_v(1:4,igr) = 1.0e-10 ! initial vacancy loop size
        end do
        end do

        write(10,*) 'Initialization of state variables complete.'

      END IF   ! END OF ioption=1

C *******************************************************************************

      IF (IOPTION.EQ.2) THEN    ! Time increment
        dpa_dose = dpa_dose + dpa_rate*TIME_INCR
        t_count = t_count + TIME_INCR

        do IPH=IPHBOT,IPHTOP
        IPHEL=IPH-IPHBOT+1
        do igr=NGR(IPH-1)+1,NGR(IPH)

C Reduce dislocation densities from 18 hcp systems to 4 systems: a1, a2, a3, c
C Based on projection scheme by CNT
C Grouped according to Burgers vectors
          rho_reduced(1,igr) = rho_irr(1,igr) + rho_irr(4,igr) +
     &     0.25e0*(rho_irr(10,igr) + rho_irr(12,igr) + rho_irr(16,igr)
     &     + rho_irr(18,igr))
          rho_reduced(2,igr) = rho_irr(2,igr) + rho_irr(5,igr) +
     &     0.25e0*(rho_irr(8,igr) + rho_irr(11,igr) + rho_irr(13,igr)
     &     + rho_irr(17,igr))
          rho_reduced(3,igr) = rho_irr(3,igr) + rho_irr(6,igr) +
     &     0.25e0*(rho_irr(7,igr) + rho_irr(9,igr) + rho_irr(14,igr)
     &     + rho_irr(15,igr))
          rho_reduced(4,igr) = 0.75e0*sum(rho_irr(7:18,igr))

C These values override VPSC dislocation densities to test Golubov's model parameters
          rho_tot = sum(rho_reduced(1:4,igr))
          do isl = 1,3
            rho_tot = rho_tot + 2.0e0*pi*xN_i(isl,igr)*xr_i(isl,igr)
            rho_temp(isl) = rho_reduced(isl,igr) + 2.0e0*pi*
     &                     (xN_i(isl,igr)*xr_i(isl,igr))
          end do

          rho_tot = rho_tot + 2.0e0*pi*(xN_v(4,igr)*xr_v(4,igr)) ! eq. 1
          rho_temp(4) = rho_reduced(4,igr) + 2.0e0*pi*
     &       (xN_v(4,igr)*xr_v(4,igr))

          ! Determine GB sink strength
          ! AG: crystal to 'old' sample axes
          ! AXISPH: ellipsoid to sample axes
          ! grad_e: grain radius in ellipsoid axis
          grad_e(:) = 0.0e0
          do i = 1,3
            grad_e(i) = AXISPH(0,i,IPH)*rad_mean(IPH)
            grad_e(i) = abs(grad_e(i))
          end do

          tot_wgt = 0.0e0
          do i = 1,3
            tot_wgt = tot_wgt + 1.0e0/(grad_e(i)**2.0e0)
            ctemp1(i) = sqrt(rho_tot)*grad_e(i)
            coth(i) = 1.0e1
            if (abs(ctemp1(i)) .ge. 1.0e-10) coth(i) = 1.0e0
     &       /tanh(ctemp1(i))
            if (coth(i) .gt. 10.0e0) coth(i) = 10.0e0
            ! GB sink strength in ellipsoid axes
            xk2Ge(i) = (1.0e0/(grad_e(i)**2.0e0))*3.0e0*
     &       (ctemp1(i)**2.0e0)*(ctemp1(i)*coth(i) - 1.0e0)
     &       /((ctemp1(i)**2.0e0) - 3.0e0*(ctemp1(i)*coth(i)
     &       - 1.0e0))
          end do
          do i = 1,3
            xk2Ge(i) = (1.0e0/(grad_e(i)**2.0e0))*xk2Ge(i)/tot_wgt
          end do

          ! GB sink strength in sample axes
          ! Ellipsoid to sample axes
          ! xk2Gs(i) = AXISPH(i,j,iph)xk2Ge(j)
          xk2Gs(:) = 0.0e0
          do i = 1,3
            do j = 1,3
              xk2Gs(i) = xk2Gs(i) + AXISPH(i,j,IPH)*xk2Ge(j)
            end do
            xk2Gs(i) = abs(xk2Gs(i))
          end do

          ! GB sink strength in crystal axes
          ! Sample to crystal axes
          ! xk2Gc(i) = AG_tr(i,j,igr)xk2Gs(j) = AG(j,i,igr)xk2Gs(j)
          xk2Gc(:) = 0.0e0
          do i = 1,3
            do j = 1,3
              xk2Gc(i) = xk2Gc(i) + AG(j,i,igr)*xk2Gs(j)
            end do
            xk2Gc(i) = abs(xk2Gc(i))
          end do
          xk2GBtot = xk2Gc(1) + xk2Gc(2) + xk2Gc(3)

C *************************************************************************
C Update loop densities/ empirical Golubov model

          rate_Ni(:) = 0.0e0
          rate_Nv(:) = 0.0e0
          do isl = 1,3   ! eq. 24 -> Check the factor 3. This is inconsistent with Fig. 2
            if (xN_i(isl,igr) .le. (xNmax_a/3.0e0)) then
             rate_Ni(isl) = (xNmax_a/3.0e0/dose_max_a)*dpa_rate
             xN_i(isl,igr)= xN_i(isl,igr) + rate_Ni(isl)*TIME_INCR
            end if
!             if (xN_v(isl,igr) .le. (xNmax_a/3.0e0)) then
!              rate_Nv(isl) = (xNmax_a/3.0e0/dose_max_a)*dpa_rate
!              xN_v(isl,igr)= xN_v(isl,igr) + rate_Nv(isl)*TIME_INCR
!             end if
          end do

          if ((dpa_dose .ge. dose_crit_c) .and. (xN_v(4,igr) .lt.
     &     xNmax_c) .and. (dpa_dose .lt. dose_max_c)) then   ! eq. 25
            rate_Nv(4) = xNmax_c*(A_*exp(A_*(dpa_dose -
     &       dose_crit_c)/(dose_max_c - dose_crit_c))/(dose_max_c -
     &       dose_crit_c)/(exp(A_) - 1.0e0))*dpa_rate
            xN_v(4,igr) = xN_v(4,igr) + rate_Nv(4)*TIME_INCR
          end if

C Update loop radii
C AP: Based on personal communication with the authors, the evolution equations
C     for loop radii are modified here. An errata regarding this modification
C     to the original equations is to appear in JNM

          do isl = 1,3 ! eq. 11, 12
!             xr_v(isl,igr) = xr_v(isl,igr) + (1.0e0 - f_recomb)*r
!      &       f_cl*(1.0e0/(rho_tot + xk2GBtot) - 1.0e0/3.0e0/
!      &       rho_temp(isl))*dpa_rate*TIME_INCR/bmag_a -
!      &       xr_v(isl,igr)*rate_Nv(isl)*TIME_INCR/2./xN_v(isl,igr)
!
!             if (xr_v(isl,igr) .lt. 1.e-10) xr_v(isl,igr) = 1.e-10

            xr_i(isl,igr) = xr_i(isl,igr) - (1.0e0 - f_recomb)*
     &       f_cl*(1.0e0/(rho_tot + xk2GBtot) - 1.0e0/3.0e0/
     &       rho_temp(isl))*dpa_rate*TIME_INCR/bmag_a -
     &       xr_i(isl,igr)*rate_Ni(isl)*TIME_INCR/2./xN_i(isl,igr)

            if (xr_i(isl,igr) .lt. 1.e-10) xr_i(isl,igr) = 1.e-10
          end do

          if (dpa_dose .ge. dose_crit_c) then ! eq. 13
            xr_v(4,igr) = xr_v(4,igr) + (1.0e0 - f_recomb)*f_cl*
     &       dpa_rate*TIME_INCR/(rho_tot + xk2GBtot)/bmag_c -
     &       xr_v(4,igr)*rate_Nv(isl)*TIME_INCR/2./xN_v(4,igr)

            if (xr_v(4,igr) .lt. 1.e-10) xr_v(4,igr) = 1.e-10
          end if

C *********************************************************************
C Crystallographic strain rates

          EpsRate_red(1:4) = 0.0e0
          do isl = 1,3 ! eq. 14
            EpsRate_red(isl) = dpa_rate*(1.0e0 - f_recomb)*f_cl
     &       *(1.0e0/3.0e0 - (rho_temp(isl))/
     &       (rho_tot + xk2GBtot))
          end do
          EpsRate_red(4) = - dpa_rate*(1.0e0 - f_recomb)*f_cl*
     &     ((rho_temp(4))/(rho_tot + xk2GBtot)) ! eq. 15
          EpsRate_grain(:,:) = 0.0e0

C AP: This is the strain rate according to the paper
C     This holds only if the dislocation density is the same on all prismatic systems
!          EpsRate_grain(1,1) = EpsRate_red(1) + (EpsRate_red(2)
!      &    + EpsRate_red(3))*cos(pi/3.0e0)**2.0e0 ! eq. 16
!          EpsRate_grain(2,2) = (EpsRate_red(2) + EpsRate_red(3))
!      &    *cos(pi/6.0e0)**2.0e0 ! eq. 17
!          EpsRate_grain(3,3) = EpsRate_red(4)  ! eq. 18

! C Based on derivation by CNT, the grain level strain rate for prismatic systems with different dislocation densities
          EpsRate_grain(1,1) = 0.25e0*(EpsRate_red(1) + EpsRate_red(2))
     &     + EpsRate_red(3)
          EpsRate_grain(2,2) = 0.75e0*(EpsRate_red(2) + EpsRate_red(3))
          EpsRate_grain(3,3) = EpsRate_red(4)   ! eq. 18

C Add GB strain rate (expressed in x, y, z crystal axes)
          do i = 1,3
            EpsRate_grain(i,i) = EpsRate_grain(i,i) - dpa_rate*(1.0e0 -
     &       f_recomb)*f_cl*xk2Gc(i)/(rho_tot + xk2GBtot)
          end do

C Rotate strain rate tensor into sample axes using grain transformation matrix AG
          do i=1,3
            do j=1,3
              EpsRate_global(i,j) = 0.0e0
              do m=1,3
                do n=1,3
                  EpsRate_global(i,j) = EpsRate_global(i,j) +
     &             AG(i,m,igr)*AG(j,n,igr)*EpsRate_grain(m,n)
                end do
              end do
            end do
          end do

          do i=1,3
            do j=1,3
              EpsGlobal(i,j) = EpsGlobal(i,j) + EpsRate_global(i,j)
     &         *TIME_INCR
            end do
          end do

C Convert epsrate to b-basis --> THIS IS THE IRRAD GROWTH RATE OF THE GRAIN
          call chg_basis(aux5,EpsRate_global,aux55,aux3333,2,5)
          dg_trans(1:5,igr) = aux5(1:5)

C Empirical irradiation creep model
          do isl=1,18
            gamd0s(isl,igr) = B_*(rho_irr(isl,igr)/rho_ref)*dpa_rate
          end do

        end do      ! end of DO IGR
        end do      ! end of DO IPH

      END IF    ! END OF ioption=2

      RETURN
      END
     
C *****************************************************************************
C     SUBROUTINE UPDATE_GROWTH_ZIRC2      -->      version 20/JUL/2017
C
C     AP: THIS SUBROUTINE IMPLEMENTS AN IRRADIATION GROWTH MODEL FOR ZIRCALOY-2
C     The call is controlled by IHARDLAW=31

C     Reference: A. Patra, C.N. Tome, S.I. Golubov, "Crystal plasticity modeling 
C     of irradiation growth in Zircaloy-2", Philosophical Magazine (2017).
C     Refer to the paper for eq. # corresponding to various parts of the code
*******************************************************************************

      SUBROUTINE UPDATE_GROWTH_ZIRC2 (IOPTION,IPH)

      USE VPSC8DIM
      USE IRR_VARS

      DIMENSION rho_reduced(4,ngrmx),rho_temp(4),EpsRate_red(4),
     & EpsRate_grain(3,3),EpsRate_global(3,3),tau_uncoupled(18),
     & grad_e(3),xk2Ge(3),xk2Gs(3),xk2Gc(3),coth(3),ctemp1(3)   ! related to GB sinks
      DIMENSION rate_Ni(4),rate_Nv(4)  
      DIMENSION rate_ri(4),rate_rv(4),rho_dot(4)
      DIMENSION str_grn(5),str33(3,3),str33g(3,3)
	  
      iprint=0   ! controls diagnostic print-out

C ****************************************************************************
      IF (IOPTION.EQ.0) THEN    
        
 ! These files are used to write selected-grain information of irradiation growth
        OPEN(67,FILE='IRR_GROWTH_gr1.OUT',STATUS='UNKNOWN')
        OPEN(68,FILE='IRR_GROWTH_gr2.OUT',STATUS='UNKNOWN')
        OPEN(69,FILE='IRR_GROWTH_gr3.OUT',STATUS='UNKNOWN')
C        OPEN(70,FILE='IRR_GROWTH_gr4.OUT',STATUS='UNKNOWN')

C    Initialize parameters for Golubov's irradiation growth model
C    Read irradiation growth model parameters at the end of SX file

        write(10,*) 'Start reading growth model parameters:'

C Hardening model parameters
C     The growth model in its present form does not have coupled hardening
C     Parameters related to hardening are not read in at the moment

        READ(UR1,*) rad_mean(IPH) ! Mean grain radius

        ISYS=1
        NTWMODX=0
        DO IM=1,NMODES(IPH)
          READ(UR1,*)        ! reminder/identifier of system
          READ(UR1,*) NRSX   ! rate sensitivity parameter
          DO IS=1,NSM(IM,IPH)
            NRS(ISYS,IPH)=NRSX
            ISYS=ISYS+1
          END DO      ! end of DO IS
        END DO

        ! Read irradiation growth parameters
        READ(UR1,*) ! reminder/identifier of parameter type
        READ(UR1,*) dpa_rate
        READ(UR1,*) B_
        READ(UR1,*) rho_ref
        READ(UR1,*) f_recomb
        READ(UR1,*) f_cl

        READ(UR1,*) bmag_a
        READ(UR1,*) bmag_c
        READ(UR1,*) ! reminder/identifier
        READ(UR1,*) rho0(1:3)
        READ(UR1,*) rho0(4:6)
        READ(UR1,*) rho0(7:18)

!       Parameter meaning:
!       dpa_rate     ! Defect production rate NRT
!       B_           ! Irradiation creep compliance
!       rho_ref      ! Reference dislocation density for irradiation creep
!       f_recomb     ! Fraction of defects recombined in cascades
!       f_cl         ! Fraction of SIAs clustered in cascades
!       bmag_a       ! Burgers vector magnitude on prismatic slip systems
!       bmag_c       ! Burgers vector magnitude on basal slip systems
!       rho0_a       ! Initial line dislocation density on prismatic planes
!       rho0_c       ! Initial line dislocation density on basal planes

!       Parameter values from Table 1 in the paper:
!       f_recomb = 0.9e0
!       f_cl = 0.2e0
!       xNmax_a = 1.0e22
!       xNmax_c = 1.0e21
!       dose_max_a = 3.84e0
!       dose_crit_c = 3.0e0
!       dose_max_c = 23.0e0
!       A_ = 5.0e0
!       bmag_a = 3.0e-10
!       bmag_c = 5.0e-10

        write(10,*) 'Reading growth model parameters complete.'

      END IF   ! END OF ioption=0
C ****************************************************************************

      IF (IOPTION.EQ.1) THEN    ! Initialize state variables
        write(10,*) 'Initializing state variables for growth model:'

        dpa_dose = 0.0e0
        t_count = 0.0e0
        EpsGlobal(:,:) = 0.0e0

        do IPH=IPHBOT,IPHTOP
        IPHEL=IPH-IPHBOT+1
        do igr=NGR(IPH-1)+1,NGR(IPH)
          do isl = 1,18
            rho_irr(isl,igr) = rho0(isl)   ! initial line dislocation density

C AP: initialize CRSS to a dummy value to avoid NaN in grain_rate_and_moduli()
C AP: this has to be replaced by an appropriate hardening law considering defect densities
            crss(isl,igr) = 1.0e0
          end do
          xN_i(1:4,igr) = 5.0e18 ! initial interstitial loop density
          xr_i(1:4,igr) = 1.e-9  ! initial interstitial loop size
          xN_v(1:4,igr) = 5.0e18 ! initial vacancy loop density
          xr_v(1:4,igr) = 1.e-9  ! initial vacancy loop size
        end do
        end do

        write(10,*) 'Initialization of state variables complete.'

      END IF   ! END OF ioption=1
C ****************************************************************************

      IF (IOPTION.EQ.2) THEN    ! Time increment
        dpa_dose = dpa_dose + dpa_rate*TIME_INCR
        t_count = t_count + TIME_INCR

        do IPH=IPHBOT,IPHTOP
        IPHEL=IPH-IPHBOT+1
        do igr=NGR(IPH-1)+1,NGR(IPH)

C Reduce dislocation densities from 18 hcp systems to 4 systems: a1, a2, a3, c
C Based on projection scheme by CNT
C Grouped according to Burgers vectors
          ! Eq. (29)
          rho_reduced(1,igr) = rho_irr(1,igr) + rho_irr(4,igr) +
     &     0.25e0*(rho_irr(10,igr) + rho_irr(12,igr) + rho_irr(16,igr)
     &     + rho_irr(18,igr))
          rho_reduced(2,igr) = rho_irr(2,igr) + rho_irr(5,igr) +
     &     0.25e0*(rho_irr(8,igr) + rho_irr(11,igr) + rho_irr(13,igr)
     &     + rho_irr(17,igr))
          rho_reduced(3,igr) = rho_irr(3,igr) + rho_irr(6,igr) +
     &     0.25e0*(rho_irr(7,igr) + rho_irr(9,igr) + rho_irr(14,igr)
     &     + rho_irr(15,igr))
          rho_reduced(4,igr) = 0.75e0*sum(rho_irr(7:18,igr))

          rho_tot = sum(rho_reduced(1:4,igr))
          do isl = 1,3
            rho_tot = rho_tot + 2.0e0*pi*(xN_i(isl,igr)*xr_i(isl,igr)
     &       + xN_v(isl,igr)*xr_v(isl,igr)) ! eq. 1
            rho_temp(isl) = rho_reduced(isl,igr) + 2.0e0*pi*(
     &       xN_i(isl,igr)*xr_i(isl,igr) + xN_v(isl,igr)*xr_v(isl,igr))
          end do

          rho_tot = rho_tot + 2.0e0*pi*(xN_v(4,igr)*xr_v(4,igr)) ! Eq. (2)
          rho_temp(4) = rho_reduced(4,igr) + 2.0e0*pi*(
     &       xN_v(4,igr)*xr_v(4,igr))

          ! Determine GB sink strength
          ! grad_e: grad in ellipsoid axes
          ! AG: crystal to 'old' sample axes
          ! AXISPH: ellipsoid to sample axes

          ! grad_e: grain radius in ellipsoid axis
          grad_e(:) = 0.0e0
          do i = 1,3
            grad_e(i) = AXISPH(0,i,IPH)*rad_mean(IPH)
            grad_e(i) = abs(grad_e(i))
          end do

          tot_wgt = 0.0e0
          do i = 1,3
            tot_wgt = tot_wgt + 1.0e0/(grad_e(i)**2.0e0)
            ctemp1(i) = sqrt(rho_tot)*grad_e(i)
            coth(i) = 1.0e1
            if (abs(ctemp1(i)) .ge. 1.0e-10) coth(i) = 1.0e0
     &       /tanh(ctemp1(i))
            if (coth(i) .gt. 10.0e0) coth(i) = 10.0e0
            ! GB sink strength in ellipsoid axes
            xk2Ge(i) = (1.0e0/(grad_e(i)**2.0e0))*3.0e0*
     &       (ctemp1(i)**2.0e0)*(ctemp1(i)*coth(i) - 1.0e0)
     &       /((ctemp1(i)**2.0e0) - 3.0e0*(ctemp1(i)*coth(i)
     &       - 1.0e0))
          end do
          do i = 1,3
            xk2Ge(i) = (1.0e0/(grad_e(i)**2.0e0))*xk2Ge(i)/tot_wgt
          end do

          ! GB sink strength in sample axes
          ! Ellipsoid to sample axes
          ! xk2Gs(i) = AXISPH(i,j,iph)*xk2Ge(j)
          xk2Gs(:) = 0.0e0
          do i = 1,3
            do j = 1,3
              xk2Gs(i) = xk2Gs(i) + AXISPH(i,j,IPH)*xk2Ge(j)
            end do
            xk2Gs(i) = abs(xk2Gs(i))
          end do

          ! GB sink strength in crystal axes
          ! Sample to crystal axes
          ! xk2Gc(i) = AG_tr(i,j,igr)*xk2Gs(j) = AG(j,i,igr)*xk2Gs(j)
          xk2Gc(:) = 0.0e0
          do i = 1,3
            do j = 1,3
              xk2Gc(i) = xk2Gc(i) + AG(j,i,igr)*xk2Gs(j)
            end do
            xk2Gc(i) = abs(xk2Gc(i))
          end do
          xk2GBtot = xk2Gc(1) + xk2Gc(2) + xk2Gc(3)

C *************************************************************************
C Update loop densities and radii/ empirical model
C The nucleation model is modified for Zircaloy-2
C Fits experimental data from Holt et al., ASTM STP 1295 (1996) 623-637.

          rate_Ni(:) = 0.0e0
          rate_Nv(:) = 0.0e0
          rate_ri(:) = 0.0e0
          rate_rv(:) = 0.0e0
          rho_dot(:) = 0.0e0

          ! rho_dot from Holt et al.
          ! This is obtained from polynomial fit to the exp. data
          ! Eq. (35a)
          do isl = 1,3
            rho_dot(isl) = dpa_rate*(3.e0*0.0006627371*
     &       (dpa_dose**2.e0) - 2.e0*0.0311591329*dpa_dose +
     &      0.4757423555)*1.e14
          end do

          ! Eq. (35b)
          rho_dot(4) = dpa_rate*(3.e0*0.0006226342*(dpa_dose**2.e0)
     &     - 2.e0*0.0228959447*dpa_dose + 0.3257986076)*1.e14

          ! Eq. (A.2)
          rate_rv(4) = 2.e0*(dpa_rate*(1.e0 - f_recomb)*f_cl/
     &     bmag_c/(rho_tot + xk2GBtot) - rho_dot(4)/4.00/(22.0/7.0)
     &     /xN_v(4,igr))

          if (t_count .eq. TIME_INCR) rate_rv(4) = 0.0e0

          ! Eq. (A.1)
          rate_Nv(4) = (rho_dot(4)/2.0/(22.0/7.0) - rate_rv(4)*
     &     xN_v(4,igr))/xr_v(4,igr)

          do isl = 1,3
            ! Eq. (A.8)
            t1 = dpa_rate*(1.e0 - f_recomb)*f_cl*(1.e0/bmag_a)*
     &       (1.e0/3.e0/rho_temp(isl) - 1.e0/(rho_tot + xk2GBtot))

            ! Eq. (A.6)
            t2 = (rho_dot(isl)/2.0e0/(22.0/7.0)/xN_i(isl,igr) - t1)/
     &       (xr_i(isl,igr)/xN_i(isl,igr) + xr_v(isl,igr)/2.0e0/
     &       xN_i(isl,igr))

            rate_Ni(isl) = t2
            rate_Nv(isl) = t2

            ! Eq. (A.7)
            rate_ri(isl) = t1 - xr_v(isl,igr)*rate_Ni(isl)/2.e0/
     &       xN_i(isl,igr)

            if (t_count .le. TIME_INCR) then
              rate_ri(isl) = 0.0e0
!               rate_Ni(isl) = 0.5e0*rho_dot(isl)/2.0e0/(22.0/7.0)/
!      &         xr_i(isl,igr)
!               rate_Nv(isl) = rate_Ni(isl)
            end if
            rate_rv(isl) = 0.0e0
          end do

C Updates loop density and raddi in every grain
          do isl = 1,4
            xN_i(isl,igr) = xN_i(isl,igr) + rate_Ni(isl)*TIME_INCR
            xN_v(isl,igr) = xN_v(isl,igr) + rate_Nv(isl)*TIME_INCR

            xr_i(isl,igr) = xr_i(isl,igr) + rate_ri(isl)*TIME_INCR
            xr_v(isl,igr) = xr_v(isl,igr) + rate_rv(isl)*TIME_INCR
          end do

          do isl = 1,3
            if (xN_i(isl,igr) .lt. 5.0e18) xN_i(isl,igr) = 5.0e18
            if (xN_v(isl,igr) .lt. 5.0e18) xN_v(isl,igr) = 5.0e18

            if (xr_i(isl,igr) .lt. 1.e-9) xr_i(isl,igr) = 1.e-9
            if (xr_v(isl,igr) .lt. 1.e-9) xr_v(isl,igr) = 1.e-9
          end do
          if (xN_v(4,igr) .lt. 5.0e18) xN_v(4,igr) = 5.0e18
          if (xr_v(4,igr) .lt. 1.e-9)  xr_v(4,igr) = 1.e-9

C *********************************************************************
C Crystallographic strain rates

          EpsRate_red(1:4) = 0.0e0
          ! Eq. (14a), (15)
          do isl = 1,3
            EpsRate_red(isl) = dpa_rate*(1.0e0 - f_recomb)*f_cl
     &       *(1.0e0/3.0e0 - (rho_temp(isl))/
     &       (rho_tot + xk2GBtot))
          end do
          ! Eq. (14b), (15)
          EpsRate_red(4) = - dpa_rate*(1.0e0 - f_recomb)*f_cl*
     &     ((rho_temp(4))/(rho_tot + xk2GBtot)) 

          EpsRate_grain(:,:) = 0.0e0

! C Based on derivation by CNT, the grain level strain rate for prismatic systems with different dislocation densities
          ! Eq. (B.1a)
          EpsRate_grain(1,1) = 0.25e0*(EpsRate_red(1) + EpsRate_red(2))
     &     + EpsRate_red(3)
          ! Eq. (B.1b)
          EpsRate_grain(2,2) = 0.75e0*(EpsRate_red(2) + EpsRate_red(3))
          ! Eq. (B.1c)
          EpsRate_grain(3,3) = EpsRate_red(4)

C Add GB strain rate (expressed in x, y, z crystal axes)
          ! Eq. (17)
          do i = 1,3
            EpsRate_grain(i,i) = EpsRate_grain(i,i) - dpa_rate*(1.0e0 -
     &       f_recomb)*f_cl*xk2Gc(i)/(rho_tot + xk2GBtot)
          end do

C Rotate strain rate tensor into sample axes using grain transformation matrix AG
          do i=1,3
            do j=1,3
              EpsRate_global(i,j) = 0.0e0
              do m=1,3
                do n=1,3
                  EpsRate_global(i,j) = EpsRate_global(i,j) +
     &             AG(i,m,igr)*AG(j,n,igr)*EpsRate_grain(m,n)
                end do
              end do
            end do
          end do

          do i=1,3
            do j=1,3
              EpsGlobal(i,j) = EpsGlobal(i,j) + EpsRate_global(i,j)
     &         *TIME_INCR
            end do
          end do

C Convert epsrate to b-basis --> THIS IS THE IRRAD GROWTH RATE OF THE GRAIN
          call chg_basis(aux5,EpsRate_global,aux55,aux3333,2,5)
          dg_trans(1:5,igr) = aux5(1:5)

C Output data to UNIT=67-70 (IRR_GROWTH.OUT)
C AP: These files may become very large if the number of grains is large
C     These files may not be necessary and may be commented

          if (igr .le. 3) then
            iunitx=66+igr

          if (t_count .eq. TIME_INCR) then  ! Write column names for the 1st step
          write(iunitx,'(''    dpa_dose        time  igr'',3x,
     &     ''Exx    Eyy    Ezz'',24x,''Sxx    Syy    Szz'',24x,  
     &     ''Ni_a        Nv_c        ri_a        rv_a        rv_c'',
     &     13x,''rho_tot'')')
          end if

          str_grn(1:5) = SG(1:5,igr)
          CALL CHG_BASIS(str_grn,str33,AUX55,AUX3333,1,5)
          do i=1,3
            do j=1,3
              str33g(i,j) = 0.0e0
              do m=1,3
                do n=1,3
                  str33g(i,j) = str33g(i,j) +
     &             AG(m,i,igr)*AG(n,j,igr)*str33(m,n)
                end do
              end do
            end do
          end do

        write(iunitx,'(2E12.3,I5,3E12.3,5X,3E12.3,5X,5E12.3,5X,E12.3)')
     &     dpa_dose,t_count,igr,EpsGlobal(1,1),
     &     EpsGlobal(2,2),EpsGlobal(3,3),str33g(1,1),str33g(2,2),
     &     str33g(3,3),xN_i(1,igr),xN_v(4,igr),xr_i(1,igr),
     &     xr_v(1,igr),xr_v(4,igr),rho_tot

          end if

C Empirical irradiation creep model
          do isl=1,18          ! Eq. (21)
            gamd0s(isl,igr) = B_*(rho_irr(isl,igr)/rho_ref)*dpa_rate
          end do

        end do      ! end of DO igr
        end do      ! end of DO iph

      END IF    ! END OF ioption=2

      RETURN
      END

C *****************************************************************************
C     SUBROUTINE UPDATE_ORIENTATION    --->   VERSION 22/FEB/2022
C
C     USES SHEAR RATES AND GRAIN SHAPE TO CALCULATE GRAIN VELOCITY GRADIENT.
C     UPDATES GRAIN ORIENTATION DUE TO CRYSTALLOGRAPHIC SHEAR IN SLIP &
C     TWINNING SYSTEMS
C     TWIN REORIENTATION IS DONE INSIDE SUBR UPDATE_TWINNING, NOT HERE.
C
C     CNT: modified to use only phase-associated tensors when ISHAPE.LE.1
C     RAL: split in 3 DO loops to deal with co-rotations (17/02/00)
C     CNT: moved last loop to subroutine UPDATE_TWINNING (05/06/02)
C     CNT: fixed bug in LIJ calculation when FREEZE_SHAPE=1     (21/JUN/2008)
C     CNT: replaced subr REORIENT GRAIN with RODRIGUES   (OCT/2009)
C *****************************************************************************

      SUBROUTINE UPDATE_ORIENTATION

      USE VPSC8DIM

      dimension dnsa(3),dbsa(3),as(3,3,3,3),east(3,3)
      dimension EDOT0(3,3),WDOT0(3,3),rotloc(3,3),corot(3,3)
      DIMENSION ROT(3,3,NGRPEL),AGX(3,3),EXP_OMEGA(3,3)
      DIMENSION LIJGRX(3,3),LIJGR0(3,3)
	  REAL*8    LIJGRX, LIJGR0

      KGX=1
      DO 2000 IPH=IPHBOT,IPHTOP
        IPHEL=IPH-IPHBOT+1

        LIJPH(:,:,IPH)=0.
        IF(ISHAPE(IPH).LE.1) AS(:,:,:,:)=ASPH(:,:,:,:,IPH)

C *** THE FOLLOWING TO AVOID EMPTY 'CHILD' PHASE
        IF(WPH(IPH).LT.1.E-6) GO TO 2000

      DO 1000 KKK=NGR(IPH-1)+1,NGR(IPH)

C *** CALCULATES LOCAL ROTATION ROTLOC=PI*S**(-1)*(DG-DAV) FOR EVERY GRAIN.
C *** ROTLOC IS ZERO FOR TAYLOR CALCULATION.

        IF(INTERACTION.EQ.0) THEN

          ROTLOC(:,:)=0.

        ELSE IF(INTERACTION.GT.0) THEN     

          IF(ISHAPE(IPH).GE.2) THEN
            AS(:,:,:,:)=ASGR(:,:,:,:,KKK)
          ENDIF

          DO I=1,5
            AUX5(I)=DG(I,KGX)-DAV(I)
          ENDDO
          CALL CHG_BASIS(AUX5,EAST,AUX55,AUX3333,1,5)      

          DO I=1,3
          DO J=1,3
          ROTLOC(I,J)=0.
            DO K=1,3
            DO L=1,3
              ROTLOC(I,J)=ROTLOC(I,J)+AS(I,J,K,L)*EAST(K,L)
            ENDDO
            ENDDO
          ENDDO
          ENDDO

        ENDIF

C *** CALCULATES VELOCITY GRADIENT IN SAMPLE AXES IN EACH PHASE AND EACH GRAIN

        AGX(:,:)=AG(:,:,KKK)
        LIJGR0(:,:)=0.

        DO IS=1,NSYST(IPHEL)
          DNSA(:)=0.
          DBSA(:)=0.
          DO I=1,3
          DO J=1,3
            DNSA(I)=DNSA(I)+AGX(I,J)*DNCA(J,IS,IPHEL)
            DBSA(I)=DBSA(I)+AGX(I,J)*DBCA(J,IS,IPHEL)
          ENDDO
          ENDDO
          DO I=1,3
          DO J=1,3
            LIJGR0(I,J)=LIJGR0(I,J)+DBSA(I)*DNSA(J)*GAMDOT(IS,KGX)
          ENDDO
          ENDDO
        ENDDO

        DO I=1,3
        DO J=1,3
          EDOT0(I,J)=(LIJGR0(I,J)+LIJGR0(J,I))/2.
          WDOT0(I,J)=(LIJGR0(I,J)-LIJGR0(J,I))/2.
        ENDDO
        ENDDO

        DO I=1,3
        DO J=1,3
          LIJGRX(I,J)=WBARc(I,J)+ROTLOC(I,J)+EDOT0(I,J)
          LIJPH(I,J,IPH)=LIJPH(I,J,IPH)+LIJGRX(I,J)*WGT(KKK)/WPH(IPH)
        ENDDO
        ENDDO

        IF(ISHAPE(IPH).GT.0) LIJGR(:,:,KKK)=LIJGRX(:,:)

C *** CRYSTALLOGRAPHIC GRAIN INCREMENTAL ROTATION (RIGID minus PLASTIC)
        DO I=1,3
        DO J=1,3
          ROT(I,J,KGX)=(WBARc(I,J)+ROTLOC(I,J)-WDOT0(I,J)) *TIME_INCR
        ENDDO
        ENDDO

        KGX=KGX+1
1000  CONTINUE      ! END OF DO 1000 LOOP OVER GRAINS
2000  CONTINUE      ! END OF DO 2000 LOOP OVER PHASES

      KGX=1
      DO IPH=IPHBOT,IPHTOP
        IPHEL=IPH-IPHBOT+1
        DO KKK=NGR(IPH-1)+1,NGR(IPH)

        DO I=1,3
        DO J=1,3
          COROT(I,J)=0.
          DO IN=0,NNEIGH
            COROT(I,J)=COROT(I,J)+ROT(I,J,NEIGH(IN,KGX))*WNEIGH(IN,KGX)
          ENDDO
        ENDDO
        ENDDO

C *** CALCULATE THE NEW TRASFORMATION MATRIX AND UPDATE

        CALL RODRIGUES (COROT,EXP_OMEGA)

        AGX(:,:)=AG(:,:,KKK)
        DO I=1,3
        DO J=1,3
          AG(I,J,KKK)=0.
          DO K=1,3
            AG(I,J,KKK)=AG(I,J,KKK)+EXP_OMEGA(I,K)*AGX(K,J)
          ENDDO
        ENDDO
        ENDDO

        KGX=KGX+1
        ENDDO      ! END OF DO LOOP OVER GRAINS
      ENDDO      ! END OF DO LOOP OVER PHASES


      RETURN
      END
C
C **************************************************************************
C     SUBROUTINE UPDATE_SCHMID
C
C     ROTATES SCHMID TENSORS OF EACH GRAIN FROM CRYSTAL TO SAMPLE AXES
C **************************************************************************

      SUBROUTINE UPDATE_SCHMID

      USE VPSC8DIM

      DIMENSION aux33r(3,3)

      KGX=1
      DO IPH=IPHBOT,IPHTOP
        IPHEL=IPH-IPHBOT+1
      DO KKK=NGR(IPH-1)+1,NGR(IPH)
        DO IS=1,NSYST(IPHEL)

          DO J=1,5
            AUX5(J)=SCHCA(J,IS,IPHEL)
          ENDDO

          CALL CHG_BASIS(AUX5,AUX33,AUX55,AUX3333,1,5)

          DO I=1,3
          DO J=1,3
            AUX33R(I,J)=0.
            DO I1=1,3
            DO J1=1,3
              AUX33R(I,J)=AUX33R(I,J)+AG(I,I1,KKK)*AG(J,J1,KKK)*
     #                                AUX33(I1,J1)
            ENDDO
            ENDDO
          ENDDO
          ENDDO

          CALL CHG_BASIS(AUX5,AUX33R,AUX55,AUX3333,2,5)

          DO J=1,5
            SCH(J,IS,KGX)=AUX5(J)
          ENDDO
        ENDDO

      KGX=KGX+1
      ENDDO
      ENDDO

      RETURN
      END

C
C **************************************************************************
C     SUBROUTINE UPDATE_SHAPE      --->      VERSION 20/APR/2021
C
C     CALCULATES THE DIRECTION AND LENGTH (EIGENVECTORS AND EIGENVALUES OF
C     DEFORMATION TENSOR FIJ) OF THE AXES OF THE ELLIPSOID ASSOCIATED WITH
C     AVERAGE, PHASE, AND GRAIN ACCUMULATED DEFORMATION.
C     CALCULATES THE EULER ANGLES OF ELLIPSOID AXES WRT SAMPLE AXES.

C *** IF IPH=0 ELLIPSOID REPRESENTS AVERAGE DEFORMATION IN ELEMENT
C *** IF IPH>0 ELLIPSOID REPRESENTS AVERAGE DEFORMATION IN PHASE 'IPH'

C *** AXISPH(0,1:3,IPH) IS THE LENGTH OF THE 3 ELLIPSOID AXES ASSOCIATED
C     WITH PHASE 'IPH' (IPH=0 IS USED FOR THE MACROSCOPIC ELLPSOID (FULL PX).
C *** AXISPH(I,J,IPH) & B(I,J) ARE TRANSFORMATION MATRIX FROM ELLIPSOID TO
C     SAMPLE AXES (COLUMNS ARE THE COMPONENTS OF THE ELLIPSOID AXES IN SA).
C *** EULERPH(I,IPH) ARE EULER ANGLES OF ELLIPSOID AXES REFERRED TO SA.

C --> REPLACES PREVIOUS SUBR GRAXES_AVERAGE & SUBR GRAXES_LOCAL (SEP/04)
C --> GRAIN-FRAGMENTATION SCHEME ADDED ON 12/NOV/07 (CNT)
C --> SHAPE 'FREEZE' SCHEME MODIFIED ON 21/JUN/08 (CNT): NOW IT APPLIES TO 
C     THE DIAGONAL ELEMENTS OF FIJ (AXES' LENGTHS ARE 'FROZEN'). 
C --> WORST-GRAIN ASPECT RATIO CRITERION TO DECIDE FRAGMENTATION REPLACED BY
C     ASPECT RATIO OF THE PHASE.
C *****************************************************************************

      SUBROUTINE UPDATE_SHAPE (IPH)

      USE VPSC8DIM

      DIMENSION W(3),BX(3,3),B(3,3),BT(3,3),FIJEA(3,3)

C *** uses LIJBARc to prevent numerical instability in SO procedure
      IPHX=IPH
      IF(INTERACTION.EQ.5) IPHX=0

C ********************************************************************
C *** UPDATES ELLIPSOIDS (EIGENVALUES & EIGENVECTORS) OF PHASE AXES.
C *** 'AXISPH(I,J,IPH)' TRANSFORMS FROM ELLIPSOID TO SAMPLE AXES.
C ********************************************************************

      iskip=1
      if(iskip.eq.0) then
      write(10,'(''  DEFORMATION GRADIENT OF PHASE #'',I3)') IPH
      write(10,'(''  FIJPH  '',3F10.4)') ((FIJPH(I,J,IPH),J=1,3),I=1,3)
      endif

      DO I=1,3
      DO J=1,3
        BX(I,J)=0.
        DO K=1,3
          BX(I,J)=BX(I,J)+FIJPH(I,K,IPHX)*FIJPH(J,K,IPHX)
        ENDDO
      ENDDO
      ENDDO

      CALL EIGEN_VAL (BX,3,3,W,B,NROT,IER)
      CALL EIGEN_SORT(W,B,3,3)
      IF (IER.EQ.1) THEN
        WRITE(*,*) 'ERROR IN UPDATE_SHAPE FOR PHASE ELLIPSOID',IPH
        STOP
      ENDIF

C *** EIGENVALUES IN SUBROUTINE eigen_sort (AND ASSOC EIGENVECTORS) ARE ORDERED 
C     FROM LARGER TO SMALLER.
C *** REDEFINES AXIS(2) TO BE THE LARGEST IN ORDER TO IMPROVE ACCURACY IN THE
C     CALCULATION OF THE ESHELBY TENSOR AND RECALCULATE ELLIPSOID EULER ANGLES. 
C *** IF DET(B)<0 MEANS THAT THE SYSTEM IS LEFT HANDED. IT IS MADE RIGHT
C     HANDED BY EXCHANGING AXES 1 AND 2.
C *** B(I,J) TRANSFORMS FROM ELLIPSOID AXES TO SAMPLE AXES

      SIGN=-1.
      IF(DET(B).LE.0.) SIGN=1.
      DO I=1,3
        EXCHANGE=B(I,1)
        B(I,1)=B(I,2)
        B(I,2)=EXCHANGE*SIGN
      ENDDO
      EXCHANGE=W(1)
      W(1)=W(2)
      W(2)=EXCHANGE

      RATMAX=SQRT(W(2)/W(3))
      RATMIN=SQRT(W(1)/W(3))

      DO I=1,3
        IF(FREEZE_SHAPE(IPH).EQ.0) AXISPH(0,I,IPH)=SQRT(W(I))
        IF(FREEZE_SHAPE(IPH).EQ.1) AXISPH(0,I,IPH)=AXISPH(0,I,IPH)   ! freeze axes length
        DO J=1,3
          AXISPH(J,I,IPH)=B(J,I)
          BT(I,J)        =B(J,I)
        ENDDO
      ENDDO

      CALL EULER(1,ANG1,ANG2,ANG3,BT)
      EULERPH(1,IPH)=ANG1
      EULERPH(2,IPH)=ANG2
      EULERPH(3,IPH)=ANG3

C *** EXPRESSES AVERAGE DEFORM GRADIENT OF PHASE 'IPH' IN SAMPLE AXES.
C *** THE MAIN ELLIPSOID AXES 'AXISPH(0,J,IPH)' ARE INVARIANT IF 
C     CRITICAL ASPECT RATIO IS REACHED.

        DO I=1,3
          DO J=1,3
            FIJEA(I,J)=XID3(I,J)*AXISPH(0,I,IPH)
          ENDDO
        ENDDO
        DO I=1,3
        DO J=1,3
          FIJPH(I,J,IPH)=0.
          DO K=1,3
          DO L=1,3
            FIJPH(I,J,IPH)=FIJPH(I,J,IPH)+B(I,K)*B(J,L)*FIJEA(K,L)
          ENDDO
          ENDDO
        ENDDO
        ENDDO

C *** STOPS UPDATING THE ELLIPSOID SHAPE FOR A LIMIT ASPECT RATIO 

      IF (IPH.NE.0 .AND. IFRAG(IPH).EQ.0) THEN
        IF(RATMAX.GE.CRIT_SHP(IPH) .AND. FREEZE_SHAPE(IPH).EQ.0) THEN

          write(10,*) '*********************************************'
          write(10,*) '  STOPS UPDATING GRAIN SHAPE FOR PHASE ',iph
          write(10,*) '*********************************************'
          write(*, *) '*********************************************'
          write(*, *) '  STOPS UPDATING GRAIN SHAPE FOR PHASE ',iph
          write(*, *) '*********************************************'

          FREEZE_SHAPE(IPH)=1
        ENDIF
      ENDIF

C *************************************************************************

C *** SPLITS AVERAGE GRAIN WHEN IT EXCEEDS A CRITICAL ASPECT RATIO.
C *** REDEFINES LENGTH OF ELLIPSOID AXES 'AXISPH(0,I,IPH)' BUT NOT THE
C     ELLIPSOID ORIENTATION 'AXIS(J,I,IPH)'
C *** RECALCULATES DEFORMATION GRADIENT OF THE PHASE 'FIJPH' IN SAMPLE AXES

      IF (IPH.NE.0 .AND. IFRAG(IPH).EQ.1) THEN
        IF(RATMAX.GE.CRIT_SHP(IPH)) THEN

          W(2)=W(2)/4.
          IF(RATMIN.GE.CRIT_SHP(IPH)/2.) W(1)=W(1)/4.

ccc       w(1)=1.    ! resetting shape to equiaxed
ccc       w(2)=1.
ccc       w(3)=1.

          DO I=1,3
            AXISPH(0,I,IPH)=SQRT(W(I))
            DO J=1,3
              FIJEA(I,J)=XID3(I,J)*AXISPH(0,I,IPH)
            ENDDO
          ENDDO
          DO I=1,3
          DO J=1,3
            FIJPH(I,J,IPH)=0.
            DO K=1,3
            DO L=1,3
              FIJPH(I,J,IPH)=FIJPH(I,J,IPH)+B(I,K)*B(J,L)*FIJEA(K,L)
            ENDDO
            ENDDO
          ENDDO
          ENDDO

        ENDIF      ! END OF IF RATMAX>CRIT_RAT
      ENDIF      ! END OF IF IPH.NE.0 .AND. IFRAG.EQ.1

      iskip=1
      if(iskip.eq.0) then
      write(10,'(''  INSIDE UPD_SHAPE - ELLIPSOID OF PHASE #'',I3)') IPH
      write(10,'(''  EULER ANGLES'',3F9.2)') EULERPH(:,IPH)
      write(10,'(''  EIGNVAL'',3(F9.4,18x))') AXISPH(0,:,IPH)
      write(10,'(''  EIGNVEC'',9F9.4)') ((AXISPH(J,I,IPH),J=1,3),I=1,3)
      endif

C ************************************************************************
C *** UPDATES ELLIPSOIDS (EIGENVALUES & EIGENVECTORS) OF INDIVIDUAL GRAINS
C *** 'AXISGR(I,J)' TRANSFORMS FROM ELLIPSOID TO SAMPLE AXES.

      IF(ISHAPE(IPH).GT.0) THEN
C ************************************************************************

      DO KKK=NGR(IPH-1)+1,NGR(IPH)

        DO I=1,3
        DO J=1,3
          BX(I,J)=0.
          DO K=1,3
            BX(I,J)=BX(I,J)+FIJGR(I,K,KKK)*FIJGR(J,K,KKK)
          ENDDO
        ENDDO
        ENDDO

        CALL EIGEN_VAL (BX,3,3,W,B,NROT,IER)
        CALL EIGEN_SORT(W,B,3,3)
        IF (IER.EQ.1) THEN
          WRITE(*,'(''ERROR IN UPDATE_SHAPE FOR GRAIN'',I5,
     #            ''  IN PHASE'',I3)') KKK,IPH
          STOP
        ENDIF

C *** EIGENVALUES (AND ASSOC EIGENVECTORS) ARE ORDERED FROM LARGER TO SMALLER.
C *** REDEFINE AXIS(2) TO BE THE LARGEST IN ORDER TO IMPROVE ACCURACY IN THE
C     CALCULATION OF THE ESHELBY TENSOR.
C *** IF DET(B)<0 MEANS THAT THE SYSTEM IS LEFT HANDED. IT IS MADE RIGHT
C     HANDED BY EXCHANGING 1 AND 2.

        SIGN=-1.
        IF(DET(B).LE.0.) SIGN=1.
        DO I=1,3
          EXCHANGE=B(I,1)
          B(I,1)=B(I,2)
          B(I,2)=EXCHANGE*SIGN
        ENDDO
        EXCHANGE=W(1)
        W(1)=W(2)
        W(2)=EXCHANGE

        RATMAX=SQRT(W(2)/W(3))
        RATMIN=SQRT(W(1)/W(3))

        DO I=1,3
          IF(FREEZE_SHAPE(IPH).EQ.0) AXISGR(0,I,KKK)=SQRT(W(I))
          IF(FREEZE_SHAPE(IPH).EQ.1) AXISGR(0,I,KKK)=AXISGR(0,I,KKK)   ! freeze axes length
          DO J=1,3
            AXISGR(J,I,KKK)=B(J,I)
            BT(I,J)        =B(J,I)
          ENDDO
        ENDDO

      iskip=1
      if(iskip.eq.0) then
      write(10,'(''  INSIDE UPD_SHAPE - ELLIPSOID OF GRAIN #'',I5)') KKK
      write(10,'(''  EIGNVAL'',3(F9.4,18x))') (AXISGR(0,I,KKK),I=1,3)
      write(10,'(''  EIGNVEC'',9F9.4)') ((AXISGR(J,I,KKK),J=1,3),I=1,3)
      endif

C *** updates def gradient (keeps ellipsoid axes invariant if FREEZE_SHAPE=1)
        DO I=1,3
        DO J=1,3
          FIJEA(I,J)=XID3(I,J)*AXISGR(0,I,KKK)
        ENDDO
        ENDDO

        DO I=1,3
        DO J=1,3
          FIJGR(I,J,KKK)=0.
          DO K=1,3
          DO L=1,3
            FIJGR(I,J,KKK)=FIJGR(I,J,KKK)+B(I,K)*B(J,L)*FIJEA(K,L)
          ENDDO
          ENDDO
        ENDDO
        ENDDO

C ********************************************************************
        IF (IFRAG(IPH).EQ.1) THEN

C *** SPLITS IN HALF AXES OF THE GRAIN WHEN LIMIT ASPECT RATIOS ARE REACHED

          IF(RATMAX.GE.CRIT_SHP(IPH)) THEN
            W(2)=W(2)/4.
            IF(RATMIN.GE.CRIT_SHP(IPH)/2.) W(1)=W(1)/4.

            DO I=1,3
              AXISGR(0,I,KKK)=SQRT(W(I))
              DO J=1,3
                FIJEA(I,J)=XID3(I,J)*AXISGR(0,I,KKK)
              ENDDO
            ENDDO
            DO I=1,3
            DO J=1,3
              FIJGR(I,J,KKK)=0.
              DO K=1,3
              DO L=1,3
                FIJGR(I,J,KKK)=FIJGR(I,J,KKK)+B(I,K)*B(J,L)*FIJEA(K,L)
              ENDDO
              ENDDO
            ENDDO
            ENDDO
          ENDIF     ! END OF IF RATMAX

        ENDIF     ! END OF IF IFRAG.EQ.1
C ********************************************************************
      ENDDO     ! END OF DO KKK
C ********************************************************************
      ENDIF      ! END OF IF ISHAPE(IPH).GT.0
C ********************************************************************

      RETURN
      END
	  
C
C ************************************************************************
C    SUBROUTINE UPDATE_TWINNING_PTR   --->   VERSION OF 09/AUG/2017
C
C --> SEP/2005 (CNT): NOW THE STATS ON TWINNING ARE DONE INSIDE
C     1st DO LOOP. THE TAGGING FOR REORIENTATION IS DONE INSIDE 2nd LOOP.
C --> 14/NOV/2010: OPTION FOR TYPE I or TYPE II TWIN REORIENTATION ADDED 
C --> 26/MAR/2017: A MONTE-CARLO SELECTION (BASED ON PROBABILITY OF 
C     RANDOMLY PICKING THE PTS) WAS ADDED OLD PTS SELECTION (BASED ON 
C     MAXIMUM TWIN VOL FRACTION) 
C --> 09/AUG/2017: PREDOMINANT TWIN SYSTM OR MONTE-CARLO SELECTION
C     CONTROLED BY ITWINLAW=1 OR 2, RESPECTIVELY 
C ************************************************************************
C
C    IPTSGR(KGX): INDEX OF PREDOMINANT TWIN SYSTEM IN THE GRAIN (absolute)
C    KPTMGR(KGX): INDEX OF PREDOMINANT TWIN MODE IN THE GRAIN   (shifted)
C    TWFRGR(TWM,KKK):ACCUMULATED TWIN FRACTION IN EACH MODE IN GRAIN KKK
C                    IN THE STEP (relative to grain volume).
C
C    TWFRSY(TWS,KKK):ACCUMMULATED TWIN VOLUME FRACTION IN SYSTEM TWS
C                    IN GRAIN KKK (relative to grain volume).
C    TWFRPH(TWM,IPH):ACCUMULATED TWIN FRACT IN EACH MODE OF EACH PHASE
C                    (relative to phase volume).
C    EFTWFR(TWM,IPH):ACCUM TWIN FRACTION GIVEN BY THE REORIENTED GRAINS
C                    IN EACH MODE OF EACH PHASE (relative to phase volume).
C    KTWSMX(KKK):INDEX OF THE PREDOMINANT TWIN SYSTEM (MAX TWIN VOLUME)
C                IN THE GRAIN. USED FOR REORIENTING GRAIN.
C                   =0 IF THE GRAIN IS NOT TAGGED FOR REORIENTATION.
C                   >0 IF THE GRAIN IS TAGGED FOR REORIENTATION.
C                   <0 IF THE GRAIN HAS BEEN REORIENTED AND IS NOT TO BE
C                      RETAGGED FOR SECONDARY TWINNING
C
C    PRITW(IPH) :ACCUMMULATED TWIN FRACTION IN EACH PHASE ASSOCIATED WITH
C                GRAINS THAT WERE TWIN REORIENTED ONCE.
C    SECTW(IPH) :ACCUMMULATED TWIN FRACTION IN EACH PHASE ASSOCIATED WITH
C                GRAINS THAT WERE TWIN REORIENTED MORE THAN ONCE.
C    NTWEVENTS(KKK) :NUMBER OF TIMES THAT GRAIN 'KKK' HAS BEEN REORIENTED.
C                    0: NONE, 1: PRIMARY, 2:SECONDARY, 3:TERTIARY
C ************************************************************************

      SUBROUTINE UPDATE_TWINNING_PTR (IPH)

      USE VPSC8DIM

      DIMENSION IMARK(NGRPEL)
      DIMENSION TWB(3),TWN(3),TWBX(3),TWNX(3),ACRYS(3,3),ATWIN(3,3)
      
      IF(WPH(IPH).EQ.0.) RETURN

      IPHEL=IPH-IPHBOT+1

C *** UPDATES THE VOLUME FRACTION ASSOCIATED WITH EACH TWIN SYSTEM IN EACH GRAIN
      KGX=1
      DO KKK=NGR(IPH-1)+1,NGR(IPH)

        KTS=0                 ! relative counter
        ITS=NSLSYS(IPHEL)     ! absolute counter
        TWFRGR(0,KGX)=0.

        DO ITM=1,NTWMOD(IPHEL)       ! relative
          KTMODE=ITM+NSLMOD(IPHEL)   ! absolute
          GMODE=0.
          TWFRGR(ITM,KGX)=0.
          TWSHX=TWSH(ITM,IPHEL)

          DO NSMX=1,NSM(KTMODE,IPHEL)
            KTS=KTS+1             ! shifted counter for twin systems
            ITS=ITS+1             ! absolute counter over all systems
            IF(GAMDOT(ITS,KGX).GT.0.) THEN
              GABS=ABS(GAMDOT(ITS,KGX))*TIME_INCR
              GMODE=GMODE+GABS
              TWFRSY(KTS,KKK)=TWFRSY(KTS,KKK)+GABS/TWSHX    ! twin frac in syst
            ENDIF
            TWFRGR(ITM,KGX)=TWFRGR(ITM,KGX)+TWFRSY(KTS,KKK)   ! twin fr in mode
          ENDDO

C * total twin fraction in grain relative to grain volume
          TWFRGR(0,KGX)=TWFRGR(0,KGX)+TWFRGR(ITM,KGX)
C * total twin fraction in aggregate for each twin mode relative to phase volume
          TWFRPH(ITM,IPH)=TWFRPH(ITM,IPH)+
     #                      (GMODE/TWSHX)*(WGT(KKK)/WPH(IPH))

        ENDDO      ! END OF LOOP OVER TWIN MODES IN THE PHASE

C *** PICKS EITHER OLD PTR SCHEME (ITWINLAW=1) OR NEW MONTE CARLO SCHEME (ITWINLAW=2)
C ----------------------------------------------------------------------
        IF (ITWINLAW.EQ.1) THEN   
C ----------------------------------------------------------------------
C *** FIND THE PREDOMINANT TWIN SYSTEM (PTS) IN EACH GRAIN (OLD SCHEME)

        IPTSGR(KGX)=0
        KPTMGR(KGX)=0
        TWFMAX=0.
        KTS=0
        ITS=NSLSYS(IPHEL)

        DO ITM=1,NTWMOD(IPHEL)    ! relative index
          KTMODE=ITM+NSLMOD(IPHEL)
          DO NSMX=1,NSM(KTMODE,IPHEL)
            KTS=KTS+1             ! shifted counter for twin systems
            ITS=ITS+1             ! absolute counter over all systems
            rand=1.0
CCC           rand=0.5+0.5*random2(jran)   ! adds stochasticity to selection
            IF (TWFRSY(KTS,KKK)*rand.GT.TWFMAX) THEN
              IPTSGR(KGX)=ITS
              KPTMGR(KGX)=ITM
              TWFMAX=TWFRSY(KTS,KKK)
            ENDIF
          ENDDO
        ENDDO      ! END OF LOOP OVER TWIN MODES IN THE PHASE
C --------------------------------------------------------------------
	ELSE IF (ITWINLAW.EQ.2) THEN
C --------------------------------------------------------------------
C *** MONTE CARLO SELECTION OF PTS IN EACH GRAIN

        IPTSGR(KGX)=0
        KPTMGR(KGX)=0
        KTS=0
        ITS=NSLSYS(IPHEL)

        TWFRADD=0.           ! normalized twin fraction in the grain
        RAND=random2(jran)   ! random number in [0,1] interval 

        DO ITM=1,NTWMOD(IPHEL)    ! relative index
          KTMODE=ITM+NSLMOD(IPHEL)
          DO NSMX=1,NSM(KTMODE,IPHEL)
            KTS=KTS+1             ! shifted counter for twin systems
            ITS=ITS+1             ! absolute counter over all systems
            TWFRADD=TWFRADD+TWFRSY(KTS,KKK)/TWFRGR(0,KGX)
            IF (TWFRADD.GT.RAND) THEN
              IPTSGR(KGX)=ITS
              KPTMGR(KGX)=ITM
              GO TO 10
            ENDIF 
          ENDDO
        ENDDO      ! END OF LOOP OVER TWIN MODES IN THE PHASE
   10   CONTINUE
C ----------------------------------------------------------------------
   	ENDIF
C --------------------------------------------------------------------
        KGX=KGX+1
      ENDDO      ! END OF DO KKK OVER GRAINS IN THE PHASE

c     do kkk=ngr(iph-1)+1,ngr(iph)
c       write(10,'(''shear rates in grain'',i5,/,(6f10.4) )')
c    #    kkk,(gamdot(is,kkk),is=1,nsyst(iphel))
c       write(10,'(''twin fractions in grain'',i5,2x,12f10.4)')
c    #    kkk,(twfrsy(its,kkk),its=1,ntwsys(iphel))
c     enddo

C *********************************************************************
C    THIS LOOP SCANS ALL GRAINS IN THE PHASE/ELEM IN A RANDOM WAY AND
C    TAGS THEM FOR CHECKING AGAINST THE TWIN-REORIENTATION CRITERION.
C    IMARK(1:NGPHASE) TAGS WITH '1' THE GRAINS AS THEY ARE PICKED.
C    IN EACH LOOP ITERAT A GRAIN 'KKK' IS RANDOMLY PICKED AND REORIENTED
C    ABOUT THE PTS IF THE 'REAL' TWIN FRACTION IS NOT EXCEEDED.
C *********************************************************************

      NGPHASE=NGR(IPH)-NGR(IPH-1)
      KLEFT  =NGPHASE
      DO I=1,NGPHASE
        IMARK(I)=0
      ENDDO

      DO INDEX=1,NGPHASE    ! tags one of the yet unchecked grains
        RAND=random2(JRAN)
        KPOINT=INT(RAND*KLEFT)+1
        KLEFT=KLEFT-1
        KACUM=0
        IGR=0
        DO WHILE(IGR.LT.NGPHASE .AND. KPOINT.NE.0)    ! KPOINT: index of tagged grain
          IGR=IGR+1
          IF(IMARK(IGR).EQ.0) THEN
            KACUM=KACUM+1
            IF(KACUM.EQ.KPOINT) THEN
              KKK=IGR+NGR(IPH-1)     ! absolute index of grain
              KGX=IGR                ! index of grain relative to element
              IMARK(IGR)=1
              IPTS=IPTSGR(KGX)       ! absolute index
              KPTM=KPTMGR(KGX)       ! relative index

c     nslsyx=nslsys(iphel)
c     if(igr.le.10) then
c     write(10,'(''grain/twmod/twsys/twfrsy'',
c    #          3i10,f12.5)') kkk,KPTM,ipts,twfrsy(ipts-nslsyx,kkk)
c     write(10,'(''all shear rates'',3e12.4,/,(6e12.4))')
c    #          (gamdot(isx,kgx),isx=1,27)
c     endif
c     write(10,'(''  flipping  --> '',3i10,4f10.5)') kkk,KPTM,ipts
c    #      ,eftwfr(KPTM,1),twfrph(KPTM,1),twfrthres,twfrgr(KPTM,kgx)

C     LABELS THE GRAIN IF THE ACCUMULATED TWIN VOLUME IN THE PREDOMINANT
C     TWINNING SYSTEM (PTS) HAS EXCEEDED A THRESHOLD VALUE.
C     THIS GRAIN WILL BE COMPLETELY REORIENTED BY TWINNING.
C     REORIENTATION IS STOPPED WHEN THE EFFECTIVE REORIENTED FRACTION
C     REACHES THE 'REAL' TWIN FRACTION CALCULATED FROM THE SHEARS.
C     THRES1: MIN ACCUM TWIN FRACTION IN GRAIN BEFORE TWIN REOR IS SWITCHED ON.
C     THRES1+THRES2: EVENTUAL ACCUM TWIN FR IN GR REQUIRED TO TWIN REORIENT IT.

              IF(KTWSMX(KKK).EQ.0 .AND. IPTS.NE.0) THEN

                IF (EFTWFR(KPTM,IPH).LT.TWFRPH(KPTM,IPH)) THEN
                  THRES1=TWTHRES(1,KPTM,IPHEL)
                  THRES2=TWTHRES(2,KPTM,IPHEL)
                  TWFRTHRES=THRES1+THRES2*EFTWFR(KPTM,IPH)/
     #                                (TWFRPH(KPTM,IPH)+1.E-6)
                  IF (TWFRGR(KPTM,KGX).GT.TWFRTHRES) THEN
                    EFTWFR(KPTM,IPH)=EFTWFR(KPTM,IPH)+WGT(KKK)/WPH(IPH)
                    KTWSMX(KKK)=IPTS
c       write(10,'(''twins grain'',i5,3x,''twfrph/eftwfr/twfrthres'',
c    #       3f10.4)') kkk,twfrph(KPTM,iph),eftwfr(KPTM,iph),twfrthres
                  ENDIF
                ENDIF    

C               KTWSMX(KKK)=0   ! activate to suppress twin reorientation always

C *** REORIENTS BY TWINNING GRAINS TAGGED WITH KTWSMX>0 AND UPDATES THE
C     ORIENTATION MATRIX. 
C *** THIS IS DONE AFTER HAVING CALLED UPDATE_ORIENTATION, WHICH HAS
C     REORIENTED THE GRAIN BECAUSE OF SLIP AND TWIN SHEAR

                IF (KTWSMX(KKK).GT.0) THEN

                  KTW=KTWSMX(KKK)
                  ITWTYPEX=ITWTYPE(KTW,IPHEL)
                  DO I=1,3
                    DO J=1,3
                      ACRYS(I,J)=AG(I,J,KKK)
                    ENDDO
                    TWN(I)=DNCA(I,KTW,IPHEL)
                    TWB(I)=DBCA(I,KTW,IPHEL)
                  ENDDO

                  CALL STAT_TWINNING(IPTS,KPTM,KKK,IPH)  ! twinning statistics

C *** CALCULATES CRYSTAL-TO-SAMPLE ORIENTATION MATRIX FOR THE TWINNED GRAIN
                  CALL TWIN_ORIENTATION (ITWTYPEX,TWB,TWN,ATWIN)
                  DO I=1,3
                  DO J=1,3
                    AG(I,J,KKK)=0.
                    DO K=1,3
                      AG(I,J,KKK)=AG(I,J,KKK)+ACRYS(I,K)*ATWIN(K,J)
                    ENDDO
                  ENDDO
                  ENDDO

C *** RESET ACCUMULATED TWINNED FRACTIONS IN THE GRAIN
                  DO IS=1,NTWSYS(IPHEL)
                    TWFRSY(IS,KKK)=0.
                  ENDDO

C *** ACCUMULATES TWINNED VOLUME FRACTION IN EACH PHASE.
                  IF(NTWEVENTS(KKK).EQ.0) THEN
                    PRITW(IPH)=PRITW(IPH)+WGT(KKK)/WPH(IPH)
                  ELSE IF(NTWEVENTS(KKK).GT.0) THEN
                    SECTW(IPH)=SECTW(IPH)+WGT(KKK)/WPH(IPH)
                  ENDIF

                  NTWEVENTS(KKK)=NTWEVENTS(KKK)+1
                  IF(ISECTW(KTW,IPHEL).EQ.0 .OR. NTWEVENTS(KKK).GT.1)
     #                   KTWSMX(KKK)=-KTWSMX(KKK)         ! SUPPRESS SECTW
                  IF(ISECTW(KTW,IPHEL).EQ.1 .AND.NTWEVENTS(KKK).LE.1)
     #                   KTWSMX(KKK)= 0                   ! ALLOW SECTW

C *** RESET CRSS IN SL & TW SYSTEMS TO INITIAL (OR OTHER) VALUES
C                 GTOTGR(KKK)=0
C                 DO IS=1,NSLSYS(IPHEL)
C                   CRSS(IS,KKK)=TAU(IS,0,IPHEL)
C                 ENDDO
C                 DO IS=NSLSYS(IPHEL)+1,NSYST(IPHEL)
C                   CRSS(IS,KKK)=5.0*CRSS(IS,KKK)
C                 ENDDO

C *** QUICK & DIRTY: RESET CRSS IN PREDOMINANT TWIN SYSTEM TO FAVOR ITS DE-TWINNING.
C     SHOULD BE REPLACED BY WRONSKI'S ALGORITHM

                  CRSS(KTW,KKK)=0.5*CRSS(KTW,KKK)

                ENDIF      ! END OF IF (KTWSMX(KKK).GT.0)
              ENDIF      ! END OF IF (KTWSMX(KKK).EQ.0)
            ENDIF      ! END OF IF (KACUM.EQ.KPOINT)
          ENDIF      ! END OF IF (IMARK(IGR).EQ.0)
        ENDDO      ! END OF DO WHILE (IGR.LT.NGPHASE)
      ENDDO      ! END OF DO INDEX=1,NGPHASE

      RETURN
      END
C 
C *************************************************************************
C    SUBROUTINE UPDATE_TWINNING_VFT   --->   VERSION OF 28/APR/2022
C
C --> THIS SUBROUTINE IS CALLED IN MAIN WHEN ITWINLAW=3 AND ONLY FOR IPH=1.  
C     GOES OVER PARENT PHASE, CREATES A CHILD AND TRANSFERS WEIGHT GRADUALLY
C     FROM PARENT TO CHILD (INSTEAD OF FULLY REORIENTING THE PARENT)
C *************************************************************************
C
C    IPTSGR(KKK):    INDEX OF PREDOMINANT TWIN SYSTEM IN THE GRAIN.
C    KPTMGR(KKK):    INDEX OF PREDOMINANT TWIN MODE IN THE GRAIN.
C          ITM & KTM  = ABSOLUTE & RELATIVE INDEX, TWIN MODE
C          ITS & KTS  = ABSOLUTE & RELATIVE INDEX, TWIN SYSTEM
C    TWFRGR(KTM,KKK):  ACCUMULATED TWIN FRACT IN EACH MODE IN GRAIN KKK
C                    IN THE STEP (relative to polycrystal).
C    TWFRSY(KTS,KKK):ACCUMMULATED TWIN VOLUME FRACTION IN SYSTEM KTS
C                    IN GRAIN KKK (relative to polycrystal).
C    TWFRPH(KTM,IPH):ACCUMULATED TWIN FRACT IN EACH MODE OF EACH PHASE
C                    (relative to polycrystal).
C    KTWSMX(KKK):INDEX OF THE PREDOMINANT TWIN SYSTEM (MAX TWIN VOLUME)
C                IN THE GRAIN. USED FOR REORIENTING GRAIN.
C                   =0 IF THE GRAIN IS NOT TAGGED FOR REORIENTATION.
C                   >0 IF THE GRAIN IS TAGGED FOR REORIENTATION.
C                   <0 IF THE GRAIN HAS BEEN TWIN REORIENTED AND IS NOT TO BE
C                      RETAGGED FOR SECONDARY TWINNING (if ITWIN=0,1)
C                   <0 IF A TWIN CHILD HAS BEEN CREATED AND VOLUME TRANSFER APPLIES
C
C    PRITW(IPH) :ACCUMMULATED TWIN FRACT IN EACH PHASE ASSOCIATED WITH
C                GRAINS THAT WERE TWIN REORIENTED ONCE.
C    SECTW(IPH) :ACCUMMULATED TWIN FRACT IN EACH PHASE ASSOCIATED WITH
C                GRAINS THAT WERE TWIN REORIENTED MORE THAN ONCE.
C    VAR_VF(IVAR,ITWM,IPH) :VOL FRACT FOR EACH VARIANT, IN EACH MODE, IN EACH PHASE
C *********************************************************************************

      SUBROUTINE UPDATE_TWINNING_VFT 

      USE VPSC8DIM

      DIMENSION TWB(3),TWN(3),ACRYS(3,3),ATWIN(3,3)

      ISKIP=1        ! controls temporary printout
	  
C *** THIS LOOP UPDATES TWIN FRACTION FOR EACH TWIN SYSTEM (TWFRSY), FOR EACH
C     GRAIN, FOR EACH MODE (TWFRGR), AND FOR EACH PHASE
C *** IDENTIFIES AND TAGS THE SYSTEM CONTRIBUTING MAXIMUM TWIN VOLUME FRACTION 

      DO IPH=1,NPH
        IPHEL=IPH

      DO KKK=NGR(IPH-1)+1,NGR(IPH)

        IF(IPH.EQ.1) THEN
          KKKPA=KKK
          KKKCH=KKK+NGR(1)           ! related twin-child in phase 2 
        ELSE IF (IPH.EQ.2) THEN
          KKKPA=KKK-NGR(1)           ! related parent in phase 1 
          KKKCH=KKK
        ENDIF

        IF(KTWSMX(KKK).EQ.0) THEN
          IPTSGR(KKK)=0
          KPTMGR(KKK)=0
        ENDIF

        TWFRMAX=0.
        ITM=NSLMOD(IPHEL)
        ITS=NSLSYS(IPHEL)
        KTS=0

        DO KTM=1,NTWMOD(IPHEL)
          ITM=ITM+1
          TWFRGR(KTM,KKK)= 0.
          TWSHX= TWSH(KTM,IPHEL)
          TWTHRESX= TWTHRES(1,KTM,IPHEL)/NGR(1)     ! to compare absolute fractions

          DO NSMX=1,NSM(ITM,IPHEL)
            KTS=KTS+1               ! shifted counter for twin systems
            ITS=ITS+1               ! absolute counter over all systems

            IF(GAMDOT(ITS,KKK). GT. 0.) THEN
              TWIN_INCR=GAMDOT(ITS,KKK)*TIME_INCR/TWSHX * WGT(KKK)
              TWFRSY(KTS,KKK)=TWFRSY(KTS,KKK)+TWIN_INCR       ! PX fraction / cummulative  
              TWFRGR(KTM,KKK)=TWFRGR(KTM,KKK)+TWIN_INCR       ! PX fraction / per step / per mode
      if(iskip.eq.0) then
      write(10,'(''A1*** kkk,wgt(kkk),its, twin_incr, twfrsy(kts,kkk)'',
     #          i10,f10.5,i10,2f10.5)') 
     #          kkk, wgt(kkk), its, twin_incr, twfrsy(kts,kkk)
      endif
            ENDIF

            IF(KTWSMX(KKK).EQ.0) THEN      ! PTS not tagged yet in this grain 
              RAND=1.0
c --> activate next line to add stochasticity to selection
ccc              RAND=0.5+0.5*random2(jran)   ! for adding stochasticity

              IF (TWFRSY(KTS,KKK)*rand .GT. TWTHRESX . AND.
     #            TWFRSY(KTS,KKK) .GT. TWFRMAX) THEN
                IPTSGR(KKK)=ITS     ! assigned here but reset at next call if KTWSMX=0
                KPTMGR(KKK)=KTM
                TWFRMAX=TWFRSY(KTS,KKK)
              ENDIF
            ENDIF

          ENDDO      ! end of DO NSMX loop over twin systems in mode

C * tags the parent for child creation if a PTS has been identified (IPTSGR>0)
          IF(KTWSMX(KKK).EQ.0 .AND. IPTSGR(KKK).NE.0) 
     #                            KTWSMX(KKK)=IPTSGR(KKK)
C * total twin fraction in grain relative to PX (it is cummulative)   
          TWFRGR(0,KKK)=TWFRGR(0,KKK)+TWFRGR(KTM,KKK)             
C * total twin fraction relative to PX for each twin mode and for each phase
          TWFRPH(KTM,IPH)=TWFRPH(KTM,IPH)+ TWFRGR(KTM,KKK)

        ENDDO      ! end of DO KTM over twin modes in parent or child phase

C *****************************************************************************
C *** GOES OVER PARENT PHASE (only) AND IF PARENT WAS TAGGED ABOVE (KTWSMX>0)
C     CREATES CHILD ORIENTATION & SWITCHES SIGN OF KTWSMX TO FREEZE THIS AS THE 
C     VARIANT FOR FUTURE WEIGHT TRANSFER
C *** THIS IS DONE AFTER HAVING CALLED UPDATE_ORIENTATION IN MAIN,
C     WHICH REORIENTS THE GRAIN BECAUSE OF SHEAR (INCLUDING TWIN SHEARS)

         IF (KTWSMX(KKK).GT.0 .AND. IPH.EQ.1) THEN      ! valid for detwinning ?

           IPTS=KTWSMX(KKK)
	       ITWTYPEX=ITWTYPE(IPTS,IPHEL)
           DO I=1,3
             DO J=1,3
	           ACRYS(I,J)=AG(I,J,KKK)
	         ENDDO
	       TWN(I)=DNCA(I,IPTS,IPHEL)
	       TWB(I)=DBCA(I,IPTS,IPHEL)
	       ENDDO

           CALL STAT_TWINNING (IPTS,KPTMGR(KKK),KKK,IPH)
	 
           CALL TWIN_ORIENTATION (ITWTYPEX,TWB,TWN,ATWIN)

           DO I=1,3
           DO J=1,3
             AG(I,J,KKKCH)=0.
	         DO K=1,3
	           AG(I,J,KKKCH)=AG(I,J,KKKCH)+ACRYS(I,K)*ATWIN(K,J)
	         ENDDO
	       ENDDO
	       ENDDO

C *** ASSIGNS WEIGHT AND CRSS OF SL & TW SYSTEMS WHEN CREATING CHILD 
C *** RESETS ACCUMULATED TWINNED FRACTIONS IN THE CHILD
C *** TAGS THE MIRROR TWIN SYSTEM IN THE CHILD (REQUIRED FOR TWIN GROWTH AND
C     FOR DE-TWINNING) AND KEEPS IT UNCHANGED THOROUGHOUT SIMULATION.

           WGT_FRAC=TWFRGR(0,KKK)      ! current fraction in twin mode of grain   
           IF(WGT_FRAC .GT. WGT(KKK)) WGT_FRAC=WGT(KKK)
           KTWSMX(KKKCH)=KTWSMX(KKK)
           KPTMGR(KKKCH)=KPTMGR(KKK)

           IVARX=ITWVFT(1,KKK)      ! variant index from CALL STAT_TWINNING above
           ITWMX=ITWVFT(2,KKK)      ! twin mode index from CALL STAT_TWINNING above 
           VAR_VF(IVARX,ITWMX,IPH)=VAR_VF(IVARX,ITWMX,IPH)+WGT_FRAC

           WGT(KKKCH)= WGT_FRAC
           WGT(KKK)= WGT(KKK) - WGT_FRAC
           WPH(1)=WPH(1) - WGT_FRAC
           WPH(2)=WPH(2) + WGT_FRAC

      if(iskip.eq.0) then
      write(10,'(''C *** kkk, ktwsmx(kkk),twfrgr(0,kkk),wph(1),wph(2))'',
     #          2i10,3f10.5)') 
     #          kkk, ktwsmx(kkk),twfrgr(0,kkk),wph(1),wph(2)
      write( 10,'('' --> kkk,,twfrph(1,1),wph(1),wph(2)'', i10,3f10.5)')
     #                         kkk,twfrph(1,1),wph(1),wph(2)
      endif

            GTOTGR(KKKCH)=GTOTGR(KKK)
            DO IS=1,NSLSYS(IPHEL)
              CRSS(IS,KKKCH)=  CRSS(IS,KKK)       ! assigns current CRSS to child-slips
            ENDDO
            DO IS=NSLSYS(IPHEL)+1,NSYST(IPHEL)
              CRSS(IS,KKKCH)=   CRSS(IS,KKK)      ! assigns current CRSS to child-twins
            ENDDO
C --> the following favors detwinning in KTWSMX of child
            CRSS(KTWSMX(KKKCH),KKKCH)= 0.5* CRSS(KTWSMX(KKK),KKK)   
            TWFRGR(:,KKKCH)=0.      ! resets twin fract in every mode of new child

          ENDIF      ! END OF IF KTWSMX(KKK)>0 and IPH=1

C *****************************************************************************
C *** TRANSFERS WEIGHT FROM PARENT TO CHILD & UPDATES TOTAL FRACTIONS.
C *** BECAUSE ONLY ONE VARIANT CONTRIBUTES TO CHILD GROWTH, AND WE ARE NOT
C     COUNTING THE OTHER TWIN SYSTEMS, THE TOTAL TWINNED VOLUME REPRESENTED BY 
C     CHILDREN WOULD UNDERESTIMATE THE TOTAL TWINNED VOLUME IN THE PXTAL
C *** TO COMPENSATE, ASSIGN ALL TWIN FRACTIONS IN MODE TO PTS --> not ideal !!

          IF (KTWSMX(KKK).LT.0) THEN      ! this grain was tagged in a previous step

            KTS=ABS(KTWSMX(KKK)) - NSLSYS(IPH)
            WGT_FRAC= TWFRGR(KPTMGR(KKK),KKK)       ! accounts for all twins in mode

ccc            WGT_FRAC= TWFRGR(KTWSMX(KKK),KKK)    ! accounts for predominant twin only

            IF(WGT_FRAC .GT. WGT(KKK)) WGT_FRAC=WGT(KKK)
            IF (IPH.EQ.1) THEN
              WGT(KKK)=WGT(KKK) - WGT_FRAC
              WGT(KKKCH)=WGT(KKKCH) + WGT_FRAC
              WPH(1)=WPH(1) - WGT_FRAC
              WPH(2)=WPH(2) + WGT_FRAC

              IVARX=ITWVFT(1,KKK)      ! variant index from CALL STAT_TWINNING above
              ITWMX=ITWVFT(2,KKK)      ! twin mode index from CALL STAT_TWINNING above 
              VAR_VF(IVARX,ITWMX,IPH)=VAR_VF(IVARX,ITWMX,IPH)+WGT_FRAC

            ELSE IF (IPH.EQ.2) THEN
              WGT(KKK)=WGT(KKK) - WGT_FRAC
              WGT(KKKPA)=WGT(KKKPA) + WGT_FRAC
              WPH(1)=WPH(1) + WGT_FRAC
              WPH(2)=WPH(2) - WGT_FRAC
            ENDIF
            PRITW(1)=WPH(2)

          ENDIF

C *** child wgt assigned at creation (ITWINLAW=1,2) or updated incrementally (ITWINLAW=3)
          IF(KTWSMX(KKK) .GT. 0 .AND. IPH.EQ.1) THEN
            KTWSMX(KKK) = -KTWSMX(KKK)
            KTWSMX(KKKCH)= KTWSMX(KKK)
          ENDIF

      ENDDO      ! END OF DO KKK

      ENDDO      ! END OF DO IPH

      do iph=1,nph
      write( *,'('' ---> iph,wph(iph),twfrph(ntwmod,iph)'',
     #   i4,f8.5,3f10.5)') iph,wph(iph),(twfrph(m,iph),m=1,ntwmod(iph))
      write(10,'('' ---> iph,wph(iph),twfrph(ntwmod,iph)'',
     #   i4,f8.5,3f10.5)') iph,wph(iph),(twfrph(m,iph),m=1,ntwmod(iph))
      enddo

      RETURN
      END

C *******************************************************************
C     SUBROUTINE VAR_VEL_GRAD      --->      VERSION 14/FEB/2022
C
C     IMPOSES A NON-UNIFORM DEFORMATION HISTORY TO THE AGGREGATE.
C     THE VELOCITY GRADIENT IS FULLY IMPOSED AND READ FROM 'FILEHIST'.
C *******************************************************************
      SUBROUTINE VAR_VEL_GRAD (IOPTION)

      USE VPSC8DIM

      DIMENSION W(3),BX(3,3),B(3,3),U_diag(3,3),U_PD(3,3),U_PD_inv(3,3)
	  DIMENSION FNEW(3,3),LCOROT(3,3),DCOROT(3,3)
      REAL*8    LCOROT,DCOROT

C ***********************************************************************
      ICOROT=0      ! hardwires the non-corotational option (default)
      ICOROT=1      ! hardwires the corotational option
C ***********************************************************************
      IF(IOPTION.EQ.0) THEN

C *** WRITES LOAD CONDITIONS FILE INTO 'RUN_LOG.OUT' FILE
      WRITE(10,*)
      WRITE(10,'(''*** LOAD HISTORY FOR THIS RUN'')')
      DO IDUM=1,50
        READ(UR6,END=100,FMT='(A)') PROSA
        WRITE(10,'(A)') PROSA
      ENDDO
  100 REWIND UR6

C *** ASSUMES ALL COMPONENTS OF VELOCITY GRADIENT AND A TIME INCR ARE IMPOSED.
C *** SINCE SOME STRAIN COMPONENTS MAY BE ZERO CALCULATES A VON MISES STRAIN RATE
C     AND USES THE TIME INCR TO ENFORCE A VON MISES STRAIN INCREMENT (ICTRL=0)

        READ(UR6,*) NSTEPS
        ICTRL=0
	    TEMP_INI =0.
	    TEMP_FIN =0.
        READ(UR6,*)      ! label of file columns

C *** GENERAL SETTINGS FOR FULLY IMPOSED VARIABLE VELOCITY GRADIENT

        ILBAR(:,:)=1
        SBARc(:,:)=0.
		
        IDBARv(1:6)=1
        ISBARv(1:6)=0

        STRAIN_CONTROL=1
        initializeF =.TRUE.
      ENDIF

C **************************************************************************
C *** NON-COROTATIONAL ALGORITHM *******************************************
C **************************************************************************
      IF(IOPTION.EQ.1 .and. ICOROT.EQ.0) THEN

C *** READS COMPONENTS OF VELOCITY GRADIENT AND TIME INCREMENT FOR THE STEP.
C *** WILL GET REEDINED INTO CO-ROTATIONAL SYSTEM BEFORE EXITING SUBROUTINE.
        READ (UR6,*) ISTEPX,((LIJBARc(I,J),J=1,3),I=1,3),TIME_INCR
        WRITE(*,'('' VARIABLE VELOCITY GRADIENT'')')
        WRITE(*,'(9E12.3)') LIJBARc

        DO I=1,3
        DO J=1,3
          DBARc(I,J)=(LIJBARc(I,J)+LIJBARc(J,I))/2.
        ENDDO
        ENDDO
        CALL CHG_BASIS (DBAR,DBARc,AUX55,AUX3333,2,5)
      ENDIF      

C ***********************************************************************
C *** CO-ROTATIONAL ALGORITHM (not fully tested) ************************
C ***********************************************************************
      IF(IOPTION.EQ.1 .and. ICOROT.EQ.1) THEN

C *** READS COMPONENTS OF VELOCITY GRADIENT AND TIME INCREMENT FOR THE STEP
        READ (UR6,*) ISTEPX,((LIJBARc(I,J),J=1,3),I=1,3),TIME_INCR
C *** FORCE TRACELESS VELOCITY GRADIENT
        TRACE_L=LIJBARc(1,1)+LIJBARc(2,2)+LIJBARc(3,3)
        DO I=1,3
          LIJBARc(i,i)=LIJBARc(i,i)-TRACE_L/3.
        ENDDO
        WRITE(*,'('' VARIABLE VELOCITY GRADIENT'')')
        WRITE(*,'(9E12.3)') LIJBARc
      write(10,'('' variable velocity gradient LIJBARc'')')
      write(10,'(3f10.5)') ((LIJBARc(i,j),j=1,3),i=1,3)	  

C *** ACCUMULATES REAL DEF GRADIENT AND APPLIES THE POLAR DECOMPOSITION F=R*U
      IF(initializeF) FIJ_PD(:,:)=XID3(:,:)
      initializeF =.FALSE.
      DO I=1,3
      DO J=1,3
        FNEW(I,J)=0.
        DO K=1,3
          FNEW(I,J)=FNEW(I,J)+(TIME_INCR*LIJBARc(I,K)+XID3(I,K))
     #                        *FIJ_PD(K,J)              
        ENDDO
      ENDDO
      ENDDO
      FIJ_PD(:,:)=FNEW(:,:)

C *** finds orientation and axis of stretch ellipsoid, needed for polar decomposition
      DO I=1,3
      DO J=1,3
        BX(I,J)=0.
        DO K=1,3
          BX(I,J)=BX(I,J)+FIJ_PD(I,K)*FIJ_PD(J,K)
        ENDDO
      ENDDO
      ENDDO

C *** AXES LENGHTS SQRT(W(i)) AND ASSOC EIGENVECTORS B(i,j), GIVEN BY SUBROUTINE  
C     EIGEN_SORT ARE ORDERED FROM LARGER TO SMALLER. 
C *** IF DET(B)<0 MEANS THAT THE SYSTEM IS LEFT HANDED. IT IS MADE RIGHT
C     HANDED BY REVERSING AXIS 1.
C *** B(I,J) TRANSFORMS FROM STRETCH ELLIPSOID AXES TO SAMPLE AXES AND IS 
C     THE POLAR DECOMPOSITION ROTATION 'ROT_PD'

      CALL EIGEN_VAL (BX,3,3,W,B,NROT,IER)
      CALL EIGEN_SORT(W,B,3,3)
      IF (IER.EQ.1) THEN
        WRITE(*,*) 'ERROR IN VAR_VEL_GRAD CALC OF STRETCH ELLIPSOID'
        STOP
      ENDIF
      IF(DET(B).LE.0.) THEN
        B(:,1)=-B(:,1)  
        write(10,'('' DET(B)<0 inside var_vel_grad'')')
      ENDIF
      DO I=1,3
      DO J=1,3
        ROT_PD(I,J)  =B(I,J)
        ROT_PD_t(I,J)=B(J,I)
      ENDDO
      ENDDO

C *** CALCULATES TENSOR U_PD OF POLAR DECOMPOSITION: F=R*U     ! not used
      U_diag(:,:)=0.                
      DO I=1,3                      ! keep temporarily for checking
        U_diag(I,I)=SQRT(W(I))
      DO J=1,3                      ! keep temporarily for checking
        U_PD(I,J)=0.
        DO K=1,3
          U_PD(I,J)=U_PD(I,J)+ROT_PD_t(I,K)*FIJ_PD(K,J)
        ENDDO
      ENDDO
      ENDDO

C *** rotate strain rate component of vel grad to corotational axes
C *** test that rotation verifies orthogonality R*R_tr=I
      DO I=1,3
      DO J=1,3
        LCOROT(I,J)=0.
        AUX33(I,J) =0.
		DO K=1,3
          AUX33(I,J)=AUX33(I,J)+ROT_PD_t(I,K)*ROT_PD(K,J)
		DO L=1,3
          LCOROT(I,J)=LCOROT(I,J)+ROT_PD_t(I,K)*LIJBARc(K,L)*ROT_PD(L,J)           
        ENDDO
        ENDDO
        LIJBARc(I,J)= LCOROT(I,J)    
      ENDDO
      ENDDO

      DO I=1,3     
      DO J=1,3
        DBARc(I,J)=(LIJBARc(I,J)+LIJBARc(J,I))/2.
      ENDDO
      ENDDO

      iskip=0
      if(iskip.eq.0) then
        write(10,'('' rotation B(i,j) from ellipsoid to sample'')')
        write(10,'(3f10.5)') ((B(i,j),j=1,3),i=1,3)
        write(10,'('' testing orhtogonality ROT_PD*ROT_PD_t=I'')')
        write(10,'(3f10.5)') ((AUX33(i,j),j=1,3),i=1,3)
        write(10,'('' tensor U_diagonal'')')
        write(10,'(3f10.5)') (U_diag(i,i),i=1,3)	  

        write(10,'('' co-rotational rotation matrix ROT_PD'')')
        write(10,'(3f10.5)') ((ROT_PD(i,j),j=1,3),i=1,3)
        write(10,'('' corotational stretch U_PD'')')
        write(10,'(3f10.5)') ((U_PD(i,j),j=1,3),i=1,3)	  
        write(10,'('' deformation gradient FIJ_PD'')')
        write(10,'(3f10.5)') ((FIJ_PD(i,j),j=1,3),i=1,3)	  
        write(10,'('' co-rotational strain rate DCOROT'')')
        write(10,'(3f10.5)') ((DBARc(i,j),j=1,3),i=1,3)
        write(10,'('' co-rotational velocity gradient LCOROT'')')
        write(10,'(3f10.5)') ((LIJBARc(i,j),j=1,3),i=1,3)
      endif

      CALL CHG_BASIS (DBAR,DBARc,AUX55,AUX3333,2,5)

      ENDIF      ! end of ioption=1
C **********************************************************************

      RETURN
      END

C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     SUBROUTINE VOIGT   ---->   VERSION OF 15/MAR/2022
C
C     TRANSFORMS 6X1 MATRIX T1 INTO SECOND ORDER TENSOR T2 IF IOPT=1
C     AND VICEVERSA IF IOPT=2.
C     TRANSFORMS 6X6 MATRIX C2 INTO FOURTH ORDER TENSOR C4 IF IOPT=3
C     AND VICEVERSA IF IOPT=4.
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE VOIGT(T1,T2,C2,C4,IOPT)

      DIMENSION T1(6),T2(3,3),C2(6,6),C4(3,3,3,3)
      DIMENSION IJV(6,2)

      ijv(1,1)=1
      ijv(1,2)=1
      ijv(2,1)=2
      ijv(2,2)=2
      ijv(3,1)=3
      ijv(3,2)=3
      ijv(4,1)=2
      ijv(4,2)=3
      ijv(5,1)=1
      ijv(5,2)=3
      ijv(6,1)=1
      ijv(6,2)=2

      IF(IOPT.EQ.1) THEN
      DO I=1,6
        I1=IJV(I,1)
        I2=IJV(I,2)
        T2(I1,I2)=T1(I)
        T2(I2,I1)=T1(I)
      ENDDO
      ENDIF

      IF(IOPT.EQ.2) THEN
      DO I=1,6
        I1=IJV(I,1)
        I2=IJV(I,2)
        T1(I)=T2(I1,I2)
      ENDDO
      ENDIF

      IF (IOPT.EQ.3) THEN
      DO I=1,6
        I1=IJV(I,1)
        I2=IJV(I,2)
        DO J=1,6
          J1=IJV(J,1)
          J2=IJV(J,2)
          C4(I1,I2,J1,J2)=C2(I,J)
          C4(I2,I1,J1,J2)=C2(I,J)
          C4(I1,I2,J2,J1)=C2(I,J)
          C4(I2,I1,J2,J1)=C2(I,J)
        ENDDO
      ENDDO
      ENDIF

      IF(IOPT.EQ.4) THEN
      DO I=1,6
        I1=IJV(I,1)
        I2=IJV(I,2)
        DO J=1,6
          J1=IJV(J,1)
          J2=IJV(J,2)
          C2(I,J)=C4(I1,I2,J1,J2)
        ENDDO
      ENDDO
      ENDIF

      RETURN
      END

C**********************************************************************
C     SUBROUTINE VPFC    --->      VERSION 29/DEC/2017
C
C     PERFORMS A FULL CONSTRAINT CALCULATION (TAYLOR) CORRESPONDING 
C     TO INTERACTION=0. 
C     USES THE VISCO PLASTIC CONSTITUTIVE LAW AT THE GRAIN LEVEL AND
C     IMPOSES THE MACROSCOPIC STRAIN 'DBAR' TO EVERY GRAIN.
C     CALCULATES STRESS 'SG' AND MODULI 'MGTG' FOR EVERY GRAIN.   
C**********************************************************************

      SUBROUTINE VPFC (ISTEP)

      USE VPSC8DIM


        KGX=1
        DO IPH=IPHBOT,IPHTOP
          IPHEL=IPH-IPHBOT+1
          DO KKK=NGR(IPH-1)+1,NGR(IPH)
            DO I=1,5
              SGTRY(I,KGX)=SG(I,KKK)
            ENDDO
              CALL GRAIN_STRESS (KGX,KKK,IPHEL,IPH)
              CALL GRAIN_RATE_AND_MODULI (0,KGX,KKK,IPHEL,IPH)
           KGX=KGX+1
          ENDDO
        ENDDO

C *** CALCULATE AVERAGE STRESS AND STRAIN-RATE

        DO I=1,5
          SAV(I)=0.
          DAV(I)=0.
          KGX=1
          DO IPH=IPHBOT,IPHTOP
            DO KKK=NGR(IPH-1)+1,NGR(IPH)
              SAV(I)=SAV(I)+SG(I,KKK)*WGT(KGX)
              DAV(I)=DAV(I)+DG(I,KGX)*WGT(KGX)
              KGX=KGX+1
            ENDDO
          ENDDO
          SBAR(I)=SAV(I)      ! ONLY REQUIRED TO CALCULATE DVM & SVM
        ENDDO
        CALL CHG_BASIS(SBAR,SBARc,AUX55,AUX3333,1,5)

      RETURN
      END

C*******************************************************************************
C     SUBROUTINE VPSC     --->      VERSION 18/AUG/2020
C
C     PERFORMS A SELF CONSISTENT CALCULATION WHEN INTERACTION>0 .
C     USES THE VISCO PLASTIC CONSTITUTIVE LAW AT THE GRAIN LEVEL AND
C     SOLVES FOR STRAIN RATE 'DG', STRESS 'SG' AND MODULI 'MGTG' FOR
C     EVERY GRAIN.   
C     SOLVES FOR MACROSCOPIC STRAIN RATE 'DBAR', STRESS 'SBAR' AND
C     MODULI 'MBARTG' FOR THE POLYCRYSTAL.
C*******************************************************************************

      SUBROUTINE VPSC (ISTEP)

      USE VPSC8DIM

      DIMENSION LBARc(3,3,3,3),ESAc(3,3,3,3),RSAc(3,3,3,3)         
      DIMENSION ESA(MDIM,MDIM),ESAINV(MDIM,MDIM),ESAINVc(3,3,3,3)
      DIMENSION MAST(MDIM,MDIM),MBARNEW(MDIM,MDIM)
      DIMENSION ImSINV(MDIM,MDIM),FS(MDIM,MDIM)
      DIMENSION BGINV(MDIM,MDIM),BGAV(MDIM,MDIM),BGX(MDIM,MDIM),
     #          BC1(MDIM,MDIM),BC2(MDIM,MDIM),SG0AV(6),
     #          MGTGBG(MDIM,MDIM),MGTGBGAV(MDIM,MDIM)
      DIMENSION BGTX(MDIM,MDIM),BGTAV(MDIM,MDIM),
     #          BGTINV(MDIM,MDIM),BGTD0AV(MDIM)

      DIMENSION DBAUX(6),DBAUXc(3,3)
      DIMENSION EIGB(3,3),AXB(3)

	  REAL*8    MAST,MBARNEW,LBARc,ImSINV
	  REAL*8    MGTGBG,MGTGBGAV

C *********************************************************************
C *** A BLOCK HERE FOR DEALING WITH COMPOSITE GRAIN MODEL WAS
C     REMOVED --> CHECK OLDER VERSIONS.
C *********************************************************************

C *********************************************************************
C *** THE FOLLOWING LOOP IS MEANT TO INCREASE THE RATE SENSITIVITY NRS
C        GRADUALLY INSIDE EACH STEP, IN ORDER TO ACHIEVE CONVERGENCE 
C        WHEN NRS IS LARGE (say nrs=70) AND/OR WHEN CONVERGENCE IS TRICKY.
C        -->  IT IS ACTIVATED WHEN 'IRSVAR=1'
C *********************************************************************
      if(irsvar.eq.0) then
        jxrsini=nrs(1,1)
        jxrsfin=nrs(1,1)
        jxrstep=nrs(1,1)
      endif
cwx
cwx   comment next line to perform rate sens loop at every step
cwx
      if(istep.gt.1) jxrsini=jxrsfin
cwx
      IRS=0
      DO JXRS=JXRSINI,JXRSFIN,JXRSTEP
      IRS=IRS+1

C *** QUADRATIC EXTRAPOLATION GUESS FOR ASO and ESO
      IF(INTERACTION.EQ.5 .AND. IRS.GE.4) THEN
        CALL SO_EXTRAPOL(1./JXRS,IRS,2)
      ENDIF

      if(jxrsini.ne.jxrsfin) then
        write(*,*)
        write(*,*) 'NRS ITERATION', IRS
        write(*,*)
      endif

C *****************************************************************
C *** THE FOLLOWING LOOP IS MEANT FOR RUNNING A SECOND-ORDER
C        SIMULATION WHERE THE VP RESPONSE IS A FUNCTION OF INTRAGRANULAR
C        STRESS FLUCTUATIONS. ACTIVATED WHEN 'INTERACTION=5'
C ***************************************************************** 

      ITSO=0
      ERRESO=2.*ERRSO
      ERRASO=2.*ERRSO
      KSO=1

      DO WHILE( KSO.EQ.1. AND.
     #  (ERRESO.GT.ERRSO.OR.ERRASO.GT.ERRSO) .AND. ITSO.LT.ITMAXSO )

      IF(INTERACTION.NE.5) KSO=0
      ITSO=ITSO+1
      IF(INTERACTION.EQ.5) then
       write(*,*)
       write(*,'(a,i4)') ' SECOND-ORDER ITERATION ',itso
       write(*,*)
      ENDIF

      RELSGR=2*ERRS
      RELS=2*ERRS
      RELD=2*ERRD

C *****************************************************************
C *** THE FOLLOWING LOOP ITERATES OVER THE GRAIN STRESS, CALCULATES
C        THE VISCO-PLASTIC COMPLIANCE, AND FREEZES IT FOR SOLVING THE SELF
C        CONSISTENT EQUATIONS GIVING THE AGGREGATE VP COMPLIANCE
C ***************************************************************** 
      IT2=0

      DO WHILE( (RELSGR.GT.ERRS.OR.RELD.GT.ERRD.OR.RELS.GT.ERRS)
     #       .AND. IT2.LT.ITMAXEXT)

      IT2=IT2+1

C *****************************************************************
C *** THE FOLLOWING LOOP USES THE FROZEN VP COMPLIANCE OF THE GRAINS
C        FOR SOLVING THE SELF CONSISTENT EQS GIVING THE PX VP COMPLIANCE
C ***************************************************************** 

      IT1TG=0
      RER1TG=2*ERRM

      DO WHILE(RER1TG.GT.ERRM.AND.IT1TG.LT.ITMAXINT)

        IT1TG=IT1TG+1

C *************************************************************************
C     CALCULATE ESHELBY TENSOR S=S(Mtg) AND FS=(I-S)^(-1)*S
C     S(Mtg) WILL BE THE SAME FOR EVERY GRAIN IN A GIVEN PHASE WHEN
C     ISHAPE.LE.1
C     S(Mtg) WILL BE DIFFERENT FOR EVERY GRAIN WHEN ISHAPE.GE.2 BECAUSE THE
C     SIZE & ORIENTATION OF THE ELLIPSOID AXES IS DIFFERENT FOR EVERY GRAIN,
C     AND Mtg or Msec HAVE TO BE ROTATED TO ELLIPSOID AXES

C     WHEN Mtg=n*Msec THEN S(Mtg)=S(Msec), INDEPENDENTLY OF ISHAPE VALUE.
C *************************************************************************

C *** SKIP EVERY OTHER CALCULATION OF MSTAR IF ISKIP=1
      ISKIP=0
      IF(ISKIP.EQ.0) THEN

C **************************************************************************
C *** LOOP #1 OVER PHASES AND GRAINS
C **************************************************************************

      CALL CHG_BASIS(AUX5,AUX33,LBARTG,LBARc,3,MDIM)        ! LBARTG=MBARTG^(-1)
      KGX=1
      DO IPH=IPHBOT,IPHTOP

        IF(ISHAPE(IPH).LE.1) THEN      ! AVERAGE GRAIN SHAPE
          AXB(1:3)     =AXISPH(0,1:3,IPH)
		  EIGB(1:3,1:3)=AXISPH(1:3,1:3,IPH)
        ENDIF

      DO KKK=NGR(IPH-1)+1,NGR(IPH)

        IF(ISHAPE(IPH).GE.2) THEN      ! INDIVIDUAL GRAIN SHAPE
          AXB(1:3)     =AXISGR(0,1:3  ,KKK)
		  EIGB(1:3,1:3)=AXISGR(1:3,1:3,KKK)
        ENDIF

C *** ESHELBY CALCULATION FOR EVERY PHASE OR FOR EVERY GRAIN.
C *** INSIDE 'ESHELBY_TENSOR' SAVES C4GA_FLU(I,J,M,N)=C4GA(I,J,M,N)
C     TO BE USED LATER (IF RUN REQUIRES) BY SUBROUTINE GET_THEFLU

      IF(ISHAPE(IPH).GE.2 .OR. KKK.EQ.NGR(IPH-1)+1) THEN

      IOPTION=2
      CALL ESHELBY_TENSOR(AXB,EIGB,LBARc,ESAc,RSAc,IOPTION)

      CALL CHG_BASIS(AUX5,AUX33,ESA,ESAc,4, MDIM)

      DO I=1,MDIM
      DO J=1,MDIM
        ESAINV(I,J)=ESA(I,J)
        ImSINV(I,J) =XID6(I,J)-ESA(I,J)      
      ENDDO
      ENDDO

c      write(10,*) 'call #1 of lu_inverse inside vpsc'
      CALL LU_INVERSE(ImSINV, MDIM)
      CALL LU_INVERSE(ESAINV, MDIM)

      CALL CHG_BASIS(AUX5,AUX33,ESAINV,ESAINVc,3, MDIM)

      DO I=1,MDIM
      DO J=1,MDIM
        FS(I,J)=0.
        DO K=1,MDIM
          FS(I,J)=FS(I,J)+ImSINV(I,K)*ESA(K,J)      ! (I-S)^(-1)*S
        ENDDO
      ENDDO
      ENDDO

C *** TENSORS USED IN SO SUBROUTINE GET_THEFLU. KEEP DIMENSION =5
      DO I=1,5
      DO J=1,5
        ImSINVPH(I,J,IPH)=ImSINV(I,J)           ! (I-S)^(-1)
        FSPH(I,J,IPH)=FS(I,J)
      ENDDO
      ENDDO

C *** CALCULATE EFFECTIVE INTERACTION TENSOR MAST = NEFF * inv(I-S) * S * Mtg
C *** Mtg IS Msec FOR the secant, tangent & n_eff INTERACTIONS
C *** Mtg FOR AFFINE AND 2nd ORDER IS THE REAL Mtg  
C *** NEFFGR=1 for secant, affine, 2nd order / =nrs(1,1)for tangent / =nrs(1,1)/2 for neff

      DO I=1,MDIM
      DO J=1,MDIM
        MAST(I,J)=0.
        DO K=1,MDIM
          MAST(I,J)=MAST(I,J)+FS(I,K)*MBARTG(K,J)      ! (I-S)^(-1)*S*Mtg
        ENDDO
        if(interaction.eq.4 .AND . irsvar.eq.1) then
          MAST(I,J)=jxrs*MAST(I,J)                   ! legacy option
        else
          MAST(I,J)=NEFFGR(KGX)*MAST(I,J)            ! Eq 5-32b
        endif
      ENDDO
      ENDDO

C *** COPIES TENSORS INTO PHASE ARRAYS OR INTO GRAIN ARRAYS.
C *** 'ASPH' OR 'ASGR' ARE USED IN SUBR. UPDATE_ORIENTATION TO CALCULATE
C     SPIN-RATE DEVIATIONS.

      IF(ISHAPE(IPH).LE.1) THEN

          MASTPH(1:MDIM,1:MDIM,IPH)=MAST(1:MDIM,1:MDIM)

        DO 140 I=1,3
        DO 140 J=1,3
        DO 140 K=1,3
        DO 140 L=1,3
          DUMMY=0.
          DO K1=1,3
          DO L1=1,3
            DUMMY=DUMMY+RSAc(I,J,K1,L1)*ESAINVc(K1,L1,K,L)
          ENDDO
          ENDDO
          ASPH(I,J,K,L,IPH)=DUMMY     ! Eq 5-29b:  P*S^(-1)
  140   CONTINUE

      ENDIF

      IF(ISHAPE(IPH).GE.2) THEN

          MASTGR(1:MDIM,1:MDIM,KKK)=MAST(1:MDIM,1:MDIM)

        DO 150 I=1,3
        DO 150 J=1,3
        DO 150 K=1,3
        DO 150 L=1,3
          DUMMY=0.
          DO K1=1,3
          DO L1=1,3
            DUMMY=DUMMY+RSAc(I,J,K1,L1)*ESAINVc(K1,L1,K,L)
          ENDDO
          ENDDO
          ASGR(I,J,K,L,KKK)=DUMMY        ! Eq 5-29b:  P*S^(-1)
  150   CONTINUE

      ENDIF

      ENDIF           !  END OF IF(ISHAPE(IPH.GE.2 .OR. KKK.EQ.NGR(IPH-1)+1)

      KGX=KGX+1
      ENDDO           !  END OF DO KKK OVER GRAINS   (LOOP #1)
      ENDDO           !  END OF DO IPH OVER PHASES   (LOOP #1)

      ENDIF           !  END OF ESHELBY CALCULATION BIG IF

C **************************************************************************
C *** LOOP #2 OVER PHASES AND GRAINS
C **************************************************************************

C *** SOLVES ITERATIVELY SC EQUATION FOR GENERAL CASE OF DIFFERENT SHAPES.
C *** REQUIRED FOR MULTI-PHASE WHEN GRAINS OF EACH PHASE HAVE DIFFERENT SHAPE,
C     OR FOR ISHAPE>1, WHEN EVERY GRAIN HAS DIFFERENT SHAPE, OR BOTH.
C *** EQUATION NUMBERS REFER TO Section 1.5 IN VPSC MANUAL 

        DO I=1,MDIM
          BGTD0AV(I)=0.
          SG0AV(I)=0.
          DO J=1,MDIM
            BGAV(I,J)=0.
            BGTAV(I,J)=0.
            MGTGBGAV(I,J)=0.
          ENDDO
        ENDDO

      KGX=1
      DO IPH=IPHBOT,IPHTOP
        IPHEL=IPH-IPHBOT+1

        IF(ISHAPE(IPH).LE.1) THEN
          DO I=1,MDIM
          DO J=1,MDIM
            MAST(I,J)=MASTPH(I,J,IPH)           !  M* (includes factor 'neff')
            BC2(I,J) =MAST(I,J)+MBARTG(I,J)     ! (M*+Mtg) -> aff,so ; (M*+Msec) -> sec,neff,tg
          ENDDO
          ENDDO
        ENDIF

        DO KKK=NGR(IPH-1)+1,NGR(IPH)

        IF(ISHAPE(IPH).GE.2) THEN
          DO I=1,MDIM
          DO J=1,MDIM
            MAST(I,J)=MASTGR(I,J,KKK)
            BC2(I,J)  =MAST(I,J)+MBARTG(I,J)
          ENDDO
          ENDDO
        ENDIF

        DO I=1,MDIM
        DO J=1,MDIM
          BC1(I,J)=MAST(I,J)+MGTG(I,J,KGX)        ! (M*+Mc,tg) or (M*+Mc,sec)
        ENDDO
        ENDDO

c        write(10,*) 'call #2 of lu_inverse inside vpsc'
        CALL LU_INVERSE(BC1,MDIM)

        DO I=1,MDIM
        DO J=1,MDIM
          BGX(I,J) =0.
          BGTX(I,J)=0.
          DO K=1,MDIM
            BGX(I,J) =BGX(I,J) +BC1(I,K)*BC2(K,J)
            BGTX(I,J)=BGTX(I,J)+BC2(I,K)*BC1(K,J)
          ENDDO
          BG(I,J,KGX) =BGX(I,J)        ! Eq (5-34)  ! used in GRAIN_STRESS_ALT
        ENDDO
        ENDDO

C *** SEPARATE 5-D CALCULATION FOR SECOND ORDER APPROACH ARRAYS
        DO I=1,5
        DO J=1,5
          CHIFLU(I,J,KGX)=BC1(I,J)                        ! (M* + Mc,tg)^-1
          DUMMY2=0.
          DO K=1,5
            DUMMY2=DUMMY2 +MGTG(I,K,KGX)*BC1(K,J)
          ENDDO
          BETFLU(I,J,KGX)=DUMMY2
        ENDDO
        ENDDO
		
        DO I=1,MDIM
          SG0(I,KGX)=0.
          DO J=1,MDIM
            SG0(I,KGX)=SG0(I,KGX)+BC1(I,J)*
     #        (DBAR0(J)-DG0(J,KGX)-DG_TRANS(J,KGX))    ! Eq 5-34 - CT aug/08/19
            BGTD0AV(I)=BGTD0AV(I)+BGTX(I,J)*
     #        (DG0(J,KGX)+DG_TRANS(J,KGX))*WGT(KGX)
          ENDDO
        ENDDO

        DO I=1,MDIM
        DO J=1,MDIM
          MGTGBG(I,J) =0.
          DO K=1,MDIM
            MGTGBG(I,J)=MGTGBG(I,J)+MGTG(I,K,KGX)*BG(K,J,KGX)     ! Eq 5-40a
          ENDDO
        ENDDO
        ENDDO

C *** option IBGINV=0 setting SG0AV(I)=0 eliminated by CNT on 23/7/2019
C      because convergence is always better using extended SC equations  

        DO I=1,MDIM
          SG0AV(I)=SG0AV(I)+SG0(I,KGX)*WGT(KKK)
        DO J=1,MDIM
          BGAV(I,J)  =BGAV(I,J)  +BGX(I,J)   *WGT(KGX)          ! <Bg> in Eqs 5-41
          BGTAV(I,J) =BGTAV(I,J) +BGTX(I,J)  *WGT(KGX) 
          MGTGBGAV(I,J)=MGTGBGAV(I,J)+MGTGBG(I,J) *WGT(KGX)     ! <Mg_tg+Bg> in Eqs 5-40a/5-41a    
        ENDDO
        ENDDO

        KGX=KGX+1
      ENDDO
      ENDDO

C ****************************************************************
C *** END OF LOOP #2 OVER GRAINS AND PHASES
C ****************************************************************

CCC        BGINV(I,J)=XID5(I,J)      ! legacy: was used when IBGINV=0

      BGINV(1:MDIM,1:MDIM) =BGAV(1:MDIM,1:MDIM)
      BGTINV(1:MDIM,1:MDIM)=BGTAV(1:MDIM,1:MDIM)

c      write(10,*) 'call #3 of lu_inverse inside vpsc'	  
      CALL LU_INVERSE(BGINV,MDIM)
      CALL LU_INVERSE(BGTINV,MDIM)

      DO I=1,MDIM
      DO J=1,MDIM
        MBARNEW(I,J)=0.
        DO K=1,MDIM
          MBARNEW(I,J)=MBARNEW(I,J)+MGTGBGAV(I,K)*BGINV(K,J)      ! 5-41a: Mbar=<Mg,tg*Bg>:<Bg>^-1
        ENDDO
      ENDDO
      ENDDO

C *************************************************************************
C     Mold-Mnew (tangent) COMPARISON
C --> RER1TG CONTROLS 'DO WHILE(RER1TG.GT.ERRM.AND.IT1TG.LT.ITMAXINT)'
C --> RERDZERO IS NOT USED TO CONTROL ANYTHING

        RER1TG=
     #   TMISMATCH(MBARTG(1:MDIM,1:MDIM), MBARNEW(1:MDIM,1:MDIM),MDIM)

c       WRITE(*,'(1H+,5X,''IT1TAN='',I3,'' --> TAN MOD REL ERROR'',
c    #          2E12.3)') IT1TG,RER1TG
c       WRITE(*,'(1H+,I5,3E11.3,I7,E11.3,I7,E11.3)')
c    #        IT2,RELSGR,RELS,RELD, IT1TG,RER1TG

        MBARTG(1:MDIM,1:MDIM)=MBARNEW(1:MDIM,1:MDIM)
        LBARTG(1:MDIM,1:MDIM)=MBARNEW(1:MDIM,1:MDIM)

c        write(10,*) 'call #4 of lu_inverse inside vpsc'
        CALL LU_INVERSE(LBARTG,MDIM)
  
        DO I=1,MDIM
          DBAR0(I)=0.
          DO J=1,MDIM
            DBAR0(I)=DBAR0(I)+BGTINV(I,J)*BGTD0AV(J)      
          ENDDO
		ENDDO

      ENDDO         !  END OF (DO..WHILE) FOR INNER ITERATION (IT1TG)

C ************************************************************************************
C     SOLVES (DBAR-DBAR0)=MBAR * SBAR ACCOUNTING FOR MIXED BOUNDARY CONDITIONS
C     IN STRESS AND STRAIN RATE COMPONENTS.

C     IF MBAR IS VISCO-PLASTIC INCOMPRESSIBILITY DICTATES THAT EITHER 3 OR 1 
C     DIAGONAL COMPONENTS OF D(I,J) CAN BE ENFORCED (0 OR 2 DIAGONAL STRESS COMPS) 
C     (SINCE D1~(d22-d11) & D2~(d33-d11)+(d33-d22) THEN EITHER BOTH OR NONE OF THE 
C     FIRST TWO b-BASIS COMPONENTS ARE ENFORCED, RESPECTIVELY).
C     --> MBAR IS 5x5 DEVIATORIC, DBAR IS 5X1 DEVIATORIC
C     --> IF 3 DIAGONAL COMPS OF DBAR ARE ENFORCED ONLY THE DEVIATORIC PART OF SBAR 
C         CAN BE OBTAINED. THE CONSTITUTIVE EQUATION CAN BE SOLVED IN THE b-BASIS 
C         REPRESENTATION AS A 5x5 SYSTEM 
C     --> IF ONE DIAGONAL COMPONENT OF SBAR IS ENFORCED THE OVERALL PRESSURE CAN BE 
C         CALCULATED (ALTHOUGH THE BEHAVIOR IS STILL INCOMPRESSIBLE AND INDEPENDENT
C         OF THE PRESSURE). 
C         THE FLAGS ON THE 'CARTESIAN COMPONENT' CONTROLING BOUNDARY CONDITIONS
C         CANNOT BE WRITTEN IN TERMS OF b-BASIS FLAGS AND THE SYSTEM HAS TO BE SOLVED
C         IN THE CARTESIAN REPRESENTATION WRITEN AS A 6x6 SYSTEM IN VOIGT NOTATION            

C     IF MBAR IS ELASTO-VISCO-PLASTIC: MBAR, SBAR & DBAR ARE 6X6 & 6X1 NON-DEVIAT.
C     THE CONSTITUTIVE EQUATION IS 6x6 IN b-BASIS REPRESENTATION.  
C     FBECAUSE THE CARTESIAN FLAGS ON IMPOSED BC's CANNOT BE TRANSLATED INTO b-BASIS  
C     FLAGS THE PROBLEM HAS TO BE SOLVED AS A 6x6 EQUATION IN TERMS OF VOIGT COMPS. 
C ************************************************************************************

C *** SUBSTRACT INDEPENDENT TERM, SOLVE LINEAR EQS, ADD IT AFTERWARDS 

      DBAUX(:)= DBAR(:)-DBAR0(:)    

      IF (IDBARv(1)*IDBARv(2)*IDBARv(3).EQ.1 .AND. MDIM.EQ.5) THEN 
      
        IDBARb(1)=1
        IDBARb(2)=1
        IDBARb(6)=1
        IDBARb(3)=IDBARv(4)
        IDBARb(4)=IDBARv(5)
        IDBARb(5)=IDBARv(6)		
        ISBARb(1)=0
        ISBARb(2)=0
        ISBARb(6)=0
        ISBARb(3)=ISBARv(4)
        ISBARb(4)=ISBARv(5)
        ISBARb(5)=ISBARv(6)

        CALL STATE_NxN (IDBARb(1:5),DBAUX,ISBARb(1:5),SBAR,MBARTG,5)      

        DBAR(:)=DBAUX(:)+DBAR0(:)
        CALL CHG_BASIS (DBAR,DBARc  ,AUX55,AUX3333,1,5)
        CALL CHG_BASIS (SBAR,SBARC  ,AUX55,AUX3333,1,5)

        DO KGX=1,NGR(NPH)
          SG(6,KGX)=0.
          DG(6,KGX)=0.
        ENDDO		   

      ELSE IF(IDBARv(1)*IDBARv(2)*IDBARv(3).EQ.0 .OR. MDIM.EQ.6) THEN  
          
        CALL CHG_BASIS (DBAUX,DBAUXc,AUX55,AUX3333,1, MDIM)
        CALL CHG_BASIS (AUX5,AUX33, MBARTG,AUX3333,3, MDIM)
		
          CALL STATE_6x6 (IDBARv,DBAUXc,ISBARv,SBARc,AUX3333)
		  
        CALL CHG_BASIS (SBAR, SBARc, AUX55,AUX3333,2,6)         ! SBAR is non-deviatoric for these conditions
        CALL CHG_BASIS (DBAUX,DBAUXc,AUX55,AUX3333,2,6)	
        DBAR(:)=DBAUX(:)+DBAR0(:)                               ! DBAR is deviatoric when NDIM=5
        CALL CHG_BASIS (DBAR,DBARc,AUX55,AUX3333,1, MDIM)    

C *** THE FOLLOWING IS A CONSISTENCY CHECK ON MACRO PRESSURE
      CALL CHG_BASIS (SBAR, AUX33, AUX55,AUX3333,1,5)      ! deviatoric comps of SBAR  
      PBAR1=0.
      PBAR2=0.
      DO I=1,3
	    PBAR1=PBAR1+ISBARv(I)*(SBARc(I,I)-AUX33(I,I))
		PBAR2=PBAR2+SBARc(I,I) 
      ENDDO
      PBAR1=PBAR1 / (ISBARv(1)+ISBARv(2)+ISBARv(3))
	  PBAR2=PBAR2 /3.
	  
C      write(10,'('' CONSISTENCY CHECK ON MACRO PRESSURE'')')
C      write(10,'('' SBARc ='',9f12.4)') SBARc
C      write(10,'('' PBAR1  &  PBAR2 ='',9f12.4)') PBAR1,PBAR2
C      write(10,'('' DBARc ='',9f12.4)') DBARc

C *** CALCULATES HYDROSTATIC STRESS IN GRAINS USING THE ELASTIC LOCALIZATION
C      EQUATION AND THE 'RECONSTRUCTED' MACROSCOPIC CAUCHY STRESS.
C *** ASSUMES THAT THE DEVIATORIC STRESS COMPONENTS GIVEN BY THE VISCO-
C      PLASTIC FORMALISM ARE CLOSE TO THE ONES THAT AN ELASTO-VISCO-PLASTIC
C      FORMALISM WOULD PREDICT.

         IF(MDIM.EQ.5) THEN
           DO KGX=1,NGR(NPH)      
             SG(6,KGX)=0.
             DG(6,KGX)=0.
             DO J=1,6
              SG(6,KGX)=SG(6,KGX)+BG_EL(6,J,KGX)*SBAR(J)
             ENDDO
           ENDDO	
         ENDIF		 

      ENDIF

      SVM=SQRT(3./2.)*VNORM(SBAR,5)      ! Svm is defined for the deviatoric stress tensor

C *************************************************************************
C     S* & D* REQUIRED FOR DIFFERENT GRAIN SHAPES
C     SASTAV(I) IS USED INSIDE GRAIN_STRESS
C     SASTBAR(I) IS USED BELOW TO CALCULATE DASTBAR(I)

      SASTMIX=0.
      IF (INTERACTION.NE.3.AND.INTERACTION.NE.4) SASTMIX=0.75      ! why neff & tg only? 

      DO I=1,MDIM
        SASTAV(I) =0.
        SASTBAR(I)=0.
        DO J=1,MDIM
          SASTAV(I) =SASTAV(I) +BGINV(I,J)*(SAV(J) -SG0AV(J))      ! Eq 5-34a
          SASTBAR(I)=SASTBAR(I)+BGINV(I,J)*(SBAR(J)-SG0AV(J))
        ENDDO
        SASTAVX  =SASTAV(I)
        SASTBARX=SASTBAR(I)
        SASTAV(I)  =(1.-SASTMIX)*SASTAVX+    SASTMIX*SASTBARX
        SASTBAR(I)=     SASTMIX*SASTAVX+(1.-SASTMIX)*SASTBARX
      ENDDO

      DO I=1,MDIM
        DASTBAR(I)=0.
        DO J=1,MDIM
         DASTBAR(I)=DASTBAR(I)+MBARTG(I,J)*SASTBAR(J)     ! DBAR0 added to DBAR or DASTBAR inside N-R
        ENDDO
      ENDDO

C *** SAVE CURRENT STRESS IN GRAIN AS 'SGTRY'.
C *** STARTING FROM CONVERGED 'MBARTG' AND 'SGTRY' RECALCULATES 'SG' USING
C     THE INTERACTION EQUATION.
C *** USES NEW 'SG' TO CALCULATE STRAIN-RATE 'DG', SHEAR RATES 'GAMDOT',
C     AND MODULI 'MGTG' FOR EVERY GRAIN.

      KGX=1
      DO IPH=IPHBOT,IPHTOP
        IPHEL=IPH-IPHBOT+1
      DO KKK=NGR(IPH-1)+1,NGR(IPH)

        SGTRY(1:6,KGX)=SG(1:6,KKK)

        IF(INTERACTION.NE.5) THEN          ! these routines need revision if MDIM=6
          CALL GRAIN_STRESS (KGX,KKK,IPHEL,IPH)
          CALL GRAIN_RATE_AND_MODULI (1,KGX,KKK,IPHEL,IPH)
        ENDIF

        IF(INTERACTION.EQ.5) THEN          ! this case is not applicable when MDIM=6
          if(irs.eq.1.and.itso.le.2.and.istep.eq.1) THEN
            CALL GRAIN_STRESS (KGX,KKK,IPHEL,IPH)
            CALL GRAIN_RATE_AND_MODULI (1,KGX,KKK,IPHEL,IPH)
          else
            CALL SO_GRAIN_STRESS_ALT (KGX,KKK,IPHEL,IPH)
          endif
        ENDIF

      KGX=KGX+1
      ENDDO
      ENDDO

C *** UPDATE AVERAGE STRESS AND AVERAGE STRAIN-RATE.

      SAV(:) =0.
      DAV(:) =0.
      DO I=1,6
        KGX=1
        DO IPH=IPHBOT,IPHTOP
          DO KKK=NGR(IPH-1)+1,NGR(IPH)
            SAV(I)=SAV(I)+SG(I,KKK)*WGT(KGX)
            DAV(I)=DAV(I)+DG(I,KGX)*WGT(KGX)
            KGX=KGX+1
          ENDDO
        ENDDO
      ENDDO
	  
      SDEV =0.
      DO I=1,MDIM
        KGX=1
        DO IPH=IPHBOT,IPHTOP
          DO KKK=NGR(IPH-1)+1,NGR(IPH)
            SDEV=SDEV+(SG(I,KKK)-SGTRY(I,KGX))**2*WGT(KGX)
            KGX=KGX+1
          ENDDO
        ENDDO
      ENDDO	  

C *** ROUNDS-UP TENSOR COMPONENTS 	  
      SNORM=TNORM(SBARc,3) *1.E-6
      DNORM=TNORM(DBARc,3) *1.E-6
      DO I=1,3
      DO J=1,3
        IF(ABS(SBARc(I,J)) .LT. SNORM) SBARc(I,J)=0.
        IF(ABS(DBARc(I,J)) .LT. DNORM) DBARc(I,J)=0.
      ENDDO
      ENDDO

C *** CHECK CONSISTENCY OF LOCAL, MACROSCOPIC AND AVERAGE MAGNITUDES
C      BY COMPARING RELATIVE DEVIATIONS AGAINST TOLERANCES.
C      a) <SG-SGTRY> .LT.ERR
C      b) /DBAR-DAV/   .LT.ERR
C      c) /SBAR-SAV/   .LT.ERR

      RELSGR=SQRT(SDEV)/VNORM(SBAR,5)
      RELD    =VMISMATCH(DBAR,DAV,5)
      RELS    =VMISMATCH(SBAR,SAV,5)

C *** DEFINES EMPIRICAL MIX COEF FOR SASTAV & SASTBAR BASED ON CONVERGENCE
C     CONVNEW=RELSGR+RELD+RELS
C     CONVDIF=CONVNEW-CONVOLD
C     SASTMIX=0.
C     IF(IT2.GT.10 .AND. CONVDIF.GT.0.) SASTMIX=0.25
C     IF(IT2.GT.20 .AND. CONVDIF.GT.0.) SASTMIX=0.50
C     IF(IT2.GT.30 .AND. CONVDIF.GT.0.) SASTMIX=0.75
C     CONVOLD=CONVNEW

C      WRITE(*,'(1H+,I6,3E11.3,I7,E11.3,I7,E11.3)')
C     #      IT2,RELSGR,RELS,RELD
C      WRITE(UW1,'(  I6,3E12.3,I7,E12.3,I7,E12.3)')
C     #      IT2,RELSGR,RELS,RELD

      ENDDO     ! END OF (DO..WHILE) FOR OUTER ITERATION (IT2) ON 'SG'

      WRITE(  *,'(I5,3E11.3,I7,E11.3,I7,E11.3)')
     #      IT2,RELSGR,RELS,RELD
      WRITE(UW1,'(I5,3E12.3,I7,E12.3,I7,E12.3)')
     #      IT2,RELSGR,RELS,RELD
      WRITE( 10,'(I5,3E12.3,I7,E12.3,I7,E12.3)')
     #      IT2,RELSGR,RELS,RELD
      IF(RELSGR.GT.ERRS .OR. RELD.GT.ERRD .OR. RELS.GT.ERRS) THEN
        CONVERGED=.FALSE.
        WRITE( *,'('' CONVERGENCE FAILED AT STEP'',I4,
     #   '' RELSGR,RELS,RELD ='',3E12.3)') ISTEP,RELSGR,RELS,RELD
        WRITE(10,'('' CONVERGENCE FAILED AT STEP'',I4,
     #   '' RELSGR,RELS,RELD ='',3E12.3)') ISTEP,RELSGR,RELS,RELD
      ENDIF

C	  write(10,'(''  SAV '',6F10.4)') SAV
C	  write(10,'(''  SBAR'',6F10.4)') SBAR
C     call chg_basis(aux5,aux33,MBARTG,aux3333,3,5)
C     call voigt    (aux6,aux33,AUX66 ,aux3333,4)
C	  write(10,'(''  MBARTG for step'',i5,/,(6F12.4))') istep,AUX66
C	  write(*  ,'(''  SAV '',6F10.4)') SAV
C	  write(*  ,'(''  SBAR'',6F10.4)') SBAR
	 
      SVM=0.
      DVM=0.
      DO I=1,5
        SVM=SVM+SBAR(I)*SBAR(I)      ! Von Mises is defined for deviatoric tensor
        DVM=DVM+SBAR(I)*DBAR(I)
      ENDDO
      SVM=SQRT(SVM*3./2.)
      DVM=DVM/SVM

C *** CALLS N_EFFECTIVE WHEN N_EFF(kgr) IS GIVEN BY THE RELATIVE DIRECTIONAL
C     COMPLIANCE CRITERION 
C *** THIS REDEFINITION IS DONE AFTER STEP=1 HAS CONVERGED USING THE MGSEC 
C     OBTAINED WITH N_EFF(kgr)=cte

      IF(INTERACTION.EQ.3 .AND. RDCrun) CALL N_EFFECTIVE (ISTEP)

      IF (IFLU.EQ.1) CALL SO_FLUCTUATIONS(ERRESO,ERRASO)

      if(interaction.eq.5.and.itso.ge.2) then
        write(97,'(2i5,2x,3e12.3)') irs,itso,erraso,erreso
        call so_mod
      endif

      ENDDO     ! END OF (DO..WHILE) FOR SO ITERATION (ITSO)

      IF (IFLU.EQ.1) THEN
       write(83,'(a,i7)') ' NRS = ',jxrs
       if(icubcom.eq.1) write(84,'(a,i7)') ' NRS = ',jxrs
       call so_sdpx
       if(interaction.eq.5) then
        write(97,'(a,i7,4e11.4)')
     #   ' NRS,SVM,SDS,SDD,UTILDE = '
     #   ,jxrs,svm,sdseqintra,sddeqintra,utilde
        write(97,*)
       endif
      ENDIF
	  
C *** UPDATES VALUES FOR PERFORMING EXTRAPOLATION ON ESO and ASO
      IF(INTERACTION.EQ.5) THEN
        CALL SO_EXTRAPOL(1./JXRS,IRS,1)
      ENDIF

      ENDDO     ! END OF DO JXRS FOR NRS ITERATION

      IF(INTERACTION.EQ.5) CALL SO_GET_GAMDOT

      RETURN
      END

C ********************************************************************************
C     SUBROUTINE VPSC_INPUT      --->      VERSION 08/MAR/2022
C     READS CHARACTERISTICS OF THE RUN: # OF PHASES, NAMES OF INPUT FILES,
C     DEFORMATION TO BE IMPOSED, CONVERGENCE PARAMETERS, ETC.
C     READS SINGLE CRYSTAL PROPERTIES: DEFORMATION MODES, CRSS, HARDENING
C     READS CRYSTAL AND MORPHOLOGIC TEXTURES.
C     INITIALIZES ARRAYS REQUIRED TO RUN VPSC.
C     OPEN & CLOSE INPUT FILES  ***  OPEN OUTPUT FILES.
C
C *******************************************************************************

      SUBROUTINE VPSC_INPUT

      USE VPSC8DIM

      DIMENSION AUX3(3)

C *********   INITIALIZATION BLOCK   ***************************

      PI=4.*ATAN(1.)

      XID3(:,:)=0.
      DO I=1,3
        XID3(I,I)=1.
      ENDDO

      XID5(:,:)=0.
      DO I=1,5
        XID5(I,I)=1.
      ENDDO

      XID6(:,:)=0.
      DO I=1,6
        XID6(I,I)=1.
      ENDDO

C     CALCULATES TENSORS OF THE SYMMETRIC BASIS 'B(3,3,6)'
      CALL CHG_BASIS(AUX5,AUX33,AUX55,AUX3333,0,5)

C     SEED FOR RANDOM NUMBER GENERATOR (random2) (USED FOR TWINNING AND RX)
      JRAN=-1

C *** INITIALIZE GAUSS-LEGENDRE OR GAUSS-LOBATO COORDINATES AND ARRAYS USED
C     IN THE DOUBLE INTEGRALS GIVING THE ESHELBY TENSORS

      CALL ESHELBY(AUX3,AUX3333,0.,AUX3333,AUX3333,AUX33,AUX33,PDIL,
     #                  AUX3333,AUX3333,0)

C ***********************************************************************
C     WRITES 'VPSC.IN' FILE INTO 'RUN_LOG.OUT' FILE
      WRITE(10,*)
      WRITE(10,'('' **** INPUT FILE VPSC.IN FOR THIS RUN ****'')')
      DO IDUM=1,100
        READ(UNIT=UR0,END=10,FMT='(A)') PROSA
        WRITE(10,'(A)') PROSA
      ENDDO
   10 REWIND UR0
      WRITE(10,*)
C ***************************************************************************
C     READS REGIME IDENTIFIER: IREGIME=1 (VISCOPLAST); IREGIME=-1 (TH-ELAST)
C     HARDWIRES THE # OF ELEMENTS: NELEM=1
C     READS NUMBER OF PHASES AND RELATIVE VOLUME OF PHASE IN EACH ELEMENT.
C     PHASES ARE LABELED SEQUENTIALLY FOR A MULTIELEMENT RUN (i.e. IF ELEMENT
C     #1 CONSISTS OF phase1=fcc AND phase2=bcc, THEN phase3 AND phase4 ARE
C     THE fcc AND bcc PHASES IN ELEMENT #2...AND SO ON)

      NELEM=1
      READ(UR0,*) IREGIME
      IF(IREGIME.EQ.-1) INTERACTION=-1      ! ELASTIC CASE
      READ(UR0,*) NPH
      READ(UR0,*) (WPH(I),I=1,NPH)
      IPHBOT=1
      IPHTOP=NPH

      write(*,*) '***************************************************'
      if(nelem.eq.1) write(*,'('' SINGLE ELEMENT CALCULATION'')')
      if(nelem.gt.1) write(*,'('' MULTIELEMENT CALCULATION'')')
      if(nph.eq.1)   write(*,'('' SINGLE PHASE CALCULATION'')')
      if(nph.gt.1)   write(*,'('' MULTIPHASE CALCULATION'')')
      write(*,*) '***************************************************'

      if(nelem.gt.1) then
        write(*,'('' VPSC8 is programmed to implement multiple'')')
        write(*,'('' elements run but current version does not'')')
        write(*,'('' include such option ! '')')
        stop
      endif

      if(nph*nelem .gt. nphmx) then
        write(*,'('' number of phases exceeds multiphase dimens !!'')')
        write(*,'('' --> increase parameter NPHMX to'',i5)') nph*nelem
        stop
      endif

      if(nelem.gt.nelmx) then
        write(*,'('' number of elements exceeds multielement dim !!'')')
        write(*,'('' --> increase parameter NELMX to'',i5)') nelem
        stop
      endif

C     CHECKS PHASE FRACTIONS AND NORMALIZES

      wphtot=0.
      do iph=1,nph
        wphtot=wphtot+wph(iph)
      enddo
      if(abs(wphtot-1.) .gt. 1.e-5) then
        write(*,*) ' --> vol fraction of phases should add to 1 !!'
        write(*,*) ' --> will normalize phase fractions to 1 !!'
        print *, 'enter c to continue'
        read  *
        do iph=1,nph
          wph(iph)=wph(iph)/wphtot
        enddo
      endif

C *** INITIALIZE DEFORMATION GRADIENT (i.e. INITIAL ELLIPSOID) OF THE ELEMENT
      NGR(0)=0
      ISHAPE(0)=0
      IFRAG(0) =0
      FREEZE_SHAPE(0)=0
      FIJPH(:,:,0)=XID3(:,:)
	  
      CALL UPDATE_SHAPE (0)      ! AVERAGE ELLIPSOID SHAPE INITIALIZATION

C ****************************************************************************
C     LOOPS OVER PHASES AND FOR EACH PHASE READS ELASTIC & THERMAL MODULI,
C     TEXTURE, GRAIN SHAPE, SLIP SYSTEMS & HARDENING PARAMETERS. 
C     THE LAST TWO ARE NOT REQUIRED FOR A THERMO-ELASTIC SIMULATION.
C ****************************************************************************

      DO IPH=1,NPH

        READ(UR0,'(a)') PROSA
          ISHAPE(IPH)  =0
          IFRAG(IPH)   =0
          CRIT_SHP(IPH)=25
          FREEZE_SHAPE(IPH)=0
        READ(UR0,*)   ISHAPE(IPH),IFRAG(IPH),CRIT_SHP(IPH)
        READ(UR0,*)   (AXISPH(0,I,IPH),I=1,3)
        READ(UR0,*)   (EULERPH(I,IPH) ,I=1,3)

C *** CHECKS IF ELLIPSOID SHAPE EXCEEDS THE CRITICAL ASPECT RATIO		
        AXMAX=0.
		AXMIN=1.E10
		DO I=1,3
		  IF(AXISPH(0,I,IPH).GT.AXMAX) AXMAX=AXISPH(0,I,IPH)
		  IF(AXISPH(0,I,IPH).LT.AXMIN) AXMIN=AXISPH(0,I,IPH)
	    ENDDO
		IF(AXMAX/AXMIN.GE.CRIT_SHP(IPH)) THEN
          print *,' MAXIMUM ASPECT RATIO EXCEEDED FOR PHASE # ',IPH
          print *,' RESET ELLIPSOID ASPECT RAT =< ',CRIT_SHP(IPH)
	      STOP
		ENDIF

        READ(UR0,'(a)') PROSA
        READ(UR0,'(a)') FILETEXT
		
        READ(UR0,'(a)') PROSA
        READ(UR0,'(a)') FILECRYS

        READ(UR0,'(a)') PROSA
        READ(UR0,'(a)') FILEAXES      ! DUMMY IF ISHAPE.LE.2

        READ(UR0,'(a)') PROSA
        READ(UR0,*)     IDIFF(IPH)
        READ(UR0,'(a)') FILEDIFF       ! DUMMY IF IDIFF=0

C **********************************************************************
C *** READS ELASTIC AND THERMAL TENSORS Cij AND Ai FOR THE PHASE
C *** READS SLIP AND TWINNING SYSTEMS FOR THE PHASE

        OPEN (UNIT=UR1,FILE=FILECRYS,STATUS='OLD')

        CALL DATA_CRYSTAL (IPH)

C ****************************************************************************
        IF (INTERACTION.NE.-1) THEN         ! visco-plastic case
C ****************************************************************************

C *** READS HARDENING PARAMETERS OF THE HARDENING LAW FOR THIS PHASE
C       IHARDLAW=0  --> EXTENDED VOCE & PREDOMINANT TWIN REORIENTATION
C       IHARDLAW=1  --> MTS MODEL
C       IHARDLAW=2x --> DISLOCATION DENSITY & DD reversal HARDENING
C       IHARDLAW=3x --> IRRADIATION GROWTH PLUS IRRADIATION CREEP   

        READ(UR1,*)
        READ(UR1,*) IHARDLAW
        READ(UR1,*) IRATESENS
        IF(IHARDLAW.EQ.1 .AND. IRATESENS.EQ.1) THEN
          WRITE(*,*) ' FOR MTS MODEL RESETS IRATESENS=0 IN CONST LAW'
          IRATESENS=0
          print *, 'enter c to continue'
          read  *
        ENDIF
        IF(IHARDLAW/10.EQ.2 .AND. IRATESENS.EQ.1) THEN
          WRITE(*,*) ' FOR DD MODEL RESETS IRATESENS=0 IN CONST LAW'
          IRATESENS=0
          print *, 'enter c to continue'
          read  *
        ENDIF

        IF(IHARDLAW.EQ.0)  CALL UPDATE_CRSS_VOCE (0,IPH)
        IF(IHARDLAW.EQ.1)  CALL UPDATE_CRSS_MTS (EDOTX,0)
        IF(IHARDLAW.EQ.20) CALL UPDATE_CRSS_DD (0)
        IF(IHARDLAW.EQ.23) CALL UPDATE_CRSS_DD_REV (0)

        IF(IHARDLAW.EQ.30) CALL UPDATE_CRSS_VOCE(0,IPH)
        IF(IHARDLAW.EQ.30) CALL UPDATE_GROWTH_RATE(0)         ! CT
        IF(IHARDLAW.EQ.31) CALL UPDATE_GROWTH_ZIRC2(0,IPH)    ! AP
        IF(IHARDLAW.EQ.32) CALL UPDATE_GROWTH_ZR(0,IPH)       ! AP

        CLOSE(UNIT=UR1)

        IF(ITWINLAW.EQ.3) THEN
          IF(NPH.NE.2 .OR. WPH(2).NE.0.) THEN 
          WRITE(*,*) ' VFT TWIN SCHEME REQUIRES A PHASE2 WITH WPH(2)=0'
          STOP
          ENDIF
        ENDIF

C ****************************************************************************
       ENDIF     ! end of IF (INTERACTION.NE.-1)
C ****************************************************************************

C *** READS INITIAL TEXTURE (BUNGE.ROE,KOCKS CONVENTIONS) FROM 'FILETEXT'
C *** INITIALIZES DEFORMATION GRADIENT AND ELLIPSOID FOR EACH PHASE.
C *** IF ISHAPE=3 READS INITIAL GRAIN AXES AND ORIENTATION FROM FILEAXES
C     AND INITIALIZES DEF GRADIENT & ELLIPSOID FOR INDIVIDUAL GRAINS.

        OPEN(UNIT=UR2,FILE=FILETEXT,STATUS='OLD')
        IF(ISHAPE(IPH).EQ.3) OPEN(UNIT=UR3,FILE=FILEAXES,STATUS='OLD')
          CALL DATA_GRAIN(IPH)
        CLOSE(UNIT=UR2)
        IF(ISHAPE(IPH).EQ.3) CLOSE(UNIT=UR3)

        FREEZE_SHAPE(IPH)=0   ! DEFAULT SETTING IS TO UPDATE PHASE SHAPE
        CALL UPDATE_SHAPE (IPH)

C **********************************************************************
C *** IF IDIFF=1 READS DIFFRACTION DIRECTIONS AND CRYSTAL PLANES ALONG 
C     THEM FOR CALCULATING AVERAGE ELASTIC STRAIN ON THOSE PLANES.

        IF(IDIFF(IPH).EQ.1) THEN 
          OPEN(UNIT=UR7,FILE=FILEDIFF,STATUS='OLD')
          CALL DIFF_PLANES (IPH,ICRYSYMPH(IPH),0)
          CLOSE(UNIT=UR7)
        ENDIF

      ENDDO     ! END OF DATA INPUT LOOP DO IPH=1,NPH

C ***************************************************************************

      IF(NGR(NPH).GT.NGRPEL) THEN
        WRITE(*,*)
        WRITE(*,'('' --> NUMBR OF GRAINS PER ELEMENT IS'',I6)') NGR(NPH)
        WRITE(*,'('' --> INCREASE PARAMETER NGRPEL IN VPSC8.DIM'')')
        STOP
      ENDIF

      NGTOT=NELEM*NGR(NPH)
      WRITE(*,'('' --> TOTAL NUMBER OF GRAINS IS'',I6)') NGTOT
      IF(NGTOT.GT.NGRMX) THEN
        WRITE(*,'('' --> INCREASE PARAMETER NGRMX IN VPSC8.DIM'')')
        STOP
      ENDIF

C *****************************************************************************
C *** IF INTERACTION=-1 (THERMO-ELASTIC) THE PLASTICITY RELATED PARAMETERS ARE
C     NOT USED.  READS NEXT 14 LINES IN VPSC.IN AS DUMMY AND RETURNS TO MAIN.

      IF (INTERACTION.EQ.-1) THEN
        DO IDUMMY=1,14
          READ(UR0,*)
        ENDDO
        RETURN
      ENDIF

C ****************************************************************************
C *** READS SETTINGS FOR CONVERGENCE PROCEDURES.

      READ(UR0,'(A)') PROSA
      READ(UR0,*) ERRS,ERRD,ERRM,ERRSO
      READ(UR0,*) ITMAXEXT,ITMAXINT,ITMAXSO
      IRSVAR=0      ! this is the default 
      READ(UR0,*) IRSVAR,JXRSINI,JXRSFIN,JXRSTEP
CCC      READ(UR0,*) IBGINV      ! eliminated 23/07/2019

C ****************************************************************************
C *** READS INPUT/OUTPUT SETTINGS FOR THE RUN

      IRECOVER=0      ! default values
      ISAVE   =0
      ICUBCOM =0
      NWRITE  =0

      READ(UR0,'(A)') PROSA
      READ(UR0,*) IRECOVER
      READ(UR0,*) ISAVE
      READ(UR0,*) ICUBCOM
      READ(UR0,*) NWRITE

C ****************************************************************************
C *** READS MODELING CONDITIONS

      IUPDORI=1      ! default values
      IUPDSHP=1
      IUPDHAR=1
      NNEIGH =0
      IFLU   =0

      READ(UR0,'(A)') PROSA
      READ(UR0,*) INTERACTION, NEFFGRX
      READ(UR0,*) IUPDORI,IUPDSHP,IUPDHAR
      READ(UR0,*) NNEIGH
      READ(UR0,*) IFLU
	  
      IF(IUPDORI.EQ.0 .AND. IUPDSHP.EQ.1) THEN
        WRITE(*,'('' --> IUPDSHP=1 AND IUPDORI=0'')')
        WRITE(*,'('' --> CANNOT UPDATE SHAPE WITHOUT ALSO UPDATING'',
     #            '' ORIENTATION'')')
        STOP
      ENDIF

      IF(NNEIGH.GT.NNEIMX) THEN
        WRITE(*,'('' --> NNEIGH='',I3,'' > NNEIMX='',I3)') NNEIGH,NNEIMX
        WRITE(*,'('' --> UPDATE NNEIMX IN VPSC8DIM'')')
        STOP
      ENDIF

C *** CHECKS PRECISION CONTROL FOR SO CALCULATIONS
      IF(INTERACTION.EQ.5) THEN
        IF(ERRS.GT.1.E-7 .OR. ERRD.GT.1.E-7 .OR.ERRM.GT.1.E-7) THEN
          WRITE(*,'('' SO CALCULATIONS MAY NOT CONVERGE IF TOLERANCES'',
     #              '' ERRS,ERRD,ERRM ARE LARGER THAN 1.E-7'')')
          print *, 'enter c to continue'
          read  *
        ENDIF
        IF(ERRSO.LT.1.E-2) THEN
          WRITE(*,'('' SO CALCULATIONS MAY NOT CONVERGE IF TOLERANCE'',
     #              '' ERRSO IS SMALLER THAN 1.E-2'')')
          print *, 'enter c to continue'
          read  *
        ENDIF
      ENDIF

C **********************************************************************
C *** CHECKS IF RATE SENSITIVITY IS THE SAME FOR ALL SYSTEMS.
C *** SEARCH FOR NRSMIN (NEEDED TO GET TAUMAX INSIDE NR SUBROUTINE)

      NUNIQUE=1
      NCOMPA=NRS(1,1)
      NRSMIN=NRS(1,1)
      DO IPH=1,NPH
        DO IS=1,NSYST(IPH)
          IF(NRS(IS,IPH).NE.NCOMPA) NUNIQUE=0
          IF(NRS(IS,IPH).LT.NRSMIN) NRSMIN=NRS(IS,IPH)
        ENDDO
      ENDDO
      IF(NUNIQUE.EQ.0) THEN
        WRITE(*,'('' --> THIS RUN CONTAINS MIXED nrs EXPONENTS !!'')')
        IF(INTERACTION.NE.1) THEN
          WRITE(*,'('' --> WILL FORCE AFFINE (INTERACTION=1)'')')
          INTERACTION=1
          print *, 'enter c to continue'
          read  *
        ENDIF
      ENDIF

C **********************************************************************
C *** DEFINES PARAMETER 'NEFF' THAT ENTERS IN THE INTERACTION EQUATION &
C     CONTROLS THE STRENGTH OF THE GRAIN INTERACTION WITH THE EFF MEDIUM.
C *** WHEN RUNNING RelDirCompl USES NRS/2 FOR ALL GRAINS ONLY IN 1ST STEP.

      IF(INTERACTION.EQ.3) THEN
        RDCrun=.FALSE.
        IF(NEFFGRX.EQ.0.) THEN
          WRITE(*,'('' RUNNING RelDirCompl '')')
          NEFFGRX=NRS(1,1)/2.      ! the 1st step is done with uniform neff
          RDCrun=.TRUE.
        ENDIF
      ENDIF
      WRITE(10,'(/,''*** RUNNING INTERACTION NEFF='',F8.2,/)') NEFFGRX
      DO IPH=1,NPH
      DO KKK=NGR(IPH-1)+1,NGR(IPH)
        IF(INTERACTION.EQ.1) NEFFGR(KKK)=1.            !affine
        IF(INTERACTION.EQ.2) NEFFGR(KKK)=1.            !secant
        IF(INTERACTION.EQ.3) NEFFGR(KKK)=NEFFGRX       !neff & RDC
        IF(INTERACTION.EQ.4) NEFFGR(KKK)=NRS(1,1)      !tangent
        IF(INTERACTION.EQ.5) NEFFGR(KKK)=1.            !2nd order
      ENDDO
      ENDDO

C **********************************************************************

      if(interaction.eq.5.and.iflu.eq.0) then
        write(*,*) 'FLUCTUATIONS ARE NEEDED FOR 2nd-ORDER CALCULATION'
        write(*,*) 'IFLU IS RESET TO 1'
        iflu=1
      endif

      if(interaction.eq.0.and.iflu.eq.1) then
        write(*,*)
     #  'INTRAGRANULAR FLUCTS CANNOT BE OBTAINED IN FC CALCULATION'
        write(*,*) 'IFLU IS RESET TO 0'
        iflu=0
      endif

      if(iflu.eq.1) then
        open(83,file='FLUCT.OUT',status='unknown')
        if(icubcom.eq.1) open(84,file='FLCUB.OUT',status='unknown')
      endif

      if(interaction.eq.5) then
        open(97,file='SO.OUT',status='unknown')
      endif

C ****************************************************************************
C *** INITIALIZE NEIGHBOURS TABLE IF NNEIGH IS NOT ZERO.
C     IOPTION=0 --> RANDOM PAIRING / IOPTION=1 --> SEQUENTIAL PAIRING

      IOPTION=0    ! random pairing option hardwired
      CALL NEIGHBOURS (IOPTION)

C *******************************************************************
C *** INITIALIZE ARRAYS ASSOCIATED WITH GRAINS:
C     HARDENING, ACCUMULATED SHEAR, POWER, TWINNING PARAMETERS, etc

      DO IPH=1,NPH

        PRITW(IPH)=0.
        SECTW(IPH)=0.
        SCH_STAT(:,:,IPH)=0.
	    VAR_STAT(:,:,IPH)=0.
        DO KKK=NGR(IPH-1)+1,NGR(IPH)
          EPSVMGR(KKK)=0.
          WORKGR(KKK) =0.
        ENDDO

        DO KKK=NGR(IPH-1)+1,NGR(IPH)
          GTOTGR(KKK)=0.
          DO IS=1,NSYST(IPH)
            GAMD0S(IS,KKK)=1.0
          ENDDO

          KTWSMX(KKK)=0
          NTWEVENTS(KKK)=0
          DO ITS=1,NTWSYS(IPH)
            TWFRSY(ITS,KKK)=0.
          ENDDO
        ENDDO

        DO ITM=1,NTWMOD(IPH)
          EFTWFR(ITM,IPH)=0.
          TWFRPH(ITM,IPH)=0.
        ENDDO

C *** UNLESS IHARDLAW=30 TRANSFORMATION/GROWTH RATE DG_TRANS(I,KKK)=0.
        DO KKK=NGR(IPH-1)+1,NGR(IPH)
          DO J=1,5
            DG_TRANS(J,KKK)=0.d0
          ENDDO
        ENDDO

      ENDDO      ! END OF DO IPH=1,NPH

C *** INITIALIZE CRSS IN EACH SYSTEM OF EACH GRAIN USING PARAMETERS READ
C     FROM SX FILE.
C *** FOR 'MTS', 'DD' AND 'DD_REV' ONLY INITIALIZE STRUCTURE-RELATED ARRAYS
C     NOT DEPENDING ON PROCESS VARIABLES (TEMP & RATE). 

      IF(IHARDLAW.EQ.0)  CALL UPDATE_CRSS_VOCE (1,IDUMMY)
      IF(IHARDLAW.EQ.1)  CALL UPDATE_CRSS_MTS (EDOTX,1)
C --> initialize DD parameters after the process and its temperature is read 
C      IF(IHARDLAW.EQ.20) CALL UPDATE_CRSS_DD (1) 
      IF(IHARDLAW.EQ.23) CALL UPDATE_CRSS_DD_REV (1)

      IF(IHARDLAW.EQ.30) CALL UPDATE_CRSS_VOCE (1,IDUMMY)
      IF(IHARDLAW.EQ.30) CALL UPDATE_GROWTH_RATE(1)
      IF(IHARDLAW.EQ.31) CALL UPDATE_GROWTH_ZIRC2(1,IDUMMY)
      IF(IHARDLAW.EQ.32) CALL UPDATE_GROWTH_ZR(1,IDUMMY)

C ***************************************************************************
C *** WHEN THERE IS MORE THAN ONE ELEMENT COPIES CRYSTAL, GRAINS & TEXTURE
C     PARAMETERS FROM ELEMENT #1 INTO THE OTHER ELEMENTS.
C *** THE NUMERICAL STRATEGY IS TO USE THE 'PHASE' DIMENSION TO ASSIGN ARRAYS
C     TO THE ELEMENTS:
C     * FOR SINGLE PHASE: GRAINS IN PHASE(N) ARE ASSIGNED TO ELEMENT(N)
C     * FOR 2-PHASE: GRAINS IN PHASES (2N-1) & (2N) ARE ASSIGNED TO ELEMEN(N)
C     * FOR 3-PHASE: GRAINS IN PHASES (3N-2),(3N-1) & (3N) ARE ASSIGNED TO
C       ELEMENT(N)
C ***************************************************************************

      IF(NELEM.GE.2) THEN

      DO IELEM=2,NELEM
        DO IPH=1,NPH
          IPHAC=(IELEM-1)*NPH+IPH

C *** PARAMETERS OF EACH PHASE IN THE ELEMENT

          WPH(IPHAC)=WPH(IPH)
          ISHAPE(IPHAC)=ISHAPE(IPH)
          FREEZE_SHAPE (IPHAC)=FREEZE_SHAPE (IPH)
          IFRAG (IPHAC)=IFRAG (IPH)          ! REDUNDANT DIMENSIONING
          CRIT_SHP(IPHAC)=CRIT_SHP(IPH)      ! REDUNDANT DIMENSIONING
          DO I=1,3
            AXISPH(0,I,IPHAC)=AXISPH(0,I,IPH)
            EULERPH(I,IPHAC) =EULERPH(I,IPH)
            DO J=1,3
              FIJPH(I,J,IPHAC) =FIJPH(I,J,IPH)
              AXISPH(I,J,IPHAC)=AXISPH(I,J,IPH)
            ENDDO
          ENDDO
          NGR(IPHAC)=NGR(IPHAC-1)+(NGR(IPH)-NGR(IPH-1))

          PRITW(IPHAC)=0.
          SECTW(IPHAC)=0.
          DO ITM=1,NTWMOD(IPH)
            EFTWFR(ITM,IPHAC)=0.
            TWFRPH(ITM,IPHAC)=0.
          ENDDO

C *** PARAMETERS OF EACH GRAIN IN THE PHASE

          DO KGX=NGR(IPH-1)+1,NGR(IPH)

            KKK=NGR(NPH)*(IELEM-1)+KGX
            WGT(KKK)   =WGT(KGX)
            GTOTGR(KKK)=0.
            DO IS=1,NSYST(IPH)
              GAMD0S(IS,KKK)=1.0
            ENDDO

            DO I=1,5
              DG_TRANS(I,KKK)=DG_TRANS(I,KGX)   ! cnt
            ENDDO

            DO I=1,3
            DO J=1,3
              AG(I,J,KKK)=AG(I,J,KGX)
            ENDDO
            ENDDO

            DO IS=1,NSYST(IPH)
              CRSS(IS,KKK)=CRSS(IS,KGX)
              TAUE(IS,KKK)=TAUE(IS,KGX)
            ENDDO

            KTWSMX(KKK)=0
            NTWEVENTS(KKK)=0
            DO ITS=1,NTWSYS(IPH)
              TWFRSY(ITS,KKK)=0.
            ENDDO

            EPSVMGR(KKK)=0.
            WORKGR(KKK) =0.

            IF(ISHAPE(IPH).NE.0) THEN
              DO I=1,3
                AXISGR(0,I,KKK)=AXISGR(0,I,KGX)
                DO J=1,3
                  FIJGR(I,J,KKK) =FIJGR(I,J,KGX)
                  AXISGR(I,J,KKK)=AXISGR(I,J,KGX)
                ENDDO
              ENDDO
            ENDIF

          ENDDO      ! END OF LOOP OVER KGX
        ENDDO      ! END OF LOOP OVER IPH
      ENDDO      ! END OF LOOP OVER IELEM
	  
      ENDIF      ! END OF IF(NELEM.GE.2)

      RETURN
      END

C **********************************************************************
C     SUBROUTINE WRITE_SHEAR_ACTIVITY    --->    VERSION 30/APR/2022
C
C     WRITES ACTIVITY STATISTICS OVER DEFORMATION MODES IN 'ACT_PHn.OUT'
C     ACTIVITIES ARE CALCULATED INSIDE SUBR STAT_SHEAR_ACTIVITY 
C **********************************************************************

      SUBROUTINE WRITE_SHEAR_ACTIVITY (istep)

      USE VPSC8DIM

      CHARACTER*7 MODE(NMODMX),TWFR(NMODMX),EFFR(NMODMX),VFPH(NMODMX)

      DO IPH=IPHBOT,IPHTOP
        IPHEL=IPH-IPHBOT+1
        IUNIT=50+IPHEL

        DO I=1,NMODES(IPHEL)
          MODE(I)='  MODE'//CHAR(48+I)
          TWFR(I)='  TWFR'//CHAR(48+I)
          EFFR(I)='  EFFR'//CHAR(48+I)
          VFPH(I)='  VFPH'//CHAR(48+I)
        ENDDO

C *** for a PCYS or LANKFORD calculation

        IF(IVGVAR.GE.2) THEN
          if(istep.eq.1) write(IUNIT,'(''    STEP   AVACS'',
     #                    10A8)') (MODE(I),I=1,NMODES(IPHEL))
          write(IUNIT,'(i8,12f8.3)') istep,ACTPH(iph),
     #                   (gavmod(kmo,iph),kmo=1,nmodes(IPHEL))
        ENDIF

C *** for loading histories

        IF(IVGVAR.LE.1) THEN

          IF (NTWMOD(IPHEL).EQ.0) THEN
            if(istep.eq.1) write(IUNIT,'(''   STRAIN    AVACS'',10a9)')
     #                     (MODE(I),I=1,NMODES(IPHEL))
            write(IUNIT,'(12f9.4)') epsacu,ACTPH(iph),
     #                   (gavmod(kmo,iph),kmo=1,nmodes(iphel))

          ELSEIF (NTWMOD(IPHEL).NE.0) THEN

            itm1=nmodes(iphel)-ntwmod(iphel)+1
            itm2=nmodes(iphel)

            IF(IHARDLAW.EQ.0) THEN

            IF(ITWINLAW.LE.2) THEN
              IF(istep.eq.1)
     #        write(IUNIT,'(''   STRAIN    AVACS    PRITW    SECTW'',
     #                     30a9)') (MODE(I),I=1,NMODES(IPHEL)),
     #                     (TWFR(J),EFFR(J),j=itm1,itm2)
              write(IUNIT,'(34f9.4)') epsacu,ACTPH(iph),
     #                               pritw(iph),sectw(iph),
     #                   (gavmod(kmo,iph),kmo=1,nmodes(iphel)),
     #         (twfrph(jmo,iph),eftwfr(jmo,iph),jmo=1,ntwmod(iphel))
            ENDIF

            IF(ITWINLAW.EQ.3) THEN
              IF(istep.eq.1)
     #        write(IUNIT,'(''   STRAIN    AVACS    PRITW    SECTW'',
     #                     30a9)') (MODE(I),I=1,NMODES(IPHEL)),
     #                     (TWFR(J),j=itm1,itm2),(VFPH(K),K=1,2)
              write(IUNIT,'(34f9.4)') epsacu,ACTPH(iph),
     #                               pritw(iph),sectw(iph),
     #                   (gavmod(kmo,iph),kmo=1,nmodes(iphel)),
     #         (twfrph(jmo,iph),jmo=1,ntwmod(iphel)),WPH(1),WPH(2)
            ENDIF

            ENDIF

            IF(IHARDLAW.EQ.20 .or. IHARDLAW.EQ.23) THEN                                     
            if(istep.eq.1)    
     #      write(IUNIT,'(''   STRAIN    AVACS    PRITW    SECTW'',   
     #                   30a9)') (MODE(I),I=1,NMODES(IPHEL)),                    
     #                   (TWFR(J),EFFR(J),j=itm1,itm2)
            write(IUNIT,'(34f9.4)') epsacu,ACTPH(iph),      
     #                              pritw(iph),sectw(iph),          
     #                   (gavmod(kmo,iph),kmo=1,nmodes(iphel)),      
     #         (twfrph(jmo,iph),eftwfr(jmo,iph),jmo=1,ntwmod(iphel))   
            ENDIF       
			
          ENDIF

        ENDIF      ! END OF IF(IVGVAR.LE.1)

      ENDDO

      RETURN
      END
C
C *****************************************************************
C     SUBROUTINE WRITE_STRESS_STRAIN    -->   VERSION 26/APR/2022
C *****************************************************************

      SUBROUTINE WRITE_STRESS_STRAIN (ISTEP)

      USE VPSC8DIM

      IU=13      ! OUTPUT UNIT (file STR_STR.OUT)
	  
      IDILAT=IDBARv(1)*IDBARv(2)*IDBARv(3) 

      IF( (ISTEP.EQ.1 .AND. INTERACTION.NE.-1) .OR.  
     #    (ISTEP.EQ.0 .AND. INTERACTION.EQ.-1) ) THEN
	  
        IF(IDILAT.EQ.1) THEN                      ! CAUCHY STRESS NOT AVAILABLE
c *** for loading path		
        IF(IVGVAR.LE.1) WRITE (IU,'(10x,''Evm'',10x,''Svm'',
     #    10x,''E11'',10x,''E22'',10x,''E33'',
     #    10x,''E23'',10x,''E13'',10x,''E12'',
     #    13x,''SDEV11'',7x,''SDEV22'',7x,''SDEV33'',
     #     7x,''SDEV23'',7x,''SDEV13'',7x,''SDEV12'',
     #     ''      TEMP      TIME    TAYLAV      WORKAV'',
     #     ''      WRATE1      WRATE2 '')')
c *** for PCYS & Lankford		
        IF(IVGVAR.EQ.2) WRITE (IU,'(10x,''Dvm'',10x,''Svm'',
     #    10x,''D11'',10x,''D22'',10x,''D33'',
     #    10x,''D23'',10x,''D13'',10x,''D12'',
     #    13x,''SDEV11'',7x,''SDEV22'',7x,''SDEV33'',
     #     7x,''SDEV23'',7x,''SDEV13'',7x,''SDEV12'',
     #     ''      TEMP    TAYLAV      WORKAV'',
     #     ''      WRATE1      WRATE2 '')')
        ELSE IF(IDILAT.EQ.0) THEN                  ! CAUCHY STRESS AVAILABLE
c *** for loading path
        IF(IVGVAR.LE.1.and.STRAIN_CONTROL.EQ.1) WRITE (IU,'(10x,''Evm'',
     #    10x,''Svm'',10x,''E11'',10x,''E22'',10x,''E33'',
     #                10x,''E23'',10x,''E13'',10x,''E12'',
     #                13x,''SCAU11'',7x,''SCAU22'',7x,''SCAU33'', 
     #                 7x,''SCAU23'',7x,''SCAU13'',7x,''SCAU12'',
     #     ''      TEMP      TIME    TAYLAV      WORKAV'', 
     #     ''      WRATE1      WRATE2 '')')
        IF(IVGVAR.LE.1.and.STRAIN_CONTROL.EQ.0) WRITE (IU,'(10x,''Evm'',
     #    10x,''Dvm'',10x,''E11'',10x,''E22'',10x,''E33'', 
     #                10x,''E23'',10x,''E13'',10x,''E12'',
     #                13x,''DBAR11'',7x,''DBAR22'',7x,''DBAR33'',
     #                 7x,''DBAR23'',7x,''DBAR13'',7x,''DBAR12'',
     #     ''      TEMP      TIME    TAYLAV      WORKAV'', 
     #     ''      WRATE1      WRATE2 '')')
c *** for PCYS & Lankford
        IF(IVGVAR.GE.2) WRITE (IU,'(10x,''Dvm'',10x,''Svm'',
     #    10x,''D11'',10x,''D22'',10x,''D33'',
     #    10x,''D23'',10x,''D13'',10x,''D12'',
     #    13x,''SCAU11'',7x,''SCAU22'',7x,''SCAU33'',
     #     7x,''SCAU23'',7x,''SCAU13'',7x,''SCAU12'',
     #     ''      TEMP    TAYLAV      WORKAV'', 
     #     ''      WRATE1      WRATE2 '')')

        ENDIF

      ENDIF

      IF(IVGVAR.LE.1.and.STRAIN_CONTROL.EQ.1) WRITE(IU,'(2e13.5,
     #   3x,6e13.5,3x,6e13.5,f10.0,e10.3,f10.3,3e12.4)')
     #   EPSVM,SVM,EPSTOTv(:),SBARv(:),
     #   TEMP,TIME,TAYLAV,WORKAV,WRATE1,WRATE2
      IF(IVGVAR.LE.1.and.STRAIN_CONTROL.EQ.0) WRITE(IU,'(2e13.5,
     #   3x,6e13.5,3x,6e13.5,f10.0,e10.3,f10.3,3e12.4)')
     #   EPSVM,DVM,EPSTOTv(:),DBARv(:),
     #   TEMP,TIME,TAYLAV,WORKAV,WRATE1,WRATE2

      IF(IVGVAR.GE.2) WRITE(IU,'(2e13.5,3x,6e13.5,3x,6e13.5,f10.0,
     #   f10.3,3e12.4)') DVM,SVM,DBARv(:),SBARv(:),     
     #                  TEMP,TAYLAV,WORKAV,WRATE1,WRATE2

C *****************************************************************************
C *** THIS OUTPUT WAS PREVIOUSLY COMING FROM SUBROUTINE WRITE_STAT (deprecated)
C       WRITES IN FILE STR_STR_STAT.OUT THE STANDARD DEVIATION OF STRESS
C       AND STRAIN-RATE COMPONENTS IN b-BASIS, CALCULATED INSIDE SUBROUTINE
C       STAT_STRESS_STRAIN.
C     FOR A THERMO-ELASTIC CASE WRITES CAUCHY AND ELASTIC STRAIN COMPONENTS
C *****************************************************************************

      IU=11      ! OUTPUT UNIT (file STR_STR_STAT.OUT)

      IF(ISTEP.EQ.1 .AND. INTERACTION.NE.-1) WRITE(IU,'(
     #   ''        EPSvm          Svm'', 4x,
     #   ''    SDEV1    SDEV2    SDEV3    SDEV4    SDEV5'',
     #   ''    SDEV6'', 4x,
     #   ''    DDEV1    DDEV2    DDEV3    DDEV4    DDEV5'',
     #   ''    DDEV6'')')

      IF(ISTEP.EQ.0 .AND. INTERACTION.EQ.-1) WRITE(IU,'(
     #   ''        EPSvm          Svm'',4x,
     #   ''    SDEV1    SDEV2    SDEV3    SDEV4    SDEV5'',
     #   ''    SDEV6'', 4x,
     #   ''    EDEV1    EDEV2    EDEV3    EDEV4    EDEV5'',
     #   ''    EDEV6'')')

      IF(IVGVAR.LE.1 .AND. INTERACTION.NE.-1)
     #          WRITE(IU,'(2E13.5,4x,6f9.3,4x,6f9.3)') 
     #          EPSVM,SVM,STDEVS(:),STDEVD(:)
      IF(IVGVAR.LE.1 .AND. INTERACTION.EQ.-1)
     #          WRITE(IU,'(2E13.5,4x,6f9.3,4x,6f9.3)') 
     #          EPSVM,SVM,STDEVS(:),STDEVE(:)
      IF(IVGVAR.GE.2) WRITE(IU,'(2E13.5,4x,6f9.3,4x,6f9.3)') 
     #          DVM,   SVM,STDEVS(:),STDEVD(:)
	 
      RETURN
      END

C **************************************************************************
C     SUBROUTINE WRITE_TEXTURE      --->     VERSION 05/JUNE/2021
C
C     WRITES CRYSTALLOGRAPHIC TEXTURE FOR EACH PHASE.
C     IF ISHAPE(IPH) >0 ALSO WRITES MORPHOLOGIC TEXTURE FOR EACH PHASE.
C **************************************************************************

      SUBROUTINE WRITE_TEXTURE

      USE VPSC8DIM

      DIMENSION AGT(3,3),EULX(3,2*NGRPEL),WGTX(2*NGRPEL)

      DO IPH=IPHBOT,IPHTOP
        IPHEL=IPH-IPHBOT+1
		NGRX=1
        DO KKK=NGR(IPH-1)+1,NGR(IPH)
          DO I=1,3
          DO J=1,3
            AGT(I,J)=AG(J,I,KKK)
          ENDDO
          ENDDO
          CALL EULER(1,GDA,GDB,GDC,AGT)
          EULX(1,NGRX)=GDA
          EULX(2,NGRX)=GDB
          EULX(3,NGRX)=GDC
          WGTX(NGRX)  =WGT(KKK)

          NGRX=NGRX+1
        ENDDO      ! END OF DO KKK
        NGRX=NGRX-1

C *************************************************************************
C     WRITE CRYSTALLOGRAPHIC TEXTURE (CRYSTAL ORIENTATION AND WEIGHT)
C *** STANDARD OUTPUT: WRITE TEXTURE OF PHASE1,PHASE2,etc IN UNITS=31,32,etc
C *** FOR THE VFT SCHEME WRITES PARENT & CHILD PHASES IN A SINGLE FILE.
C     ITWINLAW= 1, 2, 3 for: PTR , MC , VFT scheme, respectively
C *************************************************************************

        IF(ITWINLAW.NE.3) THEN
          IUNIT=30+IPHEL                ! WRITES ONE TEXTURE FILE PER PHASE

          WRITE(IUNIT,'(A,F10.4)') 'TEXTURE AT STRAIN =',EPSACU
          WRITE(IUNIT,'(3f8.3,A)') (AXISPH(0,I,IPHEL),I=1,3),
     #       '  <-- length of average ellipsoid axes for phase'
          WRITE(IUNIT,'(3f8.2,A)') (EULERPH(I,IPHEL),I=1,3),
     #       '  <-- Euler angles of average ellipsoid axes for phase'
          WRITE(IUNIT,'(''B'',I10)') NGRX
          WRITE(IUNIT,'(3F8.2,F12.7)')
     #              ((EULX(I,KG),I=1,3),WGTX(KG),KG=1,NGRX)

        ELSE IF (ITWINLAW.EQ.3) THEN
         IUNIT=31           

         IF(IPHEL.EQ.1) then 
           WRITE(IUNIT,'(A,F10.4)') 'TEXTURE AT STRAIN =',EPSACU
           WRITE(IUNIT,'(3f8.3,A)') (AXISPH(0,I,IPHEL),I=1,3),
     #       '  <-- length of average ellipsoid axes for phase'
           WRITE(IUNIT,'(3f8.2,A)') (EULERPH(I,IPHEL),I=1,3),
     #       '  <-- Euler angles of average ellipsoid axes for phase'
           WRITE(IUNIT,'(''B'',I10)') 2*ngrx
         ENDIF 

           WRITE(IUNIT,'(3F8.2,1F12.7)')
     #               ((EULX(I,KG),I=1,3),WGTX(KG),KG=1,NGRX)

c --> outputs orientations, weight & elastic energy (for quartz)
c          WRITE(IUNIT,'(3F8.2,2F12.7)')
c     #   ((EULX(I,KG),I=1,3),WGTX(KG),W_EL_GR(kg+NGR(IPH-1)),KG=1,NGRX)

      ENDIF       

      ENDDO      ! END OF DO IPH

C *************************************************************************
C     WRITE MORPHOLOGICAL TEXTURE (ELLIPSOID EULER ANGLES AND AXES' LENGTH)
C     --> AXISGR(I,J) TRANSFORMS FROM ELLIPSOID TO SAMPLE AXES
C     --> GDA,GDB,GDC: DESCRIBE ELLIPSOID AXES WITH RESPECT TO SAMPLE AXES
C     --> AVERAGE LENGTH OF AXES AND STANDARD DEVIATION ARE CALCULATED 
C         IN SUBROUTINE STAT_GRAIN_SHAPE
C **************************************************************************

      DO IPH=IPHBOT,IPHTOP

        IF(ISHAPE(IPH).NE.0 .AND. WPH(IPH).NE.0) THEN
          IUNIT=40+IPH
      iskip=0
      if(iskip.eq.0) then
          write(iunit,'(a,f10.4,a,i3)') 
     #           'MORPHOLOGY AT STRAIN =',EPSACU,'  FOR PHASE',IPH
          write(iunit,'(3f8.3,a)') (AXISPH(0,I,IPH),I=1,3) ,
     #           '  <-- AXES OF PHASE ELLIPSOID'
          write(iunit,'(3f8.2,a)') (EULERPH(i,iph),i=1,3),
     #           '  <-- EULER ANGLES OF PHASE ELLIPSOID (deg)'
          write(iunit,'(3f8.3,a)') AVAX(:,IPH),
     #           '  <-- AVERAGE OF GRAIN AXES FOR PHASE'
          write(iunit,'(3f8.3,a)') SDAX(:,IPH),
     #           '  <-- STANDARD DEV OF AVG GRAIN AXES FOR PHASE'
          write(iunit,'(''B'',2i10)') (ngr(iph)-ngr(iph-1))
      endif
          DO KKK=NGR(IPH-1)+1,NGR(IPH)
            DO I=1,3
            DO J=1,3
              AGT(I,J)=AXISGR(J,I,KKK)
            ENDDO
            ENDDO
            CALL EULER (1,GDA,GDB,GDC,AGT)
            write(iunit,'(3f8.2,3x,f10.6,3x,3f10.4)') 
     #            GDA,GDB,GDC,WGT(KKK),(AXISGR(0,I,KKK),I=1,3)
          ENDDO      ! END OF DO KKK

        ENDIF       ! END OF IF(ISHAPE(IPH).NE.0)

      ENDDO      ! END OF DO IPH=IPHBOT,IPHTOP

      RETURN
      END

C *****************************************************************************
C     SUBROUTINE WRITE_TWIN_STATS      --->    VERSION 13/AUG/2018
C
C     WRITES TWINNING STATISTICS CALCULATED IN SUBROUTINE 'STAT_TWINNING'.
C     --> DISTRIBUTION OF SCHMID FACTORS OF THE TWIN SYSTEMS IN EACH MODE  
C         USED TO REORIENT GRAINS BY TWINNING
C     --> DISTRIBUTION OF VARIANTS (BY DECREASING ORDER OF SCHMID FACTORS)
C         ASSOCIATED WITH THE TWIN SYSTEMS USED TO REORIENT GRAINS IN EACH MODE 
C *****************************************************************************

      SUBROUTINE WRITE_TWIN_STATS(ISTEP,IPH)

      USE VPSC8DIM

      INTEGER :: ITWMOD, NGR_STAT
      
C * OPEN OUTPUT STATISTIC FILES IN THE FIRST STEP 
      IF(ISTEP.EQ.1) THEN
        DO ITWMOD= 1,NTWMOD(IPH)

          IU1=100+10*IPH+ITWMOD
          OPEN(IU1,FILE='TW_SCH_STATS_PH'//CHAR(48+IPH)//'_MODE'
     #        //CHAR(48+ITWMOD)//'.OUT', STATUS='UNKNOWN')
          WRITE(IU1,'(A,I2)')
     #    '   SCHMID FACTOR DISTRIBUTION FOR TWIN MODE:',ITWMOD
     	  WRITE(IU1,'(3A,21F8.3)')
     #    'ISTEP',' STRAIN','  #TWGR', ((0.025+0.05*(J-1)),J=-5,15)

          IU2=200+10*IPH+ITWMOD
          OPEN(IU2,FILE='TW_VAR_STATS_PH'//CHAR(48+IPH)//'_MODE'
     #        //CHAR(48+ITWMOD)//'.OUT', STATUS='UNKNOWN')
     	  WRITE(IU2,'(A,I2)') 
     #    '   VARIANT FREQ & VOL FRACT DISTRIB FOR TWIN MODE:',ITWMOD 
     	  WRITE(IU2,'(3A,12I8)') 
     #    'ISTEP',' STRAIN','  #TWGR',((J),J=1,6),((J),J=1,6)

        ENDDO
      ENDIF

        IPHEL=IPH-IPHBOT+1
        DO ITWMOD= 1,NTWMOD(IPHEL)

          IU1=100+10*IPHEL+ITWMOD
          IU2=200+10*IPHEL+ITWMOD
	   
          NGR_STAT = SUM(SCH_STAT(:,ITWMOD,IPHEL))
          NGR_STATX=NGR_STAT
          IF(NGR_STAT.EQ.0) NGR_STATX=1      ! AVOID DIVIDING BY ZERO WHEN NO VARIANT IS ACTIVE

          WRITE(IU1, '(I4,F8.3,I7,30F8.4)')
     #             ISTEP,EPSACU,NGR_STAT,
     #	  ((SCH_STAT(J,ITWMOD,IPHEL)/NGR_STATX),J=-5,15)

          NGR_STAT = SUM(VAR_STAT(:,ITWMOD,IPHEL))
          NGR_STATX=NGR_STAT
          IF(NGR_STAT.EQ.0) NGR_STATX=1      ! AVOID DIVIDING BY ZERO WHEN NO VARIANT IS ACTIVE

          WRITE(IU2, '(I4,F8.3,I7,10F8.4,10F8.4)')
     #             ISTEP,EPSACU,NGR_STAT,
     #	  ((VAR_STAT(J,ITWMOD,IPHEL)/NGR_STATX),J=1,6),
     #	  ((VAR_VF(J,ITWMOD,IPHEL)),J=1,6)
    	   
        ENDDO	  ! END OF DO ITWMOD

      RETURN
      END
