      PROGRAM READ_RES
C----------------------------------------------------------------------
C-
C-   Package:
C-
C-   Program:            READ_RES
C-
C-   References:         None
C-
C-   Libraries Required: None
C-
C-   Description:
C-
C-
C----------------------------------------------------------------------
C-   Created:   30-NOV-1997         M.V. Yurkov
C----------------------------------------------------------------------
C-   Comments:
C-
C-
C-
C-
C-
C-
C-
C-
C-
C-
C-
C----------------------------------------------------------------------
C-
* Global sections
C__________________
C-
      IMPLICIT REAL*8 ( A-H, O-Z )
      PARAMETER ( NUGDIM = 10000 )
C-
	COMMON /CIDATA /
     &  IFBB,         ! Start No. of simulation run
     &  MODC,         ! Number of simulation runs
     &  MODPRF,       ! Mode for transverse profile of the electron beam (1-3)
     &  NREC,         ! Number of output records
     &  MODETAP,      ! Mode of undulator tapering
     &  NAZ,       ! Azimuthal harmonics: -NAZ ... NAZ
     &  NHARM,        ! Number of odd frequency harmonics, nharm = 1, 2, 3 (k=2*nharm-1)
     &  NRING,        ! Number of radial nodes within electron beam
     &  NNEAR,        ! Number of radial nodes outside electron beam
     &  NPSI,         ! Number of divisions in phase psi
     &  NPSIP,        ! Number of macroparticles for beamlet simulations
     &  MODESS,       ! Flag for steady-state/time-dependent mode (1 - steady state, 0 - time-dependent)
     &  MSEED,        ! Mode of initial longitudinal distribution of particles: smooth / shot noise (0 - smooth /1 - shot noise)
     &  MODEPSI,      ! Mode of smooth start in phase (0 - even, 1 - Sobol)
     &  MODEACT,      ! Mode of shot noise simulations ( -1 - McNeil, 0 - Fawley, 1 - actual number of particles)
     &  KAZEXT,       ! Azimuthal harmonic of external seed
     &  KEXT,         ! Frequency harmonic of external seed = 2*khext-1
     &  MODESRLOSS,   ! Flag for calculation of SR losses
     &  MODEQFL,      ! Flag for calculation of quantum diffusion
     &  MODEVAR,      ! Flag to read tapering pattern from external file
     &  MASTRA,       ! Flag to read bunch slice properties from external file
     &  MODEUNDSTR,   ! Flag to read undulator structure from external file
     &  MODEPHASEJUMP, ! Flag to read phase shifts from external file
     &  MODETRF,      ! Flag to read radiation filter properties
     &  MODER56,      ! Flag for exernal file with R56 FAST.R56
     &  MODEWFRESIST,   	!  Flag to read wakefield profile (resistive) FAST.WAKERES
     &  MODEWFROUGHNESS,	! 	Flag to read wakefield profile (roughness) FAST.WAKEROUGH
     &  MODEWFGEOMET,	  !  Flag to read wakefield profile (geometrical) FAST.WAKEGEOM
     &  MODEWFLSC,    ! Flag to read wakefield profile (LSC) FAST.WAKELSC
     &  IRESTORE,     ! Flag for recovering simulation run      30
     &  NPJUMP, NPR56, NCF
	COMMON /CRDATA /
     &  RLAM,         ! Radiation wavelength, cm
     &  RLAMU,        ! Undulator period, cm
     &  EBEAM,        ! Electron beam energy, MeV
     &  RIBEAM,       ! Peak beam current, A
     &  SIGMAZ,       ! rms length of the electron beam, cm
     &  BLENGTH,      ! length of the electron beam, cm
     &  RBEAM,        ! rms transverse size of electron beam, cm
     &  PPR,           ! Parameter profile for parabolic and gaussian beam profile
     &  BROUND,       ! Diffraction parameter
     &  SCP,          ! Space charge parameter
     &  ESP,          ! Energy spread parameter
     &  BGAM,         ! Parameter fo betatron oscillations
     &  RHO,          ! FEL efficiency parameter (3D)
     &  CDELIN,       ! Detuning parameter
     &  CHIRP,        ! Energy chirp parameter
     &  ZFIN,         ! Undulator length (reduced)
     &  STEPZ,        ! Integration step (reduced)
     &  STEPZMIN,     ! Minimum integration step = 4 PI RHO
     &  BOUT,         ! Start of output
     &  STOUT,        ! Step of output
     &  ZBVAR,        ! Beginning of tapering
     &  AVAR1,        ! Linear term of tapering
     &  AVAR2,        ! Quadratic term of tapering
     &  RNEAR,        ! Step of radial mesh outside electron beam
     &  PEXT,         ! External radiation power (normalized)
     &  ZPER,         ! Position of the laser beam waist
     &  W0,           ! Spot size of the laser beam
     &  A1EXT,        ! Initial modulation of electron beam (prescribed, for test only)
     &  TSR,          ! Intervals for restore point (hour)
     &  RKUNDEXT     ! Undulator parameter for AJJ calculations (used for tests) 30
C
	CHARACTER*256 NAME
	CHARACTER*128 NAMGIN
	CHARACTER*128 NAMG
	CHARACTER*32 NAMT
C
C **** Input arrays
C
	REAL*8, DIMENSION(:) :: RAZB, PROFIL
	REAL*8, DIMENSION(:,:) :: PRAD
	REAL*8, DIMENSION(:) :: PAMPR
	COMPLEX*16, DIMENSION(:,:,:) :: FIELDT
	ALLOCATABLE RAZB, PROFIL, PRAD, PAMPR, FIELDT
      COMPLEX*16, DIMENSION(:,:,:) ::  A1
      ALLOCATABLE A1
      REAL*8, DIMENSION(:) ::  A1AVER2
      ALLOCATABLE A1AVER2
C
	REAL*8, DIMENSION(:) :: EXY
	ALLOCATABLE EXY
C
      REAL*8 ZPHASE(NUGDIM), PHASESHIFT(NUGDIM)
      Integer NZSHIFT(NUGDIM)
	INTEGER NUG(NUGDIM)          !  0/1 - gap/undulator
	REAL*8 ANGU(NUGDIM)
	REAL*8 CDELZ(NUGDIM)
	INTEGER NOUT(NUGDIM)
	REAL*8 ZOUT(NUGDIM)
      REAL*8 ZR56(NUGDIM), DPSIDP(NUGDIM)
      INTEGER NZR56(NUGDIM)
C
	REAL*8 TRF(10,NUGDIM), TRFTEM(NUGDIM), ZTRF(NUGDIM)
	INTEGER NZTRF(NUGDIM)
C
      INTEGER KAZINT(100), KHARMINT (100)
C
      REAL*8, DIMENSION(:,:,:) :: PNZXY
      REAL*8, DIMENSION(:,:,:) :: PNZXYAV
 	ALLOCATABLE PNZXY
 	ALLOCATABLE PNZXYAV
	REAL*8, DIMENSION(:) :: XQ, YQ
	REAL*8, DIMENSION(:,:) :: RQ, FIQ
	ALLOCATABLE XQ, YQ, RQ, FIQ
C-
* Local variables
C__________________
C-
      COMPLEX*16 ICOM
      COMPLEX*16, DIMENSION(:,:,:) :: CFIELD
	ALLOCATABLE CFIELD
C
      COMPLEX*16 CPRAD
C
C-
* External functions and subroutines
C_____________________________________
C-
	integer nprc, nfork, OMP_GET_MAX_THREADS, OMP_GET_THREAD_NUM
	integer OMP_GET_num_THREADS
	external OMP_GET_num_THREADS
	external OMP_GET_MAX_THREADS, OMP_GET_THREAD_NUM
C-
* Equivalence
C______________
C-
C-
* Data
C_______
C-
      DATA NTRES / 4 /
      DATA NNRES / 1 /
      DATA NFRES / 2 /
C
      DATA PI    / 3.14159265358979322702D0 /
      DATA PI2   / 6.28318530717958645404D0 /
      DATA CLIGHT / 3.D10 /
C
      DATA ICOM   / (0.D0,1.D0) /
C-
* Executable statements
C________________________
C-
	CALL CPU_TIME(TSC)
      OPEN(UNIT=1,FILE='FAST2XY_2013.DAT',STATUS='OLD')
      READ(1,*) TRD1
      READ(1,*) TRD2
      READ(1,*) NXY
      READ(1,*) NSKPOINT
      READ(1,*) IFB
      READ(1,*) NZC
	READ(1,*) NAMGIN
      CLOSE(UNIT=1)
      IF ( NSKPOINT .LT. 1 ) NSKPOINT = 1
C
	NAMT = 'T'
	NAMG = NAMGIN
	CALL SETNAME2013(NAME,NAMG,NAMT,IFB,NZC)
	WRITE(6,*) 'Searching file ', TRIM(NAME)
      OPEN(UNIT=NTRES,FORM='UNFORMATTED',
     &  ACCESS='SEQUENTIAL',STATUS='OLD',ERR=20)
	GO TO 10
20	CONTINUE
	WRITE(6,*) 'No files with requested names'
	STOP
10	CONTINUE
	WRITE(6,*) 'Reading file header of reference file and allocating arrays '
      CALL RDHEADER(NTRES)
C
	MODEOUT = 1
      NSLIP = NINT(STEPZ/(4.D0*PI*RHO))
      RNORMT = PI2*RHO*NSLIP
      IF ( (MODEOUT-1)*(MODEOUT-2) .EQ. 0 ) RNORMT = NSLIP*RLAM/CLIGHT*1.D15
      RNORMR = RBEAM
      IF ( MODPRF .EQ. 3 ) RNORMR = RNORMR*DSQRT(2.D0)
C
      NR = NRING + NNEAR
C
	NAZTEM = NAZ-1
	KRING = NRING+NNEAR
C
	ALLOCATE(RAZB(KRING))
	ALLOCATE(PROFIL(KRING))
	ALLOCATE(PRAD(-NAZTEM:NAZTEM,NHARM))
	ALLOCATE(PAMPR(NHARM))
	READ(NTRES) (KAZINT(I),I=-NAZTEM,NAZTEM)
	READ(NTRES) (KHARMINT(I),I=1,NHARM)
      READ(NTRES) (RAZB(I),I=1,NRING+NNEAR)
      READ(NTRES) (PROFIL(I),I=1,NRING)
      READ(NTRES) RADEVE, RADEVP, NFINPR
	CLOSE(UNIT=NTRES)
C
	ALLOCATE(FIELDT(KRING,-NAZTEM:NAZTEM,NHARM))
      ALLOCATE(A1(NRING,-NAZTEM:NAZTEM,NHARM))
      ALLOCATE(A1AVER2(NHARM))
C
      ALLOCATE(CFIELD(-NXY:NXY,-NXY:NXY,NHARM))
      ALLOCATE(PNZXY(-NXY:NXY,-NXY:NXY,NHARM))
      ALLOCATE(PNZXYAV(-NXY:NXY,-NXY:NXY,NHARM))
C
      ALLOCATE(XQ(-NXY:NXY))
      ALLOCATE(YQ(-NXY:NXY))
      ALLOCATE(RQ(-NXY:NXY,-NXY:NXY))
      ALLOCATE(FIQ(-NXY:NXY,-NXY:NXY))
C
	ALLOCATE(EXY(NHARM))
C
      NR = NRING + NNEAR
	DR1 = RAZB(2) - RAZB(1)
	DR2 = RAZB(NRING+1) - RAZB(NRING)
      XYMAX = RAZB(NR-1)/DSQRT(2.D0)
      write(6,*) 'XYMAX: ', XYMAX
      DXY = XYMAX/NXY
      DO IX = -NXY, NXY
         XQ(IX) = DXY*IX
         YQ(IX) = DXY*IX
      ENDDO
      DO IY = -NXY, NXY
      DO IX = -NXY, NXY
         RQ(IX,IY) = SQRT(XQ(IX)**2+YQ(IY)**2)
         FIQ(IX,IY) = 0.D0
      IF ( RQ(IX,IY) .NE. 0.D0 ) FIQ(IX,IY) = ATAN2(YQ(IY),XQ(IX))
      ENDDO
      ENDDO
	WRITE(6,*) 'Input data:'
      WRITE(6,*) 'Time step of input data  (fs) = ', RNORMT
      WRITE(6,*) 'Mesh step (1) of input data (cm) = ', DR1*RNORMR
      WRITE(6,*) 'Mesh step (2) of input data (cm) = ', DR2*RNORMR
	WRITE(6,*) 'Output data:'
      WRITE(6,*) 'Time step of output data (fs) = ', RNORMT*NSKPOINT, '  Ratio: ', NSKPOINT
      WRITE(6,*) 'Number of nodes: ', NXY*2+1
      WRITE(6,*) 'Mesh step (cm) = ', DXY*RNORMR, '  Ratio: ', DXY/DR1
      WRITE(6,*)  'XYMAX, XYMAX*RNORMR, DXY*RNORMR: ', XYMAX, XYMAX*RNORMR, DXY*RNORMR
C
        c = clight
        ri = ribeam*3.d9
        ria = 17043*3.d9
        gam = ebeam/0.511
        step = dxy*rnormr
        rk1 = 8*pi*ri*ria*rho**2*gam/rlam/rlamu/c*step**2*1d-7
        write(6,*) rk1
C
	NAMT = 'T'
	NAMG = NAMGIN
	CALL SETNAME2013(NAME,NAMG,NAMT,IFB,NZC)
      OPEN(UNIT=NTRES,FORM='UNFORMATTED',
     &  ACCESS='SEQUENTIAL',STATUS='OLD',ERR=30)
	GO TO 40
30	CONTINUE
	WRITE(6,*) 'No file ', TRIM(NAME)
	STOP
40	CONTINUE
	WRITE(6,*) 'Input file: ', TRIM(NAME)
	DO KH = 1 , NHARM
	NH = 2*KH-1
	WRITE(NAMT,'(A3,I1,A1)') 'FXY', NH, '_'
	NAMG = NAMGIN
	CALL SETNAME2013(NAME,NAMG,NAMT,IFB,NZC)
      OPEN(UNIT=10+KH,FILE=NAME,ACCESS='SEQUENTIAL',STATUS='UNKNOWN')
	WRITE(6,*) 'Output file: ', TRIM(NAME), '  Nh = ', NH
C      WRITE(10+KH,'(3(G13.6,1X),I5)') RLAM*NH, RNORMT*NSKPOINT*1.D-15, DXY*RNORMR, NXY*2+1
      WRITE(10+KH,'(3(G13.6,1X),I5,(1XG13.6))') RLAM*NH, RNORMT*NSKPOINT*1.D-15, DXY*RNORMR, NXY*2+1, RK1
      WRITE(10+KH,'(2(G13.6,1X))') ICOM
	ENDDO
	NAMT = 'E1'
	NAMG = NAMGIN
	CALL SETNAME2013(NAME,NAMG,NAMT,IFB,NZC)
      OPEN(UNIT=20,FILE=NAME,ACCESS='SEQUENTIAL',STATUS='UNKNOWN')
	WRITE(6,*) 'Output file: ', TRIM(NAME)
      CALL RDHEADER(NTRES)
	READ(NTRES) (KAZINT(I),I=-NAZTEM,NAZTEM)
	READ(NTRES) (KHARMINT(I),I=1,NHARM)
      READ(NTRES) (RAZB(I),I=1,NRING+NNEAR)
      READ(NTRES) (PROFIL(I),I=1,NRING)
      READ(NTRES) RADEVE, RADEVP, NFINPR
	READ(NTRES) (NUG(NRKST),ANGU(NRKST),CDELZ(NRKST),NRKST=1,NFINPR)
	READ(NTRES) (NOUT(I),I=1,NREC)
	IF ( NPJUMP .GE. 1 ) READ(NTRES)(nzshift(i), phaseshift(i),I=1,NPJUMP)
	IF ( NPR56 .GE. 1 ) READ(NTRES)(nzR56(i),DPSIDP(i),I=1,NPR56)
	IF ( NCF .GE. 1 ) THEN
	DO NF = 1 , NCF
	READ(NTRES) ZTRF(NF), NZTRF(NF),(TRF(KH,NF),KH=1,NHARM)
	ENDDO
	ENDIF
C
	PNZXYAV = 0.
C
	NPC = 0
	NSLOUT = 0
      DO NB = 1 , 1000000
C
      DO IJK = 1 , NSKPOINT
      READ(NTRES,END=70,ERR=70) PROFLONG, ENSHIFT, ESPADD, BRTEM
      READ(NTRES,END=70,ERR=70) ESLAV, SIGMAESL
      READ(NTRES,END=70,ERR=70) GR, G
	DO KH = 1 , NHARM
      READ(NTRES,END=70,ERR=70) PAMPR(KH), (PRAD(KAZ,KH),KAZ=-NAZTEM,NAZTEM)
      DO KAZ = -NAZTEM, NAZTEM
      READ(NTRES,END=70,ERR=70) (FIELDT(I,KAZ,KH),I=1,KRING)
      ENDDO  ! KAZ
      READ(NTRES,END=70,ERR=70) A1AVER2(KH)
      DO KAZ = -NAZTEM, NAZTEM
      READ(NTRES,END=70,ERR=70) (A1(I,KAZ,KH),I=1,NRING)
      ENDDO    ! KAZ
      ENDDO    ! KH
	NPC = NPC + 1
      ENDDO	! IJK
C
      TRD = NPC*RNORMT
      IF ( TRD .LT. TRD1 ) GO TO 80
      IF ( TRD .GT. TRD2 ) GO TO 70
C
      PNZXY = 0.D0
	EXY(1:NHARM) = 0.
	DO KH = 1 , NHARM
	NH = KH*2-1
      DO IY = -NXY,NXY
      DO IX = -NXY,NXY
      DO N = 1 , NR-1
      IF ( RAZB(N) .GT. RQ(IX,IY) ) GO TO 201
      ENDDO
 201  CONTINUE
	N = N-1
	IF ( N .EQ. 0 ) N = 1
	R1 = RAZB(N)
	R2 = RAZB(N+1)
	DR = R2-R1
      CFIELD(IX,IY,KH) = (0.D0,0.D0)
      DO KAZ = -NAZTEM, NAZTEM
      CPRAD = FIELDT(N,KAZ,KH) + (FIELDT(N+1,KAZ,KH)-FIELDT(N,KAZ,KH))/DR*(RQ(IX,IY)-R1)
      CFIELD(IX,IY,KH) = CFIELD(IX,IY,KH) + CPRAD*CDEXP(ICOM*KAZ*FIQ(IX,IY))
      ENDDO
      PNZXY(IX,IY,KH) = CDABS(CFIELD(IX,IY,KH))**2
      PNZXYAV(IX,IY,KH) = PNZXYAV(IX,IY,KH) + CDABS(CFIELD(IX,IY,KH))**2
	EXY(KH) = EXY(KH) + PNZXY(IX,IY,KH)
      ENDDO
      ENDDO
      ENDDO
C
	WRITE(20,'(100(G13.6,1X))') TRD, (EXY(KH)*rk1,KH=1,NHARM)
C
	DO KH = 1 , NHARM
      DO IX = -NXY, NXY
      WRITE(10+KH,'(10000(G13.6,1X))') (CFIELD(IX,IY,KH),IY=-NXY,NXY)
      ENDDO
      ENDDO
	NSLOUT = NSLOUT + 1
80	CONTINUE
      ENDDO
70    CONTINUE
C
      CLOSE(UNIT=NTRES)
      CLOSE(UNIT=20)
	DO KH = 1 , NHARM
      CLOSE(UNIT=10+KH)
      ENDDO
      LENGB = NB-1
	WRITE(6,*) 'Number of slices: ', NSLOUT
	WRITE(6,*) 'Pulse length: ', NSLOUT*NSKPOINT*RNORMT
C
	DO KH = 1 , NHARM
	NH = 2*KH-1
	WRITE(NAMT,'(A3,I1,A1)') 'PXY', NH, '_'
	NAMG = NAMGIN
	CALL SETNAME2013(NAME,NAMG,NAMT,IFB,NZC)
      OPEN(UNIT=10+KH,FILE=NAME,ACCESS='SEQUENTIAL',STATUS='UNKNOWN')
	WRITE(6,*) 'Output file: ', TRIM(NAME), '  Nh = ', NH
      DO IX = -NXY, NXY
      WRITE(10+KH,'(1000(G13.6,1X))') (PNZXY(IX,IY,KH),IY=-NXY,NXY) ! (PNZXY(IX,IY,KH),IY=-NXY,NXY)
      ENDDO
      CLOSE(UNIT=10+KH)
      ENDDO
C
 50   CONTINUE
	CALL CPU_TIME(TEC)
	TCPU = TEC - TSC
	WRITE(6,*) 'CPU time: ', TCPU
      END
C-
C-
C-
	SUBROUTINE SETNAME2013(NAME,NAMG,NAMT,IFB,NZ)
	IMPLICIT REAL*8 ( A-H , O-Z )
	CHARACTER*256 NAME
	CHARACTER*128 NAMG
	CHARACTER*32 NAMT
	CHARACTER*4 CIFB
	CHARACTER*3 CNZ
	WRITE(CIFB,'(I4.4)'), IFB
	WRITE(CNZ,'(I3.3)'), NZ
	NAME = ''//TRIM(NAMG)//TRIM(NAMT)//CIFB//CNZ//'.RES'
	RETURN
	END
C-
C-
C-
      SUBROUTINE RDHEADER(NTRES)
	IMPLICIT REAL*8 ( A-H , O-Z )
	COMMON /CIDATA / ID(33)
	COMMON /CRDATA / RD(30)
	READ(NTRES) (ID(I),I=1,33)
	READ(NTRES) (RD(I),I=1,30)
	RETURN
	END
