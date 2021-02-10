PROGRAM DENSITYTOOL
!###                            ###!
!###       Lucas Lodeiro        ###!
!###  Departamento de Química   ###!
!### Lucas.Lodeiro@ug.uchile.cl ###!
!###   Universidad de Chile     ###!
!###                            ###!
!###        Tomáš Rauch         ###!
!###  Institut für Festkörper-  ###!
!###     theorie und -Optik     ###!
!###  tomas.rauch@uni-jena.de   ###!
!###     FSU Jena, Germany      ###!
!###                            ###!
!###          Dic/2020          ###!
!###         Version 0.2        ###!
!###                            ###!
!###     READ THE MANUAL :)     ###!
!###                            ###!
IMPLICIT NONE
INTEGER, PARAMETER :: RSP=KIND(1.0),RDP=KIND(1.0D0)
REAL(RSP), PARAMETER :: PI=3.14159265,EE=2.71828183,ANGTOA0=1.88972613
REAL(RSP), ALLOCATABLE, DIMENSION(:) :: W,PARCHGDATAAVG,PARCHGDATAAVG2,PARCHGDATAAVG3,PARCHGDATAAVG4,L1,L2,L3
REAL(RSP), ALLOCATABLE, DIMENSION(:,:) :: EN,M2
REAL(RSP), ALLOCATABLE, DIMENSION(:,:,:) :: ENERGY,PARCHGDATA,PARCHGDATA2,PARCHGDATA3,PARCHGDATA4
REAL(RSP) :: SIGMA,EMIN,EMAX,DE,E,SMEARING,ETHR,VOLUME,ALAT,ETHRCUTOFF,NPOINTS,TRASHR,MAGMOM,ABSMAGMOM,START,FINISH
INTEGER(4) :: NEN,NATOMS,NKP,NBANDS,NG1,NG2,NG3,NGAVG,DIRECTION,SPINCASE,OPEN_ERROR,TRASHI,I,J,K,L,M,N
CHARACTER(132) :: INAME,ONAME,ONAME2,ONAME3,ONAME4,NUM,BAND,KPOINT,DIR,NGS,READ_WRITE,DATE
CHARACTER(4) :: NGA,NGB,NGC
LOGICAL :: FILE_EXIST,DOPCA,DOPSA,DOPARCHGSPIN,DOCHGCARSPIN,DOCHGCARAVG
LOGICAL :: DOLDOS,DOLDOSFULL,DOLSDOS,DOLSDOSFULL,DOLDOSMAG,DOLDOSFULLMAG,DOLSDOSMAG,DOLSDOSFULLMAG

WRITE(*,*) 'Welcome to DensityTool v0.2 !' ; WRITE(*,*) 'Created by Lucas Lodeiro and Tomáš Rauch.' ; WRITE(*,*)
WRITE(*,*) 'About questions, suggestions or errors within the program or the manual:'
WRITE(*,*) 'lucas.lodeiro@ug.uchile.cl' ; WRITE(*,*) 'tomas.rauch@uni-jena.de' ; WRITE(*,*)
WRITE(*,*) 'The FORTRAN code, input parameter file, manual and application examples can be'
WRITE(*,*) 'found at:  https://github.com/llodeiro/DensityTool' ; WRITE(*,*) ; WRITE(*,*)
WRITE(*,*) 'Enter SIGMA value in eV, for gaussian broadening (0.02 Ry = 0.272113961 eV):'
READ(*,*) SIGMA ; WRITE(*,*)
WRITE(*,*) 'Enter EMIN value in eV, minimum energy grid value to computate local functions' ; WRITE(*,*) '(EMAX > EMIN):'
READ(*,*) EMIN ; WRITE(*,*)
WRITE(*,*) 'Enter EMAX value in eV, maximum energy grid value to computate local functions' ; WRITE(*,*) '(EMAX > EMIN):'
READ(*,*) EMAX ; WRITE(*,*)
WRITE(*,*) 'Enter NEN value (Integer), number of energy grid points between the minimum'
WRITE(*,*) 'and maximum energy (10 to 100 per eV):'
READ(*,*) NEN ; WRITE(*,*)
WRITE(*,*) 'Enter ETHR value, energy threshold for smearing contribution cutoff' ; WRITE(*,*) '(0.001 is the standard):'
READ(*,*) ETHR ; WRITE(*,*)
WRITE(*,*) 'Enter SPINCASE value (1,2 or 3), non-magnetic, collinear magnetic and non-' ; WRITE(*,*) 'collinear cases:'
READ(*,*) SPINCASE ; WRITE(*,*)
WRITE(*,*) 'Enter DIRECTION value (1,2 or 3), for r1,r2 or r3 (r_i) averaging direction:'
READ(*,*) DIRECTION ; WRITE(*,*)
WRITE(*,*) 'Enter DOPCA logical value (T or F), to do (or not) the partial charge average'
WRITE(*,*) 'over the plane (r_j,r_k) (PCA):'
READ(*,*) DOPCA ; WRITE(*,*)
WRITE(*,*) 'Enter DOPSA logical value (T or F), to do (or not) the partial spin average'
WRITE(*,*) 'over the plane (r_j,r_k) (PSA):'
READ(*,*) DOPSA ; WRITE(*,*)
WRITE(*,*) 'Enter DOLDOS logical value (T or F), to do (or not) the LDOS(E,r_i), previous'
WRITE(*,*) '"PCA" calculation is mandatory:'
READ(*,*) DOLDOS ; WRITE(*,*)
WRITE(*,*) 'Enter DOLDOSFULL logical value (T or F), to do (or not) the LDOS(E,r),'
WRITE(*,*) 'previous "PCA" calculation is not mandatory:'
READ(*,*) DOLDOSFULL ; WRITE(*,*)
WRITE(*,*) 'Enter DOLSDOS logical value (T or F), to do (or not) the LSDOS(E,r_i),'
WRITE(*,*) 'previous "PSA" calculation is mandatory:'
READ(*,*) DOLSDOS ; WRITE(*,*)
WRITE(*,*) 'Enter DOLSDOSFULL logical value (T or F), to do (or not) the LSDOS(E,r),'
WRITE(*,*) 'previous "PSA" calculation is not mandatory:'
READ(*,*) DOLSDOSFULL ; WRITE(*,*)
WRITE(*,*) 'Enter DOPARCHGSPIN logical value (T or F), to do (or not) the Alpha and Beta'
WRITE(*,*) 'partial densities and "PCA" for each spin density:'
READ(*,*) DOPARCHGSPIN ; WRITE(*,*)
WRITE(*,*) 'Enter DOCHGCARSPIN logical value (T or F), to do (or not) the Alpha and Beta'
WRITE(*,*) 'total densities and "PCA" for each spin density:'
READ(*,*) DOCHGCARSPIN ; WRITE(*,*)
WRITE(*,*) 'Enter DOCHGCARAVG logical value (T or F), to do (or not) charge and spin'
WRITE(*,*) 'magnetization average:'
READ(*,*) DOCHGCARAVG ; WRITE(*,*)
WRITE(*,*) 'Enter DOLDOSMAG logical value (T or F), to do (or not) the LDOS(E,r_i),'
WRITE(*,*) 'previous "PARCHGSPIN" calculation is mandatory:'
READ(*,*) DOLDOSMAG ; WRITE(*,*)
WRITE(*,*) 'Enter DOLDOSFULLMAG logical value (T or F), to do (or not) the LDOS(E,r),'
WRITE(*,*) 'previous "PARCHGSPIN" calculation is mandatory:'
READ(*,*) DOLDOSFULLMAG ; WRITE(*,*)
WRITE(*,*) 'Enter DOLSDOSMAG logical value (T or F), to do (or not) the LSDOS(E,r_i),'
WRITE(*,*) 'previous "PARCHGSPIN" calculation is mandatory:'
READ(*,*) DOLSDOSMAG ; WRITE(*,*)
WRITE(*,*) 'Enter DOLSDOSFULLMAG logical value (T or F), to do (or not) the LSDOS(E,r),'
WRITE(*,*) 'previous "PARCHGSPIN" calculation is mandatory:'
READ(*,*) DOLSDOSFULLMAG ; WRITE(*,*)
WRITE(*,*) 
WRITE(*,*) 'The selected variable values are:'
WRITE(*,*) 'SIGMA          = ',SIGMA ; WRITE(*,*) 'EMIN           = ',EMIN ; WRITE(*,*) 'EMAX           = ',EMAX
WRITE(*,*) 'NEN            = ',NEN ; WRITE(*,*) 'ETHR           = ',ETHR ; WRITE(*,*) 'SPINCASE       = ',SPINCASE
WRITE(*,*) 'DIRECTION      = ',DIRECTION ; WRITE(*,*) 'DOPCA          = ',DOPCA ; WRITE(*,*) 'DOPSA          = ',DOPSA
WRITE(*,*) 'DOLDOS         = ',DOLDOS ; WRITE(*,*) 'DOLDOSFULL     = ',DOLDOSFULL ; WRITE(*,*) 'DOLSDOS        = ',DOLSDOS
WRITE(*,*) 'DOLSDOSFULL    = ',DOLSDOSFULL ; WRITE(*,*) 'DOPARCHGSPIN   = ',DOPARCHGSPIN
WRITE(*,*) 'DOCHGCARSPIN   = ',DOCHGCARSPIN ; WRITE(*,*) 'DOCHGCARAVG    = ',DOCHGCARAVG
WRITE(*,*) 'DOLDOSMAG      = ',DOLDOSMAG ; WRITE(*,*) 'DOLDOSFULLMAG  = ',DOLDOSFULLMAG
WRITE(*,*) 'DOLSDOSMAG     = ',DOLSDOSMAG ; WRITE(*,*) 'DOLSDOSFULLMAG = ',DOLSDOSFULLMAG ; WRITE(*,*) ; WRITE(*,*) 

CALL CPU_TIME(START) ; CALL FDATE(DATE) 

OPEN(12,FILE='EIGENVAL',FORM='FORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=OPEN_ERROR)
IF (OPEN_ERROR > 0) STOP 'An error occurred while opening the file EIGENVAL'
WRITE(*,*) 'Reading EIGENVAL'
READ(12,*) NATOMS
DO I=1,4
  READ(12,*)
ENDDO
READ(12,*) TRASHI,NKP,NBANDS ; READ(12,*)
ALLOCATE(EN(NKP,NBANDS),M2(NKP,NBANDS),ENERGY(NKP,NBANDS,2),W(NKP))
DO J=1,NKP
  READ(12,*) TRASHR,TRASHR,TRASHR,W(J)
  DO I=1,NBANDS
    READ(12,*) TRASHI,ENERGY(J,I,1),ENERGY(J,I,2)
  ENDDO
  IF (J.LT.NKP) THEN
    READ(12,*)
  ENDIF
ENDDO
CLOSE(12)
IF (SPINCASE .EQ. 2) THEN
  WRITE(*,*) 'SPINCASE = 2. If you are computing LDOS or LSDOS, remember that an average'
  WRITE(*,*) 'energy between alpha and beta orbitals is used. Check the standard deviation'
  WRITE(*,*) '(between alpha/beta energy and average value), and the maximum error in energy'
  WRITE(*,*) 'due to average values (generally in the higher energy unoccupied orbitals). If'
  WRITE(*,*) 'these parameters are high for you, check the LDOS and LSDOS subroutines using'
  WRITE(*,*) 'alpha and beta data (orbitals and energies): DOLDOSMAG, DOLDOSFULLMAG,'
  WRITE(*,*) 'DOLSDOSMAG and DOLSDOSFULLMAG.'
  EN(:,:)=0.5*ENERGY(:,:,1)+0.5*ENERGY(:,:,2) ; M2(:,:)=ABS(0.5*ENERGY(:,:,1)-0.5*ENERGY(:,:,2))
  WRITE(*,*) 'Maximum Error      = ', MAXVAL(M2), ' eV'
  TRASHR=0.0
  DO J=1,NKP
    DO I=1,NBANDS
      TRASHR=TRASHR+(2/REAL((NKP*NBANDS),RSP))*(ENERGY(J,I,1)-EN(J,I))**2.0
    ENDDO
  ENDDO
  WRITE(*,*) 'Standard Deviation = ', SQRT(TRASHR), ' eV'
ELSE 
  EN(:,:)=ENERGY(:,:,1)
ENDIF
WRITE(*,*) ; WRITE(*,*) 

ALLOCATE(L1(3),L2(3),L3(3))
OPEN(12,FILE='CHGCAR',FORM='FORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=OPEN_ERROR)
IF (OPEN_ERROR > 0) STOP 'An error occurred while opening the file CHGCAR'
WRITE(*,*) 'Reading CHGCAR'
READ(12,*) ; READ(12,*) ALAT ; READ(12,*) L1 ; READ(12,*) L2 ; READ(12,*) L3
DO I=1,NATOMS+4
  READ(12,*)
ENDDO
READ(12,*) NG1,NG2,NG3
CLOSE(12)
VOLUME=(ALAT**3)*(ANGTOA0**3)*ABS(L3(1)*(L1(2)*L2(3)-L1(3)*L2(2))+L3(2)*(L1(3)*L2(1)-L1(1)*L2(3))+L3(3)*(L1(1)*L2(2)-L1(2)*L2(1)))
WRITE(*,*) 'VOLUME = ', VOLUME, 'bohr³'
DEALLOCATE(L1,L2,L3)
WRITE(NGA,'(I4)') NG1 ; WRITE(NGB,'(I4)') NG2 ; WRITE(NGC,'(I4)') NG3
NGS=' '//TRIM(NGA)//' '//TRIM(NGB)//' '//TRIM(NGC)
NPOINTS=REAL(NG1*NG2*NG3,RSP)
WRITE(*,*) ; WRITE(*,*) 

IF (DOPCA) THEN
WRITE(*,*) 'PCA activated.'
ALLOCATE(PARCHGDATA(NG1,NG2,NG3))
DO K=1,NKP
  DO J=1,NBANDS
    WRITE(BAND,'(I4.4)') J ; WRITE(KPOINT,'(I4.4)') K
    INAME = 'PARCHG.'//TRIM(BAND)//'.'//TRIM(KPOINT)
    INQUIRE(FILE=TRIM(INAME), EXIST=FILE_EXIST)
    IF (FILE_EXIST) THEN
      WRITE(*,*) 'Reading: ',TRIM(INAME)
      OPEN(12,FILE=TRIM(INAME),FORM='FORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=OPEN_ERROR)
      IF (OPEN_ERROR > 0) STOP 'An error occurred while opening the file '//TRIM(INAME)
      DO
        READ(12,'(A132)') READ_WRITE ; IF ((TRIM(READ_WRITE)).EQ.(TRIM(NGS))) EXIT
      ENDDO
      READ(12,*) PARCHGDATA
      CLOSE(12)
      SELECT CASE (DIRECTION)
        CASE(1)
        ALLOCATE(PARCHGDATAAVG(NG1))
        PARCHGDATAAVG=0.0
        DO L=1,NG2
          DO M=1,NG3
            PARCHGDATAAVG(:) = PARCHGDATAAVG(:) + PARCHGDATA(:,L,M)/(NG2*NG3)
          ENDDO
        ENDDO
        ONAME = TRIM(INAME)//'.R1'
        CASE(2)
        ALLOCATE(PARCHGDATAAVG(NG2))
        PARCHGDATAAVG=0.0
        DO L=1,NG1
          DO M=1,NG3
            PARCHGDATAAVG(:) = PARCHGDATAAVG(:) + PARCHGDATA(L,:,M)/(NG1*NG3)
          ENDDO
        ENDDO
        ONAME = TRIM(INAME)//'.R2'
        CASE(3)
        ALLOCATE(PARCHGDATAAVG(NG3))
        PARCHGDATAAVG=0.0
        DO L=1,NG1
          DO M=1,NG2
            PARCHGDATAAVG(:) = PARCHGDATAAVG(:) + PARCHGDATA(L,M,:)/(NG1*NG2)
          ENDDO
        ENDDO
        ONAME = TRIM(INAME)//'.R3'
      END SELECT
      OPEN(12,FILE=TRIM(ONAME),FORM='FORMATTED',STATUS='NEW',ACTION='WRITE',IOSTAT=OPEN_ERROR)
      IF (OPEN_ERROR > 0) STOP 'An error occurred while creating the file '//TRIM(ONAME)//', the file already exists'
      WRITE(12,*) PARCHGDATAAVG
      CLOSE(12) 
      DEALLOCATE(PARCHGDATAAVG)
    ENDIF
  ENDDO
ENDDO
DEALLOCATE(PARCHGDATA) ; WRITE(*,*) ; WRITE(*,*) 
ENDIF

IF ((DOPSA).AND.(SPINCASE .EQ. 1)) THEN
WRITE(*,*) 'PSA is not activated. PSA in collinear non-magnetic case has no sense... The'
WRITE(*,*) 'spin density is zero everywhere by definition. :P' ; WRITE(*,*) ; WRITE(*,*) 

ELSEIF ((DOPSA).AND.(SPINCASE .EQ. 2)) THEN
WRITE(*,*) 'PSA activated, doing PSA in collinear magnetic case.'
ALLOCATE(PARCHGDATA(NG1,NG2,NG3))
DO K=1,NKP
  DO J=1,NBANDS
    WRITE(BAND,'(I4.4)') J ; WRITE(KPOINT,'(I4.4)') K
    INAME = 'PARCHG.'//TRIM(BAND)//'.'//TRIM(KPOINT)
    INQUIRE(FILE=TRIM(INAME), EXIST=FILE_EXIST)
    IF (FILE_EXIST) THEN
      WRITE(*,*) 'Reading: ',TRIM(INAME)
      OPEN(12,FILE=TRIM(INAME),FORM='FORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=OPEN_ERROR)
      IF (OPEN_ERROR > 0) STOP 'An error occurred while opening the file '//TRIM(INAME)
      DO
        READ(12,'(A132)') READ_WRITE ; IF ((TRIM(READ_WRITE)).EQ.(TRIM(NGS))) EXIT
      ENDDO
      DO
        READ(12,'(A132)') READ_WRITE ; IF ((TRIM(READ_WRITE)).EQ.(TRIM(NGS))) EXIT
      ENDDO
      READ(12,*) PARCHGDATA
      CLOSE(12)
      SELECT CASE (DIRECTION)
        CASE(1)
        ALLOCATE(PARCHGDATAAVG(NG1))
        PARCHGDATAAVG=0.0
        DO L=1,NG2
          DO M=1,NG3
            PARCHGDATAAVG(:) = PARCHGDATAAVG(:) + PARCHGDATA(:,L,M)/(NG2*NG3)
          ENDDO
        ENDDO
        ONAME = TRIM(INAME)//'.SR1'
        CASE(2)
        ALLOCATE(PARCHGDATAAVG(NG2))
        PARCHGDATAAVG=0.0
        DO L=1,NG1
          DO M=1,NG3
            PARCHGDATAAVG(:) = PARCHGDATAAVG(:) + PARCHGDATA(L,:,M)/(NG1*NG3)
          ENDDO
        ENDDO
        ONAME = TRIM(INAME)//'.SR2'
        CASE(3)
        ALLOCATE(PARCHGDATAAVG(NG3))
        PARCHGDATAAVG=0.0
        DO L=1,NG1
          DO M=1,NG2
            PARCHGDATAAVG(:) = PARCHGDATAAVG(:) + PARCHGDATA(L,M,:)/(NG1*NG2)
          ENDDO
        ENDDO
        ONAME = TRIM(INAME)//'.SR3'
      END SELECT
      OPEN(12,FILE=TRIM(ONAME),FORM='FORMATTED',STATUS='NEW',ACTION='WRITE',IOSTAT=OPEN_ERROR)
      IF (OPEN_ERROR > 0) STOP 'An error occurred while creating the file '//TRIM(ONAME)//', the file already exists'
      WRITE(12,*) PARCHGDATAAVG
      CLOSE(12) 
      DEALLOCATE(PARCHGDATAAVG)
    ENDIF
  ENDDO
ENDDO
DEALLOCATE(PARCHGDATA) ; WRITE(*,*) ; WRITE(*,*) 

ELSEIF ((DOPSA).AND.(SPINCASE .EQ. 3)) THEN
WRITE(*,*) 'PSA is not activated. PSA in non-collinear case is not implemented yet, due to'
WRITE(*,*) 'VASP (6.0.8 or latter) does not print the spin densities in PARCHG, only in CHG'
WRITE(*,*) 'and CHGCAR. If VASP developers implement the spin densities in PARCHG in the'
WRITE(*,*) 'future, we can implement it here :) . Ask about it to VASP developers.' ; WRITE(*,*) ; WRITE(*,*) 
ENDIF

IF (DOLDOS) THEN
WRITE(*,*) 'LDOS activated.'
DE=(EMAX-EMIN)/NEN ; ETHRCUTOFF=ETHR/(SIGMA*SQRT(2*PI))
WRITE(*,*) 'Threshold Cutoff=', ETHRCUTOFF ; WRITE(*,*) 'Format: Reading "file name"'
WRITE(*,*) '        "Smearing value", "Eigen-energy", "Energy"'
WRITE(*,*) '        Is smearing value > Threshold Cutoff? YES or NO (not printed)'
SELECT CASE (DIRECTION)
  CASE(1) ; NGAVG=NG1 ; CASE(2); NGAVG=NG2 ; CASE(3) ; NGAVG=NG3
END SELECT
ALLOCATE(PARCHGDATAAVG2(NGAVG))
DO I=1,NEN+1
  E = DE*(I-1)+EMIN
  PARCHGDATAAVG2 = 0.0
  ALLOCATE(PARCHGDATAAVG(NGAVG))
  DO K=1,NKP
    DO J=1,NBANDS
      WRITE(BAND,'(I4.4)') J ; WRITE(KPOINT,'(I4.4)') K ; WRITE(DIR,'(I1.1)') DIRECTION
      INAME = 'PARCHG.'//TRIM(BAND)//'.'//TRIM(KPOINT)//'.R'//TRIM(DIR)
      INQUIRE(FILE=TRIM(INAME), EXIST=FILE_EXIST)
      IF (FILE_EXIST) THEN
        SMEARING = (1.0/SQRT(2.0*PI*SIGMA**2.0))*(EE**(-((E-EN(K,J))**2.0)/(2.0*SIGMA**2.0)))
        WRITE(*,*) 'Reading: ',TRIM(INAME) ; WRITE(*,*) SMEARING, EN(K,J), E
        IF (SMEARING.GT.(ETHRCUTOFF)) THEN
          WRITE(*,*) 'YES'
          OPEN(12,FILE=TRIM(INAME),FORM='FORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=OPEN_ERROR)
          IF (OPEN_ERROR > 0) STOP 'An error occurred while opening the file '//TRIM(INAME)
          READ(12,*) PARCHGDATAAVG
          CLOSE(12)
          PARCHGDATAAVG2 = PARCHGDATAAVG2 + SMEARING*W(K)*PARCHGDATAAVG/VOLUME
        ENDIF
      ENDIF
    ENDDO
  ENDDO
  DEALLOCATE(PARCHGDATAAVG)
  WRITE (NUM,'(I4.4)') I-1 ; WRITE(DIR,'(I1.1)') DIRECTION
  ONAME = 'LDOS.R'//TRIM(DIR)//'.'//TRIM(NUM)//'.dat'
  OPEN(12,FILE=TRIM(ONAME),FORM='FORMATTED',STATUS='NEW',ACTION='WRITE',IOSTAT=OPEN_ERROR)
  IF (OPEN_ERROR > 0) STOP 'An error occurred while creating the file '//TRIM(ONAME)//', the file already exists'
  DO J=1,NGAVG
    WRITE(12,*) REAL(J-1),'	',E,'	',PARCHGDATAAVG2(J)
  ENDDO
  CLOSE(12)
ENDDO
DEALLOCATE(PARCHGDATAAVG2) ; WRITE(*,*) ; WRITE(*,*) 
ENDIF

IF ((DOLSDOS).AND.(SPINCASE .EQ. 1)) THEN
WRITE(*,*) 'LSDOS is not activated. LSDOS in collinear non-magnetic case has no sense...'
WRITE(*,*) 'The spin density is zero everywhere by definition. :P' ; WRITE(*,*) ; WRITE(*,*) 

ELSEIF ((DOLSDOS).AND.(SPINCASE .EQ. 2)) THEN
WRITE(*,*) 'LSDOS activated, computing LSDOS in collinear magnetic case.'
DE=(EMAX-EMIN)/NEN ; ETHRCUTOFF=ETHR/(SIGMA*SQRT(2*PI))
WRITE(*,*) 'Threshold Cutoff=', ETHRCUTOFF ; WRITE(*,*) 'Format: Reading "file name"'
WRITE(*,*) '        "Smearing value", "Eigen-energy", "Energy"'
WRITE(*,*) '        Is smearing value > Threshold Cutoff? YES or NO (not printed)'
SELECT CASE (DIRECTION)
  CASE(1) ; NGAVG=NG1 ; CASE(2) ; NGAVG=NG2 ; CASE(3) ; NGAVG=NG3
END SELECT
ALLOCATE(PARCHGDATAAVG2(NGAVG))
DO I=1,NEN+1
  E = DE*(I-1)+EMIN
  PARCHGDATAAVG2 = 0.0
  ALLOCATE(PARCHGDATAAVG(NGAVG))
  DO K=1,NKP
    DO J=1,NBANDS
      WRITE(BAND,'(I4.4)') J ; WRITE(KPOINT,'(I4.4)') K ; WRITE(DIR,'(I1.1)') DIRECTION
      INAME = 'PARCHG.'//TRIM(BAND)//'.'//TRIM(KPOINT)//'.SR'//TRIM(DIR)
      INQUIRE(FILE=TRIM(INAME), EXIST=FILE_EXIST)
      IF (FILE_EXIST) THEN
        SMEARING = (1.0/SQRT(2.0*PI*SIGMA**2.0))*(EE**(-((E-EN(K,J))**2.0)/(2.0*SIGMA**2.0)))
        WRITE(*,*) 'Reading: ',TRIM(INAME) ; WRITE(*,*) SMEARING, EN(K,J), E
        IF (SMEARING.GT.(ETHRCUTOFF)) THEN
          WRITE(*,*) 'YES'
          OPEN(12,FILE=TRIM(INAME),FORM='FORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=OPEN_ERROR)
          IF (OPEN_ERROR > 0) STOP 'An error occurred while opening the file '//TRIM(INAME)
          READ(12,*) PARCHGDATAAVG
          CLOSE(12)
          PARCHGDATAAVG2 = PARCHGDATAAVG2 + SMEARING*W(K)*PARCHGDATAAVG/VOLUME
        ENDIF
      ENDIF
    ENDDO
  ENDDO
  DEALLOCATE(PARCHGDATAAVG)
  WRITE (NUM,'(I4.4)') I-1 ; WRITE(DIR,'(I1.1)') DIRECTION
  ONAME = 'LSDOS.R'//TRIM(DIR)//'.'//TRIM(NUM)//'.dat'
  OPEN(12,FILE=TRIM(ONAME),FORM='FORMATTED',STATUS='NEW',ACTION='WRITE',IOSTAT=OPEN_ERROR)
  IF (OPEN_ERROR > 0) STOP 'An error occurred while creating the file '//TRIM(ONAME)//', the file already exists'
  DO J=1,NGAVG
    WRITE(12,*) REAL(J-1),'	',E,'	',PARCHGDATAAVG2(J)
  ENDDO
  CLOSE(12)
ENDDO
DEALLOCATE(PARCHGDATAAVG2) ; WRITE(*,*) ; WRITE(*,*) 

ELSEIF ((DOLSDOS).AND.(SPINCASE .EQ. 3)) THEN
WRITE(*,*) 'LSDOS is not activated. LSDOS in non-collinear case is not implemented yet, due'
WRITE(*,*) 'to VASP (6.0.8 or later) does not print the spin densities in PARCHG, only in'
WRITE(*,*) 'CHG and CHGCAR. If VASP developers implement the spin densities in PARCHG in'
WRITE(*,*) 'the future, we can implement it here :) . Ask the VASP developers about it.' ; WRITE(*,*) ; WRITE(*,*) 
ENDIF

IF (DOLDOSFULL) THEN
WRITE(*,*) 'LDOSFULL activated, computing full LDOS'
DE=(EMAX-EMIN)/NEN ; ETHRCUTOFF=ETHR/(SIGMA*SQRT(2*PI))
WRITE(*,*) 'Threshold Cutoff=', ETHRCUTOFF ; WRITE(*,*) 'Format: Reading "file name"'
WRITE(*,*) '        "Smearing value", "Eigen-energy", "Energy"'
WRITE(*,*) '        Is smearing value > Threshold Cutoff? YES or NO (not printed)'
ALLOCATE(PARCHGDATA2(NG1,NG2,NG3))
DO I=1,NEN+1
  E = DE*(I-1)+EMIN
  PARCHGDATA2 = 0.0
  ALLOCATE(PARCHGDATA(NG1,NG2,NG3))
  DO K=1,NKP
    DO J=1,NBANDS
      WRITE(BAND,'(I4.4)') J ; WRITE(KPOINT,'(I4.4)') K
      INAME = 'PARCHG.'//TRIM(BAND)//'.'//TRIM(KPOINT)
      INQUIRE(FILE=TRIM(INAME), EXIST=FILE_EXIST)
      IF (FILE_EXIST) THEN
        SMEARING = (1.0/SQRT(2.0*PI*SIGMA**2.0))*(EE**(-((E-EN(K,J))**2.0)/(2.0*SIGMA**2.0)))
        WRITE(*,*) 'Reading: ',TRIM(INAME) ; WRITE(*,*) SMEARING, EN(K,J), E
        IF (SMEARING.GT.(ETHRCUTOFF)) THEN
          WRITE(*,*) 'YES'
          OPEN(12,FILE=TRIM(INAME),FORM='FORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=OPEN_ERROR)
          IF (OPEN_ERROR > 0) STOP 'An error occurred while opening the file '//TRIM(INAME)
          DO
            READ(12,'(A132)') READ_WRITE ; IF ((TRIM(READ_WRITE)).EQ.(TRIM(NGS))) EXIT
          ENDDO
          READ(12,*) PARCHGDATA
          CLOSE(12)
          PARCHGDATA2 = PARCHGDATA2 + SMEARING*W(K)*PARCHGDATA/VOLUME
        ENDIF
      ENDIF
    ENDDO
  ENDDO
  DEALLOCATE(PARCHGDATA)
  WRITE (NUM,'(I4.4)') I-1
  ONAME = 'LDOS.FULL.'//TRIM(NUM)//'.dat'
  OPEN(12,FILE=TRIM(ONAME),FORM='FORMATTED',STATUS='NEW',ACTION='WRITE',IOSTAT=OPEN_ERROR)
  IF (OPEN_ERROR > 0) STOP 'An error occurred while creating the file '//TRIM(ONAME)//', the file already exists'
  OPEN(13,FILE='CHGCAR',FORM='FORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=OPEN_ERROR)
  IF (OPEN_ERROR > 0) STOP 'An error occurred while opening the file CHGCAR'
  READ(13,'(A132)') READ_WRITE ; WRITE(12,*) TRIM(READ_WRITE),'     ENERGY = ',E
  DO
    READ(13,'(A132)') READ_WRITE ; WRITE(12,*) TRIM(READ_WRITE) ; IF ((TRIM(READ_WRITE)).EQ.(TRIM(NGS))) EXIT
  ENDDO
  WRITE(12,*) PARCHGDATA2
  CLOSE(13) ; CLOSE(12) 
ENDDO
DEALLOCATE(PARCHGDATA2) ; WRITE(*,*) ; WRITE(*,*) 
ENDIF

IF ((DOLSDOSFULL).AND.(SPINCASE .EQ. 1)) THEN
WRITE(*,*) 'LSDOSFULL is not activated. Full LSDOS in collinear non-magnetic case has no'
WRITE(*,*) 'sense... The spin density is zero everywhere by definition. :P' ; WRITE(*,*) ; WRITE(*,*) 

ELSEIF ((DOLSDOSFULL).AND.(SPINCASE .EQ. 2)) THEN
WRITE(*,*) 'LSDOSFULL activated, computing full LSDOS in collinear magnetic case.'
DE=(EMAX-EMIN)/NEN ; ETHRCUTOFF=ETHR/(SIGMA*SQRT(2*PI))
WRITE(*,*) 'Threshold Cutoff=', ETHRCUTOFF ; WRITE(*,*) 'Format: Reading "file name"'
WRITE(*,*) '        "Smearing value", "Eigen-energy", "Energy"'
WRITE(*,*) '        Is smearing value > Threshold Cutoff? YES or NO (not printed)'
ALLOCATE(PARCHGDATA2(NG1,NG2,NG3))
DO I=1,NEN+1
  E = DE*(I-1)+EMIN
  PARCHGDATA2 = 0.0
  ALLOCATE(PARCHGDATA(NG1,NG2,NG3))
  DO K=1,NKP
    DO J=1,NBANDS
      WRITE(BAND,'(I4.4)') J ; WRITE(KPOINT,'(I4.4)') K
      INAME = 'PARCHG.'//TRIM(BAND)//'.'//TRIM(KPOINT)
      INQUIRE(FILE=TRIM(INAME), EXIST=FILE_EXIST)
      IF (FILE_EXIST) THEN
        SMEARING = (1.0/SQRT(2.0*PI*SIGMA**2.0))*(EE**(-((E-EN(K,J))**2.0)/(2.0*SIGMA**2.0)))
        WRITE(*,*) 'Reading: ',TRIM(INAME) ; WRITE(*,*) SMEARING, EN(K,J), E
        IF (SMEARING.GT.(ETHRCUTOFF)) THEN
          WRITE(*,*) 'YES'
          OPEN(12,FILE=TRIM(INAME),FORM='FORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=OPEN_ERROR)
          IF (OPEN_ERROR > 0) STOP 'An error occurred while opening the file '//TRIM(INAME)
          DO
            READ(12,'(A132)') READ_WRITE ; IF ((TRIM(READ_WRITE)).EQ.(TRIM(NGS))) EXIT
          ENDDO
          DO
            READ(12,'(A132)') READ_WRITE ; IF ((TRIM(READ_WRITE)).EQ.(TRIM(NGS))) EXIT
          ENDDO
          READ(12,*) PARCHGDATA
          CLOSE(12)
          PARCHGDATA2 = PARCHGDATA2 + SMEARING*W(K)*PARCHGDATA/VOLUME
        ENDIF
      ENDIF
    ENDDO
  ENDDO
  DEALLOCATE(PARCHGDATA)
  WRITE (NUM,'(I4.4)') I-1
  ONAME = 'LSDOS.FULL.'//TRIM(NUM)//'.dat'
  OPEN(12,FILE=TRIM(ONAME),FORM='FORMATTED',STATUS='NEW',ACTION='WRITE',IOSTAT=OPEN_ERROR)
  IF (OPEN_ERROR > 0) STOP 'An error occurred while creating the file '//TRIM(ONAME)//', the file already exists'
  OPEN(13,FILE='CHGCAR',FORM='FORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=OPEN_ERROR)
  IF (OPEN_ERROR > 0) STOP 'An error occurred while opening the file CHGCAR'
  READ(13,'(A132)') READ_WRITE ; WRITE(12,*) TRIM(READ_WRITE),'     ENERGY = ',E
  DO
    READ(13,'(A132)') READ_WRITE ; WRITE(12,*) TRIM(READ_WRITE) ; IF ((TRIM(READ_WRITE)).EQ.(TRIM(NGS))) EXIT
  ENDDO
  WRITE(12,*) PARCHGDATA2
  CLOSE(13) ; CLOSE(12)
ENDDO
DEALLOCATE(PARCHGDATA2) ; WRITE(*,*) ; WRITE(*,*) 

ELSEIF ((DOLSDOSFULL).AND.(SPINCASE .EQ. 3)) THEN
WRITE(*,*) 'LSDOSFULL is not activated. Full LSDOS in non-collinear case is not implemented'
WRITE(*,*) 'yet, due to VASP (6.0.8 or latter) does not print the spin densities in PARCHG,'
WRITE(*,*) 'only in CHG and CHGCAR. If VASP developers implement the spin densities in'
WRITE(*,*) 'PARCHG in the future, we can implement it here :) . Ask about it to VASP'
WRITE(*,*) 'developers.' ; WRITE(*,*) ; WRITE(*,*) 
ENDIF

IF ((DOPARCHGSPIN).AND.(SPINCASE .EQ. 1)) THEN
WRITE(*,*) 'PARCHGSPIN is not activated. PARCHGSPIN in collinear non-magnetic case has no'
WRITE(*,*) 'sense... The Alpha and Beta densities are equal by definition, and their value'
WRITE(*,*) 'is a half of the partial density... :P' ; WRITE(*,*) ; WRITE(*,*) 

ELSEIF ((DOPARCHGSPIN).AND.(SPINCASE .EQ. 2)) THEN
WRITE(*,*) 'PARCHGSPIN is activated. Computing the Alpha and Beta Partial Densities, and'
WRITE(*,*) 'PCA over each data.'
ALLOCATE(PARCHGDATA(NG1,NG2,NG3),PARCHGDATA2(NG1,NG2,NG3),PARCHGDATA3(NG1,NG2,NG3),PARCHGDATA4(NG1,NG2,NG3))
DO K=1,NKP
  DO J=1,NBANDS
    WRITE(BAND,'(I4.4)') J ; WRITE(KPOINT,'(I4.4)') K
    INAME = 'PARCHG.'//TRIM(BAND)//'.'//TRIM(KPOINT)
    INQUIRE(FILE=TRIM(INAME), EXIST=FILE_EXIST)
    IF (FILE_EXIST) THEN
      WRITE(*,*) 'Reading: ',TRIM(INAME)
      OPEN(12,FILE=TRIM(INAME),FORM='FORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=OPEN_ERROR)
      IF (OPEN_ERROR > 0) STOP 'An error occurred while opening the file '//TRIM(INAME)
      DO
        READ(12,'(A132)') READ_WRITE ; IF ((TRIM(READ_WRITE)).EQ.(TRIM(NGS))) EXIT
      ENDDO
      READ(12,*) PARCHGDATA
      DO
        READ(12,'(A132)') READ_WRITE ; IF ((TRIM(READ_WRITE)).EQ.(TRIM(NGS))) EXIT
      ENDDO
      READ(12,*) PARCHGDATA2
      CLOSE(12)
      PARCHGDATA3=0.5*(PARCHGDATA+PARCHGDATA2) ; PARCHGDATA4=0.5*(PARCHGDATA-PARCHGDATA2)
      ONAME = TRIM(INAME)//'.ALPHA'
      OPEN(12,FILE=TRIM(ONAME),FORM='FORMATTED',STATUS='NEW',ACTION='WRITE',IOSTAT=OPEN_ERROR)
      IF (OPEN_ERROR > 0) STOP 'An error occurred while creating the file '//TRIM(ONAME)//', the file already exists'
      OPEN(13,FILE='CHGCAR',FORM='FORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=OPEN_ERROR)
      IF (OPEN_ERROR > 0) STOP 'An error occurred while opening the file CHGCAR'
      DO
        READ(13,'(A132)') READ_WRITE ; WRITE(12,*) TRIM(READ_WRITE) ; IF ((TRIM(READ_WRITE)).EQ.(TRIM(NGS))) EXIT
      ENDDO
      WRITE(12,*) PARCHGDATA3
      CLOSE(13) ; CLOSE(12)
      ONAME = TRIM(INAME)//'.BETA'
      OPEN(12,FILE=TRIM(ONAME),FORM='FORMATTED',STATUS='NEW',ACTION='WRITE',IOSTAT=OPEN_ERROR)
      IF (OPEN_ERROR > 0) STOP 'An error occurred while creating the file '//TRIM(ONAME)//', the file already exists'
      OPEN(13,FILE='CHGCAR',FORM='FORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=OPEN_ERROR)
      IF (OPEN_ERROR > 0) STOP 'An error occurred while opening the file CHGCAR'
      DO
        READ(13,'(A132)') READ_WRITE ; WRITE(12,*) TRIM(READ_WRITE) ; IF ((TRIM(READ_WRITE)).EQ.(TRIM(NGS))) EXIT
      ENDDO
      WRITE(12,*) PARCHGDATA4
      CLOSE(13) ; CLOSE(12)
      SELECT CASE (DIRECTION)
        CASE(1)
        ALLOCATE(PARCHGDATAAVG(NG1),PARCHGDATAAVG2(NG1))
        PARCHGDATAAVG=0.0 ; PARCHGDATAAVG2=0.0
        DO L=1,NG2
          DO M=1,NG3
            PARCHGDATAAVG(:) = PARCHGDATAAVG(:) + PARCHGDATA3(:,L,M)/(NG2*NG3)
            PARCHGDATAAVG2(:) = PARCHGDATAAVG2(:) + PARCHGDATA4(:,L,M)/(NG2*NG3)
          ENDDO
        ENDDO
        ONAME = TRIM(INAME)//'.ALPHA.R1' ; ONAME2 = TRIM(INAME)//'.BETA.R1'
        CASE(2)
        ALLOCATE(PARCHGDATAAVG(NG2),PARCHGDATAAVG2(NG2))
        PARCHGDATAAVG=0.0 ; PARCHGDATAAVG2=0.0
        DO L=1,NG1
          DO M=1,NG3
            PARCHGDATAAVG(:) = PARCHGDATAAVG(:) + PARCHGDATA3(L,:,M)/(NG1*NG3)
            PARCHGDATAAVG2(:) = PARCHGDATAAVG2(:) + PARCHGDATA4(L,:,M)/(NG1*NG3)
          ENDDO
        ENDDO
        ONAME = TRIM(INAME)//'.ALPHA.R2' ; ONAME2 = TRIM(INAME)//'.BETA.R2'
        CASE(3)
        ALLOCATE(PARCHGDATAAVG(NG3),PARCHGDATAAVG2(NG3))
        PARCHGDATAAVG=0.0 ; PARCHGDATAAVG2=0.0
        DO L=1,NG1
          DO M=1,NG2
            PARCHGDATAAVG(:) = PARCHGDATAAVG(:) + PARCHGDATA3(L,M,:)/(NG1*NG2)
            PARCHGDATAAVG2(:) = PARCHGDATAAVG2(:) + PARCHGDATA4(L,M,:)/(NG1*NG2)
          ENDDO
        ENDDO
        ONAME = TRIM(INAME)//'.ALPHA.R3' ; ONAME2 = TRIM(INAME)//'.BETA.R3'
      END SELECT
      OPEN(12,FILE=TRIM(ONAME),FORM='FORMATTED',STATUS='NEW',ACTION='WRITE',IOSTAT=OPEN_ERROR)
      IF (OPEN_ERROR > 0) STOP 'An error occurred while creating the file '//TRIM(ONAME)//', the file already exists'
      WRITE(12,*) PARCHGDATAAVG
      CLOSE(12) 
      OPEN(12,FILE=TRIM(ONAME2),FORM='FORMATTED',STATUS='NEW',ACTION='WRITE',IOSTAT=OPEN_ERROR)
      IF (OPEN_ERROR > 0) STOP 'An error occurred while creating the file '//TRIM(ONAME)//', the file already exists'
      WRITE(12,*) PARCHGDATAAVG2
      CLOSE(12) 
      DEALLOCATE(PARCHGDATAAVG,PARCHGDATAAVG2)
    ENDIF
  ENDDO
ENDDO
DEALLOCATE(PARCHGDATA,PARCHGDATA2,PARCHGDATA3,PARCHGDATA4) ; WRITE(*,*) ; WRITE(*,*) 

ELSEIF ((DOPARCHGSPIN).AND.(SPINCASE .EQ. 3)) THEN
WRITE(*,*) 'PARCHGSPIN is not activated. PARCHGSPIN in non-collinear case is not posible,'
WRITE(*,*) 'due to the spin polarization is not well defined, the spinors are no longer'
WRITE(*,*) 'spin-wavefunctions, i.e. pure Alpha or Beta spin functions. The magnetization'
WRITE(*,*) 'is not projected on the z direction.' ; WRITE(*,*) ; WRITE(*,*) 
ENDIF

IF ((DOCHGCARSPIN).AND.(SPINCASE .EQ. 1)) THEN
WRITE(*,*) 'CHGCARSPIN is not activated. CHGCARSPIN in collinear non-magnetic case has no'
WRITE(*,*) 'sense... The Alpha and Beta densities are equal by definition, and their value'
WRITE(*,*) 'is a half of the total density... :P' ; WRITE(*,*) ; WRITE(*,*) 

ELSEIF ((DOCHGCARSPIN).AND.(SPINCASE .EQ. 2)) THEN
WRITE(*,*) 'CHGCARSPIN is activated. Computing the Alpha and Beta Total Densities, and PCA'
WRITE(*,*) 'over each data.'
ALLOCATE(PARCHGDATA(NG1,NG2,NG3),PARCHGDATA2(NG1,NG2,NG3),PARCHGDATA3(NG1,NG2,NG3),PARCHGDATA4(NG1,NG2,NG3))
INQUIRE(FILE='CHGCAR', EXIST=FILE_EXIST)
IF (FILE_EXIST) THEN
  WRITE(*,*) 'Reading CHGCAR'
  OPEN(12,FILE='CHGCAR',FORM='FORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=OPEN_ERROR)
  IF (OPEN_ERROR > 0) STOP 'An error occurred while opening the file CHGCAR'
  DO
    READ(12,'(A132)') READ_WRITE ; IF ((TRIM(READ_WRITE)).EQ.(TRIM(NGS))) EXIT
  ENDDO
  READ(12,*) PARCHGDATA
  DO 
    READ(12,'(A132)') READ_WRITE ; IF ((TRIM(READ_WRITE)).EQ.(TRIM(NGS))) EXIT
  ENDDO
  READ(12,*) PARCHGDATA2
  CLOSE(12)
  PARCHGDATA3=0.5*(PARCHGDATA+PARCHGDATA2) ; PARCHGDATA4=0.5*(PARCHGDATA-PARCHGDATA2)
  OPEN(12,FILE='CHGCAR.ALPHA',FORM='FORMATTED',STATUS='NEW',ACTION='WRITE',IOSTAT=OPEN_ERROR)
  IF (OPEN_ERROR > 0) STOP 'An error occurred while creating the file CHGCAR.ALPHA, the file already exists'
  OPEN(13,FILE='CHGCAR',FORM='FORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=OPEN_ERROR)
  IF (OPEN_ERROR > 0) STOP 'An error occurred while opening the file CHGCAR'
  DO
    READ(13,'(A132)') READ_WRITE ; WRITE(12,*) TRIM(READ_WRITE) ; IF ((TRIM(READ_WRITE)).EQ.(TRIM(NGS))) EXIT
  ENDDO
  WRITE(12,*) PARCHGDATA3
  CLOSE(13) ; CLOSE(12)
  ONAME = 'CHGCAR.BETA'
  OPEN(12,FILE='CHGCAR.BETA',FORM='FORMATTED',STATUS='NEW',ACTION='WRITE',IOSTAT=OPEN_ERROR)
  IF (OPEN_ERROR > 0) STOP 'An error occurred while creating the file CHGCAR.BETA, the file already exists'
  OPEN(13,FILE='CHGCAR',FORM='FORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=OPEN_ERROR)
  IF (OPEN_ERROR > 0) STOP 'An error occurred while opening the file CHGCAR'
  DO
    READ(13,'(A132)') READ_WRITE ; WRITE(12,*) TRIM(READ_WRITE) ; IF ((TRIM(READ_WRITE)).EQ.(TRIM(NGS))) EXIT
  ENDDO
  WRITE(12,*) PARCHGDATA4
  CLOSE(13) ; CLOSE(12)
  SELECT CASE (DIRECTION)
    CASE(1)
    ALLOCATE(PARCHGDATAAVG(NG1),PARCHGDATAAVG2(NG1))
    PARCHGDATAAVG=0.0 ; PARCHGDATAAVG2=0.0
    DO L=1,NG2
      DO M=1,NG3
        PARCHGDATAAVG(:) = PARCHGDATAAVG(:) + PARCHGDATA3(:,L,M)/(NG2*NG3)
        PARCHGDATAAVG2(:) = PARCHGDATAAVG2(:) + PARCHGDATA4(:,L,M)/(NG2*NG3)
      ENDDO
    ENDDO
    ONAME = 'CHGCAR.ALPHA.R1' ; ONAME2 = 'CHGCAR.BETA.R1'
    CASE(2)
    ALLOCATE(PARCHGDATAAVG(NG2),PARCHGDATAAVG2(NG2))
    PARCHGDATAAVG=0.0 ; PARCHGDATAAVG2=0.0
    DO L=1,NG1
      DO M=1,NG3
        PARCHGDATAAVG(:) = PARCHGDATAAVG(:) + PARCHGDATA3(L,:,M)/(NG1*NG3)
        PARCHGDATAAVG2(:) = PARCHGDATAAVG2(:) + PARCHGDATA4(L,:,M)/(NG1*NG3)
      ENDDO
    ENDDO
    ONAME = 'CHGCAR.ALPHA.R2' ; ONAME2 = 'CHGCAR.BETA.R2'
    CASE(3)
    ALLOCATE(PARCHGDATAAVG(NG3),PARCHGDATAAVG2(NG3))
    PARCHGDATAAVG=0.0 ; PARCHGDATAAVG2=0.0
    DO L=1,NG1
      DO M=1,NG2
        PARCHGDATAAVG(:) = PARCHGDATAAVG(:) + PARCHGDATA3(L,M,:)/(NG1*NG2)
        PARCHGDATAAVG2(:) = PARCHGDATAAVG2(:) + PARCHGDATA4(L,M,:)/(NG1*NG2)
      ENDDO
    ENDDO
    ONAME = 'CHGCAR.ALPHA.R3' ; ONAME2 = 'CHGCAR.BETA.R3'
  END SELECT
  OPEN(12,FILE=TRIM(ONAME),FORM='FORMATTED',STATUS='NEW',ACTION='WRITE',IOSTAT=OPEN_ERROR)
  IF (OPEN_ERROR > 0) STOP 'An error occurred while creating the file '//TRIM(ONAME)//', the file already exists'
  WRITE(12,*) PARCHGDATAAVG
  CLOSE(12) 
  OPEN(12,FILE=TRIM(ONAME2),FORM='FORMATTED',STATUS='NEW',ACTION='WRITE',IOSTAT=OPEN_ERROR)
  IF (OPEN_ERROR > 0) STOP 'An error occurred while creating the file '//TRIM(ONAME)//', the file already exists'
  WRITE(12,*) PARCHGDATAAVG2
  CLOSE(12) 
  DEALLOCATE(PARCHGDATAAVG,PARCHGDATAAVG2)
ENDIF
DEALLOCATE(PARCHGDATA,PARCHGDATA2,PARCHGDATA3,PARCHGDATA4) ; WRITE(*,*) ; WRITE(*,*) 

ELSE IF ((DOCHGCARSPIN).AND.(SPINCASE .EQ. 3)) THEN
WRITE(*,*) 'CHGCARSPIN is not activated. CHGCARSPIN in non-collinear case is not posible,'
WRITE(*,*) 'due to the spin polarization is not well defined, the spinors are no longer'
WRITE(*,*) 'spin-wavefunctions, i.e. no pure Alpha or Beta spin functions. The spin is'
WRITE(*,*) 'not projected on the z direction.' ; WRITE(*,*) ; WRITE(*,*) 
ENDIF

IF (DOCHGCARAVG) THEN
WRITE(*,*) 'CHGCARAVG is activated. Computing planar average for CHGCAR density.'
ALLOCATE(PARCHGDATA(NG1,NG2,NG3),PARCHGDATA2(NG1,NG2,NG3),PARCHGDATA3(NG1,NG2,NG3),PARCHGDATA4(NG1,NG2,NG3))
INQUIRE(FILE='CHGCAR', EXIST=FILE_EXIST)
IF (FILE_EXIST) THEN
  WRITE(*,*) 'Reading CHGCAR'
  OPEN(12,FILE='CHGCAR',FORM='FORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=OPEN_ERROR)
  IF (OPEN_ERROR > 0) STOP 'An error occurred while opening the file CHGCAR'
  DO
    READ(12,'(A132)') READ_WRITE ; IF ((TRIM(READ_WRITE)).EQ.(TRIM(NGS))) EXIT
  ENDDO
  READ(12,*) PARCHGDATA
  IF (SPINCASE .EQ. 1) THEN
    WRITE(*,*) 'There is no spin density matrix in SPINCASE = 1 (non-magnetic case).'
    CLOSE(12)
  ELSE IF (SPINCASE .EQ. 2) THEN
    WRITE(*,*) 'SPINCASE = 2 (magnetic case), searching spin density matrix.'
    DO
      READ(12,'(A132)') READ_WRITE ; IF ((TRIM(READ_WRITE)).EQ.(TRIM(NGS))) EXIT
    ENDDO
    READ(12,*) PARCHGDATA2
    CLOSE(12)
    OPEN(14,FILE='CHGCAR.S',FORM='FORMATTED',STATUS='NEW',ACTION='WRITE',IOSTAT=OPEN_ERROR)
    IF (OPEN_ERROR > 0) STOP 'An error occurred while creating the file CHGCAR.S, the file already exists'
    OPEN(13,FILE='CHGCAR',FORM='FORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=OPEN_ERROR)
    IF (OPEN_ERROR > 0) STOP 'An error occurred while opening the file CHGCAR'
    DO
      READ(13,'(A132)') READ_WRITE ; WRITE(14,*) TRIM(READ_WRITE) ; IF ((TRIM(READ_WRITE)).EQ.(TRIM(NGS))) EXIT
    ENDDO
    WRITE(14,*) PARCHGDATA2
    CLOSE(13) ; CLOSE(14)
  ELSE IF (SPINCASE .EQ. 3) THEN
    WRITE(*,*) 'SPINCASE = 3 (non-collinear case), searching spin density matrices (3).'
    DO
      READ(12,'(A132)') READ_WRITE ; IF ((TRIM(READ_WRITE)).EQ.(TRIM(NGS))) EXIT
    ENDDO
    READ(12,*) PARCHGDATA2
    DO
      READ(12,'(A132)') READ_WRITE ; IF ((TRIM(READ_WRITE)).EQ.(TRIM(NGS))) EXIT
    ENDDO
    READ(12,*) PARCHGDATA3
    DO
      READ(12,'(A132)') READ_WRITE ; IF ((TRIM(READ_WRITE)).EQ.(TRIM(NGS))) EXIT
    ENDDO
    READ(12,*) PARCHGDATA4
    CLOSE(12)
    OPEN(14,FILE='CHGCAR.S1',FORM='FORMATTED',STATUS='NEW',ACTION='WRITE',IOSTAT=OPEN_ERROR)
    IF (OPEN_ERROR > 0) STOP 'An error occurred while creating the file CHGCAR.S1, the file already exists'
    OPEN(13,FILE='CHGCAR',FORM='FORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=OPEN_ERROR)
    IF (OPEN_ERROR > 0) STOP 'An error occurred while opening the file CHGCAR'
    DO
      READ(13,'(A132)') READ_WRITE ; WRITE(14,*) TRIM(READ_WRITE) ; IF ((TRIM(READ_WRITE)).EQ.(TRIM(NGS))) EXIT
    ENDDO
    WRITE(14,*) PARCHGDATA2
    CLOSE(13) ; CLOSE(14)
    OPEN(14,FILE='CHGCAR.S2',FORM='FORMATTED',STATUS='NEW',ACTION='WRITE',IOSTAT=OPEN_ERROR)
    IF (OPEN_ERROR > 0) STOP 'An error occurred while creating the file CHGCAR.S2, the file already exists'
    OPEN(13,FILE='CHGCAR',FORM='FORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=OPEN_ERROR)
    IF (OPEN_ERROR > 0) STOP 'An error occurred while opening the file CHGCAR'
    DO
      READ(13,'(A132)') READ_WRITE ; WRITE(14,*) TRIM(READ_WRITE) ; IF ((TRIM(READ_WRITE)).EQ.(TRIM(NGS))) EXIT
    ENDDO
    WRITE(14,*) PARCHGDATA3
    CLOSE(13) ; CLOSE(14)
    OPEN(14,FILE='CHGCAR.S3',FORM='FORMATTED',STATUS='NEW',ACTION='WRITE',IOSTAT=OPEN_ERROR)
    IF (OPEN_ERROR > 0) STOP 'An error occurred while creating the file CHGCAR.S2, the file already exists'
    OPEN(13,FILE='CHGCAR',FORM='FORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=OPEN_ERROR)
    IF (OPEN_ERROR > 0) STOP 'An error occurred while opening the file CHGCAR'
    DO
      READ(13,'(A132)') READ_WRITE ; WRITE(14,*) TRIM(READ_WRITE) ; IF ((TRIM(READ_WRITE)).EQ.(TRIM(NGS))) EXIT
    ENDDO
    WRITE(14,*) PARCHGDATA4
    CLOSE(13) ; CLOSE(14)
  ENDIF
  SELECT CASE (DIRECTION)
    CASE(1)
    ALLOCATE(PARCHGDATAAVG(NG1),PARCHGDATAAVG2(NG1),PARCHGDATAAVG3(NG1),PARCHGDATAAVG4(NG1))
    PARCHGDATAAVG=0.0 ; PARCHGDATAAVG2=0.0 ; PARCHGDATAAVG3=0.0 ; PARCHGDATAAVG4=0.0
    IF (SPINCASE .EQ. 1) THEN
      ONAME = 'CHGCAR.R1'
      DO L=1,NG2
        DO M=1,NG3
          PARCHGDATAAVG(:) = PARCHGDATAAVG(:) + PARCHGDATA(:,L,M)/(NG2*NG3)
        ENDDO
      ENDDO
    ELSE IF (SPINCASE .EQ. 2) THEN
      ONAME = 'CHGCAR.R1' ; ONAME2 = 'CHGCAR.SR1'
      DO L=1,NG2
        DO M=1,NG3
          PARCHGDATAAVG(:) = PARCHGDATAAVG(:) + PARCHGDATA(:,L,M)/(NG2*NG3)
          PARCHGDATAAVG2(:) = PARCHGDATAAVG2(:) + PARCHGDATA2(:,L,M)/(NG2*NG3)
        ENDDO
      ENDDO
    ELSE IF (SPINCASE .EQ. 3) THEN
      ONAME = 'CHGCAR.R1' ; ONAME2 = 'CHGCAR.S1R1' ; ONAME3 = 'CHGCAR.S2R1' ; ONAME4 = 'CHGCAR.S3R1'
      DO L=1,NG2
        DO M=1,NG3
          PARCHGDATAAVG(:) = PARCHGDATAAVG(:) + PARCHGDATA(:,L,M)/(NG2*NG3)
          PARCHGDATAAVG2(:) = PARCHGDATAAVG2(:) + PARCHGDATA2(:,L,M)/(NG2*NG3)
          PARCHGDATAAVG3(:) = PARCHGDATAAVG3(:) + PARCHGDATA3(:,L,M)/(NG2*NG3)
          PARCHGDATAAVG4(:) = PARCHGDATAAVG4(:) + PARCHGDATA4(:,L,M)/(NG2*NG3)
        ENDDO
      ENDDO
    ENDIF
    CASE(2)
    ALLOCATE(PARCHGDATAAVG(NG2),PARCHGDATAAVG2(NG2),PARCHGDATAAVG3(NG2),PARCHGDATAAVG4(NG2))
    PARCHGDATAAVG=0.0 ; PARCHGDATAAVG2=0.0 ; PARCHGDATAAVG3=0.0 ; PARCHGDATAAVG4=0.0
    IF (SPINCASE .EQ. 1) THEN
      ONAME = 'CHGCAR.R2'
      DO L=1,NG1
        DO M=1,NG3
          PARCHGDATAAVG(:) = PARCHGDATAAVG(:) + PARCHGDATA(L,:,M)/(NG1*NG3)
        ENDDO
      ENDDO
    ELSE IF (SPINCASE .EQ. 2) THEN
      ONAME = 'CHGCAR.R2' ; ONAME2 = 'CHGCAR.SR2'
      DO L=1,NG1
        DO M=1,NG3
          PARCHGDATAAVG(:) = PARCHGDATAAVG(:) + PARCHGDATA(L,:,M)/(NG1*NG3)
          PARCHGDATAAVG2(:) = PARCHGDATAAVG2(:) + PARCHGDATA2(L,:,M)/(NG1*NG3)
        ENDDO
      ENDDO
    ELSE IF (SPINCASE .EQ. 3) THEN
      ONAME = 'CHGCAR.R2' ; ONAME2 = 'CHGCAR.S1R2' ; ONAME3 = 'CHGCAR.S2R2' ; ONAME4 = 'CHGCAR.S3R2'
      DO L=1,NG1
        DO M=1,NG3
          PARCHGDATAAVG(:) = PARCHGDATAAVG(:) + PARCHGDATA(L,:,M)/(NG1*NG3)
          PARCHGDATAAVG2(:) = PARCHGDATAAVG2(:) + PARCHGDATA2(L,:,M)/(NG1*NG3)
          PARCHGDATAAVG3(:) = PARCHGDATAAVG3(:) + PARCHGDATA3(L,:,M)/(NG1*NG3)
          PARCHGDATAAVG4(:) = PARCHGDATAAVG4(:) + PARCHGDATA4(L,:,M)/(NG1*NG3)
        ENDDO
      ENDDO
    ENDIF
    CASE(3)
    ALLOCATE(PARCHGDATAAVG(NG3),PARCHGDATAAVG2(NG3),PARCHGDATAAVG3(NG3),PARCHGDATAAVG4(NG3))
    PARCHGDATAAVG=0.0 ; PARCHGDATAAVG2=0.0 ; PARCHGDATAAVG3=0.0 ; PARCHGDATAAVG4=0.0
    IF (SPINCASE .EQ. 1) THEN
      ONAME = 'CHGCAR.R3'
      DO L=1,NG1
        DO M=1,NG2
          PARCHGDATAAVG(:) = PARCHGDATAAVG(:) + PARCHGDATA(L,M,:)/(NG1*NG2)
        ENDDO
      ENDDO
    ELSE IF (SPINCASE .EQ. 2) THEN
      ONAME = 'CHGCAR.R3' ; ONAME2 = 'CHGCAR.SR3'
      DO L=1,NG1
        DO M=1,NG2
          PARCHGDATAAVG(:) = PARCHGDATAAVG(:) + PARCHGDATA(L,M,:)/(NG1*NG2)
          PARCHGDATAAVG2(:) = PARCHGDATAAVG2(:) + PARCHGDATA2(L,M,:)/(NG1*NG2)
        ENDDO
      ENDDO
    ELSE IF (SPINCASE .EQ. 3) THEN
      ONAME = 'CHGCAR.R3' ; ONAME2 = 'CHGCAR.S1R3' ; ONAME3 = 'CHGCAR.S2R3' ; ONAME4 = 'CHGCAR.S3R3'
      DO L=1,NG1
        DO M=1,NG2
          PARCHGDATAAVG(:) = PARCHGDATAAVG(:) + PARCHGDATA(L,M,:)/(NG1*NG2)
          PARCHGDATAAVG2(:) = PARCHGDATAAVG2(:) + PARCHGDATA2(L,M,:)/(NG1*NG2)
          PARCHGDATAAVG3(:) = PARCHGDATAAVG3(:) + PARCHGDATA3(L,M,:)/(NG1*NG2)
          PARCHGDATAAVG4(:) = PARCHGDATAAVG4(:) + PARCHGDATA4(L,M,:)/(NG1*NG2)
        ENDDO
      ENDDO
    ENDIF
  END SELECT
  OPEN(12,FILE=TRIM(ONAME),FORM='FORMATTED',STATUS='NEW',ACTION='WRITE',IOSTAT=OPEN_ERROR)
  IF (OPEN_ERROR > 0) STOP 'An error occurred while creating the file '//TRIM(ONAME)//', the file already exists'
  WRITE(12,*) PARCHGDATAAVG
  CLOSE(12) 
  IF ((SPINCASE .EQ. 2).OR.(SPINCASE .EQ. 3)) THEN
    OPEN(12,FILE=TRIM(ONAME2),FORM='FORMATTED',STATUS='NEW',ACTION='WRITE',IOSTAT=OPEN_ERROR)
    IF (OPEN_ERROR > 0) STOP 'An error occurred while creating the file '//TRIM(ONAME)//', the file already exists'
    WRITE(12,*) PARCHGDATAAVG2
    CLOSE(12)
    IF (SPINCASE .EQ. 3) THEN
      OPEN(12,FILE=TRIM(ONAME3),FORM='FORMATTED',STATUS='NEW',ACTION='WRITE',IOSTAT=OPEN_ERROR)
      IF (OPEN_ERROR > 0) STOP 'An error occurred while creating the file '//TRIM(ONAME)//', the file already exists'
      WRITE(12,*) PARCHGDATAAVG3
      CLOSE(12)
      OPEN(12,FILE=TRIM(ONAME4),FORM='FORMATTED',STATUS='NEW',ACTION='WRITE',IOSTAT=OPEN_ERROR)
      IF (OPEN_ERROR > 0) STOP 'An error occurred while creating the file '//TRIM(ONAME)//', the file already exists'
      WRITE(12,*) PARCHGDATAAVG4
      CLOSE(12)
    ENDIF
  ENDIF
  DEALLOCATE(PARCHGDATAAVG,PARCHGDATAAVG2,PARCHGDATAAVG3,PARCHGDATAAVG4)
  IF (SPINCASE .EQ. 2) THEN
    WRITE(*,*) 'Calculating Total Magnetization (M) and Absolut Magnetization (|M|) magnitudes,'
    WRITE(*,*) 'in Bohr magneton unit'
    MAGMOM=0.0D0 ; ABSMAGMOM=0.0D0
    DO I=1,NG1
      DO J=1,NG2
        DO K=1,NG3
          MAGMOM = MAGMOM + PARCHGDATA2(I,J,K)/NPOINTS ; ABSMAGMOM = ABSMAGMOM + ABS(PARCHGDATA2(I,J,K))/NPOINTS
        ENDDO
      ENDDO
    ENDDO
    WRITE(*,*) 'MAGNETIZATION = ',MAGMOM ; WRITE(*,*) 'ABSOLUTE_MAGNETIZATION = ',ABSMAGMOM
  ELSEIF (SPINCASE .EQ. 3) THEN
    WRITE(*,*) 'Calculating Total Magnetization (M) and Absolut Magnetization (|M|) magnitudes,'
    WRITE(*,*) 'for each spin projection, in Bohr magneton unit'
    MAGMOM=0.0D0 ; ABSMAGMOM=0.0D0
    DO I=1,NG1
      DO J=1,NG2
        DO K=1,NG3
          MAGMOM = MAGMOM + PARCHGDATA2(I,J,K)/NPOINTS ; ABSMAGMOM = ABSMAGMOM + ABS(PARCHGDATA2(I,J,K))/NPOINTS
        ENDDO
      ENDDO
    ENDDO
    WRITE(*,*) 'MAGNETIZATION(1) = ',MAGMOM ; WRITE(*,*) 'ABSOLUTE_MAGNETIZATION(1) = ',ABSMAGMOM
    MAGMOM=0.0D0 ; ABSMAGMOM=0.0D0
    DO I=1,NG1
      DO J=1,NG2
        DO K=1,NG3
          MAGMOM = MAGMOM + PARCHGDATA3(I,J,K)/NPOINTS ; ABSMAGMOM = ABSMAGMOM + ABS(PARCHGDATA3(I,J,K))/NPOINTS
        ENDDO
      ENDDO
    ENDDO
    WRITE(*,*) 'MAGNETIZATION(2) = ',MAGMOM ; WRITE(*,*) 'ABSOLUTE_MAGNETIZATION(2) = ',ABSMAGMOM
    MAGMOM=0.0D0 ; ABSMAGMOM=0.0D0
    DO I=1,NG1
      DO J=1,NG2
        DO K=1,NG3
          MAGMOM = MAGMOM + PARCHGDATA4(I,J,K)/NPOINTS ; ABSMAGMOM = ABSMAGMOM + ABS(PARCHGDATA4(I,J,K))/NPOINTS
        ENDDO
      ENDDO
    ENDDO
    WRITE(*,*) 'MAGNETIZATION(3) = ',MAGMOM ; WRITE(*,*) 'ABSOLUTE_MAGNETIZATION(3) = ',ABSMAGMOM
  ENDIF
ENDIF
DEALLOCATE(PARCHGDATA,PARCHGDATA2,PARCHGDATA3,PARCHGDATA4) ; WRITE(*,*) ; WRITE(*,*)
ENDIF

IF ((DOLDOSMAG).AND.((SPINCASE .EQ. 1).OR.(SPINCASE .EQ. 3))) THEN
WRITE(*,*) 'LDOSMAG is not activated. LDOSMAG subroutine is only for SPINCASE = 2' ; WRITE(*,*) ; WRITE(*,*)

ELSEIF ((DOLDOSMAG).AND.(SPINCASE .EQ. 2)) THEN
WRITE(*,*) 'LDOSMAG activated'
DE=(EMAX-EMIN)/NEN ; ETHRCUTOFF=ETHR/(SIGMA*SQRT(2*PI))
WRITE(*,*) 'Threshold Cutoff=', ETHRCUTOFF ; WRITE(*,*) 'Format: Reading "file name"'
WRITE(*,*) '        "Smearing value", "Eigen-energy", "Energy"'
WRITE(*,*) '        Is smearing value > Threshold Cutoff? YES or NO (not printed)'
SELECT CASE (DIRECTION)
  CASE(1) ; NGAVG=NG1 ; CASE(2) ; NGAVG=NG2 ; CASE(3) ; NGAVG=NG3
END SELECT
ALLOCATE(PARCHGDATAAVG2(NGAVG))
DO I=1,NEN+1
  E = DE*(I-1)+EMIN
  PARCHGDATAAVG2 = 0.0
  ALLOCATE(PARCHGDATAAVG(NGAVG))
  DO K=1,NKP
    DO J=1,NBANDS
      WRITE(BAND,'(I4.4)') J ; WRITE(KPOINT,'(I4.4)') K ; WRITE(DIR,'(I1.1)') DIRECTION
      INAME = 'PARCHG.'//TRIM(BAND)//'.'//TRIM(KPOINT)//'.ALPHA.R'//TRIM(DIR)
      INQUIRE(FILE=TRIM(INAME), EXIST=FILE_EXIST)
      IF (FILE_EXIST) THEN
        SMEARING = (1.0/SQRT(2.0*PI*SIGMA**2.0))*(EE**(-((E-ENERGY(K,J,1))**2.0)/(2.0*SIGMA**2.0)))
        WRITE(*,*) 'Reading: ',TRIM(INAME) ; WRITE(*,*) SMEARING, ENERGY(K,J,1), E
        IF (SMEARING.GT.(ETHRCUTOFF)) THEN
          WRITE(*,*) 'YES'
          OPEN(12,FILE=TRIM(INAME),FORM='FORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=OPEN_ERROR)
          IF (OPEN_ERROR > 0) STOP 'An error occurred while opening the file '//TRIM(INAME)
          READ(12,*) PARCHGDATAAVG
          CLOSE(12)
          PARCHGDATAAVG2 = PARCHGDATAAVG2 + SMEARING*W(K)*PARCHGDATAAVG/VOLUME
        ENDIF
      ENDIF
      INAME = 'PARCHG.'//TRIM(BAND)//'.'//TRIM(KPOINT)//'.BETA.R'//TRIM(DIR)
      INQUIRE(FILE=TRIM(INAME), EXIST=FILE_EXIST)
      IF (FILE_EXIST) THEN
        SMEARING = (1.0/SQRT(2.0*PI*SIGMA**2.0))*(EE**(-((E-ENERGY(K,J,2))**2.0)/(2.0*SIGMA**2.0)))
        WRITE(*,*) 'Reading: ',TRIM(INAME) ; WRITE(*,*) SMEARING, ENERGY(K,J,2), E
        IF (SMEARING.GT.(ETHRCUTOFF)) THEN
          WRITE(*,*) 'YES'
          OPEN(12,FILE=TRIM(INAME),FORM='FORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=OPEN_ERROR)
          IF (OPEN_ERROR > 0) STOP 'An error occurred while opening the file '//TRIM(INAME)
          READ(12,*) PARCHGDATAAVG
          CLOSE(12)
          PARCHGDATAAVG2 = PARCHGDATAAVG2 + SMEARING*W(K)*PARCHGDATAAVG/VOLUME
        ENDIF
      ENDIF
    ENDDO
  ENDDO
  DEALLOCATE(PARCHGDATAAVG)
  WRITE (NUM,'(I4.4)') I-1 ; WRITE(DIR,'(I1.1)') DIRECTION
  ONAME = 'LDOSMAG.R'//TRIM(DIR)//'.'//TRIM(NUM)//'.dat'
  OPEN(12,FILE=TRIM(ONAME),FORM='FORMATTED',STATUS='NEW',ACTION='WRITE',IOSTAT=OPEN_ERROR)
  IF (OPEN_ERROR > 0) STOP 'An error occurred while creating the file '//TRIM(ONAME)//', the file already exists'
  DO J=1,NGAVG
    WRITE(12,*) REAL(J-1),'	',E,'	',PARCHGDATAAVG2(J)
  ENDDO
  CLOSE(12)
ENDDO
DEALLOCATE(PARCHGDATAAVG2) ; WRITE(*,*) ; WRITE(*,*) 
ENDIF

IF ((DOLSDOSMAG).AND.((SPINCASE .EQ. 1).OR.(SPINCASE .EQ. 3))) THEN
WRITE(*,*) 'LSDOSMAG is not activated. LSDOSMAG subroutine is only for SPINCASE = 2' ; WRITE(*,*) ; WRITE(*,*) 

ELSEIF ((DOLSDOSMAG).AND.(SPINCASE .EQ. 2)) THEN
WRITE(*,*) 'LSDOSMAG activated, computing LSDOSMAG in collinear magnetic case.'
DE=(EMAX-EMIN)/NEN ; ETHRCUTOFF=ETHR/(SIGMA*SQRT(2*PI))
WRITE(*,*) 'Threshold Cutoff=', ETHRCUTOFF ; WRITE(*,*) 'Format: Reading "file name"'
WRITE(*,*) '        "Smearing value", "Eigen-energy", "Energy"'
WRITE(*,*) '        Is smearing value > Threshold Cutoff? YES or NO (not printed)'
SELECT CASE (DIRECTION)
  CASE(1) ; NGAVG=NG1 ; CASE(2) ; NGAVG=NG2 ; CASE(3) ; NGAVG=NG3
END SELECT
ALLOCATE(PARCHGDATAAVG2(NGAVG))
DO I=1,NEN+1
  E = DE*(I-1)+EMIN
  PARCHGDATAAVG2 = 0.0
  ALLOCATE(PARCHGDATAAVG(NGAVG))
  DO K=1,NKP
    DO J=1,NBANDS
      WRITE(BAND,'(I4.4)') J ; WRITE(KPOINT,'(I4.4)') K ; WRITE(DIR,'(I1.1)') DIRECTION
      INAME = 'PARCHG.'//TRIM(BAND)//'.'//TRIM(KPOINT)//'.ALPHA.R'//TRIM(DIR)
      INQUIRE(FILE=TRIM(INAME), EXIST=FILE_EXIST)
      IF (FILE_EXIST) THEN
        SMEARING = (1.0/SQRT(2.0*PI*SIGMA**2.0))*(EE**(-((E-ENERGY(K,J,1))**2.0)/(2.0*SIGMA**2.0)))
        WRITE(*,*) 'Reading: ',TRIM(INAME) ; WRITE(*,*) SMEARING, ENERGY(K,J,1), E
        IF (SMEARING.GT.(ETHRCUTOFF)) THEN
          WRITE(*,*) 'YES'
          OPEN(12,FILE=TRIM(INAME),FORM='FORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=OPEN_ERROR)
          IF (OPEN_ERROR > 0) STOP 'An error occurred while opening the file '//TRIM(INAME)
          READ(12,*) PARCHGDATAAVG
          CLOSE(12)
          PARCHGDATAAVG2 = PARCHGDATAAVG2 + SMEARING*W(K)*PARCHGDATAAVG/VOLUME
        ENDIF
      ENDIF
      INAME = 'PARCHG.'//TRIM(BAND)//'.'//TRIM(KPOINT)//'.BETA.R'//TRIM(DIR)
      INQUIRE(FILE=TRIM(INAME), EXIST=FILE_EXIST)
      IF (FILE_EXIST) THEN
        SMEARING = (1.0/SQRT(2.0*PI*SIGMA**2.0))*(EE**(-((E-ENERGY(K,J,2))**2.0)/(2.0*SIGMA**2.0)))
        WRITE(*,*) 'Reading: ',TRIM(INAME) ; WRITE(*,*) SMEARING, ENERGY(K,J,2), E
        IF (SMEARING.GT.(ETHRCUTOFF)) THEN
          WRITE(*,*) 'YES'
          OPEN(12,FILE=TRIM(INAME),FORM='FORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=OPEN_ERROR)
          IF (OPEN_ERROR > 0) STOP 'An error occurred while opening the file '//TRIM(INAME)
          READ(12,*) PARCHGDATAAVG
          CLOSE(12)
          PARCHGDATAAVG2 = PARCHGDATAAVG2 - SMEARING*W(K)*PARCHGDATAAVG/VOLUME
        ENDIF
      ENDIF
    ENDDO
  ENDDO
  DEALLOCATE(PARCHGDATAAVG)
  WRITE (NUM,'(I4.4)') I-1 ; WRITE(DIR,'(I1.1)') DIRECTION
  ONAME = 'LSDOSMAG.R'//TRIM(DIR)//'.'//TRIM(NUM)//'.dat'
  OPEN(12,FILE=TRIM(ONAME),FORM='FORMATTED',STATUS='NEW',ACTION='WRITE',IOSTAT=OPEN_ERROR)
  IF (OPEN_ERROR > 0) STOP 'An error occurred while creating the file '//TRIM(ONAME)//', the file already exists'
  DO J=1,NGAVG
    WRITE(12,*) REAL(J-1),'	',E,'	',PARCHGDATAAVG2(J)
  ENDDO
  CLOSE(12)
ENDDO
DEALLOCATE(PARCHGDATAAVG2) ; WRITE(*,*) ; WRITE(*,*) 
ENDIF

IF ((DOLDOSFULLMAG).AND.((SPINCASE .EQ. 1).OR.(SPINCASE .EQ. 3))) THEN
WRITE(*,*) 'LDOSFULLMAG is not activated. LDOSFULLMAG subroutine is only for SPINCASE = 2' ; WRITE(*,*) ; WRITE(*,*) 

ELSEIF ((DOLDOSFULLMAG).AND.(SPINCASE .EQ. 2)) THEN
WRITE(*,*) 'LDOSFULLMAG activated'
DE=(EMAX-EMIN)/NEN ; ETHRCUTOFF=ETHR/(SIGMA*SQRT(2*PI))
WRITE(*,*) 'Threshold Cutoff=', ETHRCUTOFF ; WRITE(*,*) 'Format: Reading "file name"'
WRITE(*,*) '        "Smearing value", "Eigen-energy", "Energy"'
WRITE(*,*) '        Is smearing value > Threshold Cutoff? YES or NO (not printed)'
ALLOCATE(PARCHGDATA2(NG1,NG2,NG3))
DO I=1,NEN+1
  E = DE*(I-1)+EMIN
  PARCHGDATA2 = 0.0
  ALLOCATE(PARCHGDATA(NG1,NG2,NG3))
  DO K=1,NKP
    DO J=1,NBANDS
      WRITE(BAND,'(I4.4)') J ; WRITE(KPOINT,'(I4.4)') K
      INAME = 'PARCHG.'//TRIM(BAND)//'.'//TRIM(KPOINT)//'.ALPHA'
      INQUIRE(FILE=TRIM(INAME), EXIST=FILE_EXIST)
      IF (FILE_EXIST) THEN
        SMEARING = (1.0/SQRT(2.0*PI*SIGMA**2.0))*(EE**(-((E-ENERGY(K,J,1))**2.0)/(2.0*SIGMA**2.0)))
        WRITE(*,*) 'Reading: ',TRIM(INAME) ; WRITE(*,*) SMEARING, ENERGY(K,J,1), E
        IF (SMEARING.GT.(ETHRCUTOFF)) THEN
          WRITE(*,*) 'YES'
          OPEN(12,FILE=TRIM(INAME),FORM='FORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=OPEN_ERROR)
          IF (OPEN_ERROR > 0) STOP 'An error occurred while opening the file '//TRIM(INAME)
          DO
            READ(12,'(A132)') READ_WRITE ; IF ((TRIM(READ_WRITE)).EQ.(' '//TRIM(NGS))) EXIT
          ENDDO
          READ(12,*) PARCHGDATA
          CLOSE(12)
          PARCHGDATA2 = PARCHGDATA2 + SMEARING*W(K)*PARCHGDATA/VOLUME
        ENDIF
      ENDIF
      INAME = 'PARCHG.'//TRIM(BAND)//'.'//TRIM(KPOINT)//'.BETA'
      INQUIRE(FILE=TRIM(INAME), EXIST=FILE_EXIST)
      IF (FILE_EXIST) THEN
        SMEARING = (1.0/SQRT(2.0*PI*SIGMA**2.0))*(EE**(-((E-ENERGY(K,J,2))**2.0)/(2.0*SIGMA**2.0)))
        WRITE(*,*) 'Reading: ',TRIM(INAME) ; WRITE(*,*) SMEARING, ENERGY(K,J,2), E
        IF (SMEARING.GT.(ETHRCUTOFF)) THEN
          WRITE(*,*) 'YES'
          OPEN(12,FILE=TRIM(INAME),FORM='FORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=OPEN_ERROR)
          IF (OPEN_ERROR > 0) STOP 'An error occurred while opening the file '//TRIM(INAME)
          DO
            READ(12,'(A132)') READ_WRITE ; IF ((TRIM(READ_WRITE)).EQ.(' '//TRIM(NGS))) EXIT
          ENDDO
          READ(12,*) PARCHGDATA
          CLOSE(12)
          PARCHGDATA2 = PARCHGDATA2 + SMEARING*W(K)*PARCHGDATA/VOLUME
        ENDIF
      ENDIF
    ENDDO
  ENDDO
  DEALLOCATE(PARCHGDATA)
  WRITE (NUM,'(I4.4)') I-1
  ONAME = 'LDOSMAG.FULL.'//TRIM(NUM)//'.dat'
  OPEN(12,FILE=TRIM(ONAME),FORM='FORMATTED',STATUS='NEW',ACTION='WRITE',IOSTAT=OPEN_ERROR)
  IF (OPEN_ERROR > 0) STOP 'An error occurred while creating the file '//TRIM(ONAME)//', the file already exists'
  OPEN(13,FILE='CHGCAR',FORM='FORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=OPEN_ERROR)
  IF (OPEN_ERROR > 0) STOP 'An error occurred while opening the file CHGCAR'
  READ(13,'(A132)') READ_WRITE ; WRITE(12,*) TRIM(READ_WRITE),'     ENERGY = ',E
  DO
    READ(13,'(A132)') READ_WRITE ; WRITE(12,*) TRIM(READ_WRITE) ; IF ((TRIM(READ_WRITE)).EQ.(TRIM(NGS))) EXIT
  ENDDO
  WRITE(12,*) PARCHGDATA2
  CLOSE(13) ; CLOSE(12)
ENDDO
DEALLOCATE(PARCHGDATA2) ; WRITE(*,*) ; WRITE(*,*) 
ENDIF

IF ((DOLSDOSFULLMAG).AND.((SPINCASE .EQ. 1).OR.(SPINCASE .EQ. 3))) THEN
WRITE(*,*) 'LSDOSFULLMAG is not activated. LSDOSFULLMAG subroutine is only for SPINCASE = 2' ; WRITE(*,*) ; WRITE(*,*) 

ELSEIF ((DOLSDOSFULLMAG).AND.(SPINCASE .EQ. 2)) THEN
WRITE(*,*) 'LSDOSFULLMAG activated, computing full LSDOS in collinear magnetic case.'
DE=(EMAX-EMIN)/NEN ; ETHRCUTOFF=ETHR/(SIGMA*SQRT(2*PI))
WRITE(*,*) 'Threshold Cutoff=', ETHRCUTOFF ; WRITE(*,*) 'Format: Reading "file name"'
WRITE(*,*) '        "Smearing value", "Eigen-energy", "Energy"'
WRITE(*,*) '        Is smearing value > Threshold Cutoff? YES or NO (not printed)'
ALLOCATE(PARCHGDATA2(NG1,NG2,NG3))
DO I=1,NEN+1
  E = DE*(I-1)+EMIN
  PARCHGDATA2 = 0.0
  ALLOCATE(PARCHGDATA(NG1,NG2,NG3))
  DO K=1,NKP
    DO J=1,NBANDS
      WRITE(BAND,'(I4.4)') J ; WRITE(KPOINT,'(I4.4)') K
      INAME = 'PARCHG.'//TRIM(BAND)//'.'//TRIM(KPOINT)//'.ALPHA'
      INQUIRE(FILE=TRIM(INAME), EXIST=FILE_EXIST)
      IF (FILE_EXIST) THEN
        SMEARING = (1.0/SQRT(2.0*PI*SIGMA**2.0))*(EE**(-((E-ENERGY(K,J,1))**2.0)/(2.0*SIGMA**2.0)))
        WRITE(*,*) 'Reading: ',TRIM(INAME) ; WRITE(*,*) SMEARING, ENERGY(K,J,1), E
        IF (SMEARING.GT.(ETHRCUTOFF)) THEN
          WRITE(*,*) 'YES'
          OPEN(12,FILE=TRIM(INAME),FORM='FORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=OPEN_ERROR)
          IF (OPEN_ERROR > 0) STOP 'An error occurred while opening the file '//TRIM(INAME)
          DO
            READ(12,'(A132)') READ_WRITE ; IF ((TRIM(READ_WRITE)).EQ.(' '//TRIM(NGS))) EXIT
          ENDDO
          READ(12,*) PARCHGDATA
          CLOSE(12)
          PARCHGDATA2 = PARCHGDATA2 + SMEARING*W(K)*PARCHGDATA/VOLUME
        ENDIF
      ENDIF
      INAME = 'PARCHG.'//TRIM(BAND)//'.'//TRIM(KPOINT)//'.BETA'
      INQUIRE(FILE=TRIM(INAME), EXIST=FILE_EXIST)
      IF (FILE_EXIST) THEN
        SMEARING = (1.0/SQRT(2.0*PI*SIGMA**2.0))*(EE**(-((E-ENERGY(K,J,2))**2.0)/(2.0*SIGMA**2.0)))
        WRITE(*,*) 'Reading: ',TRIM(INAME) ; WRITE(*,*) SMEARING, ENERGY(K,J,2), E
        IF (SMEARING.GT.(ETHRCUTOFF)) THEN
          WRITE(*,*) 'YES'
          OPEN(12,FILE=TRIM(INAME),FORM='FORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=OPEN_ERROR)
          IF (OPEN_ERROR > 0) STOP 'An error occurred while opening the file '//TRIM(INAME)
          DO
            READ(12,'(A132)') READ_WRITE ; IF ((TRIM(READ_WRITE)).EQ.(' '//TRIM(NGS))) EXIT
          ENDDO
          READ(12,*) PARCHGDATA
          CLOSE(12)
          PARCHGDATA2 = PARCHGDATA2 - SMEARING*W(K)*PARCHGDATA/VOLUME
        ENDIF
      ENDIF
    ENDDO
  ENDDO
  DEALLOCATE(PARCHGDATA)
  WRITE (NUM,'(I4.4)') I-1
  ONAME = 'LSDOSMAG.FULL.'//TRIM(NUM)//'.dat'
  OPEN(12,FILE=TRIM(ONAME),FORM='FORMATTED',STATUS='NEW',ACTION='WRITE',IOSTAT=OPEN_ERROR)
  IF (OPEN_ERROR > 0) STOP 'An error occurred while creating the file '//TRIM(ONAME)//', the file already exists'
  OPEN(13,FILE='CHGCAR',FORM='FORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=OPEN_ERROR)
  IF (OPEN_ERROR > 0) STOP 'An error occurred while opening the file CHGCAR'
  READ(13,'(A132)') READ_WRITE ; WRITE(12,*) TRIM(READ_WRITE),'     ENERGY = ',E
  DO
    READ(13,'(A132)') READ_WRITE ; WRITE(12,*) TRIM(READ_WRITE) ; IF ((TRIM(READ_WRITE)).EQ.(TRIM(NGS))) EXIT
  ENDDO
  WRITE(12,*) PARCHGDATA2
  CLOSE(13) ; CLOSE(12)
ENDDO
DEALLOCATE(PARCHGDATA2) ; WRITE(*,*) ; WRITE(*,*) 
ENDIF

CALL CPU_TIME(FINISH)

WRITE(*,*) 'Please cite the following article:'; WRITE(*,*)

WRITE(*,*) 'Normal termination of DensityTool.'
WRITE(*,*) 'Job starts at : ',TRIM(DATE)
WRITE(*,*) 'Job Time = ',FINISH-START,' seconds.'

END PROGRAM DENSITYTOOL
