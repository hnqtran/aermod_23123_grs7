      SUBROUTINE RXNRATES_GRS7
C***********************************************************************
C        RXNRATES_GRS7 Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE: Calculate the reaction rates for GRS chemistry
C
C        PROGRAMMER: CERC
C
C        DATE: November 2020
C
C        INPUTS:  
C
C        OUTPUTS: 
C      
C        NOTES:
C      
C        R1, R2, R3, R4, R5, R6, and R7 are reaction rates for:
C                 
C                 RPG + hv -> RP + ROC    (1)
C                 RP + NO -> NO2          (2)
C                 NO2 + hv -> NO + O3     (3)
C                 NO + O3 -> NO2          (4)
C                 RP + RP -> RP           (5)
C                 RP + NO2 -> SGN         (6)
C                 RP + NO2 -> SNGN        (7)
C
C        CALLED FROM: GRS7_CALC
C***********************************************************************

C     Variable Declarations
      USE MAIN1
      USE GRS7MOD
      IMPLICIT NONE
      CHARACTER MODNAM*12
      
      DOUBLE PRECISION, PARAMETER:: SCRNHT=1.2D0
      DOUBLE PRECISION, PARAMETER:: EXP1=4710.0D0
      DOUBLE PRECISION, PARAMETER:: CONST1=3.16D-3 ! 1/316
      DOUBLE PRECISION, PARAMETER:: CONST2=3.58D6
      DOUBLE PRECISION, PARAMETER:: CONST3A=8.0D-4
      DOUBLE PRECISION, PARAMETER:: CONST3B=7.4D-6
      DOUBLE PRECISION, PARAMETER:: EXP3=10.0D0
      DOUBLE PRECISION:: CONST4, EXP4
      DOUBLE PRECISION:: CHEMTK, CHEMPT, TEMPI
      DOUBLE PRECISION:: KP ! Activity coefficient for NO2
      DOUBLE PRECISION:: sinelev   ! sind(solar elevation angle)
      
      DOUBLE PRECISION:: XLON, LSIGN ! LONGITUDE of Met site and its sign (+/-)
      INTEGER::          EORW,WORE ! East or West Hemisphere
      
      INTEGER::NSCRNHT

C --- Variable Initializations
      MODNAM = 'RXNRATES_GRS7'
      
C --- Determine temperature (K) at screen height
      IF(TREFHT==SCRNHT)THEN
        !Reference height is screen height
        CHEMTK=TA
      ELSE
        !Interpolate from potential temperature
        CALL LOCATE(GRIDHT, 1, MXGLVL, SCRNHT, NSCRNHT)    
        CALL GINTRP(GRIDHT(NSCRNHT),GRIDPT(NSCRNHT),GRIDHT(NSCRNHT+1),
     &              GRIDPT(NSCRNHT+1),SCRNHT,CHEMPT)    
        !Convert to temperature
        CHEMTK=CHEMPT-GOVRCP*(SCRNHT+ZBASE)
      END IF
C     IF(INLNDEBUG)write(*,*)"HTdbg: RXNRATES_GRS7; CHEMTK",CHEMTK

C     Decode Longitude string read from SFC file
C---- Determine if the letter 'N' or 'n' is in the latitude field
c     IF(INLNDEBUG)write(*,*)"ALAT=",ALAT,"ALON=",ALON
      EORW = INDEX(ALON,'E') + INDEX(ALON,'e')
      WORE = INDEX(ALON,'W') + INDEX(ALAT,'w')
      READ( ALON, '(F9.1)',ERR=1000 ) XLON

      IF( EORW .NE. 0 )THEN
C        The longitude is in the eastern hemisphere; decode the longitude
         LSIGN = 1.0D0
      ELSE
C        The longitude may be in the western hemisphere
         IF( WORE .NE. 0 )THEN
            LSIGN = -1.0D0
         ELSE
C           Write a warning to the user - error decoding the latitude
            CALL ERRHDL(PATH,MODNAM,'W','835',ALON)
c           Assuming western hemisphere            
            LSIGN = -1.0D0
         END IF
      END IF
      XLON = LSIGN * XLON
c     IF(INLNDEBUG)write(*,*)"XLON=",XLON

C --- Calculate reaction rates
      TEMPI = 1.0/CHEMTK
      
C     R4: NO + O3 -> NO2          (4)
      IF (R4REF==1) THEN ! GRSM
         IF (QSW==0) THEN
            R3=0.0
         ELSE
            R3=8.0D-4*EXP(-10.0D0/QSW)+(7.4D-6*QSW) ! Same in GRSM; sec-1
         END IF
         R4 =4.405D-2*EXP(-1370/CHEMTK)  ! (ppb-1 sec-1); @273K R4 = 0.291D-3 ppb-1 sec-1
         R4 = MAX(R4,1.0D-6)             ! 1e-6 corresponds to temperature of around -150 degC; same as GRSM
      ELSE IF (R4REF==2) THEN         ! PNZR
c        Calculate solar elevation angles    
         CALL SUNAE( IYEAR, JDAY, IHOUR, XLAT, XLON)
c        IF(INLNDEBUG)write(*,*)"HTdbg:",
c    &           "Solar Elevation Angle",SOLELV,
c    &           "sin(SOLELV)=",sin(SOLELV*PI/180)
         sinelev = SIN(SOLELV*(PI/180))
         R3 = EXP(-0.575/(sinelev))      ! min-1
         R3 = R3/6.0D1                   ! sec-1         
         R4 = (9.24D5/CHEMTK)*EXP(-1.450D3/CHEMTK)! @ 273 K R4 = 16.703 ppm-1 min-1
         R4 = R4/6.0D4                   ! ppb-1 sec-1         
      END IF

      R1 = 0.0067*R3*EXP(-4.7*1000*((1/CHEMTK)-CONST1)); ! sec-1
      !R1=10000*exp(-4710/CHEMTK) * R3

      R2 = 5482.0*EXP(242.0/CHEMTK)      ! (ppm-1 min-1); @273K R2 = 13302 ppm-1 min-1
      R2 = R2/6.0D4                      ! ppb-1 sec-1         
      R5 = 10200.0                       ! (ppm-1 min-1)
      R5 = R5/6.0D4                      ! ppb-1 sec-1         
      R6 = 125.0                         ! (ppm-1 min-1)
      R6 = R6/6.0D4                      ! ppb-1 sec-1         
      R7 = 125.0                         ! (ppm-1 min-1)
      R7 = R7/6.0D4                      ! ppb-1 sec-1         

1000  CONTINUE
C     Write a warning to the user - error decoding the latitude
      CALL ERRHDL(PATH,MODNAM,'E','834',ALON)

      RETURN
      
      END SUBROUTINE RXNRATES_GRS7
      
      
      SUBROUTINE GETBGDEQUIL_GRS7
C***********************************************************************
C        GETBDGEQUIL Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE: Calculate the equilibrium background for GRS7 chemistry
C
C        PROGRAMMER: CERC
C
C        DATE: November 2020
C
C        INPUTS:  
C
C        OUTPUTS:
C
C        CALLED FROM: GRS7_CALC
C***********************************************************************
C     Variable Declarations
      USE MAIN1
      USE GRS7MOD
      IMPLICIT NONE
      CHARACTER MODNAM*12

C --- Variable Initializations
      MODNAM = 'GETBGDEQUIL_GRS7'  

      L_NO2pseudo = .False.

C     write(*,*)"HTdbg: GRSM_CALC; GETBGDEQUIL_GRS7; Inputs: ",
C    &          "|L_CalcNOxFromNO2=",L_CalcNOxFromNO2,
C    &          "|L_CalcNoxFromNO2=",L_CalcNoxFromNO2,
C    &          "|L_CalcNOXFromNO2=",L_CalcNOXFromNO2,
C    &          "|NOXMISS=", NOXMISS,
C    &          "|NOXCONC_BG=",NOXCONC_BG,
C    &          "|NO2CONC_BG=",NO2CONC_BG,
C    &          "|NOCONC_BG=",NOCONC_BG,
C    &          "|O3CONC_BG=",O3CONC_BG,
C    &          "VOCCONC_BG=",VOCCONC_BG,
C    &          "|ROCCONC_BG=",ROCCONC_BG,
C    &          "|RPCONC_BG=",RPCONC_BG

C --- Check for day or night
      IF (.NOT.L_NIGHTHOUR) THEN ! Daytime
C       IF(L_CalcNOxFromNO2.OR.NOXMISS.OR.L_NO2pseudo)THEN ! NOx background is missing
        IF(L_CalcNOxFromNO2.OR.NOXMISS)THEN                ! NOx background is missing
            IF(O3CONC_BG/=0.0D0)THEN        ! Background O3 is available
               IF(ROCCONC_BG==0.0D0)THEN    ! VOC/ROC background is missing
C ---             Calculate NOx bgd from NO2 same as for GRSM
C ---             Solve equilibrium eqns for NOx
C                 write(*,*)"HTdbg: GETBGDEQUIL_GRS7; Daytime ",
C    &                   "with NOX from NO2CONC_BG and ",
C    &                   "zero ROC background"
                  NOXCONC_BG=NO2CONC_BG/(O3CONC_BG)*(O3CONC_BG+R3/R4)
                  NOCONC_BG=NOXCONC_BG-NO2CONC_BG
               ELSE IF(ROCCONC_BG/=0.0D0)THEN  ! VOC/ROC background is available
C                 write(*,*)"HTdbg: GETBGDEQUIL_GRS7; ",
C    &                   "with NOX from NO2CONC_BG and ",
C    &                   "ROC background"
                  NO2CONC_BG = NOXCONC_BG
                  CALL QuadraticEquil_GRS7(NO2CONC_BG,NOCONC_BG,
     &                O3CONC_BG,ROCCONC_BG,RPCONC_BG)
               END IF
            ELSE  IF(O3CONC_BG==0.0D0)THEN  ! Background O3 is zero
               IF(ROCCONC_BG/=0.0D0)THEN    ! Backgound ROC is non-zero
                  RPCONC_BG = (R1*ROCCONC_BG - R3*NO2CONC_BG)/R5
                  NOCONC_BG = (R3*NO2CONC_BG)/(R2*RPCONC_BG)
                  NO2CONC_BG = NO2CONC_BG - NOCONC_BG
C                 NOXCONC_BG = NO2CONC_BG + NOCONC_BG
                  O3CONC_BG  = (R3*NO2CONC_BG)/(R4*NOCONC_BG)
               END IF
            END IF
        ELSE
C       IF NOx background is avaliable
            IF (NOXCONC_BG<NO2CONC_BG) THEN
               CALL ERRHDL(PATH, MODNAM, 'W','613','GRS7')
C              write(*,*)"HTdbg: GETBGDEQUIL_GRS7; NO2 > NOX"
               NOXCONC_BG = NO2CONC_BG 
            END IF
            NOCONC_BG = NOXCONC_BG - NO2CONC_BG
            IF(O3CONC_BG/=0.0D0)THEN        ! Background O3 is available
               IF(ROCCONC_BG==0.0D0)THEN    ! VOC/ROC background is missing
C ---             Use GRSM routine to solve NOX and O3
C                 write(*,*)"HTdbg: GETBGDEQUIL_GRS7; ",
C    &                   "with NOX from NOXCONC_BG and ",
C    &                   "zero ROC background"
                  CALL QuadraticEquil_GRSM(NO2CONC_BG,NOCONC_BG,
     &                                    O3CONC_BG)
               ELSE IF(ROCCONC_BG/=0.0D0)THEN  ! ROC background is available
C                 write(*,*)"HTdbg: GETBGDEQUIL_GRS7; ",
C    &                   "with NOX from NOXCONC_BG and ",
C    &                   "ROC background"
                  CALL QuadraticEquil_GRS7(NOCONC_BG,NO2CONC_BG,
     &                O3CONC_BG,ROCCONC_BG,RPCONC_BG)
               END IF
            ELSE  IF(O3CONC_BG==0.0D0)THEN  ! Background O3 is zero
               IF(ROCCONC_BG/=0.0D0)THEN    ! Backgound ROC is non-zero
                  RPCONC_BG = (R1*ROCCONC_BG - R3*NO2CONC_BG)/R5
                  NOCONC_BG = (R3*NO2CONC_BG)/(R2*RPCONC_BG)
                  NO2CONC_BG = NOXCONC_BG - NOCONC_BG
C                 NOXCONC_BG = NO2CONC_BG + NOCONC_BG
                  O3CONC_BG  = (R3*NO2CONC_BG)/(R4*NOCONC_BG)
               END IF
            END IF            
        END IF 

      ELSE ! Nighttime
C       write(*,*)"HTdbg: GETBGDEQUIL_GRS7; Nighttime NOX"
        IF(L_CalcNOxFromNO2.OR.NOXMISS)THEN
          !NOx is missing but don't need it
          NOXCONC_BG=-999.0D0
          NOCONC_BG=0.0D0
        ELSE  
          ! No adjustments made to background at night
          NOCONC_BG = NOXCONC_BG - NO2CONC_BG
        END IF
        RPCONC_BG = 0.0D0
      END IF
C     write(*,*)"HTdbg: GRSM_CALC; GETBGDEQUIL_GRS7; Outputs: ",
C    &          "|L_CalcNOxFromNO2=",L_CalcNOxFromNO2,
C    &          "|L_CalcNoxFromNO2=",L_CalcNoxFromNO2,
C    &          "|L_CalcNOXFromNO2=",L_CalcNOXFromNO2,
C    &          "|NOXMISS=", NOXMISS,
C    &          "|NOXCONC_BG=",NOXCONC_BG,
C    &          "|NO2CONC_BG=",NO2CONC_BG,
C    &          "|NOCONC_BG=",NOCONC_BG,
C    &          "|O3CONC_BG=",O3CONC_BG,
C    &          "VOCCONC_BG=",VOCCONC_BG,
C    &          "|ROCCONC_BG=",ROCCONC_BG,
C    &          "|RPCONC_BG=",RPCONC_BG
     
      RETURN
      END SUBROUTINE GETBGDEQUIL_GRS7

      SUBROUTINE GETBGDEQUIL_GRSM
C***********************************************************************
C        GETBDGEQUIL Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE: Calculate the equilibrium background for GRS chemistry
C
C        PROGRAMMER: CERC
C
C        DATE: November 2020
C
C        INPUTS:  
C
C        OUTPUTS:
C
C        CALLED FROM: GRSM_CALC
C***********************************************************************
C     Variable Declarations
      USE MAIN1
      USE GRS7MOD
      IMPLICIT NONE
      CHARACTER MODNAM*12

C --- Variable Initializations
      MODNAM = 'GETBGDEQUIL_GRSM'  
C --- Check for day or night
      IF (.NOT.L_NIGHTHOUR) THEN      
        IF(L_CalcNOxFromNO2.OR.NOXMISS)THEN
C ---     Calculate NOx bgd from NO2
C ---     Solve equilibrium eqns for NOx
          IF(O3CONC_BG/=0.0D0)THEN
            NOXCONC_BG=NO2CONC_BG/(O3CONC_BG)*(O3CONC_BG+R2/R1)
            NOCONC_BG=NOXCONC_BG-NO2CONC_BG
          END IF
        ELSE
!Calculate NO2 background from NOX values
C ---     If the concentration of NOx < NO2, then complain
          IF (NOXCONC_BG<NO2CONC_BG) THEN
            CALL ERRHDL(PATH, MODNAM, 'W','613','GRSM') 
            NOXCONC_BG = NO2CONC_BG 
          END IF
          NOCONC_BG = NOXCONC_BG - NO2CONC_BG
c         IF(INLNDEBUG)write(*,*)"HTdbg: GETBGDEQUIL_GRSM Inputs: ",
c    &          "|NOXCONC_BG=",NOXCONC_BG,
c    &          "|NO2CONC_BG=",NO2CONC_BG,
c    &          "|NOCONC_BG=",NOCONC_BG
          CALL QuadraticEquil_GRSM(NO2CONC_BG,NOCONC_BG,O3CONC_BG)
c         IF(INLNDEBUG)write(*,*)"HTdbg: GETBGDEQUIL_GRSM Outputs: ",
c    &          "|NOXCONC_BG=",NOXCONC_BG,
c    &          "|NO2CONC_BG=",NO2CONC_BG,
c    &          "|NOCONC_BG=",NOCONC_BG
        END IF
      ELSE ! Nighttime
        IF(L_CalcNOxFromNO2.OR.NOXMISS)THEN
          !NOx is missing but don't need it
          NOXCONC_BG=-999.0D0
          NOCONC_BG=0.0D0
        ELSE  
          ! No adjustments made to background at night
          NOCONC_BG = NOXCONC_BG - NO2CONC_BG
        END IF
      END IF
          
      RETURN
      END SUBROUTINE GETBGDEQUIL_GRSM      


      SUBROUTINE QuadraticEquil_GRS7(NO2,NO,O3,ROC,RP)
C***********************************************************************
C        QuadraticEquil_GRS7 Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE: Puts NO2, NO and O3 concentrations into photo-stationary
C                 equilibrium, conserving mass of NOx and NO2 + O3
C
C        PROGRAMMER: CERC
C
C        DATE: November 2020
C
C        INPUTS:  NO2 - Input NO2 concentration in volumetric units
C                 NO - Input NO concentraiton in volumetric units
C                 O3 - Input O3 concentration involumetric units
C
C        OUTPUTS: NO2 - Equilibrium NO2 concentration in the same volumetric units
C                 NO - Equilibrium NO concentraiton in the same volumetric units
C                 O3 - Equilibrium O3 concentration in the same volumetric units
C
C        CALLED FROM: GRS7_CALC, GETBGDEQUIL_GRS7
C***********************************************************************
      USE GRS7MOD
      IMPLICIT NONE
C --- Arguments      
      DOUBLE PRECISION, INTENT(INOUT)::NO2,NO,O3,ROC,RP
C --- Local variables
      DOUBLE PRECISION :: A, B, C, NOX, NOXplusO3

C --- Example for comparing reaction rates
C     Daytime QSW: 1572.9153279943321 ; R1:    1.7097255653683020E-005 ; R2:    12345.070472985180      
C                                       R3:    1.2434503463645231E-002 ; R4:    4.4471245115039996E-004 
C                                       R5:    10200.000000000000      ; R6:    125.00000000000000      ; R7:    125.00000000000000
C     Nighttime: QSW:  0.0000000      ; R1:    0.0000000000000000      ; R2:    12496.831943543026      
C                                     ; R3:    0.0000000000000000      ; R4:    4.1499147421085751E-004 
C                                     ; R5:    10200.000000000000      ; R6:    125.00000000000000      ; R7:    125.00000000000000   
C     The example show that R2 and R5 are much faster than R6+R7
C     To solve equilibrium of NOx and RP, we have to assume R6 = R7 = 0 because they led to termination of RP and NO2
C     We have 4 unknown and 4 equations:
C
C     1. NO @ equilibrium; balanced by R2, R3 and R4
C     R3[NO2] - R2[RP][NO] - R4[NO][O3] = 0 --> [NO] = R3[NO2]/(R2[RP]+R4[O3])
C     [NO] = [NO2](R3-R2[RP])/R4[O3]
C
C     2. RP @ equilibrium, balanced by R1; R2 and R5
C     R1[ROC] - R2[RP][NO] - R5[RP] = 0     --> R1[ROC] -R2[RP]{R3[NO2]/(R2[RP]+R4[O3])} -R5[RP]
C     R1R4[ROC][O3] + R1R2[ROC][RP] -R2R3[RP][NO2] - R2R5[RP]^2 - R4R5[RP][O3] = 0
C     R2R5[RP]^2 + [RP]{R2R3[NO2] + R4R5[O3] - R1R2[ROC]} - R1R4[ROC][O3] = 0
C     A = R2R5
C     B = R2R3[NO2] + R4R5[O3] - R1R2[ROC]
C     C = - R1R4[ROC][O3]
C     [RP]  = {-B +- SQRT(B^-4AC)}/2A
C     [NO]  = {R1[ROC] - R5[RP]}/R2[RP]
C     [NO2] = {R2[RP][NO]+R4[O3]}/R3 
C
C     2. Conservation of NOx
C      [NO2] + [NO] = [NO2]initial + [NO]initial = [NOX]
C
C     3. Conservation of [NOX] + [O3]
C      [NOX] + [O3] = [NOX] + [O3]

      NOX = NO2 + NO
      NOXplusO3 = NOX + O3

C     Quadratic in final A*[RP]^2 + B*[RP] + C = 0      
C     where:
C     A = R2 * R5
C     B = (R2*R3*NO2) + (R4*R5*O3) - (R1*R2*ROC)
C     C = -(R1*R4*ROC*O3)
      A = R1 * R2 * R5
      B = (R1*R3*R4*R5*(NOXplusO3-NOX)) +(R2*R3*NOX)-(R1*R2*ROC)
      C = -(R1*R3*R4*ROC*(NOXplusO3-NOX))

C --- Take negative root, as positive root will be bigger than NO2plusO3 and NOX
C     RP = 0.5*(-B - SQRT(B*B - 4.*C))/A
      RP = 0.5*(-B + SQRT(B*B - 4.*C))/A

C --- Use equations 2 and 3 above to calculate [NO] and [O3]
C     NO  = (R1*ROC - R5*RP)/(R2*RP)
C     NO  = (R3*NO2)/(R4*O3 + R2*RP)
      NO2 = NOX*((R2*RP)+(R4*O3))/(1+R3)
      NO  = NOX - NO2
C     NO  = (NOX*R3)/((R2*RP) + (R4*(NOXplusO3 - NOX)) + R3)
C     NO2 = NOX - NO
      O3  = NOXplusO3 - NO2 - NO
C     O3 = NO2plusO3 - NO2
      
      RETURN
      END SUBROUTINE QuadraticEquil_GRS7

      SUBROUTINE QuadraticEquil_GRSM(NO2,NO,O3)
C***********************************************************************
C        QuadraticEquil Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE: Puts NO2, NO and O3 concentrations into photo-stationary
C                 equilibrium, conserving mass of NOx and NO2 + O3
C
C        PROGRAMMER: CERC
C
C        DATE: November 2020
C
C        INPUTS:  NO2 - Input NO2 concentration in volumetric units
C                 NO - Input NO concentraiton in volumetric units
C                 O3 - Input O3 concentration involumetric units
C
C        OUTPUTS: NO2 - Equilibrium NO2 concentration in the same volumetric units
C                 NO - Equilibrium NO concentraiton in the same volumetric units
C                 O3 - Equilibrium O3 concentration in the same volumetric units
C
C        CALLED FROM: GRSM_CALC, GETBGDEQUIL
C***********************************************************************
      USE GRS7MOD
      IMPLICIT NONE
C --- Arguments      
      DOUBLE PRECISION, INTENT(INOUT)::NO2,NO,O3
C --- Local variables
      DOUBLE PRECISION :: B, C, NOX, NO2plusO3
      
C --- We have 3 equations in 3 unknowns:
C     1. Photo-stationary equilibrium
C      [NO2] = (R1/R2)*[NO]*[O3]
C     2. Conservation of NOx
C      [NO2] + [NO] = [NO2]initial + [NO]initial = [NOX]
C     3. Conservation of [NO2] + [O3]
C      [NO2] + [O3] = [NO2]initial + [O3]initial
C     These 3 lead to quadratic in [NO2]

C --- Set conserved quantities:
      NOX = NO2 + NO
      NO2plusO3 = NO2 + O3

C --- Quadratic in final [NO2]:  [NO2]*[NO2] + B*[NO2] + C = 0
C     where:
      B = -(NOX + NO2plusO3 + (R4/R3))
      C = NOX*NO2plusO3

C --- Take negative root, as positive root will be bigger than NO2plusO3 and NOX
      NO2 = 0.5*(-B - SQRT(B*B - 4.*C))

C --- Use equations 2 and 3 above to calculate [NO] and [O3]
      NO = NOX - NO2
      O3 = NO2plusO3 - NO2
      
      RETURN
      END SUBROUTINE QuadraticEquil_GRSM      

      FUNCTION FindInitDt_GRS7()
C***********************************************************************
C        FindInitDt_GRS7 Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE: Find initial GRSM chemistry time step
C
C        PROGRAMMER: CERC
C
C        DATE: November 2020
C
C        INPUTS:  
C
C        OUTPUTS:
C
C        CALLED FROM: DoGRS7Chem
C***********************************************************************
      USE GRS7MOD
      USE MAIN1, ONLY: DEBUG
      IMPLICIT NONE
C --- Local variables
      INTEGER::ConcLoop
      DOUBLE PRECISION:: dCdt(nPolsGRS7), TimeScale, TimeLoop 
      DOUBLE PRECISION:: FindInitDt_GRS7
      CHARACTER::MODNAM*12

C --- Initialisations
      MODNAM='FindInitDt_GRS7'
      
C --- Calculate the first order dreivatives
c     IF(INLNDEBUG)write(*,*)"HTdbg: GRS7OPT=",GRS7OPT
      IF(GRS7OPT==1)CALL DConcsDt_GRS7(CONCTEMP,dCdt)
      IF(GRS7OPT==2)CALL DConcsDt_ADM(CONCTEMP,dCdt)

C --- Loop through each non zero concentration and estimate the
C     'time scale' for each pollutant. We want to find the minimum non-
C     zero time period so that all pollutants change by no more than
C     CFrac*Initial Concentration.

      TimeScale = 1. ! initial value
      DO ConcLoop = 1, nPolsGRS7
c       IF(INLNDEBUG)write(*,*)"HTdbg: FindInitDt_GRS7 in ",
c    &            " | Pols=",GRS7PNAME(ConcLoop),
c    &            " | ConcTemp",ConcTemp(ConcLoop ),
c    &            " | dCdt",dCdt(ConcLoop)
        IF ( ConcTemp(ConcLoop ).GT.0.) THEN
          IF(ABS(dCdt(ConcLoop))>CFrac*ConcTemp(ConcLoop))THEN
            TimeLoop = CFrac*ABS(ConcTemp(ConcLoop)/dCdt(ConcLoop))
            IF ( TimeLoop.LT.TimeScale ) TimeScale = TimeLoop
          END IF
        END IF
c       write(*,*)"HTdbg: FindInitDt_GRS7 out ",
c    &            " | nPols=",ConcLoop,
c    &            " | ConcTemp",ConcTemp(ConcLoop ),
c    &            " | dCdt",dCdt(ConcLoop)
      END DO

      FindInitDt_GRS7 = TimeScale

      RETURN
      END FUNCTION FindInitDt_GRS7

      SUBROUTINE DConcsDt_ADM(CONC, dCdt)
C***********************************************************************
C        DConcsDt_GRS7 Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE: Calculate first-order time derivatives of each chemistry
C                 pollutant from the reaction equations
C
C        PROGRAMMER: CERC
C
C        DATE: November 2020
C
C        INPUTS:  CONC - input concentrations
C
C        OUTPUTS: dCdt - first-order time derivatives
C
C        CALLED FROM: DoGRS7Chem, FindInitDt_GRS7, CashKarpRK_GRS7
C***********************************************************************
      USE GRS7MOD
      IMPLICIT NONE
C --- Arguments
      DOUBLE PRECISION, INTENT(OUT):: dCdt(nPolsGRS7)
      DOUBLE PRECISION, INTENT(INOUT):: CONC(nPolsGRS7)
      
C --- Local variables
      CHARACTER::MODNAM*12
CRCO 3/4/2021 Removing unused variable
C      INTEGER::I
     
C --- Adopted from ADM code
      DOUBLE PRECISION :: bquad, root2 
 
C --- Initialisations
      MODNAM='DConcsDt_ADM'

C EVALUATE RATE FUNCTIONS

C calculate concentrations of RP
 
      IF(.NOT.L_NightHour) THEN
        bquad = R2*CONC(nNO) + R6*CONC(nNO2) + R7*CONC(nNO2)
        root2 = bquad*bquad+(4.*R5*R1*CONC(nROC))
        IF(root2 > 0.0D0)THEN
            CONC(nRP) = -bquad+sqrt(root2)
            CONC(nRP) = CONC(nRP)*0.5/R5
        ELSE
            CONC(nRP) = 0.0D0
        END IF
      ELSE
        CONC(nRP) = 0.0D0
      END IF

C --- Calculate VALUE OF dC/dt

C --- NO2
      dCdt(nNO2)= R2*CONC(nNO)*CONC(nRP) 
     &          + R4*CONC(nNO)*CONC(nO3)
     &          - R3*CONC(nNO2) 
     &          - (R6+R7)*CONC(nRP)*CONC(nNO2)
      !dCdt(nNO2)=R4*CONC(nNO)*CONC(nO3)-R3*CONC(nNO2)+R2*CONC(nNO)

C --- NO
      dCdt(nNO) = R3*CONC(nNO2)
     &          - R4*CONC(nNO)*CONC(nO3) 
     &          - R2*CONC(nNO)*CONC(nRP)

C --- O3
      dCdt(nO3) = R3*CONC(nNO2)
     &          - R4*CONC(nNO)*CONC(nO3)

C --- ROC
      dCdt(nROC)=0.0D0

C --- RP
      dCdt(nRP) = R1*CONC(nROC)
     &          - R2*CONC(nRP)*CONC(nNO)
     &          - R5*(CONC(nRP)*CONC(nRP))
     &          - CONC(nNO2)*CONC(nRP)*(R6+R7)
c     dCdt(nRP)=0.

C --- SGN
      dCdt(nSGN)=R6*CONC(nRP)*CONC(nNO2)
C --- SNGN
      dCdt(nSNGN)=R6*CONC(nRP)*CONC(nNO2)

      RETURN
      END SUBROUTINE DConcsDt_ADM

      SUBROUTINE DConcsDt_GRS7(CONC, dCdt)
C***********************************************************************
C        DConcsDt_GRS7 Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE: Calculate first-order time derivatives of each chemistry
C                 pollutant from the reaction equations
C
C        PROGRAMMER: CERC
C
C        DATE: November 2020
C
C        INPUTS:  CONC - input concentrations
C
C        OUTPUTS: dCdt - first-order time derivatives
C
C        CALLED FROM: DoGRS7Chem, FindInitDt_GRS7, CashKarpRK_GRS7
C***********************************************************************
      USE GRS7MOD
      IMPLICIT NONE
C --- Arguments
      DOUBLE PRECISION, INTENT(OUT):: dCdt(nPolsGRS7)
      DOUBLE PRECISION, INTENT(INOUT):: CONC(nPolsGRS7)
      
C --- Local variables
      CHARACTER::MODNAM*12
CRCO 3/4/2021 Removing unused variable
C      INTEGER::I
     
C --- Initialisations
      MODNAM='DConcsDt_GRS7'

C --- Calculate VALUE OF dC/dt

C --- NO2
      dCdt(nNO2)= R2*CONC(nNO)*CONC(nRP) 
     &          + R4*CONC(nNO)*CONC(nO3)
     &          - R3*CONC(nNO2) 
     &          - (R6+R7)*CONC(nRP)*CONC(nNO2)
      !dCdt(nNO2)=R4*CONC(nNO)*CONC(nO3)-R3*CONC(nNO2)+R2*CONC(nNO)

C --- NO
      dCdt(nNO) = R3*CONC(nNO2) 
     &          - R4*CONC(nNO)*CONC(nO3) 
     &          - R2*CONC(nNO)*CONC(nRP)

C --- O3
      dCdt(nO3) = R3*CONC(nNO2)
     &          - R4*CONC(nNO)*CONC(nO3)

C --- ROC
      dCdt(nROC)=0.

C --- RP
      dCdt(nRP) = R1*CONC(nROC) 
     &          - R2*CONC(nRP)*CONC(nNO)
     &          - R5*CONC(nRP)*CONC(nRP) 
     &          - CONC(nNO2)*CONC(nRP)*(R6+R7)
c     dCdt(nRP)=0.

C --- SGN
      dCdt(nSGN)=R6*CONC(nRP)*CONC(nNO2)

C --- SNGN
      dCdt(nSNGN)=R6*CONC(nRP)*CONC(nNO2)

      RETURN
      END SUBROUTINE DConcsDt_GRS7


      SUBROUTINE DoGRS7Chem
C***********************************************************************
C        DoGRS7Chem Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE: Calculate the final concentrations after applying
C                 GRS7 chemistry scheme
C
C        PROGRAMMER: CERC
C
C        DATE: November 2020
C
C        INPUTS:  
C
C        OUTPUTS:
C
C        CALLED FROM: GRS7_CALC
C***********************************************************************
      USE MAIN1
      USE GRS7MOD
      IMPLICIT NONE
       
C --- Local variables
      CHARACTER::MODNAM*12
      INTEGER, PARAMETER::MAXSTP=1000000
      INTEGER:: i,nstp
      DOUBLE PRECISION:: dt, dtnext, tLocal, 
     &         dCdt(nPolsGRS7),Cscal(nPolsGRS7)
      DOUBLE PRECISION, PARAMETER::dtmin=1.0D-6
C --- functions
      DOUBLE PRECISION::FindInitDt_GRS7
!     
C --- Initialisations
      MODNAM='DoGRS7Chem'
      tLocal=0.D0
      
C --- Check for negative concentrations (should never start off
C     negative but check for safety)
      DO  I=1,nPolsGRS7
        IF (CONCTEMP(I).LT.0.) CONCTEMP(I) = 0.
      END DO

c     CONCTEMP = CONCTEMP * 1.0D-03 ! Convert ppb to ppm
      
C --- Estimate initial time step
      dt = FindInitDt_GRS7()
c     IF(INLNDEBUG)write(*,*)"HTdbg: FindInitDt_GRS7 init dt = ",dt

C --- Increment in timesteps until end time reached or we hit max No. of steps
      DO nstp=1,MAXSTP

C       Calculate the derivatives
        IF(GRS7OPT==1)CALL DConcsDt_GRS7(CONCTEMP,dCdt)
        IF(GRS7OPT==2)CALL DConcsDt_ADM(CONCTEMP,dCdt)

        DO i=1,nPolsGRS7
          Cscal(i)=CONCTEMP(i)+abs(dt*dCdt(i))
c         IF(INLNDEBUG.AND.IREC==3)write(*,*)"HTdbg: DoGRS7Chem ",
c    &            "nstp=",nstp,     
c    &            "dt=",dt,     
c    &            " | Pols=",GRS7PNAME(i),
c    &            " | ConcTemp=",ConcTemp(i),
c    &            " | dCdt=",dCdt(i),
c    &            " | Cscal=",Cscal(i)
        END DO
        
C       Ensure we don't go past the travel time
        IF(tLocal+dt.gt.TTRAVCHM(IREC)) dt=TTRAVCHM(IREC)-tLocal
        
C       Perform Runge-Kutta integration using adaptive timestep
        CALL AdaptiveRK_GRS7(CONCTEMP,dCdt,tLocal,dt,Cscal,dtnext)
        
C       Check for negative concentrations
        DO  I=1,nPolsGRS7
c         IF(INLNDEBUG)write(*,*)"HTdbg: DoGRS7Chem ",
c    &            " | Pols=",GRS7PNAME(i),
c    &            " | ConcTemp=",ConcTemp(i)
          IF (CONCTEMP(I).LT.0.) CONCTEMP(I) = 0.
        END DO

C       If we've reached the travel time, return
        IF(tLocal.ge.TTRAVCHM(IREC)) THEN
c           CONCTEMP = CONCTEMP * 1.0D+03 ! Convert ppm to ppb
            RETURN
        ENDIF
C       Otherwise set the timestep for the next iteration
        IF(dtnext.lt.dtmin) dtnext=dtmin
        dt=dtnext
      END DO
C     Reached maximum No. of time steps - return with latest concentrations
c     CONCTEMP = CONCTEMP * 1.0D+03 ! Convert ppm to ppb
      RETURN
      END SUBROUTINE DoGRS7Chem

      SUBROUTINE AdaptiveRK_GRS7(CLocal,dCdt,tLocal,dttry,Cscal,dtnext)
C***********************************************************************
C        AdaptiveRK_GRS7 Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE: Performs Runge-Kutta integration using input timestep fisrt.
C                 If the error is too high, reduces the timestep and
C                 tries again until error is below threshold. Also tries
C                 to set a sensible timestep for the next time this
C                 routine is called based on the magnitude of the error
C
C        PROGRAMMER: CERC
C
C        DATE: November 2020
C
C        INPUTS:  CLocal - Input concentrations
C                 dCdt - Time derivatives of input concentrations
C                 tLocal - Total time so far
C                 dttry - First timestep to try
C                 Cscal - CLocal(i)+abs(dttry*dCdt(i))
C
C        OUTPUTS: tLocal - Total time so far (updated)
C                 CLocal - Output concentrations (updated)
C                 dtnext - Timestep to use next time around
C
C        CALLED FROM: DoGRS7Chem
C***********************************************************************
      USE MAIN1
      USE GRS7MOD
      IMPLICIT NONE
!      
C --- Arguments
      DOUBLE PRECISION, INTENT(INOUT)::CLocal(nPolsGRS7), tLocal
      DOUBLE PRECISION, INTENT(IN)::dttry, dCdt(nPolsGRS7),
     &                              Cscal(nPolsGRS7)
      DOUBLE PRECISION, INTENT(OUT)::dtnext

C --- Local variables
      CHARACTER::MODNAM*12
      INTEGER::I
      DOUBLE PRECISION::dt, CTemp(nPolsGRS7), GrowthFac,
     &               CErr(nPolsGRS7), ErrMax, ttemp, tNew
      !P_Shrink - Determines how much we shrink dt by for the next 
      !iteration if the error is too high
      DOUBLE PRECISION, PARAMETER::P_Shrink=-0.25D0
      !P_Grow - Determines how much we increase dt by for the next
      !call to this routine
      DOUBLE PRECISION, PARAMETER::P_Grow=-0.2D0
      !MaxGrow - Maximum growth factor for dt
      DOUBLE PRECISION, PARAMETER::MaxGrow=5.D0
      !SAFETY - reduces dt slightly for the next iteration
      DOUBLE PRECISION, PARAMETER::SAFETY=0.9D0 
          
C --- Initialisations
      MODNAM = 'AdaptiveRK_GRS7'
      dt=dttry
1     CALL CashKarpRK_GRS7(CLocal,dCdt,dt,CTemp,CErr)
      errmax=0.

      DO i=1,nPolsGRS7
        IF (Cscal(i).gt.0.)THEN
          errmax=max(errmax,abs(CErr(i)/Cscal(i)))
        END IF
      END DO 

      errmax=errmax/CFrac
      IF(errmax.gt.1.)THEN
        !Error too large - shrink the timestep and try again
        ttemp=SAFETY*dt*(errmax**P_Shrink)
        dt=max(ttemp,0.1*dt)
        tNew=tLocal+dt
        GOTO 1
      ELSE
        !Error below threshold - calculate initial timestep for
        !next call to this routine and return
        GrowthFac = SAFETY*(errmax**P_Grow)
        IF(GrowthFac.gt.MaxGrow)THEN
          dtnext=MaxGrow*dt
        ELSE
          dtnext=GrowthFac*dt
        END IF
        tLocal=tLocal+dt
        DO i=1,nPolsGRS7
          CLocal(i)=CTemp(i)
        END DO 
        RETURN
      END IF
      END SUBROUTINE AdaptiveRK_GRS7 

      SUBROUTINE CashKarpRK_GRS7(C,dCdt,dt,COut,CErr)
C***********************************************************************
C        CashKarpRK_GRS7 Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE: Performs Runge-Kutta integration (Cash-Karp method) once
C                 using input timestep. Also calculates error.
C
C        PROGRAMMER: CERC
C
C        DATE: November 2020
C
C        INPUTS:  C - Input concentrations
C                 dCdt - Time derivatives of input concentrations
C                 dt - Timestep to use
C
C        OUTPUTS: COut - Output concentrations
C                 CErr - Error per pollutant
C
C        CALLED FROM: AdaptiveRK_GRS7
C***********************************************************************
      USE GRS7MOD
      IMPLICIT NONE

C --- Arguments      
      DOUBLE PRECISION, INTENT(IN)::dt,dCdt(nPolsGRS7),
     &                              C(nPolsGRS7)
      DOUBLE PRECISION, INTENT(OUT)::CErr(nPolsGRS7),
     &                               COut(nPolsGRS7)

C --- Local variables
      CHARACTER::MODNAM*12
      INTEGER::i
      DOUBLE PRECISION:: ak2(nPolsGRS7),ak3(nPolsGRS7),
     &                   ak4(nPolsGRS7),ak5(nPolsGRS7),
     &                   ak6(nPolsGRS7),CTemp(nPolsGRS7)
CRCO 3/4/2021 removing unused variables
C      DOUBLE PRECISION::A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51,
      DOUBLE PRECISION::B21,B31,B32,B41,B42,B43,B51,
     &                  B52,B53,B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,
     &                  DC1,DC3,DC4,DC5,DC6
C      PARAMETER (A2=.2D0,A3=.3D0,A4=.6D0,A5=1.D0,A6=.875D0,B21=.2D0,
      PARAMETER (B21=.2D0,
     &B31=3.D0/40.D0,B32=9.D0/40.D0,B41=.3D0,B42=-.9D0,B43=1.2D0,
     &B51=-11.D0/54.D0,B52=2.5D0,B53=-70.D0/27.D0,B54=35.D0/27.D0,
     &B61=1631.D0/55296.D0,B62=175.D0/512.D0,B63=575.D0/13824.D0,
     &B64=44275.D0/110592.D0,B65=253.D0/4096.D0,C1=37.D0/378.D0,
     &C3=250.D0/621.D0,C4=125.D0/594.D0,C6=512.D0/1771.D0,
     &DC1=C1-2825.D0/27648.D0,DC3=C3-18575.D0/48384.D0,
     &DC4=C4-13525.D0/55296.D0,DC5=-277.D0/14336.D0,DC6=C6-.25D0)
      
C --- Initialisations
      MODNAM = 'CashKarpRK_GRS7'
      
      DO i=1,nPolsGRS7
        CTemp(i)=C(i)+B21*dt*dCdt(i)
        IF(CTemp(i).LT.0.) CTemp(i) = 0. !Prevent negative concs
      END DO
      IF(GRS7OPT==1)CALL DConcsDt_GRS7(CTemp,ak2)
      IF(GRS7OPT==2)CALL DConcsDt_ADM(CTemp,ak2)
      DO i=1,nPolsGRS7
        CTemp(i)=C(i)+dt*(B31*dCdt(i)+B32*ak2(i))
        IF(CTemp(i).LT.0.) CTemp(i) = 0. !Prevent negative concs
      END DO
      IF(GRS7OPT==1)CALL DConcsDt_GRS7(CTemp,ak3)
      IF(GRS7OPT==2)CALL DConcsDt_ADM(CTemp,ak3)
      DO i=1,nPolsGRS7
        CTemp(i)=C(i)+dt*(B41*dCdt(i)+B42*ak2(i)+B43*ak3(i))
        IF(CTemp(i).LT.0.) CTemp(i) = 0. !Prevent negative concs
      END DO
      IF(GRS7OPT==1)CALL DConcsDt_GRS7(CTemp,ak4)
      IF(GRS7OPT==2)CALL DConcsDt_ADM(CTemp,ak4)
      DO i=1,nPolsGRS7
        CTemp(i)=C(i)+dt*(B51*dCdt(i)+B52*ak2(i)+B53*ak3(i)+B54*ak4(i))
        IF(CTemp(i).LT.0.) CTemp(i) = 0. !Prevent negative concs
      END DO
      IF(GRS7OPT==1)CALL DConcsDt_GRS7(CTemp,ak5)
      IF(GRS7OPT==2)CALL DConcsDt_ADM(CTemp,ak5)
      DO i=1,nPolsGRS7
        CTemp(i)=C(i)+dt*(B61*dCdt(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+
     &           B65*ak5(i))
        IF(CTemp(i).LT.0.) CTemp(i) = 0. !Prevent negative concs
      END DO
      IF(GRS7OPT==1)CALL DConcsDt_GRS7(CTemp,ak6)
      IF(GRS7OPT==2)CALL DConcsDt_ADM(CTemp,ak6)
      DO i=1,nPolsGRS7
        COut(i)=C(i)+dt*(C1*dCdt(i)+C3*ak3(i)+C4*ak4(i)+C6*ak6(i))
      END DO
      DO i=1,nPolsGRS7
        CErr(i)=dt*(DC1*dCdt(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6*
     &          ak6(i))
      END DO
      RETURN
      END SUBROUTINE CashKarpRK_GRS7

      SUBROUTINE ROCSTK
C***********************************************************************
C        ROCSTK Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE: Processes VOC to ROC Default In-stack Ratio Value for GRS7
C                 based on the ROCSTACK keyword
C
C        PROGRAMMER: Huy Tran - UNC, based on NO2STACK (coset.f) by Roger W. Brode
C
C        DATE:     September 24, 2023
C        MODIFIED: November 28, 2023
C                  Since NO2STK in aermod_23123 is identical to aermod_22112,
C                  no update made to ROCSTK which was based on NO2STK
C
C        INPUTS:  Input Runstream Image Parameters
C
C        OUTPUTS:
C
C        CALLED FROM:   COCARD
C***********************************************************************
          
C     Variable Declarations
         USE MAIN1
         IMPLICIT NONE
         INTEGER I
         CHARACTER MODNAM*12
          
C        Variable Initializations
         MODNAM = 'ROCSTK'
          
C        Check The Number Of The Fields
         IF (IFC .LE. 2) THEN
C           Error Message: No Parameters; same error code 200
            CALL ERRHDL(PATH,MODNAM,'E','200',KEYWRD)
            GO TO 999
         ELSE IF (IFC .GT. 3) THEN
C        Error Message: Too Many Parameters; same error code 202
            CALL ERRHDL(PATH,MODNAM,'E','202',KEYWRD)
            GO TO 999
         END IF
          
C        Gets Double Precision of Real Number from a character string
         CALL STODBL(FIELD(3),ILEN_FLD,DNUM,IMIT)
C        Check The Numerical Field; same error code 208
         IF (IMIT .NE. 1) THEN
            CALL ERRHDL(PATH,MODNAM,'E','208',KEYWRD)
            GO TO 999
         END IF
          
C        Assign value to ROCStack variable
         ROCStack = DNUM
C        write(*,*)"HTdbg: sub ROCSTK: ROCStack = ", ROCStack
          
C        Check range of value; same error code 380
         IF (ROCStack .LT. 0.0D0) THEN
            CALL ERRHDL(PATH,MODNAM,'E','380','ROCStack')
            GO TO 999
         END IF
          
         DO I = 1, NSRC
            AROC_RATIO(I) = ROCStack
         END DO
          
 999     RETURN
      END SUBROUTINE ROCSTK

      SUBROUTINE SEMVARY
C***********************************************************************
C         EMVARY Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE: Processes Secondary Variable Emission Rate Factors for VOC
C
C        PROGRAMMER:  Huy Tran, based on EMVARY (soset.f) by Jeff Wang, Roger Brode
C
C        DATE:    September 24, 2023
C
C        MODIFIED: November 28, 2023
C                  Since EMVARY in aermod_23123 is identical in aermod_22112,
C                  no update made to SEMVARY which was based on EMVARY
C
C
C        INPUTS:  Input Runstream Image Parameters
C
C        OUTPUTS: Variable Emmission Rate Factors
C
C        CALLED FROM:   SOCARD
C***********************************************************************

C     Variable Declarations
      USE MAIN1
      IMPLICIT NONE
      CHARACTER MODNAM*12

      INTEGER :: I, IH, IL, ISDX, IQMAX
      CHARACTER (LEN=12) :: LID, HID, LID1, LID2, HID1, HID2
      CHARACTER (LEN=ILEN_FLD) :: SOID
      LOGICAL FOUND, INGRP, RMARK

C     Variable Initializations
      FOUND  = .FALSE.
      INGRP  = .FALSE.
      MODNAM = 'SEMVARY'

C     Check The Number Of The Fields
      IF (IFC .LE. 2) THEN
C        Error Message: No Parameters
         CALL ERRHDL(PATH,MODNAM,'E','200',KEYWRD)
         GO TO 999
      ELSE IF (IFC .EQ. 3) THEN
C        Error Message: No Numerical Parameters
         CALL ERRHDL(PATH,MODNAM,'E','201',KEYWRD)
         GO TO 999
      ELSE IF (IFC .LT. 5) THEN
C        Error Message: Not Enough Parameters
         CALL ERRHDL(PATH,MODNAM,'E','201',KEYWRD)
         GO TO 999
      END IF

C     Get The Source ID(s)
      SOID = FIELD(3)
      CALL FSPLIT(PATH,KEYWRD,SOID,ILEN_FLD,'-',RMARK,LID,HID)

C     Verify The Effective Srcid
      IF (LID .EQ. HID) THEN
C        Search For The Index
         CALL SINDEX(SRCID,NSRC,SOID,ISDX,FOUND)
         IF (FOUND) THEN
            SQFLAG(ISDX) = FIELD(4)
            IF (SQFLAG(ISDX) .EQ. 'SEASON') THEN
               IQMAX = 4
            ELSE IF (SQFLAG(ISDX) .EQ. 'MONTH') THEN
               IQMAX = 12
            ELSE IF (SQFLAG(ISDX) .EQ. 'HROFDY') THEN
               IQMAX = 24
            ELSE IF (SQFLAG(ISDX) .EQ. 'WSPEED') THEN
               IQMAX = 6
            ELSE IF (SQFLAG(ISDX) .EQ. 'SEASHR') THEN
               IQMAX = 96
            ELSE IF (SQFLAG(ISDX) .EQ. 'HRDOW') THEN
               IQMAX = 72
               L_DayOfWeekOpts = .TRUE.
            ELSE IF (SQFLAG(ISDX) .EQ. 'HRDOW7') THEN
               IQMAX = 168
               L_DayOfWeekOpts = .TRUE.
            ELSE IF (SQFLAG(ISDX) .EQ. 'SHRDOW') THEN
               IQMAX = 288
               L_DayOfWeekOpts = .TRUE.
            ELSE IF (SQFLAG(ISDX) .EQ. 'SHRDOW7') THEN
               IQMAX = 672
               L_DayOfWeekOpts = .TRUE.
            ELSE IF (SQFLAG(ISDX) .EQ. 'MHRDOW') THEN
               IQMAX = 864
               L_DayOfWeekOpts = .TRUE.
            ELSE IF (SQFLAG(ISDX) .EQ. 'MHRDOW7') THEN
               IQMAX = 2016
               L_DayOfWeekOpts = .TRUE.
            ELSE
C              WRITE Error Message    ! Invalid SQFLAG Field Entered
               CALL ERRHDL(PATH,MODNAM,'E','203','SQFLAG')
            END IF
            IF (IQMAX .LE. SNQF) THEN
               CALL SEFFILL(ISDX,IQMAX)
            ELSE
C              WRITE Error Message     ! SNQF Parameter Not Large Enough
               WRITE(DUMMY,'(''SNQF ='',I6)') SNQF
               CALL ERRHDL(PATH,MODNAM,'E','260',DUMMY)
            END IF
         ELSE
C           WRITE Error Message     ! Source Location Has Not Been Identified
            CALL ERRHDL(PATH,MODNAM,'E','300',KEYWRD)
         END IF
      ELSE
C        First Check Range for Upper Value < Lower Value
         CALL SETIDG(LID,LID1,IL,LID2)
         CALL SETIDG(HID,HID1,IH,HID2)
         IF ((HID1.LT.LID1) .OR. (IH.LT.IL) .OR. (HID2.LT.LID2)) THEN
C           WRITE Error Message:  Invalid Range,  Upper < Lower
            CALL ERRHDL(PATH,MODNAM,'E','203','SRCRANGE')
            GO TO 999
         END IF
         DO I = 1, NUMSRC
C           See Whether It's In The Group
            CALL ASNGRP(SRCID(I),LID,HID,INGRP)
            IF (INGRP) THEN
               ISDX = I
               SQFLAG(ISDX) = FIELD(4)
               IF (SQFLAG(ISDX) .EQ. 'SEASON') THEN
                  IQMAX = 4
               ELSE IF (SQFLAG(ISDX) .EQ. 'MONTH') THEN
                  IQMAX = 12
               ELSE IF (SQFLAG(ISDX) .EQ. 'HROFDY') THEN
                  IQMAX = 24
               ELSE IF (SQFLAG(ISDX) .EQ. 'WSPEED') THEN
                  IQMAX = 6
               ELSE IF (SQFLAG(ISDX) .EQ. 'SEASHR') THEN
                  IQMAX = 96
               ELSE IF (SQFLAG(ISDX) .EQ. 'HRDOW') THEN
                  IQMAX = 72
                  L_DayOfWeekOpts = .TRUE.
               ELSE IF (SQFLAG(ISDX) .EQ. 'HRDOW7') THEN
                  IQMAX = 168
                  L_DayOfWeekOpts = .TRUE.
               ELSE IF (SQFLAG(ISDX) .EQ. 'SHRDOW') THEN
                  IQMAX = 288
                  L_DayOfWeekOpts = .TRUE.
               ELSE IF (SQFLAG(ISDX) .EQ. 'SHRDOW7') THEN
                  IQMAX = 672
                  L_DayOfWeekOpts = .TRUE.
               ELSE IF (SQFLAG(ISDX) .EQ. 'MHRDOW') THEN
                  IQMAX = 864
                  L_DayOfWeekOpts = .TRUE.
               ELSE IF (SQFLAG(ISDX) .EQ. 'MHRDOW7') THEN
                  IQMAX = 2016
                  L_DayOfWeekOpts = .TRUE.
               ELSE
C                 WRITE Error Message    ! Invalid SQFLAG Field Entered
                  CALL ERRHDL(PATH,MODNAM,'E','203','SQFLAG')
               END IF
               IF (IQMAX .LE. SNQF) THEN
                  CALL SEFFILL(ISDX,IQMAX)
               ELSE
C                 WRITE Error Message    ! NQF Parameter Not Large Enough
                  WRITE(DUMMY,'(''SNQF ='',I6)') SNQF
                  CALL ERRHDL(PATH,MODNAM,'E','260',DUMMY)
               END IF
            END IF
         END DO
      END IF

 999  RETURN
      END

      SUBROUTINE SEFFILL(ISDX,IQMAX)
C***********************************************************************
C        SEFFILL Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE: Fill Secondary Variable Emission Rate Array (for VOC)
C
C        PROGRAMMER:  Huy Tran, based on EFFILL (soset.f) by Roger Brode, Jeff Wang
C
C        DATE:    September 24, 2023
C
C        MODIFIED: November 28, 2023
C                  Since EFFILL in aermod_23123 is identical in aermod_22112,
C                  no update made to SEFFILL which was based on EFFILL
C
C        INPUTS:  Input Runstream Image Parameters
C
C        OUTPUTS: Direction Specific Building Directions
C
C        CALLED FROM:   SEMVARY
C***********************************************************************

C     Variable Declarations
      USE MAIN1
      IMPLICIT NONE
      CHARACTER MODNAM*12

      INTEGER :: J, K, ISDX, IQMAX

C     Variable Initializations
      MODNAM = 'SEFFILL'

      ISET = IWRK2(ISDX,4)

      DO K = 5, IFC
C        Change Fields To Numbers
         CALL STODBL(FIELD(K),ILEN_FLD,DNUM,IMIT)
C        Check The Numerical Field
         IF (IMIT .EQ. -1) THEN
            CALL ERRHDL(PATH,MODNAM,'E','208',KEYWRD)
            CYCLE
         END IF
         DO J = 1, IMIT
            ISET = ISET + 1
C           Assign The Field
            IF (ISET .LE. IQMAX) THEN
               SQFACT(ISET,ISDX) = DNUM
               IF (DNUM .LT. 0.0D0) THEN
C                 WRITE Error Message:  Negative Value for SQFACT
                  CALL ERRHDL(PATH,MODNAM,'E','209',KEYWRD)
               END IF
            ELSE
C              WRITE Error Message    ! Too Many SQFACT Values Input
               IF (ISDX .LE. 999) THEN
                  WRITE(DUMMY,'(''SQFACT Src'',I3.3)') ISDX
               ELSE IF (ISDX .LE. 99999) THEN
                  WRITE(DUMMY,'(''QF Src'',I5.5)') ISDX
               ELSE
                  WRITE(DUMMY,'(''QF Src>99999'')')
               END IF
               CALL ERRHDL(PATH,MODNAM,'E','231',DUMMY)
            END IF
         END DO
      END DO

      IWRK2(ISDX,4) = ISET

      RETURN
      END     
      

C***********************************************************************
      SUBROUTINE VOC2ROC
C***********************************************************************
C        VOC2ROC Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE: Processes In-stack VOC to ROC Ratios by Source for
C                 GRS7 and MGRS Options
C
C        PROGRAMMER: Huy Tran; based on NO2RAT (soset.f) by Roger W. Brode, PES, Inc.
C
C        DATE:    September 25, 2023
C
C        MODIFIED: November 28, 2023
C                  Since NO2RAT in aermod_23123 is identical in aermod_22112,
C                  no update made to VOC2ROC which was based on NO2RAT
C
C        INPUTS:  Input Runstream Image Parameters
C
C        OUTPUTS: Array of in-stack VOC/ROC ratios
C
C        CALLED FROM:   SOCARD
C***********************************************************************

C     Variable Declarations
      USE MAIN1
      IMPLICIT NONE
      CHARACTER MODNAM*12

      INTEGER :: I, IH, IL, ISDX
      CHARACTER (LEN=12) :: LID, HID, LID1, LID2, HID1, HID2
      CHARACTER (LEN=ILEN_FLD) :: SOID
      LOGICAL FOUND, INGRP, RMARK

C     Variable Initializations
      FOUND  = .FALSE.
      INGRP  = .FALSE.
      MODNAM = 'VOC2ROC'

C     Check the Number of Fields
      IF (IFC .LE. 2) THEN
C        Error Message: No Parameters
         CALL ERRHDL(PATH,MODNAM,'E','200',KEYWRD)
         GO TO 999
      ELSE IF (IFC .LT. 4) THEN
C        Error Message: Not Enough Parameters
         CALL ERRHDL(PATH,MODNAM,'E','201',KEYWRD)
         GO TO 999
      ELSE IF (IFC .GT. 4) THEN
C        Error Message: Too Many Parameters
         CALL ERRHDL(PATH,MODNAM,'E','202',KEYWRD)
         GO TO 999
      END IF

C     Get The Source ID(s)
      SOID = FIELD(3)
      CALL FSPLIT(PATH,KEYWRD,SOID,ILEN_FLD,'-',RMARK,LID,HID)

      IF (LID .EQ. HID) THEN
C        Search For The Index
         CALL SINDEX(SRCID,NSRC,SOID,ISDX,FOUND)
         IF (FOUND) THEN
C           Read NO2/NOX Ratio and Convert to Real
            CALL STODBL(FIELD(4),ILEN_FLD,DNUM,IMIT)
C           Check The Numerical Field
            IF (IMIT .NE. 1) THEN
               CALL ERRHDL(PATH,MODNAM,'E','208',KEYWRD)
               GO TO 999
            END IF
            IF (DNUM .LT. 0.0D0) THEN
C              WRITE Error Message: VOC2ROC Ratio Out-of-Range
               CALL ERRHDL(PATH,MODNAM,'E','804',SRCID(ISDX))
            END IF
C           Assign The Field
            AROC_RATIO(ISDX) = DNUM
         ELSE
C           WRITE Error Message     ! Source Location Has Not Been Identified
            CALL ERRHDL(PATH,MODNAM,'E','300',KEYWRD)
         END IF
      ELSE
C        First Check Range for Upper Value < Lower Value
         CALL SETIDG(LID,LID1,IL,LID2)
         CALL SETIDG(HID,HID1,IH,HID2)
         IF ((HID1.LT.LID1) .OR. (IH.LT.IL) .OR. (HID2.LT.LID2)) THEN
C           WRITE Error Message:  Invalid Range,  Upper < Lower
            CALL ERRHDL(PATH,MODNAM,'E','203','SRCRANGE')
            GO TO 999
         END IF
         DO I = 1, NUMSRC
C           See Whether It's In The Group
            CALL ASNGRP(SRCID(I),LID,HID,INGRP)
            IF (INGRP) THEN
C              Read NO2/NOX Ratio and Convert to Real
               CALL STODBL(FIELD(4),ILEN_FLD,DNUM,IMIT)
C              Check The Numerical Field
               IF (IMIT .NE. 1) THEN
                  CALL ERRHDL(PATH,MODNAM,'E','208',KEYWRD)
                  GO TO 999
               END IF
               IF (DNUM .LT. 0.0D0) THEN
C                 WRITE Error Message: VOC2ROC Out-of-Range
                  CALL ERRHDL(PATH,MODNAM,'E','804',SRCID(I))
               END IF
C              Assign The Field
               AROC_RATIO(I) = DNUM
            END IF
         END DO
      END IF

 999  RETURN
      END     
      

C***********************************************************************
      SUBROUTINE SHREMIS
C***********************************************************************
C        SHREMIS Module of AERMOD
C
C        PURPOSE: To process Secondary Hourly Emissions Data
C
C        PROGRAMMER: Huy Tran; based on HREMIS (soset.f) by Jayant Hardikar, Roger Brode
C
C        DATE:    September 25, 2023
C
C        MODIFIED: November 28, 2023
C                  Since HREMIS in aermod_23123 is identical in aermod_22112,
C                  no update made to SHREMIS which was based on HREMIS
C
C        INPUTS:  Pathway (SO) and Keyword (SHOURLY)
C
C        OUTPUTS: Source SQFLAG Array
C
C        CALLED FROM:   SOCARD
C***********************************************************************

C     Variable Declarations
      USE MAIN1
C     USE BUOYANT_LINE, ONLY: L_BLHOURLY   ! Multiple_BuoyLines_D41_Wood; deactivated here
      IMPLICIT NONE
      CHARACTER MODNAM*12

      INTEGER :: I, K, IH, IL

      LOGICAL FOPEN, INGRP
      LOGICAL RMARK

      CHARACTER (LEN=12) :: LOWID, HIGID, LID1, LID2, HID1, HID2, TEMPID

C     Variable Initializations
      MODNAM = 'SHREMIS'

      FOPEN  = .FALSE.

      IF (IFC .GE. 4) THEN
C        Retrieve Hourly Emissions Data Filename as Character Substring to
C        Maintain Case
         IF ((LOCE(3)-LOCB(3)) .LE. (ILEN_FLD - 1) ) THEN
C           Retrieve Filename as Character Substring to Maintain Original Case
C           Also Check for Filename Larger Than ILEN_FLD Characters
            SHRFILE = RUNST1(LOCB(3):LOCE(3))
         ELSE
C           WRITE Error Message:  SHRFILE Field is Too Long; same error code 291
            WRITE(DUMMY,'(I8)') ILEN_FLD
            CALL ERRHDL(PATH,MODNAM,'E','291',DUMMY)
            GO TO 999
         END IF

C        Open Hourly Emissions Data File If Not Already Open
         INQUIRE (FILE=SHRFILE,OPENED=FOPEN)

         IF (.NOT. FOPEN) THEN
C           Open Hourly Emissions Data File If Not Already Open
C           Open with ACTION='READ' to prevent overwrite and allow multiple access
            INQUIRE (UNIT=ISHREMI,OPENED=FOPEN)
            IF (.NOT. FOPEN) THEN
               OPEN(UNIT=ISHREMI,ERR=998,FILE=SHRFILE,IOSTAT=IOERRN,
     &              ACTION='READ',STATUS='OLD')
            END IF
         END IF

      ELSE
C        WRITE Error Message         ! Not Enough Parameters Specified; same error code 201
         CALL ERRHDL(PATH,MODNAM,'E','201',KEYWRD)
         GO TO 999
      END IF

      TEMPID = FIELD(4)

C     Set Up The Source Group Array
      IF (TEMPID .EQ. 'ALL') THEN
         DO K = 1, NUMSRC
            SQFLAG(K) = 'HOURLY'

C        Multiple_BuoyLines_D41_Wood: Begin; deactivated here
C           Set the flag indicating that there is one
C            (or more) buoyant line sources in the hourly emissions file
C           IF (SRCTYP(K) .EQ. 'BUOYLINE' .AND. .NOT. L_BLHOURLY) THEN
C              L_BLHOURLY = .TRUE.
C           END IF
C        Multiple_BuoyLines_D41_Wood: End

         END DO
      ELSE
C        Loop Through Fields
         DO I = 4, IFC
            CALL FSPLIT(PATH,KEYWRD,FIELD(I),ILEN_FLD,'-',RMARK,
     &                  LOWID,HIGID)
C           First Check Range for Upper Value < Lower Value
            CALL SETIDG(LOWID,LID1,IL,LID2)
            CALL SETIDG(HIGID,HID1,IH,HID2)
            IF ((HID1.LT.LID1) .OR. (IH.LT.IL) .OR. (HID2.LT.LID2)) THEN
C              WRITE Error Message:  Invalid Range,  Upper < Lower
               CALL ERRHDL(PATH,MODNAM,'E','203','SRCRANGE')
               CYCLE
            END IF
            DO K = 1, NUMSRC
               CALL ASNGRP(SRCID(K),LOWID,HIGID,INGRP)
               IF (INGRP) THEN
                  SQFLAG(K) = 'HOURLY'
C              Multiple_BuoyLines_D41_Wood: Begin; deactivated here
C                 IF ((SRCTYP(K) .EQ. 'BUOYLINE' .AND. 
C    &                           .NOT. L_BLHOURLY)) THEN
C                    L_BLHOURLY = .TRUE.
C                 END IF
C              Multiple_BuoyLines_D41_Wood: End
               END IF
            END DO
         END DO
      END IF

      GO TO 999

C     Process Error Messages; same error code 500
998   CALL ERRHDL(PATH,MODNAM,'E','500',KEYWRD)

999   RETURN
      END     
      

C***********************************************************************
      SUBROUTINE HRSQREAD (IS)
C***********************************************************************
C*         HRSQREAD Module of AERMOD
C*
C*         PURPOSE: To Assign Secondary Hourly Emission for GRS7/MGRS
C*
C*         PROGRAMMER:  Huy Tran; based on HRQREAD (aermod.f) by Jayant Hardikar, Roger Brode
C*
C*         DATE:    September 25, 2023
C*
C*         MODIFIED: November 28, 2023
C*                  Update with new subroutine CENT_DATE in aermod_22112; and
C*                  fix error by changing leftover HRQS to HRSQS
C*
C*         INPUTS:  Current Source Number Being Processed
C*
C*         OUTPUTS: Source Arrays
C*
C*         Revision History:
C*
C*         CALLED FROM:  HRLOOP
C************************************************************************
C*
C*    Variable Declarations
      USE MAIN1
      IMPLICIT NONE
      CHARACTER MODNAM*12

      INTEGER :: I, IS
      INTEGER :: IHYEAR, IHMON, IHDAY, IHHOUR, IHYEAR2
      INTEGER :: ILSAVE
      CHARACTER (LEN=20) :: RDFRM

      CHARACTER (LEN=12) :: HRSOID

C*    Variable Initializations
      MODNAM = 'HRSQREAD'

C*    Assign IQLINE counter to ILINE for passing to ERRHDL if needed, save as ILSAVE first
      
C     write(*,*)"HTdbg: HRSQREAD; ILINE: ",ILINE,
C    &          "IQLINE: ",IQLINE," ISQLINE: ",ISQLINE
C     ILSAVE = ILINE
C     ILINE  = ISQLINE

C*    READ Record to Buffers, A'num' and 'num'A1, where num=ISTRG
C*    Length of ISTRG is Set in PARAMETER Statement in MAIN1
C     Setup READ format and ECHO format for runstream record,
C     based on the ISTRG PARAMETER (set in MAIN1)
      WRITE(RDFRM,9100) ISTRG, ISTRG
 9100 FORMAT('(A',I4.4,',T1,',I4.4,'A1)')
      READ (ISHREMI,RDFRM,END=888,ERR=99) RUNST1, (RUNST(I), I=1, ISTRG)
C*
C*    Convert Lower Case to Upper Case Letters              ---   CALL LWRUPR
      CALL LWRUPR
C*
C*    Define Fields on Card                                 ---   CALL DEFINE
      CALL DEFINE
C*
C*    Get the Contents of the Fields                        ---   CALL GETFLD
      CALL GETFLD
C*
C*    Check for number of fields - error if less than 7.
C     write(*,*)"HTdbg: HRSQREAD; IFC: ",IFC
      IF (IFC .LT. 7) THEN
         WRITE(DUMMY,'(I8)') KURDAT
         CALL ERRHDL(PATH,MODNAM,'E','805',DUMMY)
         RUNERR = .TRUE.
         GO TO 999
      END IF
C*
C*    Assign the Fields to Local Varables and Check The Numerical Field
C*
C*    Date and time variables common to all source types
C*
      CALL STONUM(FIELD(3), ILEN_FLD, FNUM, IMIT)
      IHYEAR = NINT(FNUM)
      IF (IMIT .NE. 1) THEN
         CALL ERRHDL(PATH,MODNAM,'E','208','VOCEHRLY')
         RUNERR = .TRUE.
         GO TO 999
      END IF
C     write(*,*)"HTdbg: HRSQREAD; IHYEAR: ",IHYEAR

      CALL STONUM(FIELD(4), ILEN_FLD, FNUM, IMIT)
      IHMON = NINT(FNUM)
      IF (IMIT .NE. 1) THEN
         CALL ERRHDL(PATH,MODNAM,'E','208','VOCEHRLY')
         RUNERR = .TRUE.
         GO TO 999
      END IF
C     write(*,*)"HTdbg: HRSQREAD; IHMON: ",IHMON

      CALL STONUM(FIELD(5), ILEN_FLD, FNUM, IMIT)
      IHDAY = NINT(FNUM)
      IF (IMIT .NE. 1) THEN
         CALL ERRHDL(PATH,MODNAM,'E','208','VOCEHRLY')
         RUNERR = .TRUE.
         GO TO 999
      END IF
C     write(*,*)"HTdbg: HRSQREAD; IHDAY: ",IHDAY

      CALL STONUM(FIELD(6), ILEN_FLD, FNUM, IMIT)
      IHHOUR = NINT(FNUM)
      IF (IMIT .NE. 1) THEN
         CALL ERRHDL(PATH,MODNAM,'E','208','VOCEHRLY')
         RUNERR = .TRUE.
         GO TO 999
      END IF
C     write(*,*)"HTdbg: HRSQREAD; IHHOUR: ",IHHOUR

C      D001 Call CENT_DATE to determine the current Julian Day and Calculate Current Gregorian Date First Convert Year to 4-Digit Value Wood 9/15/22
       IF (IHYEAR .LE. 99) THEN
         CALL CENT_DATE(IHYEAR2,IHYEAR) 
       END IF
C ---  D001 remove original calculation of 4-Digit year Wood 9/15/22
C --- Check for use of 2-digit year in HOUREMIS file, adjust to 4-digit
C     year for comparison with FULLDATE based on met data file
C     IF (IHYEAR .LE. 99) THEN
C        IHYEAR2 = IHYEAR
C        IF (IHYEAR2 .GE. ISTRT_WIND .AND.
C    &                        IHYEAR2 .LE. 99) THEN
C           IHYEAR = ISTRT_CENT*100 + IHYEAR2
C        ELSE IF (IHYEAR2 .LT. ISTRT_WIND) THEN
C           IHYEAR = (ISTRT_CENT+1)*100 + IHYEAR2
C        END IF
C     END IF
C     write(*,*)"HTdbg: HRSQREAD; 2nd IHYEAR: ",IHYEAR

C --- Calculate current date (YYYYMMDDHH) from HOUREMIS file record, FULLHRQ
      FULLHRSQ = IHYEAR*1000000 + IHMON*10000 + IHDAY*100 + IHHOUR
C     write(*,*)"HTdbg: HRSQREAD; FULLHRQ: ",FULLHRQ,
C    &          " FULLHRSQ: ",FULLHRSQ

C --- Assign source ID but check for field length > 12 first
      IF( LEN_TRIM(FIELD(7)) .LE. 12 ) THEN
         HRSOID = FIELD(7)
      ELSE
         HRSOID = FIELD(7)(1:12)
      END IF
C     write(*,*)"HTdbg: HRSQREAD; HRSOID: ",HRSOID,
C    &          "SRCID: ",SRCID(IS)

C*    Check for Source ID Consistency ; If Failed Issue Error
      IF ( HRSOID .NE. SRCID(IS) ) THEN
         WRITE(DUMMY,'(A12)') SRCID(IS)
         CALL ERRHDL(PATH,MODNAM,'E','342',SRCID(IS))
         RUNERR = .TRUE.
         GO TO 999
      END IF

      IF (IFC .EQ. 7) THEN
C*       All parameters missing for this hour/source - WRITE Warning Message
C*       Assign zeros to all parameters
         WRITE(DUMMY,'(I10.10)') FULLHRSQ
         CALL ERRHDL(PATH,MODNAM,'W','808',DUMMY)
         HRSQS = 0.0D0
         GO TO 999
      ELSE IF (IFC.GT.8) THEN
C*       Too many parameters - WRITE Error Message
C*       Assign zeros to all parameter
         WRITE(DUMMY,'(I10.10)') FULLHRSQ
         CALL ERRHDL(PATH,MODNAM,'E','807',DUMMY)
         HRSQS = 0.0D0
         RUNERR = .TRUE.
         GO TO 999
      ELSE IF (IFC.LT.7) THEN
C*       Too few parameters - WRITE Error Message
C*       Assign zeros to all parameter
         WRITE(DUMMY,'(I10.10)') FULLHRSQ
         CALL ERRHDL(PATH,MODNAM,'E','807',DUMMY)
         HRSQS = 0.0D0
         RUNERR = .TRUE.
         GO TO 999
      ELSE IF (IFC.EQ.8) THEN
C        FIELD(8) contain hourly emission values
         CALL STODBL(FIELD(8), ILEN_FLD, HRSQS, IMIT)
         IF (IMIT .NE. 1) THEN
            CALL ERRHDL(PATH,MODNAM,'E','208','VOCEHRLY')
            RUNERR = .TRUE.
C ---    Check for large negative values, could be a missing indicator
         ELSE IF ( HRSQS .LE. -90.0D0 ) THEN
C*          Assume emissions are missing; assign value of 0.0 and issue Warning;
C*          Note that AERMOD User's Guide (p. 3-38, 3-39) implies that blanks
C*          should be used for missing values.
            HRSQS = 0.0D0
            WRITE(DUMMY,'(I10.10)') FULLHRSQ
            CALL ERRHDL(PATH,MODNAM,'W','341',DUMMY)
         END IF
         GO TO 999
      END IF

C*    Write Error Message for Error Reading Hourly Emissions File
 99   CALL ERRHDL(PATH,MODNAM,'E','510','VOCEHRLY')
      RUNERR = .TRUE.
      GO TO 999

888   CONTINUE

      EOF = .TRUE.

999   RETURN
      END      


C***********************************************************************
      SUBROUTINE SEMFACT (SQARG)
C***********************************************************************
C                 SEMFACT Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE: Applies Variable Emission Rate and
C                 Unit Conversion Factors of Secondary Emission Stream
C
C        PROGRAMMER: Huy Tran, based on EMFACT (calc2.f) by Roger Brode, Jeff Wang
C        MODIFIED  : for handling OpenPit Source Type - PES, 7/26/94
C
C        DATE:    Septemer 15, 2023
C
C        MODIFIED: November 28, 2023
C                  Since EMFACT in aermod_23123 is identical in aermod_22112,
C                  no update made to SEMFACT which was based on EMFACT
C
C        INPUTS:  Arrays of Source Parameters
C                 Date and Hour
C                 Meteorological Variables for One Hour
C                 Variable Emission Rate Flags and Factors
C                 Unit Conversion Rate Factors
C
C        OUTPUTS: Adjusted Emission Rate, SQTK
C
C        CALLED FROM:   PCALC
C                       VCALC
C                       ACALC
C***********************************************************************

C     Variable Declarations
      USE MAIN1
      IMPLICIT NONE
      CHARACTER MODNAM*12

      DOUBLE PRECISION :: SQARG

C     Variable Initializations
      MODNAM = 'SEMFACT'

C --- Apply Variable Emission Rate Factor, Based on Value of SQFLAG
C     Emission unit factor is applied later since it varies by
C     output type
      IF (SQFLAG(ISRC) .EQ. ' ') THEN
         SQTK = SQARG

C*----   ISCSTM Modification: To handle hourly emissions - jah 11/4/94
      ELSE IF (SQFLAG(ISRC) .EQ. 'HOURLY') THEN
         SQTK = SQARG
C*----
C*#

      ELSE IF (SQFLAG(ISRC) .EQ. 'MONTH') THEN
         SQTK = SQARG * SQFACT(IMONTH,ISRC)

      ELSE IF (SQFLAG(ISRC) .EQ. 'HROFDY') THEN
         SQTK = SQARG * SQFACT(IHOUR,ISRC)

      ELSE IF (SQFLAG(ISRC) .EQ. 'WSPEED') THEN
         SQTK = SQARG * SQFACT(IUCAT,ISRC)

      ELSE IF (SQFLAG(ISRC) .EQ. 'SEASON') THEN
         SQTK = SQARG * SQFACT(ISEAS,ISRC)

      ELSE IF (SQFLAG(ISRC) .EQ. 'SEASHR') THEN
         SQTK = SQARG * SQFACT((IHOUR+(ISEAS-1)*24),ISRC)

      ELSE IF (SQFLAG(ISRC) .EQ. 'HRDOW') THEN
         SQTK = SQARG * SQFACT((IHOUR +
     &        (IDAY_OF_WEEK-1)*24),ISRC)

      ELSE IF (SQFLAG(ISRC) .EQ. 'HRDOW7') THEN
         SQTK = SQARG * SQFACT((IHOUR +
     &        (IDAY_OF_WEEK7-1)*24),ISRC)

      ELSE IF (SQFLAG(ISRC) .EQ. 'SHRDOW') THEN
         SQTK = SQARG * SQFACT((IHOUR+(ISEAS-1)*24+
     &        (IDAY_OF_WEEK-1)*96),ISRC)

      ELSE IF (SQFLAG(ISRC) .EQ. 'SHRDOW7') THEN
         SQTK = SQARG * SQFACT((IHOUR+(ISEAS-1)*24+
     &        (IDAY_OF_WEEK7-1)*96),ISRC)

      ELSE IF (SQFLAG(ISRC) .EQ. 'MHRDOW') THEN
         SQTK = SQARG * SQFACT((IHOUR+(IMONTH-1)*24+
     &        (IDAY_OF_WEEK-1)*288),ISRC)

      ELSE IF (SQFLAG(ISRC) .EQ. 'MHRDOW7') THEN
         SQTK = SQARG * SQFACT((IHOUR+(IMONTH-1)*24+
     &        (IDAY_OF_WEEK7-1)*288),ISRC)

      END IF

      RETURN
      END      


C***********************************************************************
      SUBROUTINE ROCVAL
C***********************************************************************
C        ROCVAL Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE: Processes Non-temporally-varying ROC Value Option,
C                 CO ROCVAL
C
C        PROGRAMMER: Huy Tran; based on O3VAL (coset.f) by Roger W. Brode
C
C        DATE:    Sep 15, 2023
C
C        MODIFIED: November 28, 2023
C                  Since O3VAL in aermod_23123 is identical in aermod_22112,
C                  no update made to ROCVALS which was based on O3VALS
C                  Corrected ROCVAL to ROCVAL for non-temporally varying input
C
C        INPUTS:  Input Runstream Image Parameters
C
C        OUTPUTS:
C
C        CALLED FROM:   COCARD
C***********************************************************************

C     Variable Declarations
      USE MAIN1
      IMPLICIT NONE

      INTEGER :: I

      CHARACTER MODNAM*12

C     Variable Initializations
      MODNAM = 'ROCVALS'

C     GO TO 999 ! Skip everything, i.e., do nothing right now

C --- Check The Number Of The Fields, accounting for sector-varying values
      IF (.NOT.L_ROCSector) THEN
         IF (IFC .LE. 2) THEN
C           Error Message: No Parameters
            CALL ERRHDL(PATH,MODNAM,'E','200',KEYWRD)
            GO TO 999
         ELSE IF (IFC .GT. 4) THEN
C           Error Message: Too Many Parameters
            CALL ERRHDL(PATH,MODNAM,'E','202',KEYWRD)
            GO TO 999
         END IF
C ---    Check for SECT ID in field 3 in case ROCSECTOR keyword was omitted
         IF (FIELD(3)(1:4) .EQ. 'SECT') THEN
C           Error Message: SECT ID without O3SECTOR keyword
            CALL ERRHDL(PATH,MODNAM,'E','171',KEYWRD)
            GO TO 999
         END IF
C ---    Assign sector ID to 1 since sector-varying values not being used;
C        also set field index for the user-specified O3VALUES option and
C        assign the option to O3FLAG variable
         IROCSECT = 1
         I = 3
         L_ROCVAL(IROCSECT) = .TRUE.

      ELSE
C ---    Process inputs based on ROCSECTOR option
         IF (IFC .LE. 2) THEN
C           Error Message: No Parameters
            CALL ERRHDL(PATH,MODNAM,'E','200',KEYWRD)
            GO TO 999
         ELSE IF (IFC .EQ. 4) THEN
            IF (FIELD(3)(1:4) .NE. 'SECT') THEN
C              Error Message: Invalid sector field
               CALL ERRHDL(PATH,MODNAM,'E','203','O3SECTOR ID')
               GO TO 999
            ELSE
C              Error Message: No Numerical Parameters
               CALL ERRHDL(PATH,MODNAM,'E','201',KEYWRD)
               GO TO 999
            END IF
         ELSE IF (IFC .GT. 5) THEN
C           Error Message: Too Many Parameters
            CALL ERRHDL(PATH,MODNAM,'E','202',KEYWRD)
            GO TO 999
         END IF
C ---    Determine user-specified sector
         IF (FIELD(3) .EQ. 'SECT1') THEN
            IROCSECT = 1
         ELSE IF (FIELD(3) .EQ. 'SECT2') THEN
            IROCSECT = 2
         ELSE IF (FIELD(3) .EQ. 'SECT3') THEN
            IROCSECT = 3
         ELSE IF (FIELD(3) .EQ. 'SECT4') THEN
            IROCSECT = 4
         ELSE IF (FIELD(3) .EQ. 'SECT5') THEN
            IROCSECT = 5
         ELSE IF (FIELD(3) .EQ. 'SECT6') THEN
            IROCSECT = 6
         ELSE
C           Error Message: Invalid sector definition
            CALL ERRHDL(PATH,MODNAM,'E','203','ROCSECTOR ID')
            GO TO 999
         END IF

C ---    Set field index for the user-specified Ozone Value
         I = 4
         L_ROCVAL(IROCSECT) = .TRUE.

      END IF

C     Get ROC Value, ROCBACK, for applicable sector
      CALL STODBL(FIELD(I),ILEN_FLD,DNUM,IMIT)
C     Check The Numerical Field
      IF (IMIT .NE. 1) THEN
         CALL ERRHDL(PATH,MODNAM,'E','208',KEYWRD)
         GO TO 999
      END IF

C     Assign value to ROCBACK variable for this sector
      ROCBACK(IROCSECT) = DNUM
C     write(*,*)"HTdbg: ROCVALS, IROCSECT: ",IROCSECT,
C    &          " ROCBACK(IROCSECT): ",ROCBACK(IROCSECT)

C     Check for units of ozone value
      IF (IFC .EQ. I+1) THEN
         IF (FIELD(I+1).EQ.'PPM' .OR. FIELD(I+1).EQ.'PPB' .OR.
     &       FIELD(I+1).EQ.'UG/M3') THEN
            ROCVALUNITS = FIELD(I+1)
         ELSE
C           Write Error Message:  Invalid units for ozone value
            CALL ERRHDL(PATH,MODNAM,'E','203',' ROCUNITS')
         END IF
      END IF

C     Convert unit to ug/m3 to make it easier to convert back to ppb unit later
      IF (ROCVALUNITS .EQ. 'PPB') THEN
         ROCBACK(IROCSECT) = ROCBACK(IROCSECT) / VOC_PPB
      ELSE IF (ROCVALUNITS .EQ. 'PPM') then
         ROCBACK(IROCSECT) = ROCBACK(IROCSECT) / VOC_PPM
      END IF

C     Check range of value
      IF (ROCBACK(IROCSECT) .LE. 0.0D0 .OR.
     &    ROCBACK(IROCSECT) .GT. 500.0D0)THEN
         CALL ERRHDL(PATH,MODNAM,'W','320',' ROCBACK ')
      END IF

 999  RETURN
      END


C***********************************************************************
      SUBROUTINE ROC_UNIT
C***********************************************************************
C                 ROC_UNIT Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE: Processes user-specified units for ROCVAL keyword
C                 ROC concentrations, based on the ROCUNIT keyword
C
C        PROGRAMMER:  Huy Tran; based on OZON_UNIT (coset.f) by Roger Brode
C
C        DATE:     Sep 15, 2023
C        MODIFIED: November 28, 2023
C                  Since OZONE_UNIT in aermod_23123 is identical in aermod_22112,
C                  no update made to ROC_UNIT which was based on OZON_UNIT
C
C        INPUTS:  Input Runstream Image Parameters
C
C        OUTPUTS: Variable Emmission Rate Factors
C
C        CALLED FROM:   COCARD
C***********************************************************************

C     Variable Declarations
      USE MAIN1
      IMPLICIT NONE
      CHARACTER MODNAM*12

C     Variable Initializations
      MODNAM = 'ROC_UNIT'

C     Check The Number Of The Fields
      IF (IFC .LE. 2) THEN
C        Error Message: No Parameters
         CALL ERRHDL(PATH,MODNAM,'E','200',KEYWRD)
         GO TO 999
      ELSE IF (IFC .GT. 3) THEN
C        Error Message: Too Many Parameters
         CALL ERRHDL(PATH,MODNAM,'E','202',KEYWRD)
         GO TO 999
      END IF

C     Check for units of background values
      IF (FIELD(3).EQ.'PPM' .OR. FIELD(3).EQ.'PPB' .OR.
     &    FIELD(3).EQ.'UG/M3') THEN
          ROCUnits = FIELD(3)
      ELSE
C        Write Error Message:  Invalid units for O3VALUES
         CALL ERRHDL(PATH,MODNAM,'E','203','ROCUnits')
      END IF

 999  RETURN
      END
      
      SUBROUTINE SEMUNIT
C***********************************************************************
C                 SEMUNIT Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE: Processes Emission Rate Unit Conversion Factors for secondary emissions
C
C        PROGRAMMER: Huy Tran; based on EMUNIT (soset.f) by Jeff Wang, Roger Brode
C
C        DATE:    Sep 15, 2023
C        MODIFIED: November 28, 2023
C                  no update needed for aermod 23132
C
C        INPUTS:  Input Runstream Image Parameters
C
C        OUTPUTS: Emission Rate Unit Conversion Factors
C
C        CALLED FROM:   SOCARD
C***********************************************************************

C     Variable Declarations
      USE MAIN1
      IMPLICIT NONE
      CHARACTER MODNAM*12

C     Variable Initializations
      MODNAM = 'SEMUNIT'

C     Check The Number Of The Fields
      IF (IFC .LE. 2) THEN
C        Error Message: No Parameters
         CALL ERRHDL(PATH,MODNAM,'E','200',KEYWRD)
         GO TO 999
      ELSE IF (IFC .LT. 5) THEN
C        Error Message: Not Enough Parameters
         CALL ERRHDL(PATH,MODNAM,'E','201',KEYWRD)
         GO TO 999
      ELSE IF (IFC .GT. 5) THEN
C        Error Message: Too Many Parameters
         CALL ERRHDL(PATH,MODNAM,'E','202',KEYWRD)
         GO TO 999
      END IF

C     Fetch Each Field
      CALL STODBL(FIELD(3),ILEN_FLD,DNUM,IMIT)
C     Check The Numerical Field
      IF (IMIT .NE. 1) THEN
         CALL ERRHDL(PATH,MODNAM,'E','208',KEYWRD)
         GO TO 999
      END IF

      SEMIFAC(1) = DNUM
      SEMILBL(1) = FIELD(4)
      SOUTLBL(1) = FIELD(5)
      SPERLBL(1) = FIELD(5)

 999  RETURN
      END

C***********************************************************************
      SUBROUTINE PRM_PLUME2 (ZARG, UEFFPH, COUT1, COUT2)                   ! ORD (EMM) change
C***********************************************************************
C             PRM_PLUME Module of the AMS/EPA Regulatory Model - AERMOD
C
C   PURPOSE: Calculate the contribution to the concentration due to
C            PRIME downwash component
C
C   PROGRAMMER: Roger Brode, PES, Inc.
C
C   DATE:    July 5, 2001
C
C   INPUTS:  Receptor height, ZARG
C
C   OUTPUTS: Contribution due to PRIME, COUT
C
C   CALLED FROM:   PRM_PCHI
C
C***********************************************************************

C     Variable Declarations
      USE MAIN1
      IMPLICIT NONE
      CHARACTER MODNAM*12
      DOUBLE PRECISION :: ZARG, COUT1, COUT2
      DOUBLE PRECISION :: UEFFPH

      MODNAM = 'PRM_PLUME2'

C     Assign receptor height for vertical term calculations
      ZR = ZARG

      IF (STABLE) THEN
         CALL VRTSBN (SZ, HE)
      ELSE IF (UNSTAB .AND. HE.LE.ZI) THEN
         CALL VRTSBL (SZ, HE, ZI)
      ELSE
         FSUBZ = 0.0D0
      END IF

C --- Calculate the WRAP term for a stable atmosphere
C      UEFFPH is defined properly in PRM_PCHI (the calling program)
C      for regulatory AERMOD, AWMA_DOWNWASH, and ORD_DOWNWASH
C      COUT = (QTK / US) * ( FSUBY * FSUBZ )                    ! Original AERMOD
      COUT1 = (QTK / UEFFPH) * ( FSUBY * FSUBZ )
      COUT2 = (SQTK / UEFFPH) * ( FSUBY * FSUBZ )
c      IF (AWMADWDBG) THEN
c         WRITE(AWMADWDBUNT,*) '     Contribution due to PRIME:'
c         WRITE(AWMADWDBUNT,*) '     SY=',SY,' SZ=',SZ,
c     &                        ' UEFFPH=',UEFFPH
c         WRITE(AWMADWDBUNT,*) '     FSUBY=',FSUBY,  ', FSUBZ=',FSUBZ
c         WRITE(AWMADWDBUNT,*) '     QTK / UEFFPH)*( FSUBY*FSUBZ )=',COUT
c      END IF

      RETURN
      END


C***********************************************************************
      SUBROUTINE AER_SPCHI( XARG, ADJ, VDINP, JIN, AEROUT, SAEROUT )
C***********************************************************************
C        AER_SPCHI Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE: Calculates Hourly Concentration for POINT Sources
C                 Using Gaussian Plume Equation for both primary and secondary emissions
C
C        PROGRAMMER: Huy Tran; based on AER_PCHI (calc1.f) by Roger Brode, PES, Inc.
C
C        DATE:    Sep 25, 2023
C
C        MODIFIED: November 28, 2023
C                  modified for AER_SPCHI to process both primary and 2ndary streams 
C                  instead of just 2ndary like September 2023 version
C
C        INPUTS:  Distance, XARG (downwind for plume; radial for pancake)
C                 Crosswind Distance
C                 Plume Height
C                 Stack Top Wind Speed
C                 Lateral Dispersion Parameter
C                 Vertical Dispersion Parameter
C                 Stability Class
C                 Mixing Height
C                 Receptor Height Above Ground
C                 Emission Rate and Units Scaling Factor
C                 Source Parameter Arrays
C
C        OUTPUTS: AEROUT, AERMOD Concentration for Particular
C                 Source/Receptor Combination from primary emission
C                 SAEROUT, AERMOD Concentration for Particular
C                 Source/Receptor Combination from secondary emission
C
C        CALLED FROM:   AERCALC, VOLCALC, ACALC
C***********************************************************************

C     Variable Declarations
      USE MAIN1
      IMPLICIT NONE
      INTEGER :: JIN
      DOUBLE PRECISION :: AEROUT(NUMTYP), XARG, ADJ, VDINP, DRYFLUX,
     &                    WETFLUX,
     &                    SAEROUT(NUMTYP)
      CHARACTER MODNAM*12

C     Variable Initializations
      MODNAM = 'AER_SPCHI'
      DRYFLUX = 0.0D0
      WETFLUX = 0.0D0

C---- Calculate the contribution due to horizontal plume, CWRAP
      IF (FOPT .EQ. 0.0D0) THEN
         CWRAP = 0.0D0
         SCWRAP = 0.0D0
      ELSE
         CALL SCPLUME (ZRT, CWRAP, SCWRAP) ! Call SCPLUME from grs7_supl.f
      END IF

C---- Calculate the contribution due to terrain-following plume, CLIFT
      IF (ZRT .EQ. ZFLAG) THEN
C----    Effective receptor heights are equal, therefore CLIFT = CWRAP
         CLIFT = CWRAP
         SCLIFT = SCWRAP
      ELSE IF (FOPT .EQ. 1.0D0) THEN
         CLIFT = 0.0D0
         SCLIFT = 0.0D0
      ELSE
         CALL SCPLUME (ZFLAG, CLIFT, SCLIFT)
      END IF

C---- Calculate the exponential decay term, D               ---   CALL DECAY
      Call DECAY (XARG)

C---- Calculate the hourly concentration and deposition values
      ITYP = 0
      IF (CONC) THEN
         ITYP = 1
         AEROUT(ITYP) = ADJ * EMIFAC(ITYP) *
     &                 (FOPT * CWRAP + (1.0D0 - FOPT) * CLIFT) * D
         SAEROUT(ITYP) = ADJ * SEMIFAC(ITYP) *
     &                 (FOPT * SCWRAP + (1.0D0 - FOPT) * SCLIFT) * D

C   ENHANCEMENT TO DEBUG OUTPUT BASED ON ENSR
         IF (DEBUG) THEN
            WRITE(DBGUNT,10) ITYP, ADJ, FOPT, CWRAP, CLIFT,
     &           SCWRAP, SCLIFT, D, AEROUT(ITYP), SAEROUT(ITYP) 
10          FORMAT(/,'ITYP = ',I2,' - CONC:',
     &             /,'AEROUT(ITYP) = ADJ * EMIFAC(ITYP) * (FOPT * ',
     &               'CWRAP + (1.0 -FOPT) * CLIFT) * D',
     &             /,'SAEROUT(ITYP) = ADJ * SEMIFAC(ITYP) * (FOPT * ',
     &               'SCWRAP + (1.0 -FOPT) * SCLIFT) * D',
     &             /,' ADJ   = ',G16.8,
     &             /,' FOPT  = ',G16.8,
     &             /,' CWRAP = ',G16.8,
     &             /,' CLIFT = ',G16.8,
     &             /,' SCWRAP = ',G16.8,
     &             /,' SCLIFT = ',G16.8,
     &             /,' D     = ',G16.8,
     &             /,' AEROUT(ITYP) = ',G16.8,
     &             /,' SAEROUT(ITYP) = ',G16.8,/)
         ENDIF

      END IF

      IF (DEPOS .OR. DDEP) THEN
C        Calculate DRYFLUX, vertical term for wet deposition
C----    Calculate the contribution due to horizontal plume, CWRAP
         IF (FOPT .EQ. 0.0D0) THEN
            CWRAP = 0.0D0
         ELSE
            CALL CPLUME (ZRT-ZFLAG+ZRDEP, CWRAP)
         END IF

C----    Calculate the contribution due to terrain-following plume, CLIFT
         IF (ZRT .EQ. ZFLAG) THEN
C----       Effective receptor heights are equal, therefore CLIFT = CWRAP
            CLIFT = CWRAP
         ELSE IF (FOPT .EQ. 1.0D0) THEN
            CLIFT = 0.0D0
         ELSE
            CALL CPLUME (ZRDEP, CLIFT)
         END IF

         DRYFLUX = (FOPT * CWRAP + (1.0D0 - FOPT) * CLIFT) * D
     &
      END IF

      IF (DEPOS .OR. WDEP) THEN
C        Calculate WETFLUX, vertical term for wet deposition.
C        Note that the SRT2PI for the integrated vertical term
C        has been removed since it should be divided by SRT2PI.
C        Additional factor of 3600. has been added to denominator
C        to account for conversion from seconds to hours when
C        divided by wind speed below.
         IF (PRATE .GT. 0.0D0) THEN
            IF (NPD .EQ. 0) THEN
               WETFLUX = (ADJ*FRACSAT*PRATE*1.0D6*RGAS*TA)/
     &                   (ZSUBP*HENRY(ISRC)*1.0D9*DENOM*3600.0D0)
            ELSE
               WETFLUX = 1.0D-3*ADJ*WASHOUT(JIN)*PRATE/
     &                   (ZSUBP*3600.0D0)
            END IF
         ELSE
            WETFLUX = 0.0D0
         END IF
      END IF

      IF (DEPOS) THEN
         ITYP = ITYP + 1
         IF( STABLE  .OR.  (UNSTAB .AND. (HS .GE. ZI) ) ) THEN
            AEROUT(ITYP) = ADJ * VDINP * EMIFAC(ITYP) * DRYFLUX +
     &                     QTK * WETFLUX * EMIFAC(ITYP) * FSUBY/UEFF
         ELSE IF (UNSTAB) THEN
            AEROUT(ITYP) = ADJ * VDINP * EMIFAC(ITYP) * DRYFLUX +
     &                     QTK * WETFLUX * EMIFAC(ITYP) *
     &                     (PPF*FSUBY3/UEFF3+(1.0D0-PPF)*FSUBY/UEFFD)
         END IF

         IF (DEBUG) THEN
            WRITE(DBGUNT,11) ITYP, ADJ, VDINP, DRYFLUX, WETFLUX,
     &                       AEROUT(ITYP)
11          FORMAT(/,'ITYP = ',I2,' - DEPOS:',
     &             /,' ADJ     = ',G16.8,
     &             /,' VDINP   = ',G16.8,
     &             /,' DRYFLUX = ',G16.8,
     &             /,' WETFLUX = ',G16.8,
     &             /,' AEROUT(ITYP) = ',G16.8,/)
         END IF

      END IF

      IF (DDEP) THEN
         ITYP = ITYP + 1
         AEROUT(ITYP) = ADJ * VDINP * EMIFAC(ITYP) * DRYFLUX

         IF (DEBUG) THEN
            WRITE(DBGUNT,12) ITYP, ADJ, VDINP, DRYFLUX,
     &                       AEROUT(ITYP)
12          FORMAT(/,'ITYP = ',I2,' - DDEP:',
     &             /,' ADJ     = ',G16.8,
     &             /,' VDINP   = ',G16.8,
     &             /,' DRYFLUX = ',G16.8,
     &             /,' AEROUT(ITYP) = ',G16.8,/)
         END IF

      END IF

      IF (WDEP) THEN
         ITYP = ITYP + 1
         IF( STABLE  .OR.  (UNSTAB .AND. (HS .GE. ZI) ) ) THEN
            AEROUT(ITYP) = QTK * WETFLUX * EMIFAC(ITYP) * FSUBY/UEFF
         ELSE IF (UNSTAB) THEN
            AEROUT(ITYP) = QTK * WETFLUX * EMIFAC(ITYP) *
     &                     (PPF*FSUBY3/UEFF3+(1.0D0-PPF)*FSUBY/UEFFD)
         END IF

         IF (DEBUG) THEN
            WRITE(DBGUNT,13) ITYP, ADJ, ZSUBP, PRATE, WETFLUX,
     &                       AEROUT(ITYP)
13         FORMAT(/,'ITYP = ',I2,' - WDEP:',
     &             /,' ADJ     = ',G16.8,
     &             /,' ZSUBP   = ',G16.8,
     &             /,' PRATE   = ',G16.8,
     &             /,' WETFLUX = ',G16.8,
     &             /,' AEROUT(ITYP) = ',G16.8,/)
         END IF

      END IF


CCRFL Call to METDEB was moved here from METEXT on 9/26/94, R.F. Lee.
CCRFL Print meteorological debug output.                   ---   CALL METDEB
      IF (METEORDBG) CALL METDEB

      IF ( DEBUG ) THEN
C        Print Out Debugging Information                    ---   CALL DEBOUT
         CALL DEBOUT
      END IF

      RETURN
      END

C***********************************************************************
      SUBROUTINE SCPLUME (ZARG, COUT, SCOUT)
C***********************************************************************
C             CPLUME Module of the AMS/EPA Regulatory Model - AERMOD
C
C   PURPOSE: Calculate the contribution to the concentration due to
C            plume component, either horizontal or terrain-following,
C            depending on the input receptor height, ZARG
C
C   PROGRAMMER: Roger Brode, PES, Inc.
C
C   DATE:    September 30, 1993
C
C   REVISIONS:
C               Make stable plume reflections dependent on the
C               developmental option switch, OPTG1 & OPTG2,
C               R. Brode, PES, 1/6/95
C
C               Remove stable plume reflections off of ZI for
C               Base Case model.  R. Brode, PES - 12/7/94
C
C               Revised emission rates for each plume to QTK*(1.-PPF)
C               for the direct and indirect plumes, and to QTK*PPF
C               for the penetrated plume.  Ref:  P.D.F. Model for
C               Dispersion in the Convective Boundary Layer,
C               J.C. Weil, 6/27/94. Changes made 7/19/94, R.F. Lee.
C
C               Added true centerline concentration calculations.
C               Changes made 7/25/94, R.F. Lee.
C
C   INPUTS:  Stability, STABLE/UNSTAB
C            Fraction of plume vertical flux remaining in the CBL, FOPT
C            Mixing height, ZI
C            Plume heights, HE/HED1/HED2/HEN1/HEN2
C            sigma_Z's: SZ, SZD1, SZD2, SZN1, SZN2, SZ3
C            Receptor height, ZARG
C
C   OUTPUTS: Contribution due to WRAP, CWRAP
C
C   CALLED FROM:   PCHI
C
C   Assumptions:  For receptor height (ZR) above the mixing height (ZI)
C                 for unstable conditions, the direct and indirect plume
C                 impacts are set to zero.
C
C   References:   "A Dispersion Model for the Convective Boundary
C                 Layer", J. Weil, 8/17/93
C
C***********************************************************************

C     Variable Declarations
      USE MAIN1
      IMPLICIT NONE
      CHARACTER MODNAM*12
      DOUBLE PRECISION :: ZARG, COUT, SCOUT

      MODNAM = 'SCPLUME'

C     Assign receptor height for vertical term calculations
      ZR = ZARG

      IF( STABLE  .OR.  (UNSTAB .AND. (HS .GE. ZI) ) ) THEN
C        Calculate the vertical term, FSUBZ              ---   CALL VRTSBL
C        With stable plume reflections and effective Zi
         IF (ZR .LE. HSBL) THEN
            CALL VRTSBL (SZ, MAX( 0.0D0, HE-HV ), HSBL)
         ELSE
            CALL VRTSBN (SZ, MAX( 0.0D0, HE-HV ))
         END IF

C        Calculate the concentration for a stable atmosphere
         COUT  = (QTK / UEFF) * ( FSUBY * FSUBZ )
         SCOUT = (SQTK / UEFF) * ( FSUBY * FSUBZ )

      ELSEIF( UNSTAB )THEN
         IF (PPF .LT. 1.0D0) THEN
C           Calculate the vertical term for the direct plume, FSUBZD
            IF (ZR .LE. ZI) THEN
C              Calculation for Receptor below Zi      ---   CALL VRTCBL
               CALL VRTCBL ( HED1-HV, HED2-HV, SZD1, SZD2, 1.0D0)
               FSUBZD = FSUBZ
            ELSE
C              Set FSUBZ = 0.0 for "receptor height" (ZR) > ZI
               FSUBZD = 0.0D0
            END IF

C           Calculate the vertical term for the indirect plume, FSUBZN
            IF (ZR .LE. ZI) THEN
C              Calculation for Receptor below Zi      ---   CALL VRTCBL
               CALL VRTCBL ( HEN1-HV, HEN2-HV, SZN1, SZN2, -1.0D0 )
               FSUBZN = FSUBZ
            ELSE
C              Set FSUBZ = 0.0 for "receptor height" (ZR) > ZI
               FSUBZN = 0.0D0
            END IF
         ELSE
            FSUBZD = 0.0D0
            FSUBZN = 0.0D0

         END IF

C        Note that UEFF and UEFF3 can never be zero, since they get
C        set to a minimum value earlier on.

         IF( PPF .GT. 0.0D0 )THEN
C           Calculate the vertical term for the penetrated
C           plume, FSUBZ3                                ---   CALL VRTSBL
           IF (ZR .LE. HPEN) THEN
              CALL VRTSBL (SZ3, MAX(0.0D0,HE3-HV), HPEN)
           ELSE
              CALL VRTSBN (SZ3, MAX(0.0D0,HE3-HV))
           END IF
           FSUBZ3 = FSUBZ
! Added for HIGHLY BUOYANT PLUME (HBP), JAN 2023
! Apply HBP only for identified sources
           IF ((HBPLUME) .AND. (HBPSRC(ISRC).EQ.'Y')) THEN
             IF (PPF .LT. 1.0D0) THEN
C  ************************************* modified code JAN 2023  --kja
               COUT = (QTK * (1.0D0-PPF) / UEFFD) * ( FSUBYD*FSUBZD ) +
     &                (QTK * (1.0D0-PPF) / UEFFN) * ( FSUBYN*FSUBZN ) +
     &                (QTK * PPF*PPFN / UEFF3) * ( FSUBY3*FSUBZ3 )
               SCOUT= (SQTK * (1.0D0-PPF) / UEFFD) * (FSUBYD*FSUBZD) +
     &                (SQTK * (1.0D0-PPF) / UEFFN) * (FSUBYN*FSUBZN) +
     &                (SQTK * PPF*PPFN / UEFF3) * ( FSUBY3*FSUBZ3 )
             ELSE
               COUT  = (QTK * PPF*PPFN / UEFF3) * ( FSUBY3*FSUBZ3 )
               SCOUT = (SQTK * PPF*PPFN / UEFF3) * ( FSUBY3*FSUBZ3 )
C  ************************************* modified code end  --kja
             ENDIF
           ELSE ! For DEFAULT plume bouyancy
             IF (PPF .LT. 1.0D0) THEN
               COUT = (QTK * (1.0D0-PPF) / UEFFD) * ( FSUBYD*FSUBZD ) +
     &                (QTK * (1.0D0-PPF) / UEFFN) * ( FSUBYN*FSUBZN ) +
     &                (QTK * PPF / UEFF3) * ( FSUBY3*FSUBZ3 )
               SCOUT= (SQTK * (1.0D0-PPF) / UEFFD) * (FSUBYD*FSUBZD) +
     &                (SQTK * (1.0D0-PPF) / UEFFN) * (FSUBYN*FSUBZN) +
     &                (SQTK * PPF / UEFF3) * ( FSUBY3*FSUBZ3 )

             ELSE
               COUT  = (QTK * PPF / UEFF3) * ( FSUBY3*FSUBZ3 )
               SCOUT = (SQTK * PPF / UEFF3) * ( FSUBY3*FSUBZ3 )
             ENDIF             
           ENDIF
! End HBP insert
         ELSE
            FSUBZ3 = 0.0D0
            HPEN   = 0.0D0
            COUT = (QTK / UEFFD) * ( FSUBYD * FSUBZD ) +
     &             (QTK / UEFFN) * ( FSUBYN * FSUBZN )
            SCOUT = (SQTK / UEFFD) * ( FSUBYD * FSUBZD ) +
     &             (SQTK / UEFFN) * ( FSUBYN * FSUBZN )

         END IF

      END IF

      RETURN
      END

C***********************************************************************
      SUBROUTINE SAERCALC( XARG, L_PLUME, AEROUT, SAEROUT )
C***********************************************************************
C             SAERCALC Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE: Calculates the AERMOD concentration without downwash
C
C        PROGRAMMER: Huy Tran; based on AERCALC (cal1.f) by Roger Brode, PES, Inc.
C
C        DATE:     September 15, 2000
C
C        INPUTS:   XARG         - Real - Distance (m), downwind for coherent
C                                        plume component (X) and radial for
C                                        random component (DISTR)
C                  L_PLUME      - Log  - Specifies coherent plume calculation
C                                        if TRUE, otherwise random component
C
C        OUTPUTS:  AEROUT(NTYP) - Real - AERMOD component of concentration
C                                        without building downwash for either
C                                        coherent plume component or for
C                                        random component, depending on
C                                        L_PLUME.
C
C        CALLED FROM:   PCALC
C
C***********************************************************************
C     Variable Declarations
      USE MAIN1
      IMPLICIT NONE
      CHARACTER MODNAM*12
      INTEGER :: J
      DOUBLE PRECISION :: AEROUT(NUMTYP), AERTMP(NUMTYP), FYOUT, XARG,
     &                    ADJ, FRAN, FRAN3, 
     &                    SAEROUT(NUMTYP),SAERTMP(NUMTYP)
      LOGICAL :: L_PLUME
C Unused:       INTEGER :: NDXZMID
C Unused:       DOUBLE PRECISION :: ZMID, VALABV, VBELOW

CMGS  Wood 8/5/2021 D063 Added to make plume rise adjustment due to PLATFORM
CMGS  External Functions:
      DOUBLE PRECISION, EXTERNAL  :: PLAT_GRADPLUMADJ

C     Variable Initializations
      MODNAM = 'SAERCALC'

C     Initialize AEROUT(NUMTYP) and AERTMP(NUMTYP) arrays
      AEROUT(:) = 0.0D0
      SAEROUT(:) = 0.0D0
      AERTMP(:) = 0.0D0
      SAERTMP(:) = 0.0D0

C --- Initialize FRAN and FRAN3
      FRAN  = 0.0D0
      FRAN3 = 0.0D0

      IF (XARG .LT. 1.0D0) THEN
C        Receptor is too close to source for calculation
C        or upwind of source for coherent plume component
         AEROUT(:) = 0.0D0
         SAEROUT(:) = 0.0D0
         RETURN

      END IF

C     Determine Deposition Correction Factors
      IF (NPD .EQ. 0 .AND. (LDGAS .OR. LWGAS)) THEN
         CALL PDEPG (XARG)
      ELSE
         DQCORG = 1.0D0
         WQCORG = 1.0D0
      END IF
      IF (NPD .GT. 0 .AND. (LDPART .OR. LWPART)) THEN
         CALL PDEP (XARG)
      ELSE IF (NPD .GT. 0) THEN
C        Set DQCOR(NPD) and WQCOR(NPD) arrays to 1.0
         DQCOR = 1.0D0
         WQCOR = 1.0D0
      END IF

C     Set initial effective parameters
      UEFF  = US
      SVEFF = SVS
      SWEFF = SWS
      TGEFF = TGS
      IF ( UNSTAB  .AND.  (HS .LT. ZI) ) THEN
         UEFFD  = US
         SVEFFD = SVS
         SWEFFD = SWS
         UEFFN  = US
         SVEFFN = SVS
         SWEFFN = SWS
         UEFF3  = US
         SVEFF3 = SVS
         SWEFF3 = SWS
         TGEFF3 = TGS
      END IF

CRJP  Add temporary debugging statement here.

C     ENSR ENHANCEMENT OF WRITE STATEMENT TO IDENTIFY COMPONENT CONCENTRATION
      IF(DEBUG) THEN
        WRITE(DBGUNT, 6014) SRCID(ISRC)
6014    FORMAT(//,'SRCID: ', A8)
        IF(L_PLUME)THEN
           WRITE(DBGUNT, 6015) UEFF, SVEFF, SWEFF
6015       FORMAT(//,'COHERENT PLUME COMPONENT',/,5X,
     &       'Initial effective parameters for ',
     &       'stable or direct convective ',
     &       'plume:',//,5x,'Ueff = ',F7.2,' m/s; ',
     &       'SVeff = ',F7.2,
     &       ' m/s; SWeff = ',F7.2,' m/s.',/)
        ELSE
           WRITE(DBGUNT, 6016) UEFF, SVEFF, SWEFF
6016       FORMAT(//,'MEANDER COMPONENT',/,5X,
     &       'Initial effective parameters for ',
     &       'stable or direct convective ',
     &       'plume:',//,5x,'Ueff = ',F7.2,' m/s; ',
     &       'SVeff = ',F7.2,
     &       ' m/s; SWeff = ',F7.2,' m/s.',/)
        END IF
      END IF

C     Define plume centroid height (CENTER) for use in
C     inhomogeneity calculations
      CALL CENTROID ( XARG )

C     Calculate the plume rise                     ---   CALL DELTAH
      CALL DELTAH ( XARG )

C     If the atmosphere is unstable and the stack
C     top is below the mixing height, calculate
C     the CBL PDF coefficients                     ---   CALL PDF
      IF( UNSTAB  .AND.  (HS .LT. ZI) ) THEN
         CALL PDF
      END IF

C     Determine Effective Plume Height             ---   CALL HEFF
      CALL HEFF ( XARG )
C ----------------------------------------------------------------------
C ---------- BEGIN PLATFORM DOWNWASH PLUME RISE ADJUSTMENT -------------
C ----------------------------------------------------------------------
C        Adjust plume rise if PLATFORM is present. This will decrease the
C        effective height when there is a PLATFORM due to downwash effects
C        of the structure on top of the platform. Calls PLAT_GRADPLUMEADJ,
C        PLAT_DOWNWASH, and PLAT_CUBIC (all contained in calc1.f).
C        Michelle G. Snyder, WOOD, 8/5/2021
CCRT     D063
         IF (OSPLAT(ISRC) .AND. (XARG > 0.0D0)) THEN
            IF (STABLE) THEN
               DHP = DHP - PLAT_GRADPLUMADJ( XARG )
            ELSE !modify direct plume downdraft (DHP1) only
               DHP1 = DHP1 - PLAT_GRADPLUMADJ( XARG )
            END IF

C           Get new effective heights based on adjusted plume rise from platform
            CALL HEFF ( XARG )

         END IF

C ----------------------------------------------------------------------
C --------- END OF PLATFORM DOWNWASH PLUME RISE ADJUSTMENT -------------
C----------             RESUME NORMAL PROCESSING           -------------
C ----------------------------------------------------------------------

C     Compute effective parameters using an
C     average through plume layer
      CALL IBLVAL ( XARG )

C     Call PDF & HEFF again for final CBL plume heights
      IF (UNSTAB .AND. (HS.LT.ZI) ) THEN
         CALL PDF
         CALL HEFF ( XARG )
      END IF

C     Determine Dispersion Parameters              ---   CALL PDIS
      CALL PDIS ( XARG )

C     Calculate the 'y-term' contribution to
C     dispersion, FSUBY
      IF (L_PLUME) THEN
         IF (L_EFFSIGY) THEN
C ---       Calculate fraction of random kinetic energy to total kinetic energy
C           for FASTALL option to optimize meander using effective sigma-y.
C           Note that these effective parameters are based on the downwind distance,
C           rather than the radial distance used in the standard meander approach
            IF (STABLE .OR. (UNSTAB.AND.(HS.GE.ZI))) THEN
               CALL MEANDR( UEFF, SVEFF, FRAN )
            ELSE IF (UNSTAB) THEN
               CALL MEANDR( UEFFD, SVEFFD, FRAN )
               IF (PPF .GT. 0.0D0) THEN
C                 For penetrated source calculate weighted average of
C                 direct/indirect plume component and penetrated component
                  CALL MEANDR( UEFF3, SVEFF3, FRAN3 )
                  FRAN = PPF*FRAN3 + (1.0D0-PPF)*FRAN
               END IF
            END IF

C           Calculate effective sigma-y for non-DFAULT FASTALL option (L_EFFSIGY)
            SYEFF = 1.0D0/((FRAN/(SRT2PI*XARG)) + (1.0D0-FRAN)/SY)
            IF (L_EFFSIGY) THEN
               IF (DABS(Y) .GT. NUMSYEFF*SYEFF) THEN
C                 Plume is more than 4 sigmas off centerline, skip calculation
                  FYOUT = 0.0D0
               ELSE
C                 Calculate FSUBY for coherent plume        ---   CALL FYPLM
                  CALL FYPLM(SYEFF,FYOUT)
               END IF
            ELSE
               IF (X.LT.1.0D0 .OR. DABS(Y) .GT. NUMSYEFF*SYEFF) THEN
C                 Receptor is upwind of source or more than 4 sigmas off
C                 ceenterline, skip calculation
                  FYOUT = 0.0D0
               ELSE
C                 Calculate FSUBY for coherent plume        ---   CALL FYPLM
                  CALL FYPLM(SYEFF,FYOUT)
               END IF
            END IF

         ELSE
C           Calculate FSUBY for coherent plume        ---   CALL FYPLM
            CALL FYPLM(SY,FYOUT)

         END IF

      ELSE
C        Calculate FSUBY for random component      ---   CALL FYPAN
         CALL FYPAN(FYOUT)

      END IF

      FSUBY  = FYOUT
      FSUBYD = FSUBY
      FSUBYN = FSUBYD

C     Calculate the 'y-term' contribution to dispersion
C     for the penetrated plume, FSUBY3
      IF( UNSTAB  .AND.  (HS .LT. ZI)  .AND. (PPF .GT. 0.0D0) )THEN
C        Compute meander fraction of horizontal distribution function
C        from Venky's memo of 6/24/98.
         IF (L_PLUME) THEN
            IF (L_EFFSIGY) THEN
C ---          Calculate fraction of random kinetic energy to total kinetic energy
C              for FASTALL option to optimize meander using effective sigma-y.
C              Note that these effective parameters are based on the downwind distance,
C              rather than the radial distance used in the standard meander approach
               SYEFF = 1.0D0/((FRAN/(SRT2PI*XARG)) + (1.0D0-FRAN)/SY3)

               IF (DABS(Y) .GT. NUMSYEFF*SYEFF) THEN
C                 Plume is more than 4 sigmas off centerline, skip calculation
                  FYOUT = 0.0D0
               ELSE
C                 Calculate FSUBY for coherent plume        ---   CALL FYPLM
                  CALL FYPLM(SYEFF,FYOUT)
               END IF

            ELSE
C              Calculate FSUBY for coherent plume        ---   CALL FYPLM
               CALL FYPLM(SY3,FYOUT)
            END IF
         ELSE
C           Calculate FSUBY for random component   ---   CALL FYPAN
            CALL FYPAN(FYOUT)
         END IF

         FSUBY3 = FYOUT
      ELSE
         FSUBY3 = 0.0D0
      END IF

C     Check for zero "y-terms"; if zero then skip calculations
C     and go to next receptor.
      IF( FSUBY.EQ.0.0D0 .AND. FSUBY3.EQ.0.0D0 )THEN
C        Set AEROUT(NUMTYP) array to 0.0
         AEROUT(:) = 0.0D0
         SAEROUT(:) = 0.0D0

      ELSE

         IF (NPD .EQ. 0) THEN
C           Perform calculations for gases
C           Assign plume tilt, HV = 0.0
            HV = 0.0D0

            ADJ = DQCORG * WQCORG

            IF (STABLE .OR. (UNSTAB.AND.(HS.GE.ZI))) THEN
C              Calculate height of the "effective reflecting surface"
               CALL REFL_HT (HE, XARG, SZB, 0.0D0, HSBL)
            ELSEIF ( UNSTAB ) THEN
               HSBL = 0.0D0
            END IF

            IF (UNSTAB .AND. (HS.LT.ZI) .AND. (PPF.GT.0.0D0)) THEN
C              Calculate height of the "effective reflecting surface"
               CALL REFL_HT (HE3, XARG, SZB3, 0.0D0, HPEN)
            ELSE
               HPEN = 0.0D0
            END IF

C           Determine the CRITical Dividing Streamline---   CALL CRITDS
            CALL CRITDS (HE)

C           Calculate the fraction of plume below
C           HCRIT, PHEE                               ---   CALL PFRACT
            CALL PFRACT (HE)

C           Calculate FOPT = f(PHEE)                  ---   CALL FTERM
            CALL FTERM

C           Calculate AERMOD Concentration     ---   CALL AER_PCHI
            CALL AER_SPCHI( XARG, ADJ, VDEPG, 0, AEROUT(:), SAEROUT(:) )

         ELSE
C           Perform calculations for particles, loop through particle sizes

C           Begin loop over particle sizes
            DO J = 1, NPD

C              Calculate Plume Tilt Due to Settling, HV
               HV = (XARG/US) * VGRAV(J)

C              Adjust Jth contribution by mass fraction and source
C              depletion
               ADJ = PHI(J) * DQCOR(J) * WQCOR(J)

               IF (STABLE .OR. (UNSTAB.AND.(HS.GE.ZI))) THEN
C                 Calculate height of the "effective reflecting surface"
C                 Calculate Settled Plume Height(s), HESETL
                  HESETL = MAX( 0.0D0, HE - HV )
                  CALL REFL_HT (HESETL, XARG, SZB, 0.0D0, HSBL)
               ELSEIF ( UNSTAB ) THEN
                  HESETL = MAX( 0.0D0, 0.5D0*(HED1+HED2) - HV )
                  HSBL = 0.0D0
               END IF

               IF (UNSTAB .AND. (HS.LT.ZI) .AND. (PPF.GT.0.0D0)) THEN
C                 Calculate height of the "effective reflecting surface"
C                 Calculate Settled Plume Height(s), HE3SETL
                  HE3SETL = MAX( 0.0D0, HE3 - HV )
                  CALL REFL_HT (HE3SETL, XARG, SZB3, 0.0D0, HPEN)
                  HPEN = MAX( HPEN, ZI )
               ELSE
                  HPEN = 0.0D0
               END IF

C              Determine the CRITical Dividing Streamline---   CALL CRITDS
               CALL CRITDS (HESETL)

C              Calculate the fraction of plume below
C              HCRIT, PHEE                               ---   CALL PFRACT
               CALL PFRACT (HESETL)

C              Calculate FOPT = f(PHEE)                  ---   CALL FTERM
               CALL FTERM

C              Calculate AERMOD Concentration            ---   CALL AER_PCHI
               CALL AER_SPCHI(XARG,ADJ,VDEP(J),J,AERTMP(:),SAERTMP(:))
               AEROUT(:) = AEROUT(:) + AERTMP(:)
               SAEROUT(:) = SAEROUT(:) + SAERTMP(:)

            END DO
C           End loop over particle sizes

         END IF

      END IF

      RETURN
      END
      
C***********************************************************************
      SUBROUTINE SUNAE( YEAR, DAY, lHOUR, LAT, LONG )
C***********************************************************************

c     Calculates the local solar azimuth and elevation angles, and
c     the distance to and angle subtended by the Sun, at a specific 
c     location and time using approximate formulas in The Astronomical 
c     Almanac.  Accuracy of angles is 0.01 deg or better (the angular 
c     width of the Sun is about 0.5 deg, so 0.01 deg is more than
c     sufficient for most applications).

c     Unlike many GCM (and other) sun angle routines, this
c     one gives slightly different sun angles depending on
c     the year.  The difference is usually down in the 4th
c     significant digit but can slowly creep up to the 3rd
c     significant digit after several decades to a century.

c     A refraction correction appropriate for the "US Standard
c     Atmosphere" is added, so that the returned sun position is
c     the APPARENT one.  The correction is below 0.1 deg for solar
c     elevations above 9 deg.  To remove it, comment out the code
c     block where variable REFRAC is set, and replace it with
c     REFRAC = 0.0.

c   References:

c     Michalsky, J., 1988: The Astronomical Almanac's algorithm for
c        approximate solar position (1950-2050), Solar Energy 40,
c        227-235 (but the version of this program in the Appendix
c        contains errors and should not be used)

c     The Astronomical Almanac, U.S. Gov't Printing Office, Washington,
c        D.C. (published every year): the formulas used from the 1995
c        version are as follows:
c        p. A12: approximation to sunrise/set times
c        p. B61: solar elevation ("altitude") and azimuth
c        p. B62: refraction correction
c        p. C24: mean longitude, mean anomaly, ecliptic longitude,
c                obliquity of ecliptic, right ascension, declination,
c                Earth-Sun distance, angular diameter of Sun
c        p. L2:  Greenwich mean sidereal time (ignoring T^2, T^3 terms)


c   Authors:  Dr. Joe Michalsky (joe@asrc.albany.edu)
c             Dr. Lee Harrison (lee@asrc.albany.edu)
c             Atmospheric Sciences Research Center
c             State University of New York
c             Albany, New York

c   Modified by:  Dr. Warren Wiscombe (wiscombe@climate.gsfc.nasa.gov)
c                 NASA Goddard Space Flight Center
c                 Code 913
c                 Greenbelt, MD 20771


c   WARNING:  Do not use this routine outside the year range
c             1950-2050.  The approximations become increasingly
c             worse, and the calculation of Julian date becomes
c             more involved.

c   Input:

c      YEAR     year (INTEGER; range 1950 to 2050)

c      DAY      day of year at LAT-LONG location (INTEGER; range 1-366)

c      HOUR     hour of DAY [GMT or UT] (REAL; range -13.0 to 36.0)
c               = (local hour) + (time zone number)
c                 + (Daylight Savings Time correction; -1 or 0)
C               where (local hour) range is 0 to 24,
c               (time zone number) range is -12 to +12, and
c               (Daylight Time correction) is -1 if on Daylight Time
c               (summer half of year), 0 otherwise;  
c               Example: 8:30 am Eastern Daylight Time would be
c
c                           HOUR = 8.5 + 5 - 1 = 12.5

c      LAT      latitude [degrees]
c               (REAL; range -90.0 to 90.0; north is positive)

c      LONG     longitude [degrees]
c               (REAL; range -180.0 to 180.0; east is positive)


c   Output:

c      AZ       solar azimuth angle (measured east from north,
c               0 to 360 degs)

c      EL       solar elevation angle [-90 to 90 degs]; 
c               solar zenith angle = 90 - EL

c      SOLDIA   solar diameter [degs]

c      SOLDST   distance to sun [Astronomical Units, AU]
c               (1 AU = mean Earth-sun distance = 1.49597871E+11 m
c                in IAU 1976 System of Astronomical Constants)


c   Local Variables:

c     DEC       Declination (radians)

c     ECLONG    Ecliptic longitude (radians)

c     GMST      Greenwich mean sidereal time (hours)

c     HA        Hour angle (radians, -pi to pi)

c     JD        Modified Julian date (number of days, including 
c               fractions thereof, from Julian year J2000.0);
c               actual Julian date is JD + 2451545.0

c     LMST      Local mean sidereal time (radians)

c     MNANOM    Mean anomaly (radians, normalized to 0 to 2*pi)

c     MNLONG    Mean longitude of Sun, corrected for aberration 
c               (deg; normalized to 0-360)

c     OBLQEC    Obliquity of the ecliptic (radians)

c     RA        Right ascension  (radians)

c     REFRAC    Refraction correction for US Standard Atmosphere (degs)

c --------------------------------------------------------------------
c   Uses double precision for safety and because Julian dates can
c   have a large number of digits in their full form (but in practice
c   this version seems to work fine in single precision if you only
c   need about 3 significant digits and aren't doing precise climate
c   change or solar tracking work).
c --------------------------------------------------------------------

c   Why does this routine require time input as Greenwich Mean Time 
c   (GMT; also called Universal Time, UT) rather than "local time"?
c   Because "local time" (e.g. Eastern Standard Time) can be off by
c   up to half an hour from the actual local time (called "local mean
c   solar time").  For society's convenience, "local time" is held 
c   constant across each of 24 time zones (each technically 15 longitude
c   degrees wide although some are distorted, again for societal 
c   convenience).  Local mean solar time varies continuously around a 
c   longitude belt;  it is not a step function with 24 steps.  
c   Thus it is far simpler to calculate local mean solar time from GMT,
c   by adding 4 min for each degree of longitude the location is
c   east of the Greenwich meridian or subtracting 4 min for each degree
c   west of it.  

c --------------------------------------------------------------------

c   TIME
c   
c   The measurement of time has become a complicated topic.  A few
c   basic facts are:
c   
c   (1) The Gregorian calendar was introduced in 1582 to replace 
c   Julian calendar; in it, every year divisible by four is a leap 
c   year just as in the Julian calendar except for centurial years
c   which must be exactly divisible by 400 to be leap years.  Thus 
c   2000 is a leap year, but not 1900 or 2100.

c   (2) The Julian day begins at Greenwich noon whereas the calendar 
c   day begins at the preceding midnight;  and Julian years begin on
c   "Jan 0" which is really Greenwich noon on Dec 31.  True Julian 
c   dates are a continous count of day numbers beginning with JD 0 on 
c   1 Jan 4713 B.C.  The term "Julian date" is widely misused and few
c   people understand it; it is best avoided unless you want to study
c   the Astronomical Almanac and learn to use it correctly.

c   (3) Universal Time (UT), the basis of civil timekeeping, is
c   defined by a formula relating UT to GMST (Greenwich mean sidereal
c   time).  UTC (Coordinated Universal Time) is the time scale 
c   distributed by most broadcast time services.  UT, UTC, and other
c   related time measures are within a few sec of each other and are
c   frequently used interchangeably.

c   (4) Beginning in 1984, the "standard epoch" of the astronomical
c   coordinate system is Jan 1, 2000, 12 hr TDB (Julian date 
c   2,451,545.0, denoted J2000.0).  The fact that this routine uses
c   1949 as a point of reference is merely for numerical convenience.
c --------------------------------------------------------------------

C     Variable Declarations
      USE MAIN1
      USE GRS7MOD
      IMPLICIT NONE
      CHARACTER MODNAM*12

c     .. Scalar Arguments ..

      INTEGER::   YEAR, DAY, lHOUR
      DOUBLE PRECISION:: AZ, EL, HOUR, LAT, LONG,
     &                   toffset
c     ..
c     .. Local Scalars ..

      LOGICAL   PASS1
      INTEGER   DELTA, LEAP
      DOUBLE PRECISION  DEC, DEN, ECLONG, GMST, HA, JD, LMST,
     &                  MNANOM, MNLONG, NUM, OBLQEC, RA,
     &                  RPD, REFRAC, TIME
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC AINT, ASIN, ATAN, COS, MOD, SIN, TAN
c     ..
c     .. Data statements ..

      SAVE     PASS1, RPD
      DATA     PASS1 /.True./
c     ..
C     Variable Initializations
      MODNAM = 'SUNAE'

c     write(*,*)"YEAR=",YEAR,"JDAY=",JDAY,
c    &          "HR=",lHOUR,
c    &          "LAT=",LAT,"LONG=",LONG
      IF (YEAR < 2000)YEAR = YEAR + 2000
      IF( YEAR.LT.1950 .OR. YEAR.GT.2050 ) 
     &    STOP 'SUNAE--bad input variable YEAR'
      IF( DAY.LT.1 .OR. DAY.GT.366 ) 
     &    STOP 'SUNAE--bad input variable DAY'
      IF( HOUR.LT.-13.0 .OR. HOUR.GT.36.0 ) 
     &    STOP 'SUNAE--bad input variable HOUR'
      IF( LAT.LT.-90.0 .OR. LAT.GT.90.0 ) 
     &    STOP 'SUNAE--bad input variable LAT'
      IF( LONG.LT.-180.0 .OR. LONG.GT.180.0 ) 
     &    STOP 'SUNAE--bad input variable LONG'

      IF(PASS1) THEN
c        PI  = 2.*ASIN( 1.0 )
c        TWOPI  = 2.*PI
         RPD    = PI / 180.
         PASS1 = .False.
      ENDIF

c     Calculate Time Zone offset to convert locatime to UTC
      toffset = floor(LONG/15)
      HOUR = lHOUR*1.0D0 - toffset
      IF (HOUR > 24) THEN
            DAY = DAY + 1
            HOUR = HOUR -24
      ENDIF

c     ** current Julian date (actually add 2,400,000 
c     ** for true JD);  LEAP = leap days since 1949;
c     ** 32916.5 is midnite 0 jan 1949 minus 2.4e6

      DELTA  = YEAR - 1949
      LEAP   = DELTA / 4
      JD     = 32916.5 + (DELTA*365 + LEAP + DAY) + HOUR / 24.

c     ** last yr of century not leap yr unless divisible
c     ** by 400 (not executed for the allowed YEAR range,
c     ** but left in so our successors can adapt this for 
c                    ** the following 100 years)

      IF( MOD( YEAR, 100 ).EQ.0 .AND.
     &    MOD( YEAR, 400 ).NE.0 ) JD = JD - 1.

c                     ** ecliptic coordinates
c                     ** 51545.0 + 2.4e6 = noon 1 jan 2000

      TIME  = JD - 51545.0

c                    ** force mean longitude between 0 and 360 degs

      MNLONG = 280.460 + 0.9856474*TIME
      MNLONG = MOD( MNLONG, 360. )
      IF( MNLONG.LT.0. ) MNLONG = MNLONG + 360.

c                    ** mean anomaly in radians between 0 and 2*pi

      MNANOM = 357.528 + 0.9856003*TIME
      MNANOM = MOD( MNANOM, 360. )
      IF( MNANOM.LT.0.) MNANOM = MNANOM + 360.

      MNANOM = MNANOM*RPD

c                    ** ecliptic longitude and obliquity 
c                    ** of ecliptic in radians

      ECLONG = MNLONG + 1.915*SIN( MNANOM ) + 0.020*SIN( 2.*MNANOM )
      ECLONG = MOD( ECLONG, 360. )
      IF( ECLONG.LT.0. ) ECLONG = ECLONG + 360.

      OBLQEC = 23.439 - 0.0000004*TIME
      ECLONG = ECLONG*RPD
      OBLQEC = OBLQEC*RPD

c                    ** right ascension

      NUM  = COS( OBLQEC )*SIN( ECLONG )
      DEN  = COS( ECLONG )
      RA   = ATAN( NUM / DEN )

c                    ** Force right ascension between 0 and 2*pi

      IF( DEN.LT.0.0 ) THEN
         RA  = RA + PI
      ELSE IF( NUM.LT.0.0 ) THEN
         RA  = RA + TWOPI
      END IF

c                    ** declination

      DEC  = ASIN( SIN( OBLQEC )*SIN( ECLONG ) )

c                    ** Greenwich mean sidereal time in hours

      GMST = 6.697375 + 0.0657098242*TIME + HOUR

c                    ** Hour not changed to sidereal time since 
c                    ** 'time' includes the fractional day

      GMST  = MOD( GMST, 24. )
      IF( GMST.LT.0. ) GMST   = GMST + 24.

c                    ** local mean sidereal time in radians

      LMST  = GMST + LONG / 15.
      LMST  = MOD( LMST, 24. )
      IF( LMST.LT.0. ) LMST   = LMST + 24.

      LMST   = LMST*15.*RPD

c                    ** hour angle in radians between -pi and pi

      HA  = LMST - RA

      IF( HA.LT.- PI ) HA  = HA + TWOPI
      IF( HA.GT.PI )   HA  = HA - TWOPI

c                    ** solar azimuth and elevation

      EL  = ASIN( SIN( DEC )*SIN( LAT*RPD ) +
     &            COS( DEC )*COS( LAT*RPD )*COS( HA ) )

      AZ  = ASIN( - COS( DEC )*SIN( HA ) / COS( EL ) )

c                    ** Put azimuth between 0 and 2*pi radians

      IF( SIN( DEC ) - SIN( EL )*SIN( LAT*RPD ).GE.0. ) THEN

         IF( SIN(AZ).LT.0.) AZ  = AZ + TWOPI

      ELSE

         AZ  = PI - AZ

      END IF

c                     ** Convert elevation and azimuth to degrees
      EL  = EL / RPD
      AZ  = AZ / RPD

c  ======== Refraction correction for U.S. Standard Atmos. ==========
c      (assumes elevation in degs) (3.51823=1013.25 mb/288 K)

      IF( EL.GE.19.225 ) THEN

         REFRAC = 0.00452*3.51823 / TAN( EL*RPD )

      ELSE IF( EL.GT.-0.766 .AND. EL.LT.19.225 ) THEN

         REFRAC = 3.51823 * ( 0.1594 + EL*(0.0196 + 0.00002*EL) ) /
     &            ( 1. + EL*(0.505 + 0.0845*EL) )

      ELSE IF( EL.LE.-0.766 ) THEN

         REFRAC = 0.0

      END IF

      EL  = EL + REFRAC
c ===================================================================

c                   ** distance to sun in A.U. & diameter in degs

      SOLDST = 1.00014 - 0.01671*COS(MNANOM) - 0.00014*COS( 2.*MNANOM )
      SOLDIA = 0.5332 / SOLDST
      SOLELV = EL

      IF( EL.LT.-90.0 .OR. EL.GT.90.0 )
     &    STOP 'SUNAE--output argument EL out of range'
      IF( AZ.LT.0.0 .OR. AZ.GT.360.0 )
     &    STOP 'SUNAE--output argument AZ out of range'

      RETURN

      END