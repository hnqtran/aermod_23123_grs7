      SUBROUTINE GRS7_CALC
C***********************************************************************
C        GRSM_CALC Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE: Processes Hourly Results for GRSM Option
C
C        PROGRAMMER: CERC
C
C        DATE:    November 2020
C
C        MODIFIED:  CERC adjustments to calculation of instantaneous
C        plume spread (PLUMESIZES), building effects on
C        plume spread, multiple plume effects on plume spread, and
C        setting 0.1 m minimum distance for downstream receptors and
C        dispersion time for going between pancake and gaussian plumes.
C        See CERC technical memo dated January 17, 2023.
C
C        INPUTS:
C
C
C        OUTPUTS:
C
C
C        CALLED FROM:   HRLOOP, EV_LOOP, MAXDCALC
C***********************************************************************
      USE MAIN1
      USE GRS7MOD
      IMPLICIT NONE
      
C --- Local variables
      INTEGER::NUMGRP_USE, IGRP_USE, nStat
      CHARACTER MODNAM*12, NightStr*3
      DOUBLE PRECISION, ALLOCATABLE:: PLUMEBKGNO2(:), PLUMEBKGNO(:),!Plume background concentrations
     &                                PLUMEBKGO3(:),                ! Plume background Ozone concentration
     &                                PLUMEBKGROC(:),               ! Plume background ROC concentration
     &                                PLUMEBKGRP(:)                 ! Plume background RP concentration
      DOUBLE PRECISION, ALLOCATABLE:: NOXPERSOURCE(:),              ! NOX concentration for each source
     &                                VOCPERSOURCE(:),              ! VOC concentration for each source 
     &                                ROCPERSOURCE(:),              ! ROC concentration for each source 
     &                                SRCNO2(:),SRCNO(:)            ! NO2 and NO concentration for each source
      DOUBLE PRECISION, ALLOCATABLE:: SIGYSIGZ_I_X(:), SIGYSIGZ_E_X(:) !Instantaneous and ensemble plume size at X
      DOUBLE PRECISION:: SIGYSIGZ_I_0  !Instantaneous plume size at origin
      DOUBLE PRECISION:: MINBKGNO2, MINBKGNO, MINBKGNOX, MINBKGO3    !Minimum background concentrations
      DOUBLE PRECISION:: MINBKGROC, MINBKGRP, MINBKGSGN, MINBKGSNGN  !Minimum background concentrations of ROC, RP, SGN, SNG
      DOUBLE PRECISION:: CONCMAXNOX, INSTCONCNOX                     !Maximum and total NOX concentrations
      DOUBLE PRECISION:: CONCMAXVOC, INSTCONCVOC                     !Maximum and total VOC concentrations
      DOUBLE PRECISION:: CONCSOURCENO2, CONCSOURCENO                 !Concentrations NO2 and NO for a single source
      DOUBLE PRECISION:: CONCSOURCEVOC, CONCSOURCEROC                !Concentrations VOC and ROC for a single source
      DOUBLE PRECISION:: ENSCONCNO2, ENSCONCNO, INSTCONCNO2, INSTCONCNO !Total ensemble and instantaneous concentrations
      DOUBLE PRECISION:: ENSCONCROC,  INSTCONCROC                    !Total ensemble and instantaneous concentrations of ROC
      DOUBLE PRECISION:: ENSCONCRP,   INSTCONCRP                     !Total ensemble and instantaneous concentrations of RP
      DOUBLE PRECISION:: ENSCONCSGN,  INSTCONSGN                     !Total ensemble and instantaneous concentrations of SGN
      DOUBLE PRECISION:: ENSCONCSNGN, INSTCONSNGN                    !Total ensemble and instantaneous concentrations of SG
      DOUBLE PRECISION:: ROCALLSRCS, VOCALLSRCS                      !Total ensemble ROC and VOC concentration from all sources
      DOUBLE PRECISION:: NO2Frac                                     !NO2 fraction
      DOUBLE PRECISION:: ROCFrac                                     ! VOC to ROC ratio
      DOUBLE PRECISION:: BGCONC_IN, BGCONC_OUT, BGCONC_PPB
      LOGICAL::L_DoFullConv
      DOUBLE PRECISION:: fDum, MaxNOx, TotalNOx, fMultPlmFac
      DOUBLE PRECISION:: MaxROC, TotalROC
      DOUBLE PRECISION:: fBldFacMaxMinusMin,fMultPlmFacMaxMinusMin
      
      
C --- Constants
      DOUBLE PRECISION, PARAMETER::fDownstrmStart = 0.1D0
      DOUBLE PRECISION, PARAMETER::fBldFacMin = 0.0D0
      DOUBLE PRECISION, PARAMETER::fBldFacMax = 0.5D0
      DOUBLE PRECISION, PARAMETER::fMultPlmFacMin = 0.25D0
      DOUBLE PRECISION, PARAMETER::fMultPlmFacMax = 0.75D0
      
C --- Initialise
      MODNAM        = 'GRS7_CALC'
      L_DoFullConv  = .FALSE.
      HRVAL         = 0.0D0
      SIGYSIGZ_I_0  = 0.0D0
      CONCSOURCENO2 = 0.0
      CONCSOURCENO  = 0.0
      CONCSOURCEVOC = 0.0
      CONCSOURCEROC = 0.0
      CONCMAXNOX    = 0.0
      CONCMAXVOC    = 0.0
      MINBKGNO2     = 0.0
      MINBKGNO      = 0.0
      MINBKGNOX     = 0.0      
      MINBKGO3      = 0.0
      MINBKGROC     = 0.0
      MINBKGRP      = 0.0
      MINBKGSGN     = 0.0
      MINBKGSNGN    = 0.0
      fBldFacMaxMinusMin = fBldFacMax - fBldFacMin
      fMultPlmFacMaxMinusMin = fMultPlmFacMax - fMultPlmFacMin
          
      IF(ALLOCATED(PLUMEBKGNO2)) DEALLOCATE(PLUMEBKGNO2,STAT=nStat)
      ALLOCATE(PLUMEBKGNO2(NUMSRC),STAT=nStat)
      IF(nStat/=0)THEN
        CALL ERRHDL(PATH, MODNAM, 'E','615','GRS7') 
      END IF
      IF(ALLOCATED(PLUMEBKGNO)) DEALLOCATE(PLUMEBKGNO,STAT=nStat)
      ALLOCATE(PLUMEBKGNO(NUMSRC),STAT=nStat)
      IF(nStat/=0)THEN
        CALL ERRHDL(PATH, MODNAM, 'E','615','GRS7') 
      END IF
      IF(ALLOCATED(PLUMEBKGO3)) DEALLOCATE(PLUMEBKGO3,STAT=nStat)
      ALLOCATE(PLUMEBKGO3(NUMSRC),STAT=nStat)
      IF(nStat/=0)THEN
        CALL ERRHDL(PATH, MODNAM, 'E','615','GRS7') 
      END IF
      IF(ALLOCATED(NOXPERSOURCE)) DEALLOCATE(NOXPERSOURCE,STAT=nStat)
      ALLOCATE(NOXPERSOURCE(NUMSRC),STAT=nStat)
      IF(nStat/=0)THEN
        CALL ERRHDL(PATH, MODNAM, 'E','615','GRSM') 
      END IF
      IF(ALLOCATED(PLUMEBKGROC)) DEALLOCATE(PLUMEBKGROC,STAT=nStat)
      ALLOCATE(PLUMEBKGROC(NUMSRC),STAT=nStat)
      IF(nStat/=0)THEN
        CALL ERRHDL(PATH, MODNAM, 'E','615','GRS7') 
      END IF
      IF(ALLOCATED(PLUMEBKGRP)) DEALLOCATE(PLUMEBKGRP,STAT=nStat)
      ALLOCATE(PLUMEBKGRP(NUMSRC),STAT=nStat)
      IF(nStat/=0)THEN
        CALL ERRHDL(PATH, MODNAM, 'E','615','GRS7') 
      END IF
      IF(ALLOCATED(VOCPERSOURCE)) DEALLOCATE(VOCPERSOURCE,STAT=nStat)
      ALLOCATE(VOCPERSOURCE(NUMSRC),STAT=nStat)
      IF(nStat/=0)THEN
        CALL ERRHDL(PATH, MODNAM, 'E','615','GRS7') 
      END IF
      IF(ALLOCATED(ROCPERSOURCE)) DEALLOCATE(ROCPERSOURCE,STAT=nStat)
      ALLOCATE(ROCPERSOURCE(NUMSRC),STAT=nStat)
      IF(nStat/=0)THEN
        CALL ERRHDL(PATH, MODNAM, 'E','615','GRS7') 
      END IF
      IF(ALLOCATED(SRCNO2)) DEALLOCATE(SRCNO2,STAT=nStat)
      ALLOCATE(SRCNO2(NUMSRC),STAT=nStat)
      IF(nStat/=0)THEN
        CALL ERRHDL(PATH, MODNAM, 'E','615','GRS7') 
      END IF
      IF(ALLOCATED(SRCNO)) DEALLOCATE(SRCNO,STAT=nStat)
      ALLOCATE(SRCNO(NUMSRC),STAT=nStat)
      IF(nStat/=0)THEN
        CALL ERRHDL(PATH, MODNAM, 'E','615','GRS7') 
      END IF
      
      IF(INLNDEBUG)write(*,*)'<------------------------->'
      write(*,'(A,I10,5A,F9.3)')"TIME= ",FULLDATE,
     &     " LAT=",ALAT," LON=",ALON,
     &     "|QSW: ",QSW

      PLUMEBKGNO2   = 0.0
      PLUMEBKGNO    = 0.0
      PLUMEBKGO3    = 0.0
      PLUMEBKGROC   = 0.0
      NOXPERSOURCE  = 0.0
      VOCPERSOURCE  = 0.0
      
      !Calculate NOx reaction rates
      CALL RXNRATES_GRS7
      IF(INLNDEBUG)write(*,'(A,7(A,D9.3))')"RXNRATES_GRS7:",
     &     " | R1: ",R1,
     &     " | R2: ",R2,
     &     " | R3: ",R3,
     &     " | R4: ",R4,
     &     " | R5: ",R5,
     &     " | R6: ",R6,
     &     " | R7: ",R7

      !Determine if its a night-time hour
      IF(QSW>0.0D0)THEN  
        L_NightHour = .FALSE.
        NightStr    = 'NO'
      ELSE
        L_NightHour = .TRUE.
        NightStr    = 'YES'
      END IF

C     This is background NO2      
      IF(EVONLY)THEN
        BGCONC_IN = EV_BGCONC(IHOUR)
      ELSE
        BGCONC_IN = BGCONC
      END IF

C --- Initialise concs by converting from ug m-3 to ppb.  Don't want to alter input bgd concs
      NO2CONC_BG = BGCONC_IN*NO2_PPB     !This is NO2 bgd
      IF(NOXMISS.OR.L_CalcNOxFromNO2)THEN
        !No NOx bgd
      ELSE
        NOXCONC_BG = NOXBGCONC*NO2_PPB   !This is NOx bgd (NOx as NO2); converted from ug/m3 back to ppb         
      END IF

      ROCCONC_BG = ROCBGCONC*VOC_PPB     !This is ROC bgd; converted back from ug/m3 to ppb
C     RPCONC_BG  = 0.0D0                 ! Keep background RP concentration zero for now      
C --- Calculate the upwind concentrations using the given background
C     values. Assumes some type of equilibrium
      IF(.NOT.O3MISS)THEN
        !If O3 bgd is missing assume full conversion later
        O3CONC_BG  = O3CONC/O3_PPB      !This is O3 bgd
        IF(INLNDEBUG)write(*,'(A,6(A,F6.3))')"pre-GETBGDEQUIL:",
     &          " | NOXCONC_BG=",NOXCONC_BG,
     &          " | NO2CONC_BG=",NO2CONC_BG,
     &          " | NOCONC_BG=",NOCONC_BG,
     &          " | O3CONC_BG=",O3CONC_BG,
     &          " | ROCCONC_BG=",ROCCONC_BG,
     &          " | RPCONC_BG=",RPCONC_BG
        IF (ROCCONC_BG > 0.0D0) THEN
            CALL GETBGDEQUIL_GRS7
        ELSE    
            CALL GETBGDEQUIL_GRSM
        END IF
        IF(INLNDEBUG)write(*,'(A,6(A,F6.3))')"pst-GETBGDEQUIL:",
     &          " | NOXCONC_BG=",NOXCONC_BG,
     &          " | NO2CONC_BG=",NO2CONC_BG,
     &          " | NOCONC_BG=",NOCONC_BG,
     &          " | O3CONC_BG=",O3CONC_BG,
     &          " | ROCCONC_BG=",ROCCONC_BG,
     &          " | RPCONC_BG=",RPCONC_BG
      ELSE
        O3CONC_BG  = -999.
        L_DoFullConv=.TRUE.    
      END IF

      IF(EVONLY)THEN
        !Event only - 1 grp
        NUMGRP_USE=1
      ELSE
        !Normal calc - use actual number of groups
        NUMGRP_USE=NUMGRP
      END IF
      
C --- Begin Source Group LOOP to sum values
      DO IGRP = 1, NUMGRP_USE
        
        IF(EVONLY)THEN
          !Event only
          IGRP_USE=IDXEV(IEVENT)
        ELSE
          !Normal calculation use actual group index
          IGRP_USE=IGRP
        END IF          

        DO IREC = 1, NUMREC            
          IF(.NOT.L_DoFullConv)THEN
            !Not doing full conversion, doing proper chemistry          
            INSTCONCNO2 = 0.0
            INSTCONCNO  = 0.0         
            INSTCONCNOX = 0.0
            CONCMAXNOX  = 0.0
            VOCALLSRCS  = 0.0
            ROCALLSRCS  = 0.0
            INSTCONCROC = 0.0
            INSTCONCRP  = 0.0
            INSTCONSGN  = 0.0
            INSTCONSNGN = 0.0
            CONCMAXVOC  = 0.0
            
            !Allocate variables for dilution and entrainment
            IF(ALLOCATED(SIGYSIGZ_I_X)) DEALLOCATE(SIGYSIGZ_I_X,
     +                                                     STAT=nStat)
            IF(ALLOCATED(SIGYSIGZ_E_X)) DEALLOCATE(SIGYSIGZ_E_X,
     +                                                     STAT=nStat)
            ALLOCATE(SIGYSIGZ_I_X(NUMSRC),STAT=nStat)
            IF(nStat/=0)THEN
              CALL ERRHDL(PATH, MODNAM, 'E','615','GRS7') 
            END IF
            ALLOCATE(SIGYSIGZ_E_X(NUMSRC),STAT=nStat)
            IF(nStat/=0)THEN
              CALL ERRHDL(PATH, MODNAM, 'E','615','GRS7') 
            END IF
            SIGYSIGZ_I_X(:) = 0.0D0
            SIGYSIGZ_E_X(:) = 0.0D0
            
            !Calculate total NOx and highest single contribution to total NOx
            !Do the same think for ROC
            MaxNOx = 0.0D0
            MaxROC = 0.0D0
            TotalNOx = 0.0D0
            TotalROC = 0.0D0
            DO ISRC = 1, NUMSRC
              IF (IGROUP(ISRC,IGRP_USE) == 0) CYCLE
              TotalNOx = TotalNOx + CHI(IREC,ISRC,1)
              IF(CHI(IREC,ISRC,1) > MaxNOx) MaxNOx = CHI(IREC,ISRC,1)
              TotalROC = TotalROC + SCHI(IREC,ISRC,1)
              IF(SCHI(IREC,ISRC,1) > MaxROC) MaxROC = SCHI(IREC,ISRC,1)
            ENDDO
            IF(TotalNOx > 0.0D0)THEN
              fMultPlmFac = 1.D0 - MaxNOx/TotalNOx
            ELSE
              fMultPlmFac = 0.0D0
            ENDIF

            DO ISRC = 1, NUMSRC
              NOXPERSOURCE(ISRC) = 0.0D0
              VOCPERSOURCE(ISRC) = 0.0D0
              !Is src included in group?
              IF (IGROUP(ISRC,IGRP_USE) == 0) CYCLE  
              IF ( CHI(IREC,ISRC,1) == 0.0D0) CYCLE
C ---         Get concentrations
              CONCSOURCENO2 = 
     &              ANO2_RATIO(ISRC)*CHI(IREC,ISRC,1)*NO2_PPB         !This is the NO2 conc, ppb
              CONCSOURCENO  = 
     &             (1.0D0-ANO2_RATIO(ISRC))*CHI(IREC,ISRC,1)*NO2_PPB  !This is the NO conc, ppb
              CONCSOURCEVOC = 
     &                         1.0D0*SCHI(IREC,ISRC,1)*VOC_PPB        !This is the VOC conc, ppb
              CONCSOURCEROC = 
     &              AROC_RATIO(ISRC)*SCHI(IREC,ISRC,1)*VOC_PPB        !This is the ROC conc, ppb    
                
              !Find downwind distance to current receptor
              CALL SETSRC
              WDSIN = AWDSIN(ISRC)  ! Sine of wind. dir
              WDCOS = AWDCOS(ISRC)  ! Cosine of wind dir. from which wind is blowing
              IF (EVONLY) THEN
                CALL XYDIST(IEVENT)
              ELSE
                CALL XYDIST(IREC)   ! Call XYDIST from calc2.f, outputs are X, Y, DISTR
                                    ! X = downwind distant (m), Y = cross wind distance (m)
                                    ! DISTR = Source - receptor distance (m) = sqrt(X^2 + Y^2)
              END IF
      
              !Implement dilution and entrainment for downstream receptors
              IF(X > fDownstrmStart)THEN
                
                !Find the sizes of the instantaneous plume at origin (source)
                !Call PLUMESIZES subroutine in grms.f
                !Output: SIGYSIGZ_I,SIGYSIGZ_E - instantaneous and ensemble plume size
                CALL PLUMESIZES(0.0D0,SIGYSIGZ_I_0,fDum)                                                         
              
                !Find the sizes of the plume at X
                CALL PLUMESIZES(X,SIGYSIGZ_I_X(ISRC),SIGYSIGZ_E_X(ISRC))
                
                !Ensure I_0 <= I_X <= E_X
                IF(SIGYSIGZ_I_X(ISRC) > SIGYSIGZ_E_X(ISRC))
     +                           SIGYSIGZ_I_X(ISRC) = SIGYSIGZ_E_X(ISRC)
                IF(SIGYSIGZ_I_0 > SIGYSIGZ_I_X(ISRC))
     +                           SIGYSIGZ_I_0 = SIGYSIGZ_I_X(ISRC)
                
C***************************************************************************************
C*************** BEGIN BUILDING EFFECTS ADJUSTMENT *************************************
C*************** (Comment out this block to remove adjustment) *************************
C***************************************************************************************                
                !Tend calculated instantaneous plume spread back towards ensemble
                !plume spread in regions that are significnatly buildings-affected
C!              UNC-IE Note: BLDFAC calculated from subroutine PCAL in calc1.f
C!              IF(HRVAL(1) > 0.0D0)THEN
C!                  BLDFAC(IREC,ISRC) = 
C!    &                          GAMFACT*(PRMVAL(1)-PRMVAL_Src1)/HRVAL(1)
C!              ENDIF
C!              End UNC-IE Note
                IF(BLDFAC(IREC,ISRC) <= fBldFacMin)THEN
                  !Do nothing
                ELSEIF(BLDFAC(IREC,ISRC) >= fBldFacMax)THEN
                   SIGYSIGZ_I_X(ISRC) = SIGYSIGZ_E_X(ISRC)
                ELSE
                   SIGYSIGZ_I_X(ISRC) = SIGYSIGZ_I_X(ISRC) +
     &               (SIGYSIGZ_E_X(ISRC) - SIGYSIGZ_I_X(ISRC))*
     &               (BLDFAC(IREC,ISRC)-fBldFacMin)/fBldFacMaxMinusMin
                ENDIF
C***************************************************************************************            
C*************** END BUILDING EFFECTS ADJUSTMENT ***************************************
C***************************************************************************************              
                
                !Also tend instantaneous plume spread back towards ensemble
                !plume spread in regions with multiple significant plumes
                IF(fMultPlmFac <= fMultPlmFacMin)THEN
                  !Do nothing
                ELSEIF(fMultPlmFac >= fMultPlmFacMax)THEN
                   SIGYSIGZ_I_X(ISRC) = SIGYSIGZ_E_X(ISRC)
                ELSE
                   SIGYSIGZ_I_X(ISRC) = SIGYSIGZ_I_X(ISRC) +
     &               (SIGYSIGZ_E_X(ISRC) - SIGYSIGZ_I_X(ISRC))*
     &               (fMultPlmFac-fMultPlmFacMin)/fMultPlmFacMaxMinusMin
                ENDIF
         
                !Get the instantaneous plume concentration
                CONCSOURCENO2= CONCSOURCENO2*
     +                         SIGYSIGZ_E_X(ISRC)/(SIGYSIGZ_I_X(ISRC))
                CONCSOURCENO = CONCSOURCENO*
     +                         SIGYSIGZ_E_X(ISRC)/(SIGYSIGZ_I_X(ISRC))
                CONCSOURCEVOC= CONCSOURCEVOC*
     +                         SIGYSIGZ_E_X(ISRC)/(SIGYSIGZ_I_X(ISRC))
                CONCSOURCEROC= CONCSOURCEROC*
     +                         SIGYSIGZ_E_X(ISRC)/(SIGYSIGZ_I_X(ISRC))

                !Calculate the entrained concentration
                IF(GRP_BACK(IGRP_USE).AND..NOT.L_NIGHTHOUR)THEN
                  PLUMEBKGNO2(ISRC)= (1.0D0-
     +                     SIGYSIGZ_I_0/SIGYSIGZ_I_X(ISRC))*NO2CONC_BG
                  PLUMEBKGNO(ISRC) = (1.0D0-
     +                     SIGYSIGZ_I_0/SIGYSIGZ_I_X(ISRC))*NOCONC_BG
                  PLUMEBKGROC(ISRC) = (1.0D0-
     +                     SIGYSIGZ_I_0/SIGYSIGZ_I_X(ISRC))*ROCCONC_BG
                  PLUMEBKGRP(ISRC) = (1.0D0-
     +                     SIGYSIGZ_I_0/SIGYSIGZ_I_X(ISRC))*RPCONC_BG
                ELSE
                  PLUMEBKGNO2(ISRC) = 0.0D0
                  PLUMEBKGNO(ISRC)  = 0.0D0
                  PLUMEBKGROC(ISRC) = 0.0D0
                  PLUMEBKGRP(ISRC)  = 0.0D0
                END IF
                
                PLUMEBKGO3(ISRC)=(1.0D0-SIGYSIGZ_I_0/
     +                           SIGYSIGZ_I_X(ISRC))*O3CONC_BG
              ELSE
                !Upwind fully entrains background
                IF(GRP_BACK(IGRP_USE).AND..NOT.L_NIGHTHOUR)THEN
                  PLUMEBKGNO2(ISRC) = NO2CONC_BG
                  PLUMEBKGNO(ISRC)  = NOCONC_BG
                  PLUMEBKGROC(ISRC) = ROCCONC_BG     
                  PLUMEBKGRP(ISRC)  = RPCONC_BG      
                ELSE
                  PLUMEBKGNO2(ISRC) = 0.0D0
                  PLUMEBKGNO(ISRC)  = 0.0D0
                  PLUMEBKGROC(ISRC) = 0.0D0
                  PLUMEBKGRP(ISRC)  = 0.0D0
                END IF
                PLUMEBKGO3(ISRC)= O3CONC_BG        
              END IF
                
              INSTCONCNO2 = INSTCONCNO2 + CONCSOURCENO2
              INSTCONCNO  = INSTCONCNO  + CONCSOURCENO
              NOXPERSOURCE(ISRC) = CONCSOURCENO2 + CONCSOURCENO
              INSTCONCNOX = INSTCONCNOX + NOXPERSOURCE(ISRC)

c             INSTCONCVOC = INSTCONCVOC + CONCSOURCEVOC 
              INSTCONCROC = INSTCONCROC + CONCSOURCEROC
c             Not calcualting INSTCONCRP
              VOCPERSOURCE(ISRC) = CONCSOURCEVOC
              ROCPERSOURCE(ISRC) = CONCSOURCEROC
              VOCALLSRCS         = VOCALLSRCS + CONCSOURCEVOC 
              ROCALLSRCS         = ROCALLSRCS + CONCSOURCEROC
                
              IF (NOXPERSOURCE(ISRC) > CONCMAXNOX) THEN
                CONCMAXNOX = NOXPERSOURCE(ISRC)
              END IF
              IF (VOCPERSOURCE(ISRC) > CONCMAXVOC) THEN
                CONCMAXVOC = VOCPERSOURCE(ISRC)
              END IF
                
            END DO
           
            IF(GRP_BACK(IGRP_USE).AND..NOT.L_NIGHTHOUR)THEN
              MINBKGNO2 = NO2CONC_BG
              MINBKGNO  = NOCONC_BG
              MINBKGROC = ROCCONC_BG
              MINBKGRP  = RPCONC_BG
            ELSE
              MINBKGNO2 = 0.0D0
              MINBKGNO  = 0.0D0
              MINBKGROC = 0.0D0
              MINBKGRP  = 0.0D0
            END IF
            MINBKGO3  = O3CONC_BG
            
            !Calculate the minimum background based on O3
            DO ISRC = 1, NUMSRC
              IF (NOXPERSOURCE(ISRC)<=0.5*CONCMAXNOX) CYCLE
              IF (MINBKGO3 > PLUMEBKGO3(ISRC)) THEN          
                MINBKGNO2 = PLUMEBKGNO2(ISRC)
                MINBKGNO  = PLUMEBKGNO(ISRC)
                MINBKGO3  = PLUMEBKGO3(ISRC)
                MINBKGROC = PLUMEBKGROC(ISRC)
                MINBKGRP  = PLUMEBKGRP(ISRC)
              END IF
            END DO
            
            !Put this set of background values into equilibrium
            IF(GRP_BACK(IGRP_USE).AND..NOT.L_NIGHTHOUR)THEN
              IF (MINBKGROC > 0.0D0) THEN    
                 CALL QuadraticEquil_GRS7(MINBKGNO2,MINBKGNO,
     &                MINBKGO3,MINBKGROC,MINBKGRP)
              ELSE              
                 CALL QuadraticEquil_GRSM(MINBKGNO2,MINBKGNO,
     &                MINBKGO3)
              END IF
            ENDIF
            
            INSTCONCNO2 = INSTCONCNO2 + MINBKGNO2
            INSTCONCNO  = INSTCONCNO  + MINBKGNO
            INSTCONCROC = INSTCONCROC + MINBKGROC
            
C ---       Calculate weighted time
            IF(SUM(CHI(IREC,:,1),MASK=IGROUP(:,IGRP_USE)==1)/=0.0D0)THEN
              TTRAVCHM(IREC)=
     &        SUM(CHI_TTRAVCHM(IREC,:),MASK=IGROUP(:,IGRP_USE)==1)/
     &        SUM(CHI(IREC,:,1),MASK=IGROUP(:,IGRP_USE)==1)
            ELSE
              TTRAVCHM(IREC)=0.0D0  
            END IF
C ---       Set time to be a maximum of one hour
            TTRAVCHM(IREC)=MIN(TTRAVCHM(IREC),3600.0D0)
            
            CONCTEMP(nNO)  = INSTCONCNO
            CONCTEMP(nNO2) = INSTCONCNO2
            CONCTEMP(nO3)  = MINBKGO3  
            CONCTEMP(nROC)  = INSTCONCROC
            CONCTEMP(nRP)   = MINBKGRP
            CONCTEMP(nSGN)  = 0.
            CONCTEMP(nSNGN) = 0.        

            IF(INLNDEBUG.AND.IREC==3)write(*,'(A,(A,I1),15(A,D9.3))')
     +          "pre-DoGRS7Chem",
     &          " | IREC=",IREC,
     &          " | Travel_Time_sec=",TTRAVCHM(IREC),
     &          " | total_plume_NO2=",CONCTEMP(nNO2),
     &          " | total_plume_NO=",CONCTEMP(nNO),
     &          " | total_plume_NOx=",CONCTEMP(nNO)+CONCTEMP(nNO2),
     &          " | total_plume_ROC=",CONCTEMP(nROC),
     &          " | total_plume_RP=",CONCTEMP(nRP),
     &          " | total_plume_SGN=",CONCTEMP(nSGN),
     &          " | total_plume_SNGN=",CONCTEMP(nSNGN),
     &          " | bgnd_NO2=",MINBKGNO2,
     &          " | bgnd_NO=",MINBKGNO,
     &          " | bgnd_O3=",CONCTEMP(nO3),
     &          " | bgnd_ROC=",MINBKGROC,
     &          " | bgnd_RP=",MINBKGRP,
     &          " | orig_plume_NOX=",INSTCONCNOX,
     &          " | orig_plume_VOC=",VOCALLSRCS
     
C ---       Do GRS7 chemistry without dilution and entrainment
            CALL DoGRS7Chem 
        
            IF(INLNDEBUG.AND.IREC==3)write(*,'(A,(A,I1),15(A,D9.3))')
     +          "pst-DoGRS7Chem",
     &          " | IREC=",IREC,
     &          " | Travel_Time_sec=",TTRAVCHM(IREC),
     &          " | total_plume_NO2=",CONCTEMP(nNO2),
     &          " | total_plume_NO=",CONCTEMP(nNO),
     &          " | total_plume_NOx=",CONCTEMP(nNO)+CONCTEMP(nNO2),
     &          " | total_plume_ROC=",CONCTEMP(nROC),
     &          " | total_plume_RP=",CONCTEMP(nRP),
     &          " | total_plume_SGN=",CONCTEMP(nSGN),
     &          " | total_plume_SNGN=",CONCTEMP(nSNGN),
     &          " | bgnd_NO2=",MINBKGNO2,
     &          " | bgnd_NO=",MINBKGNO,
     &          " | bgnd_O3=",CONCTEMP(nO3),
     &          " | bgnd_ROC=",MINBKGROC,
     &          " | bgnd_RP=",MINBKGRP,
     &          " | orig_plume_NOX=",INSTCONCNOX,
     &          " | orig_plume_VOC=",VOCALLSRCS

            INSTCONCNO  = CONCTEMP(nNO)  - MINBKGNO
            INSTCONCNO2 = CONCTEMP(nNO2) - MINBKGNO2
            INSTCONCROC = CONCTEMP(nROC) - MINBKGROC
            
            ENSCONCNO   = 0.0
            ENSCONCNO2  = 0.0
            ENSCONCROC  = 0.0
            
            DO ISRC = 1, NUMSRC
              SRCNO2(isrc) = 0.0D0
              SRCNO(isrc) = 0.0D0
              IF (INSTCONCNOX==0.OR.NOXPERSOURCE(ISRC)==0) THEN
                CONCSOURCENO2 = 0.0
                CONCSOURCENO  = 0.0
                CONCSOURCEROC = 0.0
              ELSE
                CONCSOURCENO2=INSTCONCNO2*NOXPERSOURCE(ISRC)/INSTCONCNOX
                CONCSOURCENO =INSTCONCNO *NOXPERSOURCE(ISRC)/INSTCONCNOX
                CONCSOURCEROC=INSTCONCROC*ROCPERSOURCE(ISRC)/ROCALLSRCS
              
                CALL SETSRC
                WDSIN = AWDSIN(ISRC)
                WDCOS = AWDCOS(ISRC)
                IF (EVONLY) THEN
                  CALL XYDIST(IEVENT)
                ELSE
                  CALL XYDIST(IREC)
                END IF
                IF(X > fDownstrmStart)THEN
                  CONCSOURCENO2=CONCSOURCENO2*
     +                         SIGYSIGZ_I_X(ISRC)/(SIGYSIGZ_E_X(ISRC))
                  CONCSOURCENO =CONCSOURCENO*
     +                         SIGYSIGZ_I_X(ISRC)/(SIGYSIGZ_E_X(ISRC))
                  CONCSOURCEROC=CONCSOURCEROC*
     +                         SIGYSIGZ_I_X(ISRC)/(SIGYSIGZ_E_X(ISRC))
                END IF
                SRCNO2(isrc) = CONCSOURCENO2 
                SRCNO(isrc)  = CONCSOURCENO
              END IF
              
              ! Calculate ensemable total NO2 and NO, and ROC
              ENSCONCNO2 = ENSCONCNO2 + CONCSOURCENO2
              ENSCONCNO  = ENSCONCNO  + CONCSOURCENO
              ENSCONCROC = ENSCONCROC + CONCSOURCEROC
            END DO
            
            ! HTran: Why background concentration is added here?
            ! Answer: beacause ENSCONCNO2 is not the objective NO2 concentration output in GRSM
            ! GRSm objective is to estimated new NO2Frac
            ! To-do list: Think how to adopt this for GRS7
            IF(GRP_BACK(IGRP_USE).AND..NOT.L_NIGHTHOUR)THEN
              ENSCONCNO2 = ENSCONCNO2 + NO2CONC_BG
              ENSCONCNO  = ENSCONCNO  + NOCONC_BG
              ENSCONCROC = ENSCONCROC + ROCCONC_BG      
            END IF
            
C ---       Calculate secondary NO2 fraction after chemistry
            IF((ENSCONCNO2+ENSCONCNO)/=0.0D0)THEN
              NO2Frac = ENSCONCNO2 / (ENSCONCNO2 + ENSCONCNO)  
            ELSE
              NO2Frac = 0.0D0
            END IF
C ---       Calculate post chemistry concentrations from secondary NO2 fraction, partitioning background and plume
            IF(L_NIGHTHOUR)THEN
              !Night-time
              BGCONC_OUT = NO2CONC_BG/NO2_PPB  !No NO at night-time
              BGCONC_PPB = NO2CONC_BG
            ELSE
              !Day-time
              !Background NO2 (other background concs not o/p)
              BGCONC_OUT = NO2Frac * NOXCONC_BG/NO2_PPB 
              BGCONC_PPB = NO2Frac * NOXCONC_BG
            END IF
          ELSE  !Missing O3 or NOx doing full conversion
            NO2Frac=1.0D0
            !Warning has already been issued in O3EXT
            BGCONC_OUT = BGCONC_IN
            BGCONC_PPB = BGCONC_IN*NO2_PPB
          END IF

          
          SOURCE_LOOP: DO ISRC = 1, NUMSRC          

            IF( IGROUP(ISRC,IGRP_USE) .NE. 1 )THEN
C ---         This source is not included in the current SRCGROUP
              CYCLE SOURCE_LOOP            
            END IF
 
            !Plume
            DO ITYP = 1,NUMTYP
              HRVAL(ITYP) = NO2Frac * CHI(IREC,ISRC,ITYP)  
            END DO 
                   
C ---       Call SUMVAL to add current source's contribution to AVEVAL
            IF( .NOT.EVONLY )THEN
C ---         Processing for NON-EVENT application 
              CALL SUMVAL            
              
C ---         Write debug file
              IF (GRSMDEBUG) THEN
                IF (IGROUP(ISRC,IGRP) .EQ. 1) THEN
                  WRITE(GRSMDBG,9989,ERR=999) KURDAT, IREC,
     &                 GRPID(IGRP),  ISRC, SRCID(ISRC), NIGHTSTR, 
     &                 O3CONC_BG,NOXCONC_BG,NO2CONC_BG,TTRAVCHM(IREC),
     &                 CHI_TTRAVCHM(IREC,ISRC), ANO2_RATIO(isrc),
     &                 CHI(IREC,ISRC,1), NO2Frac, BGCONC_PPB, HRVAL(1),
     &                 AVEVAL(IREC,IGRP,1,1),SIGYSIGZ_I_X(ISRC),
     &                 SIGYSIGZ_E_X(ISRC), BLDFAC(IREC,ISRC),fMultPlmFac
                END IF
              END IF
9989          FORMAT(1x,7x,i8.8,2x,i6,2x,a8,2x,i6,2x,a12,2x,a3,8x,2x,7x,
     &               e12.5,2x,8x,e12.5,2x,8x,e12.5,2x,3x,e12.5,2x,20x,
     &               e12.5,2x,4x,e12.5,2x,13x,e12.5,2x,2x,e12.5,2x,7x,
     &               e12.5,2x,5x,e12.5,2x,6x,e12.5,2x,7x,e12.5,2x,10x,
     &               e12.5,2x,4x,e12.5,2x,6x,e12.5)
            ELSE
C ---         Processing for EVENT application
              CALL EV_SUMVAL             
              
C ---         Write debug file
              IF (GRSMDEBUG) THEN
                IF (IGROUP(ISRC,IDXEV(IEVENT)) .EQ. 1) THEN
                  WRITE(GRSMDBG,9990,ERR=999) KURDAT, IEVENT, 
     &                 EVNAME(IEVENT), EVAPER(IEVENT), 
     &                 GRPID(IDXEV(IEVENT)), ISRC, SRCID(ISRC), 
     &                 NIGHTSTR, O3CONC_BG,NOXCONC_BG, NO2CONC_BG,
     &                 TTRAVCHM(IREC), CHI_TTRAVCHM(IREC,ISRC), 
     &                 ANO2_RATIO(isrc),CHI(IREC,ISRC,1), NO2Frac, 
     &                 BGCONC_PPB, HRVAL(1), GRPAVE(IDXEV(IEVENT)),
     &                 SIGYSIGZ_I_X(ISRC),SIGYSIGZ_E_X(ISRC),
     &                 BLDFAC(IREC,ISRC), fMultPlmFac
                END IF
              END IF
9990          FORMAT(1x,7x,i8.8,2x,i6,2x,a10,2x,1x,i3,2x,a8,2x,i6,2x,
     &               a12,2x,a3,8x,2x,7x,e12.5,2x,8x,e12.5,2x,8x,e12.5,
     &               2x,3x,e12.5,2x,20x,e12.5,2x,4x,e12.5,2x,13x,e12.5,
     &               2x,2x,e12.5,2x,7x,e12.5,2x,5x,e12.5,2x,9x,e12.5,2x,
     &               7x,e12.5,2x,10x,e12.5,2x,4x,e12.5,2x,6x,e12.5)
            END IF
          
            IF (EVAL(ISRC)) THEN
C ---         Check ARC centerline values for EVALFILE output - this is same for EVENT/NON-EVENT
              CALL EVALCK
            END IF              

C ---       ReInitialize __VAL arrays (1:NUMTYP)
            HRVAL  = 0.0D0              
          END DO SOURCE_LOOP          
            
          IF (GRP_BACK(IGRP_USE)) THEN
            IF( .NOT.EVONLY )THEN  
C ---         Call SUMBACK_NO2 to update BACKGROUND contributions
              BGCONC = BGCONC_OUT
              CALL SUMBACK_NO2
            ELSE
C ---         Call EV_SUMBACK to update BACKGROUND contributions
              EV_BGCONC(IHOUR) = BGCONC_OUT
              CALL EV_SUMBACK
            ENDIF
          ENDIF
        END DO  !Loop over points
      END DO  !Loop over groups

      GO TO 3456
999   CALL ERRHDL(PATH,MODNAM,'E','520','GRSMDEBUG')
3456  CONTINUE

      RETURN    
      END SUBROUTINE GRS7_CALC