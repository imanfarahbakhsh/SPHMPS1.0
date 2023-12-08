MODULE MODULE1
IMPLICIT NONE
!GEOMETRICAL-------------------------------------------------------------------
INTEGER,PARAMETER            :: IM=30
INTEGER,PARAMETER            :: JM=30
INTEGER,PARAMETER            :: BBB=1           !Lower Boundary Lower Limit
INTEGER,PARAMETER            :: BBU=3           !Lower Boundary Upper Limit
INTEGER,PARAMETER            :: UBU=JM          !Upper Boundary Upper Limit
INTEGER,PARAMETER            :: UBB=JM-BBU+BBB  !Upper Boundary Lower Limit
INTEGER,PARAMETER            :: LBL=1           !Left Boundary Left Limit
INTEGER,PARAMETER            :: LBR=3           !Left Boundary Right Limit
INTEGER,PARAMETER            :: RBR=IM          !Right Boundary Right Limit
INTEGER,PARAMETER            :: RBL=IM-LBR+LBL  !Right Boundary Left Limit
INTEGER,PARAMETER            :: KPLATE=IM*JM
REAL(8),PARAMETER            :: PI=3.1415926535897932D0
REAL(8),PARAMETER            :: LX_PLATE=1.D0!0.01D0
REAL(8),PARAMETER            :: LY_PLATE=1.D0!0.01D0
REAL(8),PARAMETER            :: DX=LX_PLATE/(IM-1)
REAL(8),PARAMETER            :: DY=LY_PLATE/(JM-1)
REAL(8),PARAMETER            :: LX_PIN=1.D0!0.01D0
REAL(8),PARAMETER            :: LY_PIN=1.D0!0.01D0
REAL(8),PARAMETER            :: LX=LX_PLATE
REAL(8),PARAMETER            :: LY=LY_PLATE
REAL(8),PARAMETER            :: ORP=0.25D0!0.003175  !Outer Radius of Pin
REAL(8),PARAMETER            :: XCI=-0.05D0    !Inital x-position of pin center !Can be changed depend on the resolution
REAL(8),PARAMETER            :: YCI=0.D0                         !Initial y-position of pin center
REAL(8),PARAMETER            :: SIDE=DSQRT(PI)*ORP
REAL(8),PARAMETER            :: LSIDE=XCI-0.5D0*SIDE
REAL(8),PARAMETER            :: RSIDE=XCI+0.5D0*SIDE
REAL(8),PARAMETER            :: USIDE=YCI+0.5D0*SIDE
REAL(8),PARAMETER            :: BSIDE=YCI-0.5D0*SIDE
REAL(8),PARAMETER            :: SIDE_TA=DSQRT(3.D0)*ORP!2.D0*DSQRT(PI)*ORP/(3.D0**0.25D0)!
REAL(8),PARAMETER            :: X1=XCI-0.5D0*SIDE_TA
REAL(8),PARAMETER            :: X2=XCI+0.5D0*SIDE_TA
REAL(8),PARAMETER            :: X3=XCI
REAL(8),PARAMETER            :: Y1=YCI-0.5D0*SIDE_TA/DSQRT(3.D0)
REAL(8),PARAMETER            :: Y2=YCI-0.5D0*SIDE_TA/DSQRT(3.D0)
REAL(8),PARAMETER            :: Y3=YCI+SIDE_TA/DSQRT(3.D0)
REAL(8),PARAMETER            :: M1=(Y3-Y1)/(X3-X1)
REAL(8),PARAMETER            :: M2=(Y3-Y2)/(X3-X2)
!COMPUTATIONAL-----------------------------------------------------------------
INTEGER,PARAMETER            :: TM=100
REAL(8),PARAMETER            :: DT=1.D-6
!KINEMATIC---------------------------------------------------------------------
REAL(8),PARAMETER            :: UT=0.366D0!0.00465D0
REAL(8),PARAMETER            :: VT=0.D0
REAL(8),PARAMETER            :: OMEGA=2.D0*PI*1200/60.D0
!OUTPUT------------------------------------------------------------------------
INTEGER,PARAMETER            :: NPRINT=10
LOGICAL,PARAMETER            :: ICFLAG=.FALSE.
INTEGER,PARAMETER            :: KSTART=0
!INPUT-------------------------------------------------------------------------

!PHYSICAL----------------------------------------------------------------------
REAL(8),PARAMETER            :: H_CONV=20.D0
REAL(8),PARAMETER            :: T_INFTY=298.15D0
REAL(8),PARAMETER            :: SOUND_VELOCITY_ZERO_TOOL=5800.D0      !-------->>
REAL(8),PARAMETER            :: SOUND_VELOCITY_ZERO_PLATE=5800.D0     !-------->>
REAL(8),PARAMETER            :: RHO_ZERO_TOOL=7800.D0
REAL(8),PARAMETER            :: RHO_ZERO_PLATE=1777.D0
REAL(8),PARAMETER            :: K_COND_ZERO_TOOL=24.3D0
REAL(8),PARAMETER            :: K_COND_ZERO_PLATE=400.D0
REAL(8),PARAMETER            :: CP_ZERO_TOOL=460.D0
REAL(8),PARAMETER            :: CP_ZERO_PLATE=1164.D0
REAL(8),PARAMETER            :: P_ZERO_TOOL=RHO_ZERO_TOOL*SOUND_VELOCITY_ZERO_TOOL**2.D0
REAL(8),PARAMETER            :: P_ZERO_PLATE=RHO_ZERO_PLATE*SOUND_VELOCITY_ZERO_PLATE**2.D0
REAL(8),PARAMETER            :: Q_ACT_EN=129.D3
REAL(8),PARAMETER            :: R_GAS_CONST=8.314D0
REAL(8),PARAMETER            :: SIGMA_R=53.3D6
REAL(8),PARAMETER            :: MASS0=1.D0
REAL(8),PARAMETER            :: T_LIQUIDUS=630.D0
REAL(8),PARAMETER            :: T_SOLIDUS=575.D0
REAL(8),PARAMETER            :: MI=80.D0
REAL(8),PARAMETER            :: MS=1.D6
REAL(8),PARAMETER            :: MU_MELTING=0.15D0
REAL(8),PARAMETER            :: MU_PIN=1.D8
REAL(8),PARAMETER            :: FS_PIN=2.D0
!SPH---------------------------------------------------------------------------
REAL(8),PARAMETER            :: BETA=7.78D8
REAL(8),PARAMETER            :: M=4.36D0
REAL(8),PARAMETER            :: INVM=1.D0/M
REAL(8),PARAMETER            :: W3=1.D-5
REAL(8),PARAMETER            :: EPSILON=0.01D0
REAL(8),PARAMETER            :: H_ZERO=0.5D0*DX
!------------------------------------------------------------------------------
END MODULE MODULE1
!==============================================================================
PROGRAM MAIN
USE MODULE1
IMPLICIT NONE
INTEGER                                :: K,KB,KPIN,KTOTAL,KBTOTAL,TIME_STEP
REAL(8)                                :: XC,YC,MASS
REAL(8),DIMENSION(IM,JM)               :: XPIN,YPIN
REAL(8),DIMENSION(:),ALLOCATABLE       :: X,Y,XBOUND,YBOUND,RP,ITHETA,F_BOUND_X,F_BOUND_Y,V_MAX
REAL(8),DIMENSION(:),ALLOCATABLE       :: RHO,RHS_MASS_CONS,U,V,FS,MU,T,VIS_TERM,RHS_ENERGY_CONS
REAL(8),DIMENSION(:),ALLOCATABLE       :: EPSILON_11,EPSILON_12,EPSILON_21,EPSILON_22,RHO_ZERO,P_ZERO
REAL(8),DIMENSION(:),ALLOCATABLE       :: RHS_MOMENTUM_CONS_X,RHS_MOMENTUM_CONS_Y,CP,K_COND,PHASE
INTEGER,DIMENSION(:),ALLOCATABLE       :: BOUND_PARTICLE_COUNTER,NONBOUND_PARTICLE_COUNTER
!------------------------------------------------------------------------------
CALL PIN_PLATE_PARTICLES (XPIN,YPIN,K,KPIN,KB,KTOTAL,KBTOTAL)
ALLOCATE(X(KTOTAL),Y(KTOTAL),RP(KPIN),ITHETA(KPIN),BOUND_PARTICLE_COUNTER(KBTOTAL),V_MAX(KTOTAL))
ALLOCATE(XBOUND(KBTOTAL),YBOUND(KBTOTAL),RHO(KTOTAL),U(KTOTAL),V(KTOTAL),RHS_MASS_CONS(KTOTAL),PHASE(KTOTAL))
ALLOCATE(EPSILON_11(KTOTAL),EPSILON_12(KTOTAL),EPSILON_21(KTOTAL),EPSILON_22(KTOTAL),P_ZERO(KTOTAL))
ALLOCATE(MU(KTOTAL),FS(KTOTAL),T(KTOTAL),VIS_TERM(KTOTAL),F_BOUND_X(KTOTAL),F_BOUND_Y(KTOTAL),RHO_ZERO(KTOTAL))
ALLOCATE(RHS_MOMENTUM_CONS_X(KTOTAL),RHS_MOMENTUM_CONS_Y(KTOTAL),CP(KTOTAL),K_COND(KTOTAL),RHS_ENERGY_CONS(KTOTAL))
ALLOCATE(NONBOUND_PARTICLE_COUNTER(KTOTAL-KBTOTAL))
CALL PINPARTICLES (K,X,Y,XPIN,YPIN,XBOUND,YBOUND,KPIN,KTOTAL,KBTOTAL,BOUND_PARTICLE_COUNTER,XC,YC,RP,ITHETA)
CALL PLATEPARTICLES (K,KB,KTOTAL,KBTOTAL,X,Y,XBOUND,YBOUND,&
                     BOUND_PARTICLE_COUNTER,NONBOUND_PARTICLE_COUNTER)
CALL PARTICLE_PHASE (KPIN,KTOTAL,Y,PHASE)
!CALL WRITE_INITIAL_POSITION (KTOTAL,X,Y,PHASE)
CALL INITIAL_CONDITIONS (KPIN,KTOTAL,X,Y,U,V,T,PHASE,MASS,RHO_ZERO,K_COND,CP,&
                               P_ZERO,RHO,RHS_MASS_CONS,RHS_MOMENTUM_CONS_X,&
                               RHS_MOMENTUM_CONS_Y,RHS_ENERGY_CONS,XC,YC)
!-------------------------------------------------------------------Time Loop
DO TIME_STEP=KSTART+1,TM
PRINT*,TIME_STEP
!IF (((TIME_STEP/NPRINT)*NPRINT)==TIME_STEP)&
!CALL READ_TIME_SOLUTION (TIME_STEP,KTOTAL,Y,PHASE)
CALL RIGID_CENTER_POSITION (XC,YC)
CALL RIGID_GEOMETRY_POSITION_VELOCITY (TIME_STEP,KTOTAL,KPIN,RP,ITHETA,XC,YC,X,Y,U,V)
CALL RHS_MASS_CONSERVATION (KTOTAL,X,Y,MASS,RHO,RHO_ZERO,U,V,RHS_MASS_CONS)
!!CALL STRAIN_RATE_TENSOR (KTOTAL,MASS,RHO,RHO_ZERO,U,V,X,Y,EPSILON_11,EPSILON_12,EPSILON_21,EPSILON_22)
!!CALL VISCOSITY_1 (KTOTAL,EPSILON_11,EPSILON_12,EPSILON_21,EPSILON_22,T,MU)
!!CALL SOLID_FRACTION (KTOTAL,KPIN,T,FS)
!!CALL VISCOSITY_2 (KPIN,KTOTAL,FS,MU)
CALL VISCOSITY_3 (KPIN,KTOTAL,T,MU)
CALL VISCOSITY_TERM (KTOTAL,MU,RHO,RHO_ZERO,U,V,X,Y,VIS_TERM)
CALL MAXIMUM_VELOCITY (KBTOTAL,KTOTAL,BOUND_PARTICLE_COUNTER,NONBOUND_PARTICLE_COUNTER,U,V,V_MAX)
CALL BOUNDARY_FORCE (KBTOTAL,KTOTAL,BOUND_PARTICLE_COUNTER,NONBOUND_PARTICLE_COUNTER,MASS,RHO,RHO_ZERO,&
                           V_MAX,X,Y,F_BOUND_X,F_BOUND_Y)
CALL RHS_MOMENTUM_CONSERVATION (KTOTAL,X,Y,F_BOUND_X,F_BOUND_Y,&
                                VIS_TERM,MASS,RHO,RHO_ZERO,P_ZERO,RHS_MOMENTUM_CONS_X,&
                                RHS_MOMENTUM_CONS_Y)
CALL RHS_ENERGY_CONSERVATION (KTOTAL,MASS,RHO,RHO_ZERO,CP,T,U,V,K_COND,VIS_TERM,&
                              X,Y,RHS_ENERGY_CONS)
CALL EULER_INTEGRATION (KBTOTAL,KTOTAL,NONBOUND_PARTICLE_COUNTER,&
                              RHS_MASS_CONS,RHS_MOMENTUM_CONS_X,&
                              RHS_MOMENTUM_CONS_Y,RHS_ENERGY_CONS,RHO,U,V,T)
CALL PARTICLE_POSITION (KBTOTAL,KTOTAL,NONBOUND_PARTICLE_COUNTER,X,Y,U,V)

IF (((TIME_STEP/NPRINT)*NPRINT)==TIME_STEP) THEN
!CALL PHASE_STATISTIC (KTOTAL,TIME_STEP,PHASE,Y)
CALL WRITE_TIME_SOLUTION (TIME_STEP,KTOTAL,X,Y,RHO,U,V,T,PHASE)!CALL WRITE_TIME_SOLUTION (TIME_STEP,KTOTAL,KBTOTAL,NONBOUND_PARTICLE_COUNTER,X,Y,RHO,U,V,T)
END IF
END DO
!------------------------------------------------------------------------------
END PROGRAM MAIN
!==============================================================================
SUBROUTINE PIN_PLATE_PARTICLES (XPIN,YPIN,K,KPIN,KB,KTOTAL,KBTOTAL)
USE MODULE1
IMPLICIT NONE
INTEGER                            :: I,J,K,KPIN,KB,KTOTAL,KBTOTAL
REAL(8),DIMENSION(IM,JM)           :: XPIN,YPIN
!------------------------------------------------------------------------------
!YPIN(:,1)=-0.5D0*LY_PIN
!DO I=1,IM
!DO J=1,JM-1
!YPIN(I,J+1)=YPIN(I,J)+DY
!END DO
!END DO
!!------------------------------------------------------------------------------
!XPIN(1,:)=-0.5D0*LX_PLATE-LX_PIN
!DO J=1,JM
!DO I=1,IM-1
!XPIN(I+1,J)=XPIN(I,J)+DX
!END DO
!END DO
!!------------------------------------------------------------------------------
!K=0
!DO J=1,JM
!DO I=1,IM
!    IF ((XPIN(I,J)-XCI)**2.D0+(YPIN(I,J)-YCI)**2.D0<=ORP**2.D0) THEN
!        K=K+1
!    END IF
!END DO
!END DO
!KPIN=K
!KB=K
!KTOTAL=KPIN+KPLATE
!KBTOTAL=KPIN+2*BBU*IM
!------------------------------------------------------------------------------
YPIN(:,1)=-0.5D0*LY
DO I=1,IM
DO J=1,JM-1
YPIN(I,J+1)=YPIN(I,J)+DY
END DO
END DO
!------------------------------------------------------------------------------
XPIN(1,:)=-0.5D0*LX
DO J=1,JM
DO I=1,IM-1
XPIN(I+1,J)=XPIN(I,J)+DX
END DO
END DO
!------------------------------------------------------------------------------
K=0
DO J=1,JM
DO I=1,IM
    IF ((YPIN(I,J)>=Y1).AND.(YPIN(I,J)<=(M1*(XPIN(I,J)-X1)+Y1)).AND.(YPIN(I,J)<=(M2*(XPIN(I,J)-X2)+Y2))) THEN
        K=K+1
    END IF
END DO
END DO
KPIN=K
KB=K
KTOTAL=KPLATE
KBTOTAL=KPIN+2*BBU*IM+2*LBR*(JM-2*BBU)
!------------------------------------------------------------------------------
END SUBROUTINE PIN_PLATE_PARTICLES
!==============================================================================
SUBROUTINE PINPARTICLES (K,X,Y,XPIN,YPIN,XBOUND,YBOUND,KPIN,KTOTAL,KBTOTAL,BOUND_PARTICLE_COUNTER,XC,YC,RP,ITHETA)
USE MODULE1
IMPLICIT NONE
INTEGER                            :: I,J,K,KTOTAL,KBTOTAL,KPIN
REAL(8)                            :: XC,YC
REAL(8),DIMENSION(IM,JM)           :: XPIN,YPIN
REAL(8),DIMENSION(KTOTAL)          :: X,Y
REAL(8),DIMENSION(KBTOTAL)         :: XBOUND,YBOUND
REAL(8),DIMENSION(KPIN)            :: RP,ITHETA
INTEGER,DIMENSION(KBTOTAL)         :: BOUND_PARTICLE_COUNTER
!------------------------------------------------------------------------------
!K=0
!DO J=1,JM
!DO I=1,IM
!    IF ((XPIN(I,J)-XCI)**2.D0+(YPIN(I,J)-YCI)**2.D0<=ORP**2.D0) THEN
!        K=K+1
!        X(K)=XPIN(I,J)
!        Y(K)=YPIN(I,J)
!        RP(K)=DSQRT((X(K)-XCI)**2.D0+(Y(K)-YCI)**2.D0)
!        ITHETA(K)=DATAN2((Y(K)-YCI),(X(K)-XCI))+2.D0*PI
!        XBOUND(K)=X(K)
!        YBOUND(K)=Y(K)
!        BOUND_PARTICLE_COUNTER(K)=K
!    END IF
!END DO
!END DO
!XC=XCI
!YC=YCI
!------------------------------------------------------------------------------
!K=0
!DO J=1,JM
!DO I=1,IM
!    IF (((XPIN(I,J)<=RSIDE).AND.(XPIN(I,J)>=LSIDE)).AND.((YPIN(I,J)<=USIDE).AND.(YPIN(I,J)>=BSIDE))) THEN
!        K=K+1
!        X(K)=XPIN(I,J)
!        Y(K)=YPIN(I,J)
!        RP(K)=DSQRT((X(K)-XCI)**2.D0+(Y(K)-YCI)**2.D0)
!        ITHETA(K)=DATAN2((Y(K)-YCI),(X(K)-XCI))+2.D0*PI
!        XBOUND(K)=X(K)
!        YBOUND(K)=Y(K)
!        BOUND_PARTICLE_COUNTER(K)=K
!    END IF
!END DO
!END DO
!XC=XCI
!YC=YCI
!------------------------------------------------------------------------------
K=0
DO J=1,JM
DO I=1,IM
    IF ((YPIN(I,J)>=Y1).AND.(YPIN(I,J)<=(M1*(XPIN(I,J)-X1)+Y1)).AND.(YPIN(I,J)<=(M2*(XPIN(I,J)-X2)+Y2))) THEN
        K=K+1
        X(K)=XPIN(I,J)
        Y(K)=YPIN(I,J)
        RP(K)=DSQRT((X(K)-XCI)**2.D0+(Y(K)-YCI)**2.D0)
        ITHETA(K)=DATAN2((Y(K)-YCI),(X(K)-XCI))+2.D0*PI
        XBOUND(K)=X(K)
        YBOUND(K)=Y(K)
        BOUND_PARTICLE_COUNTER(K)=K
    END IF
END DO
END DO
XC=XCI
YC=YCI
!------------------------------------------------------------------------------
END SUBROUTINE PINPARTICLES
!==============================================================================
SUBROUTINE PLATEPARTICLES (K,KB,KTOTAL,KBTOTAL,X,Y,XBOUND,YBOUND,&
                           BOUND_PARTICLE_COUNTER,NONBOUND_PARTICLE_COUNTER)
USE MODULE1
IMPLICIT NONE
INTEGER                                 :: I,J,K,KB,KNB,KTOTAL,KBTOTAL
REAL(8),DIMENSION(IM,JM)                :: XPL,YPL
REAL(8),DIMENSION(KTOTAL)               :: X,Y
REAL(8),DIMENSION(KBTOTAL)              :: XBOUND,YBOUND
INTEGER,DIMENSION(KBTOTAL)              :: BOUND_PARTICLE_COUNTER
INTEGER,DIMENSION(KTOTAL-KBTOTAL)       :: NONBOUND_PARTICLE_COUNTER
!------------------------------------------------------------------------------
!XPL(1,:)=-0.5D0*LX_PLATE
!DO J=1,JM
!    DO I=1,IM-1
!      XPL(I+1,J)=XPL(I,J)+DX
!    END DO
!END DO
!YPL(:,1)=-0.5D0*LY_PLATE
!DO I=1,IM
!    DO J=1,JM-1
!      YPL(I,J+1)=YPL(I,J)+DY
!    END DO
!END DO
!!------------------------------------------------------------------------------
!KNB=0       !K Non-Boundary !Non-Boundary Particles Counter
!!------------------------------------------------------------------------------
!DO J=1,JM
!DO I=1,IM
!        K=K+1
!        X(K)=XPL(I,J)
!        Y(K)=YPL(I,J)
!        IF (((J>=BBB).AND.(J<=BBU)).OR.((J>=UBB).AND.(J<=UBU))) THEN
!        KB=KB+1
!        XBOUND(KB)=XPL(I,J)
!        YBOUND(KB)=YPL(I,J)
!        BOUND_PARTICLE_COUNTER(KB)=K
!        ELSE
!        KNB=KNB+1
!        NONBOUND_PARTICLE_COUNTER(KNB)=K
!        END IF
!    END DO
!END DO
!------------------------------------------------------------------------------
XPL(1,:)=-0.5D0*LX
DO J=1,JM
    DO I=1,IM-1
      XPL(I+1,J)=XPL(I,J)+DX
    END DO
END DO
YPL(:,1)=-0.5D0*LY
DO I=1,IM
    DO J=1,JM-1
      YPL(I,J+1)=YPL(I,J)+DY
    END DO
END DO
!------------------------------------------------------------------------------
KNB=0       !K Non-Boundary !Non-Boundary Particles Counter
!------------------------------------------------------------------------------
DO J=1,JM
DO I=1,IM
        IF ((YPL(I,J)>=Y1).AND.(YPL(I,J)<=(M1*(XPL(I,J)-X1)+Y1)).AND.(YPL(I,J)<=(M2*(XPL(I,J)-X2)+Y2))) THEN
        CYCLE
        ELSE
        K=K+1
        X(K)=XPL(I,J)
        Y(K)=YPL(I,J)
        IF (((J>=BBB).AND.(J<=BBU)).OR.((J>=UBB).AND.(J<=UBU)).OR.((I>=LBL).AND.(I<=LBR)).OR.((I>=RBL).AND.(I<=RBR))) THEN
        KB=KB+1
        XBOUND(KB)=XPL(I,J)
        YBOUND(KB)=YPL(I,J)
        BOUND_PARTICLE_COUNTER(KB)=K
        ELSE
        KNB=KNB+1
        NONBOUND_PARTICLE_COUNTER(KNB)=K
        END IF
        END IF
END DO
END DO

!print*,NONBOUND_PARTICLE_COUNTER
!------------------------------------------------------------------------------
END SUBROUTINE PLATEPARTICLES
!==============================================================================
SUBROUTINE RIGID_CENTER_POSITION (XC,YC)
USE MODULE1
IMPLICIT NONE
REAL(8)            :: XC,YC
!------------------------------------------------------------------------------
XC=XC+UT*DT
YC=YC+VT*DT
!------------------------------------------------------------------------------
END SUBROUTINE RIGID_CENTER_POSITION
!==============================================================================
!SUBROUTINE RIGID_GEOMETRY_POSITION (TIME_STEP,KTOTAL,KPIN,RP,ITHETA,XC,YC,X,Y)
!USE MODULE1
!IMPLICIT NONE
!INTEGER                            :: K,TIME_STEP,KTOTAL,KPIN
!REAL(8)                            :: XC,YC
!REAL(8),DIMENSION(KTOTAL)          :: X,Y
!REAL(8),DIMENSION(KPIN)            :: RP,ITHETA
!------------------------------------------------------------------------------
!DO K=1,KPIN
!X(K)=RP(K)*(DCOS(ITHETA(K))*DCOS(OMEGA*TIME_STEP*DT)-DSIN(ITHETA(K))*DSIN(OMEGA*TIME_STEP*DT))+XC
!Y(K)=RP(K)*(DSIN(ITHETA(K))*DCOS(OMEGA*TIME_STEP*DT)+DCOS(ITHETA(K))*DSIN(OMEGA*TIME_STEP*DT))+YC
!END DO
!------------------------------------------------------------------------------
!END SUBROUTINE RIGID_GEOMETRY_POSITION
!==============================================================================
SUBROUTINE RIGID_GEOMETRY_POSITION_VELOCITY (TIME_STEP,KTOTAL,KPIN,RP,ITHETA,XC,YC,X,Y,U,V)
USE MODULE1
IMPLICIT NONE
INTEGER                            :: K,TIME_STEP,KTOTAL,KPIN
REAL(8)                            :: XC,YC
REAL(8),DIMENSION(KTOTAL)          :: X,Y,U,V
REAL(8),DIMENSION(KPIN)            :: RP,ITHETA
!------------------------------------------------------------------------------
DO K=1,KPIN
X(K)=RP(K)*DCOS(ITHETA(K)+OMEGA*TIME_STEP*DT)+XC
Y(K)=RP(K)*DSIN(ITHETA(K)+OMEGA*TIME_STEP*DT)+YC
U(K)=-RP(K)*OMEGA*DSIN(ITHETA(K)+OMEGA*TIME_STEP*DT)+UT
V(K)=RP(K)*OMEGA*DCOS(ITHETA(K)+OMEGA*TIME_STEP*DT)+VT
END DO
!------------------------------------------------------------------------------
END SUBROUTINE RIGID_GEOMETRY_POSITION_VELOCITY
!==============================================================================
SUBROUTINE WRITE_INITIAL_POSITION (KTOTAL,X,Y,PHASE)
USE MODULE1
IMPLICIT NONE
INTEGER                          :: K,KTOTAL
REAL(8),DIMENSION(KTOTAL)        :: X,Y,PHASE
!------------------------------------------------------------------------------
OPEN(1,FILE="INITIAL_CONFIG.TXT")
WRITE(1,*)'VARIABLES= "X", "Y", "PHASE"'
WRITE(1,*)'ZONE, K=',KTOTAL,',F=POINT'
DO K=1,KTOTAL
WRITE(1,*) X(K),Y(K),PHASE(K)
END DO
!------------------------------------------------------------------------------
OPEN(2,FILE="INITIAL_CONFIG.PLT")
WRITE(2,*)'VARIABLES= "X", "Y", "PHASE"'
WRITE(2,*)'ZONE, K=',KTOTAL,',F=POINT'
DO K=1,KTOTAL
WRITE(2,*) X(K),Y(K),PHASE(K)
END DO
CLOSE(1)
CLOSE(2)
!------------------------------------------------------------------------------
END SUBROUTINE WRITE_INITIAL_POSITION
!==============================================================================
SUBROUTINE WRITE_TIME_SOLUTION (TIME_STEP,KTOTAL,X,Y,RHO,U,V,T,PHASE)
USE MODULE1
IMPLICIT NONE
INTEGER                             :: K,TIME_STEP,COUNTER,KTOTAL
REAL(8),DIMENSION(KTOTAL)           :: X,Y,U,V,RHO,T,PHASE
CHARACTER*7                         :: EXT
CHARACTER*12                         :: FN1
CHARACTER*27                        :: FNAME
!------------------------------------------------------------------------------
FN1='/results/FLD'
COUNTER=0
!------------------------------------------------------------------------------
COUNTER=COUNTER+1
WRITE(EXT,'(I7)') TIME_STEP
FNAME=FN1//EXT//'.DAT'
!------------------------------------------------------------------------------
OPEN(COUNTER,FILE=FNAME,POSITION='REWIND')
WRITE(COUNTER,*)'VARIABLES= "X", "Y", "RHO", "U", "V", "T", "PHASE"'
WRITE(COUNTER,*)'ZONE, K=',KTOTAL,',F=POINT'
DO K=1,KTOTAL
WRITE(COUNTER,*) X(K),Y(K),RHO(K),U(K),V(K),T(K),PHASE(K)
END DO
CLOSE(COUNTER)
!------------------------------------------------------------------------------
WRITE(*,*) '========================'
WRITE(*,*) 'PRINTING ON ',FNAME
WRITE(*,*) '========================'
!------------------------------------------------------------------------------
END SUBROUTINE WRITE_TIME_SOLUTION
!==============================================================================
!SUBROUTINE WRITE_TIME_SOLUTION (TIME_STEP,KTOTAL,KBTOTAL,NONBOUND_PARTICLE_COUNTER,X,Y,RHO,U,V,T)
!USE MODULE1
!IMPLICIT NONE
!INTEGER                             :: K,TIME_STEP,COUNTER,KTOTAL,KBTOTAL
!REAL(8),DIMENSION(KTOTAL)           :: X,Y,U,V,RHO,T
!INTEGER,DIMENSION(KTOTAL-KBTOTAL)   :: NONBOUND_PARTICLE_COUNTER
!CHARACTER*7                         :: EXT
!CHARACTER*3                         :: FN1
!CHARACTER*18                        :: FNAME
!!------------------------------------------------------------------------------
!FN1='FLD'
!COUNTER=0
!!------------------------------------------------------------------------------
!COUNTER=COUNTER+1
!WRITE(EXT,'(I7)') TIME_STEP
!FNAME=FN1//EXT//'.DAT'
!!------------------------------------------------------------------------------
!OPEN(COUNTER,FILE=FNAME,POSITION='REWIND')
!WRITE(COUNTER,*)'VARIABLES= "X", "Y", "RHO", "U", "V", "T"'
!WRITE(COUNTER,*)'ZONE, K=',KTOTAL-KBTOTAL,',F=POINT'
!DO K=1,KTOTAL-KBTOTAL
!WRITE(COUNTER,*) X(NONBOUND_PARTICLE_COUNTER(K)),Y(NONBOUND_PARTICLE_COUNTER(K)),RHO(NONBOUND_PARTICLE_COUNTER(K)),&
!                 U(NONBOUND_PARTICLE_COUNTER(K)),V(NONBOUND_PARTICLE_COUNTER(K)),T(NONBOUND_PARTICLE_COUNTER(K))
!END DO
!CLOSE(COUNTER)
!!------------------------------------------------------------------------------
!WRITE(*,*) '========================'
!WRITE(*,*) 'PRINTING ON ',FNAME
!WRITE(*,*) '========================'
!!------------------------------------------------------------------------------
!END SUBROUTINE WRITE_TIME_SOLUTION
!==============================================================================
FUNCTION WEIGHTING_FUNCTION(HBARI,DISTANCEI)
USE MODULE1
IMPLICIT NONE
REAL(8)                     :: C1,Q,WEIGHTING_FUNCTION,HBARI,DISTANCEI
!------------------------------------------------------------------------------
Q=DISTANCEI/HBARI
C1=1.D0/(PI*HBARI**3.D0)
IF ((Q>=0.D0).AND.(Q<=1.D0)) THEN
WEIGHTING_FUNCTION=C1*(1.D0-1.5D0*(Q**2.D0)+0.75D0*(Q**3.D0))
ELSE IF ((Q>1.D0).AND.(Q<=2.D0)) THEN
WEIGHTING_FUNCTION=C1*(0.25D0*(2.D0-Q)**3.D0)
ELSE
WEIGHTING_FUNCTION=0.D0
END IF
!------------------------------------------------------------------------------
END FUNCTION WEIGHTING_FUNCTION
!!==============================================================================
FUNCTION WEIGHTING_FUNCTION_X(HBARI,X_DISTANCEI,DISTANCEI)
USE MODULE1
IMPLICIT NONE
REAL(8)                     :: C2,Q,WEIGHTING_FUNCTION_X,X_DISTANCEI,HBARI,DISTANCEI
!------------------------------------------------------------------------------
!Q=DISTANCEI/HBARI
!C2=1.D0/(PI*HBARI**3.D0)
!IF ((Q>=0.D0).AND.(Q<=1.D0)) THEN
!WEIGHTING_FUNCTION_X=C2*(X_DISTANCEI/DISTANCEI)*(2.25D0*(Q**2.D0)-3.D0*Q)
!ELSE IF ((Q>1.D0).AND.(Q<=2.D0)) THEN
!WEIGHTING_FUNCTION_X=C2*(X_DISTANCEI/DISTANCEI)*(-0.75D0*(2.D0-Q)**2.D0)
!ELSE
!WEIGHTING_FUNCTION_X=0.D0
!END IF
!------------------------------------------------------------------------------
Q=DISTANCEI/HBARI
C2=15.D0/(14.D0*PI*HBARI**2.D0)
IF ((Q>=0.D0).AND.(Q<=1.D0)) THEN
WEIGHTING_FUNCTION_X=(C2/(Q*HBARI**2.D0))*(-3.D0*(2.D0-Q)**2.D0+12.D0*(1.D0-Q)**2.D0)*X_DISTANCEI
ELSE IF ((Q>1.D0).AND.(Q<=2.D0)) THEN
WEIGHTING_FUNCTION_X=(C2/(Q*HBARI**2.D0))*(-3.D0*(2.D0-Q)**2.D0)*X_DISTANCEI
ELSE
WEIGHTING_FUNCTION_X=0.D0
END IF
!------------------------------------------------------------------------------
END FUNCTION WEIGHTING_FUNCTION_X
!==============================================================================
FUNCTION WEIGHTING_FUNCTION_Y(HBARI,Y_DISTANCEI,DISTANCEI)
USE MODULE1
IMPLICIT NONE
REAL(8)                     :: C2,Q,WEIGHTING_FUNCTION_Y,Y_DISTANCEI,HBARI,DISTANCEI
!------------------------------------------------------------------------------
!Q=DISTANCEI/HBARI
!C2=1.D0/(PI*HBARI**3.D0)
!IF ((Q>=0.D0).AND.(Q<=1.D0)) THEN
!WEIGHTING_FUNCTION_Y=C2*(Y_DISTANCEI/DISTANCEI)*(2.25D0*(Q**2.D0)-3.D0*Q)
!ELSE IF ((Q>1.D0).AND.(Q<=2.D0)) THEN
!WEIGHTING_FUNCTION_Y=C2*(Y_DISTANCEI/DISTANCEI)*(-0.75D0*(2.D0-Q)**2.D0)
!ELSE
!WEIGHTING_FUNCTION_Y=0.D0
!END IF
!------------------------------------------------------------------------------
Q=DISTANCEI/HBARI
C2=15.D0/(14.D0*PI*HBARI**2.D0)
IF ((Q>=0.D0).AND.(Q<=1.D0)) THEN
WEIGHTING_FUNCTION_Y=(C2/(Q*HBARI**2.D0))*(-3.D0*(2.D0-Q)**2.D0+12.D0*(1.D0-Q)**2.D0)*Y_DISTANCEI
ELSE IF ((Q>1.D0).AND.(Q<=2.D0)) THEN
WEIGHTING_FUNCTION_Y=(C2/(Q*HBARI**2.D0))*(-3.D0*(2.D0-Q)**2.D0)*Y_DISTANCEI
ELSE
WEIGHTING_FUNCTION_Y=0.D0
END IF
!------------------------------------------------------------------------------
END FUNCTION WEIGHTING_FUNCTION_Y
!==============================================================================
SUBROUTINE RHS_MASS_CONSERVATION (KTOTAL,XPOS,YPOS,MASS,RHO,RHO_ZERO,U,V,RHS_MASS_CONS)
USE MODULE1
IMPLICIT NONE
INTEGER                     :: K,KK,KTOTAL
REAL(8)                     :: WEIGHTING_FUNCTION_Y,DISTANCE
REAL(8)                     :: WEIGHTING_FUNCTION_X,H,HBAR,MASS
REAL(8),DIMENSION(KTOTAL)   :: U,V,RHO,YPOS,XPOS,RHS_MASS_CONS,RHO_ZERO
!------------------------------------------------------------------------------
RHS_MASS_CONS=0.D0
DO K=1,KTOTAL
DO KK=1,KTOTAL
    IF (K/=KK) THEN
RHS_MASS_CONS(K)=RHS_MASS_CONS(K)+(MASS/RHO(KK))*((U(K)-U(KK))*WEIGHTING_FUNCTION_X&
                 (HBAR(H(RHO(K),RHO_ZERO(K)),H(RHO(KK),RHO_ZERO(KK))),&
                XPOS(K)-XPOS(KK),DISTANCE(XPOS(K),XPOS(KK),YPOS(K),YPOS(KK)))&
                +(V(K)-V(KK))*WEIGHTING_FUNCTION_Y(HBAR(H(RHO(K),RHO_ZERO(K)),H(RHO(KK),RHO_ZERO(KK))),&
                YPOS(K)-YPOS(KK),DISTANCE(XPOS(K),XPOS(KK),YPOS(K),YPOS(KK))))
    END IF
END DO
RHS_MASS_CONS(K)=RHO(K)*RHS_MASS_CONS(K)
END DO
!------------------------------------------------------------------------------
END SUBROUTINE RHS_MASS_CONSERVATION
!==============================================================================
SUBROUTINE RHS_MOMENTUM_CONSERVATION (KTOTAL,XPOS,YPOS,F_BOUND_X,F_BOUND_Y,&
                                      VIS_TERM,MASS,RHO,RHO_ZERO,P_ZERO,RHS_MOMENTUM_CONS_X,&
                                      RHS_MOMENTUM_CONS_Y)
USE MODULE1
IMPLICIT NONE
INTEGER                     :: K,KK,KTOTAL
REAL(8)                     :: WEIGHTING_FUNCTION_Y,DISTANCE,MASS
REAL(8)                     :: WEIGHTING_FUNCTION_X,H,HBAR,P
REAL(8),DIMENSION(KTOTAL)   :: RHO,F_BOUND_X,F_BOUND_Y,XPOS,YPOS,RHO_ZERO,P_ZERO
REAL(8),DIMENSION(KTOTAL)   :: RHS_MOMENTUM_CONS_X,RHS_MOMENTUM_CONS_Y,VIS_TERM
!------------------------------------------------------------------------------
!RHS_MOMENTUM_CONS_X=0.D0
!RHS_MOMENTUM_CONS_Y=0.D0
!DO K=1,KTOTAL
!DO KK=1,KTOTAL
!    IF (K==KK) CYCLE
!RHS_MOMENTUM_CONS_X(K)=RHS_MOMENTUM_CONS_X(K)-MASS*((P(P_ZERO(K),RHO_ZERO(K),RHO(K))&
!                      +P(P_ZERO(KK),RHO_ZERO(KK),RHO(KK)))/(RHO(K)*RHO(KK))&
!                      +VIS_TERM(K))*WEIGHTING_FUNCTION_X(HBAR(H(RHO(K),RHO_ZERO(K)),H(RHO(KK),RHO_ZERO(KK))),&
!                    XPOS(K)-XPOS(KK),DISTANCE(XPOS(K),XPOS(KK),YPOS(K),YPOS(KK)))+F_BOUND_X(K)
!RHS_MOMENTUM_CONS_Y(K)=RHS_MOMENTUM_CONS_Y(K)-MASS*((P(P_ZERO(K),RHO_ZERO(K),RHO(K))&
!                    +P(P_ZERO(KK),RHO_ZERO(KK),RHO(KK)))/(RHO(K)*RHO(KK))&
!                    +VIS_TERM(K))*WEIGHTING_FUNCTION_Y(HBAR(H(RHO(K),RHO_ZERO(K)),H(RHO(KK),RHO_ZERO(KK))),&
!                    YPOS(K)-YPOS(KK),DISTANCE(XPOS(K),XPOS(KK),YPOS(K),YPOS(KK)))+F_BOUND_Y(K)
!END DO
!END DO
!------------------------------------------------------------------------------
RHS_MOMENTUM_CONS_X=0.D0
RHS_MOMENTUM_CONS_Y=0.D0
DO K=1,KTOTAL
DO KK=1,KTOTAL
    IF (K/=KK) THEN
RHS_MOMENTUM_CONS_X(K)=RHS_MOMENTUM_CONS_X(K)+MASS*(P(P_ZERO(K),RHO_ZERO(K),RHO(K))/RHO(K)**2.D0&
                      +P(P_ZERO(KK),RHO_ZERO(KK),RHO(KK))/RHO(KK)**2.D0&
                      +VIS_TERM(K))&
                      *WEIGHTING_FUNCTION_X(HBAR(H(RHO(K),RHO_ZERO(K)),H(RHO(KK),RHO_ZERO(KK))),&
                       XPOS(K)-XPOS(KK),DISTANCE(XPOS(K),XPOS(KK),YPOS(K),YPOS(KK)))
RHS_MOMENTUM_CONS_Y(K)=RHS_MOMENTUM_CONS_Y(K)+MASS*(P(P_ZERO(K),RHO_ZERO(K),RHO(K))/RHO(K)**2.D0&
                      +P(P_ZERO(KK),RHO_ZERO(KK),RHO(KK))/RHO(KK)**2.D0&
                      +VIS_TERM(K))&
                      *WEIGHTING_FUNCTION_Y(HBAR(H(RHO(K),RHO_ZERO(K)),H(RHO(KK),RHO_ZERO(KK))),&
                       YPOS(K)-YPOS(KK),DISTANCE(XPOS(K),XPOS(KK),YPOS(K),YPOS(KK)))
     END IF
END DO
RHS_MOMENTUM_CONS_X(K)=-RHS_MOMENTUM_CONS_X(K)+F_BOUND_X(K)
RHS_MOMENTUM_CONS_Y(K)=-RHS_MOMENTUM_CONS_Y(K)+F_BOUND_Y(K)
END DO
!------------------------------------------------------------------------------
END SUBROUTINE RHS_MOMENTUM_CONSERVATION
!==============================================================================
SUBROUTINE VISCOSITY_TERM (KTOTAL,MU,RHO,RHO_ZERO,U,V,XPOS,YPOS,VIS_TERM)
USE MODULE1
IMPLICIT NONE
INTEGER                     :: K,KK,KTOTAL
REAL(8),DIMENSION(KTOTAL)   :: U,V,RHO,XPOS,YPOS,MU,VIS_TERM,RHO_ZERO
REAL(8)                     :: H,HBAR,DISTANCE
!------------------------------------------------------------------------------
DO K=1,KTOTAL
DO KK=1,KTOTAL
    IF (K/=KK) THEN
VIS_TERM(K)=-16.D0*(MU(K)*MU(KK)/(RHO(K)*RHO(KK)*(MU(K)+MU(KK))))&
                  *(((U(K)-U(KK))*(XPOS(K)-XPOS(KK))+(V(K)-V(KK))&
                  *(YPOS(K)-YPOS(KK)))/(DISTANCE(XPOS(K),XPOS(KK),YPOS(K),YPOS(KK))&
                  **2.D0+EPSILON*HBAR(H(RHO(K),RHO_ZERO(K)),H(RHO(KK),RHO_ZERO(KK)))**2.D0))
    END IF
END DO
END DO
!------------------------------------------------------------------------------
END SUBROUTINE VISCOSITY_TERM
!==============================================================================
SUBROUTINE BOUNDARY_FORCE (KBTOTAL,KTOTAL,BOUND_PARTICLE_COUNTER,NONBOUND_PARTICLE_COUNTER,MASS,RHO,RHO_ZERO,&
                           V_MAX,XPOS,YPOS,F_BOUND_X,F_BOUND_Y)
USE MODULE1
IMPLICIT NONE
INTEGER                            :: K,KK,KTOTAL,KBTOTAL
INTEGER,DIMENSION(KBTOTAL)         :: BOUND_PARTICLE_COUNTER
INTEGER,DIMENSION(KTOTAL-KBTOTAL)  :: NONBOUND_PARTICLE_COUNTER
REAL(8)                            :: WENDLAND_WEIGHTING_FUNCTION,DISTANCE,H,HBAR,MASS
REAL(8),DIMENSION(KTOTAL)          :: XPOS,YPOS,F_BOUND_X,F_BOUND_Y,RHO,RHO_ZERO,V_MAX
!------------------------------------------------------------------------------
F_BOUND_X=0.D0
F_BOUND_Y=0.D0
!F_BOUND_X(NONBOUND_PARTICLE_COUNTER(1:KTOTAL-KBTOTAL))=0.01D0
!F_BOUND_Y(NONBOUND_PARTICLE_COUNTER(1:KTOTAL-KBTOTAL))=0.01D0
DO K=1,KTOTAL-KBTOTAL
DO KK=1,KBTOTAL
F_BOUND_X(NONBOUND_PARTICLE_COUNTER(K))=F_BOUND_X(NONBOUND_PARTICLE_COUNTER(K))+(V_MAX(NONBOUND_PARTICLE_COUNTER(K))**2.D0)*&  !------>>
WENDLAND_WEIGHTING_FUNCTION&
(HBAR(H(RHO(NONBOUND_PARTICLE_COUNTER(K)),RHO_ZERO(NONBOUND_PARTICLE_COUNTER(K))),&
H(RHO(BOUND_PARTICLE_COUNTER(KK)),RHO_ZERO(BOUND_PARTICLE_COUNTER(KK)))),&
DISTANCE(XPOS(NONBOUND_PARTICLE_COUNTER(K)),XPOS(BOUND_PARTICLE_COUNTER(KK)),&
YPOS(NONBOUND_PARTICLE_COUNTER(K)),YPOS(BOUND_PARTICLE_COUNTER(KK))))&
*(XPOS(NONBOUND_PARTICLE_COUNTER(K))-XPOS(BOUND_PARTICLE_COUNTER(KK)))&
/DISTANCE(XPOS(NONBOUND_PARTICLE_COUNTER(K)),XPOS(BOUND_PARTICLE_COUNTER(KK)),&
YPOS(NONBOUND_PARTICLE_COUNTER(K)),YPOS(BOUND_PARTICLE_COUNTER(KK)))**2.D0
!
F_BOUND_Y(NONBOUND_PARTICLE_COUNTER(K))=F_BOUND_Y(NONBOUND_PARTICLE_COUNTER(K))+(V_MAX(NONBOUND_PARTICLE_COUNTER(K))**2.D0)*&    !------>>
WENDLAND_WEIGHTING_FUNCTION&
(HBAR(H(RHO(NONBOUND_PARTICLE_COUNTER(K)),RHO_ZERO(NONBOUND_PARTICLE_COUNTER(K))),&
H(RHO(BOUND_PARTICLE_COUNTER(KK)),RHO_ZERO(BOUND_PARTICLE_COUNTER(KK)))),&
DISTANCE(XPOS(NONBOUND_PARTICLE_COUNTER(K)),XPOS(BOUND_PARTICLE_COUNTER(KK)),&
YPOS(NONBOUND_PARTICLE_COUNTER(K)),YPOS(BOUND_PARTICLE_COUNTER(KK))))*&
(YPOS(NONBOUND_PARTICLE_COUNTER(K))-YPOS(BOUND_PARTICLE_COUNTER(KK)))&
/DISTANCE(XPOS(NONBOUND_PARTICLE_COUNTER(K)),XPOS(BOUND_PARTICLE_COUNTER(KK)),&
YPOS(NONBOUND_PARTICLE_COUNTER(K)),YPOS(BOUND_PARTICLE_COUNTER(KK)))**2.D0
END DO
END DO
!------------------------------------------------------------------------------
!F_BOUND_X=0.D0
!F_BOUND_Y=0.D0
!DO K=1,KTOTAL-KBTOTAL
!DO KK=1,KBTOTAL
!F_BOUND_X(NONBOUND_PARTICLE_COUNTER(K))=F_BOUND_X(NONBOUND_PARTICLE_COUNTER(K))+(V_MAX(NONBOUND_PARTICLE_COUNTER(K))**2.D0)*&  !------>>
!(2.D0*MASS(BOUND_PARTICLE_COUNTER(KK))/(MASS(NONBOUND_PARTICLE_COUNTER(K))+&
!MASS(BOUND_PARTICLE_COUNTER(KK))))*WENDLAND_WEIGHTING_FUNCTION&
!(HBAR(H(RHO(NONBOUND_PARTICLE_COUNTER(K)),RHO_ZERO(NONBOUND_PARTICLE_COUNTER(K))),&
!H(RHO(BOUND_PARTICLE_COUNTER(KK)),RHO_ZERO(BOUND_PARTICLE_COUNTER(KK)))),&
!DISTANCE(XPOS(NONBOUND_PARTICLE_COUNTER(K)),XPOS(BOUND_PARTICLE_COUNTER(KK)),&
!YPOS(NONBOUND_PARTICLE_COUNTER(K)),YPOS(BOUND_PARTICLE_COUNTER(KK))))&
!*(XPOS(NONBOUND_PARTICLE_COUNTER(K))-XPOS(BOUND_PARTICLE_COUNTER(KK)))&
!/(DISTANCE(XPOS(NONBOUND_PARTICLE_COUNTER(K)),XPOS(BOUND_PARTICLE_COUNTER(KK)),&
!YPOS(NONBOUND_PARTICLE_COUNTER(K)),YPOS(BOUND_PARTICLE_COUNTER(KK)))**2.D0-&
!DX*DISTANCE(XPOS(NONBOUND_PARTICLE_COUNTER(K)),XPOS(BOUND_PARTICLE_COUNTER(KK)),&
!YPOS(NONBOUND_PARTICLE_COUNTER(K)),YPOS(BOUND_PARTICLE_COUNTER(KK))))
!!
!F_BOUND_Y(NONBOUND_PARTICLE_COUNTER(K))=F_BOUND_Y(NONBOUND_PARTICLE_COUNTER(K))+(V_MAX(NONBOUND_PARTICLE_COUNTER(K))**2.D0)*&    !------>>
!(2.D0*MASS(BOUND_PARTICLE_COUNTER(KK))/(MASS(NONBOUND_PARTICLE_COUNTER(K))+&
!MASS(BOUND_PARTICLE_COUNTER(KK))))*WENDLAND_WEIGHTING_FUNCTION&
!(HBAR(H(RHO(NONBOUND_PARTICLE_COUNTER(K)),RHO_ZERO(NONBOUND_PARTICLE_COUNTER(K))),&
!H(RHO(BOUND_PARTICLE_COUNTER(KK)),RHO_ZERO(BOUND_PARTICLE_COUNTER(KK)))),&
!DISTANCE(XPOS(NONBOUND_PARTICLE_COUNTER(K)),XPOS(BOUND_PARTICLE_COUNTER(KK)),&
!YPOS(NONBOUND_PARTICLE_COUNTER(K)),YPOS(BOUND_PARTICLE_COUNTER(KK))))*&
!(YPOS(NONBOUND_PARTICLE_COUNTER(K))-YPOS(BOUND_PARTICLE_COUNTER(KK)))&
!/(DISTANCE(XPOS(NONBOUND_PARTICLE_COUNTER(K)),XPOS(BOUND_PARTICLE_COUNTER(KK)),&
!YPOS(NONBOUND_PARTICLE_COUNTER(K)),YPOS(BOUND_PARTICLE_COUNTER(KK)))**2.D0-&
!DX*DISTANCE(XPOS(NONBOUND_PARTICLE_COUNTER(K)),XPOS(BOUND_PARTICLE_COUNTER(KK)),&
!YPOS(NONBOUND_PARTICLE_COUNTER(K)),YPOS(BOUND_PARTICLE_COUNTER(KK))))
!END DO
!END DO
!------------------------------------------------------------------------------
END SUBROUTINE BOUNDARY_FORCE
!==============================================================================
SUBROUTINE RHS_ENERGY_CONSERVATION (KTOTAL,MASS,RHO,RHO_ZERO,CP,T,U,V,K_COND,VIS_TERM,&
                                    XPOS,YPOS,RHS_ENERGY_CONS)
USE MODULE1
IMPLICIT NONE
INTEGER                     :: K,KK,KTOTAL
REAL(8),DIMENSION(KTOTAL)   :: U,V,RHO,XPOS,YPOS,T,K_COND,CP,RHO_ZERO
REAL(8),DIMENSION(KTOTAL)   :: TERM1,TERM2,TERM3,VIS_TERM,RHS_ENERGY_CONS
REAL(8)                     :: WEIGHTING_FUNCTION_Y,WEIGHTING_FUNCTION_X,DISTANCE,H,HBAR,MASS
!------------------------------------------------------------------------------
TERM1=0.D0
TERM2=0.D0
TERM3=0.D0
DO K=1,KTOTAL
DO KK=1,KTOTAL
    IF (K/=KK) THEN
TERM1(K)=TERM1(K)-(4.D0*MASS*K_COND(K)*K_COND(KK)*(T(K)-T(KK))&
                       *((XPOS(K)-XPOS(KK))*WEIGHTING_FUNCTION_X(HBAR(H(RHO(K),RHO_ZERO(K)),H(RHO(KK),RHO_ZERO(KK))),&
                       XPOS(K)-XPOS(KK),DISTANCE(XPOS(K),XPOS(KK),YPOS(K),YPOS(KK)))&
                        +(YPOS(K)-YPOS(KK))*WEIGHTING_FUNCTION_Y(HBAR(H(RHO(K),RHO_ZERO(K)),H(RHO(KK),RHO_ZERO(KK))),&
                        YPOS(K)-YPOS(KK),DISTANCE(XPOS(K),XPOS(KK),YPOS(K),YPOS(KK)))))&
                        /(RHO(K)*RHO(KK)*(K_COND(K)+K_COND(KK))*(DISTANCE(XPOS(K),XPOS(KK),YPOS(K),YPOS(KK))&
                        **2.D0+EPSILON*HBAR(H(RHO(K),RHO_ZERO(K)),H(RHO(KK),RHO_ZERO(KK)))**2.D0))

TERM2(K)=TERM2(K)-0.5D0*MASS*VIS_TERM(K)*((U(K)-U(KK))*WEIGHTING_FUNCTION_X(HBAR(H(RHO(K),RHO_ZERO(K)),&
         H(RHO(KK),RHO_ZERO(KK))),&
         XPOS(K)-XPOS(KK),DISTANCE(XPOS(K),XPOS(KK),YPOS(K),YPOS(KK)))&
        +(V(K)-V(KK))*WEIGHTING_FUNCTION_Y(HBAR(H(RHO(K),RHO_ZERO(K)),H(RHO(KK),RHO_ZERO(KK))),&
        YPOS(K)-YPOS(KK),DISTANCE(XPOS(K),XPOS(KK),YPOS(K),YPOS(KK))))!------>>>
END IF
TERM3(K)=-H_CONV*(T(K)-T_INFTY)
END DO
END DO
RHS_ENERGY_CONS=(1.D0/CP)*(TERM1+TERM2+TERM3)
!------------------------------------------------------------------------------
END SUBROUTINE RHS_ENERGY_CONSERVATION
!==============================================================================
FUNCTION H(RHOI,RHO_ZEROI)
USE MODULE1
IMPLICIT NONE
REAL(8)              :: H,RHOI,RHO_ZEROI
!------------------------------------------------------------------------------
H=H_ZERO*(RHOI/RHO_ZEROI)**(-1.D0/3.D0)
!------------------------------------------------------------------------------
END FUNCTION H
!==============================================================================
FUNCTION HBAR (HI,HJ)
USE MODULE1
IMPLICIT NONE
REAL(8)        :: HI,HJ,HBAR
!------------------------------------------------------------------------------
HBAR=0.5D0*(HI+HJ)
!------------------------------------------------------------------------------
END FUNCTION HBAR
!==============================================================================
FUNCTION DISTANCE (XI,XJ,YI,YJ)
USE MODULE1
IMPLICIT NONE
REAL(8)        :: XI,XJ,YI,YJ,DISTANCE
!------------------------------------------------------------------------------
DISTANCE=DSQRT((XI-XJ)**2.D0+(YI-YJ)**2.D0)
!------------------------------------------------------------------------------
END FUNCTION DISTANCE
!==============================================================================
FUNCTION P(P_ZEROI,RHO_ZEROI,RHOI)
USE MODULE1
IMPLICIT NONE
REAL(8)        :: P,RHOI,P_ZEROI,RHO_ZEROI
!------------------------------------------------------------------------------
P=(P_ZEROI/RHO_ZEROI)*RHOI
!------------------------------------------------------------------------------
END FUNCTION P
!==============================================================================
SUBROUTINE STRAIN_RATE_TENSOR (KTOTAL,MASS,RHO,RHO_ZERO,U,V,XPOS,YPOS,EPSILON_11,EPSILON_12,EPSILON_21,EPSILON_22)
USE MODULE1
IMPLICIT NONE
INTEGER                    :: K,KK,KTOTAL
REAL(8)                    :: H,HBAR,DISTANCE,WEIGHTING_FUNCTION_X,WEIGHTING_FUNCTION_Y,MASS
REAL(8),DIMENSION(KTOTAL)  :: RHO,U,V,EPSILON_11,EPSILON_12,EPSILON_21,EPSILON_22
REAL(8),DIMENSION(KTOTAL)  :: XPOS,YPOS,RHO_ZERO
!------------------------------------------------------------------------------
EPSILON_11=0.D0
EPSILON_12=0.D0
EPSILON_21=0.D0
EPSILON_22=0.D0
DO K=1,KTOTAL
DO KK=1,KTOTAL
EPSILON_11(K)=EPSILON_11(K)-(MASS/RHO(KK))*(U(K)-U(KK))&
             *WEIGHTING_FUNCTION_X(HBAR(H(RHO(K),RHO_ZERO(K)),H(RHO(KK),RHO_ZERO(KK))),XPOS(K)-XPOS(KK),&
             DISTANCE(XPOS(K),XPOS(KK),YPOS(K),YPOS(KK)))
EPSILON_12(K)=EPSILON_12(K)-0.5D0*(MASS/RHO(KK))*(U(K)-U(KK))&
             *WEIGHTING_FUNCTION_Y(HBAR(H(RHO(K),RHO_ZERO(K)),H(RHO(KK),RHO_ZERO(KK))),YPOS(K)-YPOS(KK),&
             DISTANCE(XPOS(K),XPOS(KK),YPOS(K),YPOS(KK)))+(V(K)-V(KK))&
             *WEIGHTING_FUNCTION_X(HBAR(H(RHO(K),RHO_ZERO(K)),H(RHO(KK),RHO_ZERO(KK))),XPOS(K)-XPOS(KK),&
             DISTANCE(XPOS(K),XPOS(KK),YPOS(K),YPOS(KK)))
EPSILON_21(K)=EPSILON_21(K)-0.5D0*(MASS/RHO(KK))*(U(K)-U(KK))&
             *WEIGHTING_FUNCTION_Y(HBAR(H(RHO(K),RHO_ZERO(K)),H(RHO(KK),RHO_ZERO(KK))),YPOS(K)-YPOS(KK),&
             DISTANCE(XPOS(K),XPOS(KK),YPOS(K),YPOS(KK)))+(V(K)-V(KK))&
             *WEIGHTING_FUNCTION_X(HBAR(H(RHO(K),RHO_ZERO(K)),H(RHO(KK),RHO_ZERO(KK))),XPOS(K)-XPOS(KK),&
             DISTANCE(XPOS(K),XPOS(KK),YPOS(K),YPOS(KK)))
EPSILON_22(K)=EPSILON_22(K)-(MASS/RHO(KK))*(V(K)-V(KK))&
             *WEIGHTING_FUNCTION_Y(HBAR(H(RHO(K),RHO_ZERO(K)),H(RHO(KK),RHO_ZERO(KK))),XPOS(K)-XPOS(KK),&
             DISTANCE(XPOS(K),XPOS(KK),YPOS(K),YPOS(KK)))
END DO
END DO
!------------------------------------------------------------------------------
END SUBROUTINE STRAIN_RATE_TENSOR
!==============================================================================
FUNCTION EPSILON_E (EPSILONI_11,EPSILONI_12,EPSILONI_21,EPSILONI_22)
IMPLICIT NONE
REAL(8)                :: EPSILONI_11,EPSILONI_12,EPSILONI_21,EPSILONI_22,EPSILON_E
!------------------------------------------------------------------------------
EPSILON_E=DSQRT((2.D0/3.D0)*(EPSILONI_11**2.D0+2.D0*EPSILONI_12*EPSILONI_21+EPSILONI_22**2.D0))
!------------------------------------------------------------------------------
END FUNCTION EPSILON_E
!==============================================================================
FUNCTION ZENER_HOLLOMON (EPSILONI_E,TI)
USE MODULE1
IMPLICIT NONE
REAL(8)       :: ZENER_HOLLOMON,EPSILONI_E,TI
!------------------------------------------------------------------------------
ZENER_HOLLOMON=EPSILONI_E*DEXP(Q_ACT_EN/(R_GAS_CONST*TI))
!------------------------------------------------------------------------------
END FUNCTION ZENER_HOLLOMON
!==============================================================================
FUNCTION SIGMA_E (ZENER_HOLLOMONI)
USE MODULE1
IMPLICIT NONE
REAL(8)        :: SIGMA_E,ZENER_HOLLOMONI
!------------------------------------------------------------------------------
SIGMA_E=SIGMA_R*ABS(ASINH((ZENER_HOLLOMONI/BETA)**(INVM)))
!------------------------------------------------------------------------------
END FUNCTION SIGMA_E
!==============================================================================
!SUBROUTINE VISCOSITY_1 (KTOTAL,EPSILON_11,EPSILON_12,EPSILON_21,EPSILON_22,T,MU)
!USE MODULE1
!IMPLICIT NONE
!INTEGER                      :: K,KTOTAL
!REAL(8),DIMENSION(KTOTAL)    :: MU,EPSILON_11,EPSILON_12,EPSILON_21,EPSILON_22,T
!REAL(8)                      :: ZENER_HOLLOMON,SIGMA_E,EPSILON_E
!!------------------------------------------------------------------------------
!DO K=1,KTOTAL
!MU(K)=SIGMA_E(ZENER_HOLLOMON(EPSILON_E(EPSILON_11(K),EPSILON_12(K),EPSILON_21(K),EPSILON_22(K)),T(K)))/&
!      3.D0*(EPSILON_E(EPSILON_11(K),EPSILON_12(K),EPSILON_21(K),EPSILON_22(K))+EPSILON)
!END DO
!!------------------------------------------------------------------------------
!END SUBROUTINE VISCOSITY_1
!==============================================================================
SUBROUTINE VISCOSITY_2 (KPIN,KTOTAL,FS,MU)
USE MODULE1
IMPLICIT NONE
INTEGER                   :: K,KPIN,KTOTAL
REAL(8),DIMENSION(KTOTAL) :: FS,MU
!------------------------------------------------------------------------------
DO K=1,KTOTAL
    IF ((FS(K)>0.D0).AND.(FS(K)<=0.67D0)) THEN
        MU(K)=MU_MELTING*(1.D0-FS(K))**(DLOG(MI)/DLOG(0.33D0))
    ELSE IF ((FS(K)>0.67D0).AND.(FS(K)<1.D0)) THEN
        MU(K)=MS*MU_MELTING*(1.D0-FS(K))**(DLOG(MI+MS)/DLOG(0.33D0))
    ELSE
        MU(K)=MS
    END IF
END DO
!------------------------------------------------------------------------------
END SUBROUTINE VISCOSITY_2
!==============================================================================
SUBROUTINE VISCOSITY_3 (KPIN,KTOTAL,T,MU)
USE MODULE1
IMPLICIT NONE
INTEGER                   :: K,KPIN,KTOTAL
REAL(8),DIMENSION(KTOTAL) :: T,MU
!------------------------------------------------------------------------------
MU(1:KPIN)=MS
DO K=KPIN+1,KTOTAL
    IF (T(K)<T_SOLIDUS) THEN
        MU(K)=1000.D0
    ELSE IF (T(K)>T_LIQUIDUS) THEN
        MU(K)=0.001D0
    ELSE
        MU(K)=0.001D0*DEXP(-13.82D0*((T(K)-T_LIQUIDUS)/(T_LIQUIDUS-T_SOLIDUS)))
    END IF
END DO
!------------------------------------------------------------------------------
END SUBROUTINE VISCOSITY_3
!==============================================================================
SUBROUTINE SOLID_FRACTION (KTOTAL,KPIN,T,FS)
USE MODULE1
IMPLICIT NONE
INTEGER                     :: K,KPIN,KTOTAL
REAL(8),DIMENSION(KTOTAL)   :: T,FS
!------------------------------------------------------------------------------
DO K=1,KTOTAL
    IF (K<=KPIN) THEN
    FS(K)=FS_PIN
    ELSE
    FS(K)=(T_LIQUIDUS-T(K))/(T_LIQUIDUS-T_SOLIDUS)
    END IF
END DO
!------------------------------------------------------------------------------
END SUBROUTINE SOLID_FRACTION
!==============================================================================
FUNCTION WENDLAND_WEIGHTING_FUNCTION(HBARI,DISTANCEI)
USE MODULE1
IMPLICIT NONE
REAL(8)                     :: Q,WENDLAND_WEIGHTING_FUNCTION,HBARI,DISTANCEI
!------------------------------------------------------------------------------
Q=DISTANCEI/HBARI
IF ((Q>=0.D0).AND.(Q<=2.D0)) THEN
WENDLAND_WEIGHTING_FUNCTION=W3*(1.D0+1.5D0*Q)*(2.D0-Q)**3.D0
ELSE
WENDLAND_WEIGHTING_FUNCTION=0.D0
END IF
!------------------------------------------------------------------------------
END FUNCTION WENDLAND_WEIGHTING_FUNCTION
!==============================================================================
SUBROUTINE MAXIMUM_VELOCITY (KBTOTAL,KTOTAL,BOUND_PARTICLE_COUNTER,&
                             NONBOUND_PARTICLE_COUNTER,U,V,V_MAX)
USE MODULE1
IMPLICIT NONE
INTEGER                           :: K,KK,KTOTAL,KBTOTAL
INTEGER,DIMENSION(KBTOTAL)        :: BOUND_PARTICLE_COUNTER
INTEGER,DIMENSION(KTOTAL-KBTOTAL) :: NONBOUND_PARTICLE_COUNTER
REAL(8),DIMENSION(KTOTAL-KBTOTAL) :: V_MAX
REAL(8),DIMENSION(KTOTAL)         :: U,V
!------------------------------------------------------------------------------
V_MAX=0.D0
DO K=1,KTOTAL-KBTOTAL
DO KK=1,KBTOTAL
V_MAX(NONBOUND_PARTICLE_COUNTER(K))=MAX(V_MAX(NONBOUND_PARTICLE_COUNTER(K)),&
DSQRT((U(NONBOUND_PARTICLE_COUNTER(K))-U(BOUND_PARTICLE_COUNTER(KK)))**2.D0+&
(V(NONBOUND_PARTICLE_COUNTER(K))-V(BOUND_PARTICLE_COUNTER(KK)))**2.D0))
END DO
END DO
!------------------------------------------------------------------------------
END SUBROUTINE MAXIMUM_VELOCITY
!==============================================================================
SUBROUTINE INITIAL_CONDITIONS (KPIN,KTOTAL,X,Y,U,V,T,PHASE,MASS,RHO_ZERO,K_COND,CP,&
                               P_ZERO,RHO,RHS_MASS_CONS,RHS_MOMENTUM_CONS_X,&
                               RHS_MOMENTUM_CONS_Y,RHS_ENERGY_CONS,XC,YC)
USE MODULE1
IMPLICIT NONE
INTEGER                    :: K,KPIN,KTOTAL
REAL(8),DIMENSION(KTOTAL)  :: U,V,T,RHO_ZERO,K_COND,CP,P_ZERO,RHO,X,Y
REAL(8),DIMENSION(KTOTAL)  :: RHS_MASS_CONS,RHS_MOMENTUM_CONS_X
REAL(8),DIMENSION(KTOTAL)  :: RHS_MOMENTUM_CONS_Y,RHS_ENERGY_CONS,PHASE
REAL(8)                    :: MASS,XC,YC,DUMMY
!------------------------------------------------------------------------------
IF (ICFLAG) THEN
   OPEN (1001,FILE="FLD 392000.DAT")
   READ (1001,*)
   READ (1001,*)
   DO K=1,KTOTAL
       READ (1001,*)  X(K),Y(K),RHO(K),U(K),V(K),T(K),PHASE(K)
   END DO
   CLOSE(1001)
   MASS=MASS0
   XC=XCI+KSTART*DT*UT
   YC=YCI+KSTART*DT*VT
   DO K=1,KTOTAL
      IF (K<=KPIN) THEN
         RHO_ZERO(K)=RHO_ZERO_TOOL
         K_COND(K)=K_COND_ZERO_TOOL
         CP(K)=CP_ZERO_TOOL
         P_ZERO(K)=P_ZERO_TOOL
      ELSE
         RHO_ZERO(K)=RHO_ZERO_PLATE
         K_COND(K)=K_COND_ZERO_PLATE
         CP(K)=CP_ZERO_PLATE
         P_ZERO(K)=P_ZERO_PLATE
      END IF
   END DO
ELSE
U=0.D0
V=0.D0
T=T_INFTY
MASS=MASS0
DO K=1,KTOTAL
IF (K<=KPIN) THEN
RHO_ZERO(K)=RHO_ZERO_TOOL
K_COND(K)=K_COND_ZERO_TOOL
CP(K)=CP_ZERO_TOOL
P_ZERO(K)=P_ZERO_TOOL
ELSE
RHO_ZERO(K)=RHO_ZERO_PLATE
K_COND(K)=K_COND_ZERO_PLATE
CP(K)=CP_ZERO_PLATE
P_ZERO(K)=P_ZERO_PLATE
END IF
RHO(K)=RHO_ZERO(K)
END DO
END IF
RHS_MASS_CONS=0.D0
RHS_MOMENTUM_CONS_X=0.D0
RHS_MOMENTUM_CONS_Y=0.D0
RHS_ENERGY_CONS=0.D0
!------------------------------------------------------------------------------
END SUBROUTINE INITIAL_CONDITIONS
!==============================================================================
SUBROUTINE EULER_INTEGRATION (KBTOTAL,KTOTAL,NONBOUND_PARTICLE_COUNTER,&
                              RHS_MASS_CONS,RHS_MOMENTUM_CONS_X,&
                              RHS_MOMENTUM_CONS_Y,RHS_ENERGY_CONS,RHO,U,V,T)
USE MODULE1
IMPLICIT NONE
INTEGER                            :: K,KTOTAL,KBTOTAL
REAL(8),DIMENSION(KTOTAL)          :: RHO,RHO0,U,U0,V,V0,T,T0,RHS_MASS_CONS,&
                                      RHS_MOMENTUM_CONS_X,RHS_MOMENTUM_CONS_Y,RHS_ENERGY_CONS
INTEGER,DIMENSION(KTOTAL-KBTOTAL)  :: NONBOUND_PARTICLE_COUNTER
!------------------------------------------------------------------------------
!T0=T
!RHO0=RHO
!U0=U
!V0=V
!U=0.D0
!V=0.D0
!T=T0+RHS_ENERGY_CONS*DT
!RHO=RHO0+RHS_MASS_CONS*DT
!DO K=1,KTOTAL-KBTOTAL
!!RHO(NONBOUND_PARTICLE_COUNTER(K))=RHO(NONBOUND_PARTICLE_COUNTER(K))+RHS_MASS_CONS(NONBOUND_PARTICLE_COUNTER(K))*DT
!U(NONBOUND_PARTICLE_COUNTER(K))=U0(NONBOUND_PARTICLE_COUNTER(K))+RHS_MOMENTUM_CONS_X(NONBOUND_PARTICLE_COUNTER(K))*DT
!V(NONBOUND_PARTICLE_COUNTER(K))=V0(NONBOUND_PARTICLE_COUNTER(K))+RHS_MOMENTUM_CONS_Y(NONBOUND_PARTICLE_COUNTER(K))*DT
!END DO
!------------------------------------------------------------------------------
T0=T
RHO0=RHO
U0=U
V0=V
T=T0+RHS_ENERGY_CONS*DT
RHO=RHO0+RHS_MASS_CONS*DT
DO K=1,KTOTAL-KBTOTAL
IF (DABS(RHS_MOMENTUM_CONS_X(NONBOUND_PARTICLE_COUNTER(K)))<=1.D-16) THEN
    U(NONBOUND_PARTICLE_COUNTER(K))=0.D0
    ELSE
    U(NONBOUND_PARTICLE_COUNTER(K))=U0(NONBOUND_PARTICLE_COUNTER(K))+RHS_MOMENTUM_CONS_X(NONBOUND_PARTICLE_COUNTER(K))*DT
END IF
IF (DABS(RHS_MOMENTUM_CONS_Y(NONBOUND_PARTICLE_COUNTER(K)))<=1.D-16) THEN
    V(NONBOUND_PARTICLE_COUNTER(K))=0.D0
    ELSE
    V(NONBOUND_PARTICLE_COUNTER(K))=V0(NONBOUND_PARTICLE_COUNTER(K))+RHS_MOMENTUM_CONS_Y(NONBOUND_PARTICLE_COUNTER(K))*DT
END IF
END DO
!------------------------------------------------------------------------------
END SUBROUTINE EULER_INTEGRATION
!==============================================================================
SUBROUTINE PARTICLE_POSITION (KBTOTAL,KTOTAL,NONBOUND_PARTICLE_COUNTER,X,Y,U,V)
USE MODULE1
IMPLICIT NONE
INTEGER                            :: K,KBTOTAL,KTOTAL
REAL(8),DIMENSION(KTOTAL)          :: X,Y,X0,Y0,U,V
INTEGER,DIMENSION(KTOTAL-KBTOTAL)  :: NONBOUND_PARTICLE_COUNTER
!------------------------------------------------------------------------------
X0=X
Y0=Y
DO K=1,KTOTAL-KBTOTAL
X(NONBOUND_PARTICLE_COUNTER(K))=X0(NONBOUND_PARTICLE_COUNTER(K))+U(NONBOUND_PARTICLE_COUNTER(K))*DT
Y(NONBOUND_PARTICLE_COUNTER(K))=Y0(NONBOUND_PARTICLE_COUNTER(K))+V(NONBOUND_PARTICLE_COUNTER(K))*DT
END DO
!------------------------------------------------------------------------------
END SUBROUTINE PARTICLE_POSITION
!==============================================================================
SUBROUTINE PARTICLE_PHASE (KPIN,KTOTAL,Y,PHASE)
USE MODULE1
IMPLICIT NONE
INTEGER                           :: K,KPIN,KTOTAL
REAL(8),DIMENSION(KTOTAL)         :: Y,PHASE
!------------------------------------------------------------------------------
DO K=1,KTOTAL
IF (K<=KPIN) THEN
    PHASE(K)=0.D0
    ELSE IF (Y(K)<=0.D0) THEN
    PHASE(K)=1.D0
    ELSE
    PHASE(K)=2.D0
END IF
END DO
!------------------------------------------------------------------------------
END SUBROUTINE PARTICLE_PHASE
!==============================================================================
SUBROUTINE PHASE_STATISTIC (KTOTAL,TIME_STEP,PHASE,Y)
USE MODULE1
IMPLICIT NONE
INTEGER                           :: K,KTOTAL,COUNTER_PHASE1,COUNTER_PHASE2,TIME_STEP
REAL(8),DIMENSION(KTOTAL)         :: Y,PHASE
!------------------------------------------------------------------------------
OPEN (10001,FILE="PHASE1_STATISTIC.PLT")
OPEN (10002,FILE="PHASE2_STATISTIC.PLT")
COUNTER_PHASE1=0
COUNTER_PHASE2=0
DO K=1,KTOTAL
IF (((PHASE(K)>0.9D0).AND.(PHASE(K)<1.1D0)).AND.(Y(K)>0.D0))  THEN
    COUNTER_PHASE1=COUNTER_PHASE1+1
ELSEIF ((PHASE(K)>1.9D0).AND.(Y(K)<0.D0))  THEN
    COUNTER_PHASE2=COUNTER_PHASE2+1
END IF
END DO
WRITE(10001,*) TIME_STEP*DT,COUNTER_PHASE1
WRITE(10002,*) TIME_STEP*DT,COUNTER_PHASE2
!------------------------------------------------------------------------------
END SUBROUTINE PHASE_STATISTIC
!==============================================================================
SUBROUTINE READ_TIME_SOLUTION (TIME_STEP,KTOTAL,Y,PHASE)
USE MODULE1
IMPLICIT NONE
INTEGER                             :: K,TIME_STEP,COUNTER,KTOTAL
REAL(8),DIMENSION(KTOTAL)           :: Y,PHASE
REAL(8)                             :: DUMMY
CHARACTER*7                         :: EXT
CHARACTER*3                         :: FN1
CHARACTER*18                        :: FNAME
!------------------------------------------------------------------------------
FN1='FLD'
COUNTER=0
!------------------------------------------------------------------------------
COUNTER=COUNTER+1
WRITE(EXT,'(I7)') TIME_STEP
FNAME=FN1//EXT//'.DAT'
!------------------------------------------------------------------------------
OPEN(COUNTER,FILE=FNAME,POSITION='REWIND')
READ(COUNTER,*)
READ(COUNTER,*)
DO K=1,KTOTAL
READ(COUNTER,*) DUMMY,Y(K),DUMMY,DUMMY,DUMMY,DUMMY,PHASE(K)
END DO
CLOSE(COUNTER)
!------------------------------------------------------------------------------
WRITE(*,*) '========================'
WRITE(*,*) 'PRINTING ON ',FNAME
WRITE(*,*) '========================'
!------------------------------------------------------------------------------
END SUBROUTINE READ_TIME_SOLUTION
!==============================================================================
