C ==================================================================== C
C SUBROUTINE K_COORDS_TRANSFORM: 
C             FUNCTION TO COMPUTE THE COORDINATE TRANSFORMATION MATRIX
C             BET EN THE GLOBAL AND LOCAL COORDINATES
C SUBROUTINE K_MATRIX_ZERO       : 
C             MATRIX OPERATION (A = 0)
C SUBROUTINE K_MATRIX_TRANSPOSE  : 
C             MATRIX OPERATION (B = A_T)
C SUBROUTINE K_MATRIX_PLUSSCALAR : 
C             MATRIX OPERATION (A = A + C * B)
C SUBROUTINE K_MATRIX_MULTIPLY   : 
C             MATRIX OPERATION (C = A * B)
C
C-----------------------------------------------------------------------------------------C
      SUBROUTINE UEL (RHS, AMATRX, SVARS, ENERGY, NDOFEL, NRHS, NSVARS,
     & PROPS, NPROPS, COORDS, MCRD, NNODE, U, DU, V, A, JTYPE, TIME,
     & DTIME, KSTEP, KINC, JELEM, PARAMS, NDLOAD, JDLTYP, ADLMAG,
     & PREDEF, NPREDF, LFLAGS, MLVARX, DDLMAG, MDLOAD, PNEWDT, JPROPS,
     & NJPRO, PERIOD)
       
      INCLUDE 'ABA_PARAM.INC'
       
      DIMENSION RHS(MLVARX,*), AMATRX(NDOFEL,NDOFEL), PROPS(*),
     & SVARS(*), ENERGY(8), COORDS(MCRD, NNODE), U(NDOFEL),
     & DU(MLVARX,*), V(NDOFEL), A(NDOFEL), TIME(2), PARAMS(*),
     & JDLTYP(MDLOAD,*), ADLMAG(MDLOAD,*), DDLMAG(MDLOAD,*),
     & PREDEF(2, NPREDF, NNODE), LFLAGS(*), JPROPS(*)
C       
      DIMENSION TRAC(MCRD,NRHS),TRAC_JACOB(MCRD,MCRD)
      DIMENSION CO_ORI(MCRD,NNODE), CO_CURR(MCRD,NNODE)
      DIMENSION U_loc(NDOFEL) , V_loc(NDOFEL) 
C       
      DIMENSION r_Matrix(MCRD,MCRD),r_Matrix_T(MCRD,MCRD)
      DIMENSION DeltaU_Loc(6),DeltaU_Loc_GP(MCRD),DeltaV_Loc_GP(MCRD)
      DIMENSION DeltaV_Loc(6)
      DIMENSION SF(2),SHAPE_N(3,6),  SHAPE_T(6,12)
      DIMENSION B_MATRIX(MCRD,NDOFEL), B_MATRIX_T(NDOFEL,MCRD)
C       
      DIMENSION RHS_TEMP1(MCRD,NRHS),      RHS_TEMP2(NDOFEL,NRHS)
      DIMENSION AMATRX_TEMP1(MCRD,MCRD),   AMATRX_TEMP2(MCRD,NDOFEL)
      DIMENSION AMATRX_TEMP3(MCRD,NDOFEL), AMATRX_TEMP4(NDOFEL,NDOFEL)
C
      INTEGER GP_N
C
      PARAMETER(N_ELEM=100000,NSDV=120,NGP=8)
C
C      Parameters required by UMAT.f
C
      PARAMETER (NDI=3,NSTATV=40,SSE=0.D0,SCD=0.D0,
     1           RPL=0.D0,DRPLDT=0.D0,TEMP=0.D0,DTEMP=0.D0,NSHR=1,
     2           CELENT=2.D0,LAYER=1,KSPT=1,NTENS=3)
C
C     Arrays to be used in UMAT.f
C
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1    DDSDDE(NTENS,NTENS),
     2    DDSDDT(NTENS),DRPLDE(NTENS),
     3    STRAN(NTENS),DSTRAN(NTENS),
     4    DPRED(1), DROT(3, 3), DFGRD0(3, 3),DFGRD1(3, 3),
     5    UPREDEF(1), UCOORDS(3)
C
      COMMON/KUSER/USRVAR(N_ELEM,NSDV,NGP)
C     -----------------------------------------------------------------------------------------C
      GP_N = 2
C     -----------------------------------------------------------------------------------------C
C     INITIALIZE MATRICES AND VECTORS
C 
      CALL K_MATRIX_ZERO(TRAC,MCRD,NRHS)
      CALL K_MATRIX_ZERO(TRAC_JACOB,MCRD,MCRD)
      CALL K_MATRIX_ZERO(RHS,NDOFEL,NRHS)
      CALL K_MATRIX_ZERO(AMATRX,NDOFEL,NDOFEL)
C      
      CALL K_VECTOR_ZERO(U_loc,NDOFEL)
      CALL K_VECTOR_ZERO(DeltaU_Loc,6)
      CALL K_VECTOR_ZERO(DeltaU_Loc_GP,MCRD)
      CALL K_VECTOR_ZERO(SF,2)
C      
      CALL K_MATRIX_ZERO(SHAPE_T,6,12)
      CALL K_VECTOR_ZERO(DeltaU_Loc,6)
C
C     MATRIX RELATED TO UMAT
C
      CALL K_VECTOR_ZERO(STRESS,NTENS)
      CALL K_VECTOR_ZERO(STRAN,NTENS)
      CALL K_VECTOR_ZERO(RELROT,1)
      CALL K_VECTOR_ZERO(DRELROT,1)
      CALL K_VECTOR_ZERO(DSTRAN,NTENS)
      CALL K_VECTOR_ZERO(DDSDDT,NTENS)
      CALL K_VECTOR_ZERO(DRPLDE,NTENS)
      CALL K_MATRIX_ZERO(DROT,3,3)
      CALL K_MATRIX_ZERO(DFGRD0,3,3)
      CALL K_MATRIX_ZERO(DFGRD1,3,3)
      CALL K_MATRIX_ZERO(DDSDDE,NTENS,NTENS)
C     -----------------------------------------------------------------------------------------C
C
      DO I = 1, MCRD
         CO_ORI(I,1)=COORDS(I,1)  
         CO_ORI(I,2)=COORDS(I,2)
         CO_ORI(I,3)=COORDS(I,3)
         CO_ORI(I,4)=COORDS(I,4)
      END DO
C
      DO I = 1, MCRD
         CO_CURR(I,1)=COORDS(I,1) + U(I)
         CO_CURR(I,2)=COORDS(I,2) + U(3+I)
         CO_CURR(I,3)=COORDS(I,3) + U(6+I)
         CO_CURR(I,4)=COORDS(I,4) + U(9+I)
      END DO
C     -----------------------------------------------------------------------------------------C
      CALL K_BASIS(CO_ORI,r_Matrix,U,U_Loc,MCRD,NNODE,NDOFEL)
      CALL K_BASIS(CO_ORI,r_Matrix,V,V_Loc,MCRD,NNODE,NDOFEL)
C
      SHAPE_T(1,1)      =   - 1.0D0
      SHAPE_T(1,10)     =     1.0D0
      SHAPE_T(2,2)      =   - 1.0D0
      SHAPE_T(2,11)     =     1.0D0
      SHAPE_T(3,3)      =   - 1.0D0
      SHAPE_T(3,12)     =     1.0D0
C      
      SHAPE_T(4,4)      =   - 1.0D0
      SHAPE_T(4,7)      =     1.0D0
      SHAPE_T(5,5)      =   - 1.0D0
      SHAPE_T(5,8)      =     1.0D0
      SHAPE_T(6,6)      =   - 1.0D0
      SHAPE_T(6,9)      =     1.0D0
C
      DeltaU_Loc = matmul(SHAPE_T, U_Loc)
      DeltaV_Loc = matmul(SHAPE_T, V_Loc)
C
C     -----------------------------------------------------------------------------------------C
C     -----------------------------------------------------------------------------------------C
C     -----------------------------------------------------------------------------------------C
C     DO CALCULATIONS AT GAUSS POINTS £¨GP = 2)
      DO I = 1, GP_N
C     -----------------------------------------------------------------------------------------C
         CALL K_SHAPE_FUN(I,SF)
C
         CALL K_MATRIX_ZERO(SHAPE_N,3,6)
         CALL K_MATRIX_ZERO(RHS_TEMP1,MCRD,NRHS)
         CALL K_MATRIX_ZERO(RHS_TEMP2,NDOFEL,NRHS)
         CALL K_MATRIX_ZERO(AMATRX_TEMP1,MCRD,MCRD)
         CALL K_MATRIX_ZERO(AMATRX_TEMP2,MCRD,NDOFEL)
         CALL K_MATRIX_ZERO(AMATRX_TEMP3,MCRD,NDOFEL)
         CALL K_MATRIX_ZERO(AMATRX_TEMP4,NDOFEL,NDOFEL)
C         
         SHAPE_N (1,1) = SF (1)
         SHAPE_N (2,2) = SF (1)
         SHAPE_N (3,3) = SF (1)
         SHAPE_N (1,4) = SF (2)
         SHAPE_N (2,5) = SF (2)
         SHAPE_N (3,6) = SF (2)
C 
         CALL K_VECTOR_ZERO(DeltaU_Loc_GP,MCRD)
         CALL K_VECTOR_ZERO(DeltaV_Loc_GP,MCRD)
C         
C         
         DeltaU_Loc_GP = matmul(SHAPE_N,DeltaU_Loc)
         DeltaV_Loc_GP = matmul(SHAPE_N,DeltaV_Loc)
C
         DeltaU_Loc_GP = DeltaU_Loc_GP 
         DeltaV_Loc_GP = DeltaV_Loc_GP 
         
C         write(6,*) 'DeltaU_Loc_GP',DeltaU_Loc_GP(1)
C     -----------------------------------------------------------------------------------------C
C        READ HISTORY VARIABLES 
C         
         GPT = I
         DO II = 1, NSTATV
             STATEV(II) = SVARS(NSDV*(GPT-1)+II)
         END DO
C
         STRAN(1) = DeltaU_Loc_GP(1)
         STRAN(2) = DeltaU_Loc_GP(2)
         STRAN(3) = DeltaU_Loc_GP(3)
C          
         DSTRAN(1) = STRAN(1) - SVARS(NSDV*(GPT-1)+41)
         DSTRAN(2) = STRAN(2) - SVARS(NSDV*(GPT-1)+42)
         DSTRAN(3) = STRAN(3) - SVARS(NSDV*(GPT-1)+43)
C     -----------------------------------------------------------------------------------------C
C        OBTIAN THE TRACTION AND TRACTION_STIFFNESS
C
         CALL UMAT_BOND(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL,
     1        DDSDDT, DRPLDE, DRPLDT, STRAN, DSTRAN, TIME, DTIME, TEMP,
     2        DTEMP, UPREDEF, DPRED, CMNAME, NDI, NSHR, NTENS, NSTATV,
     3        PROPS, NPROPS, UCOORDS, DROT, PNEWDT, CELENT, DFGRD0,
     4        DFGRD1, JELEM, GPT, LAYER, KSPT, KSTEP, KINC)
C
         ZETA = 0.001D0 
C
         DDSDDE(1,1)  = DDSDDE(1,1)  + ZETA / DTIME
	   STRESS(1)    = STRESS(1)    + ZETA  * DeltaV_Loc_GP(1)
C
         Trac_Jacob(1,1) = DDSDDE(1,1)
         Trac_Jacob(2,2) = DDSDDE(2,2)
         Trac_Jacob(3,3) = DDSDDE(3,3)
          
         Trac(1,1) = STRESS(1)
         Trac(2,1) = STRESS(2)
         Trac(3,1) = STRESS(3)
C     -----------------------------------------------------------------------------------------C
C        STORE HISTORY VARIABLES 
C         
         DO II = 1, NSTATV
            SVARS(NSDV*(GPT-1)+II) = STATEV(II)
         END DO
         SVARS(NSDV*(GPT-1)+41) = DeltaU_Loc_GP(1)
         SVARS(NSDV*(GPT-1)+42) = DeltaU_Loc_GP(2)
         SVARS(NSDV*(GPT-1)+43) = DeltaU_Loc_GP(3)
C     -----------------------------------------------------------------------------------------C
         B_MATRIX = matmul(SHAPE_N,SHAPE_T) 
         CALL K_MATRIX_TRANSPOSE(B_MATRIX,B_MATRIX_T,MCRD,NDOFEL)
         CALL K_MATRIX_TRANSPOSE(r_MATRIX,r_MATRIX_T,MCRD,MCRD)
C      
         ele_Length_dx =  COORDS(1,2) - COORDS(1, 1) 
         ele_Length_dy =  COORDS(2,2) - COORDS(2, 1) 
         ele_Length_dz =  COORDS(3,2) - COORDS(3, 1)
C         
         ele_Length = SQRT(ele_Length_dx**2.0D0 + ele_Length_dy**2.0D0 +
     +                     ele_Length_dz**2.0D0     )
C
         ele_peri = 135.02d0/1.0d3
C 
         ele_Length = ele_Length * ele_peri / 2.0D0
C
C     -----------------------------------------------------------------------------------------C
C        COMPUTE THE STIFFNESS MATRIX
C        LOCAL STIFFNESS = B_T * TRAC_JACOB * B
C
         CALL K_MATRIX_MULTIPLY(TRAC_JACOB,r_MATRIX,AMATRX_TEMP1,
     +                          MCRD,MCRD,MCRD)
         CALL K_MATRIX_MULTIPLY(AMATRX_TEMP1,B_MATRIX,AMATRX_TEMP2,
     +                          MCRD,MCRD,NDOFEL)
C     -----------------------------------------------------------------------------------------C
C        COMPUTE GLOBAL STIFFNESS MATRIX 
C        GLOBAL_K = T' * K * T
C
         CALL K_MATRIX_MULTIPLY(r_MATRIX_T,AMATRX_TEMP2,AMATRX_TEMP3,
     &                          MCRD,MCRD,NDOFEL)
         CALL K_MATRIX_MULTIPLY(B_MATRIX_T,AMATRX_TEMP3,AMATRX_TEMP4,
     &                          NDOFEL,MCRD,NDOFEL)          
C     -----------------------------------------------------------------------------------------C
C        MULTIPLY JACOBIAN WITH THE GLOBAL STIFFNESS AND ADD CONTRIBUTION
C        FROM EACH GAUSS POINT
C
         CALL K_MATRIX_PLUS_SCALAR(AMATRX,AMATRX_TEMP4,ele_Length,
     &                             NDOFEL,NDOFEL)          
C     -----------------------------------------------------------------------------------------C
C        COMPUTE THE GLOBAL RESIDUAL VECTOR
C        LOCAL_RESIDUAL = B_T * TRAC
C        GLOBAL_RESIDUAL = T' * LOCAL_RESIDUAL
C 
         TRAC = - TRAC
         CALL K_MATRIX_MULTIPLY(r_MATRIX_T,TRAC,RHS_TEMP1,MCRD,
     &                          MCRD,NRHS)
         CALL K_MATRIX_MULTIPLY(B_MATRIX_T,RHS_TEMP1,
     &                          RHS_TEMP2,NDOFEL,MCRD,NRHS)
C
C     -----------------------------------------------------------------------------------------C
C        MULTIPLY THE GLOBAL RESIDUAL BY THE JACOBIAN AND ADD THE 
C        CONTRIBUTION FROM EACH POINT
C   
         CALL K_MATRIX_PLUS_SCALAR(RHS,RHS_TEMP2,ele_Length,
     &                             NDOFEL,NRHS)
                  
C     -----------------------------------------------------------------------------------------C
      END DO
      
      
C      if(TIME(2) .eq. 14.4d0) then
C          write(6,*) 'time(2) = ', time(2)
C          write(6,*) 'jelem = ', jelem
C          write(6,*) 'DeltaU_Loc_GP(1) = ', DeltaU_Loc_GP(1)
C          write(6,*) 'STRESS(1) = ', STRESS(1)
          
C      end if
      
     
      RETURN
      END
C
      
      
      

      
      
C-----------------------------------------------------------------------------------------C
C-----------------------------------------------------------------------------------------C
C-----------------------------------------------------------------------------------------C
C = MATRIX OPERATIONS =================================================
      SUBROUTINE K_VECTOR_ZERO(A,N)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION A(N)
C
      DO I = 1, N
         A(I)=0.D0
      END DO
C
      RETURN
      END
C-----------------------------------------------------------------------------------------C
       SUBROUTINE K_MATRIX_ZERO (A,N,M)
       INCLUDE 'ABA_PARAM.INC'
       DIMENSION A(N,M)
       DO I = 1, N
          DO J = 1, M
             A(I,J) = 0.0
          END DO
       END DO
       RETURN
       END
C-----------------------------------------------------------------------------------------C
       SUBROUTINE K_MATRIX_TRANSPOSE (A,B,N,M)
       INCLUDE 'ABA_PARAM.INC'
       DIMENSION A(N,M), B(M,N)
       CALL K_MATRIX_ZERO (B,M,N)
       DO I = 1, N
          DO J = 1, M
             B(J,I) = A(I,J)
          END DO
       END DO
       RETURN
      END
C-----------------------------------------------------------------------------------------C
       SUBROUTINE K_MATRIX_MULTIPLY (A,B,C,L,N,M)
       INCLUDE 'ABA_PARAM.INC'
       DIMENSION A(L,N), B(N,M), C(L,M)
       CALL K_MATRIX_ZERO (C,L,M)
       DO I = 1, L
          DO J = 1, M
             DO K = 1, N
               C(I,J) = C(I,J) + A(I,K) * B (K,J)
             END DO
          END DO
       END DO
       RETURN
      END
C-----------------------------------------------------------------------------------------C
      SUBROUTINE K_BASIS(CO_ORI,r_Matrix,U,U_Loc,MCRD,NNODE,NDOFEL)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION CO_DE_M(MCRD,NNODE)
      DIMENSION CO_ORI(MCRD,NNODE), CO_CURR(MCRD,NNODE)
      DIMENSION U(NDOFEL), U_Loc(NDOFEL)
C
      DIMENSION v1(MCRD), v2(MCRD), v3(MCRD), temp(MCRD)
      DIMENSION r_Matrix(MCRD,MCRD),r_Matrix_Tran(MCRD,MCRD)
C     
      DIMENSION U1_Global(3),U2_Global(3),U3_Global(3),U4_Global(3)
      DIMENSION U1_Local(3),U2_Local(3),U3_Local(3),U4_Local(3)
C      
C     -----------------------------------------------------------------------
      CALL K_MATRIX_ZERO(CO_DE_M,3,4)
C
      DO I = 1, 3
         CO_DE_M(I,1) = CO_ORI(I,1)
         CO_DE_M(I,2) = CO_ORI(I,2)
         CO_DE_M(I,3) = CO_ORI(I,3)
         CO_DE_M(I,4) = CO_ORI(I,4)
      END DO
C     -----------------------------------------------------------------------      
      DO I = 1, 3
          v1(I) = 0.5D0 *(   CO_DE_M(I,3) + CO_DE_M(I,2) - 
     +                       CO_DE_M(I,4) - CO_DE_M(I,1) )
          v2(I) = 0.5D0 *(   CO_DE_M(I,4) + CO_DE_M(I,3) - 
     +                       CO_DE_M(I,2) - CO_DE_M(I,1) )
      END DO
      
      v2(1) = v2(1) + 1.0D0
      v2(2) = v2(2) + 0.0D0
      v2(3) = v2(3) + 0.0D0
C
      v1_length = SQRT(v1(1)**2.0D0 + v1(2)**2.0D0 + v1(3)**2.0D0)
      v1_length = MAX(1.0D-6, v1_length)
C
      DO I = 1, 3
          v1(I) = v1(I) / v1_length
      END DO 
C
      alpha = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)
      temp = v1    
      DO I = 1, 3
          temp(I) = temp(I) * alpha
      END DO 
      DO I = 1, 3
          v2(I) = v2(I) - temp(I)
      END DO
C
      v2_length = SQRT(v2(1)**2.0D0 + v2(2)**2.0D0 + v2(3)**2.0D0)
      v2_length = MAX(1.0D-6, v2_length)
      DO I = 1, 3
          v2(I) = v2(I) / v2_length
      END DO       
C
      V3(1) = v1(2)*v2(3)-v1(3)*v2(2)
      V3(2) = v1(3)*v2(1)-v1(1)*v2(3)
      V3(3) = v1(1)*v2(2)-v1(2)*v2(1)
C     -----------------------------------------------------------------------
      CALL K_MATRIX_ZERO(r_Matrix_Tran,MCRD,MCRD)
      CALL K_MATRIX_ZERO(r_Matrix     ,MCRD,MCRD)
C      
      r_Matrix_Tran(1,1) = v1(1)
      r_Matrix_Tran(1,2) = v2(1)
      r_Matrix_Tran(1,3) = v3(1)
      r_Matrix_Tran(2,1) = v1(2)
      r_Matrix_Tran(2,2) = v2(2)
      r_Matrix_Tran(2,3) = v3(2)
      r_Matrix_Tran(3,1) = v1(3)
      r_Matrix_Tran(3,2) = v2(3)
      r_Matrix_Tran(3,3) = v3(3)     
C
      r_Matrix(1,1) = r_Matrix_Tran(1,1)
      r_Matrix(1,2) = r_Matrix_Tran(2,1)
      r_Matrix(1,3) = r_Matrix_Tran(3,1)
      r_Matrix(2,1) = r_Matrix_Tran(1,2)
      r_Matrix(2,2) = r_Matrix_Tran(2,2)
      r_Matrix(2,3) = r_Matrix_Tran(3,2)
      r_Matrix(3,1) = r_Matrix_Tran(1,3)
      r_Matrix(3,2) = r_Matrix_Tran(2,3)
      r_Matrix(3,3) = r_Matrix_Tran(3,3) 
C     -----------------------------------------------------------------------      
      U1_Global(1) = U(1)
      U1_Global(2) = U(2)
      U1_Global(3) = U(3)
      U2_Global(1) = U(4)
      U2_Global(2) = U(5)
      U2_Global(3) = U(6)
      U3_Global(1) = U(7)
      U3_Global(2) = U(8)
      U3_Global(3) = U(9)
      U4_Global(1) = U(10)
      U4_Global(2) = U(11)
      U4_Global(3) = U(12)
C
      U1_Local = matmul(r_Matrix , U1_Global)
      U2_Local = matmul(r_Matrix , U2_Global)
      U3_Local = matmul(r_Matrix , U3_Global)
      U4_Local = matmul(r_Matrix , U4_Global)
C
      U_Loc(1)  = U1_Local(1)
      U_Loc(2)  = U1_Local(2)
      U_Loc(3)  = U1_Local(3)
      U_Loc(4)  = U2_Local(1)
      U_Loc(5)  = U2_Local(2)
      U_Loc(6)  = U2_Local(3)
      U_Loc(7)  = U3_Local(1)
      U_Loc(8)  = U3_Local(2)
      U_Loc(9)  = U3_Local(3)
      U_Loc(10) = U4_Local(1)
      U_Loc(11) = U4_Local(2)
      U_Loc(12) = U4_Local(3)
C     -----------------------------------------------------------------------            
      RETURN
      END
C-----------------------------------------------------------------------------------------C
      SUBROUTINE K_SHAPE_FUN(I,SF)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION GAUSS2(2), SF(2)
C
      GAUSS2(1) = -sqrt(1.0D0/3.0D0)
      GAUSS2(2) =  sqrt(1.0D0/3.0D0)
C
      SF(1) = 0.5D0 * ( 1 - GAUSS2(I))
      SF(2) = 0.5D0 * ( 1 + GAUSS2(I)) 
C
      RETURN
      END
C-----------------------------------------------------------------------------------------C
      SUBROUTINE K_MATRIX_PLUS_SCALAR(A,B,C,N,M)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION A(N,M),B(N,M)
C
      DO I = 1, N
         DO J = 1, M
            A(I,J)=A(I,J)+C*B(I,J)
         END DO
      END DO
C
      RETURN
      END      
C-----------------------------------------------------------------------------------------C


      
      
      
      
      
C################################################################################# C
C																				 C	
C			  SUBROUTINE DMGEV TO CALCULATE UNIAXIAL DAMAGE	EVOLUTION			 C
C						   														 C
C################################################################################# C  
C
       SUBROUTINE UMAT_BOND(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1    RPL,DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,
     2    TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,MATERL,NDI,NSHR,NTENS,
     3    NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,CELENT,
     4    DFGRD0,DFGRD1,NOEL,NPT,KSLAY,KSPT,KSTEP,KINC)
C -------------------------------------------------------------------------------  C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),
     2 DDSDDT(NTENS),DRPLDE(NTENS), 
     3 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     4 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
C 
      DIMENSION envlpPosStress(6), envlpPosStrain(6)
      DIMENSION envlpNegStress(6), envlpNegStrain(6)
      DIMENSION envlpPosDamgdStress(6), envlpNegDamgdStress(6)
      DIMENSION state3Stress(4), state3Strain(4)
      DIMENSION state4Stress(4), state4Strain(4)
C      
      INTEGER  Tstate, Cstate, DmgCyc
C      
      PARAMETER (ONE=1.0,TWO=2.0,THREE=3.0,SIX=6.0, HALF =0.5,
     1 N_ELEM=2,NSTV=18) 
C -------------------------------------------------------------------------------  C
C     Data input & parameter initialization 
C
      call UMAT_INPUT(stress1p, stress2p, stress3p, stress4p,
     +    strain1p, strain2p, strain3p, strain4p, stress1n, stress2n,
     +    stress3n, stress4n, strain1n, strain2n, strain3n, strain4n, 
     +    rDispP, rForceP, uForceP,rDispN, rForceN, uForceN,
     +    gammaK1, gammaK2, gammaK3, gammaK4, gammaKLimit,
     +    gammaD1, gammaD2, gammaD3, gammaD4, gammaDLimit,
     +    gammaF1, gammaF2, gammaF3, gammaF4, gammaFLimit,
     +    gammaE,  DmgCyc)
C      
      energyCapacity = 0.0D0 
      sfUnload       = 0.0D0      ! replace kunload
      elasticStrainEnergy = 0.0D0
C
      gammaKUsed = 0.0D0
      gammaFUsed = 0.0D0
C      
      call K_VECTOR_ZERO(envlpPosStress, 6)
      call K_VECTOR_ZERO(envlpPosStrain, 6)
      call K_VECTOR_ZERO(envlpNegStress, 6)
      call K_VECTOR_ZERO(envlpNegStrain, 6)
C
      call K_VECTOR_ZERO(state3Stress,   4)
      call K_VECTOR_ZERO(state3Strain,   4)
      call K_VECTOR_ZERO(state4Stress,   4)
      call K_VECTOR_ZERO(state4Strain,   4)
C
      call K_VECTOR_ZERO(envlpPosDamgdStress, 6)
      call K_VECTOR_ZERO(envlpNegDamgdStress, 6)          
C      
      call SetEnvelope(stress1p, stress2p, stress3p, stress4p,
     1    strain1p, strain2p, strain3p, strain4p, stress1n, stress2n,
     2    stress3n, stress4n, strain1n, strain2n, strain3n, strain4n, 
     3    gammaE, envlpPosStress,envlpPosStrain, envlpNegStress,
     4    envlpNegStrain, ElasticPos, ElasticNeg, energyCapacity)
C      
      envlpPosDamgdStress = envlpPosStress 
      envlpNegDamgdStress = envlpNegStress 
C -------------------------------------------------------------------------------  C
      if (  ((TIME(2)-DTIME) .eq. 0.0d0) ) then 
          call revertToStart(envlpPosStress, envlpPosStrain, 
     +        envlpNegStress, envlpNegStrain, ElasticPosDamgd, 
     +        ElasticNegDamgd, gammaKUsed, gammaFUsed,
     +        ElasticPos, ElasticNeg,
     +        Cstate, Cstrain, Cstress, CstrainRate, CstateStrainLow, 
     +        CstateStressLow, CstateStrainHgh, CstateStressHgh,
     +        CminStrainDmnd, CmaxStrainDmnd, Cenergy, CgammaK, CgammaD,
     +        CgammaF, CnCycle,Ttangent, uMaxDamgd, uMinDamgd)
      else
          call revertToLastCommit(Cstate, Cstrain, Cstress, 
     +        CstateStrainLow,CstateStressLow,CstateStrainHgh,
     +        CstateStressHgh, CminStrainDmnd, CmaxStrainDmnd, Cenergy, 
     +        CgammaK, CgammaD, CgammaF,CnCycle, CstrainRate,
     +        envlpPosDamgdStress, envlpNegDamgdStress,
     +        uMaxDamgd,uMinDamgd, ElasticPosDamgd,ElasticNegDamgd,
     +        STATEV,NSTATV)
      end if
C -------------------------------------------------------------------------------  C
C     begin SetTrialStrain
C
      Tstrain = STRAN(1)
      Dstrain = DSTRAN(1)
C
      Tstate = Cstate
C
      TstrainRate     = CstrainRate
C      
      TstateStrainLow = CstateStrainLow
      TstateStrainHgh = CstateStrainHgh
      TstateStressLow = CstateStressLow
      TstateStressHgh = CstateStressHgh
      TminStrainDmnd  = CminStrainDmnd
      TmaxStrainDmnd  = CmaxStrainDmnd
      Tenergy         = Cenergy
C      
      TgammaF         = CgammaF
      TgammaK         = CgammaK 
      TgammaD         = CgammaD   
C
      TnCycle         = CnCycle
C      
      gammaKUsed = STATEV(32) 
      TgammaD    = STATEV(33) 
      gammaFUsed = STATEV(34) 
C -------------------------------------------------------------------------------  C
C     Get the state of loading
C
      call getstate(Tstrain, Dstrain, CstrainRate, TstateStrainLow,
     1    TstateStressLow, TstateStrainHgh, TstateStressHgh, Tstate,
     2    envlpPosStress,envlpPosStrain, envlpNegStress, envlpNegStrain,
     3    envlpPosDamgdStress, envlpNegDamgdStress,
     4    TmaxStrainDmnd,  TminStrainDmnd, ElasticPosDamgd, 
     5    ElasticNegDamgd, gammaKUsed, gammaFUsed,
     6    uMaxDamgd, uMinDamgd, CgammaF, CgammaK, Cstress, Cstrain,
     7    ElasticPos, ElasticNeg )
C -------------------------------------------------------------------------------  C
      select case (Tstate)
C -------------------------------------------------------------------------------  C
      case (0)
          Ttangent = envlpPosStress(0+1)/envlpPosStrain(0+1)
          Tstress  = Ttangent * Tstrain
C -------------------------------------------------------------------------------  C
      case (1)
          call posEnvlpStress (Tstrain,envlpPosDamgdStress, 
     +                         envlpPosStrain,Tstress)
          call posEnvlpTangent(Tstrain,envlpPosDamgdStress, 
     +                         envlpPosStrain,Ttangent)
C -------------------------------------------------------------------------------  C
      case (2)
          call negEnvlpStress (Tstrain,envlpNegDamgdStress, 
     +                         envlpNegStrain,Tstress)
          call negEnvlpTangent(Tstrain,envlpNegDamgdStress, 
     +                         envlpNegStrain,Ttangent)
C -------------------------------------------------------------------------------  C
      case (3)
          if (  TstateStrainHgh .lt. 0.0d0  ) then
              sfUnload = ElasticNegDamgd
          else
              sfUnload = ElasticPosDamgd
          end if
		state3Strain(0+1) = TstateStrainLow
		state3Strain(3+1) = TstateStrainHgh
		state3Stress(0+1) = TstateStressLow
		state3Stress(3+1) = TstateStressHgh
C          
          call getstate3(sfUnload, envlpNegDamgdStress, envlpNegStrain,
     +    TstateStrainLow, TstateStressLow, TstateStrainHgh, 
     +    TstateStressHgh, TminStrainDmnd, ElasticNegDamgd,
     +    state3Strain, state3Stress, rDispN, uForceN, rForceN)
C
          call Envlp3Stress (Tstrain,state3Strain, state3Stress,Tstress)
          call Envlp3Tangent(Tstrain,state3Strain,state3Stress,Ttangent)
C -------------------------------------------------------------------------------  C
      case (4)
          if (  TstateStrainLow .lt. 0.0d0  ) then
              sfUnload = ElasticNegDamgd
          else
              sfUnload = ElasticPosDamgd
          end if
		state4Strain(0+1) = TstateStrainLow
		state4Strain(3+1) = TstateStrainHgh
		state4Stress(0+1) = TstateStressLow
		state4Stress(3+1) = TstateStressHgh  
C          
          call getstate4(sfUnload, envlpPosDamgdStress, envlpPosStrain,
     +    TstateStrainLow, TstateStressLow, TstateStrainHgh, 
     +    TstateStressHgh, TmaxStrainDmnd, ElasticPosDamgd,
     +    state4Strain, state4Stress, rDispP, uForceP, rForceP)
C          
          call Envlp4Stress (Tstrain,state4Strain, state4Stress,Tstress)
          call Envlp4Tangent(Tstrain,state4Strain,state4Stress,Ttangent)
C -------------------------------------------------------------------------------  C
      end select
C -------------------------------------------------------------------------------  C
C      
      denergy = 0.5 * (Tstress+Cstress) * Dstrain
C
      if (  Tstrain .gt. 0.0  ) then 
          elasticStrainEnergy = 0.5 * Tstress/ElasticPosDamgd*Tstress
      else
          elasticStrainEnergy = 0.5 * Tstress/ElasticNegDamgd*Tstress
      end if
C
      Tenergy = Cenergy + denergy
     
      call updateDmg(Tstrain, Dstrain, CnCycle, energyCapacity,  
     +       Tenergy, ElasticPos, ElasticNeg, elasticStrainEnergy,
     +       gammaK1, gammaK2, gammaK3, gammaK4, gammaKLimit,
     +       gammaD1, gammaD2, gammaD3, gammaD4, gammaDLimit,
     +       gammaF1, gammaF2, gammaF3, gammaF4, gammaFLimit,
     +       TmaxStrainDmnd,  TminStrainDmnd, DmgCyc,
     +       envlpPosDamgdStress, envlpPosStrain, envlpNegDamgdStress, 
     +       envlpNegStrain, TnCycle, TgammaK, TgammaD, TgammaF)
C -------------------------------------------------------------------------------  C
C     GET STIFFNESS MATRIX AND STRESS
      DDSDDE(1,1) = Ttangent  
      DDSDDE(2,2) = 200.0D9   ! DO NOT CHANGE 
      DDSDDE(3,3) = 200.0D9   ! DO NOT CHANGE 
C
      STRESS(1) = Tstress 
	STRESS(2) = DDSDDE(2,2) * STRAN(2)
      STRESS(3) = DDSDDE(3,3) * STRAN(3)
C
C -------------------------------------------------------------------------------  C
C     STATE VARIABLES STORE
      call commitState(envlpPosStress, envlpNegStress,
     +  envlpPosDamgdStress, envlpNegDamgdStress,  
     +  TstateStrainLow,TstateStressLow,TstateStrainHgh,TstateStressHgh,
     +  TmaxStrainDmnd,TminStrainDmnd,ElasticPosDamgd,ElasticNegDamgd,
     +  gammaKUsed, gammaFUsed, Tstate,
     +  TstrainRate,Tenergy,Tstress,Tstrain,TgammaK,TgammaD,TgammaF,
     +  ElasticPos, ElasticNeg, TnCycle, 
     +  Cstate, CstrainRate, CstateStrainLow, CstateStressLow,
     +  CstateStrainHgh, CstateStressHgh,CminStrainDmnd, CmaxStrainDmnd,
     +  Cenergy, Cstress, Cstrain, CgammaK,  CgammaD, CgammaF,
     +  uMaxDamgd, uMinDamgd, CnCycle,
     +  DSTRAIN, STATEV,NSTATV)
C      
      STATEV(32) = gammaKUsed
      STATEV(33) = CgammaD
      STATEV(34) = gammaFUsed
      STATEV(35) = Ttangent
C -------------------------------------------------------------------------------  C
      RETURN
      END    
      
      
      
      
C################################################################################# C
C																				 C	
C			            SUBROUTINE SetEnvelope  (Checked)                        C
C						   														 C
C################################################################################# C
      SUBROUTINE SetEnvelope(stress1p, stress2p, stress3p, stress4p,
     1    strain1p, strain2p, strain3p, strain4p, stress1n, stress2n,
     2    stress3n, stress4n, strain1n, strain2n, strain3n, strain4n, 
     3    gammaE, envlpPosStress,envlpPosStrain, envlpNegStress,
     4    envlpNegStrain, ElasticPos, ElasticNeg, energyCapacity)
C      
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION envlpPosStress(6), envlpPosStrain(6)
      DIMENSION envlpNegStress(6), envlpNegStrain(6)
C -------------------------------------------------------------------------------  C
      s_kPos = stress1p / strain1p
      s_kNeg = stress1n / strain1n
C
      if ( s_kPos .gt. s_kNeg ) then
          s_k = s_kPos
      else
          s_k = s_kNeg
      end if
C
      if ( strain1p  .gt.  (- strain1n) ) then
          u =  1.0D-4 * strain1p
      else
          u = -1.0D-4 * strain1n
      end if
C
C     form the key points in stress-strain curve
C
      envlpPosStrain(0+1) = u
      envlpPosStress(0+1) = u*s_k
      envlpNegStrain(0+1) = -u
      envlpNegStress(0+1) = -u*s_k 
C      
      envlpPosStrain(1+1) = strain1p
      envlpPosStrain(2+1) = strain2p
      envlpPosStrain(3+1) = strain3p
      envlpPosStrain(4+1) = strain4p
C      
      envlpNegStrain(1+1) = strain1n
      envlpNegStrain(2+1) = strain2n
      envlpNegStrain(3+1) = strain3n
      envlpNegStrain(4+1) = strain4n
C      
      envlpPosStress(1+1) = stress1p
      envlpPosStress(2+1) = stress2p
      envlpPosStress(3+1) = stress3p
      envlpPosStress(4+1) = stress4p
C      
      envlpNegStress(1+1) = stress1n
      envlpNegStress(2+1) = stress2n
      envlpNegStress(3+1) = stress3n
      envlpNegStress(4+1) = stress4n
C -------------------------------------------------------------------------------  C
      s_k1 = (stress4p - stress3p)/(strain4p - strain3p)
      s_k2 = (stress4n - stress3n)/(strain4n - strain3n)
C
      envlpPosStrain(5+1) = 1.0D6 * strain4p
      if (s_k1 .gt. 0.0D0) then
          envlpPosStress(5+1) = stress4p + 
     +                          s_k1 * (envlpPosStrain(5+1) - strain4p)
      else
          envlpPosStress(5+1) = stress4p * 1.1D0
      end if
C
      envlpNegStrain(5+1) = 1.0D6 * strain4n
      if (s_k2 .gt. 0.0)  then
          envlpNegStress(5+1) = stress4n + 
     +                          s_k2 * (envlpNegStrain(5+1) - strain4n)
      else
          envlpNegStress(5+1) = stress4n * 1.1D0
      end if
C
C     form the key points in stress-strain curve
C     define crtical material properties
C
      ElasticPos = envlpPosStress(1+1)/envlpPosStrain(1+1)
      ElasticNeg = envlpNegStress(1+1)/envlpNegStrain(1+1)
C -------------------------------------------------------------------------------  C
      energypos = 0.5D0 * envlpPosStrain(0+1) * envlpPosStress(0+1)
      do jt = 0,3
          energypos = energypos + 0.5*(envlpPosStress(jt+1) + 
     +                envlpPosStress(jt+1+1))*(envlpPosStrain(jt+1+1) -
     +                envlpPosStrain(jt+1))
      end do
C
      energyneg = 0.5*envlpNegStrain(0+1)*envlpNegStress(0+1)
      do jy = 0,3
          energyneg = energyneg + 0.5*(envlpNegStress(jy+1) + 
     +                envlpNegStress(jy+1+1))*(envlpNegStrain(jy+1+1) -
     +                envlpNegStrain(jy+1))
      end do
C
      if (energypos .gt. energyneg)  then
          energy_max = energypos
      else
          energy_max = energyneg
      end if      
C
      energyCapacity = gammaE * energy_max
C -------------------------------------------------------------------------------  C
      RETURN
      END
C################################################################################# C
C																				 C	
C			            SUBROUTINE commitState(checked)			                 C
C						   														 C
C################################################################################# C
      SUBROUTINE commitState(envlpPosStress, envlpNegStress,
     +  envlpPosDamgdStress, envlpNegDamgdStress,  
     +  TstateStrainLow,TstateStressLow,TstateStrainHgh,TstateStressHgh,
     +  TmaxStrainDmnd,TminStrainDmnd,ElasticPosDamgd,ElasticNegDamgd,
     +  gammaKUsed, gammaFUsed, Tstate,
     +  TstrainRate,Tenergy,Tstress,Tstrain,TgammaK,TgammaD,TgammaF,
     +  ElasticPos, ElasticNeg, TnCycle, 
     +  Cstate, CstrainRate, CstateStrainLow, CstateStressLow,
     +  CstateStrainHgh, CstateStressHgh,CminStrainDmnd, CmaxStrainDmnd,
     +  Cenergy, Cstress, Cstrain, CgammaK,  CgammaD, CgammaF,
     +  uMaxDamgd, uMinDamgd, CnCycle,
     +  du, STATEV,NSTATV)
C      
      INCLUDE 'ABA_PARAM.INC'
C      
      DIMENSION STATEV(NSTATV)
      DIMENSION envlpPosStress(6), envlpPosStrain(6)
      DIMENSION envlpNegStress(6), envlpNegStrain(6)
      DIMENSION envlpPosDamgdStress(6), envlpNegDamgdStress(6)
C     
      integer Tstate,Cstate
C     --------------------------------------------------------------------
C      
      Cstate = Tstate
C
      if (  (du .gt. 1.0D-12) .or. (du .lt. -1.0D-12)  ) then
	    CstrainRate = du
      else 
	    CstrainRate = TstrainRate
      end if
C
      CstateStrainLow = TstateStrainLow
      CstateStressLow = TstateStressLow
      CstateStrainHgh = TstateStrainHgh
      CstateStressHgh = TstateStressHgh
      CminStrainDmnd  = TminStrainDmnd
      CmaxStrainDmnd  = TmaxStrainDmnd
      Cenergy         = Tenergy
C
      Cstress         = Tstress
      Cstrain         = Tstrain
C
      CgammaK         = TgammaK
      CgammaD         = TgammaD
      CgammaF         = TgammaF
C
C     define adjusted strength and stiffness parameters
C
      ElasticPosDamgd    = ElasticPos*(1.0D0 - gammaKUsed)
      ElasticNegDamgd    = ElasticNeg*(1.0D0 - gammaKUsed)
C
      uMaxDamgd          = TmaxStrainDmnd            !*(1.0D0 + CgammaD)  
      uMinDamgd          = TminStrainDmnd            !*(1.0D0 + CgammaD)
C
      envlpPosDamgdStress = envlpPosStress *(1.0D0 - gammaFUsed)
      envlpNegDamgdStress = envlpNegStress *(1.0D0 - gammaFUsed)
C
      CnCycle = TnCycle                               ! number of cycles of loading
C
C     --------------------------------------------------------------------
      STATEV(1)  = Cstate
      STATEV(2)  = CstateStrainLow
      STATEV(3)  = CstateStressLow
      STATEV(4)  = CstateStrainHgh
      STATEV(5)  = CstateStressHgh
      STATEV(6)  = CminStrainDmnd
      STATEV(7)  = CmaxStrainDmnd
      STATEV(8)  = Cenergy
      STATEV(9)  = Cstress
      STATEV(10) = Cstrain
      STATEV(11) = CgammaK
      STATEV(12) = CgammaD
      STATEV(13) = CgammaF
      STATEV(14) = ElasticPosDamgd
      STATEV(15) = ElasticNegDamgd
      STATEV(16) = uMaxDamgd
      STATEV(17) = uMinDamgd
      STATEV(18) = envlpPosDamgdStress(1)
      STATEV(19) = envlpPosDamgdStress(2)
      STATEV(20) = envlpPosDamgdStress(3)
      STATEV(21) = envlpPosDamgdStress(4)
      STATEV(22) = envlpPosDamgdStress(5)
      STATEV(23) = envlpPosDamgdStress(6)
      STATEV(24) = envlpNegDamgdStress(1)
      STATEV(25) = envlpNegDamgdStress(2)
      STATEV(26) = envlpNegDamgdStress(3)
      STATEV(27) = envlpNegDamgdStress(4)
      STATEV(28) = envlpNegDamgdStress(5)
      STATEV(29) = envlpNegDamgdStress(6)
      STATEV(30) = CnCycle
      STATEV(31) = CstrainRate
C      
C     --------------------------------------------------------------------
      RETURN
      END
C################################################################################# C
C																				 C	
C			            SUBROUTINE UMAT_INPUT (checked)			                 C
C						   														 C
C################################################################################# C
      SUBROUTINE UMAT_INPUT(stress1p, stress2p, stress3p, stress4p,
     +    strain1p, strain2p, strain3p, strain4p, stress1n, stress2n,
     +    stress3n, stress4n, strain1n, strain2n, strain3n, strain4n, 
     +    rDispP, rForceP, uForceP,rDispN, rForceN, uForceN,
     +    gammaK1, gammaK2, gammaK3, gammaK4, gammaKLimit,
     +    gammaD1, gammaD2, gammaD3, gammaD4, gammaDLimit,
     +    gammaF1, gammaF2, gammaF3, gammaF4, gammaFLimit,
     +    gammaE,  DmgCyc)
      
      INCLUDE 'ABA_PARAM.INC'
      INTEGER DmgCyc
C -------------------------------------------------------------------------------  C
      stress_coeff = 3.75D6     ! 19.6 MPa   2.5*fc^0.5/7
      strain_eff   = 0.20D0
      
      stress1p = 2.0D0 * stress_coeff
      strain1p = 0.0001D0 * strain_eff
      stress2p = 6.0D0 * stress_coeff
      strain2p = 0.0055D0 * strain_eff
      stress3p = 7.0D0 * stress_coeff
      strain3p = 0.0188D0 * strain_eff
      stress4p = 2.8D0 * stress_coeff
      strain4p = 0.08D0 * strain_eff
C   
      stress1n = -2.0D0 * stress_coeff
      strain1n = -0.0001D0 * strain_eff
      stress2n = -6.0D0 * stress_coeff
      strain2n = -0.0055D0 * strain_eff
      stress3n = -7.0D0 * stress_coeff
      strain3n = -0.0188D0 * strain_eff
      stress4n = -2.8D0 * stress_coeff
      strain4n = - 0.08D0 * strain_eff
C
      rDispP  = 0.5D0
      rForceP = 0.25D0
      uForceP = 0.05D0 
      rDispN  = 0.5D0
      rForceN = 0.25D0
      uForceN = 0.05D0 
C
      gammaK1 = 1.0D0
      gammaD1 = 0.5D0
      gammaF1 = 1.0D0
      gammaK2 = 0.2D0
      gammaD2 = 0.5D0
      gammaF2 = 0.0D0
      gammaK3 = 0.3D0
      gammaD3 = 2.0D0
      gammaF3 = 1.0D0
      gammaK4 = 0.2D0
      gammaD4 = 2.0D0
      gammaF4 = 1.0D0
      gammaKLimit = 0.9D0
      gammaDLimit = 0.5D0
      gammaFLimit = 0.9D0
C
      gammaE = 10.0D0
C
      DmgCyc = 0
C -------------------------------------------------------------------------------  C
      RETURN
      END
C################################################################################# C
C																				 C	
C			        SUBROUTINE posEnvlpStress (checked)				             C
C						   														 C
C################################################################################# C
      SUBROUTINE posEnvlpStress(u, envlpPosDamgdStress, 
     +                          envlpPosStrain, s_f)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION  envlpPosDamgdStress(6), envlpPosStrain(6)
      INTEGER i
C -------------------------------------------------------------------------------  C
      s_k = 0.0D0     !  stiffness 
      i = 0
      s_f = 0.0D0     !  stress
C
      do while ( (s_k .eq. 0.0D0) .and. (i .le. 4) )   
          if (u .le. envlpPosStrain(i+1+1)) then
            s_k = (envlpPosDamgdStress(i+1+1)-envlpPosDamgdStress(i+1))/
     +            (envlpPosStrain(i+1+1)-envlpPosStrain(i+1))
            s_f = envlpPosDamgdStress(i+1) + (u-envlpPosStrain(i+1))*s_k
          end if
          i = i+1
      end do
C      
      if (s_k .eq. 0.0D0) then
          s_k = (envlpPosDamgdStress(5+1) - envlpPosDamgdStress(4+1)) / 
     +          (envlpPosStrain(5+1) - envlpPosStrain(4+1))
          s_f =  envlpPosDamgdStress(5+1) + s_k*(u-envlpPosStrain(5+1))
      end if
C -------------------------------------------------------------------------------  C
      RETURN
      END
C################################################################################# C
C																				 C	
C			        SUBROUTINE posEnvlpTangent (checked)				         C
C						   														 C
C################################################################################# C
      SUBROUTINE posEnvlpTangent(u,envlpPosDamgdStress, 
     +                           envlpPosStrain, s_k)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION  envlpPosDamgdStress(6), envlpPosStrain(6)
      INTEGER i
C -------------------------------------------------------------------------------  C
      s_k = 0.0D0
      i = 0
C
      do while ( (s_k .eq. 0.0D0) .and. (i .le. 4) )   
          if (u .le. envlpPosStrain(i+1+1)) then
            s_k = (envlpPosDamgdStress(i+1+1)-envlpPosDamgdStress(i+1))/
     +            (envlpPosStrain(i+1+1)-envlpPosStrain(i+1))
          end if 
          i = i+1
      end do
C      
      if (s_k .eq. 0.0D0) then
          s_k = (envlpPosDamgdStress(5+1) - envlpPosDamgdStress(4+1)) / 
     +        (envlpPosStrain(5+1) - envlpPosStrain(4+1))
      end if
C -------------------------------------------------------------------------------  C
      RETURN
      END     
C################################################################################# C
C																				 C	
C			        SUBROUTINE negEnvlpStress (checked)				             C
C						   														 C
C################################################################################# C
      SUBROUTINE negEnvlpStress(u,envlpNegDamgdStress, 
     +                          envlpNegStrain, s_f)
      
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION  envlpNegDamgdStress(6), envlpNegStrain(6)
      INTEGER i
C -------------------------------------------------------------------------------  C
      s_k = 0.0D0
      i = 0
      s_f = 0.0D0
C
      do while ( (s_k .eq. 0.0D0) .and. (i .le. 4) )   
          if (u .ge. envlpNegStrain(i+1+1)) then
            s_k = (envlpNegDamgdStress(i+1)-envlpNegDamgdStress(i+1+1))/
     +            (envlpNegStrain(i+1)-envlpNegStrain(i+1+1))
            s_f =  envlpNegDamgdStress(i+1+1)+ 
     +            (u-envlpNegStrain(i+1+1))*s_k
          end if 
          i = i+1
      end do
C      
      if (s_k .eq. 0.0D0) then
          s_k = (envlpNegDamgdStress(4+1) - envlpNegDamgdStress(5+1)) / 
     +          (envlpNegStrain(4+1) - envlpNegStrain(5+1))
          s_f =  envlpNegDamgdStress(5+1) + s_k *(u-envlpNegStrain(5+1))
      end if      
C -------------------------------------------------------------------------------  C
      RETURN
      END      
C################################################################################# C
C																				 C	
C			        SUBROUTINE negEnvlpTangent (checked)				         C
C						   														 C
C################################################################################# C
      SUBROUTINE negEnvlpTangent(u,envlpNegDamgdStress, 
     +                           envlpNegStrain, s_k)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION  envlpNegDamgdStress(6), envlpNegStrain(6)
      INTEGER i
C -------------------------------------------------------------------------------  C
      s_k = 0.0D0
      i = 0
C
      do while ( (s_k .eq. 0.0D0) .and. (i .le. 4) )   
          if (u .ge. envlpNegStrain(i+1+1)) then
            s_k = (envlpNegDamgdStress(i+1)-envlpNegDamgdStress(i+1+1))/
     +            (envlpNegStrain(i+1)-envlpNegStrain(i+1+1))
          end if 
          i = i+1
      end do
C      
      if (s_k .eq. 0.0D0) then
          s_k = (envlpNegDamgdStress(4+1) - envlpNegDamgdStress(5+1)) / 
     +          (envlpNegStrain(4+1) - envlpNegStrain(5+1))
      end if      
C -------------------------------------------------------------------------------  C
      RETURN
      END
C################################################################################# C
C																				 C	
C			        SUBROUTINE Envlp3Stress (checked)				             C
C						   														 C
C################################################################################# C
      SUBROUTINE Envlp3Stress(u,state3Strain, state3Stress,s_f)
      
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION state3Stress(4), state3Strain(4)
      INTEGER i
C -------------------------------------------------------------------------------  C
      s_k = 0.0
      i = 0
      s_f = 0.0
C
      do while ( ((s_k .eq. 0.0D0).or.(i .le. 2)) .and. (i .le. 2) ) 
         if (u .ge. state3Strain(i+1)) then 
		  s_k = (state3Stress(i+1+1) - state3Stress(i+1)) / 
     +            (state3Strain(i+1+1) - state3Strain(i+1))
            s_f =  state3Stress(i+1)  + (u - state3Strain(i+1)) * s_k
         end if
         i = i+1 
      end do
C
      if (s_k .eq. 0.0) then 
	    if (u .lt. state3Strain(0+1)) then 
		    i = 0
	    else 
		    i = 2
          end if
          s_k = (state3Stress(i+1+1) - state3Stress(i+1)) / 
     +          (state3Strain(i+1+1) - state3Strain(i+1))
          s_f =  state3Stress(i+1) + (u - state3Strain(i+1)) * s_k
      end if
C -------------------------------------------------------------------------------  C
      RETURN
      END
C################################################################################# C
C																				 C	
C			        SUBROUTINE Envlp3Tangent (checked)				             C
C						   														 C
C################################################################################# C
      SUBROUTINE Envlp3Tangent(u,state3Strain, state3Stress,s_k)
      
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION state3Stress(4), state3Strain(4)
      INTEGER i
C -------------------------------------------------------------------------------  C
      s_k = 0.0
      i = 0
C
      do while ( ((s_k .eq. 0.0D0).or.(i .le. 2)) .and. (i .le. 2) ) 
         if (u .ge. state3Strain(i+1)) then 
		  s_k = (state3Stress(i+1+1) - state3Stress(i+1)) / 
     +            (state3Strain(i+1+1) - state3Strain(i+1))
         end if
         i = i+1 
      end do
C
      if (s_k .eq. 0.0) then 
	    if (u .lt. state3Strain(0+1)) then 
		    i = 0
	    else 
		    i = 2
          end if
          s_k = (state3Stress(i+1+1) - state3Stress(i+1)) / 
     +          (state3Strain(i+1+1) - state3Strain(i+1))
      end if     
C -------------------------------------------------------------------------------  C
      RETURN
      END
C################################################################################# C
C																				 C	
C			        SUBROUTINE Envlp4Stress (checked)				             C
C						   														 C
C################################################################################# C
      SUBROUTINE Envlp4Stress(u,state4Strain, state4Stress,s_f)
      
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION state4Stress(4), state4Strain(4)
      INTEGER i
C -------------------------------------------------------------------------------  C
      s_k = 0.0
      i = 0
      s_f = 0.0
C
      do while ( ((s_k .eq. 0.0D0).or.(i .le. 2)) .and. (i .le. 2) ) 
         if (u .ge. state4Strain(i+1)) then 
		  s_k = (state4Stress(i+1+1) - state4Stress(i+1)) / 
     +            (state4Strain(i+1+1) - state4Strain(i+1))
            s_f =  state4Stress(i+1)  + (u-state4Strain(i+1)) * s_k
         end if
         i = i+1 
      end do
C
      if (s_k .eq. 0.0) then 
	    if (u .lt. state4Strain(0+1)) then 
		    i = 0
	    else 
		    i = 2
          end if
          s_k = (state4Stress(i+1+1) - state4Stress(i+1)) / 
     +          (state4Strain(i+1+1) - state4Strain(i+1))
          s_f =  state4Stress(i+1) + (u - state4Strain(i+1)) * s_k
      end if
C -------------------------------------------------------------------------------  C
      RETURN
      END
C################################################################################# C
C																				 C	
C			        SUBROUTINE Envlp4Tangent (checked)				             C
C						   														 C
C################################################################################# C
      SUBROUTINE Envlp4Tangent(u,state4Strain, state4Stress,s_k)
      
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION state4Stress(4), state4Strain(4)
      INTEGER i
C -------------------------------------------------------------------------------  C
      s_k = 0.0
      i = 0
C
      do while ( ((s_k .eq. 0.0D0).or.(i .le. 2)) .and. (i .le. 2) ) 
         if (u .ge. state4Strain(i+1)) then 
		  s_k = (state4Stress(i+1+1) - state4Stress(i+1)) / 
     +            (state4Strain(i+1+1) - state4Strain(i+1))
         end if
         i = i+1 
      end do
C
      if (s_k .eq. 0.0) then 
	    if (u .lt. state4Strain(0+1)) then 
		    i = 0
	    else 
		    i = 2
          end if
          s_k = (state4Stress(i+1+1) - state4Stress(i+1)) / 
     +          (state4Strain(i+1+1) - state4Strain(i+1))
      end if
C -------------------------------------------------------------------------------  C
      RETURN
      END
C################################################################################# C
C																				 C	
C			        SUBROUTINE getstate (checked)		    		             C
C						   														 C
C################################################################################# C
      SUBROUTINE getstate(u, du, CstrainRate, TstateStrainLow,
     1    TstateStressLow, TstateStrainHgh, TstateStressHgh, Tstate,
     2    envlpPosStress,envlpPosStrain, envlpNegStress, envlpNegStrain,
     3    envlpPosDamgdStress, envlpNegDamgdStress,
     4    TmaxStrainDmnd, TminStrainDmnd, ElasticPosDamgd, 
     5    ElasticNegDamgd, gammaKUsed, gammaFUsed,
     6    uMaxDamgd, uMinDamgd, CgammaF, CgammaK, Cstress, Cstrain,
     7    ElasticPos, ElasticNeg )
      
      INCLUDE 'ABA_PARAM.INC'
C      
      DIMENSION envlpPosStress(6), envlpPosStrain(6)
      DIMENSION envlpNegStress(6), envlpNegStrain(6)
      DIMENSION envlpPosDamgdStress(6), envlpNegDamgdStress(6)
C
      INTEGER   newState, Tstate, cid, cis, i
C -------------------------------------------------------------------------------  C
      cid = 0
      cis = 0        !  change in state (CIS)
      newState = 0   !
C
      if ( (du * CstrainRate) .le. 0.0D0) then    ! means the direction of dstrain change.
          cid = 1
      end if
C -------------------------------------------------------------------------------  C
      if ( (u .lt. TstateStrainLow) .or. (u .gt. TstateStrainHgh) .or. 
     +     (cid .eq. 1) )  then
C -------------------------------------------------------------------------------  C
      if (Tstate .eq. 0)  then
          if (u .gt. TstateStrainHgh)  then
              cis = 1  
              newState = 1
			TstateStrainLow = envlpPosStrain(0+1)
			TstateStressLow = envlpPosStress(0+1)
			TstateStrainHgh = envlpPosStrain(5+1)
			TstateStressHgh = envlpPosStress(5+1)
          else if (u .lt. TstateStrainLow) then
              cis = 1
			newState = 2
			TstateStrainLow = envlpNegStrain(5+1)
			TstateStressLow = envlpNegStress(5+1)
			TstateStrainHgh = envlpNegStrain(0+1)
			TstateStressHgh = envlpNegStress(0+1)              
          end if
      end if
C -------------------------------------------------------------------------------  C
      if ( (Tstate .eq. 1) .and. (du .lt. 0.0) )  then
          cis = 1
          if (Cstrain  .gt. TmaxStrainDmnd) then 
			TmaxStrainDmnd = u - du
          end if
          if (TmaxStrainDmnd .lt. uMaxDamgd) then
			TmaxStrainDmnd = uMaxDamgd
          end if
          if (u  .lt.  uMinDamgd) then
              newState   = 2
              gammaFUsed = CgammaF
              do i = 0, 5
                  envlpNegDamgdStress(i+1) = envlpNegStress(i+1) * 
     +                                       (1.0 - gammaFUsed)
              end do
              TstateStrainLow = envlpNegStrain(5+1)
			TstateStressLow = envlpNegStress(5+1)
			TstateStrainHgh = envlpNegStrain(0+1)
			TstateStressHgh = envlpNegStress(0+1)
          else
              newState        = 3
			TstateStrainLow = uMinDamgd
			gammaFUsed      = CgammaF
              do i = 0, 5
                  envlpNegDamgdStress(i+1) = envlpNegStress(i+1) * 
     +                                       (1.0-gammaFUsed)
              end do
              call negEnvlpStress (uMinDamgd,envlpNegDamgdStress, 
     +                             envlpNegStrain,TstateStressLow)
              TstateStrainHgh = Cstrain
			TstateStressHgh = Cstress
          end if
          gammaKUsed          = CgammaK    
          ElasticPosDamgd     = ElasticPos * (1.0-gammaKUsed)
      end if
C -------------------------------------------------------------------------------  C
      if ( (Tstate .eq. 2) .and. (du .gt. 0.0) )  then
          cis = 1
          if (Cstrain  .lt. TminStrainDmnd) then 
			TminStrainDmnd = Cstrain
          end if
          if (TminStrainDmnd .gt. uMinDamgd) then
			TminStrainDmnd = uMinDamgd
          end if
          if (u  .gt.  uMaxDamgd) then
              newState   = 1
              gammaFUsed = CgammaF
              do i = 0,5
                  envlpPosDamgdStress(i+1) = envlpPosStress(i+1) * 
     +                                       (1.0-gammaFUsed)
              end do
              TstateStrainLow = envlpPosStrain(0+1)
			TstateStressLow = envlpPosStress(0+1)
			TstateStrainHgh = envlpPosStrain(5+1)
			TstateStressHgh = envlpPosStress(5+1)
          else
              newState        = 4
              TstateStrainLow = Cstrain
			TstateStressLow = Cstress
			TstateStrainHgh = uMaxDamgd
			gammaFUsed      = CgammaF  
              do i = 0,5
                  envlpPosDamgdStress(i+1) = envlpPosStress(i+1) * 
     +                                       (1.0-gammaFUsed)
              end do
              call posEnvlpStress (uMaxDamgd,envlpPosDamgdStress, 
     +                             envlpPosStrain,TstateStressHgh)
          end if
          gammaKUsed          = CgammaK    
          ElasticNegDamgd     = ElasticNeg * (1.0 - gammaKUsed)
      end if      
C -------------------------------------------------------------------------------  C
      if ( Tstate .eq. 3 ) then    
          if (u .lt. TstateStrainLow) then
              cis = 1
			newState = 2
			TstateStrainLow = envlpNegStrain(5+1)
			TstateStrainHgh = envlpNegStrain(0+1)
			TstateStressLow = envlpNegDamgdStress(5+1)
			TstateStressHgh = envlpNegDamgdStress(0+1)
          else if ( (u .gt. uMaxDamgd) .and. (du .gt. 0.0)  ) then
              cis = 1
			newState = 1
			TstateStrainLow = envlpPosStrain(0+1)
			TstateStressLow = envlpPosStress(0+1)
			TstateStrainHgh = envlpPosStrain(5+1)
			TstateStressHgh = envlpPosStress(5+1)
          else if (du  .gt. 0.0) then 
              cis = 1
			newState = 4
			TstateStrainLow = Cstrain
			TstateStressLow = Cstress
			TstateStrainHgh = uMaxDamgd
			gammaFUsed = CgammaF
              do i = 0, 5
                  envlpPosDamgdStress(i+1) = envlpPosStress(i+1) * 
     +                                       (1.0-gammaFUsed)
              end do
              call posEnvlpStress (uMaxDamgd,envlpPosDamgdStress, 
     +                             envlpPosStrain,TstateStressHgh)
			gammaKUsed       = CgammaK
			ElasticNegDamgd = ElasticNeg * (1.0 - gammaKUsed)
          end if
      end if      
C -------------------------------------------------------------------------------  C
      if ( Tstate .eq. 4 ) then          
          if (u .gt. TstateStrainHgh) then
              cis = 1
			newState = 1
			TstateStrainLow = envlpPosStrain(0+1)
			TstateStressLow = envlpPosDamgdStress(0+1)
			TstateStrainHgh = envlpPosStrain(5+1)
			TstateStressHgh = envlpPosDamgdStress(5+1)
          else if ( (u .lt. uMinDamgd) .and. (du .lt. 0.0)  ) then
              cis = 1
			newState = 2
			TstateStrainLow = envlpNegStrain(5+1)
			TstateStressLow = envlpNegDamgdStress(5+1)
			TstateStrainHgh = envlpNegStrain(0+1)
			TstateStressHgh = envlpNegDamgdStress(0+1)
          else if (du  .lt. 0.0) then 
              cis = 1
			newState = 3
			TstateStrainLow = uMinDamgd
			gammaFUsed = CgammaF  
              do i = 0,5
                envlpNegDamgdStress(i+1) = envlpNegStress(i+1) * 
     +                                     (1.0 - gammaFUsed)
              end do
              call negEnvlpStress (uMinDamgd,envlpNegDamgdStress, 
     +                             envlpNegStrain,TstateStressLow)
			TstateStrainHgh  = Cstrain
			TstateStressHgh  = Cstress
			gammaKUsed       = CgammaK
			ElasticPosDamgd = ElasticPos * (1.0 - gammaKUsed)
          end if
      end if
C -------------------------------------------------------------------------------  C
      end if  
C -------------------------------------------------------------------------------  C
C -------------------------------------------------------------------------------  C
C -------------------------------------------------------------------------------  C
      if (cis .eq. 1) then
	    Tstate = newState
      else
          Tstate = Tstate
      end if
C -------------------------------------------------------------------------------  C
      RETURN
      END
C################################################################################# C
C																				 C	
C			             SUBROUTINE getstate3 (checked)				             C
C						   														 C
C################################################################################# C
      SUBROUTINE getstate3(sfUnload,envlpNegDamgdStress, envlpNegStrain,
     +    TstateStrainLow, TstateStressLow, TstateStrainHgh, 
     +    TstateStressHgh, TminStrainDmnd, ElasticNegDamgd,
     +    state3Strain, state3Stress, rDispN, uForceN, rForceN)
      
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION envlpPosDamgdStress(6),envlpNegDamgdStress(6)
      DIMENSION envlpNegStrain(6)
      DIMENSION state3Stress(4), state3Strain(4)
      INTEGER i
C -------------------------------------------------------------------------------  C
      if (sfUnload .gt. ElasticNegDamgd) then
          s_kmax = sfUnload
      else
          s_kmax = ElasticNegDamgd
      end if
C -------------------------------------------------------------------------------  C
C -------------------------------------------------------------------------------  C
C -------------------------------------------------------------------------------  C
      if (state3Strain(0+1)*state3Strain(3+1) .lt. 0.0) then 
C -------------------------------------------------------------------------------  C
C         trilinear unload reload path expected, first define point for reloading
C
          state3Strain(1+1) = TstateStrainLow * rDispN;
          if (  (rForceN - uForceN) .gt. 1.0D-8)  then 
              state3Stress(1+1) = TstateStressLow * rForceN
          else
              if (TminStrainDmnd .lt. envlpNegStrain(3+1)) then 
                  st1 = TstateStressLow * uForceN * (1.0 + 1.0D-6)
				st2 = envlpNegDamgdStress(4+1)  * (1.0 + 1.0D-6)
                  if (st1 .lt. st2) then
                      state3Stress(1+1) = st1
                  else
                      state3Stress(1+1) = st2
                  end if
              else
                  st1= envlpNegDamgdStress(3+1) * uForceN * (1.0+1.0D-6)
                  st2= envlpNegDamgdStress(4+1) * (1.0+1.0D-6)
                  if (st1 .lt. st2) then
                      state3Stress(1+1) = st1
                  else
                      state3Stress(1+1) = st2
                  end if                
              end if
          end if
C -------------------------------------------------------------------------------  C
C         if reload stiffness exceeds unload stiffness, reduce reload stiffness to make it equal to unload stiffness
C
          condition = (state3Stress(1+1) - state3Stress(0+1)) / 
     +                (state3Strain(1+1) - state3Strain(0+1))
          if (  condition .gt. ElasticNegDamgd )  then
              state3Strain(1+1) = TstateStrainLow + (state3Stress(1+1) -
     +                            state3Stress(0+1)) / ElasticNegDamgd
          end if
C -------------------------------------------------------------------------------  C
C         check that reloading point is not behind point 4 
          if (state3Strain(1+1) .gt. state3Strain(3+1))  then
              du       = state3Strain(3+1) - state3Strain(0+1)
              df       = state3Stress(3+1) - state3Stress(0+1)
              state3Strain(1+1) = state3Strain(0+1) + 0.33*du
              state3Strain(2+1) = state3Strain(0+1) + 0.67*du
              state3Stress(1+1) = state3Stress(0+1) + 0.33*df
              state3Stress(2+1) = state3Stress(0+1) + 0.67*df
          else
C -------------------------------------------------------------------------------  C
              if ( (TminStrainDmnd) .lt. (envlpNegStrain(3+1)) ) then 
                  state3Stress(2+1) = uForceN * envlpNegDamgdStress(4+1)
              else 
                  state3Stress(2+1) = uForceN * envlpNegDamgdStress(3+1)
              end if
              state3Strain(2+1) = TstateStrainHgh - (TstateStressHgh- 
     +                            state3Stress(2+1))/sfUnload
C -------------------------------------------------------------------------------  C
              condition1 = (state3Stress(2+1) - state3Stress(1+1))/
     +                     (state3Strain(2+1) - state3Strain(1+1))
              condition2 = (state3Stress(2+1)-state3Stress(1+1)) / 
     +                     (state3Strain(2+1)-state3Strain(1+1))
C
              if (state3Strain(2+1) .gt. state3Strain(3+1)) then 
                  du       = state3Strain(3+1) - state3Strain(1+1)
                  df       = state3Stress(3+1) - state3Stress(1+1)
                  state3Strain(2+1) = state3Strain(1+1) + 0.5*du
                  state3Stress(2+1) = state3Stress(1+1) + 0.5*df
              else if (condition1 .gt. s_kmax) then
                  du          = state3Strain(3+1) - state3Strain(0+1)
                  df          = state3Stress(3+1) - state3Stress(0+1)
                  state3Strain(1+1) = state3Strain(0+1) + 0.33*du
                  state3Strain(2+1) = state3Strain(0+1) + 0.67*du
                  state3Stress(1+1) = state3Stress(0+1) + 0.33*df
                  state3Stress(2+1) = state3Stress(0+1) + 0.67*df      
              else if (  (state3Strain(2+1) .lt. state3Strain(1+1)) .or.
     +                 (condition2 .lt. 0.0d0)) then
C -------------------------------------------------------------------------------  C
                  if (state3Strain(2+1) .lt. 0.0)  then
                      du      = state3Strain(3+1) - state3Strain(1+1)
                      df      = state3Stress(3+1) - state3Stress(1+1)
                      state3Strain(2+1) = state3Strain(1+1) + 0.5*du
                      state3Stress(2+1) = state3Stress(1+1) + 0.5*df    
                  elseif (state3Strain(1+1) .gt. 0.0) then
                      du      = state3Strain(2+1)-state3Strain(0+1)
                      df      = state3Stress(2+1)-state3Stress(0+1)
                      state3Strain(1+1) = state3Strain(0+1) + 0.5*du
                      state3Stress(1+1) = state3Stress(0+1) + 0.5*df   
                  else
                      avgforce=0.5*(state3Stress(2+1)+state3Stress(1+1))
                      dfr = 0.0
                      if (avgforce .lt. 0.0) then
                          dfr = -avgforce/100.0
                      else 
                          dfr =  avgforce/100.0
                      end if
                      slope12 = (state3Stress(1+1) - state3Stress(0+1))/
     +                          (state3Strain(1+1) - state3Strain(0+1))
                      slope34 = (state3Stress(3+1) - state3Stress(2+1))/
     +                          (state3Strain(3+1) - state3Strain(2+1)) 
                      state3Stress(1+1) = avgforce - dfr
                      state3Stress(2+1) = avgforce + dfr
                      state3Strain(1+1) = state3Strain(0+1) + 
     +                    (state3Stress(1+1)-state3Stress(0+1))/slope12
                      state3Strain(2+1) = state3Strain(3+1) - 
     +                    (state3Stress(3+1)-state3Stress(2+1))/slope34
                  end if
C -------------------------------------------------------------------------------  C
              end if
C -------------------------------------------------------------------------------  C
          end if
C -------------------------------------------------------------------------------  C
      else
          du                  = state3Strain(3+1)-state3Strain(0+1)
          df                  = state3Stress(3+1)-state3Stress(0+1)
          state3Strain(1+1)   = state3Strain(0+1) + 0.33*du
          state3Strain(2+1)   = state3Strain(0+1) + 0.67*du
          state3Stress(1+1)   = state3Stress(0+1) + 0.33*df
          state3Stress(2+1)   = state3Stress(0+1) + 0.67*df
C -------------------------------------------------------------------------------  C
      end if
C -------------------------------------------------------------------------------  C
      checkSlope = state3Stress(0+1)/state3Strain(0+1)
      slope      = 0.0
      i = 0
C      
      do while  (i .lt. 3)
          du = state3Strain(i+1+1) - state3Strain(i+1)
          df = state3Stress(i+1+1) - state3Stress(i+1)
C          
          if (  (du .lt. 0.0) .or. (df.lt.0.0)   ) then
              du          = state3Strain(3+1) - state3Strain(0+1)
		    df          = state3Stress(3+1) - state3Stress(0+1)
              state3Strain(1+1) = state3Strain(0+1) + 0.33*du
		    state3Strain(2+1) = state3Strain(0+1) + 0.67*du
		    state3Stress(1+1) = state3Stress(0+1) + 0.33*df
		    state3Stress(2+1) = state3Stress(0+1) + 0.67*df
		    slope = df/du
		    i = 3
          end if
C          
          if ( (slope .gt. 1.0D-8) .and. (slope .lt. checkSlope) ) then
              state3Strain(1+1) = 0.0 
              state3Stress(1+1) = 0.0
              state3Strain(2+1) = state3Strain(3+1)/2.0 
              state3Stress(2+1) = state3Stress(3+1)/2.0 
          end if
          i = i+1           
      end do
C -------------------------------------------------------------------------------  C
      RETURN
      END
C################################################################################# C
C																				 C	
C			             SUBROUTINE getstate4 (checked)				             C
C						   														 C
C################################################################################# C
      SUBROUTINE getstate4(sfUnload,envlpPosDamgdStress, envlpPosStrain,
     +    TstateStrainLow, TstateStressLow, TstateStrainHgh, 
     +    TstateStressHgh, TmaxStrainDmnd, ElasticPosDamgd,
     +    state4Strain, state4Stress, rDispP, uForceP, rForceP)
      
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION envlpNegDamgdStress(6),envlpPosDamgdStress(6)
      DIMENSION envlpPosStrain(6)
      DIMENSION state4Stress(4), state4Strain(4)
      INTEGER i
C -------------------------------------------------------------------------------  C
      if (sfUnload .gt. ElasticPosDamgd) then
          s_kmax = sfUnload
      else
          s_kmax = ElasticPosDamgd
      end if
C -------------------------------------------------------------------------------  C
C -------------------------------------------------------------------------------  C
C -------------------------------------------------------------------------------  C
      if (state4Strain(0+1)*state4Strain(3+1) .lt. 0.0) then 
C -------------------------------------------------------------------------------  C
C         trilinear unload reload path expected, first define point for reloading
C
          state4Strain(2+1) = TstateStrainHgh * rDispP
          if (uForceP .eq. 0.0) then 
            state4Stress(2+1) = TstateStressHgh * rForceP
          end if
          if (  (rForceP-uForceP) .gt. 1.0D-8)  then 
              state4Stress(2+1) = TstateStressHgh * rForceP
          else
              if (TmaxStrainDmnd  .gt. envlpPosStrain(3+1)) then 
                  st1 = TstateStressHgh * uForceP * (1.0 + 1.0D-6)
				st2 = envlpPosDamgdStress(4+1)  * (1.0 + 1.0D-6)
                  if (st1 .gt. st2) then
                      state4Stress(2+1) = st1
                  else
                      state4Stress(2+1) = st2
                  end if
              else
                  st1= envlpPosDamgdStress(3+1) * uForceP * (1.0+1.0D-6)
                  st2= envlpPosDamgdStress(4+1) * (1.0+1.0D-6)
                  if (st1 .gt. st2) then
                      state4Stress(2+1) = st1
                  else
                      state4Stress(2+1) = st2
                  end if                
              end if
          end if
C -------------------------------------------------------------------------------  C
C         if reload stiffness exceeds unload stiffness, reduce reload stiffness to make it equal to unload stiffness
C
          condition = (state4Stress(3+1) - state4Stress(2+1)) / 
     +                (state4Strain(3+1) - state4Strain(2+1))
          if (  condition .gt. ElasticPosDamgd )  then
              state4Strain(2+1) = TstateStrainHgh - (state4Stress(3+1) -
     +                            state4Stress(2+1)) / ElasticPosDamgd
          end if
C -------------------------------------------------------------------------------  C
C         check that reloading point is not behind point 4 
          if (state4Strain(2+1) .lt. state4Strain(0+1))  then
              du       = state4Strain(3+1) - state4Strain(0+1)
              df       = state4Stress(3+1) - state4Stress(0+1)
              state4Strain(1+1) = state4Strain(0+1) + 0.33*du
              state4Strain(2+1) = state4Strain(0+1) + 0.67*du
              state4Stress(1+1) = state4Stress(0+1) + 0.33*df
              state4Stress(2+1) = state4Stress(0+1) + 0.67*df
          else
C -------------------------------------------------------------------------------  C
              if (TmaxStrainDmnd .gt. envlpPosStrain(3+1)) then 
                  state4Stress(1+1) = uForceP * envlpPosDamgdStress(4+1)
              else 
                  state4Stress(1+1) = uForceP * envlpPosDamgdStress(3+1)
              end if
              state4Strain(1+1) = TstateStrainLow + (- TstateStressLow +
     +                            state4Stress(1+1))/sfUnload
C -------------------------------------------------------------------------------  C
              condition1 = (state4Stress(2+1) - state4Stress(1+1))/
     +                     (state4Strain(2+1) - state4Strain(1+1))
              condition2 = ((state4Stress(2+1) - state4Stress(1+1)) / 
     +              (state4Strain(2+1) - state4Strain(1+1)) )
C
              if (state4Strain(1+1) .lt. state4Strain(0+1)) then 
                  du       = state4Strain(2+1) - state4Strain(0+1)
                  df       = state4Stress(2+1) - state4Stress(0+1)
                  state4Strain(1+1) = state4Strain(0+1) + 0.5*du
                  state4Stress(1+1) = state4Stress(0+1) + 0.5*df
              else if (condition1 .gt. s_kmax) then
                  du          = state4Strain(3+1) - state4Strain(0+1)
                  df          = state4Stress(3+1) - state4Stress(0+1)
                  state4Strain(1+1) = state4Strain(0+1) + 0.33*du
                  state4Strain(2+1) = state4Strain(0+1) + 0.67*du
                  state4Stress(1+1) = state4Stress(0+1) + 0.33*df
                  state4Stress(2+1) = state4Stress(0+1) + 0.67*df
              else if (  (state4Strain(2+1) .lt. state4Strain(1+1)) .or.
     +                 (condition2 .lt. 0.0d0)) then
C -------------------------------------------------------------------------------  C
                  if (state4Strain(1+1) .gt. 0.0)  then
                      du      = state4Strain(2+1) - state4Strain(0+1)
                      df      = state4Stress(2+1) - state4Stress(0+1)
                      state4Strain(1+1) = state4Strain(0+1) + 0.5*du
                      state4Stress(1+1) = state4Stress(0+1) + 0.5*df
                  else if (state4Strain(2+1) .lt. 0.0) then
                      du      = state4Strain(3+1)-state4Strain(1+1)
                      df      = state4Stress(3+1)-state4Stress(1+1)
                      state4Strain(2+1) = state4Strain(1+1) + 0.5*du
                      state4Stress(2+1) = state4Stress(1+1) + 0.5*df
                  else
                      avgforce=0.5*(state4Stress(2+1)+state4Stress(1+1))
                      dfr = 0.0
                      if (avgforce .lt. 0.0) then
                          dfr = -avgforce/100.0
                      else 
                          dfr =  avgforce/100.0
                      end if
                      slope12 = (state4Stress(1+1) - state4Stress(0+1))/
     +                          (state4Strain(1+1) - state4Strain(0+1))
                      slope34 = (state4Stress(3+1) - state4Stress(2+1))/
     +                          (state4Strain(3+1) - state4Strain(2+1)) 
                      state4Stress(1+1) = avgforce - dfr
                      state4Stress(2+1) = avgforce + dfr
                      state4Strain(1+1) = state4Strain(0+1) + 
     +                    (state4Stress(1+1)-state4Stress(0+1))/slope12
                      state4Strain(2+1) = state4Strain(3+1) - 
     +                    (state4Stress(3+1)-state4Stress(2+1))/slope34
                  end if
C -------------------------------------------------------------------------------  C
              end if
C -------------------------------------------------------------------------------  C
          end if
C -------------------------------------------------------------------------------  C
      else
          du                = state4Strain(3+1)-state4Strain(0+1)
          df                = state4Stress(3+1)-state4Stress(0+1)
          state4Strain(1+1) = state4Strain(0+1) + 0.33*du
          state4Strain(2+1) = state4Strain(0+1) + 0.67*du
          state4Stress(1+1) = state4Stress(0+1) + 0.33*df
          state4Stress(2+1) = state4Stress(0+1) + 0.67*df
C -------------------------------------------------------------------------------  C
      end if
C -------------------------------------------------------------------------------  C
      checkSlope = state4Stress(0+1)/state4Strain(0+1)
      slope      = 0.0
      i = 0
C      
      do while  (i .lt. 3)
          du = state4Strain(i+1+1)-state4Strain(i+1)
          df = state4Stress(i+1+1)-state4Stress(i+1)
C          
          if (  (du .lt. 0.0) .or. (df .lt. 0.0)   ) then
              du          = state4Strain(3+1)-state4Strain(0+1)
		    df          = state4Stress(3+1)-state4Stress(0+1)
              state4Strain(1+1) = state4Strain(0+1) + 0.33*du
		    state4Strain(2+1) = state4Strain(0+1) + 0.67*du
		    state4Stress(1+1) = state4Stress(0+1) + 0.33*df
		    state4Stress(2+1) = state4Stress(0+1) + 0.67*df
		    slope = df/du
		    i = 3
          end if
C          
          if ( (slope .gt. 1.0D-8) .and. (slope .lt. checkSlope) ) then
              state4Strain(1+1) = 0.0 
              state4Stress(1+1) = 0.0
              state4Strain(2+1) = state4Strain(3+1)/2.0 
              state4Stress(2+1) = state4Stress(3+1)/2.0 
          end if
          i = i+1           
      end do
C -------------------------------------------------------------------------------  C
      RETURN
      END
C################################################################################# C
C																				 C	
C			        SUBROUTINE revertToStart (checked)				             C
C						   														 C
C################################################################################# C
      SUBROUTINE revertToStart(envlpPosStress, envlpPosStrain, 
     +    envlpNegStress, envlpNegStrain, ElasticPosDamgd, 
     +    ElasticNegDamgd, gammaKUsed, gammaFUsed,
     +    ElasticPos, ElasticNeg,
     +    Cstate, Cstrain, Cstress, CstrainRate, CstateStrainLow, 
     +    CstateStressLow, CstateStrainHgh, CstateStressHgh,
     +    CminStrainDmnd, CmaxStrainDmnd, Cenergy, CgammaK, CgammaD,
     +    CgammaF, CnCycle,Ttangent, uMaxDamgd, uMinDamgd)
C      
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION envlpPosStress(6), envlpPosStrain(6)
      DIMENSION envlpNegStress(6), envlpNegStrain(6)
C      
      INTEGER Tstate, Cstate, DmgCyc
C -------------------------------------------------------------------------------  C
      Cstate      = 0
      Cstrain     = 0.0d0
      Cstress     = 0.0d0
      CstrainRate = 0.0d0
C
      CstateStrainLow = envlpNegStrain(0+1)
      CstateStressLow = envlpNegStress(0+1)
      CstateStrainHgh = envlpPosStrain(0+1)
      CstateStressHgh = envlpPosStress(0+1)
      CminStrainDmnd  = envlpNegStrain(1+1)
      CmaxStrainDmnd  = envlpPosStrain(1+1)
C
      Cenergy = 0.0d0
      CgammaK = 0.0d0
      CgammaD = 0.0d0
      CgammaF = 0.0d0
      CnCycle = 0.0d0
C
      Ttangent   = envlpPosStress(0+1)/envlpPosStrain(0+1)
      gammaKUsed = 0.0d0
      gammaFUsed = 0.0d0
C	
      ElasticPosDamgd = ElasticPos
      ElasticNegDamgd = ElasticNeg
      uMaxDamgd        = CmaxStrainDmnd
      uMinDamgd        = CminStrainDmnd
C -------------------------------------------------------------------------------  C
      RETURN
      END
C################################################################################# C
C																				 C	
C			  SUBROUTINE DMGEV TO CALCULATE UNIAXIAL DAMAGE	EVOLUTION			 C
C						   														 C
C################################################################################# C
      SUBROUTINE revertToLastCommit(Cstate, Cstrain, Cstress, 
     +    CstateStrainLow,CstateStressLow,CstateStrainHgh,
     +    CstateStressHgh, CminStrainDmnd, CmaxStrainDmnd, Cenergy, 
     +    CgammaK, CgammaD, CgammaF,CnCycle, CstrainRate,
     +    envlpPosDamgdStress, envlpNegDamgdStress,
     +    uMaxDamgd,uMinDamgd, ElasticPosDamgd,ElasticNegDamgd,
     +    STATEV,NSTATV)
      
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION STATEV(NSTATV)
      DIMENSION envlpPosDamgdStress(6),envlpNegDamgdStress(6)
      INTEGER   Tstate, Cstate
C -------------------------------------------------------------------------------  C
      Cstate                  = STATEV(1)
      CstateStrainLow         = STATEV(2)
      CstateStressLow         = STATEV(3)  
      CstateStrainHgh         = STATEV(4)  
      CstateStressHgh         = STATEV(5)  
      CminStrainDmnd          = STATEV(6)  
      CmaxStrainDmnd          = STATEV(7)  
      Cenergy                 = STATEV(8)  
      Cstress                 = STATEV(9) 
      Cstrain                 = STATEV(10) 
      CgammaK                 = STATEV(11) 
      CgammaD                 = STATEV(12) 
      CgammaF                 = STATEV(13) 
      ElasticPosDamgd         = STATEV(14) 
      ElasticNegDamgd         = STATEV(15) 
      uMaxDamgd               = STATEV(16) 
      uMinDamgd               = STATEV(17) 
      envlpPosDamgdStress(1)  = STATEV(18) 
      envlpPosDamgdStress(2)  = STATEV(19) 
      envlpPosDamgdStress(3)  = STATEV(20) 
      envlpPosDamgdStress(4)  = STATEV(21) 
      envlpPosDamgdStress(5)  = STATEV(22) 
      envlpPosDamgdStress(6)  = STATEV(23) 
      envlpNegDamgdStress(1)  = STATEV(24) 
      envlpNegDamgdStress(2)  = STATEV(25) 
      envlpNegDamgdStress(3)  = STATEV(26) 
      envlpNegDamgdStress(4)  = STATEV(27) 
      envlpNegDamgdStress(5)  = STATEV(28) 
      envlpNegDamgdStress(6)  = STATEV(29) 
      CnCycle                 = STATEV(30)
      CstrainRate             = STATEV(31)
C -------------------------------------------------------------------------------  C
      RETURN
      END
C################################################################################# C
C																				 C	
C			        SUBROUTINE updateDmg     (checked)				             C
C						   														 C
C################################################################################# C      
      SUBROUTINE updateDmg(strain, dstrain, CnCycle, energyCapacity,  
     +       Tenergy, ElasticPos, ElasticNeg, elasticStrainEnergy,
     +       gammaK1, gammaK2, gammaK3, gammaK4, gammaKLimit,
     +       gammaD1, gammaD2, gammaD3, gammaD4, gammaDLimit,
     +       gammaF1, gammaF2, gammaF3, gammaF4, gammaFLimit,
     +       TmaxStrainDmnd,  TminStrainDmnd, DmgCyc,
     +       envlpPosDamgdStress, envlpPosStrain, envlpNegDamgdStress, 
     +       envlpNegStrain, TnCycle, TgammaK, TgammaD, TgammaF)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION envlpPosStrain(6),      envlpNegStrain(6)
      DIMENSION envlpPosDamgdStress(6), envlpNegDamgdStress(6)
C -------------------------------------------------------------------------------  C
      tes = 0.0
C
      if (TmaxStrainDmnd .gt. -TminStrainDmnd) then
          umaxAbs =  TmaxStrainDmnd
      else
          umaxAbs = -TminStrainDmnd
      end if    
C
      if ( envlpPosStrain(4+1) .gt. -envlpNegStrain(4+1) ) then
          uultAbs = envlpPosStrain(4+1)
      else
          uultAbs = -envlpNegStrain(4+1)
      end if    
C
      TnCycle = CnCycle + abs(dstrain)/(4.0D0 * umaxAbs)
C -------------------------------------------------------------------------------  C
C
      if  (   ((strain .lt. uultAbs) .and. (strain .gt. -uultAbs)) .and.
     +   (Tenergy .lt. energyCapacity)    )    then
C -------------------------------------------------------------------------------  C
          TgammaK = gammaK1 * (umaxAbs/uultAbs) ** gammaK3
	    TgammaD = gammaD1 * (umaxAbs/uultAbs) ** gammaD3
	    TgammaF = gammaF1 * (umaxAbs/uultAbs) ** gammaF3
C
          if (   (Tenergy .gt. elasticStrainEnergy)  ) then
C             
              tes = ((Tenergy-elasticStrainEnergy)/energyCapacity)
C
              TgammaK = TgammaK + gammaK2 * tes ** gammaK4
		    TgammaD = TgammaD + gammaD2 * tes ** gammaD4
		    TgammaF = TgammaF + gammaF2 * tes ** gammaF4
C              
           end if
C -------------------------------------------------------------------------------  C
           call posEnvlpStress (TmaxStrainDmnd,envlpPosDamgdStress, 
     +                          envlpPosStrain,Term1)
           call negEnvlpStress (TminStrainDmnd,envlpNegDamgdStress, 
     +                          envlpNegStrain,Term2)
C          
           s_kminP = Term1/TmaxStrainDmnd
           s_kminN = Term2/TminStrainDmnd
C
           if ( (s_kminP/ElasticPos) .gt. (s_kminN/ElasticNeg) ) then
              s_kmin = (s_kminP/ElasticPos)
           else
              s_kmin = (s_kminN/ElasticNeg)
           end if
C
           if ( 0.0 .gt. (1.0 - s_kmin) )  then
              gammaKLimEnv = 0.0
           else
              gammaKLimEnv = (1.0 - s_kmin)
           end if
C
           if (TgammaK .lt. gammaKLimit)  then
                s_k1 = TgammaK
           else
                s_k1 = gammaKLimit
           end if
C
           if ( s_k1 .lt. gammaKLimEnv )   then
                TgammaK =  s_k1
           else
                TgammaK =  gammaKLimEnv
           end if
C
           if (TgammaD .lt. gammaDLimit)  then
                TgammaD = TgammaD
           else
                TgammaD = gammaDLimit
           end if
C
           if (TgammaF .lt. gammaFLimit)  then
                TgammaF =  TgammaF
           else
                TgammaF =  gammaFLimit
           end if
C -------------------------------------------------------------------------------  C
      else if ( (strain .lt. uultAbs) .and. (strain .gt. -uultAbs) )then
           call posEnvlpStress (TmaxStrainDmnd,envlpPosDamgdStress, 
     +                         envlpPosStrain,Term1)
           call negEnvlpStress (TminStrainDmnd,envlpNegDamgdStress, 
     +                         envlpNegStrain,Term2)
C          
           s_kminP = Term1/TmaxStrainDmnd
           s_kminN = Term2/TminStrainDmnd
C
           if ( (s_kminP/ElasticPos) .ge. (s_kminN/ElasticNeg) ) then
                s_kmin = (s_kminP/ElasticPos)
           else
                s_kmin = (s_kminN/ElasticNeg)
           end if
C
           if (0.0 .gt. (1.0-s_kmin)) then
                gammaKLimEnv =  0.0
           else
                gammaKLimEnv =  (1.0-s_kmin)
           end if
C
           if (gammaKLimit .lt. gammaKLimEnv) then
                TgammaK = gammaKLimit
           else
                TgammaK = gammaKLimEnv
           end if
C
           TgammaD = gammaDLimit
           TgammaF = gammaFLimit
C -------------------------------------------------------------------------------  C
      end if
C     TgammaK = 0.0d0
C -------------------------------------------------------------------------------  C
      RETURN
      END

C======================================================================        
      
      
      
      
      
C-----------------------------------------------------------------------------------------C
C-----------------------------------------------------------------------------------------C
C-----------------------------------------------------------------------------------------C
C !!! NOTE: N_ELEM HAS TO BE CHANGED ACCORDING TO THE UEL !!!!!
C ==============================================================
C
       SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,
     2 TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,MATERL,NDI,NSHR,NTENS,
     3 NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,CELENT,
     4 DFGRD0,DFGRD1,NOEL,NPT,KSLAY,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1    DDSDDE(NTENS,NTENS),
     2    DDSDDT(NTENS),DRPLDE(NTENS),
     3    STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     4    PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
C 
       PARAMETER (ONE=1.0,TWO=2.0,THREE=3.0,SIX=6.0, HALF =0.5,
     1 N_ELEM=100000,NSDV=120)
       DATA NEWTON,TOLER/40,1.D-2/ 
C       
       COMMON/KUSER/USRVAR(N_ELEM,NSDV,8)
C 
C ----------------------------------------------------------- 
C          MATERIAL PROPERTIES
C ----------------------------------------------------------- 
C          PROPS(1) - YOUNG'S MODULUS 
C          PROPS(2) - POISSON RATIO 
C ----------------------------------------------------------- 
C
C	ELASTIC PROPERTIES
C
       EMOD=PROPS(1)
       ENU=PROPS(2)
       EG=EMOD/(TWO*(ONE+ENU))
       EG2=EG*TWO
       ELAM=EG2*ENU/(ONE-TWO*ENU)
C
C	STIFFNESS TENSOR
C
       DO K1=1, NTENS
        DO K2=1, NTENS
         DDSDDE(K2, K1)=0.0
        END DO
       END DO
C
       DO K1=1, NDI
        DO K2=1, NDI
         DDSDDE(K2, K1)=ELAM
        END DO
        DDSDDE(K1, K1)=EG2+ELAM
       END DO 
C
       DO K1=NDI+1, NTENS
        DDSDDE(K1, K1)=EG
       END DO
C
C	CALCULATE STRESSES
C
       DO K1=1, NTENS
          DO K2=1, NTENS
              STRESS(K2)=STRESS(K2)+DDSDDE(K2, K1)*DSTRAN(K1)
          END DO
       END DO 
C
       NELEMAN=NOEL-N_ELEM
C       
       DO I=1,NSTATV
        STATEV(I)=USRVAR(NELEMAN,I,NPT)
       END DO
C       
       RETURN
       END      