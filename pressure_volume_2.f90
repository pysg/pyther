Desde el Main Program, que incluye estos commons...
      COMMON /MODEL/ NMODEL
      COMMON/TCGC/TCGC(MAXC),PCGC
      COMMON/CRIT/TC(MAXC),PC(MAXC),DC(MAXC)
    COMMON/PAEP/ TP(NA),PP(NA),ZP(NA),DLP(NA),DVP(NA),NPAEP

Se llama así a la rutina para generar las campanas de saturación de ambos componentes en un sistema binario
     CALL PvcurveNewton(nout,1)
    CALL PvcurveNewton(nout,2)

(como el arranque alternativo desde 50K), etc.

    subroutine PvcurveNewton(nout,icomp) ! S para los primeros puntos debería ser ln(Vv/Vl)
      implicit DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (nco=2, RGAS=0.08314472d0, NA=2)
    DIMENSION TCmod(nco),rn(nco),SENS(3),dold(3),dot(3)
    dimension Arn(nco),ArVn(nco),ArTn(nco)
      DIMENSION XVAR(3),F(3),dFdS(3),dXdS(3),delX(3),RJAC(3,3),XOLD(3)
      DIMENSION AJ(3,3),IPIV(3),db(3)
      DIMENSION FUGx(2),FUGy(2)
      DIMENSION FUGTx(2),FUGTy(2)
      DIMENSION FUGVx(2),FUGVy(2)
      DIMENSION DFGN(2,2)
      DIMENSION DPDNx(2),DPDNy(2)
    LOGICAL oscil, volcheck
      COMMON/CRIT/TC(nco),PC(nco),DC(nco)
      COMMON/TCGC/TCGC(nco)
    COMMON /MODEL/ NMODEL
    COMMON/Pder/ DPDN(2),DPDT,DPDV
    COMMON/PAEP/ TP(NA),PP(NA),ZP(NA),DLP(NA),DVP(NA),NPAEP
    volcheck=.false.    ! false for visual executable. true to study azeotropy
    volog=0.0d0
    NTEMP=0
    NDER=1
    XOLD=0.0D0
    TOLF=1.0D-3    ! for residuals  1.0D-3
    TOL= 1.0D-4    ! for variables  1.0D-4
    N=3    ! DLSARG CONSTANTS
    INCX=1
    LDA=3
    ldb=3
c   IPATH=1
    TCmod=TC
    if(TCGC(1).ne.0.0D0)TCmod=TCGC
    inf=2
    if(icomp.eq.2)inf=1
    rn (inf) =0.0
    rn(icomp)=1.0
    WRITE (nout,*)
    WRITE (nout,*) '   T(K)    Pv(bar)    rhoL     rhoV'
    WRITE (nout,*) 'VAP'
c    Starting from critical point    c    c
    WRITE(nout,15) TCmod(icomp),Pc(icomp),Dc(icomp),Dc(icomp)
    T = 0.9999 * TCmod(icomp)
    Vc = 1/DC(icomp)
    Vv = 1.03*Vc
    Zc = Pc(icomp)*Vc/RGAS/T
    aaa = min(0.89 + (Zc - 0.2) / 2, 0.95)
    Vl = aaa*Vc
    NS = 3
    delXS = 0.10
    go to 2
c    c    c    c    c    c    c    c    c    c    c
c
c    c    c    c    c    c    c    c    c    c    c
c
 2    XVAR = log([T, Vl, Vv])
    DFDS = 0.0D0
    DFDS(3) = 1.0D0
    RJAC = 0.0D0
    RJAC(3,NS) = 1.0D0
 1    NITER=0
    T=exp(XVAR(1))
    Vl=exp(XVAR(2))
    Vv=exp(XVAR(3))
    FMAXOLD=8.0D0
    FMAX=7.0D0
    DMAXOLD=8.0D0
    DMAX=7.0D0
    F(3)=0.0D0
    delX=0.0
    NV=0
c    Newton procedure for solving the present point
    DO WHILE (DMAX.GT.TOL.or.FMAX.GT.TOLF)
        if(NV.GE.5.and.T.LT.0.4*Tcmod(icomp))return
        IF ((FMAX.GT.FMAXOLD.and.DMAX.GT.DMAXOLD).OR.NITER.eq.10) THEN
            delXS=0.8*delXS
            if(abs(delXS).lt.0.001)return
            if(XOLD(1).NE.0.0D0)go to 7
        END IF
        NITER=NITER+1
 21        CALL XTVTERMO(2,T,Vl,Pl,rn,FUGx,FUGTx,FUGVx,DFGN)
        if(Pl.eq.0.0d0)Pl=1.0D-17
        DPDTx=DPDT
        DPDVx=DPDV

        CALL XTVTERMO(2,T,Vv,Pv,rn,FUGy,FUGTy,FUGVy,DFGN)
        DPDTy=DPDT
        DPDVy=DPDV
        if((Pl.lt.Pv/2.and.Pv-Pl.gt.1.d-8).or.
    1     (Pl.gt.1.5*Pv.and.Pl-Pv.gt.1.d-8))then    ! Pv.lt.1.0D-8.or.
            NV=NV+1
            V1=Vl+(Pv-Pl)/DPDVx    ! very important for convergence and accuracy of P at low T
            V2=(Vl+B)/2
            Vl=max(V1,V2)
            XVAR(2)=log(Vl)
            go to 21
        end if
        
        
        if(Pl.lt.0)then
            F(1)=-TOLF
        else
            F(1)=log(Pl/Pv)
        end if
        
        F(2) = FUGx(icomp) - FUGy(icomp)
        
        RJAC(1,1) = DPDTx/Pl-DPDTy/Pv
        RJAC(1,1) = T*RJAC(1,1)
        RJAC(1,2) = Vl*DPDVx/Pl
        RJAC(1,3) = -Vv*DPDVy/Pv
C
        RJAC(2,1) = FUGTx(icomp)-FUGTy(icomp)
        RJAC(2,1) = T*RJAC(2,1)
        RJAC(2,2) = Vl*FUGVx(icomp)
        RJAC(2,3) = -Vv*FUGVy(icomp)
C
        dold=delx
c       CALL DLSARG (N, RJAC, LDA, -F, IPATH, delX)
        db = -F
        AJ=RJAC
        call dgesv( N, 1, AJ, lda, ipiv, db, ldb, info )
        if (info.ne.0) write(6,*)"error with dgesv in parameter ",info
        delX = db
c        
        dot(1:3)=delx(1:3)*dold(1:3)
        dotp=sum(dot(1:3))
        dxmax=maxval(abs(delx))
c    oscillations are recognized by negative dotp
        oscil=.false.
        if(maxval(dot).le.0.0.or.dotp.lt.-1.d-5)oscil=.true.
        if (oscil.and.dxmax.ge.1.0d-4) then
            delx=delx/2
        end if
      if(dxmax.ge.1.0d-3.and.dotp.eq.0.0) delx=delx/2  !half step for first iteration (precaution)
c        if (dmax.ge.5.0d-3) then
c            delx=5.0d-3*delx/dmax
c        end if
    if(delXS.eq.0.015.and.maxval(abs(delX)).gt.7)delX=6*delX/abs(delX(3))
 17        XVAR=XVAR+delX
 18        T=exp(XVAR(1))
        Vl=exp(XVAR(2))
        IF(T.GT.TCMOD(icomp).or.Vl.ge.Vc)THEN
            delX=delX/2
            XVAR=XVAR-delX
            go to 18
        END IF
        Vv=exp(XVAR(3))
        FMAXOLD=FMAX
        FMAX=MAXVAL(ABS(F))
        DMAXOLD=DMAX
        DMAX=MAXVAL(ABS(DELX))
        if(FMAX.lt.1.0D-1.and.FMAX*DMAX.lt.1.0D-11)exit
    END DO
    RHOL=1/Vl
    RHOV=1/Vv
c    calculation of the relative volatility from fug. coef. of the diluted component
            RT = RGAS*T
            ZL=Pv*Vl/RT
            ZV=Pv*Vv/RT
        call ArVnder(NDER,NTEMP,rn,Vl,T,Ar,ArV,ArTV,ArV2,Arn,ArVn,ArTn,DFGN)
            FUGl=Arn(inf)/RT - log(ZL)  ! ln(phi) inf dil in L
        call ArVnder(NDER,NTEMP,rn,Vv,T,Ar,ArV,ArTV,ArV2,Arn,ArVn,ArTn,DFGN)
            FUGv=Arn(inf)/RT - log(ZV)  ! ln(phi) inf dil in V
            vold=volog
            volog=FUGl-FUGv
            if(vold*volog.lt.0.0)then
                NPAEP=NPAEP+1
                TP(NPAEP)=T-volog*(TOLD-T)/(vold-volog)
                PP(NPAEP)=Pv-volog*(POLD-Pv)/(vold-volog)
                ZP(NPAEP)=rn(1)
                DLP(NPAEP)=RHOL-volog*(DLOLD-RHOL)/(vold-volog)
                DVP(NPAEP)=RHOV-volog*(DVOLD-RHOV)/(vold-volog)
            end if
c    c    c    c    c    c    c    c    c    c    c    c    c    c    c    c    c
    IF(volcheck)THEN
        WRITE (nout,15) T,Pv,RHOL,RHOV,NS,NITER,volog
    ELSE
        WRITE (nout,15) T,Pv,RHOL,RHOV,NS,NITER
    END IF
c    CALL DLSARG (N, RJAC, LDA, dFdS, IPATH, dXdS)
        db = dFdS
        AJ=RJAC
        call dgesv( N, 1, AJ, lda, ipiv, db, ldb, info )
        if (info.ne.0) write(6,*)"error with dgesv in parameter ",info
        dXdS = db
c        
    NSOLD=NS
    IF(VV.GT.1.1*VC)THEN  ! .and.T.GT.0.3*Tcmod(icomp)
        NS=1
    ELSE
        NS=3
C        SENS(1:3)=ABS(DXDS(1:3)/XVAR(1:3))
C        SENS(2)=SENS(2)/3    ! OTHERWISE lnVv is always specified and too many points are required
C        SENS(1)=(T95/300)*30*SENS(1)    ! OTHERWISE lnVv is always specified and too many points are required
C        if(NMODEL.LE.3)SENS(1)=5*SENS(1)
C        NS = IDMAX (3, SENS, INCX)
    END IF
    if(NS.ne.NSOLD)then
        dSdSold=dXdS(NS)
        dXdS=dXdS/dSdSold
        RJAC(3,1:3)=0.0D0
        RJAC(3,NS)=1.0D0
    end if
    NITER=min(niter,10)
    delXS=delXS*5/NITER/dXdS(NSOLD)
    IF(delXS.LT.0)then
        if(NS.EQ.1)
            delXS=max(delXS,-0.008) ! Max lnT decrease allowed
        if(NS.eq.2)
            delXS=max(delXS,-0.01) ! Max lnVl decrease allowed: 0.01  ! DO NOT CHANGE !!!
        if(NS.eq.3)
            delXS=max(delXS,-0.10) ! Max lnVv decrease allowed: 0.20
    ELSE
        if(NS.EQ.1)
            delXS=min(delXS,0.01) ! Max lnT increase allowed: 0.5 K
        if(NS.eq.2)
            delXS=min(delXS,0.05) ! Max lnVl increase allowed: 0.01
        if(NS.eq.3)
            delXS=min(delXS,0.20) ! Max lnVv increase allowed: 0.20
    END IF
    XOLD=XVAR
    TOLD=T
    POLD=Pv
    DVOLD=RHOV
    DLOLD=RHOL
    NV=0
  7    S=XOLD(NS)+delXS
    XVAR=XOLD+dXdS*delXS/dXdS(NS)    ! Initial estimates for the 3 variables in the next point
 8    T=exp(XVAR(1))
    Vl=exp(XVAR(2))
    Vv=exp(XVAR(3))
    if(Pv.gt.1.D-20.and.(T.gt.0.25*TCmod(icomp).or.T>250)     ! modified April 6, 2016
    1                              .and.T.lt.TCmod(icomp))go to 1
    if(T.gt.TCmod(icomp))WRITE(nout,15) TCmod(icomp),Pc(icomp),Dc(icomp),Dc(icomp)
    WRITE (nout,*)
15    FORMAT(F8.3,X,E11.4E3,F9.4,X,E11.4E3,2I4,X,E11.3E2)
    end



! --------------------------------------------
CALL XTVTERMO(INDIC,T,V, P, rn, FUGLOG,DLFUGT,DLFUGV,DLFUGX)
CALL XTVTERMO(2,    T,Vl,Pl,rn, FUGx,  FUGTx, FUGVx, DFGN)       
CALL XTVTERMO(2,    T,Vv,Pv,rn, FUGy,  FUGTy, FUGVy, DFGN)

! Aquí "x" ademas de "y" se refieren a la fase de líqudio y vapor, respectivamente.
! A continuación se muestran los nombre que se le dan a las variables de salidad en XTVTERMO,
! que son confusos ademas.

FUGLOG = FUGx  = FUGy
DLFUGT = FUGTx = FUGTy
DLFUGV = FUGVx = FUGVy
DLFUGX = DFGN  = DFGN

! En el caso de sustancias puras DFGN no se utiliza, puesto que DFGN es:
! DLFUGX    comp-derivative of FUGLOG (const t & v)

! --------------------------------------------
DPDV = -ArV2 - RT * TOTN / V ** 2
DPDT = -ArTV + TOTN * RGAS / V
! --------------------------------------------
DPDT = DPDTx = DPDTy
DPDV = DPDVx = DPDVy
! --------------------------------------------
RJAC[1, 1] = T * (DPDTx / Pl - DPDTy / Pv)
RJAC[1, 2] = Vl * DPDVx / Pl
RJAC[1, 3] = -Vv * DPDVy / Pv

RJAC[2, 1] = T * (FUGTx[i] - FUGTy[i])
RJAC[2, 2] = Vl * FUGVx[i]
RJAC[2, 3] = -Vv * FUGVy[i]
! --------------------------------------------




      SUBROUTINE XTVTERMO(INDIC,T,V,P,rn,
    1                    FUGLOG,DLFUGT,DLFUGV,DLFUGX)
C
C-------parameters of XTVTERMO (crit. point, LLV and CEP calculations)
C
C       rn        mixture mole numbers                     (input)
C       t            temperature (k)                       (input)
C       v            volume        (L)                     (input)
C       p            pressure    (bar)                     (output)
C       FUGLOG    vector of log. of fugacities (x*phi*P)   (output)    INDIC < 5
C       DLFUGT    t-derivative of FUGLOG (const. vol,n)    (output)    INDIC = 2 or 4
C       DLFUGV    vol-derivative of FUGLOG (const temp,n)  (output)    INDIC < 5
C       DLFUGX    comp-derivative of FUGLOG (const t & v)  (output)    INDIC > 2
C---------------------------------------------------
C---  MODIFIED AND CORRECTED july 2005
C---
C---------------------------------------------------
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    PARAMETER (MAXC=2,nco=2,RGAS=0.08314472d0)
    DIMENSION DLFUGX(MAXC,MAXC)
    DIMENSION FUGLOG(MAXC),DLFUGT(MAXC),DLFUGV(MAXC)
    dimension rn(nco),Arn(nco),ArVn(nco),ArTn(nco),Arn2(nco,nco)
    COMMON / MODEL / NMODEL
    COMMON / NG / NGR
    COMMON / Pder / DPDN(nco),DPDT,DPDV
    NG=NGR
    NC=2
    IF(NMODEL.EQ.5.OR.NMODEL.EQ.7) CALL PARAGC(T,NC,NG,1)      
    NTEMP=0
      IGZ=0
      NDER=1
      IF (INDIC.GT.2) NDER=2
      IF (INDIC.EQ.2 .OR. INDIC.EQ.4) NTEMP=1
    TOTN = sum(rn)
      RT = RGAS*T
    call ArVnder(NDER,NTEMP,rn,V,T,Ar,ArV,ArTV,ArV2,Arn,ArVn,ArTn,Arn2)
      P = TOTN*RT/V - ArV
      DPDV = -ArV2-RT*TOTN/V**2
      IF(INDIC.GT.4)GOTO 62
c      Z = P*V/(TOTN*RT)
    DPDT = -ArTV+TOTN*RGAS/V
    DO 60 I=1,NC
    IF(RN(I).EQ.0.0)GOTO 60
C        FUGLOG(I)=-LOG(Z)+Arn(I)/RT + log(rn(I)/TOTN) + log(P)
C        FUGLOG(I)=Arn(I)/RT + log(rn(I)/TOTN) + log(P/Z) this crashes at very low T LLV when Z=P=0.000000...
        FUGLOG(I)=Arn(I)/RT + log(rn(I)) + log(RT/V)
        DPDN(I) = RT/V-ArVn(I)
        DLFUGV(I)=-DPDN(I)/RT                    ! term DPDV/P is cancelled out
        IF(NTEMP.EQ.0) GOTO 60
        DLFUGT(I)=(ArTn(I)-Arn(I)/T)/RT+1.D0/T    ! term DPDT/P is cancelled out
   60 CONTINUE
   62 IF(NDER.LT.2) GOTO 64
      DO 63 I=1,NC
      DO 61 K=I,NC
        DLFUGX(I,K)=Arn2(I,K)/RT        ! term 1/TOTN is cancelled out
   61        DLFUGX(K,I)=DLFUGX(I,K)
        DLFUGX(I,I)=DLFUGX(I,I)+1.0/rn(I)
   63 CONTINUE
   64 RETURN
      END
C


CALL XTVTERMO(2,T,Vl,Pl,rn,FUGx,FUGTx,FUGVx,DFGN)


if(Pl.eq.0.0d0)Pl=1.0D-17
DPDTx = DPDT
DPDVx = DPDV
        
CALL XTVTERMO(2,T,Vv,Pv,rn,FUGy,FUGTy,FUGVy,DFGN)
DPDTy = DPDT
DPDVy = DPDV
        


RJAC[1,1] = T * (DPDTx/Pl - DPDTy/Pv)
RJAC[1,2] = Vl * DPDVx/Pl
RJAC[1,3] = -Vv * DPDVy/Pv

RJAC[2,1] = T * (FUGTx[i] - FUGTy[i])
RJAC[2,2] = Vl * FUGVx[i]
RJAC[2,3] = -Vv * FUGVy[i]






    SUBROUTINE HelmRKPR(NDE,NTD,rn,V,T,Ar,ArV,ArTV,ArV2,Arn,ArVn,ArTn,Arn2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (nco=2,RGAS=0.08314472d0)
    dimension rn(nco),Arn(nco),ArVn(nco),ArTn(nco),Arn2(nco,nco)
    dimension dBi(nco),dBij(nco,nco),dD1i(nco),dD1ij(nco,nco)
    dimension dDi(nco),dDij(nco,nco),dDiT(nco)
    dimension aij(nco,nco),daijdT(nco,nco),daijdT2(nco,nco)
    COMMON /rule/ncomb
    NC=2
    TOTN = sum(rn)
    call DELTAnder(nc,rn,D1,dD1i,dD1ij)
    D2=(1-D1)/(1+D1)
c
c
c   Comparison to test and debug cubic mixing rules
c   rn=[0.65,0.35]
c   T=460.0d0
c       call Bnder(nc,rn,Bmix,dBi,dBij)
c       call Bcubicnder(nc,rn,Bmix,dBi,dBij)
c       call DandTnder(NTD,nc,T,rn,D,dDi,dDiT,dDij,dDdT,dDdT2)
c       call DCubicandTnder(NTD,nc,T,rn,D,dDi,dDiT,dDij,dDdT,dDdT2)
c
c
c
    if(ncomb.lt.2)then
        call Bnder(nc,rn,Bmix,dBi,dBij)
        call DandTnder(NTD,nc,T,rn,D,dDi,dDiT,dDij,dDdT,dDdT2)
    else
c       call Bcubicnder(nc,rn,Bmix,dBi,dBij)
c       call DCubicandTnder(NTD,nc,T,rn,D,dDi,dDiT,dDij,dDdT,dDdT2)
    end if
c   The f's and g's used here are for Ar, not F (reduced Ar)                    ***********
c   This requires to multiply by R all g, f and its derivatives as defined by Mollerup ****
    f=log((V+D1*Bmix)/(V+D2*Bmix))/Bmix/(D1-D2)
    g=RGAS*log(1-Bmix/V)
    fv=-1/((V+D1*Bmix)*(V+D2*Bmix))
    fB=-(f+V*fv)/Bmix
    gv=RGAS*Bmix/(V*(V-Bmix))
    fv2=(-1/(V+D1*Bmix)**2+1/(V+D2*Bmix)**2)/Bmix/(D1-D2)
    gv2=RGAS*(1/V**2-1/(V-Bmix)**2)
C   DERIVATIVES OF f WITH RESPECT TO DELTA1
    auxD2=(1+2/(1+D1)**2)
    fD1=(1/(V+D1*Bmix)+2/(V+D2*Bmix)/(1+D1)**2)-f*auxD2
    fD1=fD1/(D1-D2)
    fBD1=-(fB*auxD2+D1/(V+D1*Bmix)**2+2*D2/(V+D2*Bmix)**2/(1+D1)**2)
    fBD1=fBD1/(D1-D2)
    fVD1=-(fV*auxD2+1/(V+D1*Bmix)**2+2/(V+D2*Bmix)**2/(1+D1)**2)/(D1-D2)
    fD1D1=4*(f-1/(V+D2*Bmix))/(1+D1)**3+Bmix*(-1/(V+D1*Bmix)**2+
    1       4/(V+D2*Bmix)**2/(1+D1)**4)-2*fD1*(1+2/(1+D1)**2)
    fD1D1=fD1D1/(D1-D2)
c   Reduced Helmholtz Energy and derivatives
    Ar=-TOTN*g*T-D*f
    ArV=-TOTN*gv*T-D*fv
    ArV2=-TOTN*gv2*T-D*fv2
c
    AUX=RGAS*T/(V-Bmix)
    FFB=TOTN*AUX-D*fB
    FFBV=-TOTN*AUX/(V-Bmix)+D*(2*fv+V*fv2)/Bmix
    FFBB=TOTN*AUX/(V-Bmix)-D*(2*f+4*V*fv+V**2*fv2)/Bmix**2
    do i=1,nc
    Arn(i)=-g*T+FFB*dBi(i)-f*dDi(i)-D*fD1*dD1i(i)
    ArVn(i)=-gv*T+FFBV*dBi(i)-fv*dDi(i)-D*fVD1*dD1i(i)
    IF (NDE.EQ.2) THEN
    do j=1,i
    Arn2(i,j)=AUX*(dBi(i)+dBi(j))-fB*(dBi(i)*dDi(j)+dBi(j)*dDi(i))
     &      +FFB*dBij(i,j)+FFBB*dBi(i)*dBi(j)-f*dDij(i,j)      
      Arn2(i,j)=Arn2(i,j)-D*fBD1*(dBi(i)*dD1i(j)+dBi(j)*dD1i(i))
     &      -fD1*(dDi(i)*dD1i(j)+dDi(j)*dD1i(i))
     &      -D*fD1*dD1ij(i,j)-D*fD1D1*dD1i(i)*dD1i(j)
    Arn2(j,i)=Arn2(i,j)
    end do
    END IF
    end do
C   TEMPERATURE DERIVATIVES
    IF (NTD.EQ.1) THEN
    ArT=-TOTN*g-dDdT*f
    ArTV=-TOTN*gv-dDdT*fV
    ArTT=-dDdT2*f
    do i=1,nc
    ArTn(i)=-g+(TOTN*AUX/T-dDdT*fB)*dBi(i)-f*dDiT(i)-dDdT*fD1*dD1i(i)
    end do
    END IF
    end
c


! p      = pressure    (bar)                        (output)
! FUGLOG = vector of log. of fugacities (x*phi*P)   (output)    INDIC < 5
! DLFUGT = t-derivative of FUGLOG (const. vol,n)    (output)    INDIC = 2 or 4
! DLFUGV = vol-derivative of FUGLOG (const temp,n)  (output)    INDIC < 5
! DLFUGX = comp-derivative of FUGLOG (const t & v)  (output)    INDIC > 2

P = TOTN*RT/V - ArV
DPDV = -ArV2-RT*TOTN/V**2


DPDT = -ArTV+TOTN*RGAS/V

DO 60 I = 1 , NC

IF(RN(I) == 0.0)GOTO 60

FUGLOG(I)=Arn(I)/RT + log(rn(I)) + log(RT/V)
DLFUGV(I)=-DPDN(I)/RT                    ! term DPDV/P is cancelled out
DPDN(I) = RT/V-ArVn(I)

IF(NTEMP == 0) GOTO 60

DLFUGT(I)=(ArTn(I)-Arn(I)/T)/RT+1.D0/T    ! term DPDT/P is cancelled out
! 60 CONTINUE
! 62 IF(NDER.LT.2) GOTO 64
DO 63 I=1,NC
    DO 61 K=I,NC
        DLFUGX(I,K)=Arn2(I,K)/RT        ! term 1/TOTN is cancelled out
        61        DLFUGX(K,I)=DLFUGX(I,K)
        DLFUGX(I,I)=DLFUGX(I,I)+1.0/rn(I)










