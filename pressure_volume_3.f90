    NTEMP=0
    NDER=1
    XOLD=0.0D0
    TOLF=1.0D-3    ! for residuals  1.0D-3
    TOL= 1.0D-4    ! for variables  1.0D-4
    N=3    ! DLSARG CONSTANTS
    INCX=1
    LDA=3
    ldb=3

    TCmod=TC

    if(TCGC(1).ne.0.0D0)TCmod=TCGC
    inf=2

    if(icomp.eq.2)inf=1

    rn (inf) =0.0
    rn(icomp)=1.0

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

 2  XVAR = log([T, Vl, Vv])
    DFDS=0.0D0
    DFDS(3)=1.0D0
    RJAC=0.0D0
    RJAC(3,NS)=1.0D0
 1  NITER=0
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
 21     CALL XTVTERMO(2,T,Vl,Pl,rn,FUGx,FUGTx,FUGVx,DFGN)
        
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
        
        F(2)=FUGx(icomp)-FUGy(icomp)
        RJAC(1,1)=DPDTx/Pl-DPDTy/Pv
        RJAC(1,1)=T*RJAC(1,1)
        RJAC(1,2)=Vl*DPDVx/Pl
        RJAC(1,3)=-Vv*DPDVy/Pv
C
        RJAC(2,1)=FUGTx(icomp)-FUGTy(icomp)
        RJAC(2,1)=T*RJAC(2,1)
        RJAC(2,2)=Vl*FUGVx(icomp)
        RJAC(2,3)=-Vv*FUGVy(icomp)
C
        dold=delx

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

    if(delXS.eq.0.015.and.maxval(abs(delX)).gt.7)delX=6*delX/abs(delX(3))

 17     XVAR=XVAR+delX
 18     T=exp(XVAR(1))
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

 7  S = XOLD(NS) + delXS
    XVAR = XOLD+dXdS*delXS/dXdS(NS)    ! Initial estimates for the 3 variables in the next point
 
 8  T=exp(XVAR(1))
    Vl=exp(XVAR(2))
    Vv=exp(XVAR(3))
 
    if(Pv.gt.1.D-20.and.(T.gt.0.25*TCmod(icomp).or.T>250)     ! modified April 6, 2016
    1                              .and.T.lt.TCmod(icomp))go to 1
 



#
      SUBROUTINE TERMO(MTYP,INDIC,IC,T,P,rn,V,PHILOG,DLPHIT,DLPHIP,FUGN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXC=2,nco=2,RGAS=0.08314472d0)
      DIMENSION FUGN(MAXC,MAXC)
      DIMENSION PHILOG(MAXC),DLPHIT(MAXC),DLPHIP(MAXC),DPDN(MAXC)
    dimension rn(nco),Arn(nco),ArVn(nco),ArTn(nco),Arn2(nco,nco)
c   The output PHILOG is actually the vector ln(phi(i)*P)
    NC=2
    NTEMP=0
      IGZ=0
      NDER=1
      IF (INDIC.GT.2) NDER=2
      IF (INDIC.EQ.2 .OR. INDIC.EQ.4) NTEMP=1
    TOTN = sum(rn)
    if(P.le.0.0d0)MTYP=1
    CALL VCALC(MTYP,NC,rn,T,P,V)      
      RT = RGAS*T
      Z = V/(TOTN*RT)   ! this is Z/P
    call ArVnder(NDER,NTEMP,rn,V,T,Ar,ArV,ArTV,ArV2,Arn,ArVn,ArTn,Arn2)
      DPV = -ArV2-RT*TOTN/V**2
    DPDT = -ArTV+TOTN*RGAS/V
      DO 60 I=1,NC
      PHILOG(I)=-LOG(Z)+Arn(I)/RT
      DPDN(I) = RT/V-ArVn(I)
      DLPHIP(I)=-DPDN(I)/DPV/RT-1.D0/P
    IF(NTEMP.EQ.0) GOTO 60
    DLPHIT(I)=(ArTn(I)-Arn(I)/T)/RT+DPDN(I)*DPDT/DPV/RT+1.D0/T
   60 CONTINUE
   62 IF(NDER.LT.2) GOTO 64
      DO 63 I=1,NC
      DO 61 K=I,NC
      FUGN(I,K)=1.D0/TOTN+(Arn2(I,K)+DPDN(I)*DPDN(K)/DPV)/RT 
   61 FUGN(K,I)=FUGN(I,K)
   63 CONTINUE
   64 RETURN
      END


















