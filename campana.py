import numpy as np

#import solution_matrix as sm

#c    Starting from critical point    c    c
#    WRITE(nout,15) TCmod(icomp),Pc(icomp),Dc(icomp),Dc(icomp)
#    T = 0.9999 * TCmod(icomp)
#    Vc = 1/DC(icomp)
#    Vv = 1.03*Vc
#    Zc = Pc(icomp)*Vc/RGAS/T
#    aaa = min(0.89 + (Zc - 0.2) / 2, 0.95)
#    Vl = aaa*Vc
#    NS = 3
#    delXS = 0.10

# PARAMETER (nco=2, RGAS=0.08314472d0, NA=2)
nco = 2
RGAS = 0.08314472
RGAS = 0.80324234
NA = 2
#----------------------------------------------------------
Tc = 123
DC = 0.234
Pc = 23
#----------------------------------------------------------
temperature = 0.9999 * Tc
critical_volume = 1 / DC
vapor_volume = 1.03 * critical_volume
Zc = Pc * critical_volume / RGAS / temperature
liquid_volume = min(0.89 + (Zc - 0.2) / 2, 0.95) * critical_volume
NS = 2
delXS = 0.10
#------------------------------------------------------
XVAR = np.log([temperature, liquid_volume, vapor_volume])

#DFDS = 0.0
DFDS = np.zeros(3)
DFDS[2] = 1.0
print("DFDS = ", format(DFDS))

#RJAC=0.0
RJAC = np.zeros([3, 3])
RJAC[2, NS]=1.0
print("RJAC = ", format(RJAC))

NITER=0
T = np.exp(XVAR[0])
Vl = np.exp(XVAR[1])
Vv = np.exp(XVAR[2])


print("T = {0}, Vl = {1}, Vv = {2}".format(T, Vl, Vv))

FMAXOLD = 8.0
FMAX = 7.0
DMAXOLD = 8.0
DMAX = 7.0
F = np.zeros(3)

F[2] = 0.0
delX = 0.0
NV = 0
#------------------------------------------------------


#if(NV.GE.5.and.T.LT.0.4*Tcmod(icomp))return
#      IF ((FMAX.GT.FMAXOLD.and.DMAX.GT.DMAXOLD).OR.NITER.eq.10) THEN
#            delXS=0.8*delXS
#            if(abs(delXS).lt.0.001)return
#            if(XOLD(1).NE.0.0D0)go to 7
#        END IF
#        NITER=NITER+1



criterio_1 = NV > 0.5 and T < 0.4*Tcmod

# IF ((FMAX.GT.FMAXOLD.and.DMAX.GT.DMAXOLD).OR.NITER.eq.10) THEN
if FMAX > FMAXOLD and DMAX > DMAXOLD or NITER == 10:
	# delXS=0.8*delXS
	delXS = 0.8 * delXS
	# if(abs(delXS).lt.0.001)return
	if abs(delXS) < 0.001:
		return
		# if(XOLD(1).NE.0.0D0)go to 7
		if XOLD[0] != 0:
			pass

#---------------------------------------------------------------------------

#    NSOLD=NS

#    IF(VV.GT.1.1*VC)THEN  ! .and.T.GT.0.3*Tcmod(icomp)
#        NS=1
#    ELSE
#        NS=3
#    END IF

NSOLD = NS

if VV > 1.1 * VC:
	NS = 1
else:
	NS = 3

#----------------------------------------------------------------------------

#    if(NS.ne.NSOLD)then
#        dSdSold=dXdS(NS)
#        dXdS=dXdS/dSdSold
#        RJAC(3,1:3)=0.0D0
#        RJAC(3,NS)=1.0D0
#    end if

#    NITER=min(niter,10)
#    delXS=delXS*5/NITER/dXdS(NSOLD)



if NS != NSOLD:
	dSdSold = dXdS[NS]
	dXdS = dXdS / dSdSold
	RJAC[3, 1:3] = 0.0
	RJAC[3, NS] = 1.0

NITER = min(niter,10)
delXS = delXS * 5 / NITER / dXdS[NSOLD]

#------------------------------------------------------------------------

#IF(delXS.LT.0)then
#        if(NS.EQ.1)
#            delXS=max(delXS,-0.008) ! Max lnT decrease allowed
#        if(NS.eq.2)
#            delXS=max(delXS,-0.01) ! Max lnVl decrease allowed: 0.01  ! DO NOT CHANGE !!!
#        if(NS.eq.3)
#            delXS=max(delXS,-0.10) ! Max lnVv decrease allowed: 0.20
#    ELSE
#        if(NS.EQ.1)
#            delXS=min(delXS,0.01) ! Max lnT increase allowed: 0.5 K
#        if(NS.eq.2)
#            delXS=min(delXS,0.05) ! Max lnVl increase allowed: 0.01
#        if(NS.eq.3)
#            delXS=min(delXS,0.20) ! Max lnVv increase allowed: 0.20
#    END IF



def comprobar_delXS(delXS, NS):

	if delXS < 0.0:
		if NS == 1:
			#! Max lnT decrease allowed
			delXS=max(delXS,-0.008) 
		elif NS == 2:
			#! Max lnVl decrease allowed: 0.01  ! DO NOT CHANGE !!!
			delXS=max(delXS,-0.01) 
		elif NS == 3:
			#! Max lnVv decrease allowed: 0.20
			delXS=max(delXS,-0.10) 
	else:
		if NS == 1:
			#! Max lnT increase allowed: 0.5 K
			delXS=min(delXS,0.01)
		if NS == 2:
			#! Max lnVl increase allowed: 0.01
			delXS=min(delXS,0.05)
		if NS == 3:
			#! Max lnVv increase allowed: 0.20
			delXS=min(delXS,0.20)

	return delXS


 
XOLD = XVAR
TOLD = T
POLD = Pv
DVOLD = RHOV
DLOLD = RHOL
NV = 0

S = XOLD(NS) + delXS

#! Initial estimates for the 3 variables in the next point
XVAR = XOLD + dXdS * delXS / dXdS[NS]

T = np.exp(XVAR[1])
Vl = np.exp(XVAR[2])
Vv = np.exp(XVAR[3])


#if(Pv.gt.1.D-20.and.(T.gt.0.25*TCmod(icomp).or.T>250)     ! modified April 6, 2016
#    1                              .and.T.lt.TCmod(icomp))go to 1
 
if (Pv > 1e-20 and (T > 0.25 * TCmod or T > 250) and T < TCmod):
	print("go to 1")
 

#      SUBROUTINE XTVTERMO(INDIC,T,V,P,rn,
#    1                    FUGLOG,DLFUGT,DLFUGV,DLFUGX)
#C
#C-------parameters of XTVTERMO (crit. point, LLV and CEP calculations)
#C
#C       rn        mixture mole numbers                     (input)
#C       t            temperature (k)                       (input)
#C       v            volume        (L)                     (input)
#C       p            pressure    (bar)                     (output)
#C       FUGLOG    vector of log. of fugacities (x*phi*P)   (output)    INDIC < 5
#C       DLFUGT    t-derivative of FUGLOG (const. vol,n)    (output)    INDIC = 2 or 4
#C       DLFUGV    vol-derivative of FUGLOG (const temp,n)  (output)    INDIC < 5
#C       DLFUGX    comp-derivative of FUGLOG (const t & v)  (output)    INDIC > 2
#C---------------------------------------------------
#C---  MODIFIED AND CORRECTED july 2005
#C---
#C---------------------------------------------------
#      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#      PARAMETER (MAXC=2,nco=2,RGAS=0.08314472d0)
#      DIMENSION DLFUGX(MAXC,MAXC)
#      DIMENSION FUGLOG(MAXC),DLFUGT(MAXC),DLFUGV(MAXC)
#    dimension rn(nco),Arn(nco),ArVn(nco),ArTn(nco),Arn2(nco,nco)
#    COMMON /MODEL/ NMODEL
#    COMMON/NG/NGR
#    COMMON /Pder/ DPDN(nco),DPDT,DPDV
#    NG=NGR
#    NC=2
#    IF(NMODEL.EQ.5.OR.NMODEL.EQ.7) CALL PARAGC(T,NC,NG,1)      
#    NTEMP=0
#      IGZ=0
#      NDER=1
#      IF (INDIC.GT.2) NDER=2
#      IF (INDIC.EQ.2 .OR. INDIC.EQ.4) NTEMP=1
#    TOTN = sum(rn)
#      RT = RGAS*T
#    call ArVnder(NDER,NTEMP,rn,V,T,Ar,ArV,ArTV,ArV2,Arn,ArVn,ArTn,Arn2)
#      P = TOTN*RT/V - ArV
#      DPDV = -ArV2-RT*TOTN/V**2
#      IF(INDIC.GT.4)GOTO 62
#c      Z = P*V/(TOTN*RT)
#    DPDT = -ArTV+TOTN*RGAS/V
#    DO 60 I=1,NC
#    IF(RN(I).EQ.0.0)GOTO 60
#C        FUGLOG(I)=-LOG(Z)+Arn(I)/RT + log(rn(I)/TOTN) + log(P)
#C        FUGLOG(I)=Arn(I)/RT + log(rn(I)/TOTN) + log(P/Z) this crashes at very low T LLV when Z=P=0.000000...
#        FUGLOG(I)=Arn(I)/RT + log(rn(I)) + log(RT/V)
#        DPDN(I) = RT/V-ArVn(I)
#        DLFUGV(I)=-DPDN(I)/RT                    ! term DPDV/P is cancelled out
#        IF(NTEMP.EQ.0) GOTO 60
#        DLFUGT(I)=(ArTn(I)-Arn(I)/T)/RT+1.D0/T    ! term DPDT/P is cancelled out
#   60 CONTINUE
#   62 IF(NDER.LT.2) GOTO 64
#      DO 63 I=1,NC
#      DO 61 K=I,NC
#        DLFUGX(I,K)=Arn2(I,K)/RT        ! term 1/TOTN is cancelled out
#   61        DLFUGX(K,I)=DLFUGX(I,K)
#        DLFUGX(I,I)=DLFUGX(I,I)+1.0/rn(I)
#   63 CONTINUE
#   64 RETURN
#      END



def XTVTERMO_cal(INDIC, T, V, P, rn):

	NG = NGR
	NC = 2
	NTEMP = 0
	IGZ = 0
	NDER = 1

	if INDIC > 2:
		NDER = 2

	if INDIC > 2 or INDIC == 4:
		NTEMP = 1

	TOTN = sum(rn)
	RT = RGAS*T

	#call ArVnder(NDER, NTEMP, rn, V, T, Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2)
	ArVnder(NDER, NTEMP, rn, V, T, Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2)

	P = TOTN*RT/V - ArV
	DPDV = -ArV2 - RT * TOTN / V ** 2

	if INDIC > 4:
		print("GOTO 62")

	Z = P * V / (TOTN * RT)
	DPDT = -ArTV + TOTN * RGAS / V

	#DO 60 I=1,NC

	#IF(RN(I).EQ.0.0)GOTO 60
	if rn[I] == 0:
		print("GOTO 60")

	FUGLOG[I] = Arn[I] / RT + log(rn[I]) + log(RT / V)
	DPDN[I] = RT / V - ArVn[I]

	#! term DPDV/P is cancelled out
	DLFUGV[I] = -DPDN[I] / RT                    

	if NTEMP == 0:
		print("GOTO 60")

	#! term DPDT/P is cancelled out
	DLFUGT[I] = (ArTn[I] - Arn[I] / T) / RT + 1.0 / T

	if NDER < 2:
		#DO 63 I=1,NC
		for i in range(NC):
			pass
		#DO 61 K=I,NC
		for k in range(NC):
			pass

	#! term 1/TOTN is cancelled out
	DLFUGX[I, K] = Arn2[I, K] / RT
	DLFUGX[K, I] = DLFUGX[I, K]
	DLFUGX[I, I] = DLFUGX[I, I] + 1.0 / rn[I]

	return FUGLOG, DLFUGT, DLFUGV, DLFUGX


def SUBROUTINE ArVnder(nc,NDER,NTD,rn,V,T,Ar,ArV,ArTV,ArV2,Arn,ArVn,ArTn,Arn2):
    
    #IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	#dimension rn(nc),Arn(nc),ArVn(nc),ArTn(nc),Arn2(nc,nc)
	#COMMON /MODEL/ NMODEL
	
	#IF(NMODEL.LE.2)THEN
	if NMODEL <= 2:
		#!  SRK or PR
		#CALL HelmSRKPR(nc,NDER,NTD,rn,V,T,Ar,ArV,ArTV,ArV2,Arn,ArVn,ArTn,Arn2)
		HelmSRKPR(nc,NDER,NTD,rn,V,T,Ar,ArV,ArTV,ArV2,Arn,ArVn,ArTn,Arn2)
	#ELSE IF (NMODEL.EQ.3) THEN
	elif NMODEL == 3:
		#CALL HelmRKPR(nc,NDER,NTD,rn,V,T,Ar,ArV,ArTV,ArV2,Arn,ArVn,ArTn,Arn2)
		HelmRKPR(nc,NDER,NTD,rn,V,T,Ar,ArV,ArTV,ArV2,Arn,ArVn,ArTn,Arn2)
	#END





def SUBROUTINE HelmRKPR():
	#---------------------------------------

	#HelmRKPR

	#	SUBROUTINE HelmRKPR(nco,NDE,NTD,rn,V,T,Ar,ArV,ArTV,ArV2,Arn,ArVn,ArTn,Arn2)
	#    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	#    PARAMETER (RGAS=0.08314472d0) 
	#	dimension rn(nco),Arn(nco),ArVn(nco),ArTn(nco),Arn2(nco,nco)
	#	dimension dBi(nco),dBij(nco,nco),dD1i(nco),dD1ij(nco,nco)
	#	dimension dDi(nco),dDij(nco,nco),dDiT(nco)
		dimension aij(nco,nco),daijdT(nco,nco),daijdT2(nco,nco)
		COMMON /rule/ncomb
		nc=nco
		TOTN = sum(rn)
		call DELTAnder(nc,rn,D1,dD1i,dD1ij)
		D2=(1-D1)/(1+D1)

	#!  Comparison to test and debug cubic mixing rules
	#!  rn=[0.65,0.35]
	#!  T=460.0d0
	#!  	call Bnder(nc,rn,Bmix,dBi,dBij)
	#!  	call Bcubicnder(nc,rn,Bmix,dBi,dBij)
	#!  	call DandTnder(NTD,nc,T,rn,D,dDi,dDiT,dDij,dDdT,dDdT2)
	#!  	call DCubicandTnder(NTD,nc,T,rn,D,dDi,dDiT,dDij,dDdT,dDdT2)

		if(ncomb.lt.2)then
			call Bnder(nc,rn,Bmix,dBi,dBij)
			call DandTnder(NTD,nc,T,rn,D,dDi,dDiT,dDij,dDdT,dDdT2)
		else
	#!  	call Bcubicnder(nc,rn,Bmix,dBi,dBij)
	#!  	call DCubicandTnder(NTD,nc,T,rn,D,dDi,dDiT,dDij,dDdT,dDdT2)
		end if
	#!  The f's and g's used here are for Ar, not F (reduced Ar)					***********
	#!  This requires to multiply by R all g, f and its derivatives as defined by Mollerup ****
		f=log((V+D1*Bmix)/(V+D2*Bmix))/Bmix/(D1-D2)
		g=RGAS*log(1-Bmix/V)
		fv=-1/((V+D1*Bmix)*(V+D2*Bmix))
		fB=-(f+V*fv)/Bmix
		gv=RGAS*Bmix/(V*(V-Bmix))
		fv2=(-1/(V+D1*Bmix)**2+1/(V+D2*Bmix)**2)/Bmix/(D1-D2)
		gv2=RGAS*(1/V**2-1/(V-Bmix)**2)
	#!  DERIVATIVES OF f WITH RESPECT TO DELTA1
		auxD2=(1+2/(1+D1)**2)
		fD1=(1/(V+D1*Bmix)+2/(V+D2*Bmix)/(1+D1)**2)-f*auxD2
		fD1=fD1/(D1-D2)
		fBD1=-(fB*auxD2+D1/(V+D1*Bmix)**2+2*D2/(V+D2*Bmix)**2/(1+D1)**2)
		fBD1=fBD1/(D1-D2)
		fVD1=-(fV*auxD2+1/(V+D1*Bmix)**2+2/(V+D2*Bmix)**2/(1+D1)**2)/(D1-D2)
		fD1D1=4*(f-1/(V+D2*Bmix))/(1+D1)**3+Bmix*(-1/(V+D1*Bmix)**2+ &
				4/(V+D2*Bmix)**2/(1+D1)**4)-2*fD1*(1+2/(1+D1)**2)
		fD1D1=fD1D1/(D1-D2)
	#!  Reduced Helmholtz Energy and derivatives
		Ar=-TOTN*g*T-D*f
		ArV=-TOTN*gv*T-D*fv
		ArV2=-TOTN*gv2*T-D*fv2

		AUX=RGAS*T/(V-Bmix)
		FFB=TOTN*AUX-D*fB
		FFBV=-TOTN*AUX/(V-Bmix)+D*(2*fv+V*fv2)/Bmix
		FFBB=TOTN*AUX/(V-Bmix)-D*(2*f+4*V*fv+V**2*fv2)/Bmix**2
		do i=1,nc
		Arn(i)=-g*T+FFB*dBi(i)-f*dDi(i)-D*fD1*dD1i(i)
		ArVn(i)=-gv*T+FFBV*dBi(i)-fv*dDi(i)-D*fVD1*dD1i(i)
		IF (NDE.EQ.2) THEN
		do j=1,i
		Arn2(i,j)=AUX*(dBi(i)+dBi(j))-fB*(dBi(i)*dDi(j)+dBi(j)*dDi(i))  &
	     		+FFB*dBij(i,j)+FFBB*dBi(i)*dBi(j)-f*dDij(i,j)      
	      Arn2(i,j)=Arn2(i,j)-D*fBD1*(dBi(i)*dD1i(j)+dBi(j)*dD1i(i))    &
	     		-fD1*(dDi(i)*dD1i(j)+dDi(j)*dD1i(i))        &
	     		-D*fD1*dD1ij(i,j)-D*fD1D1*dD1i(i)*dD1i(j)
		Arn2(j,i)=Arn2(i,j)
		end do
		END IF
		end do
	#!  TEMPERATURE DERIVATIVES
		IF (NTD.EQ.1) THEN
		ArT=-TOTN*g-dDdT*f
		ArTV=-TOTN*gv-dDdT*fV
		ArTT=-dDdT2*f
		do i=1,nc
		ArTn(i)=-g+(TOTN*AUX/T-dDdT*fB)*dBi(i)-f*dDiT(i)-dDdT*fD1*dD1i(i)
		end do
		END IF
		end
	#-----------------------------------------------------------------

###########################################################
	subroutine aTder(ac,Tc,rk,T,a,dadT,dadT2)
      implicit DOUBLE PRECISION (A-H,O-Z)
!  Given ac,Tc and the k parameter of the RKPR correlation, as well as the actual T,
!  this subroutine calculates a(T) and its first and second derivatives with T.
	COMMON /MODEL/ NMODEL
	Tr=T/Tc
	IF(NMODEL.LE.2)THEN
		rm=rk
		a=ac*(1+rm*(1-sqrt(Tr)))**2
		dadT=ac*rm*(rm-(rm+1)/sqrt(Tr))/Tc
		dadT2=ac*rm*(rm+1)/(2*Tc**2*Tr**1.5D0)
	ELSE
		a=ac*(3/(2+Tr))**rk
		dadT=-rk*a/Tc/(2+Tr)
		dadT2=-(rk+1)*dadT/Tc/(2+Tr)
	END IF
	end

	subroutine aijTder(NTD,nc,T,aij,daijdT,daijdT2)
      implicit DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (nco=20)
	DOUBLE PRECISION Kinf,Kij0(nco,nco),Kij(nco,nco),Tstar(nco,nco)
	dimension ai(nc),daidT(nc),daidT2(nc)
	dimension aij(nc,nc),daijdT(nc,nc),daijdT2(nc,nc)
	dimension aux(nc,nc),ratK(nc,nc)
    COMMON/CRIT/TC(nco),PC(nco),DCeos(nco),OM(nco)
	COMMON /COMPONENTS/ ac(nco),b(nco),d1(nco),rk(nco),Kij0,NTDEP
	COMMON /bcross/bij(nco,nco)
	COMMON /rule/ncomb
	COMMON /Tdep/ Kinf,Tstar
	IF(NTDEP.GE.1)THEN
		Kij=0.0D0
	    DO i=1,nc
    		Kij(:i-1, i)=Kinf+Kij0(:i-1, i)*exp(-T/Tstar(:i-1, i))
        END DO
!       Kij(2,1)=Kij(1,2)
!  ELSE IF(NTDEP.EQ.2)THEN
!  	Kij=0.0D0
!  	Kij(1,2)=Kij0(1,2)*exp(-T/Tstar)
!  	Kij(2,1)=Kij(1,2)
	ELSE
		Kij=Kij0
	END IF
	DO i=1,nc
	    call aTder(ac(i),Tc(i),rk(i),T,ai(i),daidT(i),daidT2(i))
	    aij(i,i)=ai(i)
	    daijdT(i,i)=daidT(i)
	    daijdT2(i,i)=daidT2(i)
	    IF (i.gt.1) THEN
	    do j=1,i-1
	        aij(j,i)=sqrt(ai(i)*ai(j))*(1-Kij(j,i))
	        aij(i,j)=aij(j,i)
	        if(NTD.EQ.1)then
		        daijdT(j,i)=(1-Kij(j,i))*(sqrt(ai(i)/ai(j))*daidT(j)+sqrt(ai(j)/ai(i))*daidT(i))/2
		        daijdT(i,j)=daijdT(j,i)
		        daijdT2(j,i)=(1-Kij(j,i))*(daidT(j)*daidT(i)/sqrt(ai(i)*ai(j))  &
     		        +sqrt(ai(i)/ai(j))*(daidT2(j)-daidT(j)**2/(2*ai(j)))        &
     		        +sqrt(ai(j)/ai(i))*(daidT2(i)-daidT(i)**2/(2*ai(i))))/2
		        daijdT2(i,j)=daijdT2(j,i)
	        end if
	    end do
	    END IF
	END DO
	if (ncomb.eq.1) then
		DO i=1,nc-1
		DO j=i+1,nc
			barrgij=bij(i,j)/sqrt(b(i)*b(j))
			aij(i,j)=barrgij*aij(i,j)
			aij(j,i)=aij(i,j)
			daijdT(i,j)=barrgij*daijdT(i,j)
			daijdT(j,i)=daijdT(i,j)
			daijdT2(i,j)=barrgij*daijdT2(i,j)
			daijdT2(j,i)=daijdT2(i,j)
		END DO
		END DO
	end if
!    		Kij(:i-1, i)=Kinf+Kij0(:i-1, i)*exp(-T/Tstar(:i-1, i))
    		
	IF(NTDEP.ge.1.and.NTD.EQ.1)THEN
	    DO i=1,nc
		    aux(:i-1, i)=daijdT(:i-1, i)
		    ratK(:i-1, i)=Kij(:i-1, i)/(1-Kij(:i-1, i))/Tstar(:i-1, i)
		    daijdT(:i-1, i)=aux(:i-1, i)+aij(:i-1, i)*ratK(:i-1, i)
		    daijdT(i, :i-1)=daijdT(:i-1, i)
		    daijdT2(:i-1, i)=daijdT2(:i-1, i)+(2*aux(:i-1, i)-aij(:i-1, i)/Tstar(:i-1, i))*ratK(:i-1, i)  
		    daijdT2(i, :i-1)=daijdT2(:i-1, i)
		END DO
	END IF
	end

	subroutine DandTnder(NTD,nco,T,rn,D,dDi,dDiT,dDij,dDdT,dDdT2)
      implicit DOUBLE PRECISION (A-H,O-Z)
!      PARAMETER (nco=2)
	dimension rn(nco),dDiT(nco)
	dimension dDi(nco),dDij(nco,nco)
	dimension aij(nco,nco),daijdT(nco,nco),daijdT2(nco,nco)
	call aijTder(NTD,nc,T,aij,daijdT,daijdT2)
	nc=nco
	D=0.0D0
	dDdT=0.0D0
	dDdT2=0.0D0
	DO i=1,nc
	aux=0.0D0
	aux2=0.0D0
	dDi(i)=0.0D0
	dDiT(i)=0.0D0
	do j=1,nc
	dDi(i)=dDi(i)+2*rn(j)*aij(i,j)
	if(NTD.EQ.1)then
		dDiT(i)=dDiT(i)+2*rn(j)*daijdT(i,j)
		aux2=aux2+rn(j)*daijdT2(i,j)
	end if
	dDij(i,j)=2*aij(i,j)
	aux=aux+rn(j)*aij(i,j)
	end do
	D=D+rn(i)*aux
	if(NTD.EQ.1)then
		dDdT=dDdT+rn(i)*dDiT(i)/2
		dDdT2=dDdT2+rn(i)*aux2
	end if
	END DO
	end

	subroutine DELTAnder(nc,rn,D1m,dD1i,dD1ij)
      implicit DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (nco=20)
	DOUBLE PRECISION Kij(nco,nco)
	dimension rn(nc),dD1i(nc),dD1ij(nc,nc)
	COMMON /COMPONENTS/ ac(nco),b(nco),d1(nco),rk(nco),Kij,NTDEP
	D1m=0.0D0
	DO i=1,nc
	D1m=D1m+rn(i)*d1(i)
	END DO
	TOTN = sum(rn)
	D1m=D1m/totn
	do i=1,nc
	dD1i(i)=(d1(i)-D1m)/totn
	do j=1,nc
	dD1ij(i,j)=(2.0D0*D1m-d1(i)-d1(j))/totn**2
	end do
	end do
	end

	subroutine Bnder(nc,rn,Bmix,dBi,dBij)
      implicit DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (nco=20)
	dimension rn(nc),dBi(nc),dBij(nc,nc),aux(nc)
	COMMON /bcross/bij(nco,nco)
	TOTN = sum(rn)
	Bmix=0.0D0
	aux=0.0D0
	DO i=1,nc
		do j=1,nc
			aux(i)=aux(i)+rn(j)*bij(i,j)
		end do
		Bmix=Bmix+rn(i)*aux(i)
	END DO
	Bmix=Bmix/totn
	DO i=1,nc
		dBi(i)=(2*aux(i)-Bmix)/totn
		do j=1,i
			dBij(i,j)=(2*bij(i,j)-dBi(i)-dBi(j))/totn
			dBij(j,i)=dBij(i,j)
		end do
	END DO
	end


###########################################################



def HelmRKPR(nco,NDE,NTD,rn,V,T,Ar,ArV,ArTV,ArV2,Arn,ArVn,ArTn,Arn2):
    # IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    # PARAMETER (RGAS=0.08314472d0) 

	# dimension rn(nco),Arn(nco),ArVn(nco),ArTn(nco),Arn2(nco,nco)
	# dimension dBi(nco),dBij(nco,nco),dD1i(nco),dD1ij(nco,nco)
	# dimension dDi(nco),dDij(nco,nco),dDiT(nco)
	# dimension aij(nco,nco),daijdT(nco,nco),daijdT2(nco,nco)

	# COMMON /rule/ncomb

	nc = nco
	TOTN = sum(rn)
	call DELTAnder(nc, rn, D1, dD1i, dD1ij)
	D2 = (1 - D1) / (1 + D1)


	#if(ncomb.lt.2)then
	if ncomb < 0.2:

		call Bnder(nc,rn,Bmix,dBi,dBij)
		call DandTnder(NTD,nc,T,rn,D,dDi,dDiT,dDij,dDdT,dDdT2)

	else:

#!  	call Bcubicnder(nc,rn,Bmix,dBi,dBij)
#!  	call DCubicandTnder(NTD,nc,T,rn,D,dDi,dDiT,dDij,dDdT,dDdT2)	

#!  The f's and g's used here are for Ar, not F (reduced Ar)					***********
#!  This requires to multiply by R all g, f and its derivatives as defined by Mollerup ****

	f = log((V + D1 * Bmix) / (V + D2 * Bmix)) / Bmix / (D1 - D2)
	g = RGAS * log(1 - Bmix / V)

	fv = -1/((V + D1 * Bmix) * (V + D2 * Bmix))
	fB = -(f + V * fv) / Bmix
	gv = RGAS * Bmix / (V * (V - Bmix))

	fv2=(-1/(V+D1*Bmix)**2+1/(V+D2*Bmix)**2)/Bmix/(D1-D2)
	gv2=RGAS*(1/V**2-1/(V-Bmix)**2)

#!  DERIVATIVES OF f WITH RESPECT TO DELTA1
	auxD2 = (1 + 2 / (1 + D1) ** 2)

	fD1 = (1 / (V + D1 * Bmix) + 2 / (V + D2 * Bmix) / (1 + D1)**2) - f * auxD2
	fD1 = fD1 / (D1 - D2)

	fBD1 = -(fB * auxD2 + D1 / (V + D1 * Bmix) ** 2 + 2 * D2 / (V + D2 * Bmix) ** 2 / (1 + D1) ** 2)
	fBD1 = fBD1 / (D1 - D2)

	fVD1 = -(fV * auxD2 + 1 / (V + D1 * Bmix) ** 2 + 2 / (V + D2 * Bmix) ** 2 / (1 + D1) ** 2) / (D1 - D2)

	fD1D1 = 4 * (f - 1 / (V + D2 * Bmix)) / (1 + D1) ** 3 + Bmix * (-1 / (V + D1 * Bmix) ** 2 + &
			4 / (V + D2 * Bmix)**2 / (1 + D1)**4) - 2 * fD1 * (1 + 2 / (1 + D1)**2)
	fD1D1 = fD1D1 / (D1 - D2)

#!  Reduced Helmholtz Energy and derivatives
	Ar = -TOTN * g * T - D * f
	ArV = -TOTN * gv * T - D * fv
	ArV2 = -TOTN * gv2 * T - D * fv2

	AUX = RGAS * T / (V - Bmix)
	FFB = TOTN * AUX - D * fB
	FFBV = -TOTN * AUX / (V - Bmix) + D * (2 * fv + V * fv2) / Bmix
	FFBB = TOTN * AUX / (V - Bmix) - D*(2 * f + 4 * V * fv + V**2 * fv2) / Bmix**2

	#do i=1,nc
	for i in range(nc):

		Arn(i) = -g * T + FFB * dBi[i] - f * dDi[i] - D * fD1 * dD1i[i]
		ArVn[i] = -gv * T + FFBV * dBi[i] - fv * dDi[i] - D * fVD1 * dD1i[i]

		#IF (NDE.EQ.2) THEN
		if NDE == 2:
			
			#do j=1,i
			for j in range(j):

				Arn2(i,j)=AUX*(dBi(i)+dBi(j))-fB*(dBi(i)*dDi(j)+dBi(j)*dDi(i))  &
     				+FFB*dBij(i,j)+FFBB*dBi(i)*dBi(j)-f*dDij(i,j)      
      			Arn2(i,j)=Arn2(i,j)-D*fBD1*(dBi(i)*dD1i(j)+dBi(j)*dD1i(i))    &
     				-fD1*(dDi(i)*dD1i(j)+dDi(j)*dD1i(i))        &
     				-D*fD1*dD1ij(i,j)-D*fD1D1*dD1i(i)*dD1i(j)
				Arn2(j,i)=Arn2(i,j)

		#	end do

		#END IF

	#end do

#!  TEMPERATURE DERIVATIVES
	#IF (NTD.EQ.1) THEN
	if NTD == 0.1:
		ArT = -TOTN * g - dDdT * f
		ArTV = -TOTN * gv - dDdT * fV
		ArTT = -dDdT2 * f
	
		#do i=1,nc
		for i in range(nc):
			ArTn[i] = -g + (TOTN * AUX / T - dDdT * fB) * dBi[i] - f * dDiT[i] - dDdT * fD1 * dD1i[i]
		#end do
	
	#END IF
	
	#end
	return

#-----------------------------------------------------------------


T = 123
TOLF = 1e-3
NS = 2
DFDS = np.zeros(3)
DFDS[2] = 1
F = np.ones(3)
F[2] = 0
RJAC = np.ones([3, 3])
RJAC[2, NS] = 0

delx_s = 0
NV = 0

Pl = 0.12
Pv = 3
FUGx = 2.3
FUGy = 3.7

condition_1 = Pl < 0.5 * Pv and (Pv - Pl) > 1e-8
condition_2 = Pl > 1.5 * Pv and (Pl - Pv) > 1e-8
print(Pl)
print(condition_1)
print(condition_2)

if condition_1:
	print('condition_1 = {0}'.format(condition_1))

if condition_2:
	print('condition_2 = {0}'.format(condition_2))

if condition_1 or condition_2:
	print(condition_1, condition_2)
	NV += 1
 
if Pl < 0:
	F[0] = TOLF
else:
	F[0] = np.log(Pl/Pv)

F[1] = FUGx - FUGy


#RJAC[0, 0] = T * (DPDTx / Pl - DPDTy / Pv)
#RJAC[0, 1] = Vl * DPDVx / Pl
#RJAC[0, 2] = -Vv * DPDVy / Pv

#RJAC[1, 0] = T * (FUGTx -FUGTy)
#RJAC[1, 1] = Vl * FUGVx
#RJAC[1, 2] = -Vv * FUGVy 

RJAC[0, 0] = 23
RJAC[0, 1] = 45
RJAC[0, 2] = -23

RJAC[1, 0] = 67
RJAC[1, 1] = 0.34
RJAC[1, 2] = -34.4 

#-------------------------------------------------------
print(F)
print(RJAC)
print(RJAC[2,2])

b = F
a = RJAC

print(a)

x = np.linalg.solve(a, b)
print(x)

print('sm.A = {0}'.format(sm.A))


class Component(object):		

	def function_Ar_cal(self):
		
		self.bv = self.B / self.V
		self.f = np.log((self.V + self.s1 * self.B) / (self.V + self.s2 * self.B)) / self.B / (self.s1 - self.s2)
		self.g = self.R * np.log(1 - self.B / self.V)

		self.AUX = self.R * self.T / (self.V - self.B)
		self.fB = -(self.f + self.V * self.fv) / self.B
		self.FFB = self.nT * AUX - self.D * self.fB
		self.Di = 2 * self.nT * self.ac * self.alfa
		self.Bi = self.bc

		self.Ar = -self.nT * self.g * self.T - self.D * self.f
		'''Primera derivada de F con respecto al volumen Ecu. (68)'''
		self.gv = self.R * self.B / (self.V * (self.V - self.B))
		self.fv = - 1 / ((self.V + self.s1 * self.B) * (self.V + self.s2 * self.B))
		self.ArV = -self.nT * self.gv * self.T - self.D * self.fv
		''' Segunda derivada de F con respecto al volumen Ecu. (74) '''
		self.gv2 = self.R * (1 / self.V ** 2 - 1 / (self.V - self.B) ** 2)
		self.fv2 = (- 1 / (self.V + self.s1 * self.B) ** 2 + 1 / (self.V + self.s2 * self.B) ** 2) / self.B / (self.s1 - self.s2)
		self.ArV2 = - self.nT * self.gv2 * self.T - self.D * self.fv2
		''' pressure '''
		self.Pcal = self.nT * self.R * self.T / self.V - self.ArV
		self.dPdV = -self.ArV2 - self.R * self.T * self.nT / self.V ** 2

		if self.eq != "RKPR":
			self.Arn = -self.g * self.T + self.FFB * self.Bi - self.f * self.Di
		else:
			self.Arn = -self.g * self.T + self.FFB * self.Bi - self.f * self.Di - self.D * self.fD1 * self.dD1i

		ArT = -nT * g - dDdT * f
		ArTV = -nT * gv - dDdT * fV
		ArTn = -g + (nT * AUX/T - dDdT * fB) * dBi - f * dDiT - dDdT * fD1 * dD1i
		ArVn = - gv * T + FFBV * dBi - fv * dDi - D * fVD1 * dD1i

		return self.g, self.f, self.Ar


component = Component()
component.function_Ar_cal()


