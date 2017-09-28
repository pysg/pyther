

def SUBROUTINE_HelmRKPR():
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
	#!  Given ac,Tc and the k parameter of the RKPR correlation, as well as the actual T,
	#!  this subroutine calculates a(T) and its first and second derivatives with T.
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
	#!       Kij(2,1)=Kij(1,2)
	#!  ELSE IF(NTDEP.EQ.2)THEN
	#!  	Kij=0.0D0
	#!  	Kij(1,2)=Kij0(1,2)*exp(-T/Tstar)
	#!  	Kij(2,1)=Kij(1,2)
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
	#!    		Kij(:i-1, i)=Kinf+Kij0(:i-1, i)*exp(-T/Tstar(:i-1, i))
    		
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
	#!      PARAMETER (nco=2)
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


SUBROUTINE_HelmRKPR()