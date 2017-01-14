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

!  Comparison to test and debug cubic mixing rules
!  rn=[0.65,0.35]
!  T=460.0d0
!  	call Bnder(nc,rn,Bmix,dBi,dBij)
!  	call Bcubicnder(nc,rn,Bmix,dBi,dBij)
!  	call DandTnder(NTD,nc,T,rn,D,dDi,dDiT,dDij,dDdT,dDdT2)
!  	call DCubicandTnder(NTD,nc,T,rn,D,dDi,dDiT,dDij,dDdT,dDdT2)

	if(ncomb.lt.2)then
		call Bnder(nc,rn,Bmix,dBi,dBij)
		call DandTnder(NTD,nc,T,rn,D,dDi,dDiT,dDij,dDdT,dDdT2)
	else
!  	call Bcubicnder(nc,rn,Bmix,dBi,dBij)
!  	call DCubicandTnder(NTD,nc,T,rn,D,dDi,dDiT,dDij,dDdT,dDdT2)
	end if
!  The f's and g's used here are for Ar, not F (reduced Ar)					***********
!  This requires to multiply by R all g, f and its derivatives as defined by Mollerup ****
	f=log((V+D1*Bmix)/(V+D2*Bmix))/Bmix/(D1-D2)
	g=RGAS*log(1-Bmix/V)
	fv=-1/((V+D1*Bmix)*(V+D2*Bmix))
	fB=-(f+V*fv)/Bmix
	gv=RGAS*Bmix/(V*(V-Bmix))
	fv2=(-1/(V+D1*Bmix)**2+1/(V+D2*Bmix)**2)/Bmix/(D1-D2)
	gv2=RGAS*(1/V**2-1/(V-Bmix)**2)
!  DERIVATIVES OF f WITH RESPECT TO DELTA1
	auxD2=(1+2/(1+D1)**2)
	fD1=(1/(V+D1*Bmix)+2/(V+D2*Bmix)/(1+D1)**2)-f*auxD2
	fD1=fD1/(D1-D2)
	fBD1=-(fB*auxD2+D1/(V+D1*Bmix)**2+2*D2/(V+D2*Bmix)**2/(1+D1)**2)
	fBD1=fBD1/(D1-D2)
	fVD1=-(fV*auxD2+1/(V+D1*Bmix)**2+2/(V+D2*Bmix)**2/(1+D1)**2)/(D1-D2)
	fD1D1=4*(f-1/(V+D2*Bmix))/(1+D1)**3+Bmix*(-1/(V+D1*Bmix)**2+ &
			4/(V+D2*Bmix)**2/(1+D1)**4)-2*fD1*(1+2/(1+D1)**2)
	fD1D1=fD1D1/(D1-D2)
!  Reduced Helmholtz Energy and derivatives
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
!  TEMPERATURE DERIVATIVES
	IF (NTD.EQ.1) THEN
	ArT=-TOTN*g-dDdT*f
	ArTV=-TOTN*gv-dDdT*fV
	ArTT=-dDdT2*f
	do i=1,nc
	ArTn(i)=-g+(TOTN*AUX/T-dDdT*fB)*dBi(i)-f*dDiT(i)-dDdT*fD1*dD1i(i)
	end do
	END IF
	end













SUBROUTINE HelmRKPR(nco,NDE,NTD,rn,V,T,Ar,ArV,ArTV,ArV2,Arn,ArVn,ArTn,Arn2)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    PARAMETER (RGAS=0.08314472d0) 

	dimension rn(nco),Arn(nco),ArVn(nco),ArTn(nco),Arn2(nco,nco)
	dimension dBi(nco),dBij(nco,nco),dD1i(nco),dD1ij(nco,nco)
	dimension dDi(nco),dDij(nco,nco),dDiT(nco)
	dimension aij(nco,nco),daijdT(nco,nco),daijdT2(nco,nco)

	COMMON /rule/ncomb

	nc=nco
	TOTN = sum(rn)
	call DELTAnder(nc,rn,D1,dD1i,dD1ij)
	D2=(1-D1)/(1+D1)


	if(ncomb.lt.2)then

		call Bnder(nc,rn,Bmix,dBi,dBij)
		call DandTnder(NTD,nc,T,rn,D,dDi,dDiT,dDij,dDdT,dDdT2)

	else

!  	call Bcubicnder(nc,rn,Bmix,dBi,dBij)
!  	call DCubicandTnder(NTD,nc,T,rn,D,dDi,dDiT,dDij,dDdT,dDdT2)

	end if

!  The f's and g's used here are for Ar, not F (reduced Ar)					***********
!  This requires to multiply by R all g, f and its derivatives as defined by Mollerup ****

	f=log((V+D1*Bmix)/(V+D2*Bmix))/Bmix/(D1-D2)
	g=RGAS*log(1-Bmix/V)

	fv=-1/((V+D1*Bmix)*(V+D2*Bmix))
	fB=-(f+V*fv)/Bmix
	gv=RGAS*Bmix/(V*(V-Bmix))

	fv2=(-1/(V+D1*Bmix)**2+1/(V+D2*Bmix)**2)/Bmix/(D1-D2)
	gv2=RGAS*(1/V**2-1/(V-Bmix)**2)

!  DERIVATIVES OF f WITH RESPECT TO DELTA1
	auxD2=(1+2/(1+D1)**2)

	fD1=(1/(V+D1*Bmix)+2/(V+D2*Bmix)/(1+D1)**2)-f*auxD2
	fD1=fD1/(D1-D2)

	fBD1=-(fB*auxD2+D1/(V+D1*Bmix)**2+2*D2/(V+D2*Bmix)**2/(1+D1)**2)
	fBD1=fBD1/(D1-D2)

	fVD1=-(fV*auxD2+1/(V+D1*Bmix)**2+2/(V+D2*Bmix)**2/(1+D1)**2)/(D1-D2)

	fD1D1=4*(f-1/(V+D2*Bmix))/(1+D1)**3+Bmix*(-1/(V+D1*Bmix)**2+ &
			4/(V+D2*Bmix)**2/(1+D1)**4)-2*fD1*(1+2/(1+D1)**2)
	fD1D1=fD1D1/(D1-D2)

!  Reduced Helmholtz Energy and derivatives
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

!  TEMPERATURE DERIVATIVES
	IF (NTD.EQ.1) THEN
		ArT=-TOTN*g-dDdT*f
		ArTV=-TOTN*gv-dDdT*fV
		ArTT=-dDdT2*f
	
		do i=1,nc
			ArTn(i)=-g+(TOTN*AUX/T-dDdT*fB)*dBi(i)-f*dDiT(i)-dDdT*fD1*dD1i(i)
		end do
	
	END IF
	
	end

