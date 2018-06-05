
SUBROUTINE HelmSRKPR(nc,ND,NT,rn,V,T,Ar,ArV,ArTV,ArV2,Arn,ArVn,ArTn,Arn2)

    IMPLICIT DOUBLE PRECISION (A-H,O-Z)

    PARAMETER (nco=20,RGAS=0.08314472d0) 

	dimension rn(nc),Arn(nc),ArVn(nc),ArTn(nc),Arn2(nc,nc)
	dimension dBi(nc),dBij(nc,nc)
	dimension dDi(nc),dDij(nc,nc),dDiT(nc)
	dimension aij(nc,nc),daijdT(nc,nc),daijdT2(nc,nc)

    DOUBLE PRECISION Kij(nco,nco)
	dimension ac(nco),b(nco),del1(nco),rm(nco)

	COMMON /COMPONENTS/ ac,b,del1,rm,Kij,NTdep
	COMMON /rule/ncomb

    TOTN = sum(rn)
    D1=del1(1)
    D2=(1-D1)/(1+D1)

	if(ncomb.lt.2)then
        call Bnder(nc,rn,Bmix,dBi,dBij)
        call DandTnder(NT,nc,T,rn,D,dDi,dDiT,dDij,dDdT,dDdT2)
	else

! 	call Bcubicnder(nc,rn,Bmix,dBi,dBij)
! 	call DCubicandTnder(NT,nc,T,rn,D,dDi,dDiT,dDij,dDdT,dDdT2)
	end if

!   The f's and g's used here are for Ar, not F (reduced Ar)					***********
!   This requires to multiply by R all g, f and its derivatives as defined by Mollerup ****

    f=log((V+D1*Bmix)/(V+D2*Bmix))/Bmix/(D1-D2)
    g=RGAS*log(1-Bmix/V)
    fv=-1/((V+D1*Bmix)*(V+D2*Bmix))
    fB=-(f+V*fv)/Bmix
    gv=RGAS*Bmix/(V*(V-Bmix))
    fv2=(-1/(V+D1*Bmix)**2+1/(V+D2*Bmix)**2)/Bmix/(D1-D2)
    gv2=RGAS*(1/V**2-1/(V-Bmix)**2)

!   Reduced Helmholtz Energy and derivatives
    Ar=-TOTN*g*T-D*f
    ArV=-TOTN*gv*T-D*fv
    ArV2=-TOTN*gv2*T-D*fv2

    AUX=RGAS*T/(V-Bmix)
    FFB=TOTN*AUX-D*fB
    FFBV=-TOTN*AUX/(V-Bmix)+D*(2*fv+V*fv2)/Bmix
    FFBB=TOTN*AUX/(V-Bmix)-D*(2*f+4*V*fv+V**2*fv2)/Bmix**2

	do i=1,nc
        Arn(i)=-g*T+FFB*dBi(i)-f*dDi(i)
        ArVn(i)=-gv*T+FFBV*dBi(i)-fv*dDi(i)
	    IF (ND.EQ.2) THEN
	        do j=1,i
                Arn2(i,j)=AUX*(dBi(i)+dBi(j))-fB*(dBi(i)*dDi(j)+dBi(j)*dDi(i))  &
                          +FFB*dBij(i,j)+FFBB*dBi(i)*dBi(j)-f*dDij(i,j)      
                Arn2(j,i)=Arn2(i,j)
	        end do
	    END IF
	end do

! TEMPERATURE DERIVATIVES
	IF (NT.EQ.1) THEN
        ArT=-TOTN*g-dDdT*f
        ArTV=-TOTN*gv-dDdT*fV
        ArTT=-dDdT2*f

	    do i=1,nc
            ArTn(i)=-g+(TOTN*AUX/T-dDdT*fB)*dBi(i)-f*dDiT(i)
	    end do

	END IF

	end

 

 
# SUBROUTINE HelmSRKPR(nc,ND,NT,rn,V,T,Ar,ArV,ArTV,ArV2,Arn,ArVn,ArTn,Arn2)

#    IMPLICIT DOUBLE PRECISION (A-H,O-Z)

#    PARAMETER (nco=20,RGAS=0.08314472d0) 

# 	dimension rn(nc),Arn(nc),ArVn(nc),ArTn(nc),Arn2(nc,nc)
# 	dimension dBi(nc),dBij(nc,nc)
# 	dimension dDi(nc),dDij(nc,nc),dDiT(nc)
# 	dimension aij(nc,nc),daijdT(nc,nc),daijdT2(nc,nc)

#     DOUBLE PRECISION Kij(nco,nco)
# 	dimension ac(nco),b(nco),del1(nco),rm(nco)

# 	COMMON /COMPONENTS/ ac,b,del1,rm,Kij,NTdep
# 	COMMON /rule/ncomb

#     TOTN = sum(rn)
#     D1=del1(1)
#     D2=(1-D1)/(1+D1)

# 	if(ncomb.lt.2)then
#         call Bnder(nc,rn,Bmix,dBi,dBij)
#         call DandTnder(NT,nc,T,rn,D,dDi,dDiT,dDij,dDdT,dDdT2)
# 	else

# ! 	call Bcubicnder(nc,rn,Bmix,dBi,dBij)
# ! 	call DCubicandTnder(NT,nc,T,rn,D,dDi,dDiT,dDij,dDdT,dDdT2)
# 	end if

# !   The f's and g's used here are for Ar, not F (reduced Ar)					***********
# !   This requires to multiply by R all g, f and its derivatives as defined by Mollerup ****

#     f=log((V+D1*Bmix)/(V+D2*Bmix))/Bmix/(D1-D2)
#     g=RGAS*log(1-Bmix/V)
#     fv=-1/((V+D1*Bmix)*(V+D2*Bmix))
#     fB=-(f+V*fv)/Bmix
#     gv=RGAS*Bmix/(V*(V-Bmix))
#     fv2=(-1/(V+D1*Bmix)**2+1/(V+D2*Bmix)**2)/Bmix/(D1-D2)
#     gv2=RGAS*(1/V**2-1/(V-Bmix)**2)

# !   Reduced Helmholtz Energy and derivatives
#     Ar=-TOTN*g*T-D*f
#     ArV=-TOTN*gv*T-D*fv
#     ArV2=-TOTN*gv2*T-D*fv2

#     AUX=RGAS*T/(V-Bmix)
#     FFB=TOTN*AUX-D*fB
#     FFBV=-TOTN*AUX/(V-Bmix)+D*(2*fv+V*fv2)/Bmix
#     FFBB=TOTN*AUX/(V-Bmix)-D*(2*f+4*V*fv+V**2*fv2)/Bmix**2

# 	do i=1,nc
#         Arn(i)=-g*T+FFB*dBi(i)-f*dDi(i)
#         ArVn(i)=-gv*T+FFBV*dBi(i)-fv*dDi(i)
# 	    IF (ND.EQ.2) THEN
# 	        do j=1,i
#                 Arn2(i,j)=AUX*(dBi(i)+dBi(j))-fB*(dBi(i)*dDi(j)+dBi(j)*dDi(i))  &
#                           +FFB*dBij(i,j)+FFBB*dBi(i)*dBi(j)-f*dDij(i,j)      
#                 Arn2(j,i)=Arn2(i,j)
# 	        end do
# 	    END IF
# 	end do

# ! TEMPERATURE DERIVATIVES
# 	IF (NT.EQ.1) THEN
#         ArT=-TOTN*g-dDdT*f
#         ArTV=-TOTN*gv-dDdT*fV
#         ArTT=-dDdT2*f

# 	    do i=1,nc
#             ArTn(i)=-g+(TOTN*AUX/T-dDdT*fB)*dBi(i)-f*dDiT(i)
# 	    end do

# 	END IF

# 	end

