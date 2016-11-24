! Listado de commons utilizados, que podrían pasar a un módulo:
!   COMMON /CRIT/TC(nco),PC(nco),DCeos(nco),omg(nco)
!	COMMON /COMPONENTS/ ac(nco),b(nco),delta1(nco),rk(nco),Kij_or_K0,NTDEP
!	COMMON /MODEL/ NMODEL
!	COMMON /rule/ncomb
!	COMMON /bcross/bij(nco,nco)
!	COMMON /Tdep/ Kinf,Tstar   (escalares por ahora, pero esto cambiará)


 program calc_flash

    implicit DOUBLE PRECISION (A-H,O-Z)
    OPEN(1,FILE='flashIN.txt')
    OPEN(2,FILE='flashOUT.txt')
    read(1,*) N
    call readcase(n)
END program

subroutine readcase(n)
    implicit DOUBLE PRECISION (A-H,O-Z)
    PARAMETER (nco=10)
	dimension z(n)

	DOUBLE PRECISION Kinf
        ! pure compound physical constants
        real*8, dimension(n) :: tcn
        real*8, dimension(n) :: pcn
        real*8, dimension(n) :: omgn

        ! eos parameters
        real*8, dimension(n) :: acn  ! in bar*(L/mol)**2
        real*8, dimension(n) :: bn   ! in L/mol
        real*8, dimension(n) :: delta1n  !only required for RKPR
        real*8, dimension(n) :: k_or_mn  ! k for RKPR ; m for SRK/PR

        ! interaction parameters matrices
        real*8, dimension(n,n) :: Kij_or_K0n, Lijn
!        real*8, dimension(n,n), intent(in) :: Tstarn   (scalar in this version)

    ! interaction parameters matrices
    real*8, dimension(nco,nco) :: Kij_or_K0, Lij  ! , Tstar

    real*8, dimension(nco) :: x , y  ! composition of liquid and vapour (molar fractions)

        COMMON /CRIT/TC(nco),PC(nco),DCeos(nco),omg(nco)
	    COMMON /COMPONENTS/ ac(nco),b(nco),delta1(nco),rk_or_m(nco),Kij_or_K0,NTDEP
	    COMMON /rule/ncomb
	    COMMON /bcross/bij(nco,nco)
	    COMMON /Tdep/ Kinf,Tstar    ! (escalares por ahora, pero esto cambiará)
	    COMMON /lforin/lij

        READ(1,*) (z(j),j=1,N)
        read(1,*) T
        read(1,*) P
        read(1,*) nmodel
      if(nmodel<3)then
        call read2PcubicNC(N,1,2)
      else if(nmodel==3)then
        call readRKPRNC(N,1,2)
      end if
      WRITE (2,*)
      write (2,4) (z(i),i=1,n)
 4	FORMAT('Molar fractions: ',10F6.3)
      WRITE (2,*)
      WRITE (2,*) ' T(K)=',T
      WRITE (2,*) ' P(bar)=',P
! Passing values from commons(nco) to input arguments (n)
        TCn = tc(:n)
        PCn = pc(:n)
        OMGn= omg(:n)
        acn = ac(:n)
        bn = b(:n)
        delta1n = delta1(:n)
        k_or_mn = rk_or_m(:n)
        Kij_or_K0n = Kij_or_K0(:n, :n)
        lijn = lij(:n, :n)
      call flash(nmodel, n, z, tcn, pcn, omgn, acn, bn, k_or_mn, delta1n, &
                     Kij_or_K0n, Tstarn, Lijn, t, p, x, y, rho_x, rho_y, beta)
      WRITE (2,*)
      WRITE (2,*) 'Beta (vapor phase fraction)= ',beta
      WRITE (2,*) 'Comp',(i,i=1,N)
      WRITE (2,1) (x(j),j=1,N)
      WRITE (2,2) (y(j),j=1,N)
 1	FORMAT('  x  ', 10E12.4)
 2	FORMAT('  y  ', 10E12.4)
 end subroutine readcase

subroutine flash(model, n, z, tcn, pcn, omgn, acn, bn, k_or_mn, delta1n, &
                     Kij_or_K0n, Tstarn, Lijn, t, p, x, y, rho_x, rho_y, beta)

    implicit DOUBLE PRECISION (A-H,O-Z)
    PARAMETER (nco=10)

        ! M&M means the book by Michelsen and Mollerup, 2nd Edition (2007)

        ! eos id and  number of compounds in the system
        integer, intent(in) :: model, n
	    DOUBLE PRECISION Kinf

        ! composition of the system
        real*8, dimension(n), intent(in) :: z

        ! pure compound physical constants
        real*8, dimension(n), intent(in) :: tcn
        real*8, dimension(n), intent(in) :: pcn
        real*8, dimension(n), intent(in) :: omgn

        ! eos parameters
        real*8, dimension(n), intent(in) :: acn  ! in bar*(L/mol)**2
        real*8, dimension(n), intent(in) :: bn   ! in L/mol
        real*8, dimension(n), intent(in) :: delta1n  !only required for RKPR
        real*8, dimension(n), intent(in) :: k_or_mn  ! k for RKPR ; m for SRK/PR

        ! interaction parameters matrices
        real*8, dimension(n,n), intent(in) :: Kij_or_K0n
!        real*8, dimension(n,n), intent(in) :: Tstarn   (scalar in this version)
        real*8, dimension(n,n), intent(in) :: Lijn

        ! Temperature and Pressure for the flash
        real*8, intent(in) :: t            ! Temperature for the flash (K)
        real*8, intent(in) :: p            ! Pressure for the flash (bar)

        ! Results from flash calculation
        real*8, dimension(n), intent(out) :: x  ! composition of liquid (molar fractions)
        real*8, dimension(n), intent(out) :: y  ! composition of vapour (molar fractions)
        real*8, intent(out) :: rho_x            ! density of liquid (moles/L)
        real*8, intent(out) :: rho_y            ! density of vapour (moles/L)
        real*8, intent(out) :: beta             ! total fraction of vapour (molar base)

        ! Intermediate variables during calculation process
        real*8, dimension(n) :: PHILOGy, PHILOGx
        real*8, dimension(n) :: KFACT, LOG_K, var_K, denom, DLPHIT, DLPHIP
        real*8, dimension(n, n) :: FUGN
        real*8 :: g0, g1  ! function g valuated at beta=0 and 1, based on Wilson K factors
        real*8 :: g, dg, bmin, bmax, Vy, Vx

        real*8, dimension(nco,nco) :: Kij_or_K0  ! , Tstar
        COMMON /CRIT/TC(nco),PC(nco),DCeos(nco),omg(nco)
	    COMMON /COMPONENTS/ ac(nco),b(nco),delta1(nco),rk_or_m(nco),Kij_or_K0,NTDEP
	    COMMON /MODEL/ NMODEL
	    COMMON /rule/ncomb
	    COMMON /bcross/bij(nco,nco)
	    COMMON /Tdep/ Kinf,Tstar

! Charging the commons(nco) from input arguments (n)
        NMODEL = model
        TC(:n) = tcn
        PC(:n) = pcn
        OMG(:n)= omgn
        ac(:n) = acn
        b(:n) = bn
        delta1(:n) = delta1n
        rk_or_m(:n) = k_or_mn
        Kij_or_K0(:n, :n) = Kij_or_K0n
        Kinf = 0.0d0
        ncomb = 0  ! only  vdW combining rules and quadratic mixing rules by  the moment
        Tstar = Tstarn
	! b matrix for Classical or van der Waals combining rules:
		do i=1,n
		    do j=i,n
		        bij(i,j)=(1-lijn(i,j))*(b(i)+b(j))/2
		        bij(j,i)=bij(i,j)
		    end do
		end do
!
        !-----------------------------------------------------------
        ! This algorithm assumes that the specified T and P correspond to
        ! vapor-liquid separation predicted by the provided model (0<beta<1)

        KFACT = (PCn/P) *EXP(5.373*(1+omgn)*(1-TCn/T))
        do while (g0<0.or.g1>0)
            g0 = sum(z*KFACT) - 1.D0
            g1 = 1.D0 - sum(z/KFACT)
            if (g0<0) then
                KFACT = 1.1*KFACT  ! increased volatiliy will bring the solution from subcooled liquid into VLE
            else if (g1>0) then
                KFACT = 0.9*KFACT  ! decreased volatiliy will bring the solution from superheated vapor into VLE
            end if
        end do
        LOG_K = LOG(KFACT)
        ! now we must have  g0>0 and g1<0 and therefore 0<beta<1 (M&M page 252)
        call betalimits (n,z,KFACT,bmin,bmax)
        beta = (bmin+bmax)/2  ! first guess for beta
        ! Succesive sustitution loop starts here
        var_K=1.0
        do while (maxval(abs(var_K)) > 1.d-4)
            ! Newton starts here
            g = 1.0
            do while (abs(g)>1.d-4)
                denom = 1+beta*(KFACT-1.D0)
                g = sum(z*(KFACT-1.D0) / denom)
                dg = -sum(z*(KFACT-1.D0)**2 / denom**2)
                beta = beta - g/dg
            end do
            denom = 1+beta*(KFACT-1.D0)
            y = z * KFACT / denom
            x = y / KFACT
            ! nc,MTYP,INDIC,T,P,rn,V,PHILOG,DLPHI
            call TERMO(n,-1,1,T,P,y,Vy,PHILOGy,DLPHIP,DLPHIT,FUGN)
            call TERMO(n, 1,1,T,P,x,Vx,PHILOGx,DLPHIP,DLPHIT,FUGN)
            var_K = PHILOGx - PHILOGy - LOG_K  ! variation in LOG_K = new - old
            LOG_K = PHILOGx - PHILOGy
            KFACT = exp(LOG_K)
        end do
        rho_x = 1/Vx
        rho_y = 1/Vy
        !-----------------------------------------------------------

        print *, x
 		print *, y
 		print *, rho_x
 		print *, rho_y
 		print *, beta
    end subroutine flash

	subroutine betalimits (n,z,KFACT,bmin,bmax)

      implicit none

        integer, intent(in) :: n  ! number of compounds in the system
        real*8, dimension(n), intent(in) :: z, KFACT  ! composition of the system and K factors
        real*8, intent(out) :: bmin, bmax
        real*8, dimension(n) :: vmin, vmax
        integer :: i, in=0, ix=0

        vmin=0.d0
        vmax=1.d0
        do i=1,n
            if (KFACT(i)*z(i)>1)then
                in = in+1
                vmin(in) = (KFACT(i)*z(i)-1.d0)/(KFACT(i)-1.d0)
            else if (KFACT(i)<z(i))then
                ix = ix+1
                vmax(ix) = (1.d0-z(i))/(1.d0-KFACT(i))
            end if
        end do
        bmin = maxval(vmin)
        bmax = minval(vmax)

    end subroutine betalimits

	subroutine readRKPRNC(nc,nin,nout)
      implicit DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (nco=10,RGAS=0.08314472d0)
!  Critical constants must be given in K and bar
!  b will be in L/mol and ac in bar*(L/mol)**2
	DOUBLE PRECISION Kij(nco,nco),lij(nco,nco),Kinf
!      PARAMETER (A0=0.0017,B0=1.9681,C0=-2.7238)
!      PARAMETER (A1=-2.4407,B1=7.4513,C1=12.504)
	dimension ac(nco),b(nco),del1(nco),rk(nco),diam(nc),vc(nc)
	dimension D(6),Vceos(nc)
      CHARACTER*10 fluid(nco)
!	COMMON/names/fluid
    COMMON/CRIT/TC(nco),PC(nco),DCeos(nco),OM(nco)
	COMMON /COMPONENTS/ ac,b,del1,rk,Kij,NTDEP
!	COMMON /K12/ K12
!	COMMON/COVOL/bb1(nco)
	COMMON /rule/ncomb
	COMMON /bcross/bij(nco,nco)
	COMMON /Tdep/ Kinf,Tstar
	COMMON /lforin/lij
!  D=[0.428363,18.496215,0.338426,0.660,789.723105,2.512392]
	read(NIN,*) ncomb,NTDEP
	do i=1,nc
	READ(NIN,'(A)')fluid(i)
	READ(NIN,*)Tc(i),Pc(i),OM(i),Vceos(i),Zrat
	RT=RGAS*Tc(i)
	Zc=Pc(i)*Vceos(i)/RT
	Zcin=Zc/Zrat
	Vc(i)=Vceos(i)/Zrat
	dceos(i)=1/Vceos(i)
	READ(NIN,*)ac(i),b(i),del1(i),rk(i)
! 4	bb1(i)=b(i)
	write(nout,'(A)')fluid(i)
	write(nout,1)Tc(i),Pc(i),Vc(i),OM(i)
	write(nout,3)Zcin,Zrat,Zc,Vceos(i)
	write(nout,2)ac(i),b(i),del1(i),rk(i)
	Kij(i,i)=0.0D0
	Lij(i,i)=0.0D0
	IF(i.gt.1)then
		if(ncomb.lt.2)then
			READ(NIN,*) (Kij(j,i),j=1,i-1)
			if(NTDEP.ge.1)READ(NIN,*) Kinf
			if(NTDEP.eq.1)READ(NIN,*)Tstar
			READ(NIN,*) (lij(j,i),j=1,i-1)
			lij(i,1:i-1)=lij(1:i-1,i)
		end if
	ENDIF
	end do
!	if(NTDEP.eq.1)Tstar=258.0
!  if(NTDEP.eq.2)READ(NIN,*)Tstar
	write(NOUT,*)
	write(nout,*)'Tc, Pc and Vc are given in K, bar and L/mol respectively'
!	K12=Kij(1,2)
 1	FORMAT('Tc=',F9.4,'   Pc =',F9.4,'   Vc =',F8.4,'   OM =',F7.4)
 3	FORMAT('Zc=',F9.4,' Zcrat=',F9.4,' Zceos=',F8.4,' Vceos=',F7.4)
 2	FORMAT('ac=',F9.4,'    b =',F9.4,'  del1=',F8.4,'    k =',F7.4)
	write(NOUT,*)
	if(ncomb.lt.2)then
		if(NTDEP.EQ.0)then
!  		write(NOUT,*)' K12 = ',Kij(1,2)
!  		write(NOUT,*)
		else
			write(NOUT,*)' K012 = ',Kij(1,2)
			write(NOUT,*)
			write(NOUT,*)' Kinf = ',Kinf
			write(NOUT,*)
			write(NOUT,*)'Tstar = ',Tstar
			write(NOUT,*)
		end if
		write(NOUT,*)'  KIJ MATRIX'
		DO I=1,NC
		write(NOUT,6)FLUID(I),(Kij(j,i),j=1,i-1)
		END DO
		write(NOUT,*)
		write(NOUT,*)'  LIJ MATRIX'
		DO I=1,NC
		write(NOUT,6)FLUID(I),(Lij(j,i),j=1,i-1)
		END DO
	else
		if(NTDEP.EQ.0)then
			write(NOUT,*)' Kijk:     112      122'
			write(NOUT,7)K01,K02
			write(NOUT,*)
		else
			write(NOUT,*)' K0ijk:    112      122'
			write(NOUT,7)K01,K02
			write(NOUT,*)
			write(NOUT,*)'Kinfijk:   112      122'
			write(NOUT,7)Kinf1,Kinf2
			write(NOUT,*)
			write(NOUT,*)'Tstar  :   112      122'
			write(NOUT,8)Tstar1,Tstar2
			write(NOUT,*)
		end if
		if(NTDEP.EQ.2)then
			write(NOUT,*)' Cijk:     112      122'
			write(NOUT,7)C1,C2
			write(NOUT,*)
		end if
		write(NOUT,*)' Lijk:     112      122'
!  	write(NOUT,7)Lijk(1,1,2),Lijk(1,2,2)
		write(NOUT,*)
	end if
	write(NOUT,*)
	write(NOUT,*)' Combining rules:'
	if(ncomb.eq.0)then
	write(NOUT,*)' 0: Classical or van der Waals '
		do i=1,nc
		do j=i,nc
		bij(i,j)=(1-lij(i,j))*(b(i)+b(j))/2
		bij(j,i)=bij(i,j)
		end do
		end do
	else if(ncomb.eq.3)then
	else
	write(NOUT,*)' 1: Lorentz-Berthelot'
		third=1.0d0/3
		do i=1,nc
		diam(i)=b(i)**third
		end do
		do i=1,nc
		do j=i,nc
		bij(i,j)=((1-lij(i,j))*(diam(i)+diam(j))/2)**3
		bij(j,i)=bij(i,j)
		end do
		end do
	end if
 6	FORMAT(A10,4F8.5)
 7	FORMAT(9x,F7.4,2x,F7.4)
 8	FORMAT(9x,F7.2,2x,F7.2)
	end

!  Then Kij values will be called indicating the lower index first, e.g. Kij(1,3)

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
      PARAMETER (nco=10)
	DOUBLE PRECISION Kinf,Kij0(nco,nco),Kij(nco,nco)
	dimension ai(nc),daidT(nc),daidT2(nc)
	dimension aij(nc,nc),daijdT(nc,nc),daijdT2(nc,nc)
    COMMON /CRIT/TC(nco),PC(nco),DCeos(nco),OM(nco)
	COMMON /COMPONENTS/ ac(nco),b(nco),d1(nco),rk(nco),Kij0,NTDEP
	COMMON /bcross/bij(nco,nco)
	COMMON /rule/ncomb
	COMMON /Tdep/ Kinf,Tstar
	IF(NTDEP.GE.1)THEN
		Kij=0.0D0
		Kij(1,2)=Kinf+Kij0(1,2)*exp(-T/Tstar)
		Kij(2,1)=Kij(1,2)
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
	IF(NTDEP.ge.1.and.NTD.EQ.1)THEN
		aux=daijdT(1,2)
		ratK=Kij(1,2)/(1-Kij(1,2))/Tstar
		daijdT(1,2)=aux+aij(1,2)*ratK
		daijdT(2,1)=daijdT(1,2)
		daijdT2(1,2)=daijdT2(1,2)+(2*aux-aij(1,2)/Tstar)*ratK
		daijdT2(2,1)=daijdT2(1,2)
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
    PARAMETER (nco=10)

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
      PARAMETER (nco=10)
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
	continue
end

SUBROUTINE TERMO(nc,MTYP,INDIC,T,P,rn,V,PHILOG,DLPHIP, DLPHIT,FUGN)
!       MTYP        TYPE OF ROOT DESIRED (-1 vapor, 1 liquid, 0 lower Gibbs energy phase)
!       rn		    mixture mole numbers                     (input)
!       t			temperature (k)                          (input)
!       p			pressure    (bar)                        (input)
!       v			volume	    (L)			                 (output)
!       PHILOG    vector of ln(phi(i)*P)                     (output)	INDIC < 5
!       DLPHIT    t-derivative of ln(phi(i)) (const P, n)    (output)	INDIC = 2 or 4
!       DLPHIP    P-derivative of ln(phi(i)) (const T, n)    (output)	INDIC < 5
!       FUGN      comp-derivative of ln(phi(i)) (const t & P)(output)	INDIC > 2
!  ---------------------------------------------------
    implicit none

    integer, intent(in) :: nc, indic
    integer ::  mtyp
    real*8, intent(in) :: t, p
    real*8, intent(out) :: v
    real*8, dimension(nc), intent(out) :: PHILOG, DLPHIT, DLPHIP
    real*8, dimension(nc,nc), intent(out) :: FUGN
	real*8, dimension(nc) :: rn, Arn, ArVn, ArTn,DPDN
	real*8, dimension(nc,nc) :: Arn2


	integer :: maxc, i, k, igz, nder, ntemp
	real*8 :: rgas, dpv, rt, totn, z, Ar, arv, arv2, artv, dpdt
	PARAMETER (MAXC=10,RGAS=0.08314472d0)

	! dimension rn(nc),Arn(nc),ArVn(nc),ArTn(nc),Arn2(nc,nc)

!  The output PHILOG is actually the vector ln(phi(i)*P)
	NTEMP=0
    IGZ=0
    NDER=1
    IF (INDIC.GT.2) NDER=2
    IF (INDIC.EQ.2 .OR. INDIC.EQ.4) NTEMP=1
	TOTN = sum(rn)
	if(P.le.0.0d0)MTYP=1
	CALL VCALC(MTYP,NC,rn,T,P,V)
      RT = RGAS*T
      Z = V/(TOTN*RT)	! this is Z/P

 	call ArVnder(nc,NDER,NTEMP,rn,V,T,Ar,ArV,ArTV,ArV2,Arn,ArVn,ArTn,Arn2)

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

   subroutine PUREFUG_CALC(nc,icomp,T,P,V,phi)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (RGAS=0.08314472d0) !nc=2,
	dimension rn(nc),Arn(nc),ArVn(nc),ArTn(nc),Arn2(nc,nc)
	rn=0.0
	rn(icomp)=1.0
	RT = RGAS*T
      Z = P*V/RT
	call ArVnder(nc,0,0,rn,V,T,Ar,ArV,ArTV,ArV2,Arn,ArVn,ArTn,Arn2)
      PHILOG=-LOG(Z)+Arn(icomp)/RT
	phi=exp(PHILOG)
	return
	END

      SUBROUTINE VCALC(ITYP,nc,rn,T,P,V)

!     ROUTINE FOR CALCULATION OF VOLUME, GIVEN PRESSURE

!     INPUT:

!     ITYP:        TYPE OF ROOT DESIRED (-1 vapor, 1 liquid, 0 lower Gibbs energy phase)
!     NC:          NO. OF COMPONENTS
!     rn:          FEED MOLES
!     T:           TEMPERATURE
!     P:           PRESSURE

!     OUTPUT:

!     V:           VOLUME

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (RGAS=0.08314472d0) ! nc=2,
	dimension rn(nc),dBi(nc),dBij(nc,nc)
	dimension Arn(nc),ArVn(nc),ArTn(nc),Arn2(nc,nc)
      LOGICAL FIRST_RUN
	NDER=0
      FIRST_RUN = .TRUE.
	TOTN = sum(rn)
	call Bcalc(nc,rn,T,B)
	CPV=B
      S3R = 1.D0/CPV
      ITER = 0

      ZETMIN = 0.D0
      ZETMAX = 1.D0-0.01*T/500	!.99D0  This is flexible for low T (V very close to B)
      IF (ITYP .GT. 0) THEN
         ZETA = .5D0
         ELSE
!..............IDEAL GAS ESTIMATE
         ZETA = MIN (.5D0,CPV*P/(TOTN*RGAS*T))
      ENDIF
  100 CONTINUE
      V = CPV/ZETA
      ITER = ITER + 1
	call ArVnder(nc,NDER,NTEMP,rn,V,T,Ar,ArV,ArTV,ArV2,Arn,ArVn,ArTn,Arn2)
      PCALC = TOTN*RGAS*T/V - ArV
      IF (PCALC .GT. P) THEN
         ZETMAX = ZETA
         ELSE
         ZETMIN = ZETA
      ENDIF
      AT  = (Ar + V*P) /(T*RGAS) - TOTN*LOG(V)
!  AT is something close to Gr(P,T)
      DER = (ArV2*V**2+TOTN*RGAS*T)*S3R  ! this is dPdrho/B
      DEL = -(PCALC-P)/DER
      ZETA = ZETA + MAX (MIN(DEL,0.1D0),-.1D0)
      IF (ZETA .GT. ZETMAX .OR. ZETA .LT. ZETMIN)   &
           ZETA = .5D0*(ZETMAX+ZETMIN)
      IF (ABS(PCALC-P) .LT. 1D-12) GOTO 101
      IF (ABS(DEL) .GT. 1D-10) GOTO 100
 101	IF (ITYP .EQ. 0 ) THEN

! FIRST RUN WAS VAPOUR; RERUN FOR LIQUID

         IF (FIRST_RUN) THEN
            VVAP = V
            AVAP = AT
            FIRST_RUN = .FALSE.
            ZETA = 0.5D0
	      ZETMAX = 1.D0-0.01*T/500
            GOTO 100
            ELSE
            IF (AT .GT. AVAP) V = VVAP
          ENDIF
      ENDIF
      END

	SUBROUTINE ArVnder(nc,NDER,NTD,rn,V,T,Ar,ArV,ArTV,ArV2,Arn,ArVn,ArTn,Arn2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	dimension rn(nc),Arn(nc),ArVn(nc),ArTn(nc),Arn2(nc,nc)
	COMMON /MODEL/ NMODEL
	IF(NMODEL.LE.2)THEN
!  SRK or PR
	  CALL HelmSRKPR(nc,NDER,NTD,rn,V,T,Ar,ArV,ArTV,ArV2,Arn,ArVn,ArTn,Arn2)
	ELSE IF (NMODEL.EQ.3) THEN

		CALL HelmRKPR(nc,NDER,NTD,rn,V,T,Ar,ArV,ArTV,ArV2,Arn,ArVn,ArTn,Arn2)
 	ELSE IF (NMODEL.EQ.4) THEN
!  	CALL HelmPCSAFT(NDER,NTD,rn,V,T,Ar,ArV,ArTV,ArV2,Arn,ArVn,ArTn,Arn2)
 	ELSE IF (NMODEL.EQ.6) THEN
!  	CALL HelmSPHCT(NDER,NTD,rn,V,T,Ar,ArV,ArTV,ArV2,Arn,ArVn,ArTn,Arn2)
 	ELSE IF (NMODEL.EQ.8) THEN
! 	  CALL HelmESD  (NDER,NTD,rn,V,T,Ar,ArV,ArTV,ArV2,Arn,ArVn,ArTn,Arn2)
	ELSE	!	GC-EOS 5 (or GCA 7)
!  	CALL HelmGC(NDER,NTD,rn,V,T,Ar,ArV,ArTV,ArV2,Arn,ArVn,ArTn,Arn2)
	END IF
      END
!
	SUBROUTINE Bcalc(nc,x,T,BMIX)
!  This general subroutine provides the "co-volume" for specified composition,
!  that will be used by Evalsecond or Vcalc
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXC=10,RGAS=0.08314472d0)
	DIMENSION x(nc),dBi(nc),dBij(nc,nc)
!      DIMENSION DD(0:3,MAXC),DDT(0:3,MAXC),DTT(0:3,MAXC),DIA(MAXC),DDB(0:3,MAXC)
	COMMON /MODEL/ NMODEL
!  COMMON/covol/VX(nc)  ! ESD
!      COMMON/MOL/DC(2),D(2),DT(2),HA(2),HB(2)
!      COMMON/MIXT/VCPM,CMIX,CVYM,CVYMT,dCVYM(MAXC),dCMIX(MAXC),
!     *	dCVYMT(MAXC),d2CVYM(MAXC,MAXC),d2CMIX(MAXC,MAXC)
!      COMMON/MIXRULE/NSUB
!      COMMON/BMIX/B
!  COMMON/forB/DDB
!	COMMON/NG/NGR
	COMMON /rule/ncomb
	NG=NGR
	if(NMODEL.EQ.5.OR.NMODEL.EQ.7)then
!  	CALL PARAGC(T,nc,NG,1)
!  	PI=3.1415926536D0
!  	XLAM3=0.0d0
!  	DO 3 I=1,nc
!  	DGC=D(I)
! 3	    XLAM3=XLAM3+X(I)*DGC**3
!  	B=PI/6.D0*XLAM3/1.0D3
	else if(NMODEL.EQ.4)then
!  DD=DDB
!      CALL DIAMET(nc,T,DIA,DD,DDT,DTT,NSUB)
!          B=RGAS*(DD(3,1)*X(1)+DD(3,2)*X(2))	!S3
	else if(NMODEL.EQ.6)then
!  	CALL Mixture_Param(NSUB,NC,X,T)
!  	B=VCPM
	else if(NMODEL.EQ.8)then
!  	B=x(1)*VX(1)+x(2)*VX(2)
	else
	if(ncomb<=2)then
		call Bnder(nc,x,B,dBi,dBij)	! Bmix is used in EVALSECOND
	else
!  	call Bcubicnder(2,x,B,dBi,dBij)
	end if
	end if
	BMIX=B
      END
