2
 program calc_flash

    implicit DOUBLE PRECISION (A-H,O-Z)
    OPEN(1,FILE='flashIN.txt')
    OPEN(2,FILE='flashOUT.txt')
    read(1,*) N
    call readcase(n)
END program

subroutine readcase(n)
    implicit DOUBLE PRECISION (A-H,O-Z)
    PARAMETER (nco=64)
	dimension z(n)
    CHARACTER*4 spec
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
        real*8, dimension(n,n) :: Tstarn   

    ! interaction parameters matrices
    real*8, dimension(nco,nco) :: Kij_or_K0, Lij, Tstar

    real*8, dimension(nco) :: x , y  ! composition of liquid and vapour (molar fractions)

        ! Auxiliary variables not used here, but required to call zTVTERMO
        real*8, dimension(n) :: PHILOG
        real*8, dimension(n, n) :: FUGN

        COMMON /CRIT/TC(nco),PC(nco),DCeos(nco),omg(nco)
	    COMMON /MODEL/ NMODEL
	    COMMON /Vshift/ iVshift, Vs(nco)    ! added June 2016
	    COMMON /COMPONENTS/ ac(nco),b(nco),delta1(nco),rk_or_m(nco),Kij_or_K0,NTDEP
	    COMMON /rule/ncomb
	    COMMON /bcross/bij(nco,nco)
	    COMMON /Tdep/ Kinf,Tstar    
	    COMMON /lforin/lij

        READ(1,*) (z(j),j=1,N)
        read(1,*) T
        read(1,*) P  ! bar
        v=1.0
        read(1,*) nmodel, iVshift  
      if(nmodel<3)then
        call read2PcubicNC(N,1,2)
      else if(nmodel==3)then
        call readRKPRNC(N,1,2)
      end if
      WRITE (2,*)
      write (2,4) (z(i),i=1,n)
 4	FORMAT('Molar fractions: ',20F7.4)
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
        Tstarn = Tstar(:n, :n)
        lijn = lij(:n, :n)

!   SE SUGIERE CREAR (Y LLAMAR AQUÍ) UNA RUTINA COMO LA SIGUIENTE. 
!   (AL FINAL DE ESTE ARCHIVO SE ENCUENTRA UN TEMPLATE VACÍO PARA LA MISMA)

!      call flash(nmodel, n, z, tcn, pcn, omgn, acn, bn, k_or_mn, delta1n, &
!                     Kij_or_K0n, Tstarn, Lijn, t, p, v, x, y, rho_x, rho_y, beta, iter)

        if(iVshift==1)then
            VsV = sum(y*Vs)
            VsL = sum(x*Vs)
            rho_y = 1 / (1/rho_y-VsV)
            rho_x = 1 / (1/rho_x-VsL)
        end if
        V = beta/rho_y + (1-beta)/rho_x
        WRITE (2,*) ' v(L/mol)=',v
          WRITE (2,*) 'Beta (vapor phase mol fraction)= ',beta
          WRITE (2,*) 'Comp',(i,i=1,N)
          WRITE (2,1) (x(j),j=1,N)
          WRITE (2,2) (y(j),j=1,N)
          WRITE (2,*) 'liquid density (moles/L)= ',rho_x
          WRITE (2,*) 'vapour density (moles/L)= ',rho_y
          betav = beta / (V * rho_y)
          WRITE (2,*) 'BetaVol (vapor phase volume fraction)= ',betav 
          WRITE (2,*)
          print *, betav 
          print *, P 
          print *, v 
 	      print *, xplus
     	  print *, yplus 
          read(1,*) P  ! bar
          WRITE (2,*) ' P(bar)=',P
    
 1	FORMAT('  x  ', 30E12.4)
 2	FORMAT('  y  ', 30E12.4)
 end subroutine readcase

subroutine flash(model, n, z, tcn, pcn, omgn, acn, bn, k_or_mn, delta1n, &
                     Kij_or_K0n, Tstarn, Lijn, t, p, v, x, y, rho_x, rho_y, beta, iter)



subroutine flash(model, n, z, tcn, pcn, omgn, acn, bn, k_or_mn, delta1n, &
                     Kij_or_K0n, Tstarn, Lijn, t, p, v, x, y, rho_x, rho_y, beta, iter)

    implicit DOUBLE PRECISION (A-H,O-Z)
    PARAMETER (nco=64)

        ! M&M means the book by Michelsen and Mollerup, 2nd Edition (2007)

        ! Flash specification, eos id and  number of compounds in the system
        LOGICAL stopflash
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
        real*8, dimension(n,n), intent(in) :: Tstarn 
        real*8, dimension(n,n), intent(in) :: Lijn

        ! Temperature and Pressure for the flash
        real*8, intent(in) :: t            ! Temperature for the flash (K)
        real*8 :: p            ! Pressure for the flash (bar)
        real*8 :: v          ! (L/mol) Molar volume resulting from (TP)

        ! Results from flash calculation
        real*8, dimension(n), intent(out) :: x  ! composition of liquid (molar fractions)
        real*8, dimension(n), intent(out) :: y  ! composition of vapour (molar fractions)
        real*8, intent(out) :: rho_x            ! density of liquid (moles/L)
        real*8, intent(out) :: rho_y            ! density of vapour (moles/L)
        real*8, intent(out) :: beta             ! total fraction of vapour (molar base)
        integer, intent(out) :: iter            ! number of iterations required to converge

        ! Intermediate variables during calculation process
        real*8, dimension(n) :: PHILOGy, PHILOGx, DLPHIT, DLPHIP
        real*8, dimension(n) :: KFACT, LOG_K, AUXK, var_K, denom, varKold, logKold
        real*8, dimension(n, n) :: FUGN
        real*8 :: g0, g1  ! function g valuated at beta=0 and 1, based on Wilson K factors
        real*8 :: g, dg, bmin, bmax, Vy, Vx

        real*8, dimension(nco,nco) :: Kij_or_K0, Tstar
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
        Tstar(:n, :n) = Tstarn
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


!   Cuando sea necesario calcular los (logaritmos de) coeficientes de fugacidad, pueden obtenerse mediante los siguientes
!   llamados a la rutina "TERMO", que además devolverá el volumen molar de la fase.
!                call TERMO(n,-1,1,T,P,y,Vy,PHILOGy,DLPHIP,DLPHIT,FUGN)
!                call TERMO(n, 1,1,T,P,x,Vx,PHILOGx,DLPHIP,DLPHIT,FUGN)


  end subroutine flash
