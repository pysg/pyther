      implicit DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*30 INFILE,OUTFILE
      CHARACTER*12 NAME
	COMMON /ZCRIT/ ZC
	COMMON/CONSTANTS/c(12),n(12),t(12),p(12)
      WRITE (*,*) 'ENTER INFILE'
      READ (*,'(A)') INFILE
	NUNIT= 30
	NOUT = 31
      OPEN(NUNIT,FILE=INFILE)
C
C
      WRITE (*,*) 'ENTER A NAME FOR THE OUTFILE'
      READ (*,'(A)') OUTFILE
      OPEN(NOUT,FILE=OUTFILE)
      READ(NUNIT,*) NAME,Zc,(c(j),j=1,12)
	s1=c(1)+c(2)+c(3)
      READ(NUNIT,*) Tc,Pc,dc 
      WRITE (*,*) 'ENTER 0 FOR REDUCED TEMPERATURE Tr'
      WRITE (*,*) '   OR 1 FOR ABSOLUTE TEMPERATURE T'
      READ (*,*) NOPT
	IF (NOPT.EQ.0) THEN
	    WRITE (*,*) 'ENTER REDUCED TEMPERATURE Tr'
		READ (*,*) Tr
	ELSE
	    WRITE (*,*) 'ENTER ABSOLUTE TEMPERATURE T'
		READ (*,*) TEM
		Tr=TEM/Tc
	END IF
	WRITE (NOUT,*) ' Isotherm of',NAME,
	1	'according to Span & Wagner equation (2003)'
	WRITE (NOUT,*) ' for Tr = ',Tr
	if (nopt.eq.1)WRITE (NOUT,*) ' T = ',TEM
	WRITE (NOUT,*) 
	WRITE (nout,*) 'Rho(MOL/L)   P(bar)       phi          Z'
 16	FORMAT(' RHOr   ',A12)
 17	FORMAT(2x,F8.4,2x,F10.4,2x,F10.5,2x,F10.5)
      WRITE (NOUT,17) d,Pr,1.0d0,1.0d0
	t(1)=0.25
	t(2)=1.25
	t(3)=1.5
	t(4)=0.25
	t(5)=0.875
	t(6)=2.375
	t(7)=2
	t(8)=2.125
	t(9)=3.5
	t(10)=6.5
	t(11)=4.75
	t(12)=12.5
	n(6)=1
	n(7)=2
	n(8)=5
	n(9)=1
	n(10)=1
	n(11)=4
	n(12)=2
	p(6)=1
	p(7)=1
	p(8)=1
	p(9)=2
	p(10)=2
	p(11)=2
	p(12)=3
c	if (Tr.lt.1.0) go to 10  ! 6/7/2015 cancelling (primitive) Pv part in order to have the full isotherm also for subcritical
c	Critical or Supercritical isotherms
	step=0.05
	do d=step,3.01D0,step
	call ZetaCalc(Tr,d,zeta)
	Pr=d*Tr*zeta/Zc
	call FUG_CALC(Tr,Pr,d,phi)
      WRITE (NOUT,17) d*dc,Pr*Pc,phi,zeta
	end do
	go to 11
c	Subcritical isotherms
 10	PVini=0.20 !0.20 this INITIAL VALUES for reduced variables are OK for Tr=0.8
	RHOL=2.30
	RHOV=0.08
C	RHOL=1.50
C	RHOV=0.50
	call VaporPressure(Tr,PVini,Pv,RHOL,RHOV)
c	vapour branch
	step=RHOV/10
	do d=step,RHOV,step
	call ZetaCalc(Tr,d,zeta)
	Pr=d*Tr*zeta/Zc
	call FUG_CALC(Tr,Pr,d,phi)
      WRITE (NOUT,17) d*dc,Pr*Pc,phi,zeta
	end do
      WRITE (NOUT,*)
c	liquid branch
	do d=RHOL,3.01D0,step
	call ZetaCalc(Tr,d,zeta)
	Pr=d*Tr*zeta/Zc
	call FUG_CALC(Tr,Pr,d,phi)
      WRITE (NOUT,17) d*dc,Pr*Pc,phi,zeta
	end do
 11	end
c
	SUBROUTINE FindDensr(iph,Tr,Pr,dini,d)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	COMMON /ZCRIT/ ZC
	tol=1.0d-5
	d0=dini
	call ZetaCalc(Tr,d0,zeta)
	Pr0=d0*Tr*zeta/Zc
	if(iph.gt.0)d=dini*(1-0.02*Tr)
	if(iph.lt.0)d=dini*(1+0.02*Tr)
	call ZetaCalc(Tr,d,zeta)
	Pr1=d*Tr*zeta/Zc
	d1=d
 2	RK=(d1-d0)/(Pr1-Pr0)
	d=d0+RK*(Pr-Pr0)
	Pr0=Pr1
	call ZetaCalc(Tr,d,zeta)
	Pr1=d*Tr*zeta/Zc
	if(abs(Pr1-Pr).lt.tol)goto 3
	d0=d1
	d1=d
	goto 2
 3	end
c
	SUBROUTINE ZetaCalc(Tr,d,zeta)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION frho(12)
	COMMON/CONSTANTS/c(12),n(12),t(12),p(12)
	s1=c(1)/Tr**t(1)+c(2)/Tr**t(2)+c(3)/Tr**t(3)
		do k=6,12
		frho(k)=d**n(k)*(n(k)-p(k)*d**p(k))*exp(-d**p(k))/Tr**t(k)
		end do
	zeta=1.0D0+s1*d+3*c(4)*d**3/Tr**t(4)+7*c(5)*d**7/Tr**t(5)
		do k=6,12
		zeta=zeta+c(k)*frho(k)
		end do
	end
c
	SUBROUTINE VaporPressure(Tr,PVini,Pv,RHOL,RHOV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ERRMAX=1.D-8)
	dphiold = 0.0D0
	P = PVini
	n=1
 30	dini=RHOL
	call FindDensr(1,Tr,P,dini,d)
	RHOL = d
	dini=RHOV
	call FindDensr(-1,Tr,P,dini,d) ! SOLVE for vapor density
	RHOV = d
	call FUG_CALC(Tr,P,RHOL,phi)
	phiL = phi
	call FUG_CALC(Tr,P,d,phi)
	phiV = phi
	dphiold = dphi
	dphi = phiV - phiL
	IF (ABS(dphi).gt.ERRMAX) THEN
		Pold = Plast
		Plast = P
	if(dphiold.eq.0.0D0.or.Tr.gt.0.975) then
		P = P * (phiL/phiV)
	else
		P = Plast - dphi*(Plast-Pold)/(dphi-dphiold)
	end if
c		n=n+1
		GO TO 30
	END IF
	PV = P
c      WRITE (31,*) ' n=',n
	return
	END
C
	SUBROUTINE FUG_CALC(Tr,Pr,d,phi)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	COMMON /ZCRIT/ ZC
c	P and d are reduced variables
	Z = Zc*Pr/d/Tr
	CALL redAres(Tr,d,ra)
      PHI=EXP(ra+Z-1)/Z
	return
	END
c
	SUBROUTINE redAres(Tr,d,ra)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION frho(12)
	COMMON/CONSTANTS/c(12),n(12),t(12),p(12)
	s1=c(1)/Tr**t(1)+c(2)/Tr**t(2)+c(3)/Tr**t(3)
		do k=6,12
		frho(k)=d**n(k)*exp(-d**p(k))/Tr**t(k)
		end do
	ra=s1*d+c(4)*d**3/Tr**t(4)+c(5)*d**7/Tr**t(5)
		do k=6,12
		ra=ra+c(k)*frho(k)
		end do
	end
c
