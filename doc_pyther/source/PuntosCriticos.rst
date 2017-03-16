12. Puntos Criticos
*******************
*******************



.. math:: n_1 = z_1 + s \sqrt(z_1) u_1

.. math:: n_2 = z_2 + s \sqrt(z_2) u_2

.. math:: u_1^2 + u_2^2 = 1 

.. math:: B_{ij} = \sqrt( z_iz_j) \left(\frac {\partial lnf_i}{\partial n_j} \right)_{T,V}

.. math:: b = \lambda_1 = 0

.. math:: c = \left( \frac {\partial \lambda_1} {\partial s} \right)_{s=0} = 0

.. math:: c \simeq \frac { \lambda_1 (s = \eta) - \lambda_1(s=-\eta)} {2\eta}

.. math:: \eta = 0.0001

.. math:: \lambda = \frac {-\beta \pm \sqrt{\beta^2 - 4} }{2}






subroutine bceval(z,T,V,P,b,c)
c	INPUT:  z,T,V
c	OUTPUT: P,b,c
c
C	The following subroutine must be included in a separate .for file:
C   XTVTERMO(NTYP,T,V,P,z,FUG,FUGT,FUGV,FUGN)
C   INPUT:
C     NTYP:   LEVEL OF DERIVATIVES, SEE BELOW
C     T:      TEMPERATURE (K)
C     V:      MOLAR VOLUME (L/MOL)
C     z:      COMPOSITION (MOLES, NEED NOT BE NORMALIZED)
C   OUTPUT:
C     P:      PRESSURE (bar)
C     FUG:    VECTOR OF LOG FUGACITIES (ALL NTYP)
C     FUGT:   T-DERIVATIVE OF FUG (NTYP = 2, 4 OR 5)
C     FUGV:   V-DERIVATIVE OF FUG (ALL NTYP)
C     FUGN:   MATRIX OF COMPOSITION DERIVATIVES OF FUG (NTYP >=3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (nco=2)
      DIMENSION z(nco)
      DIMENSION sqz(nco),ym(nco),u(nco),up(nco),y(nco)
	COMMON /Pder/ DPDN(nco),DPDT,DPDV
	COMMON /CAEPcond/ DPDVcri
	eps=1.0D-4
 1	sqz(1)=sqrt(z(1))
	sqz(2)=sqrt(z(2))
	call eigcalc(z,T,V,P,b,u)
	dpdvcri=dpdv
c	calculation of b at s=eps  (e)
	y(1)=z(1)+eps*u(1)*sqz(1)
	y(2)=z(2)+eps*u(2)*sqz(2)
	if(minval(y).lt.0)then
		call modifyz(z)
		go to 1
	end if
	call eigcalc(y,T,V,Pp,bpos,up)
c	calculation of b at s=-eps  (m)
	ym(1)=z(1)-eps*u(1)*sqz(1)
	ym(2)=z(2)-eps*u(2)*sqz(2)
	if(minval(ym).lt.0)then
		call modifyz(z)
		goto 1
	end if
	call eigcalc(ym,T,V,Pn,bneg,up)
c	calculation of c
	c=(bpos-bneg)/2.0/eps
	end
c
	subroutine modifyz(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION z(2)
	if(z(1).lt.z(2))then
		z(1)=2*z(1)
		z(2)=1.0d0-z(1)
	else
		z(2)=2*z(2)
		z(1)=1.0d0-z(2)
	end if
	end
c
	subroutine eigcalc(z,T,V,P,b,u)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (nco=2)
      DIMENSION z(nco),FUG(nco),FUGT(nco),FUGV(nco),FUGN(nco,nco)
      DIMENSION u(nco)
      jac=5 ! FUGN is required, but not FLT
      call XTVTERMO(jac,T,V,P,z,FUG,FUGT,FUGV,FUGN)
	bet=-z(1)*FUGN(1,1)-z(2)*FUGN(2,2)
	gam=z(1)*z(2)*(FUGN(1,1)*FUGN(2,2)-FUGN(1,2)**2)
	sq=sqrt(bet**2-4*gam)
	rlam1=(-bet+sq)/2
	rlam2=(-bet-sq)/2
	if(abs(rlam1).lt.abs(rlam2))then
		b=rlam1
	else
		b=rlam2
	end if
	u2=(b-z(1)*FUGN(1,1))/(sqrt(z(1)*z(2))*FUGN(1,2)) ! k=u2/u1=u2
	u(1)=sqrt(1/(1+u2*u2))  !normalization
	u(2)=sqrt(1-u(1)**2)
	if(u2.lt.0) u(2)=-u(2)
	end
C
C     purpose of routine CRITSTABCHECK:
C
C     To find the composition where the tangent plane distance respect to the 
C     critical composition takes on its minimum value at given T and P
C
C     Parameters:
C
C     T       (I)       Temperature
C     P       (I)       Pressure
C     Xc	    (I)       Composition of the critical point
C     W       (O)       Composition of the minimum tpd
C     tpdm    (O)       Value of the minimum tpd