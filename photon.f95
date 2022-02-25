 	program MB
!    -------------------------------------------------------------------------
!	This program is written by Prashant Kumar Chauhan
!	It calculates complex optical conductivity 'sigma(q=0, omega)' of a 
!	superconductor with arbitrary electron collision time 'tau' from 
!	from BCS-theory; q=0 means local electrodynamics, london limit; omega=
!	circular frequency; normalization: sigma(omega=0)=1, Input: x=omega/2Delta
!	y=1/(2*Delta*tau), tt=temperature/Tc. (Delta=gap, Tc=critical temp, one 
!	needs t=T/2Delta) Optput: s=sigma = sigma1 +i*sigma2[E.H. Brandt, 18.12.1989]
!	Remark: Since the integrande (varaible e=E/2Delta) diverge at the boundaries, 
!	one substitutes e=e(u) with de/du=O at the integration boundaries. The new 
!	integration variable u goes from 0 to I with step width dx=I/integer
!    -------------------------------------------------------------------------
	implicit none

	integer             :: i,j,k
	real*8    :: x,y,tt,te,v1,v2,diff
    complex(8) 			:: s

!	open(unit = 100, file = 'freq.dat', status = 'old', action = 'read')
    open(11,file='sigma.dat',status='unknown')
    
        
	write(*,*) "whats t/tc:"
	read(*,*) tt
    write(*,*) "whats y:"
    read(*,*) y
! v1 and v2 are starting and ending value of hv/2Delta and diff is their difference
    write(*,*) "specify starting and ending value of hv/2Delta:"
     read(*,*) v1, v2
!    te = tt/(3.528*sqrt(1-tt)*(0.9963+ 0.7733*tt))
!    print*, "value of t",te
! j is the number of frequency points (for now let j=60)
	j=100
    diff = (v2 -v1)/DBLE(j)
! we start the first value of frequency
	x=0
	do i = 1,j
!		read(100,*) x
		x=x+diff
	    CALL BERS(x,y,tt,s)
!	    write(*,*) s
	    write(11,*)x,DBLE(s),DIMAG(s)
	    s=0
    end do
 !   print *,"the value of t/tc is", tt
    End program
!-------------------------------------------------------------------------
! Subroutine to calculate the integrals (normalised coductivity) for each 
! frequency
!-------------------------------------------------------------------------
SUBROUTINE BERS(x,y,tt,s)
	real*8, INTENT(IN):: x,y,tt
    complex(8), INTENT(out):: s
	complex s1, s2, s3, GK
    parameter(M=40 , d1=1./M)
    double precision u
    real*8 t
    dx = 1./int(M*max(1.,sqrt(x)))
    t = tt/(3.528*sqrt(1-tt)*(0.9963+ 0.7733*tt))

    	s1=(0.,0.)
        s2=(0.,0.)
        s3=(0.,0.)
	u=dx*.5
	do 
    If (u>1) Exit
		s2 = s2 + GK(.5 + (u/(1. -u))**2, x, y, t, 2)*u/(1. -u)**3
		u=u+dx
	end do
		s = s2*dx*2
!        write(*,*) s, dx
        		if(x.lt.1) then
         u=dx*.5
    do 
      If (u>1) Exit
 		s1 = s1 +GK(.5 +x*u*u*(3. -u-u), x, y, t, 1)*u*(1. -u)
		u=u+dx
    end do
		s  = s + s1*dx*6.0*x
!        write(*,*) s, dx
        		else
	u=dx*.5
    do 
      If (u>1) Exit      
 		s3 = s3 + GK(.5 +(x-1.)*u*u*(3. -u-u), x, y, t, 3)*u*(1. -u)
        u=u+dx
    end do
  	 	u=d1*.5 	
	    do
        If (u>1) Exit  
 			s1 = s1 +GK(x-.5 +u*u*(3. -u-u), x, y, t, 1)*u*(1. -u)
            u=u+d1
        end do  
		s  = s + (s3*dx*(x-1.) +s1*d1)*6.
        		end if
        s  = s*cmplx(0.,y)*.5/x
        		end
!------------------------------------------------------------------------------
! Complex function to calculate the three different integrands. 
! gl, g2, g3 (=gk, k=1,2,3) 
!------------------------------------------------------------------------------
	COMPLEX FUNCTION GK(e, x, y, t, k)
			real*8,INTENT(IN):: e,x,y,t
    		complex cy, p4, c42
!	if(k.eq.2) p1=sqrt((e+x)**2 -.25)
   	if(k.eq.2) p1=0
    		   p2=sqrt(e*e      -.25)
	if(k.eq.3) p3=sqrt((e-x)**2 -.25)
	if(k.eq.1) p4=cmplx(0., sqrt(.25 -(e-x)**2))
               cy=cmplx(0., y)
    if(k.eq.1) c42=(.25 +e*(e-x))/(p4*p2 +1E-20)
!    if(k.eq.2) c12=(.25 +e*(e+x))/(p1*p2 +1E-20)
    if(k.eq.2) c12=0       	
    if(k.eq.3) c32=(.25 +e*(e-x))/(p3*p2 +1E-20)
      			th=tanh(e/(t+t+0.001))
    if(k.eq.1) GK= th* ((1-c42)/(p4+p2+cy) -(1+c42)/( p4-p2+cy))
!    if(k.eq.2) GK= tanh((e+x)/(t+t+0.001))*((1+c12)/( p1-p2+cy) - (1-c12)/(-p1-p2+cy)) +th*((1-c12)/(p1+p2+cy) -(1+c12)/(p1-p2+cy))
    if(k.eq.2) GK= 0
    if(k.eq.3) GK= th* ((1-c32)/(p3+p2+cy) -(1+c32)/( p3-p2+cy))
      			end
!-------------------------------------------------------------------------------