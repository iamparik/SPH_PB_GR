!****************************************************************************
!
!  Subroutine: gamma_analytical_Ferrand_2D
!
!  PURPOSE: Calculates the gamma_as, the derivative of renormalization factor 
!           split into contirbutions from multiple edges, as 
!           defined by Leroy et.al USAW 2014
!
!   CREATED:        11/13/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  11/13/2022       by  PARIKSHIT BOREGOWDA 
!****************************************************************************
subroutine gamma_analytical_leroy_2D(xa, xv, s_norm, gamma_as,h, pi)

implicit none

real(8), intent(in), dimension(2) :: xa
real(8), intent(in), dimension(2,2) :: xv 
real(8), intent(in), dimension(2) :: s_norm
real(8), intent(in) :: h, pi
real(8), intent(out) :: gamma_as
!integer(4) :: d, i
real(8) :: zeta(2), ra0(2), rav1(2), rav2(2), qa0, qa0_tilda, y1, y2, s_tangent(2), &
        & phi1, phi2

rav1(:)= xv(:,1) - xa(:) !xa(:)-xv(:,1)
rav2(:)= xv(:,2) - xa(:) ! xa(:)-xv(:,2)

if(norm2(xv(:,2)- xv(:,1)) .ne. 0.D0) then 
    s_tangent= (xv(:,2)- xv(:,1))/norm2(xv(:,2)- xv(:,1))
else
    s_tangent=0.D0
    write(*,*) 'wrong vertices',xv(:,1) ,' and', xv(:,2)
endif

! find r0
ra0=xa(:) - xv(:,1)
ra0= xv(:,1) + s_tangent(:)*dot_product(ra0(:),s_tangent(:)) 
ra0(:)= xa(:)-ra0(:) 
qa0= norm2(ra0)/h
qa0_tilda= min( qa0, 2.D0)

!The below is added to avoid NAN 
qa0= qa0_tilda

y1= dot_product(rav1,s_tangent)
y2= dot_product(rav2,s_tangent)

zeta(1)= min(0.5D0*sqrt(qa0_tilda*qa0_tilda+(y1/h)**2),1.D0)
zeta(2)= min(0.5D0*sqrt(qa0_tilda*qa0_tilda+(y2/h)**2),1.D0)

if(qa0 .eq. 0.D0) then
   phi1=pi       
   phi2=pi
    
else           

    phi1=qa0*sqrt(zeta(1)**2 -qa0*qa0/4.D0)*(-(zeta(1)**5)*4.D0/3.D0 + 7.D0*(zeta(1)**4) - (qa0*qa0*5.D0/12.D0+14.D0)*(zeta(1)**3) &
                & + (qa0*qa0+5.D0)*(zeta(1)**2)*7.D0/3.D0 - (qa0*qa0*5.D0/8.D0+21.D0)*zeta(1)*qa0*qa0/4.D0 + 7.D0*(qa0**4)/6.D0 + 35.D0*(qa0**2)/6.D0 - 7.D0) &
                & - (5.D0*(qa0**2)/8.D0 + 21.D0)*(qa0**5)*dacosh(2*zeta(1)/qa0)/16.D0 + 2.D0*datan(sqrt(4.D0*zeta(1)*zeta(1)/(qa0*qa0) - 1.D0))


    phi2=qa0*sqrt(zeta(2)**2 -qa0*qa0/4.D0)*(-(zeta(2)**5)*4.D0/3.D0 + 7.D0*(zeta(2)**4) - (qa0*qa0*5.D0/12.D0+14.D0)*(zeta(2)**3) + &
                & (qa0*qa0+5.D0)*(zeta(2)**2)*7.D0/3.D0 - (qa0*qa0*5.D0/8.D0+21.D0)*zeta(2)*qa0*qa0/4.D0 + 7.D0*(qa0**4)/6.D0 + 35.D0*(qa0**2)/6.D0 - 7.D0) &
                & - (5.D0*(qa0**2)/8.D0 + 21.D0)*(qa0**5)*dacosh(2*zeta(2)/qa0)/16.D0 + 2.D0*datan(sqrt(4.D0*zeta(2)*zeta(2)/(qa0*qa0) - 1.D0))
    

endif


gamma_as= sign(1.D0, dot_product(s_norm,ra0))*(dsign(1.D0,y2)*phi2 - dsign(1.D0,y1)*phi1)/(4.D0*pi)


    if(isNAN(gamma_as)) then
        write(*,*) "gamma_as is NAN for xa =", xa , " and xv =", xv
        pause
    endif


end