!****************************************************************************
!
!  Subroutine: dgamma_analytical_Ferrand_2D
!
!  PURPOSE: Calculates the del_gamma_as, the derivative of renormalization factor 
!           split into contirbutions from multipel edges, as 
!           defined by Ferand et.al USAW 2013
!
!   CREATED:        05/24/2021       by  PARIKSHIT BOREGOWDA
!   Last Modified:  04/25/2022       by  PARIKSHIT BOREGOWDA 
!****************************************************************************
subroutine dgamma_analytical_Ferrand_2D(xa,xei, xe_norm, dgamma_edge,h, pi, scale_k)

implicit none

real(8),  dimension(2) :: xa
real(8), dimension(2,2) :: xei
real(8), intent(in), dimension(2) :: xe_norm
real(8), intent(in) :: h,pi,scale_k
integer(4) :: d,i       
real(8) :: Q0, Q1, Q2, cosa1, cosa2, PQ1, PQ2, ip_1, ip_2, quintic_integ 
real(8), dimension(2), intent(out) :: dgamma_edge
real(8), dimension(2,2) :: raei !
real(8), dimension(2) :: s21   


! Find the vector connecting real particle a to edge particles
! e1 and e2, and the vector from e1 to e2 
1    do d=1, 2
        do i=1,2    ! The assumption is that for 1D there is 1 vertex, 2D there are two, and 3D three (However, 3D can have four vertex)
            !         ^
            !        /
            !       /
            !      ^
            raei(d,i)= xei(d,i) - xa(d)
        enddo
    !below needs to be changed for 1D and 3D
    ! o ----> o
    s21(d)= xei(d,2) - xei(d,1)
    enddo



! Q0 is the (smoothening length normalized) normal distance between the edge and particle a
Q0=abs(dot_product(raei(:,1),xe_norm(:)))

! Q1, Q2 is the normalized norm2 of rae1 and rae2 
Q1= norm2(raei(:,1))
Q2= norm2(raei(:,2))

Q0=Q0/h
Q1=Q1/h
Q2=Q2/h

! Q1cosa and Q2cosa are projection of rae1 and rae2 about s21 
cosa1= dot_product(raei(:,1),s21(:))/(norm2(raei(:,1))*norm2(s21(:))) 
cosa2= dot_product(raei(:,2),s21(:))/(norm2(raei(:,2))*norm2(s21(:))) 



! PQ1 and PQ2 are the polynomial function
PQ1= (Q1)**5.D0 * 7.D0/ 192.D0 - (Q1)**4.D0 * 21.D0/64.D0 + (Q1)**3.D0 * 35.D0/32.D0 - (Q1)**2.D0 *35.D0/24.D0 + 7.D0/4.D0 &
    & + (Q0)**2.D0 * ( (Q1)**3.D0 * 35.D0/768.D0 - (Q1)**2.D0* 7.D0/16.D0+ Q1*105.D0/64.D0 - 35.D0/12.D0 ) &
    & + (Q0)**4.D0 * ( Q1* 35.D0/512.D0 - 7.D0/8.D0 )

PQ2= (Q2)**5.D0 * 7.D0/ 192.D0 - (Q2)**4.D0 * 21.D0/64.D0 + (Q2)**3.D0 * 35.D0/32.D0 - (Q2)**2.D0 *35.D0/24.D0 + 7.D0/4.D0 &
    & + (Q0)**2.D0 * ( (Q2)**3.D0 * 35.D0/768.D0 - (Q2)**2.D0* 7.D0/16.D0+ Q2*105.D0/64.D0 - 35.D0/12.D0 ) &
    & + (Q0)**4.D0 * ( Q2* 35.D0/512.D0 - 7.D0/8.D0 )

ip_1= log((Q1 + abs(Q1*cosa1))/Q0)
ip_2= log((Q2 + abs(Q2*cosa2))/Q0)


! Part of the discretised Gradient of renormalisation factor for a particle about an edge 
quintic_integ= Q2*cosa2*PQ2/pi - Q1*cosa1 * PQ1 /pi + (Q0**4. /pi)*(105.D0/64.D0 + Q0**2. * 35.D0/512.D0) &
            & * (sign(ip_2, Q2*cosa2) - sign(ip_1, Q1*cosa1))

if (Q0.eq.0) quintic_integ= Q2*cosa2*PQ2/pi - Q1*cosa1 * PQ1 /pi 

if (Q1 .eq. 0) quintic_integ= Q2*cosa2*PQ2/pi   ! also, implies Q0=0

if (Q2 .eq. 0)  quintic_integ= - Q1*cosa1 * PQ1 /pi ! also, implies Q0=0

!if((Q1 .ge. scale_k) .and. (Q2 .ge. scale_k))   quintic_integ=0.D0   


if(isnan(quintic_integ)) pause "dgamma is NAN, check values of quintic_integ "

 dgamma_edge(:)=(quintic_integ/h)*xe_norm(:)
 
 
end