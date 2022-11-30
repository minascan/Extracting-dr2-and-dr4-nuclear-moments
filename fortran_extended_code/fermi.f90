module fermi

  use defs
  use consts
  
  implicit none

  private :: radmoments

contains

  subroutine radmoments_diffs(Z,rch, dr_2,dr_4,dr_6,dr_8,QG_dr_2,QG_dr_4,QG_dr_6,QG_dr_8)
    
    integer(ik), intent(in) :: Z
    real(rk), intent(in) :: rch(3)
    real(rk), intent(out) :: dr_2, dr_4, dr_6, dr_8
    real(rk), intent(out) :: QG_dr_2, QG_dr_4, QG_dr_6, QG_dr_8
    real(rk) :: r_2(3), r_4(3), r_6(3), r_8(3)
    integer :: i
    
    do i=1,3
       call radmoments(Z,rch(i), r_2(i),r_4(i),r_6(i),r_8(i))
    end do
    
    dr_2 = r_2(1) - r_2(2)
    dr_4 = r_4(1) - r_4(2)
    dr_6 = r_6(1) - r_6(2)
    dr_8 = r_8(1) - r_8(2)
    
    QG_dr_2 = r_2(1) - r_2(3)
    QG_dr_4 = r_4(1) - r_4(3)
    QG_dr_6 = r_6(1) - r_6(3)
    QG_dr_8 = r_8(1) - r_8(3)
    
  end subroutine radmoments_diffs

  
  subroutine radmoments(Z,trms, r2,r4,r6,r8)

    integer(ik), intent(in) :: Z
    real(rk), intent(in) :: trms
    real(rk) :: c, a, rmax, h, h2, rms, drms
    real(rk) :: dnorm, norm, dinte, inte
    real(rk) :: dr2, dr4, dr6, dr8
    real(rk) :: r(0:300), rho(0:300)
    real(rk) :: cv(0:360), th(0:360)
    integer(ik) :: m, m2, i, j, k
    real(rk), intent(out) :: r2, r4, r6, r8
    
    h = 0.1d0
    rmax = 30.0d0
    m = int(rmax/h)
    m = m + mod(m,2)

    c = c0
    a = t/(4.d0*log(3.d0))
  
    k = 1
    drms = 1.0
    do while (drms >= 1.d-10)
       if(k>1) then
          c = trms/rms*c
       end if
       do i=0, m
          r(i) = h*i;
          if(abs(b20)<0.0001d0.and.abs(b40)<0.0001d0) then
             rho(i) = 4.0d0*pi*(1.0d0 + w*r(i)**2.0d0/c**2.0d0)/(1.0d0+exp((r(i)-c)/a))
          end if
       enddo

       if(abs(b20) > 0.0d0.or.abs(b40) > 0.0d0) then
          m2 = 360
          h2 = pi/360.d0
          do i=0, m2, 1
             th(i) = i*pi/360.d0
             cv(i) = c*(1.d0 + b20*sqrt(5.d0/(16.d0*pi))*(3.d0*(cos(th(i)))**2.d0 - 1.d0) &
                  + b40*3.d0/(16.d0*sqrt(pi))*(35.d0*(cos(th(i)))**4.d0 - 30.d0*(cos(th(i)))**2.d0 + 3.d0))
          end do
        
          do i=0, m, 1
             inte = 0.0
             do j=1, m2-1, 2
                dinte = h2/3.0d0 &
                     *((1.0d0 + w*r(i)**2.0d0/c**2.0d0)/(1.d0 + exp((r(i) - cv(j-1))/a))*sin(th(j-1)) &
                     + 4.0d0*(1.0d0 + w*r(i)**2.0d0/c**2.0d0)/(1.d0 + exp((r(i) - cv(j))/a))*sin(th(j)) &
                     + (1.0d0 + w*r(i)**2.0d0/c**2.0d0)/(1.d0 + exp((r(i) - cv(j+1))/a))*sin(th(j+1)))
                inte = inte + dinte
             end do
             rho(i) = 2.0d0*pi*inte   
          end do
       end if
       
       norm = 0.0d0
       r2 = 0.0d0
       r4 = 0.0d0
       r6 = 0.0d0
       r8 = 0.0d0
       do i=1, m-1, 2
          dnorm = (h/3.0d0)*(rho(i-1)*r(i-1)**2.0d0 + 4.0d0*rho(i)*r(i)**2.0d0 + rho(i+1)*r(i+1)**2.0d0)
          dr2 = (h/3.0d0)*(rho(i-1)*r(i-1)**4.0d0 + 4.0d0*rho(i)*r(i)**4.0d0 + rho(i+1)*r(i+1)**4.0d0)
          dr4 = (h/3.0d0)*(rho(i-1)*r(i-1)**6.0d0 + 4.0d0*rho(i)*r(i)**6.0d0 + rho(i+1)*r(i+1)**6.0d0) 
          dr6 = (h/3.0d0)*(rho(i-1)*r(i-1)**8.0d0 + 4.0d0*rho(i)*r(i)**8.0d0 + rho(i+1)*r(i+1)**8.0d0) 
          dr8 = (h/3.0d0)*(rho(i-1)*r(i-1)**10.0d0 + 4.0d0*rho(i)*r(i)**10.0d0 + rho(i+1)*r(i+1)**10.0d0) 
          norm = norm + dnorm
          r2 = r2 + dr2
          r4 = r4 + dr4
          r6 = r6 + dr6
          r8 = r8 + dr8
       end do
       r2 = r2/norm
       rms = sqrt(r2)
       r4 = r4/norm
       r6 = r6/norm
       r8 = r8/norm
       norm = dble(z)*1.0d0/norm
       
       !write(*,'(a4,i4,a16,f12.6)')'k: ',k,'rms: ', rms
       
       k = k + 1
       if(trms > 0.00001d0) then
          drms = abs(trms-rms)
       else
          drms = 0.0d0
       end if
    end do
    !rho(0:300) = norm*rho(0:300)/(4.d0*pi)
    
end subroutine radmoments
  
end module fermi
