program nuclear_moments_extraction

  use defs
  use solver
  use matrices
  
  implicit none

  real(rk) :: pseudo_dnu(2), er_order, F(2,4), concat_F(2,2), C(2,2), D(2,2)
  real(rk) :: r(2), sigma_f(2,2), y(2), r_new(2), sigma_f_new(2,2)
  integer(ik) :: A1, A2
  real(rk) :: A, dr2, dr4

  print*, 'Give the mass numbers of the reference and the target isotopes:'
  read(*,*) A1, A2 
  A = (A2+A1)/2
  
  print*, 'Give the computed pseudo - field shifts for transition 1:'
  read(*,*) pseudo_dnu(1)
  print*, 'and transition 2:'
  read(*,*) pseudo_dnu(2)

  print*, 'Give the dr^2 used to compute the pseudo - field shifts:'
  read(*,'(f7.4)') dr2 
  print*, 'and the dr^4:'
  read(*,'(f8.4)') dr4 

  print*, 'Give the order of magnitude of the errors in the'
  print*, 'pseudo - field shifts, e.g., 0.1, 0.01, ...'
  read(*,*) er_order
  
  ! Setting up the matrix containing the electronic factors
  ! Transition 1
  F(1,:) = (/-1.573512465813949D+05, 2.075124708925421D+02, &
       -0.545067800809402D+00, 0.0896047878353D-02 /)
  ! Transition 2
  F(2,:) = (/2.260087431698599D+05, -2.939843311017920D+02, &
       0.772387440181695D+00, -0.1266437615842D-02 /)

  call coefficient_matrices(F,A, concat_F,C,D)
  
  call matrix_equation(concat_F,pseudo_dnu, r)

  call error_estimation(pseudo_dnu,concat_F,er_order, sigma_f)

  call matrix_equation(C,pseudo_dnu, y)

  call matrix_equation(D,y, r_new)

  call error_estimation2(pseudo_dnu,C,D,er_order, sigma_f_new)
  
  !Printing the results on the screen
  write(*, '(a)')
  write(*,'(a75)') '---------------------------------------------------------------------------'
  write(*,'(a75)') '---------------------------------------------------------------------------'
  write(*, '(a)')
  write(*,'(a31,f7.4,a2,f6.4,a12,f8.4,a2,f7.4,a)') ' original  summation: <dr^2> = ',&
       r(1), ' (', sigma_f(1,1), '), <dr^4> = ', r(2), ' (', sigma_f(2,2), ')'
  write(*, '(a)')
  write(*,'(a31,f7.4,a2,f6.4,a12,f8.4,a2,f7.4,a)') 'rearranged summation: <dr^2> = ', &
       r_new(1),' (', sigma_f_new(1,1),'), <dr^4> = ', r_new(2),' (', sigma_f_new(2,2),')'
  write(*, '(a)')
  write(*,'(a75)') '---------------------------------------------------------------------------'
  write(*, '(a)')
  write(*, '(a31,f7.4,a20,f8.4)') '  exact   results  : <dr^2> = ', dr2, &
       ',          <dr^4> = ', dr4
  write(*, '(a)')
  write(*,'(a75)') '---------------------------------------------------------------------------' !maybe better to give the differenece , i.e., the model errors
  
end program nuclear_moments_extraction
