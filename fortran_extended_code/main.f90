program nuclear_moments_extraction

  use defs
  use fermi
  use matrices
  use solver
  
  implicit none
  logical :: reference
  integer(ik) :: Z, A1, A2
  character(30) :: upper_1, lower_1, upper_2, lower_2
  character(100) :: data, er_data
  real(rk) :: A1_upper_1, A2_upper_1, A1_lower_1, A2_lower_1 
  real(rk) :: A1_upper_2, A2_upper_2, A1_lower_2, A2_lower_2
  real(rk) :: QG_upper_1, QG_upper_2, QG_lower_1, QG_lower_2, QG_dnu(2) 
  real(rk) :: pseudo_dnu(2), er_order, F_in(2,4), F(2,4), concat_F(2,2), C(2,2), D(2,2)
  real(rk) :: r(2), sigma_f(2,2), y(2), r_new(2), sigma_f_new(2,2)
  real(rk) :: A, moder_dr2, moder_dr4, moder_dr2_new, moder_dr4_new
  real(rk) :: r_ch(3), dr_2, dr_4, dr_6, dr_8, QG_dr_2, QG_dr_4, QG_dr_6, QG_dr_8

  print*, '////////////////////////////////////////////////////////////////////////////////'
  print*, "//                          PROGRAM 'EXTRACT_DR2_DR4'                         //"
  print*, '////////////////////////////////////////////////////////////////////////////////'
  print*, '//                                                                            //'
  print*, '//Extraction of the differences in the second and forth nuclear radial moments//'
  print*, '//from available Field Shift (FS) data for two transitions in an isotope pair.//'
  print*, '//Both the original and the rearranged RFS expansions are used; the latter was//'
  print*, '//firstly introduced in: A. Papoulia, Masters thesis (2015).                  //'
  print*, '//The extraction was initially implemented in matlab by Carlsson B. G. (2015).//'
  print*, '//                                                                            //'
  print*, '//The electronic FS factors F resulting from RIS4 are scaled by comparing the //'
  print*, '//RFS with results from state-of-the-art calculations. For the scaling,the rms//'
  print*, '//radius of the target isotope is considered unknown and it is thus guessed.  //'
  print*, '//For the details of the method see Sec. V C 1 in:                            //'
  print*, '//A. Papoulia, B. G. Carlsson, and J. Ekman, Phys. Rev. A 94, 042502 (2016).  //'
  print*, '//                                                                            //'
  print*, '//                                Written by Asimina Papoulia, December 2020  //'
  print*, '//                                                                            //'
  print*, '////////////////////////////////////////////////////////////////////////////////'
  print '(a)'
  print*, 'Give the atomic number of the isotopes:'
  read(*,*) Z 
  print*, 'Give the mass number of the heavier isotope:'
  read(*,*) A1
  print*, 'Give the mass number of the lighter isotope:'
  read(*,*) A2 
  A = (A1+A2)/2
  print*, 'Is the reference isotope heavier than the target ? Give .true. or .false.'
  read(*,*) reference
  print '(a)'
  print*, '==============================================================================='
  print*, '                                    STEP 1:                                    '
  print*, '-------------------------------------------------------------------------------'
  print*, '  Getting the computed values of the energy levels to construct the PSEUDO-FS  '
  print*, '  data that will be used for extracting the dr^2 and dr^4 nuclear properties.  '
  !To construct these pseudo-data tabulated rms radii are often used, e.g., Angeli compilation
  print*, '==============================================================================='
  print '(a)'
  print*, '----------  Input data for TRANSITION 1 - Six entries are expected  -----------'
  print*, '(1) Give the configuration and term for the UPPER level:'
  read(*,*) upper_1  
  write(*,'(a59,i3,a9)') '(1a) Give the computed UPPER energy level in cm-1 for the ', A1, &
       ' isotope:'
  read(*,*) A1_upper_1 
  write(*,'(a59,i3,a9)') '(1b) Give the computed UPPER energy level in cm-1 for the ', A2, &
       ' isotope:'
  read(*,*) A2_upper_1
  print '(a)'
  print*, '(2) Give the configuration and term for the LOWER level:'
  read(*,*) lower_1 
  write(*,'(a59,i3,a9)') '(2a) Give the computed LOWER energy level in cm-1 for the ', A1, &
       ' isotope:'
  read(*,*) A1_lower_1 
  write(*,'(a59,i3,a9)') '(2b) Give the computed LOWER energy level in cm-1 for the ', A2, &
       ' isotope:'
  read(*,*) A2_lower_1 
  print '(a)'
  print*, '-------  Input data for TRANSITION 2 - Four to six entries are expected  ------'
  print*, '(1) Give the configuration and term for the UPPER level:'
  read(*,*) upper_2 
  if (upper_2.eq.upper_1) then
     A1_upper_2 = A1_upper_1
     A2_upper_2 = A2_upper_1
  elseif (upper_2.eq.lower_1) then
     A1_upper_2 = A1_lower_1
     A2_upper_2 = A2_lower_1
  else
     write(*,'(a59,i3,a9)') '(1a) Give the computed UPPER energy level in cm-1 for &
          the ', A1,' isotope:'
     read(*,*) A1_upper_2 
     write(*,'(a59,i3,a9)') '(1b) Give the computed UPPER energy level in cm-1 for &
          the ', A2,' isotope:'
     read(*,*) A2_upper_2 
  end if
  print '(a)'
  print*, '(2) Give configuration and term for the LOWER level:'
  read(*,*) lower_2 
  if (lower_2.eq.lower_1 .and. upper_2.eq.upper_1) then
     print*, 'Transition 2 is the SAME with Transition 1'
     stop
  elseif (lower_2.eq.lower_1) then
     A1_lower_2 = A1_lower_1
     A2_lower_2 = A2_lower_1
  elseif (lower_2.eq.upper_1 .and. upper_2.eq.lower_1) then
     print*, 'Transition 2 is the SAME with Transition 1'
     stop
  elseif (lower_2.eq.upper_1) then
     A1_lower_2 = A1_upper_1
     A2_lower_2 = A2_upper_1
  else
     write(*,'(a59,i3,a9)') '(2a) Give the computed LOWER energy level in cm-1 for the ', A1, &
          ' isotope:'
     read(*,*) A1_lower_2 
     write(*,'(a59,i3,a9)') '(2b) Give the computed LOWER energy level in cm-1 for the ', A2, &
          ' isotope:'
     read(*,*) A2_lower_2 
  end if

  call FS_data(A1_upper_1,A2_upper_1,A1_lower_1,A2_lower_1,A1_upper_2,A2_upper_2, &
       A1_lower_2,A2_lower_2, pseudo_dnu)
  
  print '(a)'
  print*, '==============================================================================='
  print*, '                                    STEP 2:                                    '
  print*, '-------------------------------------------------------------------------------'
  print*, '  Getting the computed values of the energy levels participating in the above  '
  print*, '      transitions 1 and 2 for the TARGET isotope and a GUESSED RMS RADIUS.     '
  print*, '==============================================================================='
  print '(a)'
  print*, '----------  Input data for TRANSITION 1 - Two entries are expected  -----------'
  write(*,'(a42,a30)') '(1) Give the computed energy in cm-1 for ', upper_1
  read(*,*) QG_upper_1  
  write(*,'(a42,a30)') '(2) Give the computed energy in cm-1 for ', lower_1
  read(*,*) QG_lower_1 
  print '(a)'
  print*, '----------  Input data for TRANSITION 2 - Two entries are expected - ----------'
  write(*,'(a42,a30)') '(1) Give the computed energy in cm-1 for ', upper_2
  read(*,*) QG_upper_2  
  write(*,'(a42,a30)') '(2) Give the computed energy in cm-1 for ', lower_2
  read(*,*) QG_lower_2 

  if (reference) then
     call FS_data(A1_upper_1,QG_upper_1,A1_lower_1,QG_lower_1,A1_upper_2,QG_upper_2, &
          A1_lower_2,QG_lower_2, QG_dnu)
  else
     call FS_data(QG_upper_1,A2_upper_1,QG_lower_1,A2_lower_1,QG_upper_2,A2_upper_2, &
          QG_lower_2,A2_lower_2, QG_dnu)
  end if

  print '(a)'
  print*, 'Give the TABULATED rms radius of the REFERENCE isotope:'
  read(*,*) r_ch(1)  
  print*, 'Give the TABULATED rms radius of the TARGET isotope:'
  read(*,*) r_ch(2) 
  print*, 'Give the GUESSED rms radius for the TARGET isotope:'
  read(*,*) r_ch(3) 
  !r_ch = (/ 5.8571d0, 5.8431d0, 5.7363d0 /)
  call radmoments_diffs(Z,r_ch, dr_2,dr_4,dr_6,dr_8,QG_dr_2,QG_dr_4,QG_dr_6,QG_dr_8)
  
  print '(a)'
  print*, '==============================================================================='
  print*, '                                    STEP 3:                                    '
  print*, '-------------------------------------------------------------------------------'
  print*, ' Rescaling the electronic factors F, given for the above transitions 1 and 2,  '
  print*, ' by comparing the computed FS values with RFS results. The rms radius of the   '
  print*, '  target isotope is guessed and the higher-order moments are deduced from a    '
  print*, '                  spherical two-parameter Fermi distribution.                  '
  print*, '==============================================================================='
  print '(a)'
    
  print*, 'Give the electronic factor FO (GHz fm-2) computed with RIS4 for Trans. 1'
  read(*,*) F_in(1,1)  
  print*, 'Give the electronic factor F2 (GHz fm-4) computed with RIS4 for Trans. 1'
  read(*,*) F_in(1,2)  
  print*, 'Give the electronic factor F4 (GHz fm-6) computed with RIS4 for Trans. 1'
  read(*,*) F_in(1,3)  
  print*, 'Give the electronic factor F6 (GHz fm-8) computed with RIS4 for Trans. 1'
  read(*,*) F_in(1,4) 
  print '(a)'
  print*, 'Give the electronic factor FO (GHz fm-2) computed with RIS4 for Trans. 2'
  read(*,*) F_in(2,1) 
  print*, 'Give the electronic factor F2 (GHz fm-4) computed with RIS4 for Trans. 2'
  read(*,*) F_in(2,2) 
  print*, 'Give the electronic factor F4 (GHz fm-6) computed with RIS4 for Trans. 2'
  read(*,*) F_in(2,3) 
  print*, 'Give the electronic factor F6 (GHz fm-8) computed with RIS4 for Trans. 2'
  read(*,*) F_in(2,4) 
 
  call scaling_F(F_in,QG_dr_2,QG_dr_4,QG_dr_6,QG_dr_8,QG_dnu, F)
  
  print '(a)'
  print*, '==============================================================================='
  print*, '                                    STEP 4:                                    '
  print*, '-------------------------------------------------------------------------------'
  print*, 'Solving the matrix equations that are constructed by 1) the original and 2) the'
  print*, '  re-arranged RFS expansions. By assuming the magnitude of the errors in the   '
  print*, 'pseudo-FS data, the uncertainties in the extracted dr^2 and dr^4 are computed. '
  print*, '==============================================================================='
  print '(a)'
 
  print*, 'Give the order of magnitude of the ERRORS in the PSEUDO-FS, e.g., 0.1, 0.01,...'
  read(*,*) er_order
  
  call coefficient_matrices(F,A, concat_F,C,D)
  
  call matrix_equation(concat_F,pseudo_dnu, r)

  call error_estimation(pseudo_dnu,concat_F,er_order, sigma_f)

  call matrix_equation(C,pseudo_dnu, y)

  call matrix_equation(D,y, r_new)

  call error_estimation2(pseudo_dnu,C,D,er_order, sigma_f_new)
  
  moder_dr2 = r(1)-dr_2
  moder_dr4 = r(2)-dr_4
  moder_dr2_new = r_new(1)-dr_2
  moder_dr4_new = r_new(2)-dr_4
  
  print '(a)'
  print*, '-------------------------------------------------------------------------------'
  print*, '----------------------------------- RESULTS -----------------------------------'
  print*, '-------------------------------------------------------------------------------'
  print '(a)'
  
  data = '(a31,f7.4,a2,f6.4,a12,f8.4,a2,f7.4,a)'
  er_data = '(a31,f7.4,a20,f8.4)'
  
  write(*, FMT=data) ' original  summation: <dr^2> = ',&
       r(1), ' (', sigma_f(1,1), '), <dr^4> = ', r(2), ' (', sigma_f(2,2), ')'
  print '(a)'
  write(*, FMT=er_data) '  model   errors  : <dr^2> = ', moder_dr2, &
       ',          <dr^4> = ', moder_dr4
  print '(a)'
  print*, '-------------------------------------------------------------------------------'
  print '(a)'
  write(*, FMT=data) 'rearranged summation: <dr^2> = ', &
       r_new(1),' (', sigma_f_new(1,1),'), <dr^4> = ', r_new(2),' (', sigma_f_new(2,2),')'
  print '(a)'
  write(*, FMT=er_data) '  model   errors  : <dr^2> = ', moder_dr2_new, &
       ',          <dr^4> = ', moder_dr4_new
  print '(a)'
  print*, '-------------------------------------------------------------------------------'
  
end program nuclear_moments_extraction
