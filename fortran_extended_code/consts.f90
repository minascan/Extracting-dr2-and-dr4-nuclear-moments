module consts

  use defs

  integer(ik), parameter :: dim = 2     !Dimensions of coefficient matrix
  integer(ik), parameter :: l_dim = 2   !Leading dim of coefficient matrix
  integer(ik), parameter :: col_dnu = 1 !Columns of B --> dn matrix
  integer(ik), parameter :: row_dnu = 2 !Rows of B --> dn matrix

  real(rk), parameter :: c0 = 5.0
  real(rk), parameter :: t = 2.3
  real(rk), parameter :: w = 0.0
  
  !For a two-parameter spherical Fermi model
  real(rk), parameter :: b20 = 0.0
  real(rk), parameter :: b40 = 0.0
  
end module consts
