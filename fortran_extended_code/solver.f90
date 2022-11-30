module solver

  use defs
  use consts
  
  implicit none

  private :: inverse

contains
  
  ! Solution to a real system of linear equations A * X = B --> C * dy = dn
  subroutine matrix_equation(matrix,dnu, x)
      
    real(rk), intent(in)  :: matrix(2,2), dnu(2)
    real(rk), intent(out) :: x(2)
    real(rk) :: p_l_u(2,2)
    integer(ik) :: info, pivot(2)
    
    p_l_u = matrix
    x = dnu
    call dgesv(dim, col_dnu, p_l_u, l_dim, pivot, x, row_dnu, info) !lapack
    
    if (info /= 0) then
       print '(A)'
       print*, 'The coefficient matrix is numerically singular'
    end if
  
  end subroutine matrix_equation


  subroutine error_estimation(dnu,matrix,error_order, error)

    real(rk), intent(in)  :: dnu(2), matrix(2,2), error_order
    real(rk), intent(out) :: error(2,2)
    real(rk) :: inv_matrix(2,2), sigma_x(2,2), first_mul(2,2)
    integer :: i
    
    ! move this to matrices module
    sigma_x = 0._rk
    !forall(i = 1:row_dnu) sigma_x(i,i) = (dnu(i)*1.0D-03)**2
    forall(i = 1:row_dnu) sigma_x(i,i) = (dnu(i)*error_order)**2
    
    call inverse(matrix, inv_matrix)
    
    first_mul = matmul(inv_matrix,sigma_x)
    error = sqrt(matmul(first_mul,transpose(inv_matrix)))
    
  end subroutine error_estimation

  
  ! Error estimation for rearranged  expansion with orthonormal polynomials 
  subroutine error_estimation2(dnu,matrix1,matrix2,error_order, error_new)
    
    real(rk), intent(in)  :: dnu(2), matrix1(2,2), matrix2(2,2), error_order
    real(rk), intent(out) :: error_new(2,2)
    real(rk) :: inv_matrix1(2,2), inv_matrix2(2,2), inv_mul(2,2)
    real(rk) :: sec_mul(2,2), sigma_x(2,2)
    integer :: i
    
    ! move this to matrices module
    sigma_x = 0._rk
    !forall(i = 1:row_dnu) sigma_x(i,i) = (dnu(i)*1.0D-03)**2
    forall(i = 1:row_dnu) sigma_x(i,i) = (dnu(i)*error_order)**2
    
    call inverse(matrix1, inv_matrix1)
    call inverse(matrix2, inv_matrix2)

    inv_mul = matmul(inv_matrix2,inv_matrix1)
    sec_mul = matmul(inv_mul,sigma_x)
    error_new = sqrt(matmul(sec_mul,transpose(inv_mul)))
    
  end subroutine error_estimation2

  
  ! Private subroutine for analytically calculating the inverse of a 2x2 matrix
  subroutine inverse(A, B)
    
    real(rk), intent(in) :: A(2,2)
    real(rk), intent(out) :: B(2,2)
    real(rk) :: detinv
    
    ! Calculate the inverse determinant of the matrix
    detinv = 1/(A(1,1)*A(2,2) - A(1,2)*A(2,1))
    
    ! Calculate the inverse of the matrix
    B(1,1) = +detinv * A(2,2)
    B(2,1) = -detinv * A(2,1)
    B(1,2) = -detinv * A(1,2)
    B(2,2) = +detinv * A(1,1)
    
  end subroutine inverse
  
end module solver
