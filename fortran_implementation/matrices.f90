module matrices
  
  use defs

  implicit none

contains

  subroutine coefficient_matrices(el_factors,A,concatenated_F,C_matrix,D_matrix)

    real(rk), intent(in) :: el_factors(2,4), A
    real(rk), intent(out) :: concatenated_F(2,2), C_matrix(2,2), D_matrix(2,2)
    
    concatenated_F(1,1:2) = el_factors(1,1:2)
    concatenated_F(2,1:2) = el_factors(2,1:2)

    C_matrix(:,1) = 0.288554_rk*(A**(2._rk/3._rk))*el_factors(:,1) + &
         0.350673_rk*(A**(4._rk/3._rk))*el_factors(:,2) + &
         0.448303_rk*(A**2)*el_factors(:,3) + &
         0.592709_rk*(A**(8._rk/3._rk))*el_factors(:,4)
    C_matrix(:,2) = 0.0799258_rk*(A**(4._rk/3._rk))*el_factors(:,2) + &
         0.172916_rk*(A**2)*el_factors(:,3) + &
         0.2972_rk*(A**(8._rk/3._rk))*el_factors(:,4)

    D_matrix(1,:) = (/ 3.46556_rk/(A**(2._rk/3._rk)), 0._rk  /)
    D_matrix(2,:) = (/ -15.2051_rk/(A**(2._rk/3._rk)), &
         12.5116_rk/(A**(4._rk/3._rk)) /)
    
  end subroutine coefficient_matrices
    
end module matrices
