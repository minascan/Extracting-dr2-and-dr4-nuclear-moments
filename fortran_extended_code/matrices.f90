module matrices
  
  use defs

  implicit none

contains

  subroutine FS_data(A1_up1,A2_up1,A1_lo1,A2_lo1,A1_up2,A2_up2,A1_lo2,A2_lo2, dnu)

    real(rk), intent(in) :: A1_up1,A2_up1,A1_lo1,A2_lo1,A1_up2,A2_up2,A1_lo2,A2_lo2
    real(rk) :: DUP_1, DLO_1, DUP_2, DLO_2
    real(rk), intent(out) :: dnu(2)

    DUP_1 = (A1_up1 - A2_up1)*invcm_GHz
    DLO_1 = (A1_lo1 - A2_lo1)*invcm_GHz
    dnu(1) = DUP_1 - DLO_1
    DUP_2 = (A1_up2 - A2_up2)*invcm_GHz
    DLO_2 = (A1_lo2 - A2_lo2)*invcm_GHz
    dnu(2) = DUP_2 - DLO_2
    
  end subroutine FS_data

  subroutine scaling_F(F_init,QG_dr_2,QG_dr_4,QG_dr_6,QG_dr_8,QG_dnu, F_scaled)

    real(rk), intent(in) :: F_init(2,4), QG_dnu(2)
    real(rk), intent(in) :: QG_dr_2, QG_dr_4, QG_dr_6, QG_dr_8
    real(rk), intent(out) :: F_scaled(2,4)
    real(rk) :: RFS(2), scaling_factor(2)
    
    RFS(1) = F_init(1,1)*QG_dr_2 + F_init(1,2)*QG_dr_4 + F_init(1,3)*QG_dr_6 + F_init(1,4)*QG_dr_8
    RFS(2) = F_init(2,1)*QG_dr_2 + F_init(2,2)*QG_dr_4 + F_init(2,3)*QG_dr_6 + F_init(2,4)*QG_dr_8

    scaling_factor(:) = QG_dnu/RFS(:)
    
    F_scaled(1,:) = scaling_factor(1)*F_init(1,:)
    F_scaled(2,:) = scaling_factor(2)*F_init(2,:)
    
  end subroutine scaling_F
  
  subroutine coefficient_matrices(el_factors,A, concatenated_F,C_matrix,D_matrix)

    real(rk), intent(in) :: el_factors(2,4)
    real(rk), intent(in) :: A
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
