program nurbs_vertical_final
    use iso_fortran_env
    implicit none
    integer, parameter :: dp = real64
    integer :: i, j, ni, nj, np, p_idx, iter, max_iter
    integer :: ns, np_count, np_i, alloc_stat
    real(dp) :: u_val, v_val, s_u, s_v, dist, min_dist_2d
    real(dp) :: u_min, u_max, v_min, v_max, range_w
    real(dp), allocatable :: KU(:), KV(:), W(:,:), P(:,:,:)
    real(dp), allocatable :: targets(:,:)
    real(dp) :: current(3), residual_2d(2), p_surf(3)
    real(dp) :: J_2x2(2,2), delta(2)
    real(dp) :: first_section(3,3), last_section(3,3)  ! Store only first and last sections

    max_iter = 50

    ! 1. CARICAMENTO DATI
    open(10, file='superficie.in', status='old')
    read(10, *) ns
    close(10)

    ! Read the first and last sections
    allocate(first_section(3, np_count), last_section(3, np_count))
    open(10, file='superficie.in', position='append')
    read(10, *) np_i
    do i = 1, np_i
        read(10, *) targets(i,1), targets(i,2), targets(i,3)
    end do
    
    ! Read the first section (top)
    np_count = np_i
    allocate(first_section(3, np_count))
    do i = 1, np_count
        read(10, *) first_section(1:3, i)
    end do
    
    ! Read the last section (bottom)
    np_count = np_i
    allocate(last_section(3, np_count))
    do i = 1, np_count
        read(10, *) last_section(1:3, i)
    end do
    close(10)

    ! Read control point dimensions
    open(10, file='superficie.in', position='append')
    read(10, *) ni, nj
    close(10)

    ! Generate knot vectors for u and v directions
    allocate(KU(0:ni+3), KV(0:nj+3))
    call generate_knots(ni, 3, KU)
    call generate_knots(nj, 3, KV)

    ! Build the control point grid
    allocate(P(0:ni-1, 0:nj-1, 3), W(0:ni-1, 0:nj-1))
    allocate(targets(np_i, 3))
    
    ! Initialize the control points
    do i = 0, ni-1
        do j = 0, nj-1
            P(i,j,1) = first_section(1, j+1) + (i/real(ni-1, dp)) * (last_section(1, j+1) - first_section(1, j+1))
            P(i,j,2) = first_section(2, j+1) + (i/real(ni-1, dp)) * (last_section(2, j+1) - first_section(2, j+1))
            P(i,j,3) = first_section(3, j+1) + (i/real(ni-1, dp)) * (last_section(3, j+1) - first_section(3, j+1))
            W(i,j) = 1.0_dp
        end do
    end do

    ! 2. ELABORAZIONE PUNTI
    open(30, file='risultati.dat', status='replace')
    do p_idx = 1, np_i
        
        ! --- FASE 1: SEARCH FOR INITIAL GUESS ---
        u_min = 0.0_dp; u_max = 1.0_dp; v_min = 0.0_dp; v_max = 1.0_dp; range_w = 1.0_dp
        best_uv = [0.5_dp, 0.5_dp]

        do step = 1, 3
            min_dist_2d = 1e30_dp
            do iu = 0, 40
                s_u = u_min + (real(iu, dp) / 40.0_dp) * (range_w / real(40, dp))
                if (s_u < 0.0_dp .or. s_u > 1.0_dp) cycle
                do iv = 0, 40
                    s_v = v_min + (real(iv, dp) / 40.0_dp) * (range_w / real(40, dp))
                    if (s_v < 0.0_dp .or. s_v > 1.0_dp) cycle
                    call get_point_eval(ni, nj, KU, KV, P, W, s_u, s_v, p_surf)
                    dist = (targets(p_idx,1) - p_surf(1))**2 + (targets(p_idx,2) - p_surf(2))**2
                    if (dist < min_dist_2d) then
                        min_dist_2d = dist
                        best_uv = [s_u, s_v]
                    end if
                end do
            end do
            range_w = range_w * 0.5_dp
            u_min = max(0.0_dp, best_uv(1) - range_w)
            u_max = min(1.0_dp, best_uv(1) + range_w)
            v_min = max(0.0_dp, best_uv(2) - range_w)
            v_max = min(1.0_dp, best_uv(2) + range_w)
        end do
        
        ! --- FASE 2: NEWTON-RAPHSON (OPTIMIZATION) ---
        u_val = best_uv(1); v_val = best_uv(2)
        do iter = 1, max_iter
            call nurbs_basis(KU, 3, ni-1, u_val, NU, dNU)
            call nurbs_basis(KV, 3, nj-1, v_val, NV, dNV)
            call eval_rational_surface(P, W, NU, NV, current)
            
            residual_2d(1) = targets(p_idx,1) - current(1)
            residual_2d(2) = targets(p_idx,2) - current(2)
            
            if (sqrt(sum(residual_2d**2)) < 1e-10_dp) exit

            call eval_jacobian_2x2(P, W, NU, NV, dNU, dNV, current, J_2x2)
            call solve_linear_2x2(J_2x2, residual_2d, delta)
            
            u_val = max(0.0_dp, min(1.0_dp, u_val + delta(1)))
            v_val = max(0.0_dp, min(1.0_dp, v_val + delta(2)))
            
            if (abs(delta(1)) < 1e-15_dp .and. abs(delta(2)) < 1e-15_dp) exit
        end do
        
        write(30, '(6F16.8)') targets(p_idx,:), current
    end do
    close(30)

contains

    subroutine generate_knots(n_ctrl, deg, knots)
        integer, intent(in) :: n_ctrl, deg
        real(dp), intent(out) :: knots(0:)
        integer :: k
        k = 0
        do k = 0, n_ctrl + deg
            if (k <= deg) then
                knots(k) = 0.0_dp
            else if (k >= n_ctrl + deg - deg + 1) then
                knots(k) = 1.0_dp
            else
                knots(k) = real(k, dp) / real(n_ctrl, dp)
            end if
        end do
    end subroutine

    subroutine nurbs_basis(knots, deg, n_max, u, N, dN)
        integer, intent(in) :: deg, n_max
        real(dp), intent(in) :: knots(0:), u
        real(dp), intent(out) :: N(0:n_max), dN(0:n_max)
        real(dp) :: d1, d2
        
        ! Initialize basis functions
        N = 0.0_dp
        dN = 0.0_dp
        
        ! Set initial conditions for the recurrence
        do i = 0, n_max
            if (u >= knots(i) .and. u < knots(i+1)) then
                N(i) = 1.0_dp
            else
                N(i) = 0.0_dp
            end if
        end do
        
        ! Recurrence relation for basis functions
        do i = 0, n_max - deg
            d1 = knots(i+deg) - knots(i)
            d2 = knots(i+deg+1) - knots(i+deg)
            if (d1 > 1e-12_dp .and. d2 > 1e-12_dp) then
                N(i) = N(i) * (u - knots(i)) / d1
                N(i+deg+1) = N(i+deg+1) * (knots(i+deg+1) - u) / d2
            else if (d1 > 1e-12_dp) then
                N(i) = N(i) * (u - knots(i)) / d1
            else if (d2 > 1e-12_dp) then
                N(i+deg+1) = N(i+deg+1) * (knots(i+deg+1) - u) / d2
            end if
        end do
        
        ! Compute derivatives if needed
        if (deg > 0) then
            do i = 0, n_max - deg
                d1 = knots(i+deg) - knots(i)
                if (d1 > 1e-12_dp) then
                    dN(i) = deg * ( (u - knots(i)) / d1 * N(i) - (knots(i+deg) - u) / d1 * N(i+deg+1) )
                else
                    dN(i) = 0.0_dp
                end if
            end do
        end if
    end subroutine

    subroutine eval_rational_surface(P_in, W_in, NU_in, NV_in, S_out)
        real(dp), intent(in) :: P_in(0:, 0:, :), W_in(0:, 0:)
        real(dp), intent(in) :: NU_in(0:), NV_in(0:)
        real(dp), intent(out) :: S_out(3)
        real(dp) :: den_loc, w_loc
        integer :: i_l, j_l
        
        den_loc = 0.0_dp
        S_out = 0.0_dp
        
        do i_l = 0, size(NU_in)-1
            do j_l = 0, size(NV_in)-1
                w_loc = NU_in(i_l) * NV_in(j_l) * W_in(i_l, j_l)
                S_out = S_out + w_loc * P_in(i_l, j_l, :)
                den_loc = den_loc + w_loc
            end do
        end do
        
        ! Normalize by denominator
        if (abs(den_loc) > 1e-15_dp) then
            S_out = S_out / den_loc
        else
            S_out = 0.0_dp
        end if
    end subroutine

    subroutine eval_jacobian_2x2(P_in, W_in, NU_in, NV_in, dNU_in, dNV_in, S_curr, J_out)
        real(dp), intent(in) :: P_in(0:, 0:, :), W_in(0:, 0:)
        real(dp), intent(in) :: NU_in(0:), NV_in(0:), dNU_in(0:), dNV_in(0:)
        real(dp), intent(out) :: J_out(2,2)
        real(dp) :: den_l, du_den, dv_den, du_num, dv_num
        integer :: i_l, j_l
        
        ! Jacobian elements
        du_den = 0.0_dp
        dv_den = 0.0_dp
        du_num = 0.0_dp
        dv_num = 0.0_dp
        
        do i_l = 0, size(NU_in)-1
            do j_l = 0, size(NV_in)-1
                du_den = du_den + dNU_in(i_l) * NV_in(j_l) * W_in(i_l, j_l)
                dv_den = dv_den + dNV_in(j_l) * NV_in(i_l) * W_in(i_l, j_l)
                du_num = du_num + (dNU_in(i_l) * NV_in(j_l) * W_in(i_l, j_l)) * P_in(i_l, j_l, 1)
                dv_num = dv_num + (dNU_in(i_l) * NV_in(j_l) * W_in(i_l, j_l)) * P_in(i_l, j_l, 2)
            end do
        end do
        
        ! Compute the denominator for the Jacobian
        den_l = 0.0_dp
        do i_l = 0, size(NU_in)-1
            do j_l = 0, size(NV_in)-1
                den_l = den_l + NU_in(i_l) * NV_in(j_l) * W_in(i_l, j_l)
            end do
        end do
        
        ! Jacobian matrix
        J_out(1,1) = (du_num * den_l - S_curr(1) * du_den) / (den_l**2 + 1e-16_dp)
        J_out(1,2) = (dv_num * den_l - S_curr(1) * dv_den) / (den_l**2 + 1e-16_dp)
        J_out(2,1) = (du_num * den_l - S_curr(2) * du_den) / (den_l**2 + 1e-16_dp)
        J_out(2,2) = (dv_num * den_l - S_curr(2) * dv_den) / (den_l**2 + 2, 1e-16_dp)
    end subroutine

    subroutine solve_linear_2x2(A_in, b_in, x_out)
        real(dp), intent(in) :: A_in(2,2), b_in(2)
        real(dp), intent(out) :: x_out(2)
        real(dp) :: det_l
        det_l = A_in(1,1)*A_in(2,2) - A_in(1,2)*A_in(2,1)
        
        if (abs(det_l) < 1e-25_dp) then
            x_out = 0.0_dp
        else
            x_out(1) = (A_in(2,2)*b_in(1) - A_in(1,2)*b_in(2)) / det_l
            x_out(2) = (A_in(1,1)*b_in(2) - A_in(2,1)*b_in(1)) / det_l
        end if
    end subroutine

    subroutine get_point_eval(ni, nj, KU, KV, P, W, u, v, pt)
        integer, intent(in) :: ni, nj
        real(dp), intent(in) :: KU(0:), KV(0:), P(0:,0:,:), W(0:,0:)
        real(dp), intent(in) :: u, v
        real(dp), intent(out) :: pt(3)
        real(dp) :: l_N(0:ni-1), l_dN(0:ni-1), l_M(0:nj-1), l_dM(0:nj-1)
        integer :: i_l, j_l
        
        ! Evaluate the basis functions for u and v
        call nurbs_basis(KU, 3, ni-1, u, l_N, l_dN)
        call nurbs_basis(KV, 3, nj-1, v, l_M, l_dM)
        
        ! Evaluate the rational surface point
        call eval_rational_surface(P, W, l_N, l_M, pt)
    end subroutine

end program
