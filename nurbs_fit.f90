program nurbs_robust_fixed
    use iso_fortran_env
    implicit none

    interface
        subroutine dgels(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
            import :: real64
            character :: trans
            integer :: m, n, nrhs, lda, ldb, lwork, info
            real(real64) :: a(lda, *), b(ldb, *), work(*)
        end subroutine dgels
    end interface

    integer, parameter :: dp = real64
    integer :: ns, ni, nj, deg, total_pts, i, j, k, r, c, info, idx
    real(dp) :: z_min, z_max, x_min, x_max
    real(dp), allocatable :: x_all(:), y_all(:), z_all(:), u_all(:), v_all(:)
    real(dp), allocatable :: KU(:), KV(:), A(:,:), B_vec(:,:), work(:)
    real(dp) :: N_val(0:10), M_val(0:10), lwork_opt(1)
    real(dp) :: tx, ty, tz, s_y, u_ev, v_ev

    ! --- 1. CONFIGURAZIONE ULTRA-STABILE ---
    deg = 2
    ni  = 4   ! Pochissimi punti = superficie tesa, impossibile da aggrovigliare
    nj  = 4
    z_min = 715.0_dp
    z_max = 1010.0_dp

    open(10, file='sezioni.dat', status='old')
    read(10, *) ns
    
    ! Primo passaggio: limiti X e conteggio
    total_pts = 0
    x_min = 1e10; x_max = -1e10
    do i = 1, ns
        read(10, *) k
        do j = 1, k
            read(10, *) tx, ty, tz
            if (tz >= z_min .and. tz <= z_max) then
                total_pts = total_pts + 1
                x_min = min(x_min, tx); x_max = max(x_max, tx)
            end if
        end do
    end do
    
    rewind(10); read(10, *)
    allocate(x_all(total_pts), y_all(total_pts), z_all(total_pts), u_all(total_pts), v_all(total_pts))

    ! Caricamento con parametrizzazione normalizzata
    r = 0
    do j = 1, ns
        read(10, *) k
        do i = 1, k
            read(10, *) tx, ty, tz
            if (tz >= z_min .and. tz <= z_max) then
                r = r + 1
                x_all(r) = tx; y_all(r) = ty; z_all(r) = tz
                ! Parametrizzazione robusta
                v_all(r) = (tx - x_min) / (x_max - x_min + 1e-8)
                u_all(r) = (tz - z_min) / (z_max - z_min + 1e-8)
            end if
        end do
    end do
    close(10)

    ! --- 2. BASIS E MATRICE ---
    allocate(KU(0:ni+deg), KV(0:nj+deg))
    call gen_clamped_knots(ni, deg, KU)
    call gen_clamped_knots(nj, deg, KV)

    allocate(A(total_pts, ni*nj), B_vec(max(total_pts, ni*nj), 1))
    A = 0.0_dp
    do k = 1, total_pts
        call bspline_basis(KU, deg, ni-1, u_all(k), N_val(0:ni-1))
        call bspline_basis(KV, deg, nj-1, v_all(k), M_val(0:nj-1))
        do j = 1, nj
            do i = 1, ni
                idx = (j-1)*ni + i
                A(k, idx) = N_val(i-1) * M_val(j-1)
            end do
        end do
        B_vec(k, 1) = y_all(k)
    end do

    ! --- 3. RISOLUZIONE ---
    call dgels('N', total_pts, ni*nj, 1, A, total_pts, B_vec, max(total_pts, ni*nj), lwork_opt, -1, info)
    allocate(work(nint(lwork_opt(1))))
    call dgels('N', total_pts, ni*nj, 1, A, total_pts, B_vec, max(total_pts, ni*nj), work, size(work), info)

    ! --- 4. OUTPUT ---
    open(40, file='superficie_filt.dat')
    do j = 0, 40
        v_ev = real(j, dp) / 40.0_dp
        do i = 0, 40
            u_ev = real(i, dp) / 40.0_dp
            call bspline_basis(KU, deg, ni-1, u_ev, N_val(0:ni-1))
            call bspline_basis(KV, deg, nj-1, v_ev, M_val(0:nj-1))
            s_y = 0.0_dp
            do r = 1, nj
                do c = 1, ni
                    idx = (r-1)*ni + c
                    s_y = s_y + N_val(c-1) * M_val(r-1) * B_vec(idx, 1)
                end do
            end do
            write(40, '(3F16.8)') x_min + v_ev*(x_max-x_min), s_y, z_min + u_ev*(z_max-z_min)
        end do
        write(40, *)
    end do
    close(40)

contains

    subroutine gen_clamped_knots(n, p, knots)
        integer, intent(in) :: n, p
        real(dp), intent(out) :: knots(0:)
        integer :: k, m
        m = n + p
        do k = 0, m
            if (k <= p) then
                knots(k) = 0.0_dp
            else if (k >= m - p) then
                knots(k) = 1.0_dp
            else
                knots(k) = real(k-p, dp) / real(n-p, dp)
            end if
        end do
    end subroutine

    subroutine bspline_basis(knots, p, n_max, u, N)
        integer, intent(in) :: p, n_max
        real(dp), intent(in) :: knots(0:), u
        real(dp), intent(out) :: N(0:n_max)
        real(dp) :: left(p+1), right(p+1), saved, temp, uf
        integer :: j, r, idx
        N = 0.0_dp; uf = max(0.0_dp, min(1.0_dp, u))
        idx = p ! Default safe
        do j = p, n_max
            if (uf >= knots(j) .and. uf <= knots(j+1)) then
                idx = j
                exit
            end if
        end do
        N(idx) = 1.0_dp
        do j = 1, p
            left(j) = uf - knots(idx+1-j)
            right(j) = knots(idx+j) - uf
            saved = 0.0_dp
            do r = 0, j-1
                temp = N(idx-r) / (right(r+1) + left(j-r) + 1e-15)
                N(idx-r) = saved + right(r+1) * temp
                saved = left(j-r) * temp
            end do
            N(idx-j) = saved
        end do
    end subroutine
end program
