! This file is part of xtb.
!
! Copyright (C) 2026 Leopold M. Seidler
!
! xtb is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! xtb is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with xtb.  If not, see <https://www.gnu.org/licenses/>.

!> O1 numerical Hessian utilities
!> Ref: https://doi.org/10.1021/acs.jctc.5c01354
module xtb_o1numhess
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : autoaa, autorcm
   use xtb_mctc_blas, only : mctc_gemm, mctc_dot
   use xtb_mctc_lapack, only : lapack_syevr
   use xtb_type_environment, only : TEnvironment
   use xtb_type_molecule, only : TMolecule
   use xtb_freq_project, only : trproj
   use xtb_param_covalentrad, only : get_cov_rad
   use xtb_param_vdwradd3, only : getVanDerWaalsRadD3
   implicit none
   private

   public :: adj_list
   public :: gen_local_hessian, lr_loop, get_vdw_neighbor_list, &
      & gen_displdir, swart, find_projected_imag_modes

   type :: adj_list
      integer, allocatable :: neighbors(:)
   end type adj_list

   abstract interface
      subroutine matvec_operator(x, y, env, ctx)
         import :: wp, TEnvironment
         real(wp), intent(in) :: x(:)
         real(wp), intent(inout) :: y(:)
         type(TEnvironment), intent(inout) :: env
         class(*), optional, target, intent(inout) :: ctx
      end subroutine matvec_operator
   end interface

   type :: odlr_operator_data
      integer :: N, ndispl_final, info
      logical, allocatable :: mask(:, :)
      real(wp), pointer :: displdir(:, :) => null()
      real(wp), allocatable :: W2(:, :)

      real(wp), allocatable :: tmp(:, :)
      real(wp), allocatable :: tmp2(:, :)
      real(wp), allocatable :: f1(:, :)
      real(wp), allocatable :: f2(:, :)
   end type odlr_operator_data

   character(len=*), parameter :: source = 'xtb_o1numhess'

contains

!> Main routine to recover local Hessian
subroutine gen_local_hessian(env, ndispl_final, distmat, displdir, g, dmax, hess_out)
   !> Calculation environment
   type(TEnvironment), intent(inout) :: env
   !> Number of displacements
   integer, intent(in) :: ndispl_final
   !> Distance matrix between atoms
   real(wp), intent(in) :: distmat(:, :)
   !> Displacement directions
   real(wp), intent(in), target :: displdir(:, :)
   !> Gradient derivatives
   real(wp), intent(in) :: g(:, :)
   !> Maximum distance threshold
   real(wp), intent(in) :: dmax
   !> Output Hessian matrix
   real(wp), intent(out) :: hess_out(:, :)

   real(wp), parameter :: lam = 1.0e-2_wp, bet = 1.5_wp, ddmax = 5.0_wp
   
   ! Local work arrays
   type(odlr_operator_data) :: ctx
   real(wp), allocatable :: rhs(:, :), rhsv(:), sol(:)
   integer :: i, ndim, N, info
   logical :: terminate_run

   N = size(distmat, 1)

   ! Prepare cache
   ctx%N = N
   ctx%ndispl_final = ndispl_final
   allocate(ctx%f1(N, N), ctx%f2(N, N), ctx%tmp(N, N), ctx%tmp2(N, ndispl_final))
   ctx%displdir => displdir

   ! Calculate Regularization Term W2
   ctx%W2 = lam * max(0.0_wp, distmat(:, :) - dmax)**(2.0_wp * bet)

   ! Calculate rhs
   allocate(rhs(N, N))
   rhs = 0.0_wp
   call mctc_gemm(g(:, :ndispl_final), displdir(:, :ndispl_final), rhs, transb="t")

   ! Masks and Packing
   allocate(ctx%mask(N, N))
   ctx%mask = (distmat < (dmax + ddmax))
   do i = 2, N
      ctx%mask(i, 1:i-1) = .false.
   end do
   
   ! RHS Vector (b in Ax=b)
   rhsv = pack_sym(rhs, ctx%mask)
   ndim = size(rhsv)

   ! Solve
   allocate(sol(ndim))
   sol = 0.0_wp
   info = 0
   call cg(env, odlr_operator, ndim, rhsv, sol, info, ctx=ctx)
   call env%check(terminate_run)
   if (info == 1) then
      call env%warning("local hessian: CG failed to converge", source)
   else if (terminate_run) then
      return
   end if

   ! Recover Hessian from solution
   hess_out = unpack_sym(sol, ctx%mask, N)
end subroutine gen_local_hessian

!> Defines the linear operator for defining the optimizational
!> problem in local Hessian calculation
subroutine odlr_operator(x, y, env, ctx)
   real(wp), intent(in) :: x(:)
   real(wp), intent(inout) :: y(:)
   type(TEnvironment), intent(inout) :: env
   class(*), optional, target, intent(inout) :: ctx

   type(odlr_operator_data), pointer :: op_data

   if (present(ctx)) then
      select type(ctx)
         type is (odlr_operator_data)
            op_data => ctx
         class default
            call env%error("odlr_operator: invalid context", source)
            return
      end select
   else
      call env%error("odlr_operator: missing context", source)
      return
   end if

   op_data%tmp = unpack_sym(x, op_data%mask, op_data%N)

   call mctc_gemm(op_data%tmp, op_data%displdir(:, :op_data%ndispl_final), op_data%tmp2)
   call mctc_gemm(op_data%tmp2, op_data%displdir(:, :op_data%ndispl_final), op_data%f1, transb="t")

   op_data%f1 = (op_data%f1 + transpose(op_data%f1)) / 2.0_wp
   op_data%f2 = op_data%W2 * op_data%tmp

   y = pack_sym(op_data%f1 + op_data%f2, op_data%mask)
end subroutine odlr_operator

!> Generic Conjugate Gradient Solver
subroutine cg(env, operator, ndim, rhs, x, info, x0, ctx)
   type(TEnvironment), intent(inout) :: env
   procedure(matvec_operator) :: operator
   integer, intent(in) :: ndim
   real(wp), intent(in) :: rhs(:)
   real(wp), intent(inout) :: x(:)
   integer, intent(out) :: info
   real(wp), intent(in), optional :: x0(:)
   class(*), optional, target, intent(inout) :: ctx

   integer, parameter :: max_iter = 1000
   real(wp), parameter :: tol = 1.0e-14_wp
   real(wp), allocatable :: r(:), p(:), Ap(:)
   real(wp) :: alpha, beta, rs_old, rs_new
   integer :: k
   logical :: terminate_run

   info = 0
   allocate(r(ndim), p(ndim), Ap(ndim))

   if (present(x0)) then
      call operator(x, Ap, env, ctx)
      r = rhs - Ap
   else
      r = rhs
   end if
   p = r
   rs_old = mctc_dot(r, r)

   do k = 1, max_iter
      call operator(p, Ap, env, ctx)
      call env%check(terminate_run)
      if (terminate_run) return

      alpha = rs_old / mctc_dot(p, Ap)
      x = x + alpha * p
      r = r - alpha * Ap
      rs_new = mctc_dot(r, r)

      if (rs_new < tol) exit

      beta = rs_new / rs_old
      p = r + beta * p
      rs_old = rs_new
   end do

    if (k == max_iter) then
       info = 1
    end if
 end subroutine cg

!> Corrects Hessian hnum using a symmetric, low-rank update
!> so that g approx hnum * displdir
!> Uses column scaling (O1NumHess eqs. 17-18) to prioritize low-frequency modes
subroutine lr_loop(env, ndispl, g, hess_out, displdir, final_err)
   type(TEnvironment), intent(inout) :: env
   integer, intent(in) :: ndispl
   real(wp), intent(in) :: g(:, :)    ! Input Gradients
   real(wp), intent(inout) :: hess_out(:, :) ! Hessian to correct
   real(wp), intent(in) :: displdir(:, :)    ! Displacement directions
   real(wp), intent(out) :: final_err ! Final residual error (relative)

   real(wp), parameter :: scale_eps = 1.0e-3_wp
   real(wp), parameter :: thresh_LR = 1.0e-8_wp
   integer, parameter :: maxiter_LR = 200

   ! Local variables
   real(wp), allocatable :: resid(:, :), hcorr(:, :), tmp(:, :)
   real(wp), allocatable :: gscale(:), gscaled(:, :), xscaled(:, :)
   real(wp) :: rnorm, rnorm_prev, gnorm, relres
   integer :: it, N, j

   N = size(g, 1)

   ! Column scaling (eqs. 17-18): scale down large-gradient columns
   allocate(gscale(ndispl), gscaled(N, ndispl), xscaled(N, ndispl))
   do j = 1, ndispl
      gscale(j) = scale_eps / max(scale_eps, sqrt(sum(g(:, j)**2)))
      gscaled(:, j) = g(:, j) * gscale(j)
      xscaled(:, j) = displdir(:, j) * gscale(j)
   end do

   allocate(hcorr(N, N), tmp(N, N))
   gnorm = sqrt(sum(gscaled(:, :)**2))
   rnorm_prev = huge(1.0_wp)
   final_err = huge(1.0_wp)

   ! Iterative Correction Loop
   loop_lr: do it = 1, maxiter_LR

      call mctc_gemm(hess_out, xscaled(:, :ndispl), tmp)
      resid = gscaled(:, :ndispl) - tmp

      rnorm = sqrt(sum(resid(:, :)**2))
      if (gnorm > 0.0_wp) then
         relres = rnorm / gnorm
      else
         relres = rnorm
      end if
      final_err = relres

      if (rnorm < thresh_LR) then
         exit loop_lr
      else if (abs(rnorm - rnorm_prev) < thresh_LR) then
         exit loop_lr
      end if
      rnorm_prev = rnorm

      call mctc_gemm(resid, xscaled(:, :ndispl), hcorr, transb="t")
      hcorr = 0.5_wp * (hcorr + transpose(hcorr))
      hess_out = hess_out + hcorr

   end do loop_lr

   if (it == maxiter_LR) then
      call env%warning("LR loop failed to converge", source)
   end if

end subroutine lr_loop

!> Helper to unpack vector to Symmetric Matrix
function unpack_sym(v, mask, n) result(H)
   real(wp), intent(in) :: v(:)
   logical, intent(in) :: mask(n, n)
   integer, intent(in) :: n
   real(wp) :: H(n, n)
   integer :: i
   
   H = 0.0_wp
   H = unpack(v, mask, field=0.0_wp)
   
   ! Symmetrize
   do i = 2, n
      H(i, 1:i-1) = H(1:i-1, i)
   end do
end function unpack_sym

!> Helper to pack symmetric matrix
function pack_sym(m, mask) result(v)
   real(wp), intent(in) :: m(:, :)
   logical, intent(in) :: mask(:, :)
   real(wp), allocatable :: v(:)
   
   ! symmetrize, then pack
   v = pack((m + transpose(m)) * 0.5_wp, mask)
end function pack_sym

!> Setup a coordinate neighbor list from atom vdW radii.
subroutine get_vdw_neighbor_list(xyz, at, delta_r, nblist)
   real(wp), intent(in) :: xyz(:, :)
   integer, intent(in) :: at(:)
   real(wp), intent(in) :: delta_r
   type(adj_list), allocatable, intent(out) :: nblist(:)

   real(wp), allocatable :: rvdw(:)
   real(wp) :: rij, cutoff
   integer :: i, j, ic, jc, nat, n

   nat = size(at)
   n = 3 * nat
   allocate(nblist(n), rvdw(nat))

   do i = 1, nat
      rvdw(i) = getVanDerWaalsRadD3(at(i))
      if (rvdw(i) <= 0.0_wp) rvdw(i) = 2.0_wp
   end do

   do i = 1, nat
      do j = 1, nat
         if (i == j) then
            rij = 0.0_wp
         else
            rij = sqrt(sum((xyz(:, i) - xyz(:, j))**2))
         end if
         cutoff = rvdw(i) + rvdw(j) + delta_r
         if (rij <= cutoff) then
            do ic = 1, 3
               do jc = 1, 3
                  call add_neighbor(nblist(3 * (i - 1) + ic), 3 * (j - 1) + jc)
               end do
            end do
         end if
      end do
   end do
end subroutine get_vdw_neighbor_list

!> Helper to add neighbors to dynamic array
subroutine add_neighbor(list, val)
   type(adj_list), intent(inout) :: list
   integer, intent(in) :: val
   integer, allocatable :: tmp(:)
   integer :: sz

   if (.not. allocated(list%neighbors)) then
      allocate(list%neighbors(1))
      list%neighbors(1) = val
   else
      sz = size(list%neighbors)
      allocate(tmp(sz + 1))
      tmp(1:sz) = list%neighbors
      tmp(sz + 1) = val
      call move_alloc(tmp, list%neighbors)
   end if
end subroutine add_neighbor

!> Calculate displacement vectors for the gradient derivatives
subroutine gen_displdir(n, ndispl0, h0, max_nb, nblist, nbcounts, &
                        eps, eps2, displdir, ndispl_final, unit)
   integer, intent(in) :: n, ndispl0, max_nb
   !> Initial guess Hessian
   real(wp), intent(in) :: h0(n,n)
   !> Neighbor list
   type(adj_list), intent(in) :: nblist(:)
   !> Number of neighbors per atom
   integer, intent(in) :: nbcounts(n)
   !> Thresholds for sign determination and convergence
   real(wp), intent(in) :: eps, eps2
   !> Displacement directions
   real(wp), intent(inout) :: displdir(n, n)
   !> Final number of displacements
   integer, intent(out) :: ndispl_final
   !> Optional unit for profiling output
   integer, intent(in), optional :: unit

   ! Local variables
   integer :: i, j, k, p, q, nnb, info, n_curr, idx, locind, qn, m_eig
   integer :: nb_idx(max_nb)
   real(wp) :: ev(n), coverage(n), locev(max_nb)
   real(wp), allocatable :: eye(:, :)
   real(wp) :: loceigs(max_nb)
   real(wp) :: norm_ev1, norm_ev2, v_norm, d_dot, u2
   real(wp), allocatable :: locev_store(:, :)
   logical :: done
   real(wp), parameter :: orth_tol = 1.0e-10_wp
   integer :: prof_rate, prof_t0, prof_t1, prof_s0, prof_s1
   integer :: prof_evals, prof_sum_nnb, prof_max_nnb
   integer :: clk_a, clk_b, clk_c, clk_d, clk_e
   real(wp) :: prof_extract, prof_orth, prof_proj, prof_diag, prof_sign, prof_wall
   ! Thread-private work arrays for parallel region
   real(wp), allocatable :: submat_local(:,:), projmat_local(:,:), vec_subset_local(:,:), tmp_local(:,:), work_local(:)
   real(wp), allocatable :: eigvec_local(:, :)
   integer, allocatable :: iwork_local(:), isuppz_local(:)

   ! Initialize identity matrix
   allocate(eye(max_nb, max_nb))
   eye = 0.0_wp
   do i = 1, max_nb
      eye(i, i) = 1.0_wp
   end do

   ndispl_final = ndispl0
   call system_clock(prof_t0, prof_rate)
   prof_extract = 0.0_wp
   prof_orth = 0.0_wp
   prof_proj = 0.0_wp
   prof_diag = 0.0_wp
   prof_sign = 0.0_wp
   prof_evals = 0
   prof_sum_nnb = 0
   prof_max_nnb = 0

   ! Storage for eigenvectors computed in parallel
   allocate(locev_store(max_nb, n))

   done = .false.

   ! --- Open parallel region once ---
   !$omp parallel default(none) &
   !$omp shared(n, ndispl0, nblist, nbcounts, h0, displdir, max_nb, eye, locev_store, &
   !$omp        ev, coverage, done, ndispl_final, eps, eps2, prof_rate, prof_sign) &
   !$omp private(j, p, q, k, nnb, nb_idx, loceigs, locind, info, idx, n_curr, qn, m_eig, &
   !$omp         norm_ev1, norm_ev2, v_norm, d_dot, u2, locev, &
   !$omp         submat_local, projmat_local, vec_subset_local, tmp_local, work_local, &
   !$omp         eigvec_local, iwork_local, isuppz_local, &
   !$omp         prof_s0, prof_s1, clk_a, clk_b, clk_c, clk_d, clk_e) &
   !$omp reduction(+:prof_extract, prof_orth, prof_proj, prof_diag, prof_evals, prof_sum_nnb) &
   !$omp reduction(max:prof_max_nnb)

   allocate(submat_local(max_nb, max_nb))
   allocate(projmat_local(max_nb, n))
   allocate(vec_subset_local(max_nb, n))
   allocate(tmp_local(max_nb, max_nb))
   allocate(work_local(max(1, 200*max_nb)))
   allocate(iwork_local(max(1, 50*max_nb)))
   allocate(isuppz_local(max(1, 2*max_nb)))
   allocate(eigvec_local(max_nb, 1))

   ! --- Outer Loop: Generate new directions ---
   do n_curr = ndispl0, n - 1
      if (done) exit

      !$omp single
      ev = 0.0_wp
      coverage = 0.0_wp
      !$omp end single

      ! --- Parallel Phase: Compute local eigenvectors ---

      !$omp do schedule(runtime)
      do j = 1, n
         nnb = nbcounts(j)
         nb_idx(:nnb) = nblist(j)%neighbors(:nnb)

         ! Skip if subspace saturated
         if (nnb <= n_curr) cycle
         prof_evals = prof_evals + 1
         prof_sum_nnb = prof_sum_nnb + nnb
         prof_max_nnb = max(prof_max_nnb, nnb)
         call system_clock(clk_a)

         ! 1. Extract submatrix H0 (submat_p)
         do p = 1, nnb
            do q = 1, nnb
               submat_local(q, p) = h0(nb_idx(q), nb_idx(p))
            end do
         end do

         ! 2. Local Projection
         ! Form matrix A = displdir[neighbors, 0:n_curr]
         do p = 1, n_curr
            do q = 1, nnb
               vec_subset_local(q, p) = displdir(nb_idx(q), p)
            end do
         end do
         call system_clock(clk_b)

         if (n_curr > 0) then
            qn = 0
            do p = 1, n_curr
               locev(1:nnb) = vec_subset_local(1:nnb, p)
               do q = 1, qn
                  locev(1:nnb) = locev(1:nnb) &
                     & - dot_product(projmat_local(1:nnb, q), locev(1:nnb)) * projmat_local(1:nnb, q)
               end do
               u2 = dot_product(locev(1:nnb), locev(1:nnb))
               if (u2 < orth_tol) cycle
               qn = qn + 1
               if (qn > nnb) exit
               projmat_local(1:nnb, qn) = locev(1:nnb) / sqrt(u2)
            end do
         else
            qn = 0
         end if

         ! Build an orthonormal complement Z to Q and solve Z^T H Z.
         ! This avoids forming the dense projector P and diagonalizing null modes.
         locind = 0
         do p = 1, nnb
            locev(1:nnb) = 0.0_wp
            locev(p) = 1.0_wp
            do q = 1, qn
               locev(1:nnb) = locev(1:nnb) &
                  & - dot_product(projmat_local(1:nnb, q), locev(1:nnb)) * projmat_local(1:nnb, q)
            end do
            do q = 1, locind
               locev(1:nnb) = locev(1:nnb) &
                  & - dot_product(tmp_local(1:nnb, q), locev(1:nnb)) * tmp_local(1:nnb, q)
            end do
            u2 = dot_product(locev(1:nnb), locev(1:nnb))
            if (u2 < orth_tol) cycle
            locind = locind + 1
            tmp_local(1:nnb, locind) = locev(1:nnb) / sqrt(u2)
            if (locind >= nnb - qn) exit
         end do
         call system_clock(clk_c)

         if (locind < 1) then
            locev_store(1:nnb, j) = 0.0_wp
            cycle
         end if

         call mctc_gemm(submat_local(:nnb, :nnb), tmp_local(:nnb, :locind), vec_subset_local(:nnb, :locind))
         call mctc_gemm(tmp_local(:nnb, :locind), vec_subset_local(:nnb, :locind), &
            & projmat_local(:locind, :locind), transa="t")

         ! Symmetrize
         projmat_local(:locind, :locind) = 0.5_wp * (projmat_local(:locind, :locind) + transpose(projmat_local(:locind, :locind)))
         call system_clock(clk_d)

         ! 3. Largest-eigenpair diagonalization
         call lapack_syevr('V', 'I', 'U', locind, projmat_local, max_nb, 0.0_wp, 0.0_wp, locind, locind, &
            & 0.0_wp, m_eig, loceigs, eigvec_local, max_nb, isuppz_local, work_local, size(work_local), &
            & iwork_local, size(iwork_local), info)
         call system_clock(clk_e)
         prof_extract = prof_extract + real(clk_b - clk_a, wp) / real(prof_rate, wp)
         prof_orth = prof_orth + real(clk_c - clk_b, wp) / real(prof_rate, wp)
         prof_proj = prof_proj + real(clk_d - clk_c, wp) / real(prof_rate, wp)
         prof_diag = prof_diag + real(clk_e - clk_d, wp) / real(prof_rate, wp)
         
         ! Store the lifted eigenvector for the serial phase.
         locev_store(1:nnb, j) = 0.0_wp
         do p = 1, locind
            locev_store(1:nnb, j) = locev_store(1:nnb, j) + tmp_local(1:nnb, p) * eigvec_local(p, 1)
         end do
      end do
      !$omp end do

       ! --- Serial Phase: Sign fixing and accumulation ---
       !$omp single
       call system_clock(prof_s0)
       do j = 1, n
         nnb = nbcounts(j)
         nb_idx(:nnb) = nblist(j)%neighbors(:nnb)

         if (nnb <= n_curr) cycle

         ! Retrieve eigenvector
         locev(1:nnb) = locev_store(1:nnb, j)

         ! 4. Patching / Phase fixing
         ! Calculate norms for sign decision
         norm_ev1 = 0.0_wp
         norm_ev2 = 0.0_wp
         do p = 1, nnb
            idx = nb_idx(p)
            norm_ev1 = norm_ev1 + ((coverage(idx)*ev(idx) + locev(p))/(coverage(idx)+1.0_wp))**2
            norm_ev2 = norm_ev2 + ((coverage(idx)*ev(idx) - locev(p))/(coverage(idx)+1.0_wp))**2
         end do
         norm_ev1 = sqrt(norm_ev1)
         norm_ev2 = sqrt(norm_ev2)
         
         ! Apply update
         if (norm_ev1 > norm_ev2 + eps) then
            do p = 1, nnb
               idx = nb_idx(p)
               ev(idx) = (coverage(idx)*ev(idx) + locev(p))/(coverage(idx)+1.0_wp)
            end do
         else if (norm_ev1 < norm_ev2 - eps) then
            do p = 1, nnb
               idx = nb_idx(p)
               ev(idx) = (coverage(idx)*ev(idx) - locev(p))/(coverage(idx)+1.0_wp)
            end do
         else
            ! Deterministic sign fix based on max element
            locind = maxloc(abs(locev(1:nnb)), dim=1)
            if (locev(locind) > 0.0_wp) then
               do p = 1, nnb
                  idx = nb_idx(p)
                  ev(idx) = (coverage(idx)*ev(idx) + locev(p))/(coverage(idx)+1.0_wp)
               end do
            else
               do p = 1, nnb
                  idx = nb_idx(p)
                  ev(idx) = (coverage(idx)*ev(idx) - locev(p))/(coverage(idx)+1.0_wp)
               end do
            end if
         end if

         ! Update coverage
         do p = 1, nnb
            coverage(nb_idx(p)) = coverage(nb_idx(p)) + 1.0_wp
         end do
      end do ! End J loop (serial phase)

      ! Project out previous columns from global ev
      do k = 1, n_curr
         d_dot = mctc_dot(ev, displdir(:, k))
         ev = ev - d_dot * displdir(:, k)
      end do
      ! --- Check Norm ---
      v_norm = norm2(ev)
      
      if (v_norm < eps2) then
         done = .true.
      else
         ! Normalize and store
         ev = ev / v_norm
         displdir(:, n_curr + 1) = ev
         ndispl_final = n_curr + 1
      end if
      call system_clock(prof_s1)
      prof_sign = prof_sign + real(prof_s1 - prof_s0, wp) / real(prof_rate, wp)
      !$omp end single
   end do

   deallocate(submat_local, projmat_local, vec_subset_local, tmp_local, work_local, &
      & iwork_local, isuppz_local, eigvec_local)
   !$omp end parallel

   deallocate(eye, locev_store)
   call system_clock(prof_t1)
   prof_wall = real(prof_t1 - prof_t0, wp) / real(prof_rate, wp)
   if (present(unit)) then
      if (unit > 0) then
         write(unit,'("PROF gen_displdir: wall=",f10.3," s, evals=",i0,", avg_nnb=",f8.2,", max_nnb=",i0)') &
            & prof_wall, prof_evals, real(prof_sum_nnb, wp) / real(max(1, prof_evals), wp), prof_max_nnb
         write(unit,'("PROF gen_displdir: extract=",f10.3," s, orth=",f10.3," s, proj=",f10.3," s, diag=",f10.3," s, sign=",f10.3," s")') &
            & prof_extract, prof_orth, prof_proj, prof_diag, prof_sign
         flush(unit)
      end if
   end if
end subroutine gen_displdir

!> Calculates a modified Swart model Hessian
subroutine swart(env, xyz, at, hess_out)
   !> Computation environment
   type(TEnvironment), intent(inout) :: env
   !> coords
   real(wp), intent(in) :: xyz(:, :)
   !> ordinal numbers
   integer, intent(in) :: at(:)
   !> the full model hessian
   real(wp), intent(inout) :: hess_out(:, :)

   real(wp), parameter :: wthr = 0.3_wp, f = 0.12_wp, tolth = 0.2_wp, eps1 = wthr**2, eps2 = wthr**2 / exp(1.0_wp)

   real(wp) :: equildist, Hint, bmat6(6), bmat9(9), bmat29(2, 9), outer6(6, 6), outer9(9, 9), s_ijjk, costh, sinth, th1, scalelin
   real(wp), allocatable :: screenfunc(:, :), hess_local(:, :)
   real(wp) :: ri, rj
   integer :: i, j, k, nat, N, i1, i2, j1, j2, k1, k2
   
   nat = size(xyz, 2)
   N = 3*nat

   hess_out = 0.0_wp

   allocate(screenfunc(nat, nat))
   do i = 1, nat
      do j = i + 1, nat
         ri = get_cov_rad(at(i))
         rj = get_cov_rad(at(j))
         if (ri < 0.0_wp .or. rj < 0.0_wp) then
            call env%error("swart: covalent radii only defined for 1-103", source)
            return
         end if

         equildist = ri + rj
         screenfunc(i, j) = exp(1.0_wp - norm2(xyz(:, i) - xyz(:, j)) / equildist)
         screenfunc(j, i) = screenfunc(i, j)
      end do
   end do

   !$omp parallel default(none) &
   !$omp shared(nat, N, xyz, screenfunc, hess_out) private(i, j, k) &
   !$omp private(Hint, bmat6, bmat9, bmat29, outer6, outer9, s_ijjk, costh, sinth, th1, scalelin) &
   !$omp private(i1, i2, j1, j2, k1, k2, hess_local)
   allocate(hess_local(N, N))
   hess_local = 0.0_wp

   !$omp do schedule(runtime)
   do i = 1, nat
      do j = i + 1, nat
         Hint = 0.35_wp * screenfunc(i, j)**3
         bmat6 = bmat_bond(xyz(:, i) - xyz(:, j))
         outer6 = (spread(bmat6, dim=2, ncopies=6) * spread(bmat6, dim=1, ncopies=6))
         i1 = 3 * i - 2
         i2 = 3 * i
         j1 = 3 * j - 2
         j2 = 3 * j
         hess_local(i1:i2, i1:i2) = hess_local(i1:i2, i1:i2) + Hint * outer6(1:3, 1:3)
         hess_local(i1:i2, j1:j2) = hess_local(i1:i2, j1:j2) + Hint * outer6(1:3, 4:6)
         hess_local(j1:j2, i1:i2) = hess_local(j1:j2, i1:i2) + Hint * outer6(4:6, 1:3)
         hess_local(j1:j2, j1:j2) = hess_local(j1:j2, j1:j2) + Hint * outer6(4:6, 4:6)
      end do
   end do
   !$omp end do

   !$omp do schedule(runtime)
   do i = 1, nat
      do j = 1, nat
         if (i == j) cycle
         if (screenfunc(i, j) < eps2) cycle
         do k = i + 1, nat
            if (k == j) cycle
            s_ijjk = screenfunc(i, j) * screenfunc(j, k)
            if (s_ijjk < eps1) cycle

            costh = cosangle(xyz(:, i) - xyz(:, j), xyz(:, k) - xyz(:, j))
            sinth = sqrt(max(0.0_wp, 1.0_wp - costh**2))
            Hint = 0.075_wp * s_ijjk**2 * (f + (1 - f) * sinth)**2
            bmat9 = bmat_angle(xyz(:, i) - xyz(:, j), xyz(:, k) - xyz(:, j))

            if (costh > 1.0_wp - tolth) then
               th1 = 1.0_wp - costh
            else
               th1 = 1.0_wp + costh
            end if

            i1 = 3 * i - 2
            i2 = 3 * i
            j1 = 3 * j - 2
            j2 = 3 * j
            k1 = 3 * k - 2
            k2 = 3 * k
            if (th1 < tolth) then
               scalelin = (1.0_wp - (th1 / tolth)**2)**2
               if (costh > 1.0_wp - tolth) then
                  bmat29 = bmat_linangle(xyz(:, i) - xyz(:, j), xyz(:, k) - xyz(:, j))
                  bmat9 = scalelin * bmat29(1, :) + (1.0_wp - scalelin) * bmat9
                  outer9 = Hint * spread(bmat29(2, :), dim=2, ncopies=9) * spread(bmat29(2, :), dim=1, ncopies=9)
                  hess_local(i1:i2, i1:i2) = hess_local(i1:i2, i1:i2) + outer9(1:3, 1:3)
                  hess_local(i1:i2, j1:j2) = hess_local(i1:i2, j1:j2) + outer9(1:3, 4:6)
                  hess_local(i1:i2, k1:k2) = hess_local(i1:i2, k1:k2) + outer9(1:3, 7:9)
                  hess_local(j1:j2, i1:i2) = hess_local(j1:j2, i1:i2) + outer9(4:6, 1:3)
                  hess_local(j1:j2, j1:j2) = hess_local(j1:j2, j1:j2) + outer9(4:6, 4:6)
                  hess_local(j1:j2, k1:k2) = hess_local(j1:j2, k1:k2) + outer9(4:6, 7:9)
                  hess_local(k1:k2, i1:i2) = hess_local(k1:k2, i1:i2) + outer9(7:9, 1:3)
                  hess_local(k1:k2, j1:j2) = hess_local(k1:k2, j1:j2) + outer9(7:9, 4:6)
                  hess_local(k1:k2, k1:k2) = hess_local(k1:k2, k1:k2) + outer9(7:9, 7:9)
               else
                  bmat9 = (1.0_wp - scalelin) * bmat9
               end if
            end if

            outer9 = Hint * spread(bmat9, dim=2, ncopies=9) * spread(bmat9, dim=1, ncopies=9)
            hess_local(i1:i2, i1:i2) = hess_local(i1:i2, i1:i2) + outer9(1:3, 1:3)
            hess_local(i1:i2, j1:j2) = hess_local(i1:i2, j1:j2) + outer9(1:3, 4:6)
            hess_local(i1:i2, k1:k2) = hess_local(i1:i2, k1:k2) + outer9(1:3, 7:9)
            hess_local(j1:j2, i1:i2) = hess_local(j1:j2, i1:i2) + outer9(4:6, 1:3)
            hess_local(j1:j2, j1:j2) = hess_local(j1:j2, j1:j2) + outer9(4:6, 4:6)
            hess_local(j1:j2, k1:k2) = hess_local(j1:j2, k1:k2) + outer9(4:6, 7:9)
            hess_local(k1:k2, i1:i2) = hess_local(k1:k2, i1:i2) + outer9(7:9, 1:3)
            hess_local(k1:k2, j1:j2) = hess_local(k1:k2, j1:j2) + outer9(7:9, 4:6)
            hess_local(k1:k2, k1:k2) = hess_local(k1:k2, k1:k2) + outer9(7:9, 7:9)
         end do
      end do
   end do
   !$omp end do

   !$omp critical (swart_hess_merge)
   hess_out(:, :) = hess_out + hess_local
   !$omp end critical (swart_hess_merge)
   deallocate(hess_local)
   !$omp end parallel
end subroutine swart

!> Wilson B matrix for bonds
function bmat_bond(vec) result(bmat)
   real(wp), intent(in) :: vec(3)

   real(wp) :: l, bmat(6)
   
   bmat = 0.0_wp
   l = norm2(vec)

   bmat(1:3) = vec(:) / l
   bmat(4:6) = -vec(:) / l
end function bmat_bond

!> Wilson B matrix for non-linear angles
function bmat_angle(vec1, vec2) result(bmat)
   real(wp), intent(in) :: vec1(3), vec2(3)
   real(wp) :: bmat(9)
   real(wp) :: l1, l2, nvec1(3), nvec2(3)
   real(wp) :: dl(2, 6), dnvec(2, 3, 6), dinprod(9)
   real(wp) :: dot_n1n2
   integer :: ii

   l1 = norm2(vec1)
   l2 = norm2(vec2)
   nvec1 = vec1 / l1
   nvec2 = vec2 / l2

   dl = 0.0_wp
   dl(1, 1:3) = nvec1
   dl(1, 4:6) = -nvec1
   dl(2, 1:3) = nvec2
   dl(2, 4:6) = -nvec2

   dnvec = 0.0_wp
   do ii = 1, 6
      dnvec(1, 1:3, ii) = -nvec1 * dl(1, ii) / l1
      dnvec(2, 1:3, ii) = -nvec2 * dl(2, ii) / l2
   end do
   do ii = 1, 3
      dnvec(1, ii, ii) = dnvec(1, ii, ii) + 1.0_wp/l1
      dnvec(2, ii, ii) = dnvec(2, ii, ii) + 1.0_wp/l2
      dnvec(1, ii, ii+3) = dnvec(1, ii, ii+3) - 1.0_wp/l1
      dnvec(2, ii, ii+3) = dnvec(2, ii, ii+3) - 1.0_wp/l2
   end do

   dinprod = 0.0_wp
   do ii = 1, 3
      dinprod(ii) = mctc_dot(dnvec(1, :, ii), nvec2)
      dinprod(ii+3) = mctc_dot(dnvec(1, :, ii+3), nvec2) + mctc_dot(dnvec(2, :, ii+3), nvec1)
      dinprod(ii+6) = mctc_dot(dnvec(2, :, ii), nvec1)
   end do

   dot_n1n2 = mctc_dot(nvec1, nvec2)
   bmat = -dinprod / sqrt(max(1.0e-15_wp, 1.0_wp - dot_n1n2**2))
end function bmat_angle

!> Wirson B matrix for linear angles
function bmat_linangle(vec1, vec2) result(bmat)
   real(wp), intent(in) :: vec1(3), vec2(3)
   real(wp) :: bmat(2,9)
   real(wp) :: l1, l2, nvec1(3), nvec2(3)
   real(wp) :: vn(3), vn2(3), nvn
   real(wp), parameter :: xaxis(3) = [1.0_wp, 0.0_wp, 0.0_wp], yaxis(3) = [0.0_wp, 1.0_wp, 0.0_wp]

   l1 = norm2(vec1)
   l2 = norm2(vec2)
   nvec1 = vec1 / l1
   nvec2 = vec2 / l2

   vn(1) = vec1(2) * vec2(3) - vec1(3) * vec2(2)
   vn(2) = vec1(3) * vec2(1) - vec1(1) * vec2(3)
   vn(3) = vec1(1) * vec2(2) - vec1(2) * vec2(1)
   nvn = norm2(vn)

   if (nvn < 1.0e-15_wp) then
      vn = xaxis - mctc_dot(xaxis, vec1) / l1**2 * vec1
      nvn = norm2(vn)
      if (nvn < 1.0e-15_wp) then
         vn = yaxis - mctc_dot(yaxis, vec1) / l1**2 * vec1
         nvn = norm2(vn)
      end if
   end if
   vn = vn / nvn

   vn2(1) = (vec1(2) - vec2(2)) * vn(3) - (vec1(3) - vec2(3)) * vn(2)
   vn2(2) = (vec1(3) - vec2(3)) * vn(1) - (vec1(1) - vec2(1)) * vn(3)
   vn2(3) = (vec1(1) - vec2(1)) * vn(2) - (vec1(2) - vec2(2)) * vn(1)
   vn2 = vn2 / norm2(vn2)

   bmat = 0.0_wp
   bmat(2, 1:3) = vn / l1
   bmat(2, 7:9) = vn / l2
   bmat(2, 4:6) = -bmat(2, 1:3) - bmat(2, 7:9)
   bmat(1, 1:3) = vn2 / l1
   bmat(1, 7:9) = vn2 / l2
   bmat(1, 4:6) = -bmat(1, 1:3) - bmat(1, 7:9)
end function bmat_linangle

!> cos(Angle) between two vectors
function cosangle(vec1, vec2) result(cos_theta)
    implicit none
    real(wp), intent(in) :: vec1(3), vec2(3)
    real(wp) :: cos_theta
    
    cos_theta = mctc_dot(vec1, vec2) / (norm2(vec1) * norm2(vec2))
end function cosangle

!> Find significant imaginary modes of a Cartesian Hessian using the same
!> projection (trproj, unit-mass geometric center) and mass weighting as the
!> vibrational frequency path in src/hessian.F90. Returns Cartesian
!> displacement directions and frequencies in cm^-1 (negative = imaginary).
!> Linear molecules skip projection (matching src/hessian.F90); caller passes
!> the linear flag so the 3x3 inertia diag is not repeated here.
subroutine find_projected_imag_modes(env, mol, hess, linear, max_modes, modes, freqs, nmodes)
   type(TEnvironment), intent(inout) :: env
   type(TMolecule), intent(in) :: mol
   real(wp), intent(in) :: hess(:, :)
   logical, intent(in) :: linear
   integer, intent(in) :: max_modes
   real(wp), intent(out) :: modes(:, :)
   real(wp), intent(out) :: freqs(:)
   integer, intent(out) :: nmodes

   character(len=*), parameter :: source = "find_projected_imag_modes"
   integer :: N, nat, i, j, k, ia, ic, ii, m_eig, info, ndiag
   real(wp), allocatable :: hpack(:), Hmw(:, :), inv_sqrt_m(:), v(:), &
      & eigvec(:, :), work(:), dummy_mode(:, :)
   real(wp) :: eigval(max_modes)
   integer, allocatable :: isuppz(:), iwork(:)
   real(wp) :: freq
   real(wp), parameter :: neg_sig_min_rcm = 5.0_wp

   nat = mol%n
   N = 3 * nat
   nmodes = 0
   if (max_modes <= 0) return

   allocate(inv_sqrt_m(N), v(N))
   do ia = 1, nat
      do ic = 1, 3
         ii = 3 * (ia - 1) + ic
         inv_sqrt_m(ii) = 1.0_wp / sqrt(mol%atmass(ia))
      end do
   end do

   ! pack upper-by-column (j<=i), same layout as src/hessian.F90 trproj input
   allocate(hpack(N * (N + 1) / 2))
   k = 0
   do i = 1, N
      do j = 1, i
         k = k + 1
         hpack(k) = hess(j, i)
      end do
   end do

   ! project translations/rotations (skip for linear, matching src/hessian.F90)
   if (.not. linear) then
      allocate(dummy_mode(N, 1))
      dummy_mode = 0.0_wp
      call trproj(nat, N, mol%xyz, hpack, .false., 0, dummy_mode, 1)
      deallocate(dummy_mode)
   end if

   ! unpack and mass-weight: Hmw(i,j) = H(i,j) / sqrt(m_au(i) * m_au(j))
   allocate(Hmw(N, N))
   Hmw = 0.0_wp
   k = 0
   do i = 1, N
      do j = 1, i
         k = k + 1
         Hmw(j, i) = hpack(k)
         Hmw(i, j) = hpack(k)
      end do
   end do
   do i = 1, N
      do j = 1, N
         Hmw(i, j) = Hmw(i, j) * inv_sqrt_m(i) * inv_sqrt_m(j)
      end do
   end do

   ! lowest eigenpairs only
   ndiag = min(max_modes, N)
   allocate(eigvec(N, ndiag), isuppz(2 * ndiag))
   allocate(work(max(1, 200 * N)), iwork(max(1, 50 * N)))
   call lapack_syevr('V', 'I', 'U', N, Hmw, N, 0.0_wp, 0.0_wp, 1, ndiag, 0.0_wp, &
      & m_eig, eigval, eigvec, N, isuppz, work, size(work), iwork, size(iwork), info)
   if (info /= 0) then
      call env%warning("find_projected_imag_modes: dsyevr failed", source)
      nmodes = 0
      return
   end if

   do i = 1, m_eig
      freq = autorcm * sign(sqrt(abs(eigval(i))), eigval(i))
      if (freq < -neg_sig_min_rcm) then
         if (nmodes >= max_modes) exit
         nmodes = nmodes + 1
         freqs(nmodes) = freq
         v = eigvec(:, i) * inv_sqrt_m
         modes(:, nmodes) = v / norm2(v)
      end if
   end do
end subroutine find_projected_imag_modes

end module xtb_o1numhess
