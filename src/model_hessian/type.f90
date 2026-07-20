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

!> Common interface for model Hessian implementations.
module xtb_modelhessian_type
   use xtb_mctc_accuracy, only : wp
   use xtb_type_setvar, only : modhess_setvar
   implicit none
   private

   !> Abstract model Hessian implementation.
   type, public, abstract :: TModelHessian
   contains
      procedure(model_hessian_stretch), deferred :: stretch
      procedure(model_hessian_mode), deferred :: bend
      procedure(model_hessian_mode), deferred :: torsion
      procedure(model_hessian_mode), deferred :: outofplane
      procedure(model_hessian_charge), deferred :: add_charge
      !> Compute Hessian in packed lower-triangle storage.
      procedure :: compute_packed => compute_packed
      !> Compute dense symmetric Hessian.
      procedure :: compute_dense => compute_dense
      generic :: compute => compute_packed, compute_dense
   end type TModelHessian

   abstract interface
      subroutine model_hessian_stretch(self, xyz, n, hess, at, kr, kd, s6, &
            & lcutoff, rcut)
         import :: TModelHessian, wp
         class(TModelHessian), intent(in) :: self
         integer, intent(in) :: n
         real(wp), intent(in) :: xyz(3, n)
         real(wp), intent(inout) :: hess((3*n)*(3*n+1)/2)
         integer, intent(in) :: at(n)
         real(wp), intent(in) :: kr
         real(wp), intent(in) :: kd
         real(wp), intent(in) :: s6
         logical, intent(inout) :: lcutoff(n, n)
         real(wp), intent(in) :: rcut
      end subroutine model_hessian_stretch

      subroutine model_hessian_mode(self, xyz, n, hess, at, force_constant, &
            & kd, lcutoff)
         import :: TModelHessian, wp
         class(TModelHessian), intent(in) :: self
         integer, intent(in) :: n
         real(wp), intent(in) :: xyz(3, n)
         real(wp), intent(inout) :: hess((3*n)*(3*n+1)/2)
         integer, intent(in) :: at(n)
         real(wp), intent(in) :: force_constant
         real(wp), intent(in) :: kd
         logical, intent(in) :: lcutoff(n, n)
      end subroutine model_hessian_mode

      subroutine model_hessian_charge(self, xyz, n, hess, at, kq)
         import :: TModelHessian, wp
         class(TModelHessian), intent(in) :: self
         integer, intent(in) :: n
         real(wp), intent(in) :: xyz(3, n)
         real(wp), intent(inout) :: hess((3*n)*(3*n+1)/2)
         integer, intent(in) :: at(n)
         real(wp), intent(in) :: kq
      end subroutine model_hessian_charge
   end interface

contains

!> Compute Hessian in packed lower-triangle storage.
!> Element (i,j), j<=i, is stored at i*(i-1)/2+j.
subroutine compute_packed(self, xyz, n, hess, at, modh)
   class(TModelHessian), intent(in) :: self
   integer, intent(in) :: n
   real(wp), intent(in) :: xyz(3, n)
   real(wp), intent(out) :: hess((3*n)*(3*n+1)/2)
   integer, intent(in) :: at(n)
   type(modhess_setvar), intent(in) :: modh

   real(wp) :: kd
   logical, allocatable :: lcutoff(:, :)

   hess = 0.0_wp
   allocate(lcutoff(n, n), source=.false.)

   kd = modh%kd/modh%kr
   call self%stretch(xyz, n, hess, at, modh%kr, kd, modh%s6, lcutoff, modh%rcut)
   if (modh%kf /= 0.0_wp) then
      call self%bend(xyz, n, hess, at, modh%kf, kd, lcutoff)
   end if
   if (modh%kt /= 0.0_wp) then
      call self%torsion(xyz, n, hess, at, modh%kt, kd, lcutoff)
   end if
   if (modh%ko /= 0.0_wp) then
      call self%outofplane(xyz, n, hess, at, modh%ko, kd, lcutoff)
   end if
   if (modh%kq /= 0.0_wp) then
      call self%add_charge(xyz, n, hess, at, modh%kq)
   end if
end subroutine compute_packed

!> Compute dense symmetric Hessian from packed implementation.
subroutine compute_dense(self, xyz, n, hess, at, modh)
   class(TModelHessian), intent(in) :: self
   integer, intent(in) :: n
   real(wp), intent(in) :: xyz(3, n)
   real(wp), intent(out) :: hess(3*n, 3*n)
   integer, intent(in) :: at(n)
   type(modhess_setvar), intent(in) :: modh

   integer :: i, j, ij
   real(wp), allocatable :: packed(:)

   allocate(packed((3*n)*(3*n+1)/2))
   call self%compute_packed(xyz, n, packed, at, modh)

   ij = 0
   do i = 1, 3*n
      do j = 1, i
         ij = ij + 1
         hess(j, i) = packed(ij)
         hess(i, j) = packed(ij)
      end do
   end do
end subroutine compute_dense

end module xtb_modelhessian_type
