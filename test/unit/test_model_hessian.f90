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

!> Behavior-lock tests for model Hessian variants.
!>
!> Reference data lives in `test/unit/fixtures/model_hessian/*.dat` as one
!> packed Hessian value per line. With GEN_REFS=1, the tests write files from
!> current code output and PASS. Otherwise they compare element-wise against
!> stored references.
!>
!> Tests may encode current buggy behavior (rkl2=sum(rjk**2) in torsion,
!> outofp2(xyz,...) in out-of-plane) — marked as behavior-lock until
!> Phase 3 fixes them, at which point references are regenerated.
module test_model_hessian
   use testdrive, only : new_unittest, unittest_type, error_type, check
   use xtb_mctc_accuracy, only : wp
   use xtb_type_molecule, only : TMolecule
   use xtb_modelhessian_lindh, only : mh_lindh, mh_lindh_d2
   use xtb_modelhessian_swart, only : mh_swart
   use xtb_type_setvar, only : modhess_setvar
   use xtb_o1numhess, only : swart
   use xtb_test_molstock, only : getMolecule
   use xtb_type_environment, only : TEnvironment, init
   implicit none
   private
   public :: collect_model_hessian

   integer, parameter :: VAR_LINDH_D2 = 1
   integer, parameter :: VAR_LINDH    = 2
   integer, parameter :: VAR_SWART    = 3

   real(wp), parameter :: thr = 10*epsilon(0.0_wp)

contains

!> Collect all exported unit tests
subroutine collect_model_hessian(testsuite)
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      ! Lindh-D2 (1995 tables), default constants
      new_unittest("lindh_d2_h2o", test_lindh_d2_h2o), &
      new_unittest("lindh_d2_mindless01", test_lindh_d2_mindless01), &
      new_unittest("lindh_d2_caffeine", test_lindh_d2_caffeine), &
      new_unittest("lindh_d2_mgh2", test_lindh_d2_mgh2), &
      ! Lindh 2007, default constants
      new_unittest("lindh_h2o", test_lindh_h2o), &
      new_unittest("lindh_mindless01", test_lindh_mindless01), &
      new_unittest("lindh_caffeine", test_lindh_caffeine), &
      new_unittest("lindh_mgh2", test_lindh_mgh2), &
      ! Swart, default constants
      new_unittest("swart_h2o", test_swart_h2o), &
      new_unittest("swart_mindless01", test_swart_mindless01), &
      new_unittest("swart_caffeine", test_swart_caffeine), &
      new_unittest("swart_mgh2", test_swart_mgh2), &
      ! Modified Swart (O1NumHess variant), hardcoded H2O reference
      new_unittest("modified_swart_h2o", test_modified_swart_h2o), &
      ! Out-of-plane behavior (ko != 0)
      new_unittest("lindh_d2_caffeine_oop", test_lindh_d2_caffeine_oop), &
      new_unittest("lindh_caffeine_oop", test_lindh_caffeine_oop), &
      new_unittest("swart_caffeine_oop", test_swart_caffeine_oop) &
      ]

end subroutine collect_model_hessian


!> Default modhess_setvar matching setparam.f90 defaults
function default_modh() result(modh)
   type(modhess_setvar) :: modh
   modh = modhess_setvar(kr=0.4000_wp, kf=0.1300_wp, kt=0.0075_wp, &
      & ko=0.0_wp, kd=0.0_wp, kq=0.0_wp, rcut=70.0_wp, s6=20.0_wp)
end function default_modh


!> modhess_setvar with out-of-plane force constant active
function oop_modh() result(modh)
   type(modhess_setvar) :: modh
   modh = modhess_setvar(kr=0.4000_wp, kf=0.1300_wp, kt=0.0075_wp, &
      & ko=0.16_wp, kd=0.0_wp, kq=0.0_wp, rcut=70.0_wp, s6=20.0_wp)
end function oop_modh

!> Generic model Hessian test driver.
!> Computes Hessian, compares against reference.
subroutine test_mh(error, molname, variant, modh, label)
   type(error_type), allocatable, intent(out) :: error
   character(len=*), intent(in) :: molname
   character(len=*), intent(in) :: label
   integer, intent(in) :: variant
   type(modhess_setvar), intent(in) :: modh

   type(TMolecule) :: mol
   real(wp), allocatable :: hess_packed(:)

   call getMolecule(mol, molname)
   call compute_mh_packed(mol, variant, modh, hess_packed)
   call compare_or_write_ref(error, label, hess_packed)
end subroutine test_mh

!> Compute model Hessian, return packed array
subroutine compute_mh_packed(mol, variant, modh, hess_packed)
   type(TMolecule), intent(in) :: mol
   integer, intent(in) :: variant
   type(modhess_setvar), intent(in) :: modh
   real(wp), allocatable, intent(out) :: hess_packed(:)

   integer :: n3

   n3 = 3 * mol%n
   allocate(hess_packed(n3*(n3+1)/2))
   hess_packed = 0.0_wp
   select case(variant)
   case(VAR_LINDH_D2)
      call mh_lindh_d2(mol%xyz, mol%n, hess_packed, mol%at, modh)
   case(VAR_LINDH)
      call mh_lindh(mol%xyz, mol%n, hess_packed, mol%at, modh)
   case(VAR_SWART)
      call mh_swart(mol%xyz, mol%n, hess_packed, mol%at, modh)
   end select
end subroutine compute_mh_packed


!> Path to reference file for a given label.
function ref_path(label) result(path)
   character(len=*), intent(in) :: label

   character(len=:), allocatable :: path
   character(len=512) :: dir
   integer :: stat

   call get_environment_variable("TEST_FIXTURES_DIR", dir, status=stat)
   if (stat == 0 .and. len_trim(dir) > 0) then
      path = trim(dir) // "/model_hessian/" // trim(label) // ".dat"
   else
      path = "test/unit/fixtures/model_hessian/" // trim(label) // ".dat"
   end if
end function ref_path


!> True when GEN_REFS is set to 1 (force regeneration of files).
function gen_refs() result(mode)
   logical :: mode
   character(len=8) :: buf
   integer :: stat
   call get_environment_variable("GEN_REFS", buf, status=stat)
   mode = (stat == 0 .and. trim(buf) == "1")
end function gen_refs


!> Write packed Hessian as one value per line to a file.
subroutine write_ref(path, packed)
   character(len=*), intent(in) :: path
   real(wp), intent(in) :: packed(:)

   integer :: u, i

   open(newunit=u, file=path, status="replace", action="write")
   do i = 1, size(packed)
      write(u, '(ES25.17E3)') packed(i)
   end do
   close(u)
end subroutine write_ref


!> Read reference packed Hessian from file. Returns .false. on read error.
function read_ref(path, packed) result(ok)
   character(len=*), intent(in) :: path
   real(wp), allocatable, intent(out) :: packed(:)

   logical :: ok
   integer :: u, i, n, stat

   ok = .false.
   open(newunit=u, file=path, status="old", action="read", iostat=stat)
   if (stat /= 0) then
      return
   end if
   n = 0
   do
      read(u, *, iostat=stat)
      if (stat /= 0) exit
      n = n + 1
   end do
   rewind(u)
   allocate(packed(n))
   do i = 1, n
      read(u, *, iostat=stat) packed(i)
      if (stat /= 0) then
         close(u)
         return
      end if
   end do
   close(u)
   ok = .true.
end function read_ref


!> Compare packed Hessian against stored reference.
!> If GEN_REFS=1, writes reference from current output.
!> On mismatch, reports the first differing packed element.
subroutine compare_or_write_ref(error, label, packed)
   type(error_type), allocatable, intent(out) :: error
   character(len=*), intent(in) :: label
   real(wp), intent(in) :: packed(:)

   character(len=:), allocatable :: path
   real(wp), allocatable :: ref(:)
   integer :: i
   logical :: ok

   path = ref_path(label)
   if (gen_refs()) then
      call write_ref(path, packed)
      return
   end if
   ok = read_ref(path, ref)
   if (.not. ok) then
      call check(error, ok)
      return
   end if
   if (size(ref) /= size(packed)) then
      call check(error, size(ref), size(packed))
      return
   end if
   do i = 1, size(packed)
      call check(error, packed(i), ref(i), thr=thr)
      if (allocated(error)) exit
   end do
end subroutine compare_or_write_ref


subroutine test_lindh_d2_h2o(error)
   type(error_type), allocatable, intent(out) :: error
   call test_mh(error, "h2o", VAR_LINDH_D2, default_modh(), "lindh_d2_h2o")
end subroutine

subroutine test_lindh_d2_mindless01(error)
   type(error_type), allocatable, intent(out) :: error
   call test_mh(error, "mindless01", VAR_LINDH_D2, default_modh(), "lindh_d2_mindless01")
end subroutine

subroutine test_lindh_d2_caffeine(error)
   type(error_type), allocatable, intent(out) :: error
   call test_mh(error, "caffeine", VAR_LINDH_D2, default_modh(), "lindh_d2_caffeine")
end subroutine

subroutine test_lindh_d2_mgh2(error)
   type(error_type), allocatable, intent(out) :: error
   call test_mh(error, "mgh2", VAR_LINDH_D2, default_modh(), "lindh_d2_mgh2")
end subroutine

subroutine test_lindh_h2o(error)
   type(error_type), allocatable, intent(out) :: error
   call test_mh(error, "h2o", VAR_LINDH, default_modh(), "lindh_h2o")
end subroutine

subroutine test_lindh_mindless01(error)
   type(error_type), allocatable, intent(out) :: error
   call test_mh(error, "mindless01", VAR_LINDH, default_modh(), "lindh_mindless01")
end subroutine

subroutine test_lindh_caffeine(error)
   type(error_type), allocatable, intent(out) :: error
   call test_mh(error, "caffeine", VAR_LINDH, default_modh(), "lindh_caffeine")
end subroutine

subroutine test_lindh_mgh2(error)
   type(error_type), allocatable, intent(out) :: error
   call test_mh(error, "mgh2", VAR_LINDH, default_modh(), "lindh_mgh2")
end subroutine

subroutine test_swart_h2o(error)
   type(error_type), allocatable, intent(out) :: error
   call test_mh(error, "h2o", VAR_SWART, default_modh(), "swart_h2o")
end subroutine

subroutine test_swart_mindless01(error)
   type(error_type), allocatable, intent(out) :: error
   call test_mh(error, "mindless01", VAR_SWART, default_modh(), "swart_mindless01")
end subroutine

subroutine test_swart_caffeine(error)
   type(error_type), allocatable, intent(out) :: error
   call test_mh(error, "caffeine", VAR_SWART, default_modh(), "swart_caffeine")
end subroutine

subroutine test_swart_mgh2(error)
   type(error_type), allocatable, intent(out) :: error
   call test_mh(error, "mgh2", VAR_SWART, default_modh(), "swart_mgh2")
end subroutine

subroutine test_lindh_d2_caffeine_oop(error)
   type(error_type), allocatable, intent(out) :: error
   call test_mh(error, "caffeine", VAR_LINDH_D2, oop_modh(), "lindh_d2_caffeine_oop")
end subroutine

subroutine test_lindh_caffeine_oop(error)
   type(error_type), allocatable, intent(out) :: error
   call test_mh(error, "caffeine", VAR_LINDH, oop_modh(), "lindh_caffeine_oop")
end subroutine

subroutine test_swart_caffeine_oop(error)
   type(error_type), allocatable, intent(out) :: error
   call test_mh(error, "caffeine", VAR_SWART, oop_modh(), "swart_caffeine_oop")
end subroutine

!> Modified Swart (O1NumHess variant) on H2O against hardcoded reference.
subroutine test_modified_swart_h2o(error)
   type(error_type), allocatable, intent(out) :: error
   integer, parameter :: nat = 3
   integer, parameter :: at(nat) = [8, 1, 1]
   real(wp), parameter :: xyz(3, nat) = reshape([&
      & 0.00000000000000_wp,    0.00000000034546_wp,    0.18900383618455_wp, &
      & 0.00000000000000_wp,    1.45674735348811_wp,   -0.88650486059828_wp, &
      &-0.00000000000000_wp,   -1.45674735383357_wp,   -0.88650486086986_wp],&
      & shape(xyz))
   real(wp), parameter :: h0(9, 9) = reshape([&
      & 1.10923040379630813E-003_wp,  0.00000000000000000E+000_wp,  0.00000000000000000E+000_wp, &
      &-5.54615203637871161E-004_wp,  0.00000000000000000E+000_wp,  0.00000000000000000E+000_wp, &
      &-5.54615200158437081E-004_wp,  0.00000000000000000E+000_wp,  0.00000000000000000E+000_wp, &
      & 0.00000000000000000E+000_wp,  4.41847483616303494E-001_wp, -2.75527608981326075E-010_wp, &
      & 0.00000000000000000E+000_wp, -2.20923741962403913E-001_wp,  1.63107694353669314E-001_wp, &
      & 0.00000000000000000E+000_wp, -2.20923741653899580E-001_wp, -1.63107694078141713E-001_wp, &
      & 0.00000000000000000E+000_wp, -2.75527608981326075E-010_wp,  2.94170327797027287E-001_wp, &
      & 0.00000000000000000E+000_wp,  1.43195816033720724E-001_wp, -1.47085164039258476E-001_wp, &
      & 0.00000000000000000E+000_wp, -1.43195815758193096E-001_wp, -1.47085163757768811E-001_wp, &
      &-5.54615203637870945E-004_wp,  0.00000000000000000E+000_wp,  0.00000000000000000E+000_wp, &
      & 1.67247069737169415E-003_wp,  0.00000000000000000E+000_wp,  0.00000000000000000E+000_wp, &
      &-1.11785549373382320E-003_wp,  0.00000000000000000E+000_wp,  0.00000000000000000E+000_wp, &
      & 0.00000000000000000E+000_wp, -2.20923741962403913E-001_wp,  1.43195816033720724E-001_wp, &
      & 0.00000000000000000E+000_wp,  2.33232059827628624E-001_wp, -1.53151755195756606E-001_wp, &
      & 0.00000000000000000E+000_wp, -1.23083178652246897E-002_wp,  9.95593916203589037E-003_wp, &
      & 0.00000000000000000E+000_wp,  1.63107694353669314E-001_wp, -1.47085164039258476E-001_wp, &
      & 0.00000000000000000E+000_wp, -1.53151755195756606E-001_wp,  1.33753677076693056E-001_wp, &
      & 0.00000000000000000E+000_wp, -9.95593915791269726E-003_wp,  1.33314869625654146E-002_wp, &
      &-5.54615200158437190E-004_wp,  0.00000000000000000E+000_wp,  0.00000000000000000E+000_wp, &
      &-1.11785549373382299E-003_wp,  0.00000000000000000E+000_wp,  0.00000000000000000E+000_wp, &
      & 1.67247069389226039E-003_wp,  0.00000000000000000E+000_wp,  0.00000000000000000E+000_wp, &
      & 0.00000000000000000E+000_wp, -2.20923741653899580E-001_wp, -1.43195815758193096E-001_wp, &
      & 0.00000000000000000E+000_wp, -1.23083178652246897E-002_wp, -9.95593915791269726E-003_wp, &
      & 0.00000000000000000E+000_wp,  2.33232059519124291E-001_wp,  1.53151754916105803E-001_wp, &
      & 0.00000000000000000E+000_wp, -1.63107694078141713E-001_wp, -1.47085163757768811E-001_wp, &
      & 0.00000000000000000E+000_wp,  9.95593916203589037E-003_wp,  1.33314869625654146E-002_wp, &
      & 0.00000000000000000E+000_wp,  1.53151754916105803E-001_wp,  1.33753676795203363E-001_wp],&
      & shape(h0))

   type(TEnvironment) :: env
   integer :: i, j
   real(wp), allocatable :: hess_out(:, :)

   call init(env)
   allocate(hess_out(3*nat, 3*nat))
   hess_out = 0.0_wp
   call swart(env, xyz, at, hess_out)

   do i = 1, 3*nat
      do j = 1, 3*nat
         call check(error, hess_out(j, i), h0(j, i), thr=thr)
         if (allocated(error)) return
      end do
   end do
end subroutine test_modified_swart_h2o

end module test_model_hessian
