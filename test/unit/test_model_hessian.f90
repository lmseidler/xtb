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

!> Behavior-lock tests for model Hessian variants
!>
!> Reference data lives in `test/unit/fixtures/model_hessian/*.dat` as one
!> packed Hessian value per line. With GEN_REFS = 1, the tests write files from
!> current code output and PASS. Otherwise they compare element-wise against
!> stored references
!>
!> Tests may encode current buggy behavior such as rkl2=sum(rjk**2) in torsion
!> and otherwise lock the current implementation until deliberate changes
!> regenerate the references
module test_model_hessian
   use testdrive, only : new_unittest, unittest_type, error_type, check
   use xtb_mctc_accuracy, only : wp
   use xtb_chargemodel, only : new_charge_model_2019
   use xtb_gfnff_calculator, only : TGFFCalculator, newGFFCalculator
   use xtb_modelhessian_gff, only : mh_gff
   use xtb_type_molecule, only : TMolecule
   use xtb_modelhessian_eeq, only : add_eeq_hessian
   use xtb_modelhessian_type, only : TModelHessian
   use xtb_modelhessian_lindh, only : TLindhModelHessian, TLindhD2ModelHessian
   use xtb_modelhessian_swart, only : TSwartModelHessian
   use xtb_type_param, only : chrg_parameter
   use xtb_type_setvar, only : modhess_setvar
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
      new_unittest("model_hessian_dense", test_model_hessian_dense), &
      new_unittest("model_hessian_charge", test_model_hessian_charge), &
      new_unittest("eeq_addition", test_eeq_addition), &
      new_unittest("gff_h2o", test_gff_h2o), &
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

!> Generic model Hessian test driver
!> Computes Hessian, compares against reference
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
   class(TModelHessian), allocatable :: model_hessian

   n3 = 3 * mol%n
   allocate(hess_packed(n3*(n3+1)/2))
   call new_model_hessian(variant, model_hessian)
   call model_hessian%compute(mol%xyz, mol%n, hess_packed, mol%at, modh)
end subroutine compute_mh_packed


!> Allocate model Hessian implementation for a test variant
subroutine new_model_hessian(variant, model_hessian)
   integer, intent(in) :: variant
   class(TModelHessian), allocatable, intent(out) :: model_hessian

   select case(variant)
   case(VAR_LINDH_D2)
      allocate(TLindhD2ModelHessian :: model_hessian)
   case(VAR_LINDH)
      allocate(TLindhModelHessian :: model_hessian)
   case(VAR_SWART)
      allocate(TSwartModelHessian :: model_hessian)
   end select
end subroutine new_model_hessian


!> Check packed and dense generic model Hessian interfaces are equivalent
subroutine test_model_hessian_dense(error)
   type(error_type), allocatable, intent(out) :: error

   type(TMolecule) :: mol
   class(TModelHessian), allocatable :: model_hessian
   integer :: i, j, ij, n3, variant
   real(wp), allocatable :: hess_packed(:), hess_dense(:, :)

   call getMolecule(mol, "h2o")
   n3 = 3 * mol%n
   allocate(hess_packed(n3*(n3+1)/2), hess_dense(n3, n3))

   do variant = VAR_LINDH_D2, VAR_SWART
      call new_model_hessian(variant, model_hessian)
      call model_hessian%compute(mol%xyz, mol%n, hess_packed, mol%at, default_modh())
      call model_hessian%compute(mol%xyz, mol%n, hess_dense, mol%at, default_modh())

      ij = 0
      do i = 1, n3
         do j = 1, i
            ij = ij + 1
            call check(error, hess_dense(j, i), hess_packed(ij), thr=thr)
            if (allocated(error)) return
            call check(error, hess_dense(i, j), hess_packed(ij), thr=thr)
            if (allocated(error)) return
         end do
      end do
   end do
end subroutine test_model_hessian_dense


!> Check model implementations dispatch the additive charge contribution
subroutine test_model_hessian_charge(error)
   type(error_type), allocatable, intent(out) :: error

   type(TMolecule) :: mol
   type(chrg_parameter) :: chrgeq
   type(modhess_setvar) :: modh
   class(TModelHessian), allocatable :: model_hessian
   integer :: i, n3, variant
   real(wp), allocatable :: base(:), charged(:), contribution(:)

   call getMolecule(mol, "h2o")
   call new_charge_model_2019(chrgeq, mol%n, mol%at)
   n3 = 3 * mol%n
   allocate(base(n3*(n3+1)/2), charged(n3*(n3+1)/2), &
      & contribution(n3*(n3+1)/2))
   contribution = 0.0_wp
   call add_eeq_hessian(mol%n, mol%at, mol%xyz, 0.0_wp, chrgeq, 0.1_wp, contribution)

   do variant = VAR_LINDH_D2, VAR_SWART
      call new_model_hessian(variant, model_hessian)
      modh = default_modh()
      call model_hessian%compute(mol%xyz, mol%n, base, mol%at, modh)
      modh%kq = 0.1_wp
      call model_hessian%compute(mol%xyz, mol%n, charged, mol%at, modh)

      do i = 1, size(base)
         call check(error, charged(i), base(i) + contribution(i), &
            & thr=100*epsilon(0.0_wp))
         if (allocated(error)) return
      end do
   end do
end subroutine test_model_hessian_charge


!> Check EEQ utility adds to, rather than replaces, packed Hessian values
subroutine test_eeq_addition(error)
   type(error_type), allocatable, intent(out) :: error

   type(TMolecule) :: mol
   type(chrg_parameter) :: chrgeq
   integer :: i, n3
   real(wp), allocatable :: contribution(:), shifted(:)

   call getMolecule(mol, "h2o")
   call new_charge_model_2019(chrgeq, mol%n, mol%at)
   n3 = 3 * mol%n
   allocate(contribution(n3*(n3+1)/2), shifted(n3*(n3+1)/2))
   contribution = 0.0_wp
   shifted = 1.0_wp

   call add_eeq_hessian(mol%n, mol%at, mol%xyz, 0.0_wp, chrgeq, 0.1_wp, contribution)
   call add_eeq_hessian(mol%n, mol%at, mol%xyz, 0.0_wp, chrgeq, 0.1_wp, shifted)

   call check(error, any(abs(contribution) > epsilon(0.0_wp)))
   if (allocated(error)) return
   do i = 1, size(contribution)
      call check(error, shifted(i), 1.0_wp + contribution(i), thr=100*epsilon(0.0_wp))
      if (allocated(error)) return
   end do
end subroutine test_eeq_addition


!> Check the GFN-FF model Hessian against values from the current implementation
subroutine test_gff_h2o(error)
   type(error_type), allocatable, intent(out) :: error

   type(TMolecule) :: mol
   type(TEnvironment) :: env
   type(TGFFCalculator) :: calc
   type(modhess_setvar) :: modh
   real(wp), allocatable :: hessian(:)
   real(wp), parameter :: hessian_ref(45) = reshape([&
      & 8.773653166565035E-01_wp,  0.000000000000000E+00_wp,  4.538894511017244E-03_wp, &
      & 0.000000000000000E+00_wp,  0.000000000000000E+00_wp,  1.090816395965398E+00_wp, &
      &-4.386826583282517E-01_wp,  0.000000000000000E+00_wp, -1.225221024048295E-01_wp, &
      & 5.234760825580048E-01_wp,  0.000000000000000E+00_wp, -2.269447255508622E-03_wp, &
      & 0.000000000000000E+00_wp,  0.000000000000000E+00_wp,  1.660799248573899E-03_wp, &
      &-3.395567216330035E-01_wp,  0.000000000000000E+00_wp, -5.454081979826992E-01_wp, &
      & 2.310394120189165E-01_wp,  0.000000000000000E+00_wp,  4.053283378675142E-01_wp, &
      &-4.386826583282517E-01_wp,  0.000000000000000E+00_wp,  1.225221024048295E-01_wp, &
      &-8.479342422975308E-02_wp,  0.000000000000000E+00_wp,  1.085173096140870E-01_wp, &
      & 5.234760825580048E-01_wp,  0.000000000000000E+00_wp, -2.269447255508622E-03_wp, &
      & 0.000000000000000E+00_wp,  0.000000000000000E+00_wp,  6.086480069347232E-04_wp, &
      & 0.000000000000000E+00_wp,  0.000000000000000E+00_wp,  1.660799248573899E-03_wp, &
      & 3.395567216330035E-01_wp,  0.000000000000000E+00_wp, -5.454081979826992E-01_wp, &
      &-1.085173096140870E-01_wp,  0.000000000000000E+00_wp,  1.400798601151851E-01_wp, &
      &-2.310394120189165E-01_wp,  0.000000000000000E+00_wp,  4.053283378675142E-01_wp],&
      & shape(hessian_ref))
   integer :: i, n3
   logical :: terminate

   call init(env)
   call getMolecule(mol, "h2o")
   call newGFFCalculator(env, mol, calc, '.param_gfnff.xtb', .false.)
   call env%check(terminate)
   call check(error, terminate .eqv. .false.)
   if (allocated(error)) return

   n3 = 3 * mol%n
   allocate(hessian(n3*(n3+1)/2))
   modh = default_modh()
   call mh_gff(mol%xyz, mol%n, hessian, mol%at, modh%s6, &
      & calc%param, calc%topo, calc%neigh)

   do i = 1, size(hessian_ref)
      call check(error, hessian(i), hessian_ref(i), thr=100*epsilon(0.0_wp))
      if (allocated(error)) return
   end do
end subroutine test_gff_h2o


!> Path to reference file for a given label
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


!> True when GEN_REFS is set to 1 (force regeneration of files)
function gen_refs() result(mode)
   logical :: mode
   character(len=8) :: buf
   integer :: stat
   call get_environment_variable("GEN_REFS", buf, status=stat)
   mode = (stat == 0 .and. trim(buf) == "1")
end function gen_refs


!> Write packed Hessian as one value per line to a file
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


!> Read reference packed Hessian from file, returning .false. on read error
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


!> Compare packed Hessian against stored reference
!> If GEN_REFS = 1, writes reference from current output
!> On mismatch, reports the first differing packed element
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

end module test_model_hessian
