! This file is part of xtb.
! SPDX-Identifier: LGPL-3.0-or-later
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

module test_dipro
   use testdrive, only : new_unittest, unittest_type, error_type, check_ => check, test_failed, skip_test
   use xtb_features, only : get_xtb_feature
   implicit none
   private

   public :: collect_dipro

contains

!> Collect all exported unit tests
subroutine collect_dipro(testsuite)
   type(unittest_type),allocatable,intent(out) :: testsuite(:)
      !! array of objects(unit tests)

   testsuite = [ &
         !! function that registers a new unit test -> result unittest_type object
      new_unittest("J_ab,eff", test_dipro_jabeff), &
      new_unittest("J_ab,eff PTB", test_dipro_ptb), &
      new_unittest("J_ab,eff PTB monomer scaled", test_dipro_ptb_paper) &
      ]

end subroutine collect_dipro

!---------------------------------------------
! Unit test for dipro
!---------------------------------------------
subroutine test_dipro_jabeff(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout
   use xtb_mctc_convert, only : aatoau

   use xtb_type_environment
   use xtb_type_molecule
   use xtb_type_restart
   use xtb_type_basisset
   use xtb_type_environment
   use xtb_setparam
   
   use xtb_dipro, only: get_jab, jab_input
   use xtb_tblite_calculator, only : TTBLiteCalculator, TTBLiteInput, newTBLiteWavefunction
   use xtb_io_writer 
   use xtb_mctc_filetypes, only : filetype
   type(error_type), allocatable, intent(out) :: error
   type(TMolecule)  :: mol
   type(TEnvironment)       :: env
   type(jab_input)  :: dipro
   type(TTBLiteInput) :: tblite

   real(wp), parameter :: thr = 0.002_wp
   !> Molecular structure data
   integer, parameter :: nat = 26
   !> fragments to be determined in call get_jab
   integer, allocatable :: splitlist(:)
      !! benzene - nitrobenzene cofacial dimer
   integer, parameter :: at(nat) = [6,6,6,6,6,6,1,1,1,1,1,6,6,6,6,6,6,1,1,1,1,1,1,7,8,8]
   real(wp), parameter :: xyz(3,nat) = reshape(aatoau*[&
      &  -1.20618_wp,  0.69599_wp,  1.75211_wp, &
      &  -1.20379_wp, -0.70017_wp,  1.75215_wp, &
      &   0.00616_wp, -1.39690_wp,  1.75210_wp, &
      &   1.21824_wp, -0.70002_wp,  1.75171_wp, &
      &   1.23033_wp,  0.70989_wp,  1.75117_wp, &
      &   0.00125_wp,  1.40067_wp,  1.75172_wp, &
      &   2.14504_wp, -1.25925_wp,  1.74975_wp, &
      &  -0.02278_wp,  2.48280_wp,  1.74976_wp, &
      &  -2.14520_wp,  1.23405_wp,  1.75041_wp, &
      &  -2.14028_wp, -1.24255_wp,  1.75044_wp, &
      &   0.00554_wp, -2.47919_wp,  1.75040_wp, &
      &  -1.21072_wp,  0.69984_wp, -1.75114_wp, &
      &  -1.21074_wp, -0.69902_wp, -1.75149_wp, &
      &   0.00072_wp, -1.39843_wp, -1.75114_wp, &
      &   1.21220_wp, -0.69897_wp, -1.75037_wp, &
      &   1.21222_wp,  0.69988_wp, -1.74995_wp, &
      &   0.00077_wp,  1.39929_wp, -1.75036_wp, &
      &   2.14941_wp, -1.24008_wp, -1.74839_wp, &
      &   2.14945_wp,  1.24099_wp, -1.74755_wp, &
      &   0.00077_wp,  2.48148_wp, -1.74837_wp, &
      &  -2.14799_wp,  1.24091_wp, -1.74971_wp, &
      &  -2.14800_wp, -1.24014_wp, -1.75028_wp, &
      &   0.00066_wp, -2.48066_wp, -1.74972_wp, &
      &   2.46723_wp,  1.42685_wp,  1.74737_wp, &
      &   2.47280_wp,  2.62117_wp,  1.74468_wp, &
      &   3.50640_wp,  0.83800_wp,  1.74471_wp], shape(xyz))
   
   if (.not.get_xtb_feature('tblite')) then
      call skip_test(error, 'tblite libary not available.')
      return
   end if

   call init(env)
      !! construct calculation environment

   call init(mol, at, xyz)
      !! construct molecular structure type
      !! interface to initMoleculeNumbers
   mol%chrg = 0.0_wp
   tblite%method = 'gfn2'
   dipro%othr = 0.1_wp
   
   call get_jab(env,tblite,mol,splitlist,dipro)
   call check_(error,dipro%totjab(1),0.119_wp, thr=thr)

endsubroutine test_dipro_jabeff

!> Ethen dimer stacked at 3.5 Angstrom, PTB with projected site energies
subroutine test_dipro_ptb(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_dipro, only: jab_input
   type(error_type), allocatable, intent(out) :: error
   type(jab_input)  :: dipro

   call run_ethen_ptb(error, dipro)
   if (allocated(error)) return
   call check_(error,dipro%totjab(1),0.00785_wp, thr=1.0e-3_wp)
   if (allocated(error)) return
   call check_(error,dipro%totjab(2),0.2792_wp, thr=2.8e-3_wp)

end subroutine test_dipro_ptb

!> Ethen dimer with PTB, monomer site energies and scaled J_ab as in the PTB DIPRO paper
subroutine test_dipro_ptb_paper(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_dipro, only: jab_input
   type(error_type), allocatable, intent(out) :: error
   type(jab_input)  :: dipro

   dipro%monomer_energies = .true.
   dipro%scale = 1.921_wp
   call run_ethen_ptb(error, dipro)
   if (allocated(error)) return
   call check_(error,dipro%totjab(1),0.8421_wp, thr=8.4e-3_wp)
   if (allocated(error)) return
   call check_(error,dipro%totjab(2),0.4643_wp, thr=4.6e-3_wp)

end subroutine test_dipro_ptb_paper

!> DIPRO with PTB for the ethen dimer stacked at 3.5 Angstrom
subroutine run_ethen_ptb(error, dipro)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : aatoau
   use xtb_type_environment
   use xtb_type_molecule
   use xtb_dipro, only: get_jab, jab_input
   use xtb_tblite_calculator, only : TTBLiteInput
   type(error_type), allocatable, intent(out) :: error
   type(jab_input), intent(inout) :: dipro
   type(TMolecule)  :: mol
   type(TEnvironment)       :: env
   type(TTBLiteInput) :: tblite
   !> fragments to be determined in call get_jab
   integer, allocatable :: splitlist(:)
   integer, parameter :: nat = 12
   integer, parameter :: at(nat) = [6,6,1,1,1,1,6,6,1,1,1,1]
   real(wp), parameter :: xyz(3,nat) = reshape(aatoau*[&
      &   0.000000_wp,  0.000000_wp,  0.000000_wp, &
      &   1.332000_wp,  0.000000_wp,  0.000000_wp, &
      &  -0.574301_wp,  0.000000_wp, -0.928785_wp, &
      &  -0.574301_wp,  0.000000_wp,  0.928785_wp, &
      &   1.906301_wp,  0.000000_wp,  0.928785_wp, &
      &   1.906301_wp,  0.000000_wp, -0.928785_wp, &
      &  -0.000000_wp,  3.500000_wp,  0.000000_wp, &
      &   1.332000_wp,  3.500000_wp, -0.000000_wp, &
      &  -0.574301_wp,  3.500000_wp,  0.928785_wp, &
      &  -0.574301_wp,  3.500000_wp, -0.928785_wp, &
      &   1.906301_wp,  3.500000_wp, -0.928785_wp, &
      &   1.906301_wp,  3.500000_wp,  0.928785_wp], shape(xyz))
   logical :: exitRun

   if (.not.get_xtb_feature('tblite')) then
      call skip_test(error, 'tblite libary not available.')
      return
   end if

   call init(env)
   call init(mol, at, xyz)
   mol%chrg = 0.0_wp
   tblite%method = 'ptb'
   dipro%othr = 0.1_wp

   call get_jab(env,tblite,mol,splitlist,dipro)
   call env%check(exitRun)
   if (exitRun) call test_failed(error, "DIPRO calculation with PTB failed")

end subroutine run_ethen_ptb

end module test_dipro
