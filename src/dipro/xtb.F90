! This file is part of the dipro module in xtb.
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

#ifndef WITH_TBLITE
#define WITH_TBLITE 0
#endif

!> Helper routines to handle calculations with the tblite library
module xtb_dipro_xtb
   use mctc_env, only : error_type, fatal_error, wp
   use mctc_io, only : structure_type
#if WITH_TBLITE
   use xtb_type_environment, only : TEnvironment
   use xtb_type_molecule, only : TMolecule, assignment(=)
   use xtb_type_data, only : scc_results
   use xtb_type_restart, only : TRestart
   use xtb_setparam, only : set
   use xtb_ptb_calculator, only : TPTBCalculator, newPTBCalculator
   use xtb_ptb_integrals, only : get_integrals
   use tblite_adjlist, only : adjacency_list, new_adjacency_list
   use tblite_basis_type, only : get_cutoff, basis_type
   use tblite_context_type, only : context_type
   use tblite_cutoff, only : get_lattice_points
   use tblite_integral_native_integrals, only : get_overlap
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_gfn2, only : new_gfn2_calculator
   use tblite_xtb_gfn1, only : new_gfn1_calculator
   use tblite_xtb_ipea1, only : new_ipea1_calculator
   use tblite_xtb_singlepoint, only : xtb_singlepoint
   implicit none
   private

   public :: get_singlepoint

   character(len=*), parameter :: source = 'xtb_dipro_xtb'
   !> Conversion factor from temperature to energy (Boltzmann's constant in atomic units)
   real(wp), parameter :: ktoau = 3.166808578545117e-06_wp

contains

!> Single point for a structure with the requested method, returns basis set, wavefunction and optionally the AO overlap
subroutine get_singlepoint(env, mol, method, accuracy, verbosity, bas, wfn, overlap)
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Name of the requested method
   character(len=*), intent(in) :: method
   !> Calculation accuracy
   real(wp), intent(in) :: accuracy
   !> Output verbosity
   integer, intent(in) :: verbosity
   !> Basis set of the calculation
   type(basis_type), intent(out) :: bas
   !> Converged wavefunction
   type(wavefunction_type), intent(out) :: wfn
   !> Overlap matrix in the AO basis of the calculation
   real(wp), allocatable, intent(out), optional :: overlap(:, :)

   type(xtb_calculator) :: xcalc
   type(TPTBCalculator) :: ptb
   type(TMolecule) :: tmol
   type(TRestart) :: chk
   type(scc_results) :: results
   type(adjacency_list) :: list
   type(context_type) :: ctx
   type(error_type), allocatable :: error
   real(wp), allocatable :: trans(:, :), gradient(:, :)
   real(wp) :: energy, cutoff, hlgap, sigma(3, 3)
   logical :: exitRun

   select case(method)
   case("ptb")
      tmol = mol
      call newPTBCalculator(env, tmol, ptb, accuracy)
      call env%check(exitRun)
      if (exitRun) return
      allocate(gradient(3, mol%nat))
      call ptb%singlepoint(env, tmol, chk, verbosity, .false., energy, gradient, sigma, hlgap, results)
      call env%check(exitRun)
      if (exitRun) return
      bas = ptb%bas
      wfn = chk%tblite

      if (present(overlap)) then
         !> PTB coefficients refer to the normalized vDZP basis
         allocate(overlap(bas%nao, bas%nao))
         cutoff = get_cutoff(bas)
         call get_lattice_points(mol%periodic, mol%lattice, cutoff, trans)
         call new_adjacency_list(list, mol, trans, cutoff)
         call get_integrals(mol, bas, trans, list, overlap)
      end if
   case default
      call get_calculator(xcalc, mol, method, error)
      if (allocated(error)) then
         call env%error(error%message, source)
         return
      end if
      call new_wavefunction(wfn, mol%nat, xcalc%bas%nsh, xcalc%bas%nao, 1, set%etemp * ktoau)
      call xtb_singlepoint(ctx, mol, xcalc, wfn, accuracy, energy, verbosity=verbosity)
      if (ctx%failed()) then
         do while (ctx%failed())
            call ctx%get_error(error)
            call env%error(error%message, source)
         end do
         return
      end if
      bas = xcalc%bas

      if (present(overlap)) then
         allocate(overlap(bas%nao, bas%nao))
         cutoff = get_cutoff(bas)
         call get_lattice_points(mol%periodic, mol%lattice, cutoff, trans)
         call get_overlap(mol, trans, cutoff, bas, overlap)
      end if
   end select
end subroutine get_singlepoint

!> Create a new calculator for a given method string
subroutine get_calculator(xcalc, mol, method, error)
   !> Instance of the new calculator
   type(xtb_calculator), intent(out) :: xcalc
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Name of the requested method
   character(len=*), intent(in) :: method
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   select case(method)
   case default
      call fatal_error(error, "Unknown method '"//method//"' requested")
!      error stop
   case("gfn2")
      call new_gfn2_calculator(xcalc, mol, error)
   case("gfn1")
      call new_gfn1_calculator(xcalc, mol, error)
   case("ipea1")
      call new_ipea1_calculator(xcalc, mol, error)
   end select
end subroutine get_calculator
#endif

end module xtb_dipro_xtb
