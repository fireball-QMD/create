! copyright info:
!
!                             @Copyright 2004
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega
! Universidad de Madrid - Pavel Jelinek

! Other contributors, past and present:
! Auburn University - Jianjun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
! Brigham Young University - Hao Wang
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Motorola, Physical Sciences Research Labs - Alex Demkov
! Motorola, Physical Sciences Research Labs - Jun Wang
! Ohio University - Dave Drabold
! University of Regensburg - Juergen Fritsch

!
! fireball-qmd is a free (GPLv3) open project.

! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

! onecentervdip.f90  
! Program Description
! ===========================================================================
!
! ===========================================================================
! Code written iy:
! Dani i JOM 
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine onecentervdip (nsh_max, nspec, nspec_max, fraction, nsshxc,    &
     &                   lsshxc, drr_rho, rcutoffa_max, what, signature)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: nsh_max
        integer, intent (in) :: nspec
        integer, intent (in) :: nspec_max

        integer, intent (in), dimension (nspec_max, nsh_max) :: lsshxc
        integer, intent (in), dimension (nspec_max) :: nsshxc

        real*8, intent (in) :: fraction

        real*8, intent (in), dimension (nspec_max) :: drr_rho
        real*8, intent (in), dimension (nspec_max) :: rcutoffa_max


        character (len=70), intent (in), dimension (nspec_max) :: what

        character (len=70), intent (in) :: signature


! Output
 
! Local Parameters and Data Declaration
! ===========================================================================
        real*8, parameter :: eq2 = 14.39975d0
 
! Local Variable Declaration and Description
! ===========================================================================
        integer irho
        integer irhop
        integer issh
        integer itype
        integer jssh
        integer lalpha
        integer lmu
        integer lqn
        integer malpha
        integer mmu
        integer mqn
        integer nnrho

        real*8 coefficient
        real*8 cg1
        real*8 cg2
        real*8 cg3
        real*8 cg4
        real*8 drho
        real*8 psi1
        real*8 psi2
        real*8 psi3
        real*8 psi4
        real*8 r
        real*8 rp
        real*8 rhomax
        real*8 rhomin
        real*8 sumr
        real*8 sumrp
        real*8 aux

        real*8, dimension (:), allocatable :: factor
        real*8, dimension (:), allocatable :: rpoint

        real*8, external :: psiofr

        integer :: l,l1,l2,l3,l4,m1,m2,m3,m4
        real*8 ::gauntReal



! Procedure
! ===========================================================================
! Open the file to store the onecenter data.
        open (unit = 36, file = 'coutput/vdip_onecenter.dat',                &
     &        status = 'unknown')

 
! Set up the header for the output file.
        write (36,100)
        write (36,*) ' All one center matrix elements created by:'
        write (36,200) signature

        do itype = 1, nspec
         write (36,300) what(itype)
        end do
        write (36,100)
 
        do itype = 1, nspec
 
! Initialize the limits of integration for the radial integral.
! Set up the grid points for the rho integration.
         rhomin = 0.0d0
         rhomax = rcutoffa_max(itype)
         drho = drr_rho(itype)
         nnrho = int((rhomax - rhomin)/drho) + 1

 
         allocate (rpoint(nnrho))
         allocate (factor(nnrho))
         do irho = 1, nnrho
          rpoint(irho) = float(irho - 1)*drho
! Set up the Simpson rule factors:
          factor(irho) = 2.0d0*drho/3.0d0
          if (mod(irho,2) .eq. 0) factor(irho) = 4.0d0*drho/3.0d0
          if (irho .eq. 1 .or. irho .eq. nnrho) factor(irho) = drho/3.0d0
         end do !irho
 

         do l1 = 0 , nsshxc(itype)
          do l2 = 0 , nsshxc(itype)
           do l3 = 0 , nsshxc(itype)
            do l4 = 0 , nsshxc(itype)
             do lqn = 0, 4
              aux=0
              do m1 = -l1 , l1
               do m2 = -l2 , l2
                do m3 = -l3 , l3
                 do m4 = -l4 , l4
                  aux=aux+abs(gauntReal(lqn,l1,l2,l3,l4,m1,m2,m3,m4))
                 enddo
                enddo
               enddo 
              enddo
! Perform the radial integration. Only do this integral if the
! coefficient is
! non-zero.
              if (aux .gt. 1.0d-4) then
 
              do irho = 1, nnrho
               r = rpoint(irho)
               if (r .lt. 1.0d-04) r = 1.0d-04
                psi1 = psiofr(itype,l1,r)
                psi2 = psiofr(itype,l2,r)
             !  write(36,*)irho,r,psi1,l1,psi2,l2
             end do


! First integrate the even pieces and then the odd pieces
              sumr = 0.0d0
              do irho = 1, nnrho
               r = rpoint(irho)
               if (r .lt. 1.0d-04) r = 1.0d-04
                psi1 = psiofr(itype,l1,r)
                psi2 = psiofr(itype,l2,r)
 
! ****************************************************************************
! Perform the radial integration over r'.
! Limits from 0 to r.
                sumrp = 0.0d0
                do irhop = 1, nnrho
                 rp = rpoint(irhop)
                 if (rp .lt. 1.0d-04) rp = 1.0d-04
                 psi3 = psiofr(itype,l3,rp)
                 psi4 = psiofr(itype,l4,rp)
                 if (rp .le. r) then
                  sumrp = sumrp + factor(irhop)*psi3*psi4*rp**(lqn + 2)/r**(lqn + 1)   ! Limits from 0 to rcutoff.
                 else
                  sumrp = sumrp + factor(irhop)*psi3*psi4*r**lqn/rp**(lqn - 1)
                 end if
!                if (sumr .gt. 1.0d-04 ) then
!                    write(36,'(10F14.10)')  sumr, psi1,psi2,r,sumrp,psi3,psi3,rp
!                end if
                end do !rhop
! ****************************************************************************
                sumr = sumr + factor(irho)*sumrp*psi1*psi2*r**2  !*coefficient
               end do !irho
 
               write (36,'("itype = ",I2,", l = ",I2,", li = (",4I2,")",2x,", sumr =",F14.10)') itype,lqn,l1,l2,l3,l4,sumr
                             !answer(malpha) = answer(malpha) + (eq2/2.0d0)*fraction*sumr
              end if !(aux .gt. 1.0d-4)
             end do !lqn
            enddo !l4
           enddo !l3
          enddo !l2
         enddo !l1
 
! ****************************************************************************
! End loop over the species.
         deallocate (rpoint)
         deallocate (factor)
        end do

        write (36,*) '  '
        write (*,*) '  '
        write (*,*) ' Writing output to: coutput/vdip_onecenter.dat '
        write (*,*) '  '
        close (unit = 36)

! Format itatements
! ===========================================================================
100     format (70('='))
200     format (2x, a45)
300     format (a70)
400     format (8d20.10)
 
        return
        end
