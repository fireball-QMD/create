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
! Code written by:
! Dani y JOM 
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine onecentervdip (nspec, nspec_max, fraction, nsshxc,    &
     &                           rcutoffa_max,  what, signature,drr_rho)
        use constants
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: nspec
        integer, intent (in) :: nspec_max

        integer, intent (in), dimension (nspec_max) :: nsshxc
 
        real*8, intent (in) :: fraction
        real*8, intent (in), dimension (nspec_max) :: rcutoffa_max
        real*8, intent (in), dimension (nspec_max) :: drr_rho
        character (len=70), intent (in) :: signature
        character (len=70), intent (in), dimension (nspec_max) :: what
 
! Output
 
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer in1
        integer nnrho
        integer irho1,irho2,ind,indtot
        real*8 rho1,rho2 
        real*8 drho
        real*8 factor1
        real*8 factor2
        real*8 rcutoff
        real*8 rho
        real*8 rhomin
        real*8 rhomax
 
        integer :: index_max
        integer, dimension (:), allocatable :: index_l
        integer, dimension (:), allocatable :: index_l1
        integer, dimension (:), allocatable :: index_l2
        integer, dimension (:), allocatable :: index_l3
        integer, dimension (:), allocatable :: index_l4
        integer, dimension (:), allocatable :: index_m1
        integer, dimension (:), allocatable :: index_m2
        integer, dimension (:), allocatable :: index_m3
        integer, dimension (:), allocatable :: index_m4
        real*8, dimension (:), allocatable :: answer
        
        integer :: index_max_l
        real*8, dimension (:), allocatable :: R
 
        real*8, external :: psiofr

        integer, parameter :: lmax = 2
        integer :: l,l1,l2,l3,l4,m1,m2,m3,m4
        real*8 :: gauntReal
! Procedure
! ===========================================================================
! Open the file to store the onecenter data.
        open (unit = 36, file = 'coutput/vdip_onecenter.dat',                &
     &        status = 'unknown')
 
! Set up the header for the output file.
        write (36,100)
        write (36,*) ' All one center matrix elements '
        write (36,*) ' created by: '
        write (36,200) signature
 
        do in1 = 1, nspec
         write (36,300) what(in1)
        end do
        write (36,100)
 

! ===========================================================================
        ind=-1
        indtot=-1
        write(*,*) 'l,l1,l2,l3,l4,m1,m2,m3,m4,gauntReal'
        do l1 = 0 , lmax
          do l2 = 0 , lmax
            do l3 = 0 , lmax
              do l4 = 0 , lmax
                do m1 = -l1 , l1
                  do m2 = -l2 , l2
                    do m3 = -l3 , l3
                      do m4 = -l4 , l4
                        do l = 0, 4 
                          indtot=indtot+1
                          if (gauntReal(l,l1,l2,l3,l4,m1,m2,m3,m4)  .ne. 0.0d0) then
                            ind=ind+1
                            write(*,'(9I5,2x,F10.8)') l, l1,l2,l3,l4,m1,m2,m3,m4,gauntReal(l,l1,l2,l3,l4,m1,m2,m3,m4)
                          end if 
                        enddo
                      enddo
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
        write(*,*) 'hay ',ind, ' diferentes de ',indtot
        index_max=ind
       
        allocate (index_l (index_max)) !nos da los que son diferentes de 0 logical
        allocate (index_l1 (index_max)) 
        allocate (index_l2 (index_max)) 
        allocate (index_l3 (index_max))
        allocate (index_l4 (index_max)) 
        allocate (index_m1 (index_max))
        allocate (index_m2 (index_max)) 
        allocate (index_m3 (index_max))
        allocate (index_m4 (index_max)) 
        allocate (answer (index_max))
        ind=0
        do l1 = 0 , lmax
          do l2 = 0 , lmax
            do l3 = 0 , lmax
              do l4 = 0 , lmax
                do m1 = -l1 , l1
                  do m2 = -l2 , l2
                    do m3 = -l3 , l3
                      do m4 = -l4 , l4
                        do l = 0, 4 
                          if (gauntReal(l,l1,l2,l3,l4,m1,m2,m3,m4)  .ne. 0.0d0) then
                            index_l(ind)=l
                            index_l1(ind)=l1
                            index_l2(ind)=l2
                            index_l3(ind)=l3
                            index_l4(ind)=l4
                            index_m1(ind)=m1
                            index_m2(ind)=m2
                            index_m3(ind)=m3
                            index_m4(ind)=m4
                            !GR(l,l1,l2,l3,l4,m1,m2,m3,m4)=gauntReal(l,l1,l2,l3,l4,m1,m2,m3,m4)
                            ind=ind+1
                          end if 
                        enddo
                      enddo
                    enddo
                  enddo
                enddo
               enddo
             enddo
          enddo
        enddo

         

       allocate (R (index_max_l)) 
       do in1 = 1, nspec
        write (36,400) in1, nsshxc(in1)
        drho = drr_rho(in1) 
        rhomin = 0.0d0
        rhomax = rcutoffa_max(in1)  
        nnrho = nint((rhomax - rhomin)/drho) + 1

! Here we loop over rho.
          
!         do irho1= 1, nnrho
!           rho1 = rhomin + dfloat(irho1 - 1)*drho
!           write(36,*) irho1,rho1,psiofr(in1,1,rho1)
!         end do


        do irho1= 1, nnrho
          rho1 = rhomin + dfloat(irho1 - 1)*drho

          factor1 = 2.0d0*drho/3.0d0
          if (mod(irho1, 2) .eq. 0) factor1 = 4.0d0*drho/3.0d0
          if (irho1 .eq. 1 .or. irho1 .eq. nnrho) factor1 = drho/3.0d0

          do irho2= 1, nnrho 
            rho2 = rhomin + dfloat(irho2 - 1)*drho

            factor2 = 2.0d0*drho/3.0d0
            if (mod(irho2, 2) .eq. 0) factor2 = 4.0d0*drho/3.0d0
            if (irho2 .eq. 1 .or. irho2 .eq. irho1) factor2 = drho/3.0d0

            if (rho1 .lt. 1.0d-04) rho1 = 1.0d-04

            do ind = 0, index_max
              l=index_l(ind)
              l1=index_l1(ind)
              l2=index_l2(ind)
              l3=index_l3(ind)
              l4=index_l4(ind)
              if (irho2 .le. irho1 ) then
                R(ind)=R(ind)+                                                      &
                & + factor1*rho1**(1-l)*(psiofr(in1,l1,rho1)*psiofr(in1,l2,rho1))*  &
                &   factor2*rho2**(2+l)*(psiofr(in1,l3,rho2)*psiofr(in1,l4,rho2))
              else 
                R(ind)=R(ind)+                                                      &
                & + factor1*rho1**(2+l)*(psiofr(in1,l1,rho1)*psiofr(in1,l2,rho1))*  &
                &   factor2*rho2**(1-l)*(psiofr(in1,l3,rho2)*psiofr(in1,l4,rho2))
              end if
            end do !ind
          end do !rho1
        end do !rho2
       end do !nspec


!I(l1,l2,l3,l4,m1,m2,m3,m4) =I(li,mi) =SUM_L (4pi/(2l+ 1)) * R(l,li) * GR(l,li,mi)
! GR(l,l1,l2,l3,l4,m1,m2,m3,m4)=gauntReal(l,l1,l2,l3,l4,m1,m2,m3,m4)

       do ind = 0, index_max
        l=index_l(ind)
        l1=index_l1(ind)
        l2=index_l2(ind)
        l3=index_l3(ind)
        l4=index_l4(ind)
        m1=index_m1(ind)
        m2=index_m2(ind)
        m3=index_m3(ind)
        m4=index_m4(ind)
        answer(ind)=R(ind)*                                                        &
         & gauntReal(index_l(ind),index_l1(ind),index_l2(ind),index_l3(ind),index_l4(ind),index_m1(ind),index_m2(ind),index_m3(ind),index_m4(ind))
        end do !ind

       do ind = 0, index_max
         write(36,'(8I5,2x,2F10.8)') index_l1(ind),index_l2(ind),index_l3(ind),index_l4(ind), &
                                  & index_m1(ind),index_m2(ind),index_m3(ind),index_m4(ind), &
                                  & answer(ind),gauntReal(index_l(ind),index_l1(ind),index_l2(ind),index_l3(ind),index_l4(ind),index_m1(ind),index_m2(ind),index_m3(ind),index_m4(ind))
       end do 

       write (36,*) '  '
       write (*,*) '  '
       write (*,*) ' Writing output to: coutput/vdip_onecenter.dat '
       write (*,*) '  '

       close (unit = 36)

! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
100     format (70('='))
200     format (2x, a45)
300     format (a70)
400     format (2x, i3, 2x, i3)
500     format (8d20.10)
        return
        end
