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

! onecentervdip.f90  ....... onecenternuxc.f90
! Program Description
! ===========================================================================
!       This routine calculates the one-center integrals for the exchange-
! correlation interactions of the extended Hubbard model.
!
!   int [n(ishell)n(jshell) Nuxc(natom)]
!
! ===========================================================================
! Code written by:
! Dani y JOM 
! ===========================================================================
!
! Program 1Declaration
! ===========================================================================
        subroutine onecentervdip (nspec, nspec_max, nsh_max, wfmax_points,   &
     &                            iexc, fraction, nsshxc, rcutoffa_max,      &
     &                            xnocc, dqorb, iderorb, what, signature,    &
     &                            drr_rho)
        use constants
        implicit none


!rcutoffa_max necesitamos que sea rcutoffa_min,


 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: iexc
        integer, intent (in) :: nsh_max
        integer, intent (in) :: nspec
        integer, intent (in) :: nspec_max
        integer, intent (in) :: wfmax_points

        integer, intent (in), dimension (nspec_max) :: iderorb
        integer, intent (in), dimension (nspec_max) :: nsshxc
 
        real*8, intent (in) :: fraction
 
        real*8, intent (in), dimension (nspec_max) :: dqorb
        real*8, intent (in), dimension (nspec_max) :: drr_rho
        real*8, intent (in), dimension (nspec_max) :: rcutoffa_max
        real*8, intent (in), dimension (nsh_max, nspec_max) :: xnocc
 
        character (len=70), intent (in) :: signature

        character (len=70), intent (in), dimension (nspec_max) :: what
 
! Output
 
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer ideriv
        integer in1
        integer irho
        integer issh
        integer jssh
        integer nnrho
        integer nssh
        integer irho1,rho1,irho2,rho2,ind,indtot
 
        real*8 dnuxc
        real*8 dnuxcs
        real*8 dq
        real*8 drho
        real*8 exc
        real*8 dexc
        real*8 factor
        real*8 factor1
        real*8 factor2
        real*8 rcutoff
        real*8 rho
        real*8 rhomin
        real*8 rhomax
        real*8 rh
        real*8 rhp
        real*8 rhpp
        real*8 vxc
 
        real*8, dimension (:,:), allocatable :: answer
        real*8, dimension (:), allocatable :: answer1
        real*8, dimension (:), allocatable :: answer2
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
        real*8, dimension (:), allocatable :: rho1c
        real*8, dimension (:), allocatable :: rhop1c
        real*8, dimension (:), allocatable :: rhopp1c
        real*8, dimension (:), allocatable :: xnocc_in
 
        real*8, external :: psiofr

        integer, parameter :: lmax = 2
        integer, parameter :: loop_lmax = 5 ! loop max = loop_lmax-1
        integer :: l,l1,l2,l3,l4,m1,m2,m3,m4
        real*8 :: g,sum_l
        real*8 :: gauntReal
        real*8, dimension (0:lmax,0:lmax,0:lmax,0:lmax,-lmax:lmax,-lmax:lmax,-lmax:lmax,-lmax:lmax) :: I
        real*8, dimension (loop_lmax,0:lmax,0:lmax,0:lmax,0) :: R
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
        g=0.0d0
        l1=0
        l2=0
        l3=0
        l4=0
        m1=0
        m2=0
        m3=0
        m4=0
        l=0
        ind=0
        indtot=0
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
                             ! write(*,'(9I5,2x,1F10.8)') l, l1,l2,l3,l4,m1,m2,m3,m4,gauntReal(l,l1,l2,l3,l4,m1,m2,m3,m4)
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
         write(*,*) 'hay ',ind, ' diferentes de 0, de ',indtot
         index_max=ind-1
         
         allocate (index_l (index_max)) !nos da los que son diferentes de 0 logical
         allocate (index_l1 (index_max)) !nos da los que son diferentes de 0 logical
         allocate (index_l2 (index_max)) !nos da los que son diferentes de 0 logical
         allocate (index_l3 (index_max)) !nos da los que son diferentes de 0 logical
         allocate (index_l4 (index_max)) !nos da los que son diferentes de 0 logical
         allocate (index_m1 (index_max)) !nos da los que son diferentes de 0 logical
         allocate (index_m2 (index_max)) !nos da los que son diferentes de 0 logical
         allocate (index_m3 (index_max)) !nos da los que son diferentes de 0 logical
         allocate (index_m4 (index_max)) !nos da los que son diferentes de 0 logical
         allocate (answer1 (index_max))
         allocate (answer2 (index_max))
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

! Loop over the different charge types (0, -1, or +1).
        ideriv = 0
 
! Loop over the species
        allocate (rho1c (wfmax_points))
        allocate (rhop1c (wfmax_points))
        allocate (rhopp1c (wfmax_points))
        do in1 = 1, nspec
         nssh = nsshxc(in1)
         write (36,400) in1, nssh
 
! Needed for charge corrections:
         dq = dqorb(in1)
         jssh = iderorb(in1)
 
         drho = drr_rho(in1)
         !rcutoff = rcutoffa_min(in1) !, rcutoffa_max(in1)
         allocate (xnocc_in (nssh))
         xnocc_in(1:nssh) = xnocc(1:nssh,in1)
 
! Obtain the density and respective derivatives needed for evaluating the
! exchange-correlation interactions (LDA or GGA).
         call rho1c_store (in1, nsh_max, nssh, dq, jssh, drho, rcutoff,      &
     &                     xnocc_in, ideriv + 1, wfmax_points, rho1c, rhop1c,&
     &                     rhopp1c)
 
! Fix the endpoints and initialize the increments dz and drho.
         rhomin = 0.0d0
         rhomax = rcutoff  ! es el minimo...
 
         nnrho = nint((rhomax - rhomin)/drho) + 1
 
! Here we loop over rho.
       
         do irho1= 1, nnrho
           rho1 = rhomin + dfloat(irho1 - 1)*drho

           factor1 = 2.0d0*drho/3.0d0
           if (mod(irho1, 2) .eq. 0) factor1 = 4.0d0*drho/3.0d0
           if (irho1 .eq. 1 .or. irho1 .eq. nnrho) factor1 = drho/3.0d0

           do irho2= 1, irho1 !ojo cuando irho2=irho1
             rho2 = rhomin + dfloat(irho2 - 1)*drho
 
             factor2 = 2.0d0*drho/3.0d0
             if (mod(irho2, 2) .eq. 0) factor2 = 4.0d0*drho/3.0d0
             if (irho2 .eq. 1 .or. irho2 .eq. irho1) factor2 = drho/3.0d0

             do ind = 1, index_max
               l=index_l(ind)
               l1=index_l1(ind)
               l2=index_l2(ind)
               l3=index_l3(ind)
               l4=index_l4(ind)
               !write(36,*) 'psiofr(in1,l1,rho1)=',in1,l1,rho1,psiofr(in1,l1,rho2)
               !write(36,*) l,l1,l2,l3,l4,R(l,l1,l2,l3,l4)
               R(l,l1,l2,l3,l4)=R(l,l1,l2,l3,l4)+                                  &
               & + factor1*rho1**(1-l*0)*(psiofr(in1,l1,rho1)*psiofr(in1,l2,rho1))*  &
               &   factor2*rho2**(2+l)*(psiofr(in1,l3,rho2)*psiofr(in1,l4,rho2))
             end do !ind
           end do !rho1
         end do !rho2


         do irho1= 1, nnrho
           rho1 = rhomin + dfloat(irho1 - 1)*drho

           factor1 = 2.0d0*drho/3.0d0
           if (mod(irho1, 2) .eq. 0) factor1 = 4.0d0*drho/3.0d0
           if (irho1 .eq. 1 .or. irho1 .eq. nnrho) factor1 = drho/3.0d0

           do irho2= (irho1+1), nnrho 
             rho2 = rhomin + dfloat(irho2 - 1)*drho
 
             factor2 = 2.0d0*drho/3.0d0
             if (mod(irho2, 2) .eq. 0) factor2 = 4.0d0*drho/3.0d0
             if (irho2 .eq. (irho1+1) .or. nnrho .eq. nnrho) factor2 = drho/3.0d0

             do ind = 1, index_max
               l=index_l(ind)
               l1=index_l1(ind)
               l2=index_l2(ind)
               l3=index_l3(ind)
               l4=index_l4(ind)
               R(l,l1,l2,l3,l4)=R(l,l1,l2,l3,l4)+                                  &
               & + factor1*rho1**(2+l)*(psiofr(in1,l1,rho1)*psiofr(in1,l2,rho1))*  &
               &   factor2*rho2**(1-l)*(psiofr(in1,l3,rho2)*psiofr(in1,l4,rho2))
             end do !ind
           end do !rho1
         end do !rho2
        
         do ind = 0, index_max
           write(36,'(9I5,2x,F10.8,2x,F10.8)') index_l(ind), index_l1(ind),index_l2(ind),index_l3(ind),index_l4(ind),index_m1(ind),index_m2(ind),index_m3(ind),index_m4(ind), &
                         &  gauntReal(index_l(ind), index_l1(ind),index_l2(ind),index_l3(ind),index_l4(ind),index_m1(ind),index_m2(ind),index_m3(ind),index_m4(ind)), &
                         &  R(index_l(ind), index_l1(ind),index_l2(ind),index_l3(ind),index_l4(ind))
                          
         enddo



         do issh = 1, nssh
          write (36,500) answer(issh,1:nssh)
         end do
         deallocate (xnocc_in)
         deallocate (answer)
        end do
 

!I(l1,l2,l3,l4,m1,m2,m3,m4) =I(li,mi) =SUM_L (4pi/(2l+ 1)) * R(l,li) * GR(l,li,mi)
! GR(l,l1,l2,l3,l4,m1,m2,m3,m4)=gauntReal(l,l1,l2,l3,l4,m1,m2,m3,m4)
         do l1 = 0 , lmax
           do l2 = 0 , lmax
             do l3 = 0 , lmax
               do l4 = 0 , lmax
                 do m1 = -l1 , l1
                   do m2 = -l2 , l2
                     do m3 = -l3 , l3
                       do m4 = -l4 , l4
                         do l = 0, 4
                           I(l1,l2,l3,l4,m1,m2,m3,m4)=I(l1,l2,l3,l4,m1,m2,m3,m4) +  &
                           & 4*pi/(2*l+1)*R(l,l1,l2,l3,l4)*gauntReal(l,l1,l2,l3,l4,m1,m2,m3,m4)
                         enddo
                       enddo
                     enddo
                   enddo
                 enddo
               enddo
             enddo
           enddo
         enddo






        write (36,*) '  '
        write (*,*) '  '
        write (*,*) ' Writing output to: coutput/vdip_onecenter.dat '
        write (*,*) '  '
 
        close (unit = 36)

! Deallocate Arrays
! ===========================================================================
        deallocate (rho1c, rhop1c, rhopp1c)
 
! Format Statements
! ===========================================================================
100     format (70('='))
200     format (2x, a45)
300     format (a70)
400     format (2x, i3, 2x, i3)
500     format (8d20.10)
        return
        end
