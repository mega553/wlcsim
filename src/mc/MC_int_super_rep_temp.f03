!---------------------------------------------------------------!
!
!     This subroutine calculates the change in the self energy for
!     a super reptation move.  That is a reptation move where the 
!     chain identities change along with position so that middle
!     beads appear not to change.
!
!     Created from MC_int_rep by Quinn on 8/10/17
!---------------------------------------------------------------
subroutine MC_int_super_rep_temp(wlc_p,wlc_d,IT1,IT2,forward)
use params,only: dp,wlcsim_params,wlcsim_data
implicit none

!   iputs
TYPE(wlcsim_params), intent(in) :: wlc_p   ! <---- Contains output
TYPE(wlcsim_data), intent(inout) :: wlc_d
integer, intent(in) :: IT1  ! Test bead position 1
integer, intent(in) :: IT2  ! Test bead position 2

!   Internal variables
integer I                 ! For looping over bins
integer II                ! For looping over IB
integer IB1,IB2               ! Bead index
integer rrdr ! -1 if r, 1 if r + dr
integer IX(2),IY(2),IZ(2)
real(dp) WX(2),WY(2),WZ(2)
real(dp) WTOT       ! total weight ascribed to bin
real(dp) RBin(3)    ! bead position
integer inDBin              ! index of bin
integer ISX,ISY,ISZ
LOGICAL isA   ! The bead is of type A
real(dp), dimension(-2:2) ::  phi2
integer m_index  ! m from Ylm spherical harmonics
integer NBinX(3)
real(dp) temp    !for speeding up code
LOGICAL, intent(in) :: forward ! move forward
integer AminusB
NBinX = wlc_p%NBinX

wlc_d%NPHI = 0

! unroll the loop

!(case #2: II = 2)
do IB2 = (IT2-wlc_p%nBPM+1),IT2
      if (forward) then
          ! moving forward I2 is added
          rrdr = 1
      else
          ! moving forward I1 is removed
          rrdr = -1
      endif
  !else
     ! print*, "Error in MC_int_rep, II = {1,2}"
     ! stop 1
  !endif
   ! subract current and add new
   if (rrdr == -1) then
       RBin(1) = wlc_d%R(1,IB2)
       RBin(2) = wlc_d%R(2,IB2)
       RBin(3) = wlc_d%R(3,IB2)
       isA = wlc_d%AB(IB2).eq.1
   else
       RBin(1) = wlc_d%RP(1,IB2)
       RBin(2) = wlc_d%RP(2,IB2)
       RBin(3) = wlc_d%RP(3,IB2)
       isA = wlc_d%ABP(IB2).eq.1
   endif
   if (wlc_p%chi_l2_on .and. isA) then
       if (rrdr == -1) then
           call Y2calc(wlc_d%U(:,IB2),phi2)
       else
           call Y2calc(wlc_d%UP(:,IB2),phi2)
       endif
   else
       ! You could give some MS parameter to B as well if you wanted
       phi2=0.0_dp
   endif
   ! --------------------------------------------------
   !
   !  Interpolate beads into bins
   !
   ! --------------------------------------------------
   call interp(wlc_p,RBin,IX,IY,IZ,WX,WY,WZ)

   ! -------------------------------------------------------
   !
   ! Count beads in bins
   !
   ! ------------------------------------------------------
   !   Add or Subtract volume fraction with weighting from each bin
   !   I know that it looks bad to have this section of code twice but it
   !   makes it faster.
   if (isA) then
       do ISX = 1,2
          if ((IX(ISX).le.0).OR.(IX(ISX).ge.(NBinX(1) + 1))) CYCLE
          do ISY = 1,2
             if ((IY(ISY).le.0).OR.(IY(ISY).ge.(NBinX(2) + 1))) CYCLE
             do ISZ = 1,2
                if ((IZ(ISZ).le.0).OR.(IZ(ISZ).ge.(NBinX(3) + 1))) cycle
                WTOT = WX(ISX)*WY(ISY)*WZ(ISZ)
                inDBin = IX(ISX) + (IY(ISY)-1)*NBinX(1) + (IZ(ISZ)-1)*NBinX(1)*NBinX(2)
                ! Generate list of which phi's change and by how much
                I = wlc_d%NPHI
                do
                   if (I.eq.0) then
                      wlc_d%NPHI = wlc_d%NPHI + 1
                      wlc_d%inDPHI(wlc_d%NPHI) = inDBin
                      temp = rrdr*WTOT*wlc_p%beadVolume/wlc_d%Vol(inDBin)
                      wlc_d%DPHIA(wlc_d%NPHI) = temp
                      wlc_d%DPHIB(wlc_d%NPHI) = 0.0_dp
                      if(wlc_p%chi_l2_on) then
                          do m_index = -2,2
                              wlc_d%DPHI_l2(m_index,wlc_d%NPHI) = &
                                  + phi2(m_index)*temp
                          enddo
                      endif
                      exit
                   elseif (inDBin == wlc_d%inDPHI(I)) then
                      temp = rrdr*WTOT*wlc_p%beadVolume/wlc_d%Vol(inDBin)
                      wlc_d%DPHIA(I) = wlc_d%DPHIA(I) + temp
                      if(wlc_p%chi_l2_on) then
                          do m_index = -2,2
                              wlc_d%DPHI_l2(m_index,I) = wlc_d%DPHI_l2(m_index,I) &
                                  + phi2(m_index)*temp
                          enddo
                      endif
                      exit
                   else
                      I = I-1
                   endif
                enddo
             enddo
          enddo
       enddo
   else ! (isB)
       do ISX = 1,2
          if ((IX(ISX).le.0).OR.(IX(ISX).ge.(NBinX(1) + 1))) CYCLE
          do ISY = 1,2
             if ((IY(ISY).le.0).OR.(IY(ISY).ge.(NBinX(2) + 1))) CYCLE
             do ISZ = 1,2
                if ((IZ(ISZ).le.0).OR.(IZ(ISZ).ge.(NBinX(3) + 1))) cycle
                WTOT = WX(ISX)*WY(ISY)*WZ(ISZ)
                inDBin = IX(ISX) + (IY(ISY)-1)*NBinX(1) + (IZ(ISZ)-1)*NBinX(1)*NBinX(2)
                ! Generate list of which phi's change and by how much
                I = wlc_d%NPHI
                do
                   if (I.eq.0) then
                      wlc_d%NPHI = wlc_d%NPHI + 1
                      wlc_d%inDPHI(wlc_d%NPHI) = inDBin
                      wlc_d%DPHIA(wlc_d%NPHI) = 0.0_dp
                      wlc_d%DPHIB(wlc_d%NPHI) = rrdr*WTOT*wlc_p%beadVolume/wlc_d%Vol(inDBin)
                      if(wlc_p%chi_l2_on) then
                          do m_index = -2,2
                              ! This is somewhat wastefull, could eliminate for speedup by having another NPHI for L=2
                              wlc_d%DPHI_l2(m_index,wlc_d%NPHI) = 0.0_dp
                          enddo
                      endif
                      exit
                   elseif (inDBin == wlc_d%inDPHI(I)) then
                      wlc_d%DPHIB(I) = wlc_d%DPHIB(I) + rrdr*WTOT*wlc_p%beadVolume/wlc_d%Vol(inDBin)
                      exit
                   else
                      I = I-1
                   endif
                enddo
             enddo !ISZ
          enddo !ISY
       enddo !ISX
   endif
enddo ! loop over IB  A.k.a. beads
!---------------------------------------------------------------!

!(case # 1: II = 1)

do IB1 = IT1,(IT1 + wlc_p%nBPM - 1)
 ! if (II.eq.1) then
     ! IB = I1
      if (forward) then
          ! moving forward I1 is removed
          rrdr = -1 
      else
          ! moving backward I2 is added
          rrdr = 1
      endif
!  elseif (II.eq.2) then
 !     IB = I2
  !    if (forward) then
   !       ! moving forward I2 is added
   !       rrdr = 1
   !   else
          ! moving forward I1 is removed
    !      rrdr = -1
    !  endif
  !else
    !  print*, "Error in MC_int_rep, II = {1,2}"
     ! stop 1
  !endif
   ! subract current and add new
   if (rrdr == -1) then
       RBin(1) = wlc_d%R(1,IB1)
       RBin(2) = wlc_d%R(2,IB1)
       RBin(3) = wlc_d%R(3,IB1)
       isA = wlc_d%AB(IB1).eq.1
   else
       RBin(1) = wlc_d%RP(1,IB1)
       RBin(2) = wlc_d%RP(2,IB1)
       RBin(3) = wlc_d%RP(3,IB1)
       isA = wlc_d%ABP(IB1).eq.1
   endif
   if (wlc_p%chi_l2_on .and. isA) then
       if (rrdr == -1) then
           call Y2calc(wlc_d%U(:,IB1),phi2)
       else
           call Y2calc(wlc_d%UP(:,IB1),phi2)
       endif
   else
       ! You could give some MS parameter to B as well if you wanted
       phi2=0.0_dp
   endif
   ! --------------------------------------------------
   !
   !  Interpolate beads into bins
   !
   ! --------------------------------------------------
   call interp(wlc_p,RBin,IX,IY,IZ,WX,WY,WZ)

   ! -------------------------------------------------------
   !
   ! Count beads in bins
   !
   ! ------------------------------------------------------
   !   Add or Subtract volume fraction with weighting from each bin
   !   I know that it looks bad to have this section of code twice but it
   !   makes it faster.
   if (isA) then
       do ISX = 1,2
          if ((IX(ISX).le.0).OR.(IX(ISX).ge.(NBinX(1) + 1))) CYCLE
          do ISY = 1,2
             if ((IY(ISY).le.0).OR.(IY(ISY).ge.(NBinX(2) + 1))) CYCLE
             do ISZ = 1,2
                if ((IZ(ISZ).le.0).OR.(IZ(ISZ).ge.(NBinX(3) + 1))) cycle
                WTOT = WX(ISX)*WY(ISY)*WZ(ISZ)
                inDBin = IX(ISX) + (IY(ISY)-1)*NBinX(1) + (IZ(ISZ)-1)*NBinX(1)*NBinX(2)
                ! Generate list of which phi's change and by how much
                I = wlc_d%NPHI
                do
                   if (I.eq.0) then
                      wlc_d%NPHI = wlc_d%NPHI + 1
                      wlc_d%inDPHI(wlc_d%NPHI) = inDBin
                      temp = rrdr*WTOT*wlc_p%beadVolume/wlc_d%Vol(inDBin)
                      wlc_d%DPHIA(wlc_d%NPHI) = temp
                      wlc_d%DPHIB(wlc_d%NPHI) = 0.0_dp
                      if(wlc_p%chi_l2_on) then
                          do m_index = -2,2
                              wlc_d%DPHI_l2(m_index,wlc_d%NPHI) = &
                                  + phi2(m_index)*temp
                          enddo
                      endif
                      exit
                   elseif (inDBin == wlc_d%inDPHI(I)) then
                      temp = rrdr*WTOT*wlc_p%beadVolume/wlc_d%Vol(inDBin)
                      wlc_d%DPHIA(I) = wlc_d%DPHIA(I) + temp
                      if(wlc_p%chi_l2_on) then
                          do m_index = -2,2
                              wlc_d%DPHI_l2(m_index,I) = wlc_d%DPHI_l2(m_index,I) &
                                  + phi2(m_index)*temp
                          enddo
                      endif
                      exit
                   else
                      I = I-1
                   endif
                enddo
             enddo
          enddo
       enddo
   else ! isB
       do ISX = 1,2
          if ((IX(ISX).le.0).OR.(IX(ISX).ge.(NBinX(1) + 1))) CYCLE
          do ISY = 1,2
             if ((IY(ISY).le.0).OR.(IY(ISY).ge.(NBinX(2) + 1))) CYCLE
             do ISZ = 1,2
                if ((IZ(ISZ).le.0).OR.(IZ(ISZ).ge.(NBinX(3) + 1))) cycle
                WTOT = WX(ISX)*WY(ISY)*WZ(ISZ)
                inDBin = IX(ISX) + (IY(ISY)-1)*NBinX(1) + (IZ(ISZ)-1)*NBinX(1)*NBinX(2)
                ! Generate list of which phi's change and by how much
                I = wlc_d%NPHI
                do
                   if (I.eq.0) then
                      wlc_d%NPHI = wlc_d%NPHI + 1
                      wlc_d%inDPHI(wlc_d%NPHI) = inDBin
                      wlc_d%DPHIA(wlc_d%NPHI) = 0.0_dp
                      wlc_d%DPHIB(wlc_d%NPHI) = rrdr*WTOT*wlc_p%beadVolume/wlc_d%Vol(inDBin)
                      if(wlc_p%chi_l2_on) then
                          do m_index = -2,2
                              ! This is somewhat wastefull, could eliminate for speedup by having another NPHI for L=2
                              wlc_d%DPHI_l2(m_index,wlc_d%NPHI) = 0.0_dp
                          enddo
                      endif
                      exit
                   elseif (inDBin == wlc_d%inDPHI(I)) then
                      wlc_d%DPHIB(I) = wlc_d%DPHIB(I) + rrdr*WTOT*wlc_p%beadVolume/wlc_d%Vol(inDBin)
                      exit
                   else
                      I = I-1
                   endif
                enddo
             enddo !ISZ
          enddo !ISY
       enddo !ISX
   endif
enddo! loop over IB  A.k.a. beads
! ---------------------------------------------------------------------
!



!----------------------------------------------------------
!
!  Intermediate Beads Don't change field
!
!-----------------------------------------------------------

! ... skipping

! ---------------------------------------------------------------------
!
! Calcualte change in energy
!
!---------------------------------------------------------------------
call hamiltonian(wlc_p,wlc_d,.false.)

!RETURN
END subroutine

!---------------------------------------------------------------!
