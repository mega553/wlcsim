!--------------------------------------------------------------*
!
!           Makes Monti Carlo Moves
!
!    Quinn created this from reptation move on 8/9/17
!
!---------------------------------------------------------------

! variables that need to be allocated only on certain branches moved into MD to prevent segfaults
! please move other variables in as you see fit
subroutine MC_superReptation(wlc_p,R,U,RP,UP,AB,ABP,IP,IT1,IT2,IB1,IB2&
                  ,rand_stat,forward)
use mersenne_twister
use params, only: dp,wlcsim_params

implicit none
type(wlcsim_params), intent(in) :: wlc_p
real(dp), intent(in) :: R(3,wlc_p%NT)  ! Bead positions
real(dp), intent(in) :: U(3,wlc_p%NT)  ! Tangent vectors
integer, intent(in) :: AB(wlc_p%NT)  ! Current chemical identity
integer, intent(out) :: ABP(wlc_p%NT)  ! Proposed chemical identity
real(dp), intent(out) :: RP(3,wlc_p%NT)  ! Bead positions
real(dp), intent(out) :: UP(3,wlc_p%NT)  ! Tangent vectors
integer, intent(out) :: IP    ! Test polymer
integer, intent(out) :: IT1   ! Index of test bead 1
integer, intent(out) :: IT2   ! Index of test bead 2
integer, intent(out) :: IB1   ! Index of test bead 1
integer, intent(out) :: IB2   ! Index of test bead 2

integer I  ! Test indices
! Things for random number generator
type(random_stat), intent(inout) :: rand_stat  ! status of random number generator
real urnd(1) ! single random number
integer irnd(1)
real(dp) MAG      ! Magnitude of vector
real(dp) DR(3)    ! Displacement for slide move
real(dp) Uvec(3) ! parallel component of triad
real(dp) pDir(3) ! perp component of triad
real(dp) tDir(3) ! twist component of triad
real(dp) r_relative(3) ! r in new coordinate system
real(dp) u_relative(3) ! u in new coordinate system
logical, intent(out) :: forward
integer beadi
real(dp),  allocatable :: Rtemp(:,:)
real(dp),  allocatable  :: Utemp(:,:)

!TOdo saving RP is not actually needed, even in these cases, but Brad's code assumes that we have RP.
if (wlc_p%ring .OR. wlc_p%interp_bead_lennard_jones) then
    RP = R
    UP = U
endif

allocate(Rtemp(3,wlc_p%nb))
allocate(Utemp(3,wlc_p%nb))


! single bead reptation
call random_index(wlc_p%NP,irnd,rand_stat)
IP=irnd(1)
IT1 = wlc_p%NB*(IP-1) + 1
IT2 = wlc_p%NB*(IP-1) + wlc_p%NB


IB1 = 1
IB2 = wlc_p%NB

Rtemp(1:3,IB1:IB2) = R(1:3,IT1:IT2)
Utemp(1:3,IB1:IB2) = U(1:3,IT1:IT2)

! move forward or backward
call random_number(urnd,rand_stat)
if (urnd(1).lt.0.5_dp) then
    forward = .true.
    ! note, this is not the most efficient way of doing things
    ! it would be more efficient to move all nbpm at a time.
    do beadi =1,wlc_p%nbpm
        dR(1) = Rtemp(1,IB1 + 1)-Rtemp(1,IB1)
        dR(2) = Rtemp(2,IB1 + 1)-Rtemp(2,IB1)
        dR(3) = Rtemp(3,IB1 + 1)-Rtemp(3,IB1)

        Uvec(1) = Utemp(1,IB1)
        Uvec(2) = Utemp(2,IB1)
        Uvec(3) = Utemp(3,IB1)
        ! chose coordinate system
        call random_perp(Uvec,pDir,tDir,rand_stat)
        ! find next r and u in new coordinate system
        u_relative(1) = Uvec(1)*Utemp(1,IB1 + 1) + &
                        Uvec(2)*Utemp(2,IB1 + 1) + &
                        Uvec(3)*Utemp(3,IB1 + 1)
        u_relative(2) = pDir(1)*Utemp(1,IB1 + 1) + &
                        pDir(2)*Utemp(2,IB1 + 1) + &
                        pDir(3)*Utemp(3,IB1 + 1)
        u_relative(3) = tDir(1)*Utemp(1,IB1 + 1) + &
                        tDir(2)*Utemp(2,IB1 + 1) + &
                        tDir(3)*Utemp(3,IB1 + 1)
        r_relative(1) = Uvec(1)*dR(1) + &
                        Uvec(2)*dR(2) + &
                        Uvec(3)*dR(3)
        r_relative(2) = pDir(1)*dR(1) + &
                        pDir(2)*dR(2) + &
                        pDir(3)*dR(3)
        r_relative(3) = tDir(1)*dR(1) + &
                        tDir(2)*dR(2) + &
                        tDir(3)*dR(3)


        ! orient coordinate system with end of chain
        Uvec(1) = Utemp(1,IB2)
        Uvec(2) = Utemp(2,IB2) 
        Uvec(3) = Utemp(3,IB2)
        call random_perp(Uvec,pDir,tDir,rand_stat)
        ! update UP and RP
        UP(1,IT2) = Uvec(1)*u_relative(1) + pDir(1)*u_relative(2) + tDir(1)*u_relative(3)
        UP(2,IT2) = Uvec(2)*u_relative(1) + pDir(2)*u_relative(2) + tDir(2)*u_relative(3)
        UP(3,IT2) = Uvec(3)*u_relative(1) + pDir(3)*u_relative(2) + tDir(3)*u_relative(3)
        mag = sqrt(UP(1,IT2)**2 + UP(2,IT2)**2 + UP(3,IT2)**2)
        UP(1,IT2) = UP(1,IT2)/mag
        UP(2,IT2) = UP(2,IT2)/mag
        UP(3,IT2) = UP(3,IT2)/mag
        RP(1,IT2) = Rtemp(1,IB2) + Uvec(1)*r_relative(1) + pDir(1)*r_relative(2) + tDir(1)*r_relative(3)
        RP(2,IT2) = Rtemp(2,IB2) + Uvec(2)*r_relative(1) + pDir(2)*r_relative(2) + tDir(2)*r_relative(3)
        RP(3,IT2) = Rtemp(3,IB2) + Uvec(3)*r_relative(1) + pDir(3)*r_relative(2) + tDir(3)*r_relative(3)

        RP(1:3,IT1:(IT2-1))=Rtemp(1:3,2:IB2)
        UP(1:3,IT1:(IT2-1))=Utemp(1:3,2:IB2)

        ! pointing Rtemp to RP saves time compared to defining Rtemp as RP every loop. 
        Rtemp(1:3,IB1:IB2) = RP(1:3,IT1:IT2)
        Utemp(1:3,IB1:IB2) = UP(1:3,IT1:IT2)
    enddo
    ABP(IT1:(IT2-wlc_p%nbpm))=AB((IT1+wlc_p%nbpm):IT2)
    ABP((IT2-wlc_p%nbpm+1):IT2)=AB(IT1:(IT1+wlc_p%nbpm-1)) ! put end segment type on other end for detail balance

   ! RperpMag = sqrt(r_relative(2)**2 + r_relative(3)**2)
   ! RparaMag = r_relative(1)
   ! call test_equiv_forward(U,R,UP,RP,wlc_p%NT,IT1,IT2,RparaMag,RperpMag)

else
    forward = .false.
    do beadi = 1,wlc_p%nbpm
    dR(1) = Rtemp(1,IB2)-Rtemp(1,IB2-1)
    dR(2) = Rtemp(2,IB2)-Rtemp(2,IB2-1)
    dR(3) = Rtemp(3,IB2)-Rtemp(3,IB2-1)


    Uvec(1) = Utemp(1,IB2); Uvec(2) = Utemp(2,IB2); Uvec(3) = Utemp(3,IB2)
    ! chose coordinate system
    call random_perp(Uvec,pDir,tDir,rand_stat)
    ! find next r and u in new coordinate system
    u_relative(1) = Uvec(1)*Utemp(1,IB2-1) + &
                  Uvec(2)*Utemp(2,IB2-1) + &
                  Uvec(3)*Utemp(3,IB2-1)
    u_relative(2) = pDir(1)*Utemp(1,IB2-1) + &
                  pDir(2)*Utemp(2,IB2-1) + &
                  pDir(3)*Utemp(3,IB2-1)
    u_relative(3) = tDir(1)*Utemp(1,IB2-1) + &
                  tDir(2)*Utemp(2,IB2-1) + &
                  tDir(3)*Utemp(3,IB2-1)
    r_relative(1) = Uvec(1)*dR(1) + &
                  Uvec(2)*dR(2) + &
                  Uvec(3)*dR(3)
    r_relative(2) = pDir(1)*dR(1) + &
                  pDir(2)*dR(2) + &
                  pDir(3)*dR(3)
    r_relative(3) = tDir(1)*dR(1) + &
                  tDir(2)*dR(2) + &
                  tDir(3)*dR(3)

    ! orient coordinate system with end of chain
    Uvec(1) = Utemp(1,IB1); Uvec(2) = Utemp(2,IB1); Uvec(3) = Utemp(3,IB1)
    call random_perp(Uvec,pDir,tDir,rand_stat)
    ! update UP and RP
    UP(1,IT1) = Uvec(1)*u_relative(1) + pDir(1)*u_relative(2) + tDir(1)*u_relative(3)
    UP(2,IT1) = Uvec(2)*u_relative(1) + pDir(2)*u_relative(2) + tDir(2)*u_relative(3)
    UP(3,IT1) = Uvec(3)*u_relative(1) + pDir(3)*u_relative(2) + tDir(3)*u_relative(3)
    mag = sqrt(UP(1,IT1)**2 + UP(2,IT1)**2 + UP(3,IT1)**2)
    UP(1,IT1) = UP(1,IT1)/mag
    UP(2,IT1) = UP(2,IT1)/mag
    UP(3,IT1) = UP(3,IT1)/mag
    RP(1,IT1) = Rtemp(1,IB1)-Uvec(1)*r_relative(1)-pDir(1)*r_relative(2)-tDir(1)*r_relative(3)
    RP(2,IT1) = Rtemp(2,IB1)-Uvec(2)*r_relative(1)-pDir(2)*r_relative(2)-tDir(2)*r_relative(3)
    RP(3,IT1) = Rtemp(3,IB1)-Uvec(3)*r_relative(1)-pDir(3)*r_relative(2)-tDir(3)*r_relative(3)



    RP(1:3,(IT1+1):IT2) = Rtemp(1:3,1:(IB2-1))
    UP(1:3,(IT1+1):IT2) = Utemp(1:3,1:(IB2-1))

    Rtemp(1:3,IB1:IB2) = RP(1:3,IT1:IT2)
    Utemp(1:3,IB1:IB2) = UP(1:3,IT1:IT2)
    enddo
    !do I = IT1 + 1,IT2
      ! RP(1,I) = R(1,I-1)
      ! RP(2,I) = R(2,I-1)
      ! RP(3,I) = R(3,I-1)
       !UP(1,I) = U(1,I-1)
      ! UP(2,I) = U(2,I-1)
      ! UP(3,I) = U(3,I-1)
      ! ABP(I) = AB(I-1)
   ! enddo
    !ABP(IT1)=AB(IT1) ! extend chemical sequence

    ABP((IT1+wlc_p%nbpm):IT2)=AB(IT1:(IT2-wlc_p%nbpm))
    ABP(IT1:(IT1+wlc_p%nbpm-1))=AB((IT2-wlc_p%nbpm+1):IT2) ! put end segment type on other end for detail balance
   endif
 deallocate(Rtemp)
 deallocate(Utemp)
end subroutine
