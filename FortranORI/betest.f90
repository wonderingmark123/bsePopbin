!----------------------------------------------------------------------
!     Test the output of calc_logBE, in reference to `betest_ref.txt`
!     (`example_output.txt` from Loveridge et al. ApJ 743, 49, 2011)
!
program binding_energy
  implicit none
  !> Double-precision kind
  integer, parameter :: double = selected_real_kind(15,307)
  !> Double-precision kind
  integer, parameter :: dbl = double
  
  integer :: gb
  real(double) :: Z,Mzams,M,R,E,Er
  real(double) :: logZ,logR,logE,logEr
  character :: path*(99),gbs(2)*(3)
  integer :: Mi,Ri,Zi,i
  real(double) :: calc_logBE
  real(double) :: Mset

  gbs = (/'RGB','AGB'/)  ! Names of the two giant branches (gb=1,2)

  Mset=0.06

  Z = 0.01
  M = Mset
  Mzams = M  ! Assume no wind mass loss since ZAMS
  R = 25.0
  gb = 1

  logZ = log10(Z)

  write(6,*)
  write(6,'(A)')'1:  Star of 1Mo, Z=0.02, no wind mass loss and different radii (RGB+AGB):'
  write(6,'(A6,4A12, 2A8, 2(5x, 2A14))')'#','Z','Mzams (Mo)','M (Mo)','R (Ro)','gb','branch','log(BE/erg)','BE (erg)', &
       'log(BE_r/erg)','BE_r (erg)'
  i = 0
  do Ri = 30,300,30
     i = i+1
     R = dble(Ri)
     if(R.gt.170.d0) gb = 2  ! Assume AGB for R>170Ro
     logR = log10(R)

     if(gb.eq.1)then
        logE = calc_logBE(3, logZ, Mzams, M, logR)
     else
        logE = calc_logBE(5, logZ, Mzams, M, logR)
     end if
     ! logEr = calc_logBE_recom(logZ, M, logR)
     logEr = 0.d0
     
     E  = 10.d0**logE
     Er = 10.d0**logEr
     
     write(6,'(I6,F12.4,3F12.2, I8,A8, 2(5x, F14.3,ES14.2))')i,Z,Mzams,M,R,gb,gbs(gb),logE,E,logEr,Er
  end do

  
  
  
  Z = 0.02
  M = Mset
  Mzams = M/0.9d0
  R = 25.0
  gb = 1

  write(6,*)
  write(6,'(A)')'2:  Stars of different masses, for Z=0.02, R=25Ro (RGB):'
  write(6,'(A6,4A12, 2A8, 2(5x, 2A14))')'#','Z','Mzams (Mo)','M (Mo)','R (Ro)','gb','branch','log(BE/erg)','BE (erg)', &
       'log(BE_r/erg)','BE_r (erg)'
  i = 0
  do Mi = 1,10
     i = i+1
     M = dble(Mi)
     Mzams = M/0.9d0   ! Assume 10% mass loss since ZAMS
     logZ = log10(Z)
     logR = log10(R)

     if(gb.eq.1)then
        logE = calc_logBE(3, logZ, Mzams, M, logR)
     else
        logE = calc_logBE(5, logZ, Mzams, M, logR)
     end if
     ! logEr = calc_logBE_recom(logZ, M, logR)
     logEr = 0.d0

     E = 10.d0**logE
     Er = 10.d0**logEr

     write(6,'(I6,F12.4,3F12.2, I8,A8, 2(5x, F14.3,ES14.2))')i,Z,Mzams,M,R,gb,gbs(gb),logE,E,logEr,Er
  end do



  Z = 0.02
  M = Mset
  Mzams = M   ! Assume no mass loss since ZAMS
  R = 100.0
  gb = 1

  write(6,*)
  write(6,'(A)')'3:  Stars of 1Mo and different metallicities (code chooses nearest Z), for R=100Ro (RGB):'
  write(6,'(A6,4A12, 2A8, 2(5x, 2A14))')'#','Z','Mzams (Mo)','M (Mo)','R (Ro)','gb','branch','log(BE/erg)','BE (erg)', &
       'log(BE_r/erg)','BE_r (erg)'
  i = 0
  do Zi = 1,300,20
     i = i+1
     Z = dble(Zi)*1.d-4
     logZ = log10(Z)
     logR = log10(R)

     if(gb.eq.1)then
        logE = calc_logBE(3, logZ, Mzams, M, logR)
     else
        logE = calc_logBE(5, logZ, Mzams, M, logR)
     end if
     ! logEr = calc_logBE_recom(logZ, M, logR)
     logEr = 0.d0

     E = 10.d0**logE
     Er = 10.d0**logEr

     write(6,'(I6,F12.4,3F12.2, I8,A8, 2(5x, F14.3,ES14.2))')i,Z,Mzams,M,R,gb,gbs(gb),logE,E,logEr,Er
  end do

  write(6,*)
end program binding_energy
!
!----------------------------------------------------------------------
