!----------------------------------------------------------------------
!     Binding energy formulae from Loveridge et al., ApJ 743, 49 (2011)
!
!     Number of available metallicities
      integer nz
      parameter (nz=6)
!     List of available metallicities
      real*8 zs(nz)
      data zs / 1.d-4, 1.d-3, 0.01d0, 0.015d0, 0.02d0, 0.03d0 /

!     Number of model groups
      integer ngr
      parameter (ngr=5)
!     List of model groups
      character*5 groups(ngr)
      data groups / 'LMR1 ','LMR2 ','LMA  ','HM   ','Recom' /

!     Maximum number of coefficients
      integer ndatmax
      parameter (ndatmax=441)

!     Fitting coefficients (alpha_m,r)
      real*8 alphas(nz,ngr,ndatmax)

!     Range of powers for M for each coefficient
      integer mrange(nz,ngr,3)
!     Range of powers for R for each coefficient
      integer rrange(nz,ngr,3)

!     Low-mass/high-mass boundaries (boundaries between groups LM* and HM)
      real*8 LMHMbs(nz)

!     Number of coefficients for RGB boundaries (boundaries between
!     groups LMR1 and LMR2)
      integer nRGBbc
      parameter (nRGBbc=5)
!     Coefficients for RGB boundaries (boundaries between groups LMR1 and LMR2)
      real*8 RGBb_coef(nz,nRGBbc)

!     The binding energy was originally expressed in erg per solar mass
!     to avoid large numbers, so we need to add this constant
!     (log10[Mo/g]) to log(BE).
      real*8 logBE0
      parameter(logBE0=33.29866d0) ! log10(1.9891e33) = 33.29866

      integer i

      include "bedata/z000010.h"
      include "bedata/z000100.h"
      include "bedata/z001000.h"
      include "bedata/z001500.h"
      include "bedata/z002000.h"
      include "bedata/z003000.h"
!
!----------------------------------------------------------------------
