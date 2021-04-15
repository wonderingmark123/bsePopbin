!----------------------------------------------------------------------
!     Binding energy formulae from Loveridge et al., ApJ 743, 49 (2011)
!
      real*8 FUNCTION calc_logBE(kw,logz,mzams,m,logr)
      implicit none

      integer kw
      real*8 logz,mzams,m,logr
      include "be.h"
      integer iz,ig,im,ir
      real*8 logm,dz,dzi,logRbound,dlogBE,lambda
      
!     Find the nearest logZ in the grid:
      iz = 0
      dz = huge(dz)
      do i=1,nz
         dzi = abs(log10(zs(i))-logz)
         if(dzi.lt.dz)then
            dz = dzi
            iz = i
         endif
      enddo

!     Find the group from the mass and evolutionary state (ig = 1-4):
      logm = log10(m)
      if(m.le.LMHMbs(iz))then   ! group LM (low mass)
         if(kw.lt.4)then        ! group LMR (low-mass RGB)
            ig = 1              ! group LMR1 (low-mass RGB 1)
!     Compute boundary radius between LMR1 an LMR2:
            logRbound = 0.d0
            do i=1,nRGBbc
               logRbound = logRbound + RGBb_coef(iz,i)*logm**(i-1)
            enddo
            if(logR.gt.logRbound) ig = 2 ! group LMR2 (low-mass RGB 2)
         else
            ig = 3              ! group LMA (low-mass AGB)
         endif
      else                      ! group HM (high mass)
         ig = 4
      endif

!     Compute the 10-logarithm of the binding energy:
      calc_logBE = 0.d0
      i = 0
      do im=mrange(iz,ig,1),mrange(iz,ig,2),mrange(iz,ig,3)
         do ir=rrange(iz,ig,1),rrange(iz,ig,2),rrange(iz,ig,3)
            i = i + 1
            calc_logBE = calc_logBE + alphas(iz,ig,i)*
     &           (logm**im)*(logr**ir)
         enddo
      enddo

!     Compute and apply the mass-loss correction factor Lambda:
      lambda = 1.d0 + 0.25d0*((mzams-m)/mzams)**2

      calc_logBE = calc_logBE * lambda

!     BE was originally expressed in erg/solar mass to avoid large
!     numbers, so we need to convert to erg here:
      calc_logBE = calc_logBE + logBE0

      return
      end
!
!----------------------------------------------------------------------
