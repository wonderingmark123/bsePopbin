!----------------------------------------------------------------------
!     Followed the instrucion from
!     "Giacobbo N., Mapelli M., Spera M., 2018, MNRAS, 474, 2959"
!     section 2.1, the MOBSE2 version.
!----------------------------------------------------------------------
***
      real*8 FUNCTION mlwind(kw,lum,r,mt,mc,rl,z)
      implicit none
      integer kw
      real*8 lum,r,mt,mc,rl,z
      real*8 dml,dms,dmt,p0,x,mew,lum0,kap,neta,bwind,hewind,mxns
      parameter(lum0=7.0d+04,kap=-0.5d0)
      common /value1/ neta,bwind,hewind,mxns
!----------------------------------------------------------------------
!     ADDED
!
      real*8 Teff, Teff_jump, alfa, logMdot, w1, w2, vinf_div_vesc, dT
!
!----------------------------------------------------------------------
*
* Calculate stellar wind mass loss.
*
* Apply mass loss of Nieuwenhuijzen & de Jager, A&A, 1990, 231, 134,
* for massive stars over the entire HRD.
      dms = 0.d0
      if(lum.gt.4000.d0)then
         x = MIN(1.d0,(lum-4000.d0)/500.d0)
         dms = 9.6d-15*x*(r**0.81d0)*(lum**1.24d0)*(mt**0.16d0)
         dms = dms*(z/0.02d0)**(1.d0/2.d0)
      endif
      if(kw.ge.2.and.kw.le.9)then
* 'Reimers' mass loss
         dml = neta*4.0d-13*r*lum/mt
         if(rl.gt.0.d0) dml = dml*(1.d0 + bwind*(MIN(0.5d0,(r/rl)))**6)
* Apply mass loss of Vassiliadis & Wood, ApJ, 1993, 413, 641,
* for high pulsation periods on AGB.
         if(kw.eq.5.or.kw.eq.6)then
            p0 = -2.07d0 - 0.9d0*log10(mt) + 1.94d0*log10(r)
            p0 = 10.d0**p0
            p0 = MIN(p0,2000.d0)
            dmt = -11.4d0+0.0125d0*(p0-100.d0*MAX(mt-2.5d0,0.d0))
            dmt = 10.d0**dmt
            dmt = 1.d0*MIN(dmt,1.36d-09*lum)
            dml = MAX(dml,dmt)
         endif
         if(kw.gt.6)then
!----------------------------------------------------------------------
!     EDITED according to:
!     Belczynski K., Bulik T., Fryer C. L., Ruiter A., Valsecchi F.,
!     Vink J. S., Hurley J. R., 2010, ApJ, 714, 1217
!
            dms = 1.0d-13*lum**(3.d0/2.d0)*(z/0.02d0)**0.86d0
            dms = MAX(dml,dms)
!
!     ORIGINAL:
!           dms = MAX(dml,1.0d-13*hewind*lum**(3.d0/2.d0))
!----------------------------------------------------------------------
         else
            dms = MAX(dml,dms)
            mew = ((mt-mc)/mt)*MIN(5.d0,MAX(1.2d0,(lum/lum0)**kap))
* reduced WR-like mass loss for small H-envelope mass
            if(mew.lt.1.d0)then
               dml = 1.0d-13*lum**(3.d0/2.d0)*(1.d0 - mew)
               dms = MAX(dml,dms)
            end if
* LBV-like mass loss beyond the Humphreys-Davidson limit.
            x = 1.0d-5*r*sqrt(lum)
            if(lum.gt.6.0d+05.and.x.gt.1.d0)then
!----------------------------------------------------------------------
!     EDITED according to:
!     Belczynski K., Bulik T., Fryer C. L., Ruiter A., Valsecchi F.,
!     Vink J. S., Hurley J. R., 2010, ApJ, 714, 1217
!
               dms = 1.5d-4
               dms = MAX(dml,dms)
!
!     ORIGINAL:
!              dml = 0.1d0*(x-1.d0)**3*(lum/6.0d+05-1.d0)
!              dms = dms + dml
!----------------------------------------------------------------------
            endif
         endif
      endif
*
      mlwind = dms
!----------------------------------------------------------------------
!     ADDED according to:
!     Vink J. S., de Koter A., Lamers H. J. G. L. M., 2001, A & A,369, 574
!     For hot massive H-rich stars (B/O spectral type)
!
      Teff = 1000.d0*((1130.d0*lum/(r**2.d0))**(1.d0/4.d0))
      if ((Teff.gt.5d4).or.(Teff.lt.1.25d4).or.(kw.gt.6)) then
         return
      endif
!
!     The following lines were copied from Mesa
!     (mesa/star/private/wind.f90: eval_Vink_wind)
!
!     alfa = 1 for hot side, = 0 for cool side
      if (Teff.gt.27500d0) then
         alfa = 1
      else if (Teff.lt.22500d0) then
         alfa = 0
      else
!     use Vink et al 2001, eqns 14 and 15 to set "jump" temperature
         Teff_jump = 61.2d3+2.59d3*(-13.636d0 + 0.889d0*log10(z/0.02d0))
         dT = 100d0
         if (Teff .gt. (Teff_jump + dT)) then
            alfa = 1
         else if (Teff .lt. (Teff_jump - dT)) then
            alfa = 0
         else
            alfa = (Teff - (Teff_jump - dT)) / (2*dT)
         end if
      end if
      if (alfa .gt. 0) then     ! eval hot side wind (eqn 24)
         vinf_div_vesc = 2.6d0  ! this is the hot side galactic value
         vinf_div_vesc = vinf_div_vesc*(z/0.02d0)**0.13d0 ! corrected for Z
         logMdot =
     &        - 6.697d0
     &        + 2.194d0*log10(lum/1d5)
     &        - 1.313d0*log10(mt/30)
     &        - 1.226d0*log10(vinf_div_vesc/2d0)
     &        + 0.933d0*log10(Teff/4d4)
     &        - 10.92d0*log10(Teff/4d4)**2
     &        + 0.85d0*log10(z/0.02d0)
         w1 = 10.d0**logMdot
      else
         w1 = 0
      end if
      if (alfa .lt. 1) then     ! eval cool side wind (eqn 25)
         vinf_div_vesc = 1.3d0  ! this is the cool side galactic value
         vinf_div_vesc = vinf_div_vesc*(z/0.02d0)**0.13d0 ! corrected for Z
         logMdot =
     &        - 6.688d0
     &        + 2.210d0*log10(lum/1d5)
     &        - 1.339d0*log10(mt/30)
     &        - 1.601d0*log10(vinf_div_vesc/2d0)
     &        + 1.07d0*log10(Teff/2d4)
     &        + 0.85d0*log10(z/0.02d0)
         w2 = 10d0**logMdot
      else
         w2 = 0
      end if
      dms = alfa * w1 + (1 - alfa) * w2
      mlwind = MAX(dms,mlwind)
!
!----------------------------------------------------------------------
*
      return
      end
***
