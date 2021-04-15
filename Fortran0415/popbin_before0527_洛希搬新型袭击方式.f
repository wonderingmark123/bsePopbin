***
      PROGRAM popbin
***
*
* Evolves a population of binaries using input parameters 
* read from input file binaries.in (M1, M2, P, e, Z, Tmax). 
*
***
      implicit none
*
      INCLUDE 'const_bse.h'
*
      integer i,j,k,jj,nm1
      integer kw,kw2,kwx,kwx2,kstar(2)
      integer i1,i2,kdum
*
      real*8 m1,m2,tmax
      real*8 mass0(2),mass(2),z,zpars(20)
        real*8 mass01,mass02
        common /ini_mass/ mass01,mass02
      real*8 epoch(2),tms(2),tphys,tphysf,dtp
      real*8 rad(2),lum(2),ospin(2)
      real*8 massc(2),radc(2),menv(2),renv(2)
      real*8 sep0,tb0,tb,ecc0,ecc,aursun,yeardy,yearsc,tol
      PARAMETER(aursun=214.95d0,yeardy=365.25d0,yearsc=3.1557d+07)
      PARAMETER(tol=1.d-07)
      real*8 t1,t2,mx,mx2,tbx,eccx

        real*8 dlnm1,dlnm2,dlna,a,a_day,aa,bb,cc,rzams1,rzams2,rl
        real*8 M1max,M2max,M1min,M2min,amax,amin
***************************************************************
        PARAMETER(a_day=4.17)
        PARAMETER(M1max=80,M2max=80,M1min=0.8,M2min=0.1,
     &            amax=10000,amin=3)
      real*8 rl1,rl2,rzamsf,q2,qnow,f,uBHL,Lx,nnn,downeq,Lx_rl,mip
      real*8 mt1,mt2,uuu,e,G,Ufun,dum,pmassc1,pmassc2,MASS_FUN
      real*8 L_Eddington,b_iso,m_transferrate
      real :: aaa(10)
      real :: bbb(10)
      DATA aaa /1.1, 1.5, 2, 3, 4 ,6, 8, 10,14,20/  

!       real::nnn0(20)=(/0.2,0.6,2,6,20/),
!       real::fff0(20)=(/50,55,60,65,70,75.0,
!      &	80,85,90,95,99/)
*
************************************************************************
* BSE parameters:
*
* neta is the Reimers mass-loss coefficent (neta*4x10^-13: 0.5 normally). 
* bwind is the binary enhanced mass loss parameter (inactive for single).
* hewind is a helium star mass loss factor (1.0 normally).
* alpha1 is the common-envelope efficiency parameter (1.0).  
* lambda is the binding energy factor for common envelope evolution (0.5).
*
* ceflag > 0 activates spin-energy correction in common-envelope (0). 
* tflag > 0 activates tidal circularisation (1).
* ifflag > 0 uses WD IFMR of HPE, 1995, MNRAS, 272, 800 (0). 
* wdflag > 0 uses modified-Mestel cooling for WDs (0). 
* bhflag > 0 allows velocity kick at BH formation (0). 
* nsflag > 0 takes NS/BH mass from Belczynski et al. 2002, ApJ, 572, 407 (1). 
* mxns is the maximum NS mass (1.8, nsflag=0; 3.0, nsflag=1). 
* idum is the random number seed used by the kick routine. 
*
* Next come the parameters that determine the timesteps chosen in each
* evolution phase:
*                 pts1 - MS                  (0.05) 
*                 pts2 - GB, CHeB, AGB, HeGB (0.01)
*                 pts3 - HG, HeMS            (0.02)
* as decimal fractions of the time taken in that phase.
*
* sigma is the dispersion in the Maxwellian for the SN kick speed (190 km/s). 
* beta is wind velocity factor: proportional to vwind**2 (1/8). 
* xi is the wind accretion efficiency factor (1.0). 
* acc2 is the Bondi-Hoyle wind accretion factor (3/2). 
* epsnov is the fraction of accreted matter retained in nova eruption (0.001). 
* eddfac is Eddington limit factor for mass transfer (1.0).
* gamma is the angular momentum factor for mass lost during Roche (-1.0). 
*
      e=2.71828d0
      G=6.67d-11
      neta = 0.5
      bwind = 0.0
      hewind = 1.0
      alpha1 = 3.0
      lambda = 1.5 
      ceflag = 0
      tflag = 1
      ifflag = 0 
      wdflag = 1 
      bhflag = 0
      nsflag = 1
      mxns = 3.0
      pts1 = 0.05
      pts2 = 0.01
      pts3 = 0.02
      sigma = 190.0
      beta = 2
      xi = 1.0 
      acc2 = 1.5
      epsnov = 0.001
      eddfac = 10.0
      gamma = -1.0
      mip=-1.7

*
* Set the seed for the random number generator. 
*
      idum = 3234
      if(idum.gt.0) idum = -idum
*
* Set the collision matrix.
*
      CALL instar
*
* Open the input file - list of binary initial parameters. 
***********************************************
       nm1=1d6
      !  *******************************************
c      OPEN(10,file='binaries.in',status='unknown')
c      READ(10,*)nm1
*
* Open the output files. 
*
      OPEN(11,file='binaries.out',status='unknown')
      OPEN(12,file='search_Lx38erg_beta1.out',status='unknown')
*
      do i = 1,nm1
        if(MOD(i,10000).eq.0.0)then
          print*,i
        endif

! 141             call random_number(aa)
!              mass0(1)=log(M1min)+(log(M1max)-log(M1min))*aa
!              mass0(1)=2.71828**mass0(1)
!              m1=mass0(1)
!              mass01=mass0(1)
             

!              call random_number(bb)
!              mass0(2)=log(M2min)+(log(M2max)-log(M2min))*bb
!              mass0(2)=2.71828**mass0(2)
!              m2=mass0(2)
!              mass02=mass0(2)

             q2=mass0(2)/mass0(1)
      !   new setting  TODO:
141             call random_number(aa)

             mass0(1)=  (aa*(M1min**mip-M1max**mip)
     &+M1max**mip)**(1/mip)
            call random_number(bb)
            mass0(2)=(mass0(1)-M2min)*bb+M2min
             call random_number(cc)
             a=log(amin)+(log(amax)-log(amin))*cc
             a=2.71828d0**a
             m1=mass0(1)
             mass01=mass0(1)
             m2=mass0(2)
             mass02=mass0(2)


        if(mass0(1).gt.100..or.mass0(1).lt.0.1)goto 141  
        if(mass0(2).gt.100..or.mass0(2).lt.0.1)goto 141
        rzams1=rzamsf(mass0(1))
        rzams2=rzamsf(mass0(2))
        
        rl1=rl(mass0(1)/mass0(2))*a
        rl2=rl(mass0(2)/mass0(1))*a
        
      !   mass0(1)>mass0(2)
        if(rzams1.gt.rl1.or.rzams2.gt.rl2.or.mass0(1).lt.mass0(2))then
            ! print*,i
           goto 141 
        endif
         ecc=0.0
         z=0.01
         tmax=1200.0
           tb=a_day**(-1.5d0)*mass0(1)**(-0.5d0)
     &               *(1d0+q2)**(-0.5d0)*a**(1.5d0)

*
* Read in parameters and set coefficients which depend on metallicity. 
*
c         READ(10,*)m1,m2,tb,ecc,z,tmax
         CALL zcnsts(z,zpars)
*
         ecc0 = ecc
         tb0 = tb/yeardy
         sep0 = aursun*(tb0*tb0*(mass(1) + mass(2)))**(1.d0/3.d0)
         tb0 = tb
*
* Initialize the binary. 
*
         kstar(1) = 1
         mass0(1) = m1
         mass(1) = m1
         massc(1) = 0.0
         ospin(1) = 0.0
         epoch(1) = 0.0
*
         kstar(2) = 1
         mass0(2) = m2
         mass(2) = m2
         massc(2) = 0.0
         ospin(2) = 0.0
         epoch(2) = 0.0
*
         tphys = 0.0
         tphysf = tmax
         dtp = 0.0


*
* Evolve the binary. 
*

         CALL evolv2(kstar,mass0,mass,rad,lum,massc,radc,
     &               menv,renv,ospin,epoch,tms,
     &               tphys,tphysf,dtp,z,zpars,tb,ecc)
*
* Search the BCM array for the formation of binaries of 
* interest (data on unit 12 if detected) and also output 
* the final state of the binary (unit 11). 
*
* In this example we will search for CVs. 
*
         jj = 0
         t1 = -1.0
         t2 = -1.0
         mt1=-1.0
 30      jj = jj + 1
* 检测最后一行
         if(bcm(jj,1).lt.0.0)then
             mt1=-1.0
             goto 40
         endif
         kw = INT(bcm(jj,2))
         kw2 = INT(bcm(jj,16))
         
*
         i1 = 15
         i2 = 29
      !    判断kw2是不是中子星
         if(kw.gt.kw2)then
            kdum = kw2
            kw2 = kw
            kw = kdum
            i2 = 15
            i1 = 29
            ! ratio of radius to roche lobe radius、
         endif
      !    kw2是中子星

      if(kw2.gt.12.and.kw2.lt.15.and.kw.lt.10.and.
     &bcm(jj,i1).gt.0.0.and.kw.gt.2)then
            if(mt1.lt.0.0d0)then
                  mt1=bcm(jj,1)
            else
                  mt1=mt2
                  mt2=bcm(jj,1)
                  if(mt2.gt.(mt1+tol))then
                        
                        f= bcm(jj,i1)*100
                        pmassc1=bcm(jj,i1-11)
                        pmassc2=bcm(jj,i2-11)
                        
                        ! print*,kw,kw2,bcm(jj,i2-9),bcm(jj,31)
                        ! nnn=3.0d0
                        if(bcm(jj,31).gt.0.and.f.gt.50)then
                              ! wrong caculation 
!       nnn=(2*G*bcm(jj,i2-11)/(10**(bcm(jj,i2-9))*6.955*10**8))**0.5
!      &/(2*3.1415926*(bcm(jj,31))/(bcm(jj,30)*6.955*10**8))
                              
                              nnn=SQRT(2*6.67*pmassc1*1.9891d+19/
     &(10**bcm(jj,i1-9)*6.955d+08))
     &/(2*3.1416*(bcm(jj,31)*6.955d+8/(bcm(jj,30)*3.1557d+07)))
                        !   check whether nnn is right    ! 
                        ! if(f.gt.50.and.nnn.gt.1.0.and.nnn.lt.40.0)then
                        !       print *,nnn
                        ! endif
                        else
                              
                              goto 30
                        endif
                        qnow=bcm(jj,i1-11)/bcm(jj,i2-11)
                        
                        
                        ! BHLmodel
!                         uBHL=(1+qnow)/qnow**3/
!      &(nnn*(1-f*rl(qnow))**beta*(1+
!      &(nnn*(1+qnow)*(1-f*rl(qnow))**beta/qnow)**1.5))
                        
                        if(f.lt.100.0d0)then
                              ! 风吸积
                              if(nnn.lt.1.1.or.nnn.gt.20)then
                                    
                                    goto 30
                                    
                              endif
                              

                              ! TODO:  uuu is ???
                              uuu=Ufun(f,nnn,beta,qnow)
                              if(uuu.gt.100)goto 30
! 260                           Lx=-uuu*nnn/2*1.0d6* bcm(jj,i1-1)
!      &*bcm(jj,i1-11)/((mt2-mt1)*1d6)
!       以 10**39 erg/s 为单位

                              ! print *,bcm(jj,i1-1),bcm(jj,i2-1)
                              mx=bcm(jj,i1-11)
                              mx2=bcm(jj,i2-11)
                              ! print *,bcm(jj,1),uuu
                              ! print *,uuu,Lx,bcm(jj,i1-1),bcm(jj,i1-11)
                              if(Lx.gt.0.1d0)then
                                    tbx=bcm(jj,30)*yeardy
                                    eccx=bcm(jj,32)
            ! print *,mass0,a
            ! print *,bcm(jj,i1-1)*bcm(jj,i1-11)/((mt2-mt1)*1d06)
                                    ! print *,uuu,nnn,Lx,bcm(jj,i1-1)
                                    ! print *,bcm(jj,i1-1),bcm(jj,i1-11)
                                    WRITE(12,113)m1,m2,Lx,ecc0,a
     &  ,tb0,mt1,mt2,kw,kw2,mx,mx2,eccx,tbx,f,nnn,uuu,bcm(jj,i1-1)
 
                              endif
!                               Lx_rl=mx*(3.0d8)**2*bcm(jj,i1-1)*1.9891d30/10**7/1.0d39
!                               if(Lx_rl.gt.0.01d0)then
!                                     tbx=bcm(jj,30)*yeardy
!                                     WRITE(12,112)m1,Lx,
!      & ecc0,tb0,mt1,mt2,kwx,kwx2,mx,mx2,eccx,tbx
!                               endif
                        endif

                  endif
            endif
        endif
            goto 30


         



! *  2020/03/25 remove from   line 282
            ! spline for uuu(eta)
!             if(kw2.eq.14)then
!                   if(INT(beta).lt.1.5)then
!                         ! NS1
!                         if(nnn.lt.6)then
!                               uuu=e**(25.0837
! & -1.2853*f + 0.021269*f*f - 0.00014105*f*f*f  
! &-2.9435*nnn + 1.7668*nnn**2 - 0.18687*nnn**3 + 3.1452d-7*f**4)
!                         else
!                               uuu=0.15
!                         endif
!                   else
!                         if(nnn.lt.12)then
! uuu=e**(-6.992 + 0.10317*f-0.0029074*f*f+ 4.6993e-05*f*f*f
! &-0.25939*nnn+0.23696*nnn **2 -0.014215*nnn **3-2.3307d-07*f **4)
!                         else
!                               uuu=0.15
!                         endif

!                   endif
!             endif
!             if(kw2.eq.14)then
!                   if(beta.lt.1.5)then
!                         if(nnn.lt.11.0)then
! uuu=e**(-31.5204 + 1.5286*f-0.033299*f*f+0.00032088*f*f*f 
! & -0.46158*nnn+ 0.3337*nnn **2-0.019411*nnn**3 - 1.1388d-06*f**4)
!                         else
!                               uuu=0.15
!                         endif
!                   else
!                         if(nnn.lt.15)then
! uuu=e**(-0.14552-0.55385*f+0.015067*f*f-0.00015465*f*f*f 
! & -0.58923*nnn+0.15002*nnn**2-0.0050181*nnn**3+5.6084d-07*f**4)
!                         else
!                               uuu=0.15
!                         endif
!                   endif
!             endif

            ! 原程序
!          if(kw.le.1.and.bcm(jj,i1).ge.1.0)then
!             if(kw2.ge.10.and.kw2.le.12)then
!                if(t1.lt.0.0)then
!                   ! 第一次筛选出来(第一行)
!                   t1 = bcm(jj,1)
!                   kwx = kw
!                   kwx2 = kw2
!                   mx = bcm(jj,i1-11)
!                   mx2 = bcm(jj,i2-11)
!                   tbx = bcm(jj,30)
!                   eccx = bcm(jj,32)
!                endif
!             endif
!          endif
! *
!          if(t1.gt.0.0.and.(bcm(jj,i1).lt.1.0.or.
!      &      kw.ne.kwx.or.kw2.ne.kwx2))then
!             if(t2.lt.0.0)then
!                t2 = bcm(jj,1)
!                if(t2.gt.(t1+tol))then
!                   f= bcm(jj,i2)
!                   qnow=bcm(jj,i1-11)/bcm(jj,i2-11)
                  
!                   nnn=3.0d0
!                   downeq=(1+
!      &(nnn*(1+qnow)*(1-f*rl(qnow))**beta/qnow)**1.5)
!                   uBHL=(1+qnow)/qnow**3/
!      &(nnn*(1-f*rl(qnow))**beta*downeq)
!                   Lx=uBHL*nnn/2/1.0d8*bcm(jj,i1-1)*bcm(jj,i1-11)

!                   ! Lx_rl=bcm(jj,i1-1)*bcm(jj,i1-11)*3.00d8**2/
!                   !    *1.0d39 erg/s
!                   if(Lx.gt.0.1d0)then 
!                         WRITE(12,112)m1,m2,ecc0,tb0,t1,t2,kwx,kwx2,
!      &                         mx,mx2,eccx,tbx
!                   endif
!                endif
!                t1 = -1.0
!                t2 = -1.0
!             endif
!          endif
! *
!          goto 30
 
! *
!          if(t1.gt.0.0)then
!             if(t2.lt.0.0) t2 = tmax
!             WRITE(12,112)m1,m2,ecc0,tb0,t1,t2,kwx,kwx2,mx,mx2,eccx,tbx
!          endif
! *
!          jj = jj - 1
!          kw = INT(bcm(jj,2))
!          kw2 = INT(bcm(jj,16))
!          mx = bcm(jj,4)
!          mx2 = bcm(jj,18)
!          tbx = bcm(jj,30)*yeardy
!          eccx = bcm(jj,32)
!          WRITE(11,111)tmax,kw,kw2,mx,mx2,eccx,tbx
40      continue

      enddo
*
 111  FORMAT(f10.1,2i3,3f8.3,1p,e14.6)
 112  FORMAT(3f8.3,1p,e14.6,0p,2f10.2,2i3,3f8.3,1p,e14.6)
 113  FORMAT(8f25.5,2i5,7f25.6,e30.6)
      CLOSE(10)
      CLOSE(11)
      CLOSE(12)
*
************************************************************************
*
      STOP
      END
***

