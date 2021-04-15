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
      integer i1,i2,kdum,anotherloop
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
        PARAMETER(M1max=150,M2max=80,M1min=2,M2min=0.1,
     &            amax=10000,amin=3)
      real*8 rl1,rl2,rzamsf,q2,qnow,f,uBHL,Lx,nnn,downeq,Lx_rl,mip
      real*8 mt1,mt2,uuu,e,G,Ufun,dum,pmassc1,pmassc2,MASS_FUN
      real*8 L_Eddington,b_iso,m_transferrate,m_tr,Msun,L_app,dd
      
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
      ceflag = 1
      tflag = 0
      ifflag = 0 
      wdflag = 1 
      bhflag = 1
      nsflag = 1
      mxns = 3.0
      pts1 = 0.05
      pts2 = 0.01
      pts3 = 0.02
      sigma = 265.0
      beta = 2
      xi = 1.0 
      acc2 = 1.5
      epsnov = 0.001
      eddfac = 10.0
      gamma = -1.0
      
      Msun=1.9891d30

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
       nm1=1d5
      !  *******************************************
c      OPEN(10,file='binaries.in',status='unknown')
c      READ(10,*)nm1
*
* Open the output files. set names for them
*
      OPEN(11,file='binaries.out',status='unknown')
      OPEN(12,file='all_Lx37erg_z0.01_wind_MS2.out',status='unknown')
      OPEN(13,file='all_Lx37erg_z0.01_rb_MS2.out',status='unknown')
      OPEN(14,file='first_Lx37erg_z0.01_wind_MS2.out',status='unknown')
      OPEN(15,file='first_Lx37erg_z0.01_rb_MS2.out',status='unknown')
*
      do mass0(1) = M1min,M1max,1
        do q=0.08,1,0.1
            do cc=0,1,0.01
        if(MOD(i,10000).eq.0.0)then
          print*,i
        endif
        anotherloop=1
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
 
! 141             call random_number(aa)
! ! check!!!!!
!             mip=-1.3
!              mass0(1)=  (aa*(M1max**mip-M1min**mip)
!      &+M1min**mip)**(1/mip)
!             call random_number(bb)
!             mass0(2)=mass0(1)*bb
            
!             ! TODO: P(logP)~ (logP)^-0.55 P(e)~e^-0.42
!              call random_number(cc)
             a=log(amin)+(log(amax)-log(amin))*cc
             a=2.71828d0**a
             m1=mass0(1)
             mass01=mass0(1)
             m2=mass0(2)
             mass02=mass0(2)


        if(bb.lt.0.08)goto 141
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
         tmax=10000.0
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
      !    bcm(jj,i2-1).gt.0
      if(kw2.gt.12.and.kw2.lt.15.and.
     &bcm(jj,i1).gt.0.0.and.bcm(jj,31).gt.0)then
            if(mt1.lt.0.0d0)then
                  mt1=bcm(jj,1)
            else
                  mt1=mt2
                  mt2=bcm(jj,1)
                  
                  if(mt2.gt.(mt1+tol))then
                        !  f=洛希搬填充率百分之n
                        f= bcm(jj,i1)*100
                        qnow=bcm(jj,i1-11)/bcm(jj,i2-11)
                        L_Eddington=2.6d38*pmassc2
                        m_transferrate=bcm(jj,i1-1)
                        mx=bcm(jj,i1-11)
                        mx2=bcm(jj,i2-11)
                        m_tr=abs(bcm(jj,i1-1))
                        if(f.lt.100.0d0)then
                              ! 风吸积 未充满洛希搬
                              !  Mc^2 first
                              L_Eddington=2.6d38*mx2
      m_transferrate=m_tr/(L_Eddington/1d7/(3.0d8)**2/Msun*yearsc)
                              if(m_transferrate.gt.1)then
                              Lx=L_Eddington*(1+log(m_transferrate))
                              else
                              Lx=0.1*(3.0d8)**2*m_tr*Msun*10**7/yearsc
                              endif
                              ! ********************************
                              ! TODO: using 0.1*Mc^2 to check
                              Lx=0.1*(3.0d8)**2*m_tr*Msun*10**7/yearsc


                              ! more range of Lx
                              if(Lx.gt.1d37)then
                                    tbx=bcm(jj,30)*yeardy
                                    eccx=bcm(jj,32)
                                    if(anotherloop.eq.1)then
                                          WRITE(14,113)m1,m2,Lx,ecc0,a
     &  ,tb0,mt1,mt2,kw,kw2,mx,mx2,eccx,tbx,f,m_transferrate,uuu,m_tr
                                          anotherloop=0
                                    endif
                                    WRITE(12,113)m1,m2,Lx,ecc0,a
     &  ,tb0,mt1,mt2,kw,kw2,mx,mx2,eccx,tbx,f,m_transferrate,uuu,m_tr
 
                              endif
                        endif

                        if(f.gt.100.0d0)then
                              ! lobe吸积 充满洛希搬
                              L_Eddington=2.6d38*mx2
      m_transferrate=m_tr/(L_Eddington/1d7/(3.0d8)**2/Msun*yearsc)
      ! print*,m_tr,m_transferrate

                              
                              if(m_transferrate.gt.1)then
                              Lx=L_Eddington*(1+log(m_transferrate))
                              else
                                    Lx=L_Eddington*m_transferrate
                              endif

                              ! to get b_iso 
                              
                              if(m_transferrate.gt.8.5)then
                                    b_iso=73/(m_transferrate)**2
                              else
                                    b_iso=1
                              endif
                              L_app=Lx/b_iso

                              ! TODO:!using 0.1*Mc^2 to check
                              Lx=(3.0d8)**2*m_tr*Msun*10**7/yearsc


                              ! if(L_app.gt.1d37)then
                              if(Lx.gt.1d37)then
                                    tbx=bcm(jj,30)*yeardy
                                    eccx=bcm(jj,32)
                                    if(anotherloop.eq.1)then
                                          WRITE(15,113)m1,m2,Lx,ecc0,a
     &  ,tb0,mt1,mt2,kw,kw2,mx,mx2,eccx,tbx,f,m_transferrate,b_iso,m_tr
                                          anotherloop=0
                                    endif
                                    WRITE(13,113)m1,m2,Lx,ecc0,a
     &  ,tb0,mt1,mt2,kw,kw2,mx,mx2,eccx,tbx,f,m_transferrate,b_iso,m_tr
 
                              endif
                        endif
                  
                  endif
            endif
        endif
            goto 30


         


40      continue

      enddo
*
 111  FORMAT(f10.1,2i3,3f8.3,1p,e14.6)
 112  FORMAT(3f8.3,1p,e14.6,0p,2f10.2,2i3,3f8.3,1p,e14.6)
 113  FORMAT(2f25.5,e30.6,5f30.6,2i5,5f25.6,3e30.6)
      CLOSE(10)
      CLOSE(11)
      CLOSE(12)
*
************************************************************************
*
      STOP
      END
***

