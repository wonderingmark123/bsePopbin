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
      integer i,j,k,jj,nm1,mass0INT,q2INT,seperationINT
      integer kw,kw2,kwx,kwx2,kstar(2)
      integer i1,i2,kdum,anotherloop,isht,ALLnumMass
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
      character(len=20) :: Ctemp
      CHARACTER*8 label(14)
      CHARACTER*4 Outname
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
      eddfac = 10000.0
      gamma = -1.0
      
      Msun=1.9891d30
      label(1) = 'INITIAL '
      label(2) = 'KW CHNGE'
      label(3) = 'BEG RCHE'
      label(4) = 'END RCHE'
      label(5) = 'CONTACT '
      label(6) = 'COELESCE'
      label(7) = 'COMENV  '
      label(8) = 'GNTAGE  '
      label(9) = 'NO REMNT'
      label(10) = 'MAX TIME'
      label(11) = 'DISRUPT '
      label(12) = 'BEG SYMB'
      label(13) = 'END SYMB'
      label(14) = 'BEG BSS'

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
       ALLnumMass=1000
       Outname='1126'
      !  *******************************************
c      OPEN(10,file='binaries.in',status='unknown')
c      READ(10,*)nm1
*
* Open the output files. set names for them
*

      OPEN(14,file='hanchen_first_neutron'
     &      //Outname//'.out',status='unknown')
*     
            isht=0
            
      
      ! do mass0INT = int(M1min),int(M1max),1

      do mass0INT = 0,ALLnumMass
            
            
      do q2INT=0,93
            
      do seperationINT=0,100
            cc=seperationINT*0.01
            
            mass0(1)=exp(real(mass0INT)/real(ALLnumMass)
     &            *(log(M1max)-log(M1min))+log(M1min))

            q2=q2INT*0.1+0.08
            
            isht=isht+1
        if(MOD(isht,10000).eq.0.0)then
          print*,isht
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

            !  q2=mass0(2)/mass0(1)
 
! 141             call random_number(aa)
! ! check!!!!!
!             mip=-1.3
!              mass0(1)=  (aa*(M1max**mip-M1min**mip)
!      &+M1min**mip)**(1/mip)
!             call random_number(bb)
        
            mass0(2)=mass0(1)*q2
            
            ! TODO: P(logP)~ (logP)^-0.55 P(e)~e^-0.42
            !  call random_number(cc)
             a=log(amin)+(log(amax)-log(amin))*cc
             a=2.71828d0**a
             m1=mass0(1)
             mass01=mass0(1)
             m2=mass0(2)
             mass02=mass0(2)
             write(Ctemp,"(I20)") isht
!             OPEN(16,file='DataGedianz0.01_'//Outname//'/'
!      &      //trim(adjustl(Ctemp))//'.out',status='unknown')
!             OPEN(17,file='DataGedianz0.01_'//Outname//'/'
!      &           //trim(adjustl(Ctemp))//'bpp.out',status='unknown')
        if(bb.lt.0.08)continue
        
        rzams1=rzamsf(mass0(1))
        rzams2=rzamsf(mass0(2))
        
        rl1=rl(mass0(1)/mass0(2))*a
        rl2=rl(mass0(2)/mass0(1))*a
        
      !   mass0(1)>mass0(2)
        if(rzams1.gt.rl1.or.rzams2.gt.rl2.or.mass0(1).lt.mass0(2))then
            ! print*,i
      !      goto 141 
            WRITE(11,101)m1,q2,a
            continue
        endif
         ecc=0.0
         z=0.01
      !    tmax=200.0
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
      if(bcm(1,1).lt.0.0) goto 40
         jj = 0
         t1 = -1.0
         t2 = -1.0
         mt1=-1.0
    


 30      jj = jj + 1
* 检测最后一行
            CLOSE(17)
            if(jj > 50000)goto 40
         if(bcm(jj,1).lt.0.0)then
             mt1=-1.0
             goto 40
         endif

         kw = INT(bcm(jj,2))
         kw2 = INT(bcm(jj,16))
         

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
      if(kw2.eq.12.and.bcm(jj,31).gt.0)then
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
                        ! m_transferrate=bcm(jj,i2-1)
                        mx=bcm(jj,i1-11)
                        mx2=bcm(jj,i2-11)
                        ! m_tr=abs(bcm(jj,i2-1))
                        m_tr=(bcm(jj,i2-1))
                        tbx=bcm(jj,30)
                        
                        eccx=bcm(jj,32)
                        WRITE(14,113)m1,m2,Lx,ecc0,a
     &  ,tb0,mt1,mt2,kw,kw2,mx,mx2,eccx,tbx,
     &     f,m_transferrate,uuu,m_tr,isht
                        goto 40
                  endif
            endif
        endif
            goto 30


         


40      continue

      enddo
      enddo
      enddo
*
 111  FORMAT(f10.1,2i3,3f8.3,1p,e14.6)
 99   FORMAT(f10.4,2i3,10f10.4,5e12.4,f7.3)
 112  FORMAT(3f8.3,1p,e14.6,0p,2f10.2,2i3,3f8.3,1p,e14.6)
 113  FORMAT(2f25.5,e30.6,5f30.6,2i5,5f25.6,3e30.6,i15)
 100  FORMAT(f11.4,2f9.3,2i3,f13.3,f6.2,2f8.3,2x,a8)
 101  FORMAT(3f11.4)
      CLOSE(10)
      CLOSE(11)
      CLOSE(12)
*
************************************************************************
*
      STOP
      END
***

