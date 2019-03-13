      subroutine wv_4_abs(taulyr,gw,ngw,npp,tlevel,plevel,
	1  ulayer,u_sat)

c     This subroutine calculates the absorption coefficient
c     fk0 at each model level, for each interval in g, called gw.
c     It also calculates the amount of solar flux in each gw.      

c     output variables
      integer pt
	parameter(pt=5)
	integer ngw,npp	,gas(4),xxx
      real taulyr(ngw,npp-1),fk(4,npp-1,ngw,pt),
	1          gw(ngw)	 ,ulayer(8,npp-1),tran(ngw,npp,npp),
	1          temp(npp-1)
c     input variables 
      real tlevel(npp), plevel(npp),tku,dx,tau(npp-1,1) 	,aaa

c     standard pressure levels
      real standp(21)
      data standp / 1.0169E-01, 2.0625E-01, 4.1834E-01,
     1              1.2180, 1.8075, 2.6824, 3.9806, 5.9072, 8.7662,
     1              13.0091, 19.3054, 28.6491, 42.5151, 63.0922,
     1              93.6284, 138.9440, 206.1920, 305.9876, 454.0837,
     1              673.8573, 1000.0000 /
c     interpolat1ion points as defined by preints
      integer inps(npp)  
c     work arrays
      real coefkl1(4,5,21,ngw,pt),wg(pt)
	real x1,x2,y1,y2,z1,z2,c1,c2,dip,player(npp-1),tlayer
	real u_sat

c     first set the interpolation between the pressure layers
c     (use gcm code)


	data wg /0.45162525 ,
     1   0.73635552  ,
     1   0.56888889  ,
     1   0.22090182  ,
     1   0.02222852 /


	do i=1,npp-1
	player(i)=(plevel(i+1)+plevel(i))/2.
	enddo

      call preints(inps,plevel,npp)

      open(8,file='abs6s.dat',status='old')

c     read in g intervals

      do j = 1, ngw
       read(8,*) gw(j)
      enddo
 55   format(F7.4,'  ')

c     read in coefficients for absorption fit
	read (8,*) ngas

	do i = 1, 21
	read(8,*)
	do igas =1,ngas
	read(8,*) gas(igas)
       do k = 1, ngw
       read(8,*)
	  do p=1,pt
        read(8,*) coefkl1(igas,1,i,k,p), coefkl1(igas,2,i,k,p),
	1             coefkl1(igas,3,i,k,p),
     3             coefkl1(igas,4,i,k,p), coefkl1(igas,5,i,k,p)
	  enddo 
       enddo
      enddo
	enddo
 	
c

cc     start loop over g intervals
      do 100 j = 1, ngw
	
c     interpolate abs coef (using gcm code)
c
c     should some better way be devised of dealing with pressures
c     above 1000 mb (for instance as in the way that p < 10mb is dealt with)?

      do k = 1, npp - 1
	
        tlayer = 0.5 * (tlevel(k) + tlevel(k+1))
        m = inps(k)
        n = m + 1
        dt = (tlayer - 250.0)

	do 101 igas=1,1 !ngas
	y1=0
 	y2=0
        if (m.gt.0)               then
	dip      = (player(k) - standp(m)) / (standp(n) - standp(m))
      do p=1,pt
	
         x1 = coefkl1(igas,1,m,j,p) + dt * (coefkl1(igas,2,m,j,p) + 
     1        dt * (coefkl1(igas,3,m,j,p) + dt *(coefkl1(igas,4,m,j,p)+ 
     1        dt * (coefkl1(igas,5,m,j,p)))))
         x2 = coefkl1(igas,1,n,j,p) + dt * (coefkl1(igas,2,n,j,p) + 
     1        dt * (coefkl1(igas,3,n,j,p) + dt *(coefkl1(igas,4,n,j,p)+ 
     1        dt * (coefkl1(igas,5,n,j,p)))))
c	y1=y1+ exp(-x1*ulayer(gas(igas),k))*wg(p)/2.0

c	y2=y2+ exp(-x2*ulayer(gas(igas),k))*wg(p)/2.0
	fk(igas,k,j,p)= x1 + (x2 - x1) * dip

	enddo

c      taulyr(j,k) =taulyr(j,k) -log(y1) + (-log(y2) + log(y1)) *dip*0.6
c	taulyr(j,k) =taulyr(j,k) - log(y1 + (y2 - y1) * dip ) !*tt1(igas,j)
	    

        else
	dip      = player(k) / standp(1)
      do p=1,pt
	
         x1 = coefkl1(igas,1,1,j,p) + dt * (coefkl1(igas,2,1,j,p) + 
     1        dt * (coefkl1(igas,3,1,j,p) + dt *(coefkl1(igas,4,1,j,p)+ 
     1        dt * (coefkl1(igas,5,1,j,p)))))
	y1=y1+ exp(-x1*ulayer(gas(igas),k))*wg(p)/2.
	fk(igas,k,j,p)= x1 * dip
	enddo
  
c        taulyr(j,k) = taulyr(j,k)-log(y1)* dip *0.6     

        endif
 101	enddo

      enddo
 100  continue
c
******
	do 111 ib = 1, ngw
	
	  do 112 igas=1, ngas
	tau1=0.
	tau0=0.
	  do j=1,npp-1
          temp(j)=0.
	
         enddo

	    do 113 ip = 1,pt
	

         do i=1,npp-1

	temp(i)=temp(i) +
	1 exp(-fk(igas,i,ib,ip)* ulayer(gas(igas),i))*wg(ip)/2

         enddo

	
 113		enddo



***********

	do i=1,npp-1
	taulyr(ib,i)=taulyr(ib,i)+temp(i)
c	print*,ib,i,tau(i,1)
	enddo

 112	  enddo
 111	enddo

      return
      end

********************************************************
      SUBROUTINE PREINTS (INPS, P, NP)
C
C----------------------------------------------------------------------C
C     THIS SUBROUTINE DETERMINES THE PRESSURE INTERPRETION POINTS      C
C     FOR ABSCOES                                                      C
C----------------------------------------------------------------------C
C
      REAL P(NP), STANDP(21)
      INTEGER INPS(NP)
C
      DATA STANDP / 1.0169E-01, 2.0625E-01, 4.1834E-01,
     1              1.2180, 1.8075, 2.6824, 3.9806, 5.9072, 8.7662,
     1              13.0091, 19.3054, 28.6491, 42.5151, 63.0922,
     1              93.6284, 138.9440, 206.1920, 305.9876, 454.0837,
     1              673.8573, 1000.0000 /
C
C----------------------------------------------------------------------C
C     STANDP ARE REVERSED TO THOSE IN ABSCOES                          C
C----------------------------------------------------------------------C
C
      JENDS = 21
      DO 400 K = 1, NP
          INPS(K) =  0
C
        DO 200 J = 1, JENDS
          IF (P(K) .GE. STANDP(J))                                THEN
            INPS(K) =  INPS(K) + 1
          ENDIF
  200   CONTINUE
C
  400 CONTINUE
C
      RETURN
      END

