
c     atmospheric conditions to be filled by 'include "mls75.pro"'
c     for np different levels
      parameter (np=60)
      parameter(vstar=990., nband=1, nv=200000, dv=.001)
c      parameter(vstar=1297.9, nband=1, nv=133100, dv=.001)
c     ngw is the number of gw intervals used      
      parameter(ngw=6)
      real  gw(ngw)
c
      real xlayer1(np),xlayer2(np),xlayer3(np),xlayer4(np),
     1     xlayer6(np),xlevel1(np+1),xlevel2(np+1) 
c     defined by .pro files
      real plevel(np+1),tlayer(np)
      real wlayer(np),clayer,olayer(np),gas_lay(7,np),
	1           ulayer(8,np)	,gas_layer(np,7)
c     interpolated values
      real player(np),tlevel(np+1), flx(np+1), coolr(np),
     1     coolrg(np,ngw)
	real vm,tku ,aaa,dx,ts(np)
c
      real tflx(np+1),heatr(np)
c
c-----temporary arrays
c     taulyr = optical depth of each layer
c     tran = transmission
c     player = pressure at layers
c
      real taulyr(ngw,np),tmean(np+1),blayer(np+1),bblayer(np+1)
      real tran(np+1,np+1),tran1(np+1,np+1),fluxu(np+1),
	1           fluxd(np+1),flxu(np+1),flxd(np+1)
      real v(nv) ,srf(nv+1),id(nv),t1,t2
	real ggg(np),ssa(np) ,ee ,LWC,c_total(np,4),U_Sat,phi
	  
c     set v
     
      do iv = 1, nv
        v(iv) = vstar + dv*((iv-1)+0.5)
      enddo
c
      npp1 = np+1
	ee=1.0
	u_sat=1 !0.7033948
c
c      include "mls60.pro"
c      tsfc = 294.0
c      include "saw60.pro"
c      tsfc = 257.1
      include "tro60.pro"
      tsfc = 300.0
c	include "mls36.pro"
c	tsfc = 302.75
c

c     open input file
	open(342,file='srf_09.dat',status='old')
c	open(343,file='id.dat',status='old')

	  read(342,*)
      do iv = 1, nv+1
        read(342,*) SRF(iv)
      enddo

c-----compute mean layer pressure
c
      do i=1,np
       player(i)=0.5*(plevel(i)+plevel(i+1))
      enddo
c
c-----compute mean level temperature
c
      do i = 2,np
       tlevel(i) = 0.5*(tlayer(i-1)+tlayer(i))
      enddo
c     set boundary levels for temperature

      tlevel(1) = tlayer(1)
      tlevel(np+1) = tsfc

c-----compute the absorption coefficient at "nv" spectral points
c     mid=1, 2, 7 for water vapor, co2 and o2 absorption, respectively.
c     absgas is the absorption ceofficient in (cm**2/molecules)

c-----assign absorber amount


      do k=1,np
	
	ulayer(1,k)=	wlayer(k)*1.0204 * (plevel(k+1) - plevel(k)) 
	ulayer(2,k)=	350.* 1.E-06 * 44. / 28.97
     *   *1.0204 * (plevel(k+1) - plevel(k))
	ulayer(3,k)=	olayer(k)*1.0204 * (plevel(k+1) - plevel(k))
	ulayer(8,k)= 28.97/18.0047*wlayer(k)**2
     &   *1.0204 * (plevel(k+1) - plevel(k)) 
      enddo


c-----compute the flux transmittance between level i and j
c     for the spectral band 'iband'

        do i=1,npp1
         fluxu(i)=0.
         fluxd(i)=0.
        enddo

c-----compute planck flux.  units are w/m**2
      do i=1,np
       tmean(i)=tlayer(i)
       blayer(i)=0.
	c_total(i,1)=0.
	c_total(i,2)=0.
	c_total(i,3)=0.
	c_total(i,4)=0.
      enddo
      tmean(np+1)=tsfc

      do j=1,nv
       do i=1,np+1
	vm=v(j)
       blayer(i)=blayer(i)+3.74e-8*vm*vm*vm/(exp(1.4385*vm/tmean(i))-1.)
     *            *dv*srf(j)
c        blayer(i)=blayer(i)+1.191e-5*vm**3 /(exp(1.4385*vm/tmean(i))-1.)
c     1      *dv*0.001  * (SRF(j) + SRF(j+1))/2 
       enddo
      enddo
c	do i=1,61
c      print*,i, blayer(i)
c	enddo
c	call wv_abs(taulyr,gw,npp1,tlayer,plevel,gas_lay)
	call wv_4_abs(taulyr,gw,ngw,np+1,tlevel,plevel,ulayer,u_sat)
c	call Habs(taulyr, gw,np,ngw,gas_layer,tlayer,plevel*100)

      do 1000 ib=1, ngw

ccc-----compute net downward flux.

      call abap4 (blayer,ee,ggg,taulyr(ib,:),ssa,flxu,flxd,npp1)
c	 call abap4i (ee,c_total,taulyr(ib,:),SSA,npp1,U_Sat,blayer,
c	1 blayer(npp1),Flxu)

       do i=1,npp1
        fluxu(i)=fluxu(i)+flxu(i) *gw(ib)
        fluxd(i)=fluxd(i)+flxd(i) *gw(ib)
       enddo
      do i=1,np
       coolrg(i,ib)=(flxu(i+1)-flxu(i)+flxd(i)-flxd(i+1))*8.441874/
     *              (plevel(i+1)-plevel(i)) *gw(ib)
      enddo

1000  continue

c	do i=1,npp1
c	write (9,*) plevel(i),fluxu(i)/phi
c	enddo
      write (9,2000) (plevel(i),fluxu(i),fluxd(i),i=1,npp1)
2000  format(1x,f12.4,1x,f12.4,1x,f12.4)

c-----compute cooling rate profile

      do i=1,np
       coolr(i)=(fluxu(i+1)-fluxu(i)+fluxd(i)-fluxd(i+1))*8.441874/
     *          (plevel(i+1)-plevel(i))
      enddo
      call CPU_TIME(t2)
       print *, t2

      end

**************************************************************
 
