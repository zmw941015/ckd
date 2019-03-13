      SUBROUTINE abap4 (bf,ee,ggg,tau,omg,flxu,flxd,n)
!bf 普朗克发射率   ee地表发射率

c      include 'para.file'           
      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
      parameter (Pi=3.1415927)
!      REAL TAUA(n-1), TAUR(n-1), TAUG(n-1), TAUOMA(n-1),	 &
!     TAUC(n-1), TAUOMC(n-1),						 &
!     TAUOMGC(n-1,4), TAUOMGA(n-1,4),CLD(n-1)

      real bf(n),tt(n),ww(n),u,fd1(n),fu1(n),bts,g(n),f(n),uu,
	1  fdu(n),fuu(n),u1,u2,fdu1(n),fuu1(n),fdu2(n),fuu2(n),aa(n),bb(n),
     1  cc(n),dd(n),aa1(n),bb1(n),cc1(n),dd1(n),mm(n),mm1(n),ff(n),
     1  gg(n),hh(n),a(n),b(n),c(n),d(n),m(n),j(n),k(n),h(n),aa2(n),
     1  bb2(n),cc2(n),dd2(n), mmu(n),mmuu(n),mmu1(n),mmu2(n),nn1(n),
     1  nn2(n),nn3(n),nn4(n),nn33(n),nn22(n),nn5(n),nn6(n),nn7(n),
     1  nn8(n),nn9(n),nn10(n),nn11(n),nn12(n),nn13(n),nn14(n),nn15(n),
     1  nn16(n),nn17(n),nn18(n),nn19(n),nn20(n),nn21(n),nn23(n),nn24(n),
	1  nn25(n),nn26(n),nn27(n),nn28(n),mm13(n),mm12(n),mm11(n),mm10(n),
	1  mm9(n),mm8(n),mm6(n),mm5(n),mm4(n),mm3(n),mm2(n),mm15,mm14,mm7,
	1  mm21,start,finish,mm4u1(n),mm4u2(n),mm6u1(n),mm6u2(n),mmuu1(n),
	1  mmuu2(n),mm10u1(n),mm10u2(n),nn1u1(n),nn1u2(n),nn33u1(n),
	1  nn33u2(n),nn2u1(n),nn2u2(n),nn22u1(n),nn22u2(n),nn5u1(n),
	1  nn5u2(n), nn6u1(n),nn6u2(n),nn9u1(n),nn9u2(n),nn10u1(n),
	1  nn10u2(n),nn11u1(n),nn11u2(n),nn12u1(n),nn12u2(n),mm9u1(n),
	1  mm9u2(n),nn3u1(n),nn3u2(n),nn4u1(n),nn4u2(n),nn13u1(n),
	1  nn13u2(n),nn14u1(n),nn14u2(n),nn15u1(n),nn15u2(n),nn16u1(n),
	1  nn16u2(n),nn17u1(n),nn17u2(n),fduu1(n),fuuu1(n),fduu2(n),
	1  fuuu2(n),mm14u1,mm14u2,mm15u1,mm15u2,mm1u1(n),mm1u2(n),mm8u1(n),
	1  mm8u2(n),nn7u1(n),nn7u2(n),nn8u1(n),nn8u2(n),nn18u1(n),
	1  nn18u2(n),nn19u1(n),nn19u2(n),nn20u1(n),nn20u2(n),nn29u1(n),
	1  nn29u2(n),nn30u1(n),nn30u2(n),nn21u1(n),nn21u2(n),nn23u1(n),
	1  nn23u2(n),nn24u1(n),nn24u2(n),nn25u1(n),nn25u2(n),nn26u1(n),
	1  nn26u2(n),nn27u1(n),nn27u2(n),nn28u1(n),nn28u2(n),aa3(n),
	1  aa4(n),bb3(n),bb4(n),cc3(n),cc4(n),dd3(n),dd4(n),bs,r,t0,uii,
	1  ttall(n),u0,f0,f1(n),flxu(n),flxd(n)

      real tau(n-1),omega(n-1),ww1(n-1),ggg(n-1),omg(n-1) ,ee

!--------------------------wk
!	do i = 1, n-1
!            tt(i)         =  TAUA(I)  + TAUG(I) 
!            OMARS         =  MAX(TAUOMA(I), 10E-20) 
!            ww(i)        =  OMARS / tt(i)
!           g(i) = TAUOMGA(I,1) / OMARS
!	  IF(CLD(I).GE.CUT)	  THEN
!	        tt(i)         =  TAUA(I)  + TAUG(I) + TAUC(I) 
 !           OMARCS        =  TAUOMC(I) 
  !          ww(i)        =  OMARCS / tt(i)
   !         g(i) = TAUOMGC(I,1) / OMARCS
	!  ENDIF	    
!	      ww(i) = MIN(ww(i),0.999999) 	
!	enddo
!------------------------------wk

      do i=1,n-1
      tt(i)=tau(i)
      ww(i)=omg(i)
      g(i)=ggg(i)
      enddo

!aa(1)=tt(1)
!do i=2,n-1
!aa(i)=tt(i)-tt(i-1)
!enddo
!do i=1,n-1
!tt(i)=aa(i)
!enddo
      do i=1,n
      if (tt(i)==0)	then
      tt(i)=1e-15
      end if
      end do

      do i=1,n
      if (ww(i)==0)	then
      ww(i)=0.0000001
      end if
      end do
!ee=1.0
!ccccccccccccccccccccccccc 2Sdeta函数调整 cccccccccccccccccccccccccccccccccccccccccccccccc
      do i=1,n-1
      f(i)=g(i)*g(i)
      tt(i)=(1-ww(i)*f(i))*tt(i)
      ww(i)=((1-f(i))*ww(i))/(1-ww(i)*f(i))
      g(i)=g(i)/(1+g(i))
      end do
!ccccccccccccccccccccccccc 基于AA的1阶变分迭代参数(相函数作4S展开) ccccccccccccccccccccccccccccccccccccccccccccccc
      u1=1/4.7320500
      u2=1/1.2679492
      mm14u1=u1*u1
      mm15u1=u1*u2
      mm14u2=u2*u1
      mm15u2=u2*u2
      fd1(1)=0
      fduu1(1)=0
      fduu2(1)=0
      fdu1(1)=0
      fdu2(1)=0
      do i=1,n-1
c	print*,i,bf(i+1),bf(i),tt(i)
      mm(i)=log(bf(i+1)/bf(i))/tt(i)
      mm1u1(i)=-tt(i)/u1
      mm1u2(i)=-tt(i)/u2
      mm2(i)=1-ww(i)
      mm3(i)=mm2(i)*tt(i)
      mm4u1(i)=u1*mm(i)
      mm4u2(i)=u2*mm(i)
      mm5(i)=u1*mm(i)
      mm6(i)=u2*mm(i)
      mm8u1(i)=u1*mm2(i)
      mm8u2(i)=u2*mm2(i)
      mmu1(i)=-mm3(i)/u1
      mmu2(i)=-mm3(i)/u2
      mm9u1(i)=exp(mm1u1(i))
      mm9u2(i)=exp(mm1u2(i))
      mm11(i)=exp(mmu1(i))
      mm12(i)=exp(mmu2(i))
      mm13(i)=3*g(i)
c	print*,i, mm4u1(i),mm4u2(i)
      nn33u1(i)=mm2(i)/(-mm4u1(i)-1)
      nn33u2(i)=mm2(i)/(-mm4u2(i)-1)
      nn22u1(i)=mm2(i)/(mm4u1(i)-1)
      nn22u2(i)=mm2(i)/(mm4u2(i)-1)
      nn3(i)=mm2(i)/(-mm5(i)-mm2(i))
      nn4(i)=mm2(i)/(mm5(i)-mm2(i))
      nn5(i)=mm2(i)/(-mm6(i)-mm2(i))
      nn6(i)=mm2(i)/(mm6(i)-mm2(i))
      nn7u1(i)=u1/(2-ww(i))
      nn7u2(i)=-u2*u1/(-mm8u2(i)-u1)
      nn8u1(i)=-u1*u2/(-mm8u1(i)-u2)
      nn8u2(i)=u2/(2-ww(i))
      nn9u1(i)=u1/ww(i)
      nn9u2(i)=u2*u1/(-mm8u2(i)+u1)
      nn10u1(i)=u1*u2/(-mm8u1(i)+u2)
      nn10u2(i)=u2/ww(i)
      nn15u1(i)=-u1/(-mm4u1(i)-1)
      nn15u2(i)=-u2/(-mm4u2(i)-1)
      nn16u1(i)=u1/(mm4u1(i)-1)
      nn16u2(i)=u2/(mm4u2(i)-1)
      nn17u1(i)=1-mm9u1(i)*mm11(i)
      nn17u2(i)=1-mm9u2(i)*mm11(i)
      nn18u1(i)=1-mm9u1(i)*mm12(i)
      nn18u2(i)=1-mm9u2(i)*mm12(i)
      nn19u1(i)=mm11(i)-mm9u1(i)
      nn19u2(i)=mm11(i)-mm9u2(i)
      nn20u1(i)=mm12(i)-mm9u1(i)
      nn20u2(i)=mm12(i)-mm9u2(i)
      nn29u1(i)=bf(i+1)-bf(i)*mm9u1(i)
      nn29u2(i)=bf(i+1)-bf(i)*mm9u2(i)
      nn30u1(i)=bf(i+1)*mm9u1(i)-bf(i)
      nn30u2(i)=bf(i+1)*mm9u2(i)-bf(i)
      nn21u1(i)=nn15u1(i)*nn29u1(i)
      nn21u2(i)=nn15u2(i)*nn29u2(i)
      nn23u1(i)=nn16u1(i)*nn30u1(i)
      nn23u2(i)=nn16u2(i)*nn30u2(i)
      nn24u1(i)=ww(i)/(4*u1)
      nn24u2(i)=ww(i)/(4*u2)
      nn25u1(i)=1-mm13(i)*mm14u1
      nn25u2(i)=1-mm13(i)*mm14u2
      nn26u1(i)=1-mm13(i)*mm15u1
      nn26u2(i)=1-mm13(i)*mm15u2
      nn27u1(i)=1+mm13(i)*mm14u1
      nn27u2(i)=1+mm13(i)*mm14u2
      nn28u1(i)=1+mm13(i)*mm15u1
      nn28u2(i)=1+mm13(i)*mm15u2
      end do
!cccccccccccccccccccccccc 吸收近似（在4S展开中使用） cccccccccccccccccccccccccccccccccccccccccccc
      do i=1,n-1
      fduu1(i+1)=fduu1(i)*mm11(i)-nn3(i)*(bf(i+1)-bf(i)*mm11(i))
      end do
      do i=1,n-1
      fduu2(i+1)=fduu2(i)*mm12(i)-nn5(i)*(bf(i+1)-bf(i)*mm12(i))
      end do
      fuuu1(n)=(1-ee)*(u1*fduu1(n)+u2*fduu2(n))+ee*bf(n)
      fuuu2(n)=fuuu1(n)
      do i=1,n-1
      fuuu1(n-i)=fuuu1(n-i+1)*mm11(n-i)-nn4(n-i)*(-bf(n-i+1)*mm11(n-i)
	1           +bf(n-i))
      end do
      do i=1,n-1
      fuuu2(n-i)=fuuu2(n-i+1)*mm12(n-i)-nn6(n-i)*(-bf(n-i+1)*mm12(n-i)
	1           +bf(n-i))
      end do
!ccccccccccccccccccccccc 基于AA的变分迭代算法（4S相函数作4S展开）有调整ccccccccccccccccccccccccccccccccccc
      do i=1,n-1
      aa1(i)=fuuu1(i+1)*nn7u1(i)*nn17u1(i)
	1       -nn4(i)*(nn21u1(i)-bf(i+1)*nn7u1(i)*nn17u1(i))
      bb1(i)=fuuu2(i+1)*nn8u1(i)*nn18u1(i)
	1       -nn6(i)*(nn21u1(i)-bf(i+1)*nn8u1(i)*nn18u1(i))
      cc1(i)=fdu1(i)*nn9u1(i)*nn19u1(i)
	1       -nn3(i)*(nn21u1(i)-bf(i)*nn9u1(i)*nn19u1(i))
      dd1(i)=fdu2(i)*nn10u1(i)*nn20u1(i)
	1       -nn5(i)*(nn21u1(i)-bf(i)*nn10u1(i)*nn20u1(i))
      aa3(i)=fuuu1(i+1)*nn7u2(i)*nn17u2(i)
	1       -nn4(i)*(nn21u2(i)-bf(i+1)*nn7u2(i)*nn17u2(i))
      bb3(i)=fuuu2(i+1)*nn8u2(i)*nn18u2(i)
	1       -nn6(i)*(nn21u2(i)-bf(i+1)*nn8u2(i)*nn18u2(i))
      cc3(i)=fdu1(i)*nn9u2(i)*nn19u2(i)
	1       -nn3(i)*(nn21u2(i)-bf(i)*nn9u2(i)*nn19u2(i))
      dd3(i)=fdu2(i)*nn10u2(i)*nn20u2(i)
	1       -nn5(i)*(nn21u2(i)-bf(i)*nn10u2(i)*nn20u2(i))
      fdu1(i+1)=fdu1(i)*mm9u1(i)
	1       -nn33u1(i)*nn29u1(i)+nn24u1(i)*(nn25u1(i)*aa1(i)
     1       +nn26u1(i)*bb1(i)+nn27u1(i)*cc1(i)+nn28u1(i)*dd1(i))
      fdu2(i+1)=fdu2(i)*mm9u2(i)
	1       -nn33u2(i)*nn29u2(i)+nn24u2(i)*(nn25u2(i)*aa3(i)
     1       +nn26u2(i)*bb3(i)+nn27u2(i)*cc3(i)+nn28u2(i)*dd3(i))
      end do 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      fuu1(n)=(1-ee)*(u1*fdu1(n)+u2*fdu2(n))+ee*bf(n)
      fuu2(n)=fuu1(n)
      do i=1,n-1
      aa2(n-i)=fuu1(n-i+1)*(-nn9u1(n-i))*(-nn19u1(n-i))
	1  -nn4(n-i)*(nn23u1(n-i)-bf(n-i+1)*(-nn9u1(n-i))*(-nn19u1(n-i)))
      bb2(n-i)=fuu2(n-i+1)*(-nn10u1(n-i))*(-nn20u1(n-i))
	1  -nn6(n-i)*(nn23u1(n-i)-bf(n-i+1)*(-nn10u1(n-i))*(-nn20u1(n-i)))
      cc2(n-i)=fdu1(n-i)*(-nn7u1(n-i))*(-nn17u1(n-i))
	1  -nn3(n-i)*(nn23u1(n-i)-bf(n-i)*(-nn7u1(n-i))*(-nn17u1(n-i)))
      dd2(n-i)=fdu2(n-i)*(-nn8u1(n-i))*(-nn18u1(n-i))
	1  -nn5(n-i)*(nn23u1(n-i)-bf(n-i)*(-nn8u1(n-i))*(-nn18u1(n-i)))
      aa4(n-i)=fuu1(n-i+1)*(-nn9u2(n-i))*(-nn19u2(n-i))
	1  -nn4(n-i)*(nn23u2(n-i)-bf(n-i+1)*(-nn9u2(n-i))*(-nn19u2(n-i)))
      bb4(n-i)=fuu2(n-i+1)*(-nn10u2(n-i))*(-nn20u2(n-i))
	1  -nn6(n-i)*(nn23u2(n-i)-bf(n-i+1)*(-nn10u2(n-i))*(-nn20u2(n-i)))
      cc4(n-i)=fdu1(n-i)*(-nn7u2(n-i))*(-nn17u2(n-i))
	1  -nn3(n-i)*(nn23u2(n-i)-bf(n-i)*(-nn7u2(n-i))*(-nn17u2(n-i)))
      dd4(n-i)=fdu2(n-i)*(-nn8u2(n-i))*(-nn18u2(n-i))
	1  -nn5(n-i)*(nn23u2(n-i)-bf(n-i)*(-nn8u2(n-i))*(-nn18u2(n-i)))
      fuu1(n-i)=fuu1(n-i+1)*mm9u1(n-i)
	1  +nn22u1(n-i)*nn30u1(n-i)+nn24u1(n-i)*(nn27u1(n-i)*aa2(n-i)
     1  +nn28u1(n-i)*bb2(n-i)+nn25u1(n-i)*cc2(n-i)+nn26u1(n-i)*dd2(n-i))
      fuu2(n-i)=fuu2(n-i+1)*mm9u2(n-i)
	1  +nn22u2(n-i)*nn30u2(n-i)+nn24u2(n-i)*(nn27u2(n-i)*aa4(n-i)
     1  +nn28u2(n-i)*bb4(n-i)+nn25u2(n-i)*cc4(n-i)+nn26u2(n-i)*dd4(n-i))
      end do 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do i=1,n
      fu1(i)=u1*fuu1(i)+u2*fuu2(i)
      fd1(i)=u1*fdu1(i)+u2*fdu2(i)
      end do
      do i=1,n
      flxu(i)=fu1(i)
      flxd(i)=fd1(i)
      end do
      end