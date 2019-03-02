     module dimen
        implicit none
        integer                 ::nmp!分区数
        integer,allocatable     ::mm(:,:),me(:,:)!mm为网格点数，me为网格单元数
    end module
    
    module grid
        implicit none
        real*8,allocatable      ::fx(:,:,:),fz(:,:,:),vol(:,:,:),ari(:,:,:,:),ark(:,:,:,:)
                                !fx中第一个位置为分区位置，后两个为坐标点，fx表示坐标点
                                !位置x方向的值。vol表示单元格的体积，ari为i方向上的边矢量
                                !在x方向的分量，ark为z方向的分量。
    end module
    
     module qinf
        implicit none
        real*8                  :: amif,pif,dif,uif,wif,eif,alf,amach
                                !amif无量纲自由流声速，pif无量纲压力，dif密度
                                !uif，wif无量纲x和z方向速度，eif自由流总能量
                                !alf自由流攻角（弧度），amach马赫数
    end module

    
    module data0
        implicit none
        real*8                  :: gma,cfl,epsl!gma气体常数，cfl：CFL常数，epsl含义未确定
    end module
        
    
    module kswitch    
        implicit none
        integer                 :: kns,kSA!kns层流控制常数,湍流控制常数
    end module
        
    
    module file
        implicit none
        character               :: fn1*30,fn2*30,fn3*30,fn4*30,fn5*30,fn6*30!程序中文件名的变量
    end module
    
    
    module chara
        implicit none
        integer,allocatable     :: kb(:,:,:),kb1(:,:,:),kb2(:,:,:),kb3(:,:,:),md(:,:,:,:,:)
                                !kb边界点的类型，kb1对接区区号，kb2对接边序号
                                !kb3对接边对应单元序号
                                !md提供对接单元编号信息
    end module
    
    
    module conver
        implicit none
        real*8,allocatable      :: wp1(:,:,:),wp2(:,:,:),wp3(:,:,:),wp4(:,:,:),pre(:,:,:)
                                !wp为流量通量的值，wp1为密度，wp2为x方向质量流量
                                !wp3为k方向的质量流量，wp4为能量流量，pre为压力值
    end module
    
        
    module resi
        implicit none
        real*8,allocatable      :: wup1(:,:,:),wup2(:,:,:),wup3(:,:,:),wup4(:,:,:)!wup经过残差处理过的守恒变量
    end module
    
    
    module dtlo
        implicit none
        real*8,allocatable      :: dt(:,:,:)!dt当地时间步长
    end module
        
    module cfdiss
        implicit none
        real*8                  :: ek2,ek4
        real*8,allocatable      :: diss(:,:,:,:),diss0(:,:,:,:)    
    end module
        
   

    module turb
        implicit none
        real*8                  :: ren
                                !ren为雷诺数
    end module
            
        
    module rungekuta        !龙格－库塔法（Runge-Kutta methods）是用于非线性常微分方程的解的重要的一类隐式或显式迭代法。
        implicit none
        real*8                  :: rgkt(5),betaa(5)
    end module
        
    module geom_acv
        implicit none
        real*8,allocatable      :: acvi_vol(:,:,:),acvk_vol(:,:,:),acvi_ari(:,:,:,:),&
                                   acvi_ark(:,:,:,:),acvk_ari(:,:,:,:),acvk_ark(:,:,:,:)
                                    !辅助单元的信息，i,k方向的辅助单元的体积
                                    !acvi_ari为i向辅助单元i向边处边矢量
                                    !acvi_ark为i向辅助单元k向边处边矢量
                                    !acvk_ari为k向辅助单元i向边处边矢量
                                    !acvk_ark为k向辅助单元k向边处边矢量
    end module

    
!----------------------------------------------------------------------------------------
    ! TODO:主程序开始
!----------------------------------------------------------------------------------------
    
    
        program ns2d        !Program 表示主程序。但可以省略。
                            !对于不包含在 module 或者 subroutine，function 中的代码段，自动认为是 Program 主程序。

        use dimen
        use grid
        use qinf
		use data0
        use conver
        use resi
        use file
		use turb
		use rungekuta
		use cfdiss
		use kswitch
		use chara
		use dtlo
		use vissflux
		use viscous
		use safun
		use geom_acv
        
        implicit none
!------------------------------------------------------------------------------------------------------                
        integer                :: i,k,j0,irg,iread,icycle,istart,icymax,numb,ix,iz,nx,nz,nx0,ii,nsec 
                                !i表示x方向上的网格数，k表示z方向上的网格数
	                            !irg龙格库塔算法迭代次数，iread
	                            !icycle循环计算的次数，istart循环计算开始时候的值
	                            !icymax迭代最大次数，numb分区数，ix为x方向单元数
	                            !iz为z方向的单元数，nx为x方向上的网格点数，nz为k方向上的网格点数
	                            !nx0为两个方向上节点数的最大值，ii循环指示变量，nesc网格数
!------------------------------------------------------------------------------------------------------        
        real*8                 :: pi,pii,dsmt,dsmt0,dden,dsum,w1,w2,w3,w4,d,u,w,p,res,reso,btime,etime
                                !pi为π对应的弧度，pii为每度对应的弧度，
	                            !dsmt,dmst0,dden,dsum这四项与残差有关具体作用不详
                                !w1,w2,w3,w4为守恒变量，d,u,w,p为原始变量。res残差值
	                            !reso与残差有关的参数，btime表示cpu开始计算时已经运行的时间
	                            !etime表示程序运行结束后cpu运行的时间
!------------------------------------------------------------------------------------------------------        
        real*8,allocatable     :: sw1(:,:,:),sw2(:,:,:),sw3(:,:,:),sw4(:,:,:)
                                !sw为龙格库塔算法中的守恒变量
!------------------------------------------------------------------------------------------------------        

        parameter (numb=1,ix=369,iz=65)  !numb取分区数,ix取最大分区x方向单元数，iz取最大分区z方向单元数   
                                         !是否改成parameter ::numb=1,ix=369,iz=65
!------------------------------------------------------------------------------------------------------                
		call cpu_time(btime)             !cpu_time返回的是一个浮点型，表示当前CPU运行时间（以秒计）
!------------------------------------------------------------------------------------------------------                
        ! TODO:输出文件名称，以及常数设置
          fn2="rae2822.grd"              
	      fn3="temp.tmp"
	      fn4="resid.dat"
	      fn5="scp.dat"
	      fn6="flow.plt"
	      amach=0.729
	      alf=2.31
	      cfl=2.0
	      epsl=1.0e-12	
	      ren=2.133e7 
	      ek2=0.5
	      ek4=0.015	
	      gma=1.4
	      iread=0
	      icymax=200
	      kns=1
!------------------------------------------------------------------------------------------------------                
    !TODO:分配数组维数
	    allocate(mm(numb,2),me(numb,2))
	    allocate(fx(numb,ix+1,iz+1),fz(numb,ix+1,iz+1),vol(numb,ix,iz),ari(numb,ix+1,iz+1,2),ark(numb,ix+1,iz+1,2)) 
	    allocate(kb(numb,4,ix+1),kb1(numb,4,ix+1),kb2(numb,4,ix+1),kb3(numb,4,ix+1),md(numb,ix,4,2,4))
	    allocate(wp1(numb,0:ix+1,0:iz+1),wp2(numb,0:ix+1,0:iz+1),wp3(numb,0:ix+1,0:iz+1),wp4(numb,0:ix+1,0:iz+1),pre(numb,ix,iz))
	    allocate(wup1(numb,ix,iz),wup2(numb,ix,iz),wup3(numb,ix,iz),wup4(numb,ix,iz))
	    allocate(dt(numb,ix,iz))
	    allocate(diss(numb,ix,iz,4),diss0(numb,ix,iz,4))   
    	allocate(sw1(numb,ix,iz),sw2(numb,ix,iz),sw3(numb,ix,iz),sw4(numb,ix,iz))         
    	allocate(acvi_vol(numb,ix+1,iz),acvk_vol(numb,ix,iz+1),acvi_ari(numb,ix+2,iz,2),&
    	         acvi_ark(numb,ix+1,iz+1,2),acvk_ari(numb,ix+1,iz+1,2),acvk_ark(numb,ix,iz+2,2))           	         
!------------------------------------------------------------------------------------------------------                
    !	TODO:to read data from the grid file
	    open (2,file=fn2,form="binary")
	    read (2) nsec,nmp
	    write(*,*) "nsec=",nsec,",   nmp=",nmp  !nsec:网格; nmp:分区数
	    do ii=1,nmp
		    read(2) mm(ii,1),mm(ii,2)  !ii区网格沿x、z方向点个数
		    nx=mm(ii,1)  
		    nz=mm(ii,2)
		    me(ii,1)=nx-1   !ii区网格沿x、z方向单元数
		    me(ii,2)=nz-1
		    write(*,*)ii,me(ii,1),me(ii,2)
		    nx0=max(nx,nz)
		    read(2) ((fx(ii,i,k),i=1,nx),k=1,nz)   !ii区第（i，k）个网格点x方向坐标
		    read(2) ((fz(ii,i,k),i=1,nx),k=1,nz)   !ii区第（i，k）个网格点z方向坐标
		    read(2) ((kb(ii,k,i),k=1,4),i=1,nx0)   !boundary point     =0 far field  =1 wing body  =2 gene 4表示二维分区四个边界
		    read(2) ((kb1(ii,k,i),k=1,4),i=1,nx0)  !对接区区号             @有一定的内存空间浪费
		    read(2) ((kb2(ii,k,i),k=1,4),i=1,nx0)  !对接边序号             @有一定的内存空间浪费
		    read(2) ((kb3(ii,k,i),k=1,4),i=1,nx0)  !对接边对应单元序号     @有一定的内存空间浪费
	    enddo
	    close(2)
!------------------------------------------------------------------------------------------------------                
    !TODO:to give initial parameter
        pi=atan(1.0d0)*4.0
        pii=pi/180.
        alf=alf*pii !自由流攻角（弧度）
        pif=1.0d0   !自由流无量纲压强
        dif=1.0d0   !自由流无量纲密度
        amif=sqrt(gma*pif/dif) !无量纲自由流声速
        uif=amach*amif*cos(alf)  !无量纲自由流x方向速度
        wif=amach*amif*sin(alf)  !无量纲自由流z方向速度
        eif=pif/(gma-1.)+.5*dif*(amif*amach)**2 !无量纲自由流总能DENSITY*E
 	    cc=0.432d0         

!------------------------------------------------------------------------------------------------------                
	    !TODO:如果iread为0，则流场从最开始计算
	      istart=0
	      do ii=1,nmp
          do k=0,mm(ii,2)
	      do i=0,mm(ii,1)
		    wp1(ii,i,k)=dif  !未知矢量分量DENSITY
		    wp2(ii,i,k)=dif*uif  !未知矢量分量 DENSITY*U
		    wp3(ii,i,k)=dif*wif  !未知矢量分量 DENSITY*V
		    wp4(ii,i,k)=eif      !未知矢量分量 DENSITY*E
	      enddo
	      enddo
	      enddo

	    betaa(1)=1.
	    betaa(2)=0.
 	    betaa(3)=0.56
 	    betaa(4)=0.
	    betaa(5)=0.44
        rgkt(1)=1./4.
        rgkt(2)=1./6.
        rgkt(3)=3./8.
        rgkt(4)=1./2.
        rgkt(5)=1.	  

        call gmpr  !该程序计算网格单元面积，边矢量以及每个区边界处的单位法向量
        call auxi_cont_vol  !该程序计算网格辅助单元信息
	    call gdinf !整理分区对接信息
    	   
        call SAdistance !计算SA紊流模型中需要的网格单元到物面的最近距离
    	  	  
	    write(*,*)"......now the loop is starting......"
	    
        do 1000 icycle=istart+1,icymax

        if(icycle.eq.1) then
	      open(4,file=fn4)
	      write(4,61) amach,ren,alf/pii
	      close(4)
	    endif
61	    format(1x,'mach=',f6.2,3x,'Re=',E12.2,3x,'alpha=',f6.2)
    
		do ii=1,nmp
		do k=1,me(ii,2)
		do i=1,me(ii,1)
			sw1(ii,i,k)=wp1(ii,i,k)  !上一步中的守恒变量值，用于R-K推进
			sw2(ii,i,k)=wp2(ii,i,k)
			sw3(ii,i,k)=wp3(ii,i,k)
			sw4(ii,i,k)=wp4(ii,i,k)
		enddo
		enddo
 		enddo
 		

		if(kns.eq.1) then		
          call lam_vis   !计算层流粘性系数
		  if(kSA.eq.1) call SA_turb_vis
		  call derivatives		  
		  call vissi     !计算总粘性系数和粘性通量项
		endif  	 

		call clts

	    if(kns.eq.1)  then
		  if(kSA.eq.1) then
		    call vortex
		    call SAmodel
		  endif
	    endif


		do 900 irg=1,5
		
         if(irg.eq.3.or.irg.eq.5) then       !用来计算R-K推进时人工粘性项所需该步压力
			do ii=1,nmp
			do k=1,me(ii,2)
			do i=1,me(ii,1)
				w1=wp1(ii,i,k)
				w2=wp2(ii,i,k)
				w3=wp3(ii,i,k)
				w4=wp4(ii,i,k)
				call trans(w1,w2,w3,w4,d,u,w,p)
				pre(ii,i,k)=p
			enddo
			enddo
			enddo
          endif
			
		if(irg.eq.1.or.irg.eq.3.or.irg.eq.5) then
		  call dissi		!计算中心格式有限体积法的人工耗散项 
		  do ii=1,nmp
		  do k=1,me(ii,2)
		  do i=1,me(ii,1)
		  do j0=1,4
			diss(ii,i,k,j0)=betaa(irg)*diss(ii,i,k,j0)&
    		              +(1-betaa(irg))*diss0(ii,i,k,j0)
			diss0(ii,i,k,j0)=diss(ii,i,k,j0)
          enddo
          enddo
          enddo
          enddo
		endif

		call resid		  !计算对流通量项，进而求出欧拉/N-S方程有限体积法的残值。

	    do ii=1,nmp
          do k=1,me(ii,2)
          do i=1,me(ii,1)
     	    wup1(ii,i,k)=rgkt(irg)*dt(ii,i,k)*wup1(ii,i,k)/vol(ii,i,k)
		    wup2(ii,i,k)=rgkt(irg)*dt(ii,i,k)*wup2(ii,i,k)/vol(ii,i,k)
		    wup3(ii,i,k)=rgkt(irg)*dt(ii,i,k)*wup3(ii,i,k)/vol(ii,i,k)
     		wup4(ii,i,k)=rgkt(irg)*dt(ii,i,k)*wup4(ii,i,k)/vol(ii,i,k)
     	  enddo
     	  enddo
        enddo
		call ravimplicit(2d0)

		do ii=1,nmp
		do k=1,me(ii,2)
		do i=1,me(ii,1)
			wp1(ii,i,k)=sw1(ii,i,k)+wup1(ii,i,k)
			wp2(ii,i,k)=sw2(ii,i,k)+wup2(ii,i,k)
			wp3(ii,i,k)=sw3(ii,i,k)+wup3(ii,i,k)
			wp4(ii,i,k)=sw4(ii,i,k)+wup4(ii,i,k)
        enddo
        enddo
        enddo
		call farcor		
900		continue

!		call enthal(0.1d0) !该处传参应为当地时间步长@？？？？？？？？？

!		to put out the residual
     		
		dsmt=0.0
		dsum=0.0
		do ii=1,nmp
		  do k=1,me(ii,2)
		  do i=1,me(ii,1)
!		    dden=dabs(wp1(ii,i,k)-sw1(ii,i,k))/dt(ii,i,k)
			dden=dabs(wp1(ii,i,k)-sw1(ii,i,k))
			dsmt=dsmt+dden**2
          enddo
          enddo
		  dsum=dsum+me(ii,1)*me(ii,2)
        enddo

		if(icycle.eq.1) dsmt0=dsmt
		res=log(sqrt(dsmt/dsmt0))
		reso=sqrt(dsmt/dsum)
		write(*,*)icycle,res,dsmt
		if((dsmt.le.epsl).or.(icycle.gt.icymax))goto 1100

!		to save the conservation variables every 100 time steps
		 
		if(icycle/100*100.eq.icycle) then
		write(*,*)"...Now saving the temporary variables..."
		open (3,file=fn3,form="binary")
		write(3)icycle,amach,alf,dsmt0
	  	do ii=1,nmp
		write(3)((wp1(ii,i,k),i=0,mm(ii,1)),k=0,mm(ii,2))
		write(3)((wp2(ii,i,k),i=0,mm(ii,1)),k=0,mm(ii,2))
		write(3)((wp3(ii,i,k),i=0,mm(ii,1)),k=0,mm(ii,2))
		write(3)((wp4(ii,i,k),i=0,mm(ii,1)),k=0,mm(ii,2))
		write(3)((rmub(ii,i,k),i=0,mm(ii,1)),k=0,mm(ii,2))
		enddo
		call output(icycle,dsmt)
		close(3)
		endif
1000	continue
1100	continue
	    call output(icycle,dsmt)
		call cpu_time(etime)
		open(10,access="append",file=fn4)
		write(10,*) etime-btime
		close(10)			
	    write(*,*)"the loop finished"
	    end program




        
       
!------------------------------------------------------------------------------------------------------  	     
	    !辅助单元的几何信息
        subroutine ACV_gmpr
        use dimen
        use geom_acv
        implicit none
        integer                 ::ii,i,k,im,km,numb,ix,iz
        parameter(numb=1,ix=369,iz=65)
        real*8                  ::p(4,2),x21,y21,x14,y14,x34,y34,x23,y23,x31,y31,x24,y24
        real*8                  ::xacvi(numb,ix+2,iz+1),zacvi(numb,ix+2,iz+1),xacvk(numb,ix+1,iz+2),zacvk(numb,ix+1,iz+2)
        real*8                  ::mmacvi(numb,2),mmacvk(numb,2)
        common /dimen_acv/mmacvi,mmacvk
        common /grid_acv/xacvi,zacvi,xacvk,zacvk

!       I direction
     
        do ii=1,nmp
        im=mmacvi(ii,1)-1
        km=mmacvi(ii,2)-1
        do k=1,km
        do i=1,im
        p(1,1)=xacvi(ii,i,k)
        p(1,2)=zacvi(ii,i,k)
        p(2,1)=xacvi(ii,i,k+1)
        p(2,2)=zacvi(ii,i,k+1)
        p(3,1)=xacvi(ii,i+1,k+1)
        p(3,2)=zacvi(ii,i+1,k+1)
        p(4,1)=xacvi(ii,i+1,k)
        p(4,2)=zacvi(ii,i+1,k)
	    x21=p(2,1)-p(1,1)
	    y21=p(2,2)-p(1,2)
	    x14=p(1,1)-p(4,1)
	    y14=p(1,2)-p(4,2)      

	    if(i.eq.im) then
		x34=p(3,1)-p(4,1)
		y34=p(3,2)-p(4,2)
        end if     
	    if(k.eq.km) then
		x23=p(2,1)-p(3,1)
		y23=p(2,2)-p(3,2)
        end if  
 
        acvi_ari(ii,i,k,1)=y21  !i向辅助单元i向边处边矢量x方向分量
        acvi_ari(ii,i,k,2)=-x21 !i向辅助单元i向边处边矢量z方向分量
        acvi_ark(ii,i,k,1)=y14  !i向辅助单元k向边处边矢量x方向分量
        acvi_ark(ii,i,k,2)=-x14 !i向辅助单元k向边处边矢量z方向分量 
	    if(i.eq.im) then
          acvi_ari(ii,i+1,k,1)=y34
          acvi_ari(ii,i+1,k,2)=-x34
        end if
        if(k.eq.km) then
          acvi_ark(ii,i,k+1,1)=y23
          acvi_ark(ii,i,k+1,2)=-x23
        end if
      
!       The area of cell
        x31=p(3,1)-p(1,1)
        y31=p(3,2)-p(1,2)
        x24=p(2,1)-p(4,1)
        y24=p(2,2)-p(4,2)
        acvi_vol(ii,i,k)=0.5*(x31*y24-y31*x24)
        enddo
        enddo
        enddo
      
      

!     K direction
     
        do ii=1,nmp
        im=mmacvk(ii,1)-1
        km=mmacvk(ii,2)-1
        do i=1,im
        do k=1,km
        p(1,1)=xacvk(ii,i,k)
        p(1,2)=zacvk(ii,i,k)
        p(2,1)=xacvk(ii,i,k+1)
        p(2,2)=zacvk(ii,i,k+1)
        p(3,1)=xacvk(ii,i+1,k+1)
        p(3,2)=zacvk(ii,i+1,k+1)
        p(4,1)=xacvk(ii,i+1,k)
        p(4,2)=zacvk(ii,i+1,k)
	    x21=p(2,1)-p(1,1)
	    y21=p(2,2)-p(1,2)
	    x14=p(1,1)-p(4,1)
	    y14=p(1,2)-p(4,2)      

	    if(i.eq.im) then
		x34=p(3,1)-p(4,1)
		y34=p(3,2)-p(4,2)
        end if     
	    if(k.eq.km) then
		x23=p(2,1)-p(3,1)
		y23=p(2,2)-p(3,2)
        end if  
 
        acvk_ari(ii,i,k,1)=y21  !k向辅助单元i向边处边矢量x方向分量
        acvk_ari(ii,i,k,2)=-x21 !k向辅助单元i向边处边矢量z方向分量
        acvk_ark(ii,i,k,1)=y14  !k向辅助单元k向边处边矢量x方向分量
        acvk_ark(ii,i,k,2)=-x14 !k向辅助单元k向边处边矢量z方向分量
        
	    if(i.eq.im) then
          acvk_ari(ii,i+1,k,1)=y34
          acvk_ari(ii,i+1,k,2)=-x34
        end if
        if(k.eq.km) then
          acvk_ark(ii,i,k+1,1)=y23
          acvk_ark(ii,i,k+1,2)=-x23
        end if
      
!       The area of cell
        x31=p(3,1)-p(1,1)
        y31=p(3,2)-p(1,2)
        x24=p(2,1)-p(4,1)
        y24=p(2,2)-p(4,2)
        acvk_vol(ii,i,k)=0.5*(x31*y24-y31*x24)
        enddo
        enddo
        enddo
        return
        end
        
!------------------------------------------------------------------------------------------------------          
        subroutine gdinf
        use dimen
        use chara
     	implicit none
        integer                 ::izone,i,k,ip,iv,mbd(4,2)
        common /mmbbdd/mbd

	    do 100 izone=1,nmp

!	i=1 and i=imax
		
		do k=1,me(izone,2)
			if(kb(izone,3,k).eq.2)then  !i=1对接
				call patch(izone,3,k)  !得到对接边界上每个单元对接信息，包括对接区区号，对接边在对接区代号，对接单元及内侧单元在对接区序号
				do ip=1,2
				do iv=1,4
					md(izone,k,iv,ip,3)=mbd(iv,ip)
				enddo
				enddo
			endif

			if(kb(izone,4,k).eq.2)then !i=imax对接
				call patch(izone,4,k)
				do ip=1,2
				do iv=1,4
					md(izone,k,iv,ip,4)=mbd(iv,ip)
				enddo
				enddo
			endif
        enddo

!	k=1 and k=kmax

		do i=1,me(izone,1)
			if(kb(izone,1,i).eq.2)then  !k=1对接
				call patch(izone,1,i)
				do ip=1,2
				do iv=1,4
					md(izone,i,iv,ip,1)=mbd(iv,ip)
				enddo
				enddo
			endif

			if(kb(izone,2,i).eq.2)then  !i=imax对接
				call patch(izone,2,i)
				do ip=1,2
				do iv=1,4
					md(izone,i,iv,ip,2)=mbd(iv,ip)
				enddo
				enddo
			endif
        enddo
100	    continue
	    return
	    end
	    
	    
	    subroutine patch(izone,ior,ie0)   !ior为单元边界代号，ie0为该对接边界上对应单元的序号
	    use dimen
	    use chara
	    implicit none
        integer                 ::izone,ior,ie0,iqh,ibh,ieh,lp,mbd(4,2)
        common /mmbbdd/mbd

	    iqh=kb1(izone,ior,ie0)  !对接区的区号
	    ibh=kb2(izone,ior,ie0)  !对接边的代号
	    ieh=kb3(izone,ior,ie0)  !对接边上对应单元在对接区边界上的序号

	    mbd(1,1)=iqh    !mbd(i,k):i=1表示所属分区编号 i=2 沿i方向单元序号 i=3 沿k方向单元序号 i=4 所属对接边代号 k=1对接单元 k=2 对接单元内侧单元
	    mbd(1,2)=iqh
    	 
	    if(ibh.eq.1)then
		    lp=1
		    mbd(2,1)=ieh
		    mbd(3,1)=lp
		    mbd(4,1)=ibh
		    mbd(2,2)=ieh
		    mbd(3,2)=lp+1
		    mbd(4,2)=ibh
	    else if(ibh.eq.2) then
		    lp=me(iqh,2)  !k向单元数
		    mbd(2,1)=ieh
		    mbd(3,1)=lp
		    mbd(4,1)=ibh
		    mbd(2,2)=ieh
		    mbd(3,2)=lp-1
		    mbd(4,2)=ibh
	    else if(ibh.eq.3) then
		    lp=1
		    mbd(2,1)=lp
		    mbd(3,1)=ieh
		    mbd(4,1)=ibh
		    mbd(2,2)=lp+1
		    mbd(3,2)=ieh
		    mbd(4,2)=ibh
	    else if(ibh.eq.4) then
		    lp=me(iqh,1)  !i向单元数
		    mbd(2,1)=lp
		    mbd(3,1)=ieh
		    mbd(4,1)=ibh
		    mbd(2,2)=lp-1
		    mbd(3,2)=ieh
		    mbd(4,2)=ibh
	    endif
	    return
	    end
	    
	    
	    subroutine SAdistance
	    use dimen
	    use grid
		use chara
	    implicit none
        integer                 ::ii,i,k,ip,numb,ix,iz,ibp,im,km
        parameter (numb=1,ix=369,iz=65)
        real*8                  ::x0g(numb,ix,iz),z0g(numb,ix,iz),xbp(1000),zbp(1000),xc,zc,dis,dissa(numb,ix,iz),&
								  x1,x2,x3,x4,z1,z2,z3,z4
        common /SAdis/dissa
        common /xz0g/x0g,z0g

        do ii=1,nmp
	    do k=1,me(ii,2)
	    do i=1,me(ii,1)
	    x1=fx(ii,i,k)
	    x2=fx(ii,i+1,k)
	    x3=fx(ii,i,k+1)
	    x4=fx(ii,i+1,k+1)
		
	    z1=fz(ii,i,k)
	    z2=fz(ii,i+1,k)
	    z3=fz(ii,i,k+1)
	    z4=fz(ii,i+1,k+1)

	    x0g(ii,i,k)=(x1+x2+x3+x4)/4.
	    z0g(ii,i,k)=(z1+z2+z3+z4)/4.
        enddo	
        enddo
        enddo
            
	    ibp=0  !物面网格点数计数变量

	    do 50 ii=1,nmp
		    im=me(ii,1)
		    km=me(ii,2)

    !得出分区内物面边界上所有网格单元边界中心坐标和物面边界端点坐标（有可能该边界上只有部分为物面边界）
		    do 10 k=1,km
    
   !		i=1
    
			    if(k.eq.1.and.kb(ii,3,1).eq.1) then 
				    ibp=ibp+1
				    xbp(ibp)=fx(ii,1,1)
				    zbp(ibp)=fz(ii,1,1)
			    endif
    								
			    if(k.eq.km.and.kb(ii,3,km).eq.1) then
				    ibp=ibp+1
				    xbp(ibp)=fx(ii,1,km+1)
				    zbp(ibp)=fz(ii,1,km+1)
			    endif

			    if(k.gt.1.and.k.lt.km) then
				    if(kb(ii,3,k).eq.1.and.kb(ii,3,k-1).ne.1)then
					    ibp=ibp+1
					    xbp(ibp)=fx(ii,1,k)
					    zbp(ibp)=fz(ii,1,k)
				    endif

				    if(kb(ii,3,k).eq.1.and.kb(ii,3,k+1).ne.1)then
					    ibp=ibp+1
					    xbp(ibp)=fx(ii,1,k+1)
					    zbp(ibp)=fz(ii,1,k+1)
				    endif
			    endif

			    if(kb(ii,3,k).eq.1)then
				    ibp=ibp+1
				    xbp(ibp)=.5*(fx(ii,1,k)+fx(ii,1,k+1))
				    zbp(ibp)=.5*(fz(ii,1,k)+fz(ii,1,k+1))
			    endif

    
    !		i=imax
    
			    if(k.eq.1.and.kb(ii,4,1).eq.1) then
				    ibp=ibp+1
				    xbp(ibp)=fx(ii,im+1,1)
				    zbp(ibp)=fz(ii,im+1,1)
			    endif
    								
			    if(k.eq.km.and.kb(ii,4,km).eq.1) then
				    ibp=ibp+1
				    xbp(ibp)=fx(ii,im+1,km+1)
				    zbp(ibp)=fz(ii,im+1,km+1)
			    endif

			    if(k.gt.1.and.k.lt.km) then
				    if(kb(ii,4,k).eq.1.and.kb(ii,4,k-1).ne.1)then
					    ibp=ibp+1
					    xbp(ibp)=fx(ii,im+1,k)
					    zbp(ibp)=fz(ii,im+1,k)
				    endif

				    if(kb(ii,4,k).eq.1.and.kb(ii,4,k+1).ne.1)then
					    ibp=ibp+1
					    xbp(ibp)=fx(ii,im+1,k+1)
					    zbp(ibp)=fz(ii,im+1,k+1)
				    endif
			    endif

			    if(kb(ii,4,k).eq.1)then
				    ibp=ibp+1
				    xbp(ibp)=.5*(fx(ii,im+1,k)+fx(ii,im+1,k+1))
				    zbp(ibp)=.5*(fz(ii,im+1,k)+fz(ii,im+1,k+1))
			    endif
    10		continue

		    do 20 i=1,im
    
    !		k=1
    
			    if(i.eq.1.and.kb(ii,1,1).eq.1)then
				    ibp=ibp+1
				    xbp(ibp)=fx(ii,1,1)
				    zbp(ibp)=fz(ii,1,1)
			    endif
    					
			    if(i.eq.im.and.kb(ii,1,im).eq.1)then
				    ibp=ibp+1
				    xbp(ibp)=fx(ii,im+1,1)
				    zbp(ibp)=fz(ii,im+1,1)
			    endif

			    if(i.gt.1.and.i.lt.im)then
				    if(kb(ii,1,i).eq.1.and.kb(ii,1,i-1).ne.1)then
					    ibp=ibp+1
					    xbp(ibp)=fx(ii,i,1)
					    zbp(ibp)=fz(ii,i,1)
				    endif
				    if(kb(ii,1,i).eq.1.and.kb(ii,1,i+1).ne.1)then
					    ibp=ibp+1
					    xbp(ibp)=fx(ii,i+1,1)
					    zbp(ibp)=fz(ii,i+1,1)
				    endif
			    endif			

			    if(kb(ii,1,i).eq.1)then
				    ibp=ibp+1
				    xbp(ibp)=.5*(fx(ii,i,1)+fx(ii,i+1,1))
				    zbp(ibp)=.5*(fz(ii,i,1)+fz(ii,i+1,1))
			    endif
    
    !		k=kmax
    
    			
			    if(i.eq.1.and.kb(ii,2,1).eq.1)then
				    ibp=ibp+1
				    xbp(ibp)=fx(ii,1,km+1)
				    zbp(ibp)=fz(ii,1,km+1)
			    endif
    					
			    if(i.eq.im.and.kb(ii,2,im).eq.1)then
				    ibp=ibp+1
				    xbp(ibp)=fx(ii,im+1,km+1)
				    zbp(ibp)=fz(ii,im+1,km+1)
			    endif

			    if(i.gt.1.and.i.lt.im)then
				    if(kb(ii,2,i).eq.1.and.kb(ii,2,i-1).ne.1)then
					    ibp=ibp+1
					    xbp(ibp)=fx(ii,i,km+1)
					    zbp(ibp)=fz(ii,i,km+1)
				    endif
				    if(kb(ii,2,i).eq.1.and.kb(ii,2,i+1).ne.1)then
					    ibp=ibp+1
					    xbp(ibp)=fx(ii,i+1,km+1)
					    zbp(ibp)=fz(ii,i+1,km+1)
				    endif
			    endif			

			    if(kb(ii,2,i).eq.1)then
				    ibp=ibp+1
				    xbp(ibp)=.5*(fx(ii,i,km+1)+fx(ii,i+1,km+1))
				    zbp(ibp)=.5*(fz(ii,i,km+1)+fz(ii,i+1,km+1))
			    endif
    20		continue
    50	continue
    	
	    write(*,*)"points at wall = ",ibp
    	
        do ii=1,nmp        !排序法
	    do k=1,me(ii,2)
	    do i=1,me(ii,1)
		  xc=(x0g(ii,i,k)-xbp(1))**2
		  zc=(z0g(ii,i,k)-zbp(1))**2
		  dissa(ii,i,k)=sqrt(xc+zc)
		  do ip=2,ibp
			xc=(x0g(ii,i,k)-xbp(ip))**2
			zc=(z0g(ii,i,k)-zbp(ip))**2
			dis=sqrt(xc+zc)
			if(dissa(ii,i,k).gt.dis)dissa(ii,i,k)=dis
          enddo
		  	if(dissa(ii,i,k).eq.0) then
			write(*,*) i,k
			pause
			endif
        enddo
        enddo
        enddo       
	    return	
	    end
	    
        subroutine trans(w1,w2,w3,w4,d,u,w,p)  !守恒变量转变为原始变量
        use qinf
        use data0
        implicit none
        real*8                    ::w1,w2,w3,w4,d,u,w,p,dd,q2      
        d=w1
        dd=1./d
        u=w2*dd
        w=w3*dd
        q2=u*u+w*w
        p=(gma-1.d0)*(w4-.5*q2*d)
        return
        end
	    
	      
        subroutine clts   !当地时间步长,同时求出处理人工耗散时所需的相耗散项尺度因子
        use dimen
        use grid
        use chara
        use conver
        use dtlo
        use data0
        use kswitch
        use qinf
        use viscous
        use geom_acv
        implicit none
        integer                   ::numb,ix,iz,ii,im,km,i,k
        parameter(numb=1,ix=369,iz=65)
        real*8                    ::prl,prlt,w1,w2,w3,w4,d,u,w,p
        real*8                    ::dti(numb,ix,iz),dtk(numb,ix,iz),gmki(numb,ix,iz),gmik(numb,ix,iz)
        real*8                    ::csi2,csk2,csi,csk,vdsi,vdsk,ssp,dtic,dtkc,vis,dtiv,dtkv,sigma,phii,phik
        common /dtik/dti,dtk
        common /ggmm/gmki,gmik

        prl=0.72 !层流普朗特数
        prlt=0.9 !紊流普朗特数

        do 100 ii=1,nmp
        im=me(ii,1)
        km=me(ii,2)

        do k=1,km
        do i=1,im
        w1=wp1(ii,i,k)
        w2=wp2(ii,i,k)
        w3=wp3(ii,i,k)
        w4=wp4(ii,i,k)
        call trans(w1,w2,w3,w4,d,u,w,p)
        pre(ii,i,k)=p

        if(p/d.le.0)then
          write(*,*)"error:",ii,i,k
        endif

        csi2=acvi_ari(ii,i+1,k,1)**2+acvi_ari(ii,i+1,k,2)**2
        csk2=acvk_ark(ii,i,k+1,1)**2+acvk_ark(ii,i,k+1,2)**2
        csi=sqrt(csi2)   !单元格心处i方向长度
        csk=sqrt(csk2)   !单元格心处k方向长度
        vdsi=abs(u*acvi_ari(ii,i+1,k,1)+w*acvi_ari(ii,i+1,k,2))/csi  !沿i方向单元格心处的法向速度
        vdsk=abs(u*acvk_ark(ii,i,k+1,1)+w*acvk_ark(ii,i,k+1,2))/csk  !沿k方向单元格心处的法向速度
        ssp=sqrt(gma*p/d)     !声速平方
        dtic=(vdsi+ssp)*csi   !对流通量Jacobians矩阵在I方向的谱半径
        dtkc=(vdsk+ssp)*csk   !对流通量Jacobians矩阵在k方向的谱半径
        if(kns.eq.1) then
        vis=ul(ii,i,k)/prl+ut(ii,i,k)/prlt
        dtiv=max(4/(3*d),gma/d)*vis*csi2/vol(ii,i,k)  !粘性通量Jacobians矩阵在I方向的谱半径
        dtkv=max(4/(3*d),gma/d)*vis*csk2/vol(ii,i,k)  !粘性通量Jacobians矩阵在k方向的谱半径
        else
        dtiv=0.
        dtkv=0.
        endif
        dt(ii,i,k)=cfl*vol(ii,i,k)/((dtic+dtkc)+4*(dtiv+dtkv))    !计算当地时间步长
        gmki(ii,i,k)=dtkc/dtic  !对流通量Jacobians谱半径比值，k/i
        gmik(ii,i,k)=dtic/dtkc  !对流通量Jacobians谱半径比值，i/k
        sigma=0.55			
        phii=1+(dtkc/dtic)**sigma   !处理人工耗散项时与i方向有关的系数
        phik=1+(dtic/dtkc)**sigma   !处理人工耗散项时与k方向有关的系数
        dti(ii,i,k)=dtic*phii       !与重新定义的耗散项尺度因子有关 手册P34
        dtk(ii,i,k)=dtkc*phik       !与重新定义的耗散项尺度因子有关
        enddo
        enddo
100	    continue
        return
        end
        
        
        subroutine lam_vis
        use dimen
        use data0
        use conver
        use viscous
        use turb
        use qinf
        implicit none
        integer                   ::ii,i,k
        real*8                    ::w1,w2,w3,w4,d,u,w,p,t
        
        do ii=1,nmp
        do k=1,me(ii,2)
        do i=1,me(ii,1)
            w1=wp1(ii,i,k)
            w2=wp2(ii,i,k)
            w3=wp3(ii,i,k)
            w4=wp4(ii,i,k)
            call trans(w1,w2,w3,w4,d,u,w,p)
            t=p/d
            ul(ii,i,k)=(t**(1.5)*(1+cc)/(t+cc))*sqrt(gma)*amach/ren
	    
        enddo
        enddo
        enddo
        return
        end
        
        subroutine SA_turb_vis
		use dimen
		use safun
		use conver
		use viscous
		implicit none
		integer				::i,k,ii
		real*8				::cv1,cv13,rkai3,fv1

!      	to evaluate turbulent viscosity 

	    cv1=7.1
	    cv13=cv1**3
	    do ii=1,nmp
	    do k=1,me(ii,2)
	    do i=1,me(ii,1)
          rmu(ii,i,k)=ul(ii,i,k)/wp1(ii,i,k)
	      rkai3=(rmub(ii,i,k)/rmu(ii,i,k))**3
	      fv1=rkai3/(rkai3+cv13)	
	      ut(ii,i,k)=rmub(ii,i,k)*fv1*wp1(ii,i,k)
        enddo
        enddo
        enddo

		return
		end


	    subroutine vortex   !计算涡量与变形率
	    use dimen
	    use grid
	    use qinf
	    use conver
	    use chara
	    implicit none
        integer                   ::ii,i,k,numb,ix,iz,kis1,kis2,kks1,kks2,im,km,ia,i1,i2
        parameter(numb=1,ix=369,iz=65)    
        real*8                    ::w1,w2,w3,w4,d,u,w,p,t,dx02,ux02,wx02,px02,tx02,dx2,ux2,wx2,px2,tx2,&
                                    dx01,ux01,wx01,px01,tx01,dx1,ux1,wx1,px1,tx1,dz02,uz02,wz02,pz02,tz02,&
                                    dz2,uz2,wz2,pz2,tz2,dz01,uz01,wz01,pz01,tz01,dz1,uz1,wz1,pz1,tz1,&
                                    x11,y11,x12,y12,tt11,tt12,u11,u12,w11,w12,t11,t12,a10,d10
        real*8                    ::strate(numb,ix,iz),wow(numb,ix,iz)
	    common /str_wow/strate,wow	

	    do 200 ii=1,nmp
	        kis1=0
		    kis2=0
		    kks1=0
		    kks2=0
		    im=me(ii,1)
		    km=me(ii,2)
     
          do 100 k=1,km
          do 100 i=1,im
!         当前单元
              w1=wp1(ii,i,k)
              w2=wp2(ii,i,k)
              w3=wp3(ii,i,k)
              w4=wp4(ii,i,k)
              call trans(w1,w2,w3,w4,d,u,w,p)
              t=p/d
!         直接计算边界4上的各个流动值   
	      if(i.eq.im) then  
		    if(kb(ii,4,k).eq.0) then  !远场 
			    w1=wp1(ii,im+1,k)
			    w2=wp2(ii,im+1,k)
			    w3=wp3(ii,im+1,k)
			    w4=wp4(ii,im+1,k)
		    elseif(kb(ii,4,k).eq.1) then  !物面
                ux02=0
                wx02=0
                px02=p
                dx02=d
                goto 5
		    else                          !对接边界
			    ia=md(ii,k,1,1,4)  !ii区编号为4的边界上沿k向序号k单元对接单元所属区编号
			    i1=md(ii,k,2,1,4)  !ii区编号为4的边界上沿k向序号k单元对接单元所属区沿i向序号
			    i2=md(ii,k,3,1,4)  !ii区编号为4的边界上沿k向序号k单元对接单元所属区沿k向序号
			    w1=.5*(wp1(ii,im,k)+wp1(ia,i1,i2))
			    w2=.5*(wp2(ii,im,k)+wp2(ia,i1,i2))
			    w3=.5*(wp3(ii,im,k)+wp3(ia,i1,i2))
			    w4=.5*(wp4(ii,im,k)+wp4(ia,i1,i2))
		    endif
              
		    call trans(w1,w2,w3,w4,dx02,ux02,wx02,px02)
5           continue
              tx02=px02/dx02
              kis2=1
              goto 10
          end if
!         当前单元右侧(i+1)单元
          w1=wp1(ii,i+1,k)
          w2=wp2(ii,i+1,k)
          w3=wp3(ii,i+1,k)
          w4=wp4(ii,i+1,k)
          call trans(w1,w2,w3,w4,dx2,ux2,wx2,px2)
          tx2=px2/dx2
10        continue
    !直接计算边界3上的各个流动值 
          if(i.eq.1) then
		    if(kb(ii,3,k).eq.0) then
		        w1=wp1(ii,0,k)
			    w2=wp2(ii,0,k)
			    w3=wp3(ii,0,k)
			    w4=wp4(ii,0,k)
		    elseif(kb(ii,3,k).eq.1) then
                ux01=0
                wx01=0
                px01=p
                dx01=d
                goto 15
		    else
			    ia=md(ii,k,1,1,3)
			    i1=md(ii,k,2,1,3)
			    i2=md(ii,k,3,1,3)
			    w1=.5*(wp1(ii,1,k)+wp1(ia,i1,i2))
			    w2=.5*(wp2(ii,1,k)+wp2(ia,i1,i2))
			    w3=.5*(wp3(ii,1,k)+wp3(ia,i1,i2))
			    w4=.5*(wp4(ii,1,k)+wp4(ia,i1,i2))
		    endif

              call trans(w1,w2,w3,w4,dx01,ux01,wx01,px01)
15          continue      
              tx01=px01/dx01
              kis1=1
              goto 20
          endif
!         当前单元左侧(i-1)单元
          w1=wp1(ii,i-1,k)
          w2=wp2(ii,i-1,k)
          w3=wp3(ii,i-1,k)
          w4=wp4(ii,i-1,k)
          call trans(w1,w2,w3,w4,dx1,ux1,wx1,px1)
          tx1=px1/dx1
20        continue
    !直接计算边界2上的各个流动值
          if(k.eq.km)then
		    if(kb(ii,2,i).eq.0) then
			    w1=wp1(ii,i,km+1)
			    w2=wp2(ii,i,km+1)
			    w3=wp3(ii,i,km+1)
			    w4=wp4(ii,i,km+1)
		    elseif(kb(ii,2,i).eq.1) then
                uz02=0
                wz02=0
                pz02=p
                dz02=d
                goto 25
		    else
			    ia=md(ii,i,1,1,2)
			    i1=md(ii,i,2,1,2)
			    i2=md(ii,i,3,1,2)
			    w1=.5*(wp1(ii,i,km)+wp1(ia,i1,i2))
			    w2=.5*(wp2(ii,i,km)+wp2(ia,i1,i2))
			    w3=.5*(wp3(ii,i,km)+wp3(ia,i1,i2))
			    w4=.5*(wp4(ii,i,km)+wp4(ia,i1,i2))
		    endif

              call trans(w1,w2,w3,w4,dz02,uz02,wz02,pz02)
25          continue  
              tz02=pz02/dz02
              kks2=1
              goto 30
          endif
!         当前单元上侧(k+1)单元     
	      w1=wp1(ii,i,k+1)
          w2=wp2(ii,i,k+1)
          w3=wp3(ii,i,k+1)
          w4=wp4(ii,i,k+1)
          call trans(w1,w2,w3,w4,dz2,uz2,wz2,pz2)
          tz2=pz2/dz2
30        continue
    !直接计算边界1上的各个流动值
          if(k.eq.1)then
              if(kb(ii,1,i).eq.0)then
			    w1=wp1(ii,i,0)
			    w2=wp2(ii,i,0)
			    w3=wp3(ii,i,0)
			    w4=wp4(ii,i,0)
		    elseif(kb(ii,1,i).eq.1)then
                uz01=0
                wz01=0
                pz01=p
                dz01=d
                goto 35
		    else
			    ia=md(ii,i,1,1,1)
			    i1=md(ii,i,2,1,1)
			    i2=md(ii,i,3,1,1)
			    w1=.5*(wp1(ii,i,1)+wp1(ia,i1,i2))
			    w2=.5*(wp2(ii,i,1)+wp2(ia,i1,i2))
			    w3=.5*(wp3(ii,i,1)+wp3(ia,i1,i2))
			    w4=.5*(wp4(ii,i,1)+wp4(ia,i1,i2))
		    endif
              call trans(w1,w2,w3,w4,dz01,uz01,wz01,pz01)
35          continue    
		    tz01=pz01/dz01
              kks1=1
              goto 40
          endif
!         当前单元下侧(k-1)单元
	      w1=wp1(ii,i,k-1)
          w2=wp2(ii,i,k-1)
          w3=wp3(ii,i,k-1)
          w4=wp4(ii,i,k-1)
          call trans(w1,w2,w3,w4,dz1,uz1,wz1,pz1)
          tz1=pz1/dz1
40        continue

    !对于边界单元直接根据kis、kks控制参数来提取；对于内部单元用中心格式计算边界流动值
          if(kis1.eq.1)then
             kis1=0
          else
             ux01=(ux1+u)/2.
             wx01=(wx1+w)/2.
             tx01=(tx1+t)/2.
          endif

          if(kis2.eq.1)then
             kis2=0
          else
             ux02=(ux2+u)/2.
             wx02=(wx2+w)/2.
             tx02=(tx2+t)/2.
          endif

          if(kks1.eq.1)then
            kks1=0
          else
            uz01=(uz1+u)/2.
            wz01=(wz1+w)/2.
            tz01=(tz1+t)/2.
          endif

          if(kks2.eq.1)then
            kks2=0
          else
            uz02=(uz2+u)/2.
            wz02=(wz2+w)/2.
            tz02=(tz2+t)/2.
          endif
    !计算得到各单元边界速度后根据格林公式计算单元速度偏导
            x11=ux02*ari(ii,i+1,k,1)-ux01*ari(ii,i,k,1)&
               +uz02*ark(ii,i,k+1,1)-uz01*ark(ii,i,k,1)
            x12=ux02*ari(ii,i+1,k,2)-ux01*ari(ii,i,k,2)&
               +uz02*ark(ii,i,k+1,2)-uz01*ark(ii,i,k,2)
            y11=wx02*ari(ii,i+1,k,1)-wx01*ari(ii,i,k,1)&
               +wz02*ark(ii,i,k+1,1)-wz01*ark(ii,i,k,1)
            y12=wx02*ari(ii,i+1,k,2)-wx01*ari(ii,i,k,2)&
               +wz02*ark(ii,i,k+1,2)-wz01*ark(ii,i,k,2)
            tt11=tx02*ari(ii,i+1,k,1)-tx01*ari(ii,i,k,1)&
               +tz02*ark(ii,i,k+1,1)-tz01*ark(ii,i,k,1)
            tt12=tx02*ari(ii,i+1,k,2)-tx01*ari(ii,i,k,2)&
                +tz02*ark(ii,i,k+1,2)-tz01*ark(ii,i,k,2)
    !格林公式计算速度导数
	      u11=x11/vol(ii,i,k)
	      u12=x12/vol(ii,i,k)
	      w11=y11/vol(ii,i,k)
	      w12=y12/vol(ii,i,k)
	      t11=tt11/vol(ii,i,k)
	      t12=tt12/vol(ii,i,k)
    !注意，涡量与剪切率公式均是张量表达式
	      wow(ii,i,k)=dabs(u12-w11)  !涡量
	      a10=(u12+w11)**2
	      d10=2*(u11**2+w12**2)
	      strate(ii,i,k)=sqrt(a10+d10)  !剪切率 

100	    continue
200	    continue
	    return
	    end
	    
	    

!      **********Spalart-Allmaras one-equation model********

	    subroutine SAmodel
	    use dimen
	    use data0
	    use turb
	    use dtlo
	    use conver
	    use viscous
	    use safun
	    use rungekuta
	    implicit none
        integer                   ::ii,i,k,numb,ix,iz,irg
        parameter(numb=1,ix=369,iz=65)    
        real*8                    ::sw5(numb,ix,iz),wup5(numb,ix,iz),aa,qt(numb,ix,iz),qtim1(numb,ix,iz)
        common /resa/wup5
	    common /SAsource/qt,qtim1
	                    
     	do ii=1,nmp
	    do k=1,me(ii,2)
	    do i=1,me(ii,1)
	      rmu(ii,i,k)=ul(ii,i,k)/wp1(ii,i,k)	 ! 层流运动粘性系数
	      sw5(ii,i,k)=rmub(ii,i,k)   !SA一方程紊流模型控制方程变量,上一时间步的rmub值
        enddo
        enddo
        enddo

	    call contravelo

	    do irg=1,5       !龙格-库塔五步迭代
            
            if(irg.eq.1) then
		    call saderimub
		    call SA_source_term
		    endif
		    
		    call residsa

		    do ii=1,nmp
		    do k=1,me(ii,2)
		    do i=1,me(ii,1)
		        aa=1-qtim1(ii,i,k)*dt(ii,i,k)     
			    rmub(ii,i,k)=max(sw5(ii,i,k)+rgkt(irg)*dt(ii,i,k)*wup5(ii,i,k)/aa,ru0t*1e-10)
			    if(rmub(ii,i,k).eq.ru0t*1e-20) then
			    write(*,*) "湍流粘性系数为负"
                pause
			    endif   
            enddo
            enddo
            enddo
        enddo    	
	    return
	    end
	    

	    subroutine contravelo   !计算单元边界法向速度
	    use dimen
	    use grid
	    use chara
	    use conver
	    implicit none
        integer                   ::ii,i,k,numb,ix,iz,im,km,ia,i1,i2
        parameter(numb=1,ix=369,iz=65)
        real*8                    ::w1,w2,w3,w4,d,u,w,p  
        real*8                    ::ctrvi(numb,ix+1,iz),ctrvk(numb,ix,iz+1)
	    common /contrv/ctrvi,ctrvk

	    do 500 ii=1,nmp
	      im=me(ii,1)
	      km=me(ii,2)	

!   	  i direction	
  
          do 100 k=1,me(ii,2)
	      do 100 i=1,mm(ii,1)
	        if(i.eq.1) then
	          if(kb(ii,3,k).eq.0) then
		        w1=wp1(ii,0,k)
		        w2=wp2(ii,0,k)
		        w3=wp3(ii,0,k)
		        w4=wp4(ii,0,k)
		      else if(kb(ii,3,k).eq.1) then
		        ctrvi(ii,i,k)=0.    !法向速度为0（无穿透）
		        goto 100
		      else  
		        ia=md(ii,k,1,1,3)
		        i1=md(ii,k,2,1,3)
		        i2=md(ii,k,3,1,3)
		        w1=0.5*(wp1(ii,1,k)+wp1(ia,i1,i2))
		        w2=0.5*(wp2(ii,1,k)+wp2(ia,i1,i2))
		        w3=0.5*(wp3(ii,1,k)+wp3(ia,i1,i2))
		        w4=0.5*(wp4(ii,1,k)+wp4(ia,i1,i2))
		      endif
		    elseif(i.eq.mm(ii,1)) then
		      if(kb(ii,4,k).eq.0) then
		        w1=wp1(ii,im+1,k)		
			    w2=wp2(ii,im+1,k)		
			    w3=wp3(ii,im+1,k)		
			    w4=wp4(ii,im+1,k)
		      elseif(kb(ii,4,k).eq.1) then
			    ctrvi(ii,i,k)=0.
			    goto 100
		      else
			    ia=md(ii,k,1,1,4)
			    i1=md(ii,k,2,1,4)
			    i2=md(ii,k,3,1,4)
			    w1=0.5*(wp1(ii,im,k)+wp1(ia,i1,i2))
			    w2=0.5*(wp2(ii,im,k)+wp2(ia,i1,i2))
			    w3=0.5*(wp3(ii,im,k)+wp3(ia,i1,i2))
			    w4=0.5*(wp4(ii,im,k)+wp4(ia,i1,i2))
		      endif
		    else
                w1=0.5*(wp1(ii,i,k)+wp1(ii,i-1,k))
                w2=0.5*(wp2(ii,i,k)+wp2(ii,i-1,k))
                w3=0.5*(wp3(ii,i,k)+wp3(ii,i-1,k))
                w4=0.5*(wp4(ii,i,k)+wp4(ii,i-1,k))
		    endif
    	    
		    call trans(w1,w2,w3,w4,d,u,w,p)
	        ctrvi(ii,i,k)=u*ari(ii,i,k,1)+w*ari(ii,i,k,2)  !单元边界法向速度沿边界积分（即边界法向速度乘以边界面积（长度））
100	      continue

!   	  k direction

	      do 200 i=1,me(ii,1)
	      do 200 k=1,mm(ii,2)
	        if(k.eq.1) then
	          if(kb(ii,1,i).eq.0) then
		        w1=wp1(ii,i,0)
		        w2=wp2(ii,i,0)
		        w3=wp3(ii,i,0)
		        w4=wp4(ii,i,0)
		      elseif(kb(ii,1,i).eq.1) then
		        ctrvk(ii,i,k)=0.
		        goto 200
		      elseif(kb(ii,1,i).eq.2) then
		        ia=md(ii,i,1,1,1)
		        i1=md(ii,i,2,1,1)
		        i2=md(ii,i,3,1,1)
		        w1=0.5*(wp1(ii,i,1)+wp1(ia,i1,i2))
		        w2=0.5*(wp2(ii,i,1)+wp2(ia,i1,i2))
		        w3=0.5*(wp3(ii,i,1)+wp3(ia,i1,i2))
		        w4=0.5*(wp4(ii,i,1)+wp4(ia,i1,i2))
		      endif		
		    elseif(k.eq.mm(ii,2)) then
     		    if(kb(ii,2,i).eq.0) then
		        w1=wp1(ii,i,km+1)
			    w2=wp2(ii,i,km+1)
			    w3=wp3(ii,i,km+1)
			    w4=wp4(ii,i,km+1)
		      else if(kb(ii,2,i).eq.1) then
			    ctrvk(ii,i,k)=0.
			    goto 200
		      else
		        ia=md(ii,i,1,1,2)
			    i1=md(ii,i,2,1,2)
			    i2=md(ii,i,3,1,2)
			    w1=0.5*(wp1(ii,i,km)+wp1(ia,i1,i2))
			    w2=0.5*(wp2(ii,i,km)+wp2(ia,i1,i2))
			    w3=0.5*(wp3(ii,i,km)+wp3(ia,i1,i2))
			    w4=0.5*(wp4(ii,i,km)+wp4(ia,i1,i2))
		      endif
		    else
		      w1=0.5*(wp1(ii,i,k)+wp1(ii,i,k-1))
		      w2=0.5*(wp2(ii,i,k)+wp2(ii,i,k-1))
		      w3=0.5*(wp3(ii,i,k)+wp3(ii,i,k-1))
		      w4=0.5*(wp4(ii,i,k)+wp4(ii,i,k-1))
		    endif	        
    		
		    call trans(w1,w2,w3,w4,d,u,w,p)
		    ctrvk(ii,i,k)=u*ark(ii,i,k,1)+w*ark(ii,i,k,2)
200	    continue
500     continue
	    return
	    end
    
    
    
	    subroutine saderimub      !计算rmub的导数（包括在边界上的值和在单元中心处的值）
	    use dimen
	    use chara
	    use grid
	    use geom_acv
	    use turb
	    use safun
	    implicit none
        integer                   ::ii,i,k,numb,ix,iz,im,km,ia,i1,i2
        parameter(numb=1,ix=369,iz=65)
        real*8                    ::dmux2,dmux1,dmuz21,dmuz22,dmuz2,dmuz11,dmuz12,dmuz1,&
                                    pmubpx,pmubpz,volume,dmux21,dmux22,dmux11,dmux12
        real*8                    ::dmui1(numb,ix+1,iz),dmui2(numb,ix+1,iz),dmuk1(numb,ix,iz+1),dmuk2(numb,ix,iz+1)
	    common /saderiv/dmui1,dmui2,dmuk1,dmuk2
           	
	    do 300 ii=1,nmp
		    im=me(ii,1)
		    km=me(ii,2)

    !
    !	  to supply values of conservative variables in the dummy cells
    !
    !
    !	  i=0
    ! 
		    do k=1,km
		    if(kb(ii,3,k).eq.1) then 
		       rmub(ii,0,k)=-rmub(ii,1,k)  !物面镜像
		    elseif(kb(ii,3,k).eq.2) then  !对接
		       ia=md(ii,k,1,1,3)
		       i1=md(ii,k,2,1,3)
		       i2=md(ii,k,3,1,3)
		       rmub(ii,0,k)=rmub(ia,i1,i2)
		    endif
		    enddo
    		
    !
    !	  i=im+1
    ! 

		    do k=1,km
		    if(kb(ii,4,k).eq.1) then 
		       rmub(ii,im+1,k)=-rmub(ii,im,k)  !物面镜像
		    elseif(kb(ii,4,k).eq.2) then  !对接
		       ia=md(ii,k,1,1,4)
		       i1=md(ii,k,2,1,4)
		       i2=md(ii,k,3,1,4)
		       rmub(ii,im+1,k)=rmub(ia,i1,i2)
		    endif
		    enddo
    	
    !
    !	  k=0
    ! 
		    do i=1,im
		    if(kb(ii,1,i).eq.1) then 
		       rmub(ii,i,0)=-rmub(ii,i,1)  !物面镜像
		    elseif(kb(ii,1,i).eq.2) then  !对接
		       ia=md(ii,i,1,1,1)
		       i1=md(ii,i,2,1,1)
		       i2=md(ii,i,3,1,1)
		       rmub(ii,i,0)=rmub(ia,i1,i2)
		    endif
		    enddo

    !
    !	  k=km+1
    ! 

		    do i=1,im
		    if(kb(ii,2,i).eq.1) then 
		       rmub(ii,i,km+1)=-rmub(ii,i,km)  !物面镜像
		    elseif(kb(ii,2,i).eq.2) then  !对接
		       ia=md(ii,i,1,1,2)
		       i1=md(ii,i,2,1,2)
		       i2=md(ii,i,3,1,2)
		       rmub(ii,i,km+1)=rmub(ia,i1,i2)
		    endif
		    enddo

    !
    !	  四个角上的虚拟单元值
    ! 
            rmub(ii,0,0)=.5*(rmub(ii,0,1)+rmub(ii,1,0))
            rmub(ii,0,km+1)=.5*(rmub(ii,0,km)+rmub(ii,1,km+1))
            rmub(ii,im+1,0)=.5*(rmub(ii,im,0)+rmub(ii,im+1,1))
            rmub(ii,im+1,km+1)=.5*(rmub(ii,im,km+1)+rmub(ii,im+1,km))


    !
    !	  to evaluate the first derivatives of the turbulent variable mu in cell edges
    !

    !
    !	  I direction
    !

            do k=1,me(ii,2)
            do i=1,mm(ii,1)
              dmux2=rmub(ii,i,k)
              dmux1=rmub(ii,i-1,k)
              dmuz21=rmub(ii,i-1,k+1)
              dmuz22=rmub(ii,i,k+1)
              dmuz2=.25*(dmux2+dmux1+dmuz21+dmuz22)
              dmuz11=rmub(ii,i-1,k-1)
              dmuz12=rmub(ii,i,k-1)
              dmuz1=.25*(dmux2+dmux1+dmuz11+dmuz12)
              
              if(i.eq.1) then
                pmubpx=dmux2*ari(ii,i+1,k,1)-dmux1*ari(ii,i,k,1)&
                     +dmuz2*ark(ii,i,k+1,1)-dmuz1*ark(ii,i,k,1)
                pmubpz=dmux2*ari(ii,i+1,k,2)-dmux1*ari(ii,i,k,2)&
                     +dmuz2*ark(ii,i,k+1,2)-dmuz1*ark(ii,i,k,2)
               volume=vol(ii,1,k)
              elseif(i.eq.mm(ii,1)) then
                pmubpx=dmux2*ari(ii,i,k,1)-dmux1*ari(ii,i-1,k,1)&
                     +dmuz2*ark(ii,im,k+1,1)-dmuz1*ark(ii,im,k,1)
                pmubpz=dmux2*ari(ii,i,k,2)-dmux1*ari(ii,i-1,k,2)&
                     +dmuz2*ark(ii,im,k+1,2)-dmuz1*ark(ii,im,k,2)
               volume=vol(ii,im,k)
              else
		      pmubpx=dmux2*acvi_ari(ii,i+1,k,1)-dmux1*acvi_ari(ii,i,k,1)&
         		    +dmuz2*acvi_ark(ii,i,k+1,1)-dmuz1*acvi_ark(ii,i,k,1)
		      pmubpz=dmux2*acvi_ari(ii,i+1,k,2)-dmux1*acvi_ari(ii,i,k,2)&
         		    +dmuz2*acvi_ark(ii,i,k+1,2)-dmuz1*acvi_ark(ii,i,k,2)
		      volume=acvi_vol(ii,i,k)
		      endif
    		  
		      dmui1(ii,i,k)=pmubpx/volume
		      dmui2(ii,i,k)=pmubpz/volume
		    enddo
		    enddo


    !
    !	  K direction
    !

            do i=1,me(ii,1)
            do k=1,mm(ii,2)
              dmuz2=rmub(ii,i,k)
              dmuz1=rmub(ii,i,k-1)
              dmux21=rmub(ii,i+1,k-1)
              dmux22=rmub(ii,i+1,k)
              dmux2=.25*(dmuz2+dmuz1+dmux21+dmux22)
              dmux11=rmub(ii,i-1,k-1)
              dmux12=rmub(ii,i-1,k)
              dmux1=.25*(dmuz2+dmuz1+dmux11+dmux12)
              
              if(k.eq.1) then
                pmubpx=dmux2*ari(ii,i+1,k,1)-dmux1*ari(ii,i,k,1)&
                     +dmuz2*ark(ii,i,k+1,1)-dmuz1*ark(ii,i,k,1)
                pmubpz=dmux2*ari(ii,i+1,k,2)-dmux1*ari(ii,i,k,2)&
                     +dmuz2*ark(ii,i,k+1,2)-dmuz1*ark(ii,i,k,2)
               volume=vol(ii,i,k)
              elseif(k.eq.mm(ii,2)) then
                pmubpx=dmux2*ari(ii,i+1,km,1)-dmux1*ari(ii,i,km,1)&
                     +dmuz2*ark(ii,i,k,1)-dmuz1*ark(ii,i,k-1,1)
                pmubpz=dmux2*ari(ii,i+1,km,2)-dmux1*ari(ii,i,km,2)&
                     +dmuz2*ark(ii,i,k,2)-dmuz1*ark(ii,i,k-1,2)
               volume=vol(ii,i,km)
              else
		      pmubpx=dmux2*acvk_ari(ii,i+1,k,1)-dmux1*acvk_ari(ii,i,k,1)&
         		    +dmuz2*acvk_ark(ii,i,k+1,1)-dmuz1*acvk_ark(ii,i,k,1)
		      pmubpz=dmux2*acvk_ari(ii,i+1,k,2)-dmux1*acvk_ari(ii,i,k,2)&
         		    +dmuz2*acvk_ark(ii,i,k+1,2)-dmuz1*acvk_ark(ii,i,k,2)
		      volume=acvk_vol(ii,i,k)
		      endif
    		  
		      dmuk1(ii,i,k)=pmubpx/volume
		      dmuk2(ii,i,k)=pmubpz/volume
		    enddo
		    enddo

    300   continue
	    return
	    end



        subroutine SA_source_term
        use dimen
		use grid
        use safun
        use chara
	    implicit none
        integer                   ::ii,i,k,numb,ix,iz,ia,i1,i2,im,km,kis1,kis2,kks1,kks2
        parameter(numb=1,ix=369,iz=65)
        real*8                    ::wu,wux01,wux1,wux02,wux2,wuz01,wuz1,wuz02,wuz2,wu11,wu12
        real*8                    ::cv1,cv2,cw2,cw3,cw36,cv13,rkai0,cvor,cb1,cb2,sigma,cw1,rkai,&
                                    rkai3,fv1,fv2,fv3,rk2d2,syb,syb_org,r,g,p1,p2,fw,q1,q2,q3 
        real*8                    ::dissa(numb,ix,iz),dmu1_cent(numb,ix,iz),dmu2_cent(numb,ix,iz),&
                                    strate(numb,ix,iz),wow(numb,ix,iz),qt(numb,ix,iz),qtim1(numb,ix,iz)
        common /SAdis/dissa
	    common /str_wow/strate,wow
	    common /SAsource/qt,qtim1
	    
	    
    !
    !	  to evaluate the first derivatives of the turbulent variable mu in cells
    !
    !

	    do 400 ii=1,nmp
		    im=me(ii,1)
		    km=me(ii,2)
		    kis1=0
		    kis2=0
		    kks1=0
		    kks2=0

		    do 200 k=1,km
		    do 200 i=1,im

		    wu=rmub(ii,i,k)
    
    !       i=im
    
		    if(i.eq.im) then
			    if(kb(ii,4,k).eq.0) then 
                    wux02=rmub(ii,im+1,k)
			    elseif(kb(ii,4,k).eq.1) then  !物面
				    wux02=0.
			    else
				    ia=md(ii,k,1,1,4)
				    i1=md(ii,k,2,1,4)
				    i2=md(ii,k,3,1,4)
				    wux02=.5*(rmub(ii,im,k)+rmub(ia,i1,i2))
			    endif
			    kis2=1
			    goto 10
		    endif

		    wux2=rmub(ii,i+1,k)
10	    	continue

    
    !       i=1
    
		    if(i.eq.1) then
			    if(kb(ii,3,k).eq.0)then
                    wux01=rmub(ii,0,k)
			    elseif(kb(ii,3,k).eq.1) then
				    wux01=0.
			    else
				    ia=md(ii,k,1,1,3)
				    i1=md(ii,k,2,1,3)
				    i2=md(ii,k,3,1,3)
				    wux01=.5*(rmub(ii,1,k)+rmub(ia,i1,i2))
			    endif
			    kis1=1
			    goto 20
		    endif

		    wux1=rmub(ii,i-1,k)
20  		continue

    
    !       k=km
    
		    if(k.eq.km)then
			    if(kb(ii,2,i).eq.0)then
                    wuz02=rmub(ii,i,km+1)
			    elseif(kb(ii,2,i).eq.1) then
				    wuz02=0.
			    else
				    ia=md(ii,i,1,1,2)
				    i1=md(ii,i,2,1,2)
				    i2=md(ii,i,3,1,2)
				    wuz02=.5*(rmub(ii,i,km)+rmub(ia,i1,i2))
			    endif
			    kks2=1
			    goto 30
		    endif
          
		    wuz2=rmub(ii,i,k+1)
30  		continue

    
    !       k=1
    
	        if(k.eq.1)then
		        if(kb(ii,1,i).eq.0)then
                    wuz01=rmub(ii,i,0)	
			    elseif(kb(ii,1,i).eq.1)then
				    wuz01=0d0
			    else
				    ia=md(ii,i,1,1,1)
				    i1=md(ii,i,2,1,1)
				    i2=md(ii,i,3,1,1)
				    wuz01=.5*(rmub(ii,i,1)+rmub(ia,i1,i2))
			    endif
			    kks1=1
			    goto 40
		    endif

		    wuz1=rmub(ii,i,k-1)
40  		continue

		    if(kis1.eq.1)then
                 kis1=0
                 goto 50
              endif
              wux01=(wux1+wu)/2.
50          continue

              if(kis2.eq.1)then
		      kis2=0
                goto 60
              endif
              wux02=(wux2+wu)/2.
60          continue

              if(kks1.eq.1)then
                kks1=0
                goto 70
              endif
              wuz01=(wuz1+wu)/2.
70          continue
       
              if(kks2.eq.1)then
                kks2=0
                goto 80
              endif
              wuz02=(wuz2+wu)/2.
80          continue

		      wu11=wux02*ari(ii,i+1,k,1)-wux01*ari(ii,i,k,1)&
                  +wuz02*ark(ii,i,k+1,1)-wuz01*ark(ii,i,k,1)
              wu12=wux02*ari(ii,i+1,k,2)-wux01*ari(ii,i,k,2)&
                  +wuz02*ark(ii,i,k+1,2)-wuz01*ark(ii,i,k,2)
              dmu1_cent(ii,i,k)=wu11/vol(ii,i,k)
              dmu2_cent(ii,i,k)=wu12/vol(ii,i,k)
200 		continue
400     continue
      
        cv1=7.1
	    cv2=5.
	    cw2=0.3
	    cw3=2.
	    cw36=cw3**6
	    cv13=cv1**3
	    rkai0=0.41
	    cvor=4.0
!   	cvor=0.

        cb1=0.1355
        cb2=0.622
        sigma=2./3.
        rkai0=0.41	
        cw1=cb1/(rkai0*rkai0)+(1+cb2)/sigma
	    do 500 ii=1,nmp
	    do 500 k=1,me(ii,2)
	    do 500 i=1,me(ii,1)

        rkai=rmub(ii,i,k)/rmu(ii,i,k)
        rkai3=rkai**3
        fv1=rkai3/(rkai3+cv13)
        fv2=(1+rkai/cv2)**(-3)
        fv3=(1+rkai*fv1)*(1-fv2)/max(rkai,0.001d0)

!		to solve syb 
		
        rk2d2=(rkai0*dissa(ii,i,k))**2  !dissa(ii,i,k)是到物面的最近距离
        syb=fv3*wow(ii,i,k)&
		        +cvor*dmin1(0d0,strate(ii,i,k)-wow(ii,i,k))&  !修正项？？？？？？？？？？？？
		        +rmub(ii,i,k)*fv2/rk2d2
        syb_org=fv3*wow(ii,i,k)+rmub(ii,i,k)*fv2/rk2d2

!		to solve fw

        r=rmub(ii,i,k)/(syb_org*rk2d2)
        g=r+cw2*(r**6-r)
        p1=1+cw36
        p2=g**6+cw36
        fw=g*((p1/p2)**(1./6.))

!	   the source term

	    q1=cb1*syb*rmub(ii,i,k)
	    q2=(cb2/sigma)*(dmu1_cent(ii,i,k)**2+dmu2_cent(ii,i,k)**2)
	    q3=cw1*fw*((rmub(ii,i,k)/dissa(ii,i,k))**2)         !ft2=0，ft2用来表示转捩
	    qt(ii,i,k)=(q1+q2-q3)*vol(ii,i,k)
	    	    
!
!	    point implicit treatment for source terms 	
!
!为了稳定性，根据源项修正当地时间步长（物理时间步长）
	    qtim1(ii,i,k)=-2*cw1*fw*rmub(ii,i,k)/(dissa(ii,i,k)**2) 
        open(10,file="qtiml.txt")
        write(10,*) qtim1(ii,i,k)
500     continue    
	    return
	    end




	    subroutine residsa     !计算SA模型控制方程有限体积法残值
	    use dimen
	    use chara
	    use grid
	    use safun
	    use conver
	    use turb
	    implicit none
        integer                   ::ii,i,k,numb,ix,iz,im,km,ia,i1,i2
        parameter(numb=1,ix=369,iz=65)
        real*8                    ::sigma,volm
        real*8                    ::wup5(numb,ix,iz),dmui1(numb,ix+1,iz),dmui2(numb,ix+1,iz),dmuk1(numb,ix,iz+1),&
                                    dmuk2(numb,ix,iz+1),ctrvi(numb,ix+1,iz),ctrvk(numb,ix,iz+1),qt(numb,ix,iz),qtim1(numb,ix,iz)
        real*8                    ::fctkf(ix),fctkb,fctif,fctib,q1,vfcof,fvtif(iz),fvtkf,fvtib,fvtkb,txxt,tzzt,&
                                    fct(ix,iz),fvt(ix,iz)
        common /resa/wup5
	    common /saderiv/dmui1,dmui2,dmuk1,dmuk2
	    common /contrv/ctrvi,ctrvk
	    common /SAsource/qt,qtim1
	    
	    do 1000 ii=1,nmp
	      im=me(ii,1)
	      km=me(ii,2)
    	  
    
    !	  the convective flux
    
    !为了反映出对流特性，通过边界的通量一定是由上游传到下游，也就是rmub需要取上游的rmub
    !对流方向均以i向正向与k向正向为正
	      do i=1,im
    
    !       k=1  
    	  
	      if(ctrvk(ii,i,1).le.0) then
		      fctkf(i)=ctrvk(ii,i,1)*rmub(ii,i,1) 
		  else
		    if(kb(ii,1,i).eq.0) then
              fctkf(i)=ctrvk(ii,i,1)*rmub(ii,i,0)
            elseif(kb(ii,1,i).eq.1) then
              fctkf(i)=0d0 
		    else
		      ia=md(ii,i,1,1,1)
		      i1=md(ii,i,2,1,1)
		      i2=md(ii,i,3,1,1)
		      fctkf(i)=ctrvk(ii,i,1)*rmub(ia,i1,i2)
		    endif
		  endif
          enddo
    		  	 
            do 200 k=1,km 
    
    !       i=1
               
	      if(ctrvi(ii,1,k).le.0)then  
		      fctif=ctrvi(ii,1,k)*rmub(ii,1,k) 
		  else
		    if(kb(ii,3,k).eq.0) then
		      fctif=ctrvi(ii,1,k)*rmub(ii,0,k) 
		    elseif(kb(ii,3,k).eq.1) then
		      fctif=0
		    else
		      ia=md(ii,k,1,1,3)
		      i1=md(ii,k,2,1,3)
		      i2=md(ii,k,3,1,3)
		      fctif=ctrvi(ii,1,k)*rmub(ia,i1,i2)
		    endif
          endif
          
		    do 100 i=1,im
    
    !	      i direction
    
	         if(ctrvi(ii,i+1,k).ge.0) then
		       fctib=ctrvi(ii,i+1,k)*rmub(ii,i,k)   !fctib(i单元i+1边界i向流出).gt.0 (此时ctrvi(ii,i+1,k).gt.0)
		     else
		       if(i.eq.im) then  !i=imax边界
		          if(kb(ii,4,k).eq.0) then
                  fctib=ctrvi(ii,im+1,k)*rmub(ii,im+1,k)
                  elseif(kb(ii,4,k).eq.1) then
                  fctib=0
		          else
		          ia=md(ii,k,1,1,4)
		          i1=md(ii,k,2,1,4)
		          i2=md(ii,k,3,1,4)
		          fctib=ctrvi(ii,im+1,k)*rmub(ia,i1,i2)   !fctib(从对接单元流入im+1边界).lt.0(此时ctrvi(ii,i+1,k).lt.0)
		          endif	  
		       else
		        fctib=ctrvi(ii,i+1,k)*rmub(ii,i+1,k)   !fctib(从相邻i+1单元流入i+1边界).lt.0 (此时ctrvi(ii,i+1,k).lt.0)
		       endif 
	        endif
              q1=fctib-fctif
		      fctif=fctib   !流出i单元的i+1边界通量=流进i+1单元的i+1边界通量
    
    !	 	  k direction
    
	        if(ctrvk(ii,i,k+1).ge.0) then
		        fctkb=ctrvk(ii,i,k+1)*rmub(ii,i,k)
		    else
		      if(k.eq.km) then
		        if(kb(ii,2,i).eq.0) then
		            fctkb=ctrvk(ii,i,km+1)*rmub(ii,i,km+1)
		         elseif(kb(ii,2,i).eq.1) then
		            fctkb=0
		         else
		            ia=md(ii,i,1,1,2)
		            i1=md(ii,i,2,1,2)
		            i2=md(ii,i,3,1,2)
		            fctkb=ctrvk(ii,i,km+1)*rmub(ia,i1,i2)
		        endif	  
		      else
		        fctkb=ctrvk(ii,i,k+1)*rmub(ii,i,k+1)
		      endif  
	       endif
             fct(i,k)=q1+fctkb-fctkf(i)   !流出i单元所有边界对流通量
		     fctkf(i)=fctkb
100		   continue
200	       continue


!          the viscous flux

!          i=1

	        sigma=2./3.
	        
            do 201 k=1,km
            if(kb(ii,3,k).eq.0) then
              fvtif(k)=0.
              goto 201
            elseif(kb(ii,3,k).eq.1) then
              vfcof=rmu(ii,1,k)/sigma   !物面湍流粘性系数为0，直接近似取为层流近物面单元的运动粘性系数
            elseif(kb(ii,3,k).eq.2) then
              ia=md(ii,k,1,1,3)
		      i1=md(ii,k,2,1,3)
		      i2=md(ii,k,3,1,3)
		      vfcof=0.5*(rmu(ii,1,k)+rmub(ii,1,k)+rmu(ia,i1,i2)+rmub(ia,i1,i2))/sigma 
            endif
            
            txxt=vfcof*dmui1(ii,1,k)
            tzzt=vfcof*dmui2(ii,1,k)
            fvtif(k)=txxt*ari(ii,1,k,1)+tzzt*ari(ii,1,k,2)       !迭代开始i=1边界值
201         continue


            do 500 i=1,im        !整体向i正向推进
            
    
    !       k=1
    
            if(kb(ii,1,i).eq.0) then
              fvtkf=0.
              goto 202
            elseif(kb(ii,1,i).eq.1) then
              vfcof=rmu(ii,i,1)/sigma
            elseif(kb(ii,1,i).eq.2) then
              ia=md(ii,i,1,1,1)
		      i1=md(ii,i,2,1,1)
		      i2=md(ii,i,3,1,1)
		      vfcof=0.5*(rmu(ii,i,1)+rmub(ii,i,1)+rmu(ia,i1,i2)+rmub(ia,i1,i2))/sigma
            endif
            
            txxt=vfcof*dmuk1(ii,i,1)
            tzzt=vfcof*dmuk2(ii,i,1)
            fvtkf=txxt*ark(ii,i,1,1)+tzzt*ark(ii,i,1,2)      
202         continue

            do 500 k=1,km      !在每一个i中向k正向扫描


    
    !       I direction
    

	          if(i.eq.im) then	    	
		        if(kb(ii,4,k).eq.0) then
			      fvtib=0
			      goto 356
		        elseif(kb(ii,4,k).eq.1) then
			      vfcof=rmu(ii,im,k)/sigma
			    elseif(kb(ii,4,k).eq.2) then
                  ia=md(ii,k,1,1,4)
		          i1=md(ii,k,2,1,4)
		          i2=md(ii,k,3,1,4)
		          vfcof=0.5*(rmu(ii,im,k)+rmub(ii,im,k)+rmu(ia,i1,i2)+rmub(ia,i1,i2))/sigma      
                endif 		
		        goto 355
	          endif

		      vfcof=0.5*(rmu(ii,i,k)+rmub(ii,i,k)+rmu(ii,i+1,k)+rmub(ii,i+1,k))/sigma     
355 	      continue
		      txxt=vfcof*dmui1(ii,i+1,k)	
		      tzzt=vfcof*dmui2(ii,i+1,k)
		      fvtib=txxt*ari(ii,i+1,k,1)+tzzt*ari(ii,i+1,k,2)
356           continue
		      q1=fvtib-fvtif(k)
		      fvtif(k)=fvtib   !保存前一个单元的下游边界值作为下一个单元的上游边界值

    
    !		  k direction
    
		      if(k.eq.km) then	  	
     			if(kb(ii,2,i).eq.0) then
     			  fvtkb=0d0
     			  goto 366		
     			elseif(kb(ii,2,i).eq.1) then
     			  vfcof=rmu(ii,i,km)/sigma			
			    else
			      ia=md(ii,i,1,1,2)
			      i1=md(ii,i,2,1,2)
			      i2=md(ii,i,3,1,2)
			      vfcof=0.5*(rmu(ii,i,km)+rmub(ii,i,km)+rmu(ia,i1,i2)+rmub(ia,i1,i2))/sigma             
			    endif
			    goto 365
		      endif

              vfcof=0.5*(rmu(ii,i,k)+rmub(ii,i,k)+rmu(ii,i,k+1)+rmub(ii,i,k+1))/sigma        
365		     continue

		      txxt=vfcof*dmuk1(ii,i,k+1)	
		      tzzt=vfcof*dmuk2(ii,i,k+1)
    		  fvtkb=txxt*ark(ii,i,k+1,1)+tzzt*ark(ii,i,k+1,2)
366          continue
		      fvt(i,k)=q1+fvtkb-fvtkf
		      fvtkf=fvtkb
400	        continue	
500	        continue

    
    
    !	  to resolve the residual  
    
            do k=1,km
	        do i=1,im
		    volm=vol(ii,i,k)
	        wup5(ii,i,k)=(-fct(i,k)+fvt(i,k)+qt(ii,i,k))/volm   !湍流方程残值
            enddo
            enddo

1000    continue
        return
	    end
    
    
    
    
        subroutine derivatives
        use dimen
        use grid
        use chara
        use conver
        use geom_acv
        use data0
        use qinf
        implicit none
        integer                   ::ii,i,k,numb,ix,iz,im,km,ia,i1,i2
        parameter(numb=1,ix=369,iz=65)
        real*8                    ::w1,w2,w3,w4,d,u,w,p,d0,u0,w0,p0,dx1,dx2,ux1,ux2,wx1,wx2,px1,px2,&
                                    tx1,tx2,dz11,uz11,wz11,pz11,tz11,dz12,uz12,wz12,pz12,tz12,dz21,&
                                    uz21,wz21,pz21,tz21,dz22,uz22,wz22,pz22,tz22,dz1,uz1,wz1,pz1,tz1,&
                                    dz2,uz2,wz2,pz2,tz2,x11,x12,z11,z12,t11,t12,dx11,ux11,wx11,px11,tx11,&
                                    dx12,ux12,wx12,px12,tx12,dx21,ux21,wx21,px21,tx21,dx22,ux22,wx22,px22,tx22,&
                                    volume                 
        real*8                    ::ui1(numb,ix+1,iz),ui2(numb,ix+1,iz),wi1(numb,ix+1,iz),wi2(numb,ix+1,iz),&
                                    ti1(numb,ix+1,iz),ti2(numb,ix+1,iz),uk1(numb,ix,iz+1),uk2(numb,ix,iz+1),&
                                    wK1(numb,ix,iz+1),wK2(numb,ix,iz+1),tk1(numb,ix,iz+1),tk2(numb,ix,iz+1)
         
        common /grad/ui1,ui2,wi1,wi2,ti1,ti2,uk1,uk2,wk1,wk2,tk1,tk2    !单元边界上的速度偏导
  
      
      
      
        do 100 ii=1,nmp
        im=me(ii,1)
        km=me(ii,2)

!
!     to supply values of conservative variables in the dummy cells(i=0,im+1;k=0,km+1)

!     i=0
!
        do k=1,km
        if(kb(ii,3,k).eq.1) then
        w1=wp1(ii,1,k)
        w2=wp2(ii,1,k)
        w3=wp3(ii,1,k)
        w4=wp4(ii,1,k)
        call trans(w1,w2,w3,w4,d,u,w,p)
        d0=d
        u0=-u
        w0=-w
        p0=p
        wp1(ii,0,k)=d0
        wp2(ii,0,k)=d0*u0
        wp3(ii,0,k)=d0*w0
        wp4(ii,0,k)=p0/(gma-1.)+0.5*d0*(u0*u0+w0*w0)
        elseif(kb(ii,3,k).eq.2) then
        ia=md(ii,k,1,1,3)
        i1=md(ii,k,2,1,3) 
        i2=md(ii,k,3,1,3)
        wp1(ii,0,k)=wp1(ia,i1,i2)
        wp2(ii,0,k)=wp2(ia,i1,i2)
        wp3(ii,0,k)=wp3(ia,i1,i2)
        wp4(ii,0,k)=wp4(ia,i1,i2)
        endif
        enddo 


!
!     i=im+1
!
        do k=1,km
        if(kb(ii,4,k).eq.1) then
        w1=wp1(ii,im,k)
        w2=wp2(ii,im,k)
        w3=wp3(ii,im,k)
        w4=wp4(ii,im,k)
        call trans(w1,w2,w3,w4,d,u,w,p)
        d0=d
        u0=-u
        w0=-w
        p0=p
        wp1(ii,im+1,k)=d0
        wp2(ii,im+1,k)=d0*u0
        wp3(ii,im+1,k)=d0*w0
        wp4(ii,im+1,k)=p0/(gma-1.)+0.5*d0*(u0*u0+w0*w0)
        elseif(kb(ii,4,k).eq.2) then
        ia=md(ii,k,1,1,4)
        i1=md(ii,k,2,1,4) 
        i2=md(ii,k,3,1,4)
        wp1(ii,im+1,k)=wp1(ia,i1,i2)
        wp2(ii,im+1,k)=wp2(ia,i1,i2)
        wp3(ii,im+1,k)=wp3(ia,i1,i2)
        wp4(ii,im+1,k)=wp4(ia,i1,i2)
        endif
        enddo 


!
!     k=0
!
        do i=1,im
        if(kb(ii,1,i).eq.1) then
        w1=wp1(ii,i,1)
        w2=wp2(ii,i,1)
        w3=wp3(ii,i,1)
        w4=wp4(ii,i,1)
        call trans(w1,w2,w3,w4,d,u,w,p)
        d0=d
        u0=-u
        w0=-w
        p0=p
        wp1(ii,i,0)=d0
        wp2(ii,i,0)=d0*u0
        wp3(ii,i,0)=d0*w0
        wp4(ii,i,0)=p0/(gma-1.)+0.5*d0*(u0*u0+w0*w0)
        elseif(kb(ii,1,i).eq.2) then
        ia=md(ii,i,1,1,1)
        i1=md(ii,i,2,1,1) 
        i2=md(ii,i,3,1,1)
        wp1(ii,i,0)=wp1(ia,i1,i2)
        wp2(ii,i,0)=wp2(ia,i1,i2)
        wp3(ii,i,0)=wp3(ia,i1,i2)
        wp4(ii,i,0)=wp4(ia,i1,i2)
        endif
        enddo 


!
!     k=km+1
!
        do i=1,im
        if(kb(ii,2,i).eq.1) then
        w1=wp1(ii,i,km)
        w2=wp2(ii,i,km)
        w3=wp3(ii,i,km)
        w4=wp4(ii,i,km)
        call trans(w1,w2,w3,w4,d,u,w,p)
        d0=d
        u0=-u
        w0=-w
        p0=p
        wp1(ii,i,km+1)=d0
        wp2(ii,i,km+1)=d0*u0
        wp3(ii,i,km+1)=d0*w0
        wp4(ii,i,km+1)=p0/(gma-1.)+0.5*d0*(u0*u0+w0*w0)
        elseif(kb(ii,2,i).eq.2) then
        ia=md(ii,i,1,1,2)
        i1=md(ii,i,2,1,2) 
        i2=md(ii,i,3,1,2)
        wp1(ii,i,km+1)=wp1(ia,i1,i2)
        wp2(ii,i,km+1)=wp2(ia,i1,i2)
        wp3(ii,i,km+1)=wp3(ia,i1,i2)
        wp4(ii,i,km+1)=wp4(ia,i1,i2)
        endif
        enddo
!
!	  四个角上的虚拟单元值
! 
        wp1(ii,0,0)=.5*(wp1(ii,0,1)+wp1(ii,1,0))
        wp2(ii,0,0)=.5*(wp2(ii,0,1)+wp2(ii,1,0))
        wp3(ii,0,0)=.5*(wp3(ii,0,1)+wp3(ii,1,0))
        wp4(ii,0,0)=.5*(wp4(ii,0,1)+wp4(ii,1,0))

        wp1(ii,0,km+1)=.5*(wp1(ii,0,km)+wp1(ii,1,km+1))
        wp2(ii,0,km+1)=.5*(wp2(ii,0,km)+wp2(ii,1,km+1))
        wp3(ii,0,km+1)=.5*(wp3(ii,0,km)+wp3(ii,1,km+1))
        wp4(ii,0,km+1)=.5*(wp4(ii,0,km)+wp4(ii,1,km+1))

        wp1(ii,im+1,0)=.5*(wp1(ii,im,0)+wp1(ii,im+1,1))
        wp2(ii,im+1,0)=.5*(wp2(ii,im,0)+wp2(ii,im+1,1))
        wp3(ii,im+1,0)=.5*(wp3(ii,im,0)+wp3(ii,im+1,1))
        wp4(ii,im+1,0)=.5*(wp4(ii,im,0)+wp4(ii,im+1,1))

        wp1(ii,im+1,km+1)=.5*(wp1(ii,im,km+1)+wp1(ii,im+1,km))
        wp2(ii,im+1,km+1)=.5*(wp2(ii,im,km+1)+wp2(ii,im+1,km))   
        wp3(ii,im+1,km+1)=.5*(wp3(ii,im,km+1)+wp3(ii,im+1,km))   
        wp4(ii,im+1,km+1)=.5*(wp4(ii,im,km+1)+wp4(ii,im+1,km))        

!
!	  to evaluate the first derivatives of the velocity components and of temperature in cell edges
! 

!
!	  I direction
!
        do k=1,me(ii,2)
        do i=1,mm(ii,1)

        w1=wp1(ii,i,k)
        w2=wp2(ii,i,k)
        w3=wp3(ii,i,k)
        w4=wp4(ii,i,k)
        call trans(w1,w2,w3,w4,dx2,ux2,wx2,px2)
        tx2=px2/dx2

        w1=wp1(ii,i-1,k)
        w2=wp2(ii,i-1,k)
        w3=wp3(ii,i-1,k)
        w4=wp4(ii,i-1,k)
        call trans(w1,w2,w3,w4,dx1,ux1,wx1,px1)
        tx1=px1/dx1

        w1=wp1(ii,i-1,k-1)
        w2=wp2(ii,i-1,k-1)
        w3=wp3(ii,i-1,k-1)
        w4=wp4(ii,i-1,k-1)
        call trans(w1,w2,w3,w4,dz11,uz11,wz11,pz11)
        tz11=pz11/dz11

        w1=wp1(ii,i,k-1)
        w2=wp2(ii,i,k-1)
        w3=wp3(ii,i,k-1)
        w4=wp4(ii,i,k-1)
        call trans(w1,w2,w3,w4,dz12,uz12,wz12,pz12)
        tz12=pz12/dz12
                
        w1=wp1(ii,i-1,k+1)
        w2=wp2(ii,i-1,k+1)
        w3=wp3(ii,i-1,k+1)
        w4=wp4(ii,i-1,k+1)
        call trans(w1,w2,w3,w4,dz21,uz21,wz21,pz21)
        tz21=pz21/dz21     

        w1=wp1(ii,i,k+1)
        w2=wp2(ii,i,k+1)
        w3=wp3(ii,i,k+1)
        w4=wp4(ii,i,k+1)
        call trans(w1,w2,w3,w4,dz22,uz22,wz22,pz22)
        tz22=pz22/dz22

        uz1=0.25*(ux1+ux2+uz11+uz12)
        wz1=0.25*(wx1+wx2+wz11+wz12)
        tz1=0.25*(tx1+tx2+tz11+tz12)
        uz2=0.25*(ux1+ux2+uz21+uz22)
        wz2=0.25*(wx1+wx2+wz21+wz22)
        tz2=0.25*(tx1+tx2+tz21+tz22)

        if(i.eq.1) then
        x11=ux2*ari(ii,i+1,k,1)-ux1*ari(ii,i,k,1)&
            +uz2*ark(ii,i,k+1,1)-uz1*ark(ii,i,k,1)
        x12=ux2*ari(ii,i+1,k,2)-ux1*ari(ii,i,k,2)&
            +uz2*ark(ii,i,k+1,2)-uz1*ark(ii,i,k,2)
        z11=wx2*ari(ii,i+1,k,1)-wx1*ari(ii,i,k,1)&
            +wz2*ark(ii,i,k+1,1)-wz1*ark(ii,i,k,1)
        z12=wx2*ari(ii,i+1,k,2)-wx1*ari(ii,i,k,2)&
            +wz2*ark(ii,i,k+1,2)-wz1*ark(ii,i,k,2)
        t11=tx2*ari(ii,i+1,k,1)-tx1*ari(ii,i,k,1)&
            +tz2*ark(ii,i,k+1,1)-tz1*ark(ii,i,k,1)
        t12=tx2*ari(ii,i+1,k,2)-tx1*ari(ii,i,k,2)&
            +tz2*ark(ii,i,k+1,2)-tz1*ark(ii,i,k,2)
        volume=vol(ii,i,k)        
        elseif(i.eq.mm(ii,1)) then
        x11=ux2*ari(ii,im+1,k,1)-ux1*ari(ii,im,k,1)&
            +uz2*ark(ii,im,k+1,1)-uz1*ark(ii,im,k,1)
        x12=ux2*ari(ii,im+1,k,2)-ux1*ari(ii,im,k,2)&
            +uz2*ark(ii,im,k+1,2)-uz1*ark(ii,im,k,2)
        z11=wx2*ari(ii,im+1,k,1)-wx1*ari(ii,im,k,1)&
            +wz2*ark(ii,im,k+1,1)-wz1*ark(ii,im,k,1)
        z12=wx2*ari(ii,im+1,k,2)-wx1*ari(ii,im,k,2)&
            +wz2*ark(ii,im,k+1,2)-wz1*ark(ii,im,k,2)  
        t11=tx2*ari(ii,im+1,k,1)-tx1*ari(ii,im,k,1)&
            +tz2*ark(ii,im,k+1,1)-tz1*ark(ii,im,k,1)
        t12=tx2*ari(ii,im+1,k,2)-tx1*ari(ii,im,k,2)&
            +tz2*ark(ii,im,k+1,2)-tz1*ark(ii,im,k,2)
        volume=vol(ii,im,k) 
        else
        x11=ux2*acvi_ari(ii,i+1,k,1)-ux1*acvi_ari(ii,i,k,1)&
            +uz2*acvi_ark(ii,i,k+1,1)-uz1*acvi_ark(ii,i,k,1) 
        x12=ux2*acvi_ari(ii,i+1,k,2)-ux1*acvi_ari(ii,i,k,2)&
            +uz2*acvi_ark(ii,i,k+1,2)-uz1*acvi_ark(ii,i,k,2) 
        z11=wx2*acvi_ari(ii,i+1,k,1)-wx1*acvi_ari(ii,i,k,1)&
            +wz2*acvi_ark(ii,i,k+1,1)-wz1*acvi_ark(ii,i,k,1) 
        z12=wx2*acvi_ari(ii,i+1,k,2)-wx1*acvi_ari(ii,i,k,2)&
            +wz2*acvi_ark(ii,i,k+1,2)-wz1*acvi_ark(ii,i,k,2)     
        t11=tx2*acvi_ari(ii,i+1,k,1)-tx1*acvi_ari(ii,i,k,1)&
            +tz2*acvi_ark(ii,i,k+1,1)-tz1*acvi_ark(ii,i,k,1) 
        t12=tx2*acvi_ari(ii,i+1,k,2)-tx1*acvi_ari(ii,i,k,2)&
            +tz2*acvi_ark(ii,i,k+1,2)-tz1*acvi_ark(ii,i,k,2)
        volume=acvi_vol(ii,i,k)  
        endif

        ui1(ii,i,k)=x11/volume !i向边偏u/偏x
        ui2(ii,i,k)=x12/volume !i向边偏u/偏y
        wi1(ii,i,k)=z11/volume !i向边偏w/偏x
        wi2(ii,i,k)=z12/volume !i向边偏w/偏y 
        ti1(ii,i,k)=t11/volume !i向边偏t/偏x
        ti2(ii,i,k)=t12/volume !i向边偏t/偏y

        enddo
        enddo 

!
!	  K direction
!
        do i=1,me(ii,1)
        do k=1,mm(ii,2)

        w1=wp1(ii,i,k)
        w2=wp2(ii,i,k)
        w3=wp3(ii,i,k)
        w4=wp4(ii,i,k)
        call trans(w1,w2,w3,w4,dz2,uz2,wz2,pz2)
        tz2=pz2/dz2

        w1=wp1(ii,i,k-1)
        w2=wp2(ii,i,k-1)
        w3=wp3(ii,i,k-1)
        w4=wp4(ii,i,k-1)
        call trans(w1,w2,w3,w4,dz1,uz1,wz1,pz1)
        tz1=pz1/dz1
              
        w1=wp1(ii,i-1,k)
        w2=wp2(ii,i-1,k)
        w3=wp3(ii,i-1,k)
        w4=wp4(ii,i-1,k)
        call trans(w1,w2,w3,w4,dx12,ux12,wx12,px12)
        tx12=px12/dx12

        w1=wp1(ii,i-1,k-1)
        w2=wp2(ii,i-1,k-1)
        w3=wp3(ii,i-1,k-1)
        w4=wp4(ii,i-1,k-1)
        call trans(w1,w2,w3,w4,dx11,ux11,wx11,px11)
        tx11=px11/dx11
                
        w1=wp1(ii,i+1,k)
        w2=wp2(ii,i+1,k)
        w3=wp3(ii,i+1,k)
        w4=wp4(ii,i+1,k)
        call trans(w1,w2,w3,w4,dx22,ux22,wx22,px22)
        tx22=px22/dx22     

        w1=wp1(ii,i+1,k-1)
        w2=wp2(ii,i+1,k-1)
        w3=wp3(ii,i+1,k-1)
        w4=wp4(ii,i+1,k-1)
        call trans(w1,w2,w3,w4,dx21,ux21,wx21,px21)
        tx21=px21/dx21

        ux1=0.25*(uz1+uz2+ux11+ux12)
        wx1=0.25*(wz1+wz2+wx11+wx12)
        tx1=0.25*(tz1+tz2+tx11+tx12)
        ux2=0.25*(uz1+uz2+ux21+ux22)
        wx2=0.25*(wz1+wz2+wx21+wx22)
        tx2=0.25*(tz1+tz2+tx21+tx22)

        if(k.eq.1) then
        x11=ux2*ari(ii,i+1,k,1)-ux1*ari(ii,i,k,1)&
            +uz2*ark(ii,i,k+1,1)-uz1*ark(ii,i,k,1)
        x12=ux2*ari(ii,i+1,k,2)-ux1*ari(ii,i,k,2)&
            +uz2*ark(ii,i,k+1,2)-uz1*ark(ii,i,k,2)
        z11=wx2*ari(ii,i+1,k,1)-wx1*ari(ii,i,k,1)&
            +wz2*ark(ii,i,k+1,1)-wz1*ark(ii,i,k,1)
        z12=wx2*ari(ii,i+1,k,2)-wx1*ari(ii,i,k,2)&
            +wz2*ark(ii,i,k+1,2)-wz1*ark(ii,i,k,2)
        t11=tx2*ari(ii,i+1,k,1)-tx1*ari(ii,i,k,1)&
            +tz2*ark(ii,i,k+1,1)-tz1*ark(ii,i,k,1)
        t12=tx2*ari(ii,i+1,k,2)-tx1*ari(ii,i,k,2)&
            +tz2*ark(ii,i,k+1,2)-tz1*ark(ii,i,k,2)
        volume=vol(ii,i,k)        
        elseif(k.eq.mm(ii,2)) then
        x11=ux2*ari(ii,i+1,km,1)-ux1*ari(ii,i,km,1)&
            +uz2*ark(ii,i,km+1,1)-uz1*ark(ii,i,km,1)
        x12=ux2*ari(ii,i+1,km,2)-ux1*ari(ii,i,km,2)&
            +uz2*ark(ii,i,km+1,2)-uz1*ark(ii,i,km,2)
        z11=wx2*ari(ii,i+1,km,1)-wx1*ari(ii,i,km,1)&
            +wz2*ark(ii,i,km+1,1)-wz1*ark(ii,i,km,1)
        z12=wx2*ari(ii,i+1,km,2)-wx1*ari(ii,i,km,2)&
            +wz2*ark(ii,i,km+1,2)-wz1*ark(ii,i,km,2)  
        t11=tx2*ari(ii,i+1,km,1)-tx1*ari(ii,i,km,1)&
            +tz2*ark(ii,i,km+1,1)-tz1*ark(ii,i,km,1)
        t12=tx2*ari(ii,i+1,km,2)-tx1*ari(ii,i,km,2)&
            +tz2*ark(ii,i,km+1,2)-tz1*ark(ii,i,km,2)
        volume=vol(ii,i,km) 
        else
        x11=ux2*acvk_ari(ii,i+1,k,1)-ux1*acvk_ari(ii,i,k,1)&
            +uz2*acvk_ark(ii,i,k+1,1)-uz1*acvk_ark(ii,i,k,1) 
        x12=ux2*acvk_ari(ii,i+1,k,2)-ux1*acvk_ari(ii,i,k,2)&
            +uz2*acvk_ark(ii,i,k+1,2)-uz1*acvk_ark(ii,i,k,2) 
        z11=wx2*acvk_ari(ii,i+1,k,1)-wx1*acvk_ari(ii,i,k,1)&
            +wz2*acvk_ark(ii,i,k+1,1)-wz1*acvk_ark(ii,i,k,1) 
        z12=wx2*acvk_ari(ii,i+1,k,2)-wx1*acvk_ari(ii,i,k,2)&
            +wz2*acvk_ark(ii,i,k+1,2)-wz1*acvk_ark(ii,i,k,2)     
        t11=tx2*acvk_ari(ii,i+1,k,1)-tx1*acvk_ari(ii,i,k,1)&
            +tz2*acvk_ark(ii,i,k+1,1)-tz1*acvk_ark(ii,i,k,1) 
        t12=tx2*acvk_ari(ii,i+1,k,2)-tx1*acvk_ari(ii,i,k,2)&
            +tz2*acvk_ark(ii,i,k+1,2)-tz1*acvk_ark(ii,i,k,2)
        volume=acvk_vol(ii,i,k)  
        endif

        uk1(ii,i,k)=x11/volume !k向边偏u/偏x
        uk2(ii,i,k)=x12/volume !k向边偏u/偏y
        wk1(ii,i,k)=z11/volume !k向边偏w/偏x
        wk2(ii,i,k)=z12/volume !k向边偏w/偏y 
        tk1(ii,i,k)=t11/volume !k向边偏t/偏x
        tk2(ii,i,k)=t12/volume !k向边偏t/偏y         
        enddo
        enddo 

100     continue
        return
        end
    
    
    
    
        subroutine vissi   !计算层流(紊流)粘性系数和粘性通量项及其在单元边界上的积分
        use dimen
        use grid
        use chara
        use conver
        use qinf
        use data0
        use viscous
        use vissflux
        use geom_acv
        implicit none
        integer                   ::ii,i,k,numb,ix,iz,im,km,ia,i1,i2
        parameter(numb=1,ix=369,iz=65)
        real*8                    ::w1,w2,w3,w4,d,u,w,p,niu00,kk00,niu01,niu02,kk01,kk02,prl,prlt,&
                                    d01,u01,w01,p01,d02,u02,w02,p02,fskv11,fskv21,fskv31,fskv41,&
                                    fskv12,fskv22,fskv32,fskv42,fsiv11,fsiv21,fsiv31,fsiv41,fsiv12,fsiv22,&
                                    fsiv32,fsiv42,fvif1,fvif2,fvif3,fvif4,fvib1,fvib2,fvib3,fvib4,fvkb1,&
                                    fvkb2,fvkb3,fvkb4,q1,q2,q3,q4,deltu
        real*8                    ::ui1(numb,ix+1,iz),ui2(numb,ix+1,iz),wi1(numb,ix+1,iz),wi2(numb,ix+1,iz),&
                                    ti1(numb,ix+1,iz),ti2(numb,ix+1,iz),uk1(numb,ix,iz+1),uk2(numb,ix,iz+1),&
                                    wK1(numb,ix,iz+1),wK2(numb,ix,iz+1),tk1(numb,ix,iz+1),tk2(numb,ix,iz+1)
        real*8                    ::fvkf1(ix),fvkf2(ix),fvkf3(ix),fvkf4(ix)
        common /grad/ui1,ui2,wi1,wi2,ti1,ti2,uk1,uk2,wk1,wk2,tk1,tk2    !单元边界上的速度偏导




        prl=0.72
        prlt=0.9

        do 1000 ii=1,nmp
        im=me(ii,1)
        km=me(ii,2)

!
!       k=1
!	  
        do i=1,im

        if(kb(ii,1,i).eq.0) then
        w1=wp1(ii,i,0)
        w2=wp2(ii,i,0)
        w3=wp3(ii,i,0)
        w4=wp4(ii,i,0)
        call trans(w1,w2,w3,w4,d,u,w,p)
        niu00=ul(ii,i,1)+ut(ii,i,1)
        kk00=gma*ul(ii,i,1)/(prl*(gma-1.))+gma*ut(ii,i,1)/(prlt*(gma-1.))
        elseif(kb(ii,1,i).eq.1) then
        u=0
        w=0
        niu00=ul(ii,i,1)  !物面湍流粘性系数为0
        kk00=0.     !绝热壁面
        elseif(kb(ii,1,i).eq.2) then
        w1=wp1(ii,i,1)
        w2=wp2(ii,i,1)
        w3=wp3(ii,i,1)
        w4=wp4(ii,i,1)
        call trans(w1,w2,w3,w4,d01,u01,w01,p01)
        niu01=ul(ii,i,1)+ut(ii,i,1)
        kk01=gma*ul(ii,i,1)/(prl*(gma-1.))+gma*ut(ii,i,1)/(prlt*(gma-1.))
        ia=md(ii,i,1,1,1)
        i1=md(ii,i,2,1,1)
        i2=md(ii,i,3,1,1)
        w1=wp1(ia,i1,i2)
        w2=wp2(ia,i1,i2)
        w3=wp3(ia,i1,i2) 
        w4=wp4(ia,i1,i2)
        call trans(w1,w2,w3,w4,d02,u02,w02,p02)
        niu02=ul(ia,i1,i2)+ut(ia,i1,i2)
        kk02=gma*ul(ia,i1,i2)/(prl*(gma-1.))+gma*ut(ia,i1,i2)/(prlt*(gma-1.))
        u=0.5*(u01+u02)
        w=0.5*(w01+w02)
        niu00=0.5*(niu01+niu02)
        kk00=0.5*(kk01+kk02)
        endif

        deltu=2/3.*(uk1(ii,i,1)+wk2(ii,i,1))
        fskv11=0.
        fskv21=niu00*(2*uk1(ii,i,1)-deltu)  !Tauxx
        fskv31=niu00*(uk2(ii,i,1)+wk1(ii,i,1))  !Tauyx
        fskv41=u*fskv21+w*fskv31+kk00*tk1(ii,i,1)

        fskv12=0
        fskv22=fskv31  !Tauxy=Tauyx
        fskv32=niu00*(2*wk2(ii,i,1)-deltu)  !Tauyy
        fskv42=u*fskv22+w*fskv32+kk00*tk2(ii,i,1)

        fvkf1(i)=fskv11*ark(ii,i,1,1)+fskv12*ark(ii,i,1,2)
        fvkf2(i)=fskv21*ark(ii,i,1,1)+fskv22*ark(ii,i,1,2)
        fvkf3(i)=fskv31*ark(ii,i,1,1)+fskv32*ark(ii,i,1,2)
        fvkf4(i)=fskv41*ark(ii,i,1,1)+fskv42*ark(ii,i,1,2)

        enddo       

        do 500 k=1,km
!
!       i=1
!   
        if(kb(ii,3,k).eq.0) then
        w1=wp1(ii,0,k)
        w2=wp2(ii,0,k)
        w3=wp3(ii,0,k)
        w4=wp4(ii,0,k)
        call trans(w1,w2,w3,w4,d,u,w,p)
        niu00=ul(ii,1,k)+ut(ii,1,k)
        kk00=gma*ul(ii,1,k)/(prl*(gma-1.))+gma*ut(ii,1,k)/(prlt*(gma-1.))
        elseif(kb(ii,3,k).eq.1) then
        u=0
        w=0
        niu00=ul(ii,1,k)
        kk00=0
        elseif(kb(ii,3,k).eq.2) then
        w1=wp1(ii,1,k)
        w2=wp2(ii,1,k)
        w3=wp3(ii,1,k)
        w4=wp4(ii,1,k)
        call trans(w1,w2,w3,w4,d01,u01,w01,p01)
        niu01=ul(ii,1,k)+ut(ii,1,k)
        kk01=gma*ul(ii,1,k)/(prl*(gma-1.))+gma*ut(ii,1,k)/(prlt*(gma-1.))
        ia=md(ii,k,1,1,3)
        i1=md(ii,k,2,1,3)
        i2=md(ii,k,3,1,3)
        w1=wp1(ia,i1,i2)
        w2=wp2(ia,i1,i2)
        w3=wp3(ia,i1,i2)
        w4=wp4(ia,i1,i2)
        call trans(w1,w2,w3,w4,d02,u02,w02,p02)
        niu02=ul(ia,i1,i2)+ut(ia,i1,i2)
        kk02=gma*ul(ia,i1,i2)/(prl*(gma-1.))+gma*ut(ia,i1,i2)/(prlt*(gma-1.))

        u=0.5*(u01+u02)
        w=0.5*(w01+w02)
        niu00=0.5*(niu01+niu02)
        kk00=0.5*(kk01+kk02)
        endif

        deltu=2/3.*(ui1(ii,1,k)+wi2(ii,1,k))
        fsiv11=0.
        fsiv21=niu00*(2*ui1(ii,1,k)-deltu)  !Tauxx
        fsiv31=niu00*(ui2(ii,1,k)+wi1(ii,1,k))  !Tauyx
        fsiv41=u*fsiv21+w*fsiv31+kk00*ti1(ii,1,k)

        fsiv12=0
        fsiv22=fsiv31  !Tauxy=Tauyx
        fsiv32=niu00*(2*wi2(ii,1,k)-deltu)  !Tauyy
        fsiv42=u*fsiv22+w*fsiv32+kk00*ti2(ii,1,k)

        fvif1=fsiv11*ari(ii,1,k,1)+fsiv12*ari(ii,1,k,2)
        fvif2=fsiv21*ari(ii,1,k,1)+fsiv22*ari(ii,1,k,2)
        fvif3=fsiv31*ari(ii,1,k,1)+fsiv32*ari(ii,1,k,2)
        fvif4=fsiv41*ari(ii,1,k,1)+fsiv42*ari(ii,1,k,2)

        do 200 i=1,im

!
!      I direction
!
        if(i.eq.im) then
        if(kb(ii,4,k).eq.0) then
        w1=wp1(ii,im+1,k)
        w2=wp2(ii,im+1,k)
        w3=wp3(ii,im+1,k)
        w4=wp4(ii,im+1,k)
        call trans(w1,w2,w3,w4,d,u,w,p)
        niu00=ul(ii,im,k)+ut(ii,im,k)
        kk00=gma*ul(ii,im,k)/(prl*(gma-1.))+gma*ut(ii,im,k)/(prlt*(gma-1.))
        elseif(kb(ii,4,k).eq.1) then
        u=0
        w=0
        niu00=ul(ii,im,k)
        kk00=0
        elseif(kb(ii,4,k).eq.2) then
        w1=wp1(ii,im,k)
        w2=wp2(ii,im,k)
        w3=wp3(ii,im,k)
        w4=wp4(ii,im,k)
        call trans(w1,w2,w3,w4,d01,u01,w01,p01)
        niu01=ul(ii,im,k)+ut(ii,im,k)
        kk01=gma*ul(ii,im,k)/(prl*(gma-1.))+gma*ut(ii,im,k)/(prlt*(gma-1.))
        ia=md(ii,k,1,1,4)
        i1=md(ii,k,2,1,4)
        i2=md(ii,k,3,1,4)
        w1=wp1(ia,i1,i2)
        w2=wp2(ia,i1,i2)
        w3=wp3(ia,i1,i2)
        w4=wp4(ia,i1,i2)
        call trans(w1,w2,w3,w4,d02,u02,w02,p02)
        niu02=ul(ia,i1,i2)+ut(ia,i1,i2)
        kk02=gma*ul(ia,i1,i2)/(prl*(gma-1.))+gma*ut(ia,i1,i2)/(prlt*(gma-1.))

        u=0.5*(u01+u02)
        w=0.5*(w01+w02)
        niu00=0.5*(niu01+niu02)
        kk00=0.5*(kk01+kk02)
        endif
        goto 210
        endif

        w1=wp1(ii,i,k)
        w2=wp2(ii,i,k)
        w3=wp3(ii,i,k)
        w4=wp4(ii,i,k)
        call trans(w1,w2,w3,w4,d01,u01,w01,p01)
        niu01=ul(ii,i,k)+ut(ii,i,k)
        kk01=gma*ul(ii,i,k)/(prl*(gma-1.))+gma*ut(ii,i,k)/(prlt*(gma-1.))
        w1=wp1(ii,i+1,k)
        w2=wp2(ii,i+1,k)
        w3=wp3(ii,i+1,k)
        w4=wp4(ii,i+1,k)
        call trans(w1,w2,w3,w4,d02,u02,w02,p02)
        niu02=ul(ii,i+1,k)+ut(ii,i+1,k)
        kk02=gma*ul(ii,i+1,k)/(prl*(gma-1.))+gma*ut(ii,i+1,k)/(prlt*(gma-1.))
        u=0.5*(u01+u02)
        w=0.5*(w01+w02)
        niu00=0.5*(niu01+niu02)
        kk00=0.5*(kk01+kk02)
        
210     continue

        deltu=2/3.*(ui1(ii,i+1,k)+wi2(ii,i+1,k))
        fsiv11=0.
        fsiv21=niu00*(2*ui1(ii,i+1,k)-deltu)  !Tauxx
        fsiv31=niu00*(ui2(ii,i+1,k)+wi1(ii,i+1,k))  !Tauyx
        fsiv41=u*fsiv21+w*fsiv31+kk00*ti1(ii,i+1,k)

        fsiv12=0
        fsiv22=fsiv31  !Tauxy=Tauyx
        fsiv32=niu00*(2*wi2(ii,i+1,k)-deltu)  !Tauyy
        fsiv42=u*fsiv22+w*fsiv32+kk00*ti2(ii,i+1,k)

        fvib1=fsiv11*ari(ii,i+1,k,1)+fsiv12*ari(ii,i+1,k,2)
        fvib2=fsiv21*ari(ii,i+1,k,1)+fsiv22*ari(ii,i+1,k,2)
        fvib3=fsiv31*ari(ii,i+1,k,1)+fsiv32*ari(ii,i+1,k,2)
        fvib4=fsiv41*ari(ii,i+1,k,1)+fsiv42*ari(ii,i+1,k,2)

        q1=fvib1-fvif1
        q2=fvib2-fvif2
        q3=fvib3-fvif3
        q4=fvib4-fvif4
        fvif1=fvib1
        fvif2=fvib2
        fvif3=fvib3
        fvif4=fvib4

!
!     K direction
!
        if(k.eq.km) then
        if(kb(ii,2,i).eq.0) then
        w1=wp1(ii,i,km+1)
        w2=wp2(ii,i,km+1)
        w3=wp3(ii,i,km+1)
        w4=wp4(ii,i,km+1)
        call trans(w1,w2,w3,w4,d,u,w,p)
        niu00=ul(ii,i,km)+ut(ii,i,km)
        kk00=gma*ul(ii,i,km)/(prl*(gma-1.))+gma*ut(ii,i,km)/(prlt*(gma-1.))
        elseif(kb(ii,2,i).eq.1) then
        u=0
        w=0
        niu00=ul(ii,i,km)
        kk00=0
        elseif(kb(ii,2,i).eq.2) then 
        w1=wp1(ii,i,km)
        w2=wp2(ii,i,km)
        w3=wp3(ii,i,km)
        w4=wp4(ii,i,km)
        call trans(w1,w2,w3,w4,d01,u01,w01,p01)
        niu01=ul(ii,i,km)+ut(ii,i,km)
        kk01=gma*ul(ii,i,km)/(prl*(gma-1.))+gma*ut(ii,i,km)/(prlt*(gma-1.))
        ia=md(ii,i,1,1,2)
        i1=md(ii,i,2,1,2)
        i2=md(ii,i,3,1,2)
        w1=wp1(ia,i1,i2)
        w2=wp2(ia,i1,i2)
        w3=wp3(ia,i1,i2)
        w4=wp4(ia,i1,i2)
        call trans(w1,w2,w3,w4,d02,u02,w02,p02)
        niu02=ul(ia,i1,i2)+ut(ia,i1,i2)
        kk02=gma*ul(ia,i1,i2)/(prl*(gma-1.))+gma*ut(ia,i1,i2)/(prlt*(gma-1.))
        u=0.5*(u01+u02)
        w=0.5*(w01+w02)
        niu00=0.5*(niu01+niu02)
        kk00=0.5*(kk01+kk02)
        endif
        goto 220
        endif

        w1=wp1(ii,i,k)
        w2=wp2(ii,i,k)
        w3=wp3(ii,i,k)
        w4=wp4(ii,i,k)
        call trans(w1,w2,w3,w4,d01,u01,w01,p01)
        niu01=ul(ii,i,k)+ut(ii,i,k)
        kk01=gma*ul(ii,i,k)/(prl*(gma-1.))+gma*ut(ii,i,k)/(prlt*(gma-1.))
        w1=wp1(ii,i,k+1)
        w2=wp2(ii,i,k+1)
        w3=wp3(ii,i,k+1)
        w4=wp4(ii,i,k+1)
        call trans(w1,w2,w3,w4,d02,u02,w02,p02)
        niu02=ul(ii,i,k+1)+ut(ii,i,k+1)
        kk02=gma*ul(ii,i,k+1)/(prl*(gma-1.))+gma*ut(ii,i,k+1)/(prlt*(gma-1.))
        u=0.5*(u01+u02)
        w=0.5*(w01+w02)
        niu00=0.5*(niu01+niu02)
        kk00=0.5*(kk01+kk02)         

220     continue
        deltu=2/3.*(uk1(ii,i,k+1)+wk2(ii,i,k+1))
        fskv11=0.
        fskv21=niu00*(2*uk1(ii,i,k+1)-deltu)  !Tauxx
        fskv31=niu00*(uk2(ii,i,k+1)+wk1(ii,i,k+1))  !Tauyx
        fskv41=u*fskv21+w*fskv31+kk00*tk1(ii,i,k+1)

        fskv12=0
        fskv22=fskv31  !Tauxy=Tauyx
        fskv32=niu00*(2*wk2(ii,i,k+1)-deltu)  !Tauyy
        fskv42=u*fskv22+w*fskv32+kk00*tk2(ii,i,k+1)

        fvkb1=fskv11*ark(ii,i,k+1,1)+fskv12*ark(ii,i,k+1,2)
        fvkb2=fskv21*ark(ii,i,k+1,1)+fskv22*ark(ii,i,k+1,2)
        fvkb3=fskv31*ark(ii,i,k+1,1)+fskv32*ark(ii,i,k+1,2)
        fvkb4=fskv41*ark(ii,i,k+1,1)+fskv42*ark(ii,i,k+1,2)
        fvis(ii,i,k,1)=q1+fvkb1-fvkf1(i)
        fvis(ii,i,k,2)=q2+fvkb2-fvkf2(i)
        fvis(ii,i,k,3)=q3+fvkb3-fvkf3(i)
        fvis(ii,i,k,4)=q4+fvkb4-fvkf4(i)
        fvkf1(i)=fvkb1
        fvkf2(i)=fvkb2
        fvkf3(i)=fvkb3
        fvkf4(i)=fvkb4
200     continue
500     continue   
1000	continue
        return
        end
        
        
        subroutine dissi   !计算中心格式有限体积法的人工耗散项
        use dimen
        use grid
        use chara
        use data0
        use cfdiss
        use conver
        use kswitch
        use viscous
        implicit none
        integer                   ::izone,i,k,numb,ix,iz,im,km,ia1,i11,i21,ibh,ia2,i12,i22,&
                                    i1,i2,i3,k1,k2,k3
        parameter(numb=1,ix=369,iz=65)
        real*8                    ::w1,w2,w3,w4,d,u,w,pssp,p2,p4,p6,d11,d12,d13,d14,d21,d22,d23,d24,&
                                    d31,d32,d33,d34,dx1m,dx2m,dx3m,dx4m,dx1p,dx2p,dx3p,dx4p,vfor,vbac,&
                                    ep2m,ep4m,ep2p,ep4p,volm,dti(numb,ix,iz),dtk(numb,ix,iz)

        common /dtik/dti,dtk



        do 100 izone=1,nmp
        im=me(izone,1)
        km=me(izone,2)

!       k-direction

        do 10 i=1,im

!	    k=1

        if(kb(izone,1,i).eq.0) then   !k=1边界 远场
        dx1m=0.
        dx2m=0.
        dx3m=0.
        dx4m=0.
        w1=wp1(izone,i,0)
        w2=wp2(izone,i,0)
        w3=wp3(izone,i,0)
        w4=wp4(izone,i,0)
        call trans(w1,w2,w3,w4,d,u,w,pssp)
        p2=pssp
        p4=pre(izone,i,1)   !K=1单元
        p6=pre(izone,i,2)   !k=2单元
        vfor=dabs(p2-2.*p4+p6)/(dabs(p2)+2.*dabs(p4)+dabs(p6))   !压强梯度探测器
        d21=wp1(izone,i,2)-wp1(izone,i,1)
        d22=wp2(izone,i,2)-wp2(izone,i,1)
        d23=wp3(izone,i,2)-wp3(izone,i,1)
        d24=wp4(izone,i,2)-wp4(izone,i,1)
        d11=wp1(izone,i,1)-wp1(izone,i,0)
        d12=wp2(izone,i,1)-wp2(izone,i,0)
        d13=wp3(izone,i,1)-wp3(izone,i,0)
        d14=wp4(izone,i,1)-wp4(izone,i,0)
        
        else if(kb(izone,1,i).eq.1) then   !k=1物面 虚拟单元镜像
        dx1m=0.
        dx2m=0.
        dx3m=0.
        dx4m=0.
        if(kns.eq.1) then
        p2=pre(izone,i,1)
        else
        p2=2*pre(izone,i,1)-pre(izone,i,2)
        endif
        p4=pre(izone,i,1)
        p6=pre(izone,i,2)
        vfor=dabs(p2-2.*p4+p6)/(dabs(p2)+2.*dabs(p4)+dabs(p6))
        d21=wp1(izone,i,2)-wp1(izone,i,1)
        d22=wp2(izone,i,2)-wp2(izone,i,1)
        d23=wp3(izone,i,2)-wp3(izone,i,1)
        d24=wp4(izone,i,2)-wp4(izone,i,1)
        if(kns.eq.1) then
        d11=0.d0
        d12=2*wp2(izone,i,1)
        d13=2*wp3(izone,i,1)
        d14=0.d0
        else
        d11=d21
        d12=d22
        d13=d23
        d14=d24
        endif
        else
        ia1=md(izone,i,1,1,1)
        i11=md(izone,i,2,1,1)
        i21=md(izone,i,3,1,1)
        ibh=md(izone,i,4,1,1)
        ia2=md(izone,i,1,2,1)
        i12=md(izone,i,2,2,1)
        i22=md(izone,i,3,2,1)

        p2=pre(ia2,i12,i22)
        p4=pre(ia1,i11,i21)
        p6=pre(izone,i,1)
        vfor=dabs(p2-2.*p4+p6)/(dabs(p2)+2.*dabs(p4)+dabs(p6))  !压力梯度探测器vi
        p2=pre(ia1,i11,i21)
        p4=pre(izone,i,1)
        p6=pre(izone,i,2)
        vbac=dabs(p2-2.*p4+p6)/(dabs(p2)+2.*dabs(p4)+dabs(p6))  !压力梯度探测器v(i+1)
        d31=wp1(izone,i,2)-wp1(izone,i,1)
        d32=wp2(izone,i,2)-wp2(izone,i,1)
        d33=wp3(izone,i,2)-wp3(izone,i,1)
        d34=wp4(izone,i,2)-wp4(izone,i,1)
        d21=wp1(izone,i,1)-wp1(ia1,i11,i21)
        d22=wp2(izone,i,1)-wp2(ia1,i11,i21)
        d23=wp3(izone,i,1)-wp3(ia1,i11,i21)
        d24=wp4(izone,i,1)-wp4(ia1,i11,i21)
        d11=wp1(ia1,i11,i21)-wp1(ia2,i12,i22)
        d12=wp2(ia1,i11,i21)-wp2(ia2,i12,i22)
        d13=wp3(ia1,i11,i21)-wp3(ia2,i12,i22)
        d14=wp4(ia1,i11,i21)-wp4(ia2,i12,i22)
        ep2m=ek2*dmax1(vbac,vfor)
        ep4m=dmax1(dble(0.),(ek4-ep2m))
        if(ibh.eq.1.or.ibh.eq.2) then
        volm=.5*(dtk(izone,i,1)+dtk(ia1,i11,i21))
        else
        volm=.5*(dtk(izone,i,1)+dti(ia1,i11,i21))
        endif
        dx1m=ep2m*d21
        dx2m=ep2m*d22
        dx3m=ep2m*d23
        dx4m=ep2m*d24
        dx1m=(dx1m-ep4m*(d31-2.*d21+d11))*volm
        dx2m=(dx2m-ep4m*(d32-2.*d22+d12))*volm
        dx3m=(dx3m-ep4m*(d33-2.*d23+d13))*volm
        dx4m=(dx4m-ep4m*(d34-2.*d24+d14))*volm
        d11=d21
        d12=d22
        d13=d23
        d14=d24
        d21=d31
        d22=d32
        d23=d33
        d24=d34
        vfor=vbac
        endif

        do 10 k=1,km/2
        k1=k
        k2=k+1
        k3=k+2
        p2=pre(izone,i,k1)
        p4=pre(izone,i,k2)
        p6=pre(izone,i,k3)	  
        vbac=dabs(p2-2.*p4+p6)/(dabs(p2)+2.*dabs(p4)+dabs(p6))
        d31=wp1(izone,i,k3)-wp1(izone,i,k2)
        d32=wp2(izone,i,k3)-wp2(izone,i,k2)
        d33=wp3(izone,i,k3)-wp3(izone,i,k2)
        d34=wp4(izone,i,k3)-wp4(izone,i,k2)
        volm=.5*(dtk(izone,i,k1)+dtk(izone,i,k2))
        ep2p=ek2*dmax1(vbac,vfor)
        ep4p=dmax1(dble(0.),(ek4-ep2p))
        dx1p=ep2p*d21
        dx2p=ep2p*d22
        dx3p=ep2p*d23
        dx4p=ep2p*d24
        dx1p=(dx1p-ep4p*(d31-2.*d21+d11))*volm
        dx2p=(dx2p-ep4p*(d32-2.*d22+d12))*volm
        dx3p=(dx3p-ep4p*(d33-2.*d23+d13))*volm
        dx4p=(dx4p-ep4p*(d34-2.*d24+d14))*volm
        diss(izone,i,k,1)=dx1p-dx1m
        diss(izone,i,k,2)=dx2p-dx2m
        diss(izone,i,k,3)=dx3p-dx3m
        diss(izone,i,k,4)=dx4p-dx4m
        vfor=vbac
        d11=d21
        d12=d22
        d13=d23
        d14=d24
        d21=d31
        d22=d32
        d23=d33
        d24=d34
        dx1m=dx1p
        dx2m=dx2p
        dx3m=dx3p
        dx4m=dx4p
10      continue


        do 20 i=1,im

!   	k=km+1


        if(kb(izone,2,i).eq.0) then   !k=km+1边界 远场
        dx1m=0.
        dx2m=0.
        dx3m=0.
        dx4m=0.
        w1=wp1(izone,i,km+1)
        w2=wp2(izone,i,km+1)
        w3=wp3(izone,i,km+1)
        w4=wp4(izone,i,km+1)
        call trans(w1,w2,w3,w4,d,u,w,pssp)
        p2=pssp
        p4=pre(izone,i,km)   !K=km单元
        p6=pre(izone,i,km-1)   !k=km-1单元
        vfor=dabs(p2-2.*p4+p6)/(dabs(p2)+2.*dabs(p4)+dabs(p6))   !压强梯度探测器
        d21=wp1(izone,i,km-1)-wp1(izone,i,km)
        d22=wp2(izone,i,km-1)-wp2(izone,i,km)
        d23=wp3(izone,i,km-1)-wp3(izone,i,km)
        d24=wp4(izone,i,km-1)-wp4(izone,i,km)
        d11=d21
        d12=d22
        d13=d23
        d14=d24
        else if(kb(izone,2,i).eq.1) then   !k=km+1物面
        dx1m=0.
        dx2m=0.
        dx3m=0.
        dx4m=0.        
        p6=pre(izone,i,km-1)
        p4=pre(izone,i,km)
        if(kns.eq.1) then
        p2=pre(izone,i,km)
        else
        p2=2*p4-p6
        endif
        vfor=dabs(p2-2.*p4+p6)/(dabs(p2)+2.*dabs(p4)+dabs(p6))
        d21=wp1(izone,i,km-1)-wp1(izone,i,km)
        d22=wp2(izone,i,km-1)-wp2(izone,i,km)
        d23=wp3(izone,i,km-1)-wp3(izone,i,km)
        d24=wp4(izone,i,km-1)-wp4(izone,i,km)
        if(kns.eq.1) then
        d11=0.d0
        d12=2*wp2(izone,i,km)
        d13=2*wp3(izone,i,km)
        d14=0.d0
        else
        d11=d21
        d12=d22
        d13=d23
        d14=d24
        endif
        else
        ia1=md(izone,i,1,1,2)
        i11=md(izone,i,2,1,2)
        i21=md(izone,i,3,1,2)
        ibh=md(izone,i,4,1,2)
        ia2=md(izone,i,1,2,2)
        i12=md(izone,i,2,2,2)
        i22=md(izone,i,3,2,2)

        p2=pre(ia2,i12,i22)
        p4=pre(ia1,i11,i21)
        p6=pre(izone,i,km)
        vfor=dabs(p2-2.*p4+p6)/(dabs(p2)+2.*dabs(p4)+dabs(p6))  !压力梯度探测器vi
        p2=pre(ia1,i11,i21)
        p4=pre(izone,i,km)
        p6=pre(izone,i,km-1)
        vbac=dabs(p2-2.*p4+p6)/(dabs(p2)+2.*dabs(p4)+dabs(p6))  !压力梯度探测器v(i+1)
        d31=wp1(izone,i,km-1)-wp1(izone,i,km)
        d32=wp2(izone,i,km-1)-wp2(izone,i,km)
        d33=wp3(izone,i,km-1)-wp3(izone,i,km)
        d34=wp4(izone,i,km-1)-wp4(izone,i,km)
        d21=wp1(izone,i,km)-wp1(ia1,i11,i21)
        d22=wp2(izone,i,km)-wp2(ia1,i11,i21)
        d23=wp3(izone,i,km)-wp3(ia1,i11,i21)
        d24=wp4(izone,i,km)-wp4(ia1,i11,i21)
        d11=wp1(ia1,i11,i21)-wp1(ia2,i12,i22)
        d12=wp2(ia1,i11,i21)-wp2(ia2,i12,i22)
        d13=wp3(ia1,i11,i21)-wp3(ia2,i12,i22)
        d14=wp4(ia1,i11,i21)-wp4(ia2,i12,i22)
        ep2m=ek2*dmax1(vbac,vfor)
        ep4m=dmax1(dble(0.),(ek4-ep2m))
        if(ibh.eq.1.or.ibh.eq.2) then
        volm=.5*(dtk(izone,i,km)+dtk(ia1,i11,i21))
        else
        volm=.5*(dtk(izone,i,km)+dti(ia1,i11,i21))
        endif
        dx1m=ep2m*d21
        dx2m=ep2m*d22
        dx3m=ep2m*d23
        dx4m=ep2m*d24
        dx1m=(dx1m-ep4m*(d31-2.*d21+d11))*volm
        dx2m=(dx2m-ep4m*(d32-2.*d22+d12))*volm
        dx3m=(dx3m-ep4m*(d33-2.*d23+d13))*volm
        dx4m=(dx4m-ep4m*(d34-2.*d24+d14))*volm
        d11=d21
        d12=d22
        d13=d23
        d14=d24
        d21=d31
        d22=d32
        d23=d33
        d24=d34
        vfor=vbac
        endif
   
        do 20 k=km,km/2+1,-1
        k1=k
        k2=k-1
        k3=k-2
        p2=pre(izone,i,k1)
        p4=pre(izone,i,k2)
        p6=pre(izone,i,k3)	  
        vbac=dabs(p2-2.*p4+p6)/(dabs(p2)+2.*dabs(p4)+dabs(p6))
        d31=wp1(izone,i,k3)-wp1(izone,i,k2)
        d32=wp2(izone,i,k3)-wp2(izone,i,k2)
        d33=wp3(izone,i,k3)-wp3(izone,i,k2)
        d34=wp4(izone,i,k3)-wp4(izone,i,k2)
        volm=.5*(dtk(izone,i,k1)+dtk(izone,i,k2))
        ep2p=ek2*dmax1(vbac,vfor)
        ep4p=dmax1(dble(0.),(ek4-ep2p))
        dx1p=ep2p*d21
        dx2p=ep2p*d22
        dx3p=ep2p*d23
        dx4p=ep2p*d24
        dx1p=(dx1p-ep4p*(d31-2.*d21+d11))*volm
        dx2p=(dx2p-ep4p*(d32-2.*d22+d12))*volm
        dx3p=(dx3p-ep4p*(d33-2.*d23+d13))*volm
        dx4p=(dx4p-ep4p*(d34-2.*d24+d14))*volm
        diss(izone,i,k,1)=dx1p-dx1m
        diss(izone,i,k,2)=dx2p-dx2m
        diss(izone,i,k,3)=dx3p-dx3m
        diss(izone,i,k,4)=dx4p-dx4m
        vfor=vbac
        d11=d21
        d12=d22
        d13=d23
        d14=d24
        d21=d31
        d22=d32
        d23=d33
        d24=d34
        dx1m=dx1p
        dx2m=dx2p
        dx3m=dx3p
        dx4m=dx4p
20      continue


!       i-direction
 
        do 30 k=1,km
!   	  (i=1)
        if(kb(izone,3,k).eq.0)then
        dx1m=0.
        dx2m=0.
        dx3m=0.
        dx4m=0.
        w1=wp1(izone,0,k)
        w2=wp2(izone,0,k)
        w3=wp3(izone,0,k)
        w4=wp4(izone,0,k)
        call trans(w1,w2,w3,w4,d,u,w,pssp)
        p2=pssp
        p4=pre(izone,1,k)
        p6=pre(izone,2,k)
        vfor=dabs(p2-2.*p4+p6)/(dabs(p2)+2.*dabs(p4)+dabs(p6))
        d21=wp1(izone,2,k)-wp1(izone,1,k)
        d22=wp2(izone,2,k)-wp2(izone,1,k)
        d23=wp3(izone,2,k)-wp3(izone,1,k)
        d24=wp4(izone,2,k)-wp4(izone,1,k)
        d11=d21
        d12=d22
        d13=d23
        d14=d24
        else if(kb(izone,3,k).eq.1)then
        dx1m=0.
        dx2m=0.
        dx3m=0.
        dx4m=0.
        if(kns.eq.1) then
        p2=pre(izone,1,k)
        else
        p2=2*pre(izone,1,k)-pre(izone,2,k)
        endif
        p4=pre(izone,1,k)        
        p6=pre(izone,2,k)
        vfor=dabs(p2-2.*p4+p6)/(dabs(p2)+2.*dabs(p4)+dabs(p6))
        d21=wp1(izone,2,k)-wp1(izone,1,k)
        d22=wp2(izone,2,k)-wp2(izone,1,k)
        d23=wp3(izone,2,k)-wp3(izone,1,k)
        d24=wp4(izone,2,k)-wp4(izone,1,k)
        if(kns.eq.1) then
        d11=0.d0
        d12=2*wp2(izone,1,k)
        d13=2*wp3(izone,1,k)
        d14=0.d0
        else
        d11=d21
        d12=d22
        d13=d23
        d14=d24
        endif
        else
        ia1=md(izone,k,1,1,3)
        i11=md(izone,k,2,1,3)
        i21=md(izone,k,3,1,3)
        ibh=md(izone,k,4,1,3)
        ia2=md(izone,k,1,2,3)
        i12=md(izone,k,2,2,3)
        i22=md(izone,k,3,2,3)

        p2=pre(ia2,i12,i22)
        p4=pre(ia1,i11,i21)
        p6=pre(izone,1,k)
        vfor=dabs(p2-2.*p4+p6)/(dabs(p2)+2.*dabs(p4)+dabs(p6))
        p2=pre(ia1,i11,i21)
        p4=pre(izone,1,k)
        p6=pre(izone,2,k)
        vbac=dabs(p2-2.*p4+p6)/(dabs(p2)+2.*dabs(p4)+dabs(p6))
        d31=wp1(izone,2,k)-wp1(izone,1,k)
        d32=wp2(izone,2,k)-wp2(izone,1,k)
        d33=wp3(izone,2,k)-wp3(izone,1,k)
        d34=wp4(izone,2,k)-wp4(izone,1,k)
        d21=wp1(izone,1,k)-wp1(ia1,i11,i21)
        d22=wp2(izone,1,k)-wp2(ia1,i11,i21)
        d23=wp3(izone,1,k)-wp3(ia1,i11,i21)
        d24=wp4(izone,1,k)-wp4(ia1,i11,i21)
        d11=wp1(ia1,i11,i21)-wp1(ia2,i12,i22)
        d12=wp2(ia1,i11,i21)-wp2(ia2,i12,i22)
        d13=wp3(ia1,i11,i21)-wp3(ia2,i12,i22)
        d14=wp4(ia1,i11,i21)-wp4(ia2,i12,i22)
        ep2m=ek2*dmax1(vbac,vfor)
        ep4m=dmax1(dble(0.),(ek4-ep2m))
        if(ibh.eq.1.or.ibh.eq.2)then
        volm=.5*(dti(izone,1,k)+dtk(ia1,i11,i21))
        else
        volm=.5*(dti(izone,1,k)+dti(ia1,i11,i21))
        endif
        dx1m=ep2m*d21
        dx2m=ep2m*d22
        dx3m=ep2m*d23
        dx4m=ep2m*d24
        dx1m=(dx1m-ep4m*(d31-2.*d21+d11))*volm
        dx2m=(dx2m-ep4m*(d32-2.*d22+d12))*volm
        dx3m=(dx3m-ep4m*(d33-2.*d23+d13))*volm
        dx4m=(dx4m-ep4m*(d34-2.*d24+d14))*volm
        d11=d21
        d12=d22
        d13=d23
        d14=d24
        d21=d31
        d22=d32
        d23=d33
        d24=d34
        vfor=vbac
        end if

        do 30 i=1,im/2
        i1=i
        i2=i+1
        i3=i+2

        p2=pre(izone,i1,k)
        p4=pre(izone,i2,k)
        p6=pre(izone,i3,k)
        vbac=dabs(p2-2.*p4+p6)/(dabs(p2)+2.*dabs(p4)+dabs(p6))
        d31=wp1(izone,i3,k)-wp1(izone,i2,k)
        d32=wp2(izone,i3,k)-wp2(izone,i2,k)
        d33=wp3(izone,i3,k)-wp3(izone,i2,k)
        d34=wp4(izone,i3,k)-wp4(izone,i2,k)
        volm=.5*(dti(izone,i1,k)+dti(izone,i2,k))

        ep2p=ek2*dmax1(vbac,vfor)
        ep4p=dmax1(dble(0.),(ek4-ep2p))
        dx1p=ep2p*d21
        dx2p=ep2p*d22
        dx3p=ep2p*d23
        dx4p=ep2p*d24
        dx1p=(dx1p-ep4p*(d31-2.*d21+d11))*volm
        dx2p=(dx2p-ep4p*(d32-2.*d22+d12))*volm
        dx3p=(dx3p-ep4p*(d33-2.*d23+d13))*volm
        dx4p=(dx4p-ep4p*(d34-2.*d24+d14))*volm
        diss(izone,i,k,1)=dx1p-dx1m+diss(izone,i,k,1)
        diss(izone,i,k,2)=dx2p-dx2m+diss(izone,i,k,2)
        diss(izone,i,k,3)=dx3p-dx3m+diss(izone,i,k,3)
        diss(izone,i,k,4)=dx4p-dx4m+diss(izone,i,k,4)
        vfor=vbac
        d11=d21
        d12=d22
        d13=d23
        d14=d24
        d21=d31
        d22=d32
        d23=d33
        d24=d34
        dx1m=dx1p
        dx2m=dx2p
        dx3m=dx3p
        dx4m=dx4p
30      continue



!       i-direction

        do 40 k=1,km
!   	  (i=im+1)
        if(kb(izone,4,k).eq.0)then
        dx1m=0.
        dx2m=0.
        dx3m=0.
        dx4m=0.
        w1=wp1(izone,im+1,k)
        w2=wp2(izone,im+1,k)
        w3=wp3(izone,im+1,k)
        w4=wp4(izone,im+1,k)
        call trans(w1,w2,w3,w4,d,u,w,pssp)
        p2=pssp
        p4=pre(izone,im,k)
        p6=pre(izone,im-1,k)
        vfor=dabs(p2-2.*p4+p6)/(dabs(p2)+2.*dabs(p4)+dabs(p6))
        d21=wp1(izone,im-1,k)-wp1(izone,im,k)
        d22=wp2(izone,im-1,k)-wp2(izone,im,k)
        d23=wp3(izone,im-1,k)-wp3(izone,im,k)
        d24=wp4(izone,im-1,k)-wp4(izone,im,k)
        d11=d21
        d12=d22
        d13=d23
        d14=d24
        else if(kb(izone,4,k).eq.1)then
        dx1m=0.
        dx2m=0.
        dx3m=0.
        dx4m=0.
        if(kns.eq.1) then
        p2=pre(izone,im,k)
        else
        p2=2*pre(izone,im,k)-pre(izone,im-1,k)
        endif
        p4=pre(izone,im,k)
        p6=pre(izone,im-1,k)
        vfor=dabs(p2-2.*p4+p6)/(dabs(p2)+2.*dabs(p4)+dabs(p6))
        d21=wp1(izone,im-1,k)-wp1(izone,im,k)
        d22=wp2(izone,im-1,k)-wp2(izone,im,k)
        d23=wp3(izone,im-1,k)-wp3(izone,im,k)
        d24=wp4(izone,im-1,k)-wp4(izone,im,k)
        if(kns.eq.1) then
        d11=0.d0
        d12=2*wp2(izone,im,k)
        d13=2*wp3(izone,im,k)
        d14=0.d0
        else
        d11=d21
        d12=d22
        d13=d23
        d14=d24
        endif
        else
        ia1=md(izone,k,1,1,4)
        i11=md(izone,k,2,1,4)
        i21=md(izone,k,3,1,4)
        ibh=md(izone,k,4,1,4)
        ia2=md(izone,k,1,2,4)
        i12=md(izone,k,2,2,4)
        i22=md(izone,k,3,2,4)

        p2=pre(ia2,i12,i22)
        p4=pre(ia1,i11,i21)
        p6=pre(izone,im,k)
        vfor=dabs(p2-2.*p4+p6)/(dabs(p2)+2.*dabs(p4)+dabs(p6))
        p2=pre(ia1,i11,i21)
        p4=pre(izone,im,k)
        p6=pre(izone,im-1,k)
        vbac=dabs(p2-2.*p4+p6)/(dabs(p2)+2.*dabs(p4)+dabs(p6))
        d31=wp1(izone,im-1,k)-wp1(izone,im,k)
        d32=wp2(izone,im-1,k)-wp2(izone,im,k)
        d33=wp3(izone,im-1,k)-wp3(izone,im,k)
        d34=wp4(izone,im-1,k)-wp4(izone,im,k)
        d21=wp1(izone,im,k)-wp1(ia1,i11,i21)
        d22=wp2(izone,im,k)-wp2(ia1,i11,i21)
        d23=wp3(izone,im,k)-wp3(ia1,i11,i21)
        d24=wp4(izone,im,k)-wp4(ia1,i11,i21)
        d11=wp1(ia1,i11,i21)-wp1(ia2,i12,i22)
        d12=wp2(ia1,i11,i21)-wp2(ia2,i12,i22)
        d13=wp3(ia1,i11,i21)-wp3(ia2,i12,i22)
        d14=wp4(ia1,i11,i21)-wp4(ia2,i12,i22)
        ep2m=ek2*dmax1(vbac,vfor)
        ep4m=dmax1(dble(0.),(ek4-ep2m))
        if(ibh.eq.1.or.ibh.eq.2)then
        volm=.5*(dti(izone,im,k)+dtk(ia1,i11,i21))
        else
        volm=.5*(dti(izone,im,k)+dti(ia1,i11,i21))
        endif
        dx1m=ep2m*d21
        dx2m=ep2m*d22
        dx3m=ep2m*d23
        dx4m=ep2m*d24
        dx1m=(dx1m-ep4m*(d31-2.*d21+d11))*volm
        dx2m=(dx2m-ep4m*(d32-2.*d22+d12))*volm
        dx3m=(dx3m-ep4m*(d33-2.*d23+d13))*volm
        dx4m=(dx4m-ep4m*(d34-2.*d24+d14))*volm
        d11=d21
        d12=d22
        d13=d23
        d14=d24
        d21=d31
        d22=d32
        d23=d33
        d24=d34
        vfor=vbac
        end if
          
        do 40 i=im,im/2+1,-1
        i1=i
        i2=i-1
        i3=i-2

        p2=pre(izone,i1,k)
        p4=pre(izone,i2,k)
        p6=pre(izone,i3,k)
        vbac=dabs(p2-2.*p4+p6)/(dabs(p2)+2.*dabs(p4)+dabs(p6))
        d31=wp1(izone,i3,k)-wp1(izone,i2,k)
        d32=wp2(izone,i3,k)-wp2(izone,i2,k)
        d33=wp3(izone,i3,k)-wp3(izone,i2,k)
        d34=wp4(izone,i3,k)-wp4(izone,i2,k)
        volm=.5*(dti(izone,i1,k)+dti(izone,i2,k))

        ep2p=ek2*dmax1(vbac,vfor)
        ep4p=dmax1(dble(0.),(ek4-ep2p))
        dx1p=ep2p*d21
        dx2p=ep2p*d22
        dx3p=ep2p*d23
        dx4p=ep2p*d24
        dx1p=(dx1p-ep4p*(d31-2.*d21+d11))*volm
        dx2p=(dx2p-ep4p*(d32-2.*d22+d12))*volm
        dx3p=(dx3p-ep4p*(d33-2.*d23+d13))*volm
        dx4p=(dx4p-ep4p*(d34-2.*d24+d14))*volm

        diss(izone,i,k,1)=dx1p-dx1m+diss(izone,i,k,1)
        diss(izone,i,k,2)=dx2p-dx2m+diss(izone,i,k,2)
        diss(izone,i,k,3)=dx3p-dx3m+diss(izone,i,k,3)
        diss(izone,i,k,4)=dx4p-dx4m+diss(izone,i,k,4)
        vfor=vbac
        d11=d21
        d12=d22
        d13=d23
        d14=d24
        d21=d31
        d22=d32
        d23=d33
        d24=d34
        dx1m=dx1p
        dx2m=dx2p
        dx3m=dx3p
        dx4m=dx4p
40      continue
100	    continue
        end


        subroutine resid
        use dimen
        use grid
        use chara
        use conver
        use resi
        use cfdiss
        use kswitch
        use qinf
        use data0
        use viscous
        use vissflux
        implicit none
        integer                   ::izone,i,k,numb,ix,iz,im,km,ii,kk,ia,i1,i2
        parameter(numb=1,ix=369,iz=65)
        real*8                    ::w1,w2,w3,w4,d,u,w,p,pb,fskf1(ix),fskf2(ix),fskf3(ix),fskf4(ix),&
                                    fsif1,fsif2,fsif3,fsif4,fsib1,fsib2,fsib3,fsib4,fskb1,fskb2,fskb3,&
                                    fskb4,q11,q21,q31,q41,q1,q2,q3,q4

        do 200 izone=1,nmp	
        im=me(izone,1)
        km=me(izone,2)

        do i=1,im

!       k=1

        if(kb(izone,1,i).eq.0) then
        w1=wp1(izone,i,0)
        w2=wp2(izone,i,0)
        w3=wp3(izone,i,0)
        w4=wp4(izone,i,0)

        call trans(w1,w2,w3,w4,d,u,w,p)
        fskf1(i)=w2*ark(izone,i,1,1)+w3*ark(izone,i,1,2)  !对流通量第一项在远场k=1边界的积分
        fskf2(i)=(u*w2+p)*ark(izone,i,1,1)+u*(w3*ark(izone,i,1,2))
        fskf3(i)=w*(w2*ark(izone,i,1,1))+(w*w3+p)*ark(izone,i,1,2)
        fskf4(i)=(w4+p)*(u*ark(izone,i,1,1)+w*ark(izone,i,1,2))

        else if(kb(izone,1,i).eq.1) then
        w1=wp1(izone,i,1)
        w2=wp2(izone,i,1)
        w3=wp3(izone,i,1)
        w4=wp4(izone,i,1)

        call trans(w1,w2,w3,w4,d,u,w,pb)			  
        fskf1(i)=0.  !对流通量第一项在物面k=1边界积分
        fskf2(i)=pb*ark(izone,i,1,1)
        fskf3(i)=pb*ark(izone,i,1,2)
        fskf4(i)=0.

        else if(kb(izone,1,i).eq.2) then
        ia=md(izone,i,1,1,1)
        i1=md(izone,i,2,1,1)
        i2=md(izone,i,3,1,1)

        w1=0.5*(wp1(izone,i,1)+wp1(ia,i1,i2))
        w2=0.5*(wp2(izone,i,1)+wp2(ia,i1,i2))
        w3=0.5*(wp3(izone,i,1)+wp3(ia,i1,i2))
        w4=0.5*(wp4(izone,i,1)+wp4(ia,i1,i2))

        call trans(w1,w2,w3,w4,d,u,w,p)
        fskf1(i)=w2*ark(izone,i,1,1)+w3*ark(izone,i,1,2)   !对流通量第一项在对接k=1边界积分
        fskf2(i)=(u*w2+p)*ark(izone,i,1,1)+u*(w3*ark(izone,i,1,2))
        fskf3(i)=w*(w2*ark(izone,i,1,1))+(w*w3+p)*ark(izone,i,1,2)
        fskf4(i)=(w4+p)*(u*ark(izone,i,1,1)+w*ark(izone,i,1,2))
        endif
        enddo

        do 100 k=1,km

!  	    i=1        

        if(kb(izone,3,k).eq.0) then
        w1=wp1(izone,0,k)
        w2=wp2(izone,0,k)
        w3=wp3(izone,0,k)
        w4=wp4(izone,0,k)

        call trans(w1,w2,w3,w4,d,u,w,p)
        fsif1=w2*ari(izone,1,k,1)+w3*ari(izone,1,k,2)
        fsif2=(u*w2+p)*ari(izone,1,k,1)+u*(w3*ari(izone,1,k,2))
        fsif3=w*(w2*ari(izone,1,k,1))+(w*w3+p)*ari(izone,1,k,2)
        fsif4=(w4+p)*(u*ari(izone,1,k,1)+w*ari(izone,1,k,2))

        else if(kb(izone,3,k).eq.1) then
        w1=wp1(izone,1,k)
        w2=wp2(izone,1,k)
        w3=wp3(izone,1,k)
        w4=wp4(izone,1,k)

        call trans(w1,w2,w3,w4,d,u,w,pb)
        fsif1=0.
        fsif2=pb*ari(izone,1,k,1)
        fsif3=pb*ari(izone,1,k,2)
        fsif4=0.

        else if(kb(izone,3,k).eq.2) then 
        ia=md(izone,k,1,1,3)
        i1=md(izone,k,2,1,3)
        i2=md(izone,k,3,1,3)

        w1=0.5*(wp1(izone,1,k)+wp1(ia,i1,i2))
        w2=0.5*(wp2(izone,1,k)+wp2(ia,i1,i2))
        w3=0.5*(wp3(izone,1,k)+wp3(ia,i1,i2))
        w4=0.5*(wp4(izone,1,k)+wp4(ia,i1,i2))

        call trans(w1,w2,w3,w4,d,u,w,p)
        fsif1=w2*ari(izone,1,k,1)+w3*ari(izone,1,k,2)
        fsif2=(u*w2+p)*ari(izone,1,k,1)+u*(w3*ari(izone,1,k,2))
        fsif3=w*(w2*ari(izone,1,k,1))+(w*w3+p)*ari(izone,1,k,2)
        fsif4=(w4+p)*(u*ari(izone,1,k,1)+w*ari(izone,1,k,2))
        endif

        do 50 i=1,im

!	  I direction

        if(i.eq.im) then
        ii=im+1
        if(kb(izone,4,k).eq.0) then
	        w1=wp1(izone,ii,k)		
	        w2=wp2(izone,ii,k)		
	        w3=wp3(izone,ii,k)		
	        w4=wp4(izone,ii,k)
        		
        else if(kb(izone,4,k).eq.1) then
            w1=wp1(izone,im,k)
            w2=wp2(izone,im,k)
            w3=wp3(izone,im,k)
            w4=wp4(izone,im,k)

            call trans(w1,w2,w3,w4,d,u,w,pb)
	        fsib1=0.
	        fsib2=pb*ari(izone,ii,k,1)
	        fsib3=pb*ari(izone,ii,k,2)
	        fsib4=0.0
	        goto 47
        else if(kb(izone,4,k).eq.2) then
	        ia=md(izone,k,1,1,4)
	        i1=md(izone,k,2,1,4)
	        i2=md(izone,k,3,1,4)
        	
	        w1=0.5*(wp1(izone,im,k)+wp1(ia,i1,i2))
	        w2=0.5*(wp2(izone,im,k)+wp2(ia,i1,i2))
	        w3=0.5*(wp3(izone,im,k)+wp3(ia,i1,i2))
	        w4=0.5*(wp4(izone,im,k)+wp4(ia,i1,i2))
        endif	  
        goto 46 
        endif

        w1=.5*(wp1(izone,i,k)+wp1(izone,i+1,k))
        w2=.5*(wp2(izone,i,k)+wp2(izone,i+1,k))
        w3=.5*(wp3(izone,i,k)+wp3(izone,i+1,k))
        w4=.5*(wp4(izone,i,k)+wp4(izone,i+1,k))

46	    continue	  

        call trans(w1,w2,w3,w4,d,u,w,p)
        fsib1=w2*ari(izone,i+1,k,1)+w3*ari(izone,i+1,k,2)
        fsib2=(u*w2+p)*ari(izone,i+1,k,1)+u*(w3*ari(izone,i+1,k,2))
        fsib3=w*(w2*ari(izone,i+1,k,1))+(w*w3+p)*ari(izone,i+1,k,2)
        fsib4=(w4+p)*(u*ari(izone,i+1,k,1)+w*ari(izone,i+1,k,2))
47	    continue

        q11=fsif1-fsib1      !沿i向流入单元i的通量（始终定义对于单元i，i+1边界（fsib）流出为正，i边界（fsif）流入为正）
        q21=fsif2-fsib2
        q31=fsif3-fsib3
        q41=fsif4-fsib4

        q1=q11+diss(izone,i,k,1)
        q2=q21+diss(izone,i,k,2)
        q3=q31+diss(izone,i,k,3)
        q4=q41+diss(izone,i,k,4)

        fsif1=fsib1
        fsif2=fsib2
        fsif3=fsib3
        fsif4=fsib4

!
!       K direction
!  
        if(k.eq.km) then
        kk=km+1	  	
	        if(kb(izone,2,i).eq.0) then
	        w1=wp1(izone,i,kk)
	        w2=wp2(izone,i,kk)
	        w3=wp3(izone,i,kk)
	        w4=wp4(izone,i,kk)

        else if(kb(izone,2,i).eq.1) then
            w1=wp1(izone,i,km)
            w2=wp2(izone,i,km)
            w3=wp3(izone,i,km)
            w4=wp4(izone,i,km)

            call trans(w1,w2,w3,w4,d,u,w,pb)	
	        fskb1=0.0
	        fskb2=pb*ark(izone,i,kk,1)
	        fskb3=pb*ark(izone,i,kk,2)
	        fskb4=0.0
	        goto 49
        else if(kb(izone,2,i).eq.2) then
	        ia=md(izone,i,1,1,2)
	        i1=md(izone,i,2,1,2)
	        i2=md(izone,i,3,1,2)

	        w1=0.5*(wp1(izone,i,km)+wp1(ia,i1,i2))
	        w2=0.5*(wp2(izone,i,km)+wp2(ia,i1,i2))
	        w3=0.5*(wp3(izone,i,km)+wp3(ia,i1,i2))
	        w4=0.5*(wp4(izone,i,km)+wp4(ia,i1,i2))
        endif
        goto 48
        endif

        w1=.5*(wp1(izone,i,k)+wp1(izone,i,k+1))
        w2=.5*(wp2(izone,i,k)+wp2(izone,i,k+1))
        w3=.5*(wp3(izone,i,k)+wp3(izone,i,k+1))
        w4=.5*(wp4(izone,i,k)+wp4(izone,i,k+1))

48      continue

        call trans(w1,w2,w3,w4,d,u,w,p)
        fskb1=w2*ark(izone,i,k+1,1)+w3*ark(izone,i,k+1,2)
        fskb2=(u*w2+p)*ark(izone,i,k+1,1)+u*(w3*ark(izone,i,k+1,2))
        fskb3=w*(w2*ark(izone,i,k+1,1))+(w*w3+p)*ark(izone,i,k+1,2)
        fskb4=(w4+p)*(u*ark(izone,i,k+1,1)+w*ark(izone,i,k+1,2))      
49	    continue

        if(kns.eq.0) then
        fvis(izone,i,k,1)=0
        fvis(izone,i,k,2)=0
        fvis(izone,i,k,3)=0
        fvis(izone,i,k,4)=0
        endif
        wup1(izone,i,k)=q1+fskf1(i)-fskb1+fvis(izone,i,k,1)      !残差第一项
        wup2(izone,i,k)=q2+fskf2(i)-fskb2+fvis(izone,i,k,2)
        wup3(izone,i,k)=q3+fskf3(i)-fskb3+fvis(izone,i,k,3)
        wup4(izone,i,k)=q4+fskf4(i)-fskb4+fvis(izone,i,k,4)
        fskf1(i)=fskb1
        fskf2(i)=fskb2
        fskf3(i)=fskb3
        fskf4(i)=fskb4
50      continue
100     continue
200	    continue
        return
        end
        
        
        subroutine ravimplicit(rcfl)
        use dimen
        use grid
        use resi
        implicit none
        integer                   ::ii,i,k,numb,ix,iz,im,km
        parameter(numb=1,ix=369,iz=65)
        real*8                    ::gmki(numb,ix,iz),gmik(numb,ix,iz),gm,epi,epk,rcfl,c1
        real*8                    ::a(ix),b(ix),c(ix),d1(ix),d2(ix),d3(ix),d4(ix),d5(ix)
        common /ggmm/gmki,gmik

	    do 200 ii=1,nmp
		im=me(ii,1)
		km=me(ii,2)

!         i-direction

		do 50 k=1,km
			do i=1,im
				gm=gmki(ii,i,k)
				epi=max(0.25*((rcfl/(1.+0.125*gm))**2-1.),0d0)    !2D时系数取0.125
				a(i)=-epi
				b(i)=2.*epi+1.
				c(i)=-epi
				d1(i)=wup1(ii,i,k)
				d2(i)=wup2(ii,i,k)
				d3(i)=wup3(ii,i,k)
				d4(i)=wup4(ii,i,k)
            enddo
!化成上三角矩阵			
			do i=2,im
				c1=a(i)/b(i-1)
				b(i)=b(i)-c(i-1)*c1
				d1(i)=d1(i)-d1(i-1)*c1
				d2(i)=d2(i)-d2(i-1)*c1
				d3(i)=d3(i)-d3(i-1)*c1
				d4(i)=d4(i)-d4(i-1)*c1
            enddo
			
			wup1(ii,im,k)=d1(im)/b(im)
			wup2(ii,im,k)=d2(im)/b(im)
			wup3(ii,im,k)=d3(im)/b(im)
			wup4(ii,im,k)=d4(im)/b(im)
!化成对角矩阵后求解        
			do i=im-1,1,-1
				d1(i)=d1(i)-d1(i+1)/b(i+1)*c(i)
				d2(i)=d2(i)-d2(i+1)/b(i+1)*c(i)
				d3(i)=d3(i)-d3(i+1)/b(i+1)*c(i)
				d4(i)=d4(i)-d4(i+1)/b(i+1)*c(i)
				d5(i)=d5(i)-d5(i+1)/b(i+1)*c(i)
				wup1(ii,i,k)=d1(i)/b(i)
				wup2(ii,i,k)=d2(i)/b(i)
				wup3(ii,i,k)=d3(i)/b(i)
				wup4(ii,i,k)=d4(i)/b(i)
            enddo
50   		continue

!         k-direction

		do 150 i=1,im
       		do k=1,km
				gm=gmik(ii,i,k)
				epk=max(0.25*((rcfl/(1+0.125*gm))**2-1),0d0)
				a(k)=-epk
				b(k)=2.*epk+1.
				c(k)=-epk
				d1(k)=wup1(ii,i,k)
				d2(k)=wup2(ii,i,k)
				d3(k)=wup3(ii,i,k)
				d4(k)=wup4(ii,i,k)
                enddo
			
			do k=2,km
				c1=a(k)/b(k-1)
				b(k)=b(k)-c(k-1)*c1
				d1(k)=d1(k)-d1(k-1)*c1
				d2(k)=d2(k)-d2(k-1)*c1
				d3(k)=d3(k)-d3(k-1)*c1
				d4(k)=d4(k)-d4(k-1)*c1
                enddo

	        wup1(ii,i,km)=d1(km)/b(km)
			wup2(ii,i,km)=d2(km)/b(km)
			wup3(ii,i,km)=d3(km)/b(km)
			wup4(ii,i,km)=d4(km)/b(km)

			do k=km-1,1,-1
				d1(k)=d1(k)-d1(k+1)/b(k+1)*c(k)
				d2(k)=d2(k)-d2(k+1)/b(k+1)*c(k)
				d3(k)=d3(k)-d3(k+1)/b(k+1)*c(k)
				d4(k)=d4(k)-d4(k+1)/b(k+1)*c(k)
				wup1(ii,i,k)=d1(k)/b(k)
				wup2(ii,i,k)=d2(k)/b(k)
				wup3(ii,i,k)=d3(k)/b(k)
				wup4(ii,i,k)=d4(k)/b(k)
                enddo
150		continue
200     continue
        return
        end
        
        
        
        !远场修正也应包括湍流方程控制变量rmub的边界修正(unfinished)
        subroutine farcor
        use dimen
		use grid
        use chara
        use conver
        use turb
        use safun
        implicit none
        integer                 ::izone,i,k,numb,ix,iz,im,km,kfario
        parameter(numb=1,ix=369,iz=65)
        real*8                  ::wpr1,wpr2,wpr3,wpr4,fx0,fz0
        real*8                  ::xnbk(numb,ix,iz),znbk(numb,ix,iz),xnbi(numb,ix,iz),znbi(numb,ix,iz)
        common /bnor/xnbk,znbk,xnbi,znbi
        common /correct/wpr1,wpr2,wpr3,wpr4,kfario

        do 100 izone=1,nmp
        im=me(izone,1)
        km=me(izone,2)

        do i=1,im
!	    (k=1)
        if(kb(izone,1,i).eq.0) then
		fx0=0.5*(fx(izone,i,1)+fx(izone,i+1,1))
		fz0=0.5*(fz(izone,i,1)+fz(izone,i+1,1))
        call bcal(wp1(izone,i,1),wp2(izone,i,1),wp3(izone,i,1),wp4(izone,i,1),&
                  wp1(izone,i,2),wp2(izone,i,2),wp3(izone,i,2),wp4(izone,i,2),&
       		      xnbk(izone,i,1),znbk(izone,i,1),fx0,fz0)
        wp1(izone,i,0)=wpr1
        wp2(izone,i,0)=wpr2
        wp3(izone,i,0)=wpr3
        wp4(izone,i,0)=wpr4
          if(kfario.eq.0) then
            rmub(izone,i,0)=ru0t
          else
            rmub(izone,i,0)=1.5*rmub(izone,i,1)-0.5*rmub(izone,i,2)
          endif
        endif
!		(k=kmax)
        if(kb(izone,2,i).eq.0) then
		fx0=0.5*(fx(izone,i,km+1)+fx(izone,i+1,km+1))
		fz0=0.5*(fz(izone,i,km+1)+fz(izone,i+1,km+1))
        call bcal(wp1(izone,i,km),wp2(izone,i,km),wp3(izone,i,km),wp4(izone,i,km),&
        	      wp1(izone,i,km-1),wp2(izone,i,km-1),wp3(izone,i,km-1),wp4(izone,i,km-1),&
                  xnbk(izone,i,2),znbk(izone,i,2),fx0,fz0)
        wp1(izone,i,km+1)=wpr1
        wp2(izone,i,km+1)=wpr2
        wp3(izone,i,km+1)=wpr3
        wp4(izone,i,km+1)=wpr4
        endif
          if(kfario.eq.0) then
            rmub(izone,i,km+1)=ru0t
          else
            rmub(izone,i,km+1)=1.5*rmub(izone,i,km)-0.5*rmub(izone,i,km-1)
          endif
        enddo

        do k=1,km
!		(i=1)
        if(kb(izone,3,k).eq.0) then
		fx0=0.5*(fx(izone,1,k)+fx(izone,1,k+1))
		fz0=0.5*(fz(izone,1,k)+fz(izone,1,k+1))
        call bcal(wp1(izone,1,k),wp2(izone,1,k),wp3(izone,1,k),wp4(izone,1,k),&
       			  wp1(izone,2,k),wp2(izone,2,k),wp3(izone,2,k),wp4(izone,2,k),&
                  xnbi(izone,k,1),znbi(izone,k,1),fx0,fz0)
        wp1(izone,0,k)=wpr1
        wp2(izone,0,k)=wpr2
        wp3(izone,0,k)=wpr3
        wp4(izone,0,k)=wpr4
          if(kfario.eq.0) then
            rmub(izone,0,k)=ru0t
          else
            rmub(izone,0,k)=1.5*rmub(izone,1,k)-0.5*rmub(izone,2,k)
          endif
        endif
!		(i=imax)
        if(kb(izone,4,k).eq.0) then
		fx0=0.5*(fx(izone,im+1,k)+fx(izone,im+1,k+1))
		fz0=0.5*(fz(izone,im+1,k)+fz(izone,im+1,k+1))
        call bcal(wp1(izone,im,k),wp2(izone,im,k),wp3(izone,im,k),wp4(izone,im,k),&
       		      wp1(izone,im-1,k),wp2(izone,im-1,k),wp3(izone,im-1,k),wp4(izone,im-1,k),&
                  xnbi(izone,k,2),znbi(izone,k,2),fx0,fz0)
        wp1(izone,im+1,k)=wpr1
        wp2(izone,im+1,k)=wpr2
        wp3(izone,im+1,k)=wpr3
        wp4(izone,im+1,k)=wpr4
          if(kfario.eq.0) then
            rmub(izone,im+1,k)=ru0t
          else
            rmub(izone,im+1,k)=1.5*rmub(izone,im,k)-0.5*rmub(izone,im-1,k)
          endif
        endif
        enddo
100   continue
        return
        end
        
        
        subroutine bcal(wpout1,wpout2,wpout3,wpout4,wpin1,wpin2,wpin3,wpin4,xnb,znb,fx0,fz0)
!       .....  correction for outer boundary 
        use dimen
        use data0
        use qinf
        implicit none
        integer                 ::numb,ix,iz,kfario
        parameter(numb=1,ix=369,iz=65)
        real*8                  ::wpout1,wpout2,wpout3,wpout4,wpin1,wpin2,wpin3,wpin4,xnb,znb,&
                                  xnbk(numb,ix,iz),znbk(numb,ix,iz),xnbi(numb,ix,iz),znbi(numb,ix,iz),&
                                  wpr1,wpr2,wpr3,wpr4,sif,dint,uint,wint,pint,vn,vnif,rif,w1,w2,w3,w4,&
                                  d,u,w,p,s,vne,a2,rext,a,fx0,fz0,uif0,wif0,amif0

        common /bnor/xnbk,znbk,xnbi,znbi
        common /correct/wpr1,wpr2,wpr3,wpr4,kfario
 		common /vortexcor/uif0,wif0,amif0
		       
        sif=pif/(dif**gma)
        call trans(wpout1,wpout2,wpout3,wpout4,dint,uint,wint,pint)
        vn=uint*xnb+wint*znb   !边界法向速度
        if(dabs(vn)/amif.lt.1.) goto 10 

!       for surpersonic inflow or outflow

        if(vn.lt.0.) then
!       in-flow  !全部为自由来流值
        kfario=0   !入流/出流的标志
        wpr1=dif
        wpr2=dif*uif
        wpr3=dif*wif
        wpr4=eif
        go to 100    

        else
!       out-flow  !全部由内场外插
        kfario=1
        wpr1=1.5*wpout1-.5*wpin1
        wpr2=1.5*wpout2-.5*wpin2
        wpr3=1.5*wpout3-.5*wpin3
        wpr4=1.5*wpout4-.5*wpin4
        go to 100
        endif

10      continue

!       for subsonic in- or out-flow
		
		call vorcor(fx0,fz0)      
        vnif=uif0*xnb+wif0*znb
        rif=vnif-2.*amif0/(gma-1.)   !r-
        w1=1.5*wpout1-.5*wpin1
        w2=1.5*wpout2-.5*wpin2
        w3=1.5*wpout3-.5*wpin3
        w4=1.5*wpout4-.5*wpin4
        call trans(w1,w2,w3,w4,d,u,w,p)
        vne=u*xnb+w*znb
        a2=gma*p/d
        rext=vne+2.*sqrt(a2)/(gma-1.)  !r+
        vn=(rext+rif)*.5            !远场边界法向速度
        a=(gma-1.)*.25*(rext-rif)  !远场边界声速
        !#        if(a.lt.(.05*amif)) a=.01*amif 
        
        if(vn.lt.0.) then
!       in-flow
        kfario=0
        s=sif
        d=(a*a/(gma*s))**(1./(gma-1.))
        p=d*a*a/gma
        wpr1=d
        u=uif0+(vn-vnif)*xnb
        w=wif0+(vn-vnif)*znb
        wpr2=u*d
        wpr3=w*d
        wpr4=p/(gma-1.)+.5*d*(u*u+w*w)
        go to 100

        else
!       out-flow
        kfario=1
        s=p/d**gma
        d=(a*a/(gma*s))**(1./(gma-1.))
        p=d*a*a/gma
        wpr1=d
        u=u+(vn-vne)*xnb
        w=w+(vn-vne)*znb
        wpr2=u*d
        wpr3=w*d
        wpr4=p/(gma-1.)+.5*d*(u*u+w*w)
        endif
100     continue
        return
        end
        
		subroutine vorcor(fx0,fz0)
		use dimen
		use grid
		use data0
		use qinf
		implicit none
		real*8                 :: fx0,fz0,d,theta,vor,gma0,fl,uif0,wif0,pif0,dif0,amif0,pi,xref,zref   
		common /vortexcor/uif0,wif0,amif0
		common /lift/fl
        xref=0.25
		zref=0.d0
		d=sqrt((fx0-xref)**2+(fz0-zref)**2)
		theta=atan2((fz0-zref),(fx0-xref))
		vor=0.5*sqrt(uif**2+wif**2)*fl
		pi=atan(1.0d0)*4.0
		uif0=uif+(vor*sqrt(1-amach**2)/(2*pi*d))*sin(theta)/&
				 (1-amach**2*(sin(theta-alf))**2)
		wif0=wif-(vor*sqrt(1-amach**2)/(2*pi*d))*cos(theta)/&
				 (1-amach**2*(sin(theta-alf))**2)
		gma0=gma/(gma-1)
		pif0=(pif**(1/gma0)+1/gma0*dif*(uif**2+wif**2-uif0**2-wif0**2)/(2*pif**(1/gma)))**gma0  
		dif0=dif*(pif0/pif)**(1/gma)
		amif0=sqrt(gma*pif0/dif0)
		return
		end subroutine
				  
        
        subroutine enthal(coenth)  !焓阻尼程序仅用于能量方程
        use dimen
        use data0
        use qinf
        use conver
        implicit none
        integer                 ::ii,i,k,numb,ix,iz
        parameter(numb=1,ix=369,iz=65)
        real*8                  ::hif,d,u,w,p,coenth

        hif=(eif+pif)/dif
        do ii=1,nmp 
        do k=1,me(ii,2)
        do i=1,me(ii,1)
        call trans(wp1(ii,i,k),wp2(ii,i,k),wp3(ii,i,k),wp4(ii,i,k),d,u,w,p)
        wp4(ii,i,k)=(wp4(ii,i,k)+coenth*(d*hif-p))/(1.+coenth)
        enddo
        enddo
        enddo
        return
        end
        
        
        
        subroutine output(icycle,res)
        use dimen
        use chara
        use grid
        use file
        use data0
        use qinf
        use conver
        use turb
        use kswitch
        implicit none
        integer                 ::ii,i,k,numb,ix,iz,im,km,test,icycle
        parameter(numb=1,ix=369,iz=65)
        real*8                  ::dp,cl,cdp,cm,w1,w2,w3,w4,d,u,w,p,t,c11,c12,c21,c22,&
                                  ul0,deltu,txx,txz,tzz,fl,fdp,fm,ldratio,res,cp
        real*8                  ::x0g(numb,ix,iz),z0g(numb,ix,iz)
        real*8                  ::ui1(numb,ix+1,iz),ui2(numb,ix+1,iz),wi1(numb,ix+1,iz),wi2(numb,ix+1,iz),&
                                    ti1(numb,ix+1,iz),ti2(numb,ix+1,iz),uk1(numb,ix,iz+1),uk2(numb,ix,iz+1),&
                                    wK1(numb,ix,iz+1),wK2(numb,ix,iz+1),tk1(numb,ix,iz+1),tk2(numb,ix,iz+1)        
        common /grad/ui1,ui2,wi1,wi2,ti1,ti2,uk1,uk2,wk1,wk2,tk1,tk2 
        common /xz0g/x0g,z0g
		common /lift/fl

        dp=0.7*amach**2

        cl=0.
        cdp=0.
        cm=0.

        do 200 ii=1,nmp
        im=me(ii,1)
        km=me(ii,2)

!       k=1
	
        do i=1,me(ii,1)		
        if(kb(ii,1,i).eq.1)then  !k=1边界
          w1=wp1(ii,i,1)
          w2=wp2(ii,i,1)
          w3=wp3(ii,i,1)
          w4=wp4(ii,i,1)
          call trans(w1,w2,w3,w4,d,u,w,p)
          c11=(p-1.)/dp*ark(ii,i,1,2)    !压差引起的升力(压强在边界的y负向积分)
          c12=(p-1.)/dp*ark(ii,i,1,1)    !压差阻力(压强在边界的x负向积分)
          c21=0.
          c22=0.

          if(kns.eq.1) then
              t=p/d
	        ul0=(t**(1.5)*(1+cc)/(t+cc))*sqrt(gma)*amach/ren   !粘性系数 物面湍流粘性系数为0
            deltu=2/3.*(uk1(ii,i,1)+wk2(ii,i,1))
            txx=ul0*(2*uk1(ii,i,1)-deltu)
            tzz=ul0*(2*wk2(ii,i,1)-deltu)
            txz=ul0*(uk2(ii,i,1)+wk1(ii,i,1))
            c21=txz*ark(ii,i,1,1)+tzz*ark(ii,i,1,2)   !粘性引起的升力
            c22=txx*ark(ii,i,1,1)+txz*ark(ii,i,1,2)   !粘性阻力
            c21=-c21/dp
            c22=-c22/dp
          endif			
          cl=cl+c21+c11
          cdp=cdp+c22+c12
          cm=cm+(c11+c21)*((fx(ii,i,1)+fx(ii,i+1,1))/2.-0.25)&
                -(c12+c22)*(fz(ii,i,1)+fz(ii,i+1,1))/2.
        endif
        enddo


!        k=kmax

        do i=1,me(ii,1)		
        if(kb(ii,2,i).eq.1)then  !k=1边界
          w1=wp1(ii,i,km)
          w2=wp2(ii,i,km)
          w3=wp3(ii,i,km)
          w4=wp4(ii,i,km)
          call trans(w1,w2,w3,w4,d,u,w,p)
          c11=-(p-1.)/dp*ark(ii,i,km+1,2)    !压差引起的升力(压强在边界的y负向积分)
          c12=-(p-1.)/dp*ark(ii,i,km+1,1)    !压差阻力(压强在边界的x正向积分)
          c21=0.
          c22=0.
          if(kns.eq.1) then
            t=p/d
	        ul0=(t**(1.5)*(1+cc)/(t+cc))*sqrt(gma)*amach/ren   !粘性系数 物面湍流粘性系数为0
            deltu=2/3.*(uk1(ii,km+1,1)+wk2(ii,km+1,1))
            txx=ul0*(2*uk1(ii,km+1,1)-deltu)
            tzz=ul0*(2*wk2(ii,km+1,1)-deltu)
            txz=ul0*(uk2(ii,km+1,1)+wk1(ii,km+1,1))
            c21=txz*ark(ii,i,km+1,1)+tzz*ark(ii,i,km+1,2)   !粘性引起的升力
            c22=txx*ark(ii,i,km+1,1)+txz*ark(ii,i,km+1,2)   !粘性阻力
            c21=c21/dp
            c22=c22/dp
          endif			
          cl=cl+c21+c11
          cdp=cdp+c22+c12
          cm=cm+(c11+c21)*((fx(ii,i,km+1)+fx(ii,i+1,km+1))/2.-0.25)&
      		   -(c12+c22)*(fz(ii,i,km+1)+fz(ii,i+1,km+1))/2.
        endif
        enddo


!        i=1
	
        do k=1,me(ii,2)		
        if(kb(ii,3,k).eq.1)then  !k=1边界
          w1=wp1(ii,1,k)
          w2=wp2(ii,1,k)
          w3=wp3(ii,1,k)
          w4=wp4(ii,1,k)
          call trans(w1,w2,w3,w4,d,u,w,p)
          c11=(p-1.)/dp*ari(ii,1,k,2)    !压差引起的升力(压强在边界的y负向积分)
          c12=(p-1.)/dp*ari(ii,1,k,1)    !压差阻力(压强在边界的x正向积分)
          c21=0.
          c22=0.
          if(kns.eq.1) then
            t=p/d
	        ul0=(t**(1.5)*(1+cc)/(t+cc))*sqrt(gma)*amach/ren   !粘性系数 物面湍流粘性系数为0
            deltu=2/3.*(ui1(ii,1,k)+wi2(ii,1,k))
            txx=ul0*(2*ui1(ii,1,k)-deltu)
            tzz=ul0*(2*wi2(ii,1,k)-deltu)
            txz=ul0*(ui2(ii,1,k)+wi1(ii,1,k))
            c21=txz*ari(ii,1,k,1)+tzz*ari(ii,1,k,2)   !粘性引起的升力
            c22=txx*ari(ii,1,k,1)+txz*ari(ii,1,k,2)   !粘性阻力
            c21=-c21/dp
            c22=-c22/dp
          endif			
          cl=cl+c21+c11
          cdp=cdp+c22+c12
          cm=cm+(c11+c21)*((fx(ii,1,k)+fx(ii,1,k+1))/2.-0.25)&
             -(c12+c22)*(fz(ii,1,k)+fz(ii,1,k+1))/2.
        endif
        enddo


!        i=imax
	
        do k=1,me(ii,2)		
        if(kb(ii,4,k).eq.1)then  !k=1边界
          w1=wp1(ii,im,k)
          w2=wp2(ii,im,k)
          w3=wp3(ii,im,k)
          w4=wp4(ii,im,k)
          call trans(w1,w2,w3,w4,d,u,w,p)
          c11=(p-1.)/dp*ari(ii,im+1,k,2)    !压差引起的升力(压强在边界的y负向积分)
          c12=(p-1.)/dp*ari(ii,im+1,k,1)    !压差阻力(压强在边界的x正向积分)
          c21=0.
          c22=0.
          if(kns.eq.1) then
            t=p/d
	        ul0=(t**(1.5)*(1+cc)/(t+cc))*sqrt(gma)*amach/ren   !粘性系数 物面湍流粘性系数为0
            deltu=2/3.*(ui1(ii,im+1,k)+wi2(ii,im+1,k))
            txx=ul0*(2*ui1(ii,im+1,k)-deltu)
            tzz=ul0*(2*wi2(ii,im+1,k)-deltu)
            txz=ul0*(ui2(ii,im+1,k)+wi1(ii,im+1,k))
            c21=txz*ari(ii,im+1,k,1)+tzz*ari(ii,im+1,k,2)   !粘性引起的升力
            c22=txx*ari(ii,im+1,k,1)+txz*ari(ii,im+1,k,2)   !粘性阻力
            c21=-c21/dp
            c22=-c22/dp
          endif			
          cl=cl+c21+c11
          cdp=cdp+c22+c12
          cm=cm+(c11+c21)*((fx(ii,im+1,k)+fx(ii,im+1,k+1))/2.-0.25)&
                -(c12+c22)*(fz(ii,im+1,k)+fz(ii,im+1,k+1))/2.
        endif
        enddo

200	continue
        fl=-cl*cos(alf)+cdp*sin(alf)
        fdp=-cl*sin(alf)-cdp*cos(alf)
        fm=cm
        write(*,15)fl,fdp,fm
15	    format(1x,"cl=",f11.6,3x,"cd=",f11.6,3x,"cm=",f11.6)

        open(4,access="append",file=fn4)
        ldratio=fl/fdp
        write(4,20) icycle,fl,fdp,fm,res,ldratio
20	    format(i7,2x,5(d14.6,2x))
        close(4)

        open (55,file=fn5)
        do 30 ii=1,nmp
        write(55,*)ii
        do i=1,me(ii,1)
	        if(kb(ii,1,i).eq.1)then
		        w1=wp1(ii,i,1)
		        w2=wp2(ii,i,1)
		        w3=wp3(ii,i,1)
		        w4=wp4(ii,i,1)
		        call trans(w1,w2,w3,w4,d,u,w,p)
		        cp=(p-1.)/dp
		        write(55,*) x0g(ii,i,1),cp,d
	        endif
        enddo
30	    continue
        close(55)

!       test

        test=0
        if(test.eq.1) then
        open(17,file="test.dat")
        do ii=1,nmp
        do i=1,me(ii,1)
        if(kb(ii,1,i).eq.1) then
		        w1=wp1(ii,i,1)
		        w2=wp2(ii,i,1)
		        w3=wp3(ii,i,1)
		        w4=wp4(ii,i,1)
		        call trans(w1,w2,w3,w4,d,u,w,p)
        write(17,*) d,u,w,p
        endif
        enddo
        enddo
        close(17)
        endif
        call drawing
        return
        end
        
        
        subroutine drawing
        use dimen
        use chara
        use grid
        use file
        use data0
        use qinf
        use conver
        implicit none
        integer                 ::ii,i,k,numb,ix,iz,im,km,imax,kmax
        parameter(numb=1,ix=369,iz=65)
        real*8                  ::w1,w2,w3,w4,d,u,w,p,aa0aa,absv
        real*8                  ::dc(ix,iz),uc(ix,iz),wc(ix,iz),pc(ix,iz),rm(ix,iz)

        open(66,file=fn6)
        WRITE(66,*)'TITLE=flow.plt'
        WRITE(66,*)'VARIABLES="X","Z","D","U","W","P","macha"'
        do ii=1,nmp
        im=me(ii,1)
        km=me(ii,2)	
        IMAX=mm(ii,1)
        KMAX=mm(ii,2)
        do k=1,me(ii,2)
        do i=1,me(ii,1)
        w1=wp1(ii,i,k)
        w2=wp2(ii,i,k)
        w3=wp3(ii,i,k)
        w4=wp4(ii,i,k)
        call trans(w1,w2,w3,w4,d,u,w,p)
        dc(i,k)=d
        uc(i,k)=u
        wc(i,k)=w
        pc(i,k)=p
        aa0aa=sqrt(1.4*p/d)
        absv=sqrt(u*u+w*w)
        rm(i,k)=absv/aa0aa
        enddo
        enddo
        WRITE(66,500)IMAX,KMAX
        WRITE(66,*) 'VARLOCATION=([3-7]=CELLCENTERED)'
        WRITE(66,*)((fx(ii,i,k),i=1,IMAX),k=1,KMAX)
        WRITE(66,*)((fz(ii,i,k),i=1,IMAX),k=1,KMAX)
        WRITE(66,*)((dc(i,k),i=1,im),k=1,km)
        WRITE(66,*)((uc(i,k),i=1,im),k=1,km)
        WRITE(66,*)((wc(i,k),i=1,im),k=1,km)
        WRITE(66,*)((pc(i,k),i=1,im),k=1,km)
        WRITE(66,*)((rm(i,k),i=1,im),k=1,km)
        enddo
        close(66)
500	    format(1x,'ZONE I=',i5,',K=',i5,',F=BLOCK')
        end