	subroutine arcscatter(x0,y0,x1,y1,arc)
	include 't0.inc'
	integer ee
	double precision x0,y0,x1,y1
	double precision arcxy,arc0,arc,arcstep

	call background_grid
c	arcstep = asin(1.0)*10
	arc0=0
	np=1
	nbs=1
	px(np)=x0
	py(np)=y0
	x=px(np)
	y=py(np)
	call search_boun(x,y,ee,arcxy)
	pstep(np)=0
	arc0=arc0+arcxy
	do 
	  if((arc-arc0)<=0)then
	    arc0=arc
	    x=x1
	    y=y1
		call search_boun(x,y,ee,arcxy)
	    if((arcxy*1.5)>(arc0-pstep(np)))then
	      exit
	    endif
	  elseif((arc-arc0)<arcxy)then
	    arc0=(arc+arc0-arcxy)/2
	    call transfer_xy(arc0)
	    call search_boun(x,y,ee,arcxy)
	  else
	    call transfer_xy(arc0)
	    call search_boun(x,y,ee,arcxy)
	  endif
C	  arcxy=arcstep
	  if((arcxy*1.5)<=(arc0-pstep(np)))then
	    arc0=arc0-arcxy*0.5
	  else
	    np=np+1
	    px(np)=x
	    py(np)=y
	    pstep(np)=arc0
	    arc0=arc0+arcxy
	  endif
	enddo
	np=np+1
	px(np)=x1
	py(np)=y1
	pstep(np)=arc
	end subroutine arcscatter
!----------------------------------------------------------
	subroutine background_grid
	include 't0.inc'
	integer i,j

	open(1,file='bgnod.dat')
	open(2,file='bgnoe.dat')
	open(3,file='bggridt.dat')
	open(4,file='bgstep.dat')
	read(1,*)bgne
	read(3,*)bgnp
	do i=1,bgnp
	  read(3,*)bgpx(i),bgpy(i)
	  read(4,*)bgstep(i)
!	bgstep(i)=0.085
	enddo
	do i=1,bgne
	  read(1,*)bgnod(1,i),bgnod(2,i),bgnod(3,i)
	  read(2,*)bgnoe(1,i),bgnoe(2,i),bgnoe(3,i)
	enddo
	close(1)
	close(2)
	close(3)
	close(4)
	do i=1,bgne,1
	  do j=1,3,1
	    gx(j)=bgpx(bgnod(j,i))
	    gy(j)=bgpy(bgnod(j,i))
	  enddo
	  call circumcircle
	  bgcx(i)=gcx
	  bgcy(i)=gcy
	  bgcr(i)=gcr
	enddo
	end subroutine background_grid
!----------------------------------------------------------
      subroutine search_boun(xx,yy,result,step)
	include 't0.inc'
	double precision xx,yy,step
	integer result,i,j,k,nbg(3),emin
	double precision arr(3,3)
	double precision dcp,bulk(3),v1,v2,v4,vmin,vol(3)
	
	 vmin=1.0
	 emin=0
	 result=0
	 do i=bgne,1,-1
	   nbg(1)=bgnoe(1,i)
	   nbg(2)=bgnoe(2,i)
	   nbg(3)=bgnoe(3,i)
!	   if(nbg(1)/=0.and.nbg(2)/=0.and.nbg(3)/=0)cycle
	   gcx=bgcx(i)
	   gcy=bgcy(i)
	   gcr=bgcr(i)
	   dcp=sqrt((xx-gcx)**2+(yy-gcy)**2)
!         if((dcp-gcr)<gcr*fm2*1000)then
 	     do k=1,3,1
	       do j=1,3,1
	         arr(1,j)=1
	         arr(2,j)=bgpx(bgnod(j,i))
	         arr(3,j)=bgpy(bgnod(j,i))
	       enddo
	       arr(2,k)=xx
	       arr(3,k)=yy
	       call det(arr,3,bulk(k))
	     enddo
	     v1=bulk(1)+bulk(2)+bulk(3)  
	     v2=abs(bulk(1))+abs(bulk(2))+abs(bulk(3))
	     v4=(v2-v1)/v1      
		 if(v4<1e-1)then
	        if(vmin>v4)then
			  vmin=v4
	          emin=i
	          vol(1)=bulk(1)
	          vol(2)=bulk(2)
	          vol(3)=bulk(3)
	        endif
	     endif
!	   endif
	 enddo
	 if(emin==0)then
	   print*,'error'
	   return
	 endif
	 v1=vol(1)+vol(2)+vol(3)
	 v2=abs(vol(1))+abs(vol(2))+abs(vol(3))
	 v4=v1+v2
	 step=0
	 do k=1,3,1
	   vol(k)=(vol(k)+abs(vol(k)))/v4
	   step=step+bgstep(bgnod(k,emin))*vol(k)
	 enddo
	 result=emin
	end subroutine search_boun
!------------------------------------------------------------------
      subroutine circumcircle
	  include 't0.inc'
	  integer i
	  double precision arr(2,2),ver(2),verx(2)
	  double precision x1,x2,x3,r1,r2,r3
	  
	  do i=1,2,1
	    arr(i,1)=gx(i)-gx(3)
		arr(i,2)=gy(i)-gy(3)
		ver(i)=0.5*(arr(i,1)**2+arr(i,2)**2)
	  enddo
	  call axeqb(arr,ver,2,verx)
	  gcx=verx(1)+gx(3)
	  gcy=verx(2)+gy(3)
	          
	  x1=gx(1)-gcx
	  x2=gy(1)-gcy
	  x3=sqrt(x1*x1+x2*x2)
	  r1=x3
	  x1=gx(2)-gcx
	  x2=gy(2)-gcy
	  x3=sqrt(x1*x1+x2*x2)
	  r2=x3
	  r3=sqrt(verx(1)**2+verx(2)**2)
	  gcr=(r1+r2+r3)/3
	end subroutine circumcircle  	  
!---------------------------------------------      
	subroutine axeqb(arr1,verb1,num,verx1)
   	  implicit none
        integer num
        double precision arr1(num,num)
	  double precision arr(num,num)
	  double precision verb1(num),verx1(num)
	  double precision verb(num),verx(num)
        integer i,j,k
	  double precision real_1,real_2
	 
	  do i=1,num
	    do j=1,num
		 arr(i,j)=arr1(i,j)
	    enddo
	  verb(i)=verb1(i)
	  enddo
	  
	  do i=1,num-1,1
		 real_1=arr(i,i)
		 do j=i+1,num,1
		   	if(abs(arr(j,i))>abs(real_1))then
		       do k=i,num,1
			      real_1=arr(i,k)
			      arr(i,k)=arr(j,k)
			      arr(j,k)=real_1
			   enddo
			   real_1=verb(i)
			   verb(i)=verb(j)
			   verb(j)=real_1
			 endif
			 real_1=arr(i,i)
		 enddo
		 real_1=arr(i,i)
		 arr(i,i)=1.0
		 if(abs(real_1)<1e-8)then
		   print*,'error1:3 point in 1 line'
		   print*,i,real_1
	       read*,k
		   return
		 endif
		 do j=i+1,num,1
		    arr(i,j)=arr(i,j)/real_1
		 enddo
		 verb(i)=verb(i)/real_1
		 do j=i+1,num,1
			real_2=arr(j,i)
			arr(j,i)=0.0
			do k=i+1,num,1
			   arr(j,k)=arr(j,k)-real_2*arr(i,k)
			enddo
			verb(j)=verb(j)-real_2*verb(i)
		 enddo
	  enddo

	  if(abs(arr(num,num))<1e-8)then
	    print*,'error2:3 point in 1 line'
		print*,arr(num,num)
	    read*,k
		return
	  endif
	  verx(num)=verb(num)/arr(num,num)
	  do i=num-1,1,-1
	    verx(i)=verb(i)
		do j=num,i+1,-1
		 verx(i)=verx(i)-arr(i,j)*verx(j)
		enddo
	  enddo

	do i=1,num
	  verx1(i)=verx(i)
	enddo
	end	subroutine axeqb

!------------------------------------------------------------------
	subroutine det(matri,num,det_matri)
   	  implicit none
        integer num
        double precision matri(num,num)
	  double precision matrix(num,num)
        integer i,j,k,sign_1
	  double precision det_matri
	  double precision det_matrix,real_1,real_2
   	 
	  do i=1,num
	   do j=1,num
	    matrix(j,i)=matri(i,j)
	   enddo
	  enddo
        det_matrix=1
	  sign_1=1
	  do i=1,num-1,1
	     real_1=matrix(i,i)
		 do j=i+1,num,1
		   	if(abs(matrix(j,i))>abs(real_1))then
		       do k=i,num,1
			      real_1=matrix(i,k)
			      matrix(i,k)=matrix(j,k)
			      matrix(j,k)=real_1
			   enddo
			   real_1=matrix(i,i)
			   sign_1=sign_1*(-1)
			endif
		 enddo
		 real_1=matrix(i,i)
		 if(abs(real_1)<1e-15)then
		   det_matri=0.0
		   return
		 endif
		 matrix(i,i)=1.0
     	 det_matrix=det_matrix*real_1
		 do j=i+1,num,1
			real_2=matrix(j,i)/real_1
			matrix(j,i)=0.0
			do k=i+1,num,1
			   matrix(j,k)=matrix(j,k)-real_2*matrix(i,k)
			enddo
		 enddo
	  enddo
	  det_matrix=det_matrix*matrix(num,num)*sign_1
      det_matri=det_matrix
	end subroutine det
!----------------------------------------------------------
      subroutine loadboun(filename)
	include 't0.inc'
	character filename*32
	integer i
      open(15,file=filename)
	read(15,*)nboun
      do i=1,nboun
	  read(15,*)bx(i),by(i),bs(i)
	enddo
      close(15)
	end subroutine loadboun
!----------------------------------------------------------
      subroutine transfer_xy(arc)
	include 't0.inc'
	double precision arc,ratio
	integer i,n
!sidetype=1 离散形式的边界
	if(sidetype==1)then
        if(bs(nbs)>=arc)then
	    n=1
	    do i=nbs,1,-1
	      if(bs(i)<arc)then
		    n=i
			exit
	      endif
	    enddo
	  elseif(bs(nbs)<arc)then
	    n=nboun
	    do i=nbs,nboun,1
	      if(bs(i)>=arc)then
		    n=i-1
			exit
	      endif
	    enddo
	  endif
	  nbs=n
	  if(n==1.or.n==nboun)then
	    x=bx(n)
	    y=by(n)
	  else
	    ratio=(arc-bs(n))/(bs(n+1)-bs(n))
	    x=bx(n)+ratio*(bx(n+1)-bx(n))
	    y=by(n)+ratio*(by(n+1)-by(n))
	  endif
!sidetype=2 连续边界，函数形式
      elseif(sidetype==2)then
	  call arc_to_xy(arc,x,y)
	elseif(sidetype==3)then
	  call arcn_to_xy(arc,x,y)
	endif
	end subroutine transfer_xy