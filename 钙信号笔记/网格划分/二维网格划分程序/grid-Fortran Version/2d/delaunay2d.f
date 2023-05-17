	program main
	include 't0.inc'
	integer i,j,k,nei,nfrontele,ep
	double precision xx,yy,step

	call boun_grid
	print*,'boun_grid ok'
	do i=1,ne,1
	  ecrit(i)=0
	  do j=1,3,1
	    nei=noe(j,i)
	    if(nei==0)then
	      ecrit(i)=1
	      exit
	    endif
	  enddo
	enddo

	call background_grid
	print*,'background grid ok'

!	open(1,file='bgeop.dat')
!	open(2,file='pstep.dat')
	bgsearch=0
	do i=1,bgne
	  bgepoch(i)=0
	enddo
      do i=1,np,1
	  xx=px(i)
	  yy=py(i)
	  call search_in_bg(xx,yy,1,ep,step)
	  if(ep==0)then
	    call search_boun(xx,yy,ep,step)
	    if(ep==0)then   
		  print*,'error: search boungrid in bg no found',i
	      read*,ep
	    endif
	  endif
	  bgeop(i)=ep
	  pstep(i)=step
!	  read(1,*)bgeop(i)
!	  read(2,*)pstep(i)
!	pstep(i)=0.08
	enddo
!	close(1)
!	close(2)
      do 
	  call advance_front
	  do i=1,ne,1
	    if(ecrit(i)/=1)cycle
		ecrit(i)=-1
	  enddo
	  print*,'np=',newnp
	  do i=np+1,newnp,1
	if(i==617)then
	print*,i
	endif
	    k=eop(i)
	    if(k>ne)k=ne
	    call search_1st(i,k,j)
	if(j==0)then
	print*,'error:search_ist=0/i=',i
	read*,j
	exit
	endif
	    call search_all(i,j)
	if(ndel>200)then
	print*,'ndel=',ndel,'/i=',i
	read*,k
	endif
	    call check_up(i)
	    call clear_up(i)
	    if(mod(i,1000)==0)print*,'i=',i,ne
	  enddo
	  np=newnp
	  nfrontele=0
	  do i=1,ne,1
	    if(ecrit(i)==-1)cycle
	    nfrontele=nfrontele+1
	    ecrit(i)=0
	    do j=1,3,1
	      nei=noe(j,i)
	      if(nei==0)then
	        ecrit(i)=1
	        exit
	      elseif(ecrit(nei)==-1)then
	        ecrit(i)=1
	        exit
	      endif
	    enddo
	  enddo

	  print*,'frontele=',nfrontele	 
	  if(nfrontele==0)exit  	      
	enddo

	do 
	  i=np
	  call smooth
	  if(i==np)exit
	enddo
	call optimize

	open(1,file='gridt.dat')
	write(1,*)np
	do i=1,np,1
	  write(1,*)px(i),py(i)
	enddo
	close(1)
	open(1,file='bgeop.dat')
	do i=1,np,1
	  write(1,'(i8)')bgeop(i)
	enddo
	close(1)
	open(1,file='pstep.dat')
	do i=1,np,1
	  write(1,*)pstep(i)
	enddo
	close(1)
      call savedata
	end program
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subroutine boun_grid
	include 't0.inc'
	integer i,j,lastnp,k
	real xmax,xmin,ymax,ymin,len_x,len_y
	
      open(1,file='sidegridt.dat')
	read(1,*)np
	nps = np
	do i=1,np,1
	  read(1,*)px(i),py(i)

	  if(i==1)then
	    xmax=px(i)
	    xmin=px(i)
	    ymax=py(i)
	    ymin=py(i)
	  else
	    if(xmax<px(i))xmax=px(i)
	    if(xmin>px(i))xmin=px(i)
	    if(ymax<py(i))ymax=py(i)
	    if(ymin>py(i))ymin=py(i)
	  endif
	enddo
	close(1)
	print*,'read ok'
	centx=(xmax+xmin)/2
	centy=(ymax+ymin)/2
	len_x=xmax-xmin
	len_y=ymax-ymin
	length=max(len_x,len_y)
	length=1.2*length
	call convex_hull
	ne=2
	lastnp=0
	nsearch=0
!	call loaddata(lastnp)
	do i=lastnp+1,np,1
	  call search_1st(i,ne,j)
	  call search_all(i,j)
	  if(ndel>200)then
	    print*,'ndel=',ndel,'/i=',i
	    read*,k
	  endif
	  call check_up(i)
	  call clear_up(i)
	if(mod(i,100)==0)print*,'i=',i,ne
!	if(i==216)then
!	  print*,'checkside start'
!	  call checkside
!	  print*,'checkside finished'
!	  call delete_outside
!	endif
	enddo
	print*,'checkside start'
	  call checkside
	  print*,'checkside finished'
	  call delete_outside
	call savedata
	
      end subroutine boun_grid

!---------------------------------------------------------------
	subroutine convex_hull
	include 't0.inc'
	integer i,j

	ipx(1)=centx-0.65*length
	ipy(1)=centy-0.6*length
	ipx(2)=centx+0.55*length
	ipy(2)=centy-0.6*length
	ipx(3)=centx+0.6*length
	ipy(3)=centy+0.6*length
      ipx(4)=centx-0.6*length
	ipy(4)=centy+0.6*length
      
	nod(1,1)=-1; nod(2,1)=-2; nod(3,1)=-4
      nod(1,2)=-2; nod(2,2)=-3; nod(3,2)=-4

	noe(1,1)=2; noe(2,1)=0; noe(3,1)=0
      noe(1,2)=0; noe(2,2)=1; noe(3,2)=0
	
	do i=1,2,1
	  do j=1,3,1
	    gx(j)=ipx(-nod(j,i))
	    gy(j)=ipy(-nod(j,i))
	  enddo
	  call circumcircle
	  cx(i)=gcx
	  cy(i)=gcy
	  cr(i)=gcr
	enddo
	end subroutine convex_hull
!-----------------------------------------------------------------
	subroutine search_1st(p,start,ending)
	include 't0.inc'
	integer p,start,ending
	integer current,last,next,key,i,j,k
	double precision dcp,dep

	x=px(p)
	y=py(p)
	nsearch=nsearch+1
      last=0
	next=0
	ending=0
	current=start
	do i=1,ne,1
	  if(epoch(current)==nsearch)exit
	  gcx=cx(current)
	  gcy=cy(current)
	  gcr=cr(current)
	  dcp=sqrt((x-gcx)**2+(y-gcy)**2)
	  if(dcp<gcr)then
	    ending=current
		exit
	  endif
	  epoch(current)=nsearch
	  dep=-1
	  next=0
	  do j=1,3,1
	    k=noe(j,current)
	    if(k<1.or.k==last.or.epoch(k)==nsearch)cycle
		gcx=cx(k)
		gcy=cy(k)
		gcr=cr(k)
     	 	dcp=sqrt((x-gcx)**2+(y-gcy)**2) 
	    if(dcp<gcr)then
	      ending=current
		  exit
	    endif
	    call centroid(k)
	    dcp=sqrt((x-xce)**2+(y-yce)**2) 
	    if(dep<0.or.dep>dcp)then
	      dep=dcp
	      next=k
	    endif
	  enddo
	  if(next==0)exit
	  last=current
	  current=next
	enddo
	 
	current=ending
	ending=0
      if(current>0)then
	  call in_or_out(p,current,key)
	  if(key==-1)then
	    call search_all(p,current)
	    do i=1,ndel,1
	      j=edel(i)
	      call in_or_out(p,j,k)
	      if(k>-1)then
	        ending=j
	        exit
	      endif
	    enddo
	  else
	    ending=current
	  endif
	endif

	if(ending==0)then
	 do i=ne,1,-1
	   gcx=cx(i)
	   gcy=cy(i)
	   gcr=cr(i)
	   dcp=sqrt((x-gcx)**2+(y-gcy)**2)
         if(dcp<gcr)then
 	    call in_or_out(p,i,key)
	    if(key==-1)then
		  cycle
	    else
		  ending=i
	      exit
	    endif
	   endif
	  enddo
	endif 
	  
	end subroutine search_1st   
!-----------------------------------------------------------------
      subroutine search_all(p,ending)
	include 't0.inc'
	integer p,ending
	integer ele,nei,num,k,j,add
	double precision dcp

	ndel=1
      edel(ndel)=ending
	x=px(p)
	y=py(p)
	num=1
      do 
	  ele=edel(num)
	  do k=1,3,1
	    nei=noe(k,ele)
	    if(nei<=0)cycle
	    gcx=cx(nei)
          gcy=cy(nei)
	    gcr=cr(nei)
	    dcp=sqrt((x-gcx)**2+(y-gcy)**2)
	    if((dcp-gcr)<(fm2*gcr))then
	      add=1
	      do j=1,ndel,1
	        if(edel(j)==nei)then
	          add=0
	          exit
	        endif
	      enddo
	      if(add==1)then
	        ndel=ndel+1
	        edel(ndel)=nei
	      endif
	    endif
	  enddo
	  num=num+1
	  if(num>ndel)exit
      enddo
	
      end subroutine search_all
!----------------------------------------------------------------
      subroutine check_up(p)
	include 't0.inc'
	integer p,j,k,ele,nei,ndel1,pp(3),flag
	double precision vol1,vol2

      do j=1,ndel,1
	  epoch(edel(j))=-nsearch
	enddo
	ndel1=ndel
	do 
	  flag=0
	  do j=1,ndel1,1
	    ele=edel(j)
	    do k=1,3,1
	      pp(k)=nod(k,ele)
	    enddo
	    do k=1,3,1
	      nei=noe(k,ele)
	      if(nei>0.and.epoch(nei)==-nsearch)then
		    cycle
		  else 
	        pp(k)=p
	        call area(pp(1),pp(2),pp(3),vol1)
	        pp(k)=nod(k,ele)
	        call distance(pp(1),pp(2),pp(3),vol2,k)
	        if(vol1/vol2<1e-8)then
	          flag=1
	          exit
			endif
	      endif
	    enddo
	   	if(flag==1)then
	        epoch(ele)=0
	        edel(j)=edel(ndel1)
	        ndel1=ndel1-1
	        exit
	    endif
	  enddo
	  if(flag==0)exit
	enddo
	ndel=ndel1
	
	end subroutine check_up
!----------------------------------------------------------------
      subroutine clear_up(p)
	include 't0.inc'
	integer p,i,j,k,ele,fill,nei,nnew,arr1(3),arr2(3),n1,n2,n3
	
      nnew=ne
	do i=1,ndel,1
	  ele=edel(i)
	  do j=1,3,1
	      nei=noe(j,ele)
	      if(nei<=0)then
	        nnew=nnew+1
	        nod(j,nnew)=p
	        do k=1,3,1
	         if(k==j)cycle
			 nod(k,nnew)=nod(k,ele)
			enddo
			noe(j,nnew)=nei
			neibn(nnew-ne)=neibn(nnew-ne)+1  
		  elseif(epoch(nei)/=-nsearch)then	
		    nnew=nnew+1
	        nod(j,nnew)=p
	        do k=1,3,1
	         if(k==j)cycle
			 nod(k,nnew)=nod(k,ele)
			enddo
			noe(j,nnew)=nei 
	        neibn(nnew-ne)=neibn(nnew-ne)+1
	        do k=1,3,1
	         if(noe(k,nei)==ele)then
	            noe(k,nei)=nnew
	            exit
	         endif
	        enddo
	      endif
	  enddo
	enddo
	
      do i=ne+1,nnew-1,1
	  if(neibn(i-ne)==3)cycle
	  do k=1,3,1
	    arr1(k)=nod(k,i)
	  enddo
	  do j=i+1,nnew,1
	    if(neibn(j-ne)==3)cycle
	    do k=1,3,1
		  arr2(k)=nod(k,j)
		enddo   
	    call common_point(arr1,arr2,n1,n2,n3)
	    if(n3/=2)cycle
	    neibn(i-ne)=neibn(i-ne)+1
	    neibn(j-ne)=neibn(j-ne)+1
	    noe(n1,i)=j
	    noe(n2,j)=i
	    if(neibn(i-ne)==3)exit
	  enddo
	enddo
	do i=1,nnew-ne,1
	  neibn(i)=0
	epoch(i+ne)=0
	ecrit(i+ne)=0
	enddo

	fill=nnew
      do i=1,ndel,1
	  ele=edel(i)
	  if(ele>nnew-ndel)cycle
	  do j=1,3,1
	    nod(j,ele)=nod(j,fill)
	    noe(j,ele)=noe(j,fill)
	  enddo
	  ecrit(ele)=ecrit(fill)
	  epoch(ele)=epoch(fill)
	  do j=1,3,1
	    nei=noe(j,ele)
	    if(nei<=0)cycle
	    do k=1,3,1
	      if(noe(k,nei)==fill)then
	        noe(k,nei)=ele
	      endif
	    enddo
	  enddo
	  fill=fill-1
	  do
	    if(epoch(fill)==-nsearch)then
		  fill=fill-1
		else
		  exit
		endif
	  enddo	      
	  do j=1,3,1
	    k=nod(j,ele)
	    if(k<0)then
		  gx(j)=ipx(-k)
	      gy(j)=ipy(-k)
	    else
	      gx(j)=px(k)
	      gy(j)=py(k)
	    endif
	  enddo
	  call circumcircle
	  cx(ele)=gcx
	  cy(ele)=gcy
	  cr(ele)=gcr
	enddo
	
      do i=ne+1,nnew-ndel,1
	  ele=i
	  do j=1,3,1
	    k=nod(j,ele)
	    if(k<0)then
		  gx(j)=ipx(-k)
	      gy(j)=ipy(-k)
	    else
	      gx(j)=px(k)
	      gy(j)=py(k)
	    endif
	  enddo
	  call circumcircle
	  cx(ele)=gcx
	  cy(ele)=gcy
	  cr(ele)=gcr
	enddo
      ne=nnew-ndel 
      end subroutine clear_up
!------------------------------------------------------------------
      subroutine common_point(arr1,arr2,n1,n2,n3)      
      implicit none
      integer arr1(3),arr2(3)
	integer n1,n2,n3,i,j
      					   
      n3=0
	n1=6
	n2=6
	do i=1,3,1
	  do j=1,3,1
	    if(arr2(j)==arr1(i))then
	      n3=n3+1
	      n1=n1-i
	      n2=n2-j
	      exit
	    endif
	  enddo
      enddo
      end subroutine common_point
!----------------------------------------------------------------      
	subroutine in_or_out(p,e,key)
	include 't0.inc'
      integer p,e,key
      double precision v1,v2,v3,v,vv
      
	call area(p,nod(2,e),nod(3,e),v1)
	call area(nod(1,e),p,nod(3,e),v2)
	call area(nod(1,e),nod(2,e),p,v3)

      v=v1+v2+v3
	vv=(abs(v1)+abs(v2)+abs(v3)-v)/v
      if(vv<1e-8)then
        key=1
	else
        key=-1
      end if
      end subroutine in_or_out
!------------------------------------------------------------------
	subroutine centroid(e)
	include 't0.inc'
	integer i,e,pp(3)

	xce=0
	yce=0
	do i=1,3,1
	  pp(i)=nod(i,e)
	  if(pp(i)>0)then
	    xce=xce+px(pp(i))/3.0
	    yce=yce+py(pp(i))/3.0
	  else
	    xce=xce+ipx(-pp(i))/3.0
	    yce=yce+ipy(-pp(i))/3.0
	  endif
	enddo
	end subroutine centroid
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
	subroutine area(p1,p2,p3,vol)
	include 't0.inc'
      integer p1,p2,p3,i,pp(3)
	double precision vol,arr(3,3),det_arr
      
	pp(1)=p1
	pp(2)=p2
	pp(3)=p3
	do i=1,3,1
	  arr(1,i)=1
	  if(pp(i)>0)then
	    arr(2,i)=px(pp(i))
	    arr(3,i)=py(pp(i))
	  else
	    arr(2,i)=ipx(-pp(i))
	    arr(3,i)=ipy(-pp(i))
	  endif
	enddo
	call det(arr,3,det_arr)
	vol=det_arr/2.0
	end subroutine area   
!------------------------------------------------------------------
      subroutine distance(p1,p2,p3,vol,n)
	include 't0.inc'
      integer p1,p2,p3,n,pp(3)
	double precision x1,y1,x2,y2,vol
      
	pp(1)=p1
	pp(2)=p2
	pp(3)=p3
	pp(n)=p3
      if(pp(1)>0)then
	  x1=px(pp(1))
	  y1=py(pp(1))
	else
	  x1=ipx(-pp(1))
	  y1=ipy(-pp(1))
	endif
	if(pp(2)>0)then
	  x2=px(pp(2))
	  y2=py(pp(2))
	else
	  x2=ipx(-pp(2))
	  y2=ipy(-pp(2))
	endif
	vol=sqrt((x1-x2)**2+(y1-y2)**2)
	end subroutine distance
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
!---------------------------------------------------------------
      subroutine loaddata(lastnp)
	  include 't0.inc' 
	  integer lastnp,i
	  
	open(1,file='nod.dat')
	open(2,file='noe.dat')
	open(3,file='cir.dat')
	read(1,*)ne
	do i=1,ne,1
	  read(1,*)nod(1,i),nod(2,i),nod(3,i)
	  read(2,*)noe(1,i),noe(2,i),noe(3,i)
	  read(3,*)cx(i),cy(i),cr(i)
	enddo
	read(1,*)lastnp
	close(1)
	close(2)
	close(3)
	print*,'loaddata ok'
       
	end subroutine loaddata       	
!----------------------------------------------------------------
      subroutine savedata
	include 't0.inc'
	integer i

	open(1,file='nod.dat')
	open(2,file='noe.dat')
	open(3,file='cir.dat')
	write(1,*)ne
	do i=1,ne,1
	write(1,*)nod(1,i),nod(2,i),nod(3,i)
	write(2,*)noe(1,i),noe(2,i),noe(3,i)
	write(3,*)cx(i),cy(i),cr(i)
	enddo
	write(1,*)np
	close(1)
	close(2)
	close(3)
	end subroutine savedata
!----------------------------------------------------------------
      subroutine checkside
	include 't0.inc'
	integer pp1(3),pp2(3)
	integer i,j,k,p,start,ending,nei,n1,n2,n3
	double precision vol
	
	open(1,file='sideface.dat')
	read(1,*)nb
	do i=1,nb,1
	  read(1,*)side(1,i),side(2,i),side(3,i)
	enddo
	close(1)

	nsearch=0
	do j=1,ne,1
	  epoch(j)=0
	enddo
      start=1
	do i=1,nb,1
	  k=mod(i,100)+1
        p=np+k
	  pp1(1)=side(1,i)
	  pp1(2)=side(2,i)
	  pp1(3)=0
	  px(p)=(px(side(1,i))+px(side(2,i)))/2.0
	  py(p)=(py(side(1,i))+py(side(2,i)))/2.0
	  call search_1st(p,start,ending)
	  pp2(1)=nod(1,ending)
	  pp2(2)=nod(2,ending)
	  pp2(3)=nod(3,ending)
	  call common_point(pp1,pp2,n1,n2,n3)
	  if(n3/=2)then
	    do  j=1,3
	      nei=noe(j,ending)
	      pp2(1)=nod(1,nei)
	      pp2(2)=nod(2,nei)
	      pp2(3)=nod(3,nei)
	      call common_point(pp1,pp2,n1,n2,n3)
	      if(n3==2.and.n1==3)then
		    ending=nei
			exit
	      endif
	    enddo
	  endif
	  if(n3/=2)then
	    call search_all(p,ending)
	    do j=1,ndel
	      nei=edel(j)
	      pp2(1)=nod(1,nei)
	      pp2(2)=nod(2,nei)
	      pp2(3)=nod(3,nei)
	      call common_point(pp1,pp2,n1,n2,n3)
	      if(n3==2.and.n1==3)then
		    ending=nei
			exit
	      endif
	    enddo
	  endif
	  if(n3==2.and.n1==3)then
	    start=ending
	    call area(pp1(1),pp1(2),pp2(n2),vol)
	    if(vol<0)then
	      eos(2,i)=ending
	      eos(1,i)=noe(n2,ending)
	      pp2(1)=nod(1,noe(n2,ending))
	      pp2(2)=nod(2,noe(n2,ending))
	      pp2(3)=nod(3,noe(n2,ending))
	      j=noe(n2,ending)
	      call common_point(pp1,pp2,n1,n2,n3)
	      call area(pp1(1),pp1(2),pp2(n2),vol)
	      if(vol<0.or.n1/=3.or.n3/=2.or.noe(n2,j)/=ending)
     &		  print*,'error'
	    else
	      eos(1,i)=ending
	      eos(2,i)=noe(n2,ending)
	      pp2(1)=nod(1,noe(n2,ending))
	      pp2(2)=nod(2,noe(n2,ending))
	      pp2(3)=nod(3,noe(n2,ending))
	      j=noe(n2,ending)
	      call common_point(pp1,pp2,n1,n2,n3)
	      call area(pp1(1),pp1(2),pp2(n2),vol)
	      if(vol>0.or.n1/=3.or.n3/=2.or.noe(n2,j)/=ending)
     &		  print*,'error'  		  
		endif
	  else
	    print*,'error:ending do not include s',i,ending
	    read*,j
	  endif
	  if(mod(i,100)==0)then
	    print*,'i=',i
	    nsearch=0
	    do j=1,ne,1
	      epoch(j)=0
	    enddo
	  endif
	enddo
	end	subroutine checkside
!--------------------------------------------------------------
	subroutine delete_outside
	include 't0.inc'
	integer i,j,k,e1,e2,nei,flag

	do i=1,nb,1
	  e1=eos(1,i)
	  e2=eos(2,i)
	  flag=0
	  do j=1,3,1
	    nei=noe(j,e1)
	    if(nei/=e2)cycle
	    noe(j,e1)=0
	    flag=flag+1
	  enddo
	  do j=1,3,1
	    nei=noe(j,e2)
	    if(nei/=e1)cycle
	    noe(j,e2)=0
	    flag=flag+1
	  enddo
	  if(flag/=2)print*,'error'
	enddo
	print*,'inside and outside separate finished!'
	
      do i=1,nb,1
	  e2=eos(2,i)
	  nod(1,e2)=0
	  nod(2,e2)=0
	  nod(3,e2)=0
	  do j=1,3,1
	    nei=noe(j,e2)
	    if(nei==0)cycle
		nod(1,nei)=0
	    noe(j,e2)=0
	  enddo
	enddo

	e1=0
      do 
	  flag=0
	  do i=1,ne,1
	    if(nod(1,i)/=0)cycle
	    if(nod(2,i)==0)cycle
	    flag=flag+1
	    nod(2,i)=0
	    nod(3,i)=0
	    do j=1,3,1
	      nei=noe(j,i)
	      if(nei==0)cycle
	      nod(1,nei)=0
	      noe(j,i)=0
	    enddo
	  enddo
	  print*,e1
	  e1=e1+1
	  if(flag==0)exit
	enddo
	print*,'delete ok, times=',e1
	
      e1=0	
	do i=1,ne,1
	  if(nod(1,i)==0)cycle
	  e1=e1+1
	  if(i==e1)cycle
	  do j=1,3,1
	    nod(j,e1)=nod(j,i)
	    noe(j,e1)=noe(j,i)
	    cx(e1)=cx(i)
		cy(e1)=cy(i)
	    cr(e1)=cr(i)
	  enddo
	  do j=1,3,1
	    nei=noe(j,i)
	    if(nei==0)cycle
	    do k=1,3,1
		  if(noe(k,nei)==i)noe(k,nei)=e1
		enddo
	  enddo
	enddo
	ne=e1   	       
	print*,'rewrite delete ok,   ne=',ne
	end subroutine delete_outside			  