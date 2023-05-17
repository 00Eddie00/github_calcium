	subroutine advance_front
	include 't0.inc'
	integer i,j,k,m,n,nei,pp(3),ee(2),ep,emx,emn,bgep,p2(3),p
	double precision sx(3),sy(3),ss(3),sr(3),f1,f3,step
	double precision xn,yn,vx,vy,lcs,lsp,xx,yy,x1,x2,y1,y2

	f1=0.7
	f3=0.3
      nsearch=0
	do i=1,ne,1
	  epoch(i)=0
	enddo    
	bgsearch=0
	do i=1,bgne,1
	  bgepoch(i)=0
	enddo
      newnp=np
      do i=1,ne,1
	  if(ecrit(i)/=1)cycle
	  pp(1)=nod(1,i)
	  pp(2)=nod(2,i)
	  pp(3)=nod(3,i)
	  do j=1,3,1
	    nei=noe(j,i) 
	    if(nei>0.and.ecrit(nei)/=-1)cycle
!该面的内法向量及外界圆心、半径	
      	m=j+1
	    n=m+1
	    if(m>3)m=m-3
	    if(n>3)n=n-3
	    ee(1)=bgeop(pp(m))
	    ee(2)=bgeop(pp(n))
	    yn=px(pp(n))-px(pp(m))
	    xn=py(pp(m))-py(pp(n))
		lsp=sqrt(xn*xn+yn*yn)
	    xn=xn/lsp
	    yn=yn/lsp
	    lcs=(cx(i)-px(pp(n)))*xn+(cy(i)-py(pp(n)))*yn
	    sx(j)=cx(i)-lcs*xn
	    sy(j)=cy(i)-lcs*yn
	    sr(j)=sqrt(cr(i)*cr(i)-lcs*lcs)
	    emx=max(ee(1),ee(2))
	    emn=min(ee(1),ee(2))
		if(emx==emn)then
	      ss(j)=(pstep(pp(1))+pstep(pp(2)))/2.0
	    else
	      call search_in_bg(sx(j),sy(j),emx,bgep,ss(j))
	      if(bgep==0)then
!	        ss(j)=(pstep(pp(1))+pstep(pp(2)))/2.0
	        if(nei/=0)then
			  print*,'error:search=0',i,j,nei
	          read*,bgep
	        elseif(nei==0)then
	          call search_boun(sx(j),sy(j),bgep,ss(j))
	          if(bgep==0)then
	            print*,'error:search=0',i,j,nei
	            read*,bgep
	          endif
			endif
		  endif
	    endif
!	    if(ss(j)<0.5*sr(j))then
!		  ss(j)=0.5*sr(j)
!		else
          if(ss(j)>3.0*sr(j))then
		  ss(j)=3.0*sr(j) 
		endif   
!推进的新点坐标及其与该面构成的四面体外心及半径
		xx=sx(j)+ss(j)*xn
		yy=sy(j)+ss(j)*yn
		gcr=0.5*ss(j)+0.5*sr(j)*sr(j)/ss(j)
		gcx=xx-gcr*xn
		gcy=yy-gcr*yn
!以下为该点是否应存在的判定
		lcs=sqrt((gcx-px(pp(j)))**2+(gcy-py(pp(j)))**2)
	    if(lcs<gcr)cycle
		if(bgep>0)emx=bgep
		call search_in_bg(xx,yy,emx,bgep,step)
		if(bgep==0)cycle
		lcs=sqrt((xx-px(pp(j)))**2+(yy-py(pp(j)))**2)
		if(lcs<f1*step)cycle
	    p=newnp+1
	if(p>pmax)then
	print*,'newnp>pmax:',newnp,pmax
	read*,p
	endif
	    px(p)=xx
	    py(p)=yy
	    call search_1st(p,i,ep)
		if(ep==0)cycle
	    p2(1)=nod(1,ep)
	    p2(2)=nod(2,ep)
	    p2(3)=nod(3,ep)
		do k=1,3,1
	      lcs=sqrt((px(p2(k))-xx)**2+(py(p2(k))-yy)**2)
	      lsp=pstep(p2(k))*f1
	      if(lcs<f1*step.and.lcs<lsp)then
	        p=0
	        exit
	      endif
	    enddo
	    if(p==0)cycle
		do k=1,3,1
	      nei=noe(k,ep)
	      if(nei/=0)cycle
	      m=k+1
	      n=m+1
	      if(m>3)m=m-3
	      if(n>3)n=n-3
	      vx=py(p2(m))-py(p2(n))
	      vy=px(p2(n))-px(p2(m))
	      lsp=sqrt(vx*vx+vy*vy)
	      vx=vx/lsp
	      vy=vy/lsp
		  lsp=abs(vx*(xx-px(p2(n)))+vy*(yy-py(p2(n))))
	      if(lsp<f3*step)then
		    p=0
	        exit
	      endif
		enddo
	    if(p==0)cycle
	    call search_all(p,ep)
	    x1=xx+step*f1
	    x2=xx-step*f1
	    y1=yy+step*f1
	    y2=yy-step*f1
		do k=2,ndel,1
	      do m=1,3,1
	        if(px(nod(m,edel(k)))>x1)cycle
	        if(px(nod(m,edel(k)))<x2)cycle
	        if(py(nod(m,edel(k)))>y1)cycle
	        if(py(nod(m,edel(k)))>y2)cycle
			vx=px(nod(m,edel(k)))
	        vy=py(nod(m,edel(k)))
	        lcs=sqrt((vx-xx)**2+(vy-yy)**2)
	        if(lcs<f1*step)then
	          p=0
	          exit
	        endif
	      enddo
	      if(p==0)exit
	    enddo
	    if(p==0)cycle
	    newnp=newnp+1
	    pstep(newnp)=step
	    eop(newnp)=ep
	    bgeop(newnp)=bgep
	  enddo
	enddo

	do i=np+1,newnp
	  if(eop(i)==0)cycle
	  do j=i+1,newnp,1
	    if(eop(j)==0)cycle
	    if(pstep(j)<pstep(i))then
	      xx=px(i);px(i)=px(j);px(j)=xx
	      yy=py(i);py(i)=py(j);py(j)=yy
	      step=pstep(i);pstep(i)=pstep(j);pstep(j)=step
	      ep=eop(i);eop(i)=eop(j);eop(j)=ep
	      bgep=bgeop(i);bgeop(i)=bgeop(j);bgeop(j)=bgep
	    endif
	  enddo
	  xx=px(i)
	  yy=py(i)
	  step=pstep(i)*f1
	  x1=xx-step
	  x2=xx+step
	  y1=yy-step
	  y2=yy+step
	  do j=i+1,newnp,1
	    if(px(j)<x1.or.px(j)>x2)cycle
	    if(py(j)<y1.or.py(J)>y2)cycle
	    lcs=sqrt((px(j)-xx)**2+(py(j)-yy)**2)
	    if(lcs<step)then
	      eop(j)=0
	    endif
	  enddo
	enddo

	p=np 
      do i=np+1,newnp,1
	  if(eop(i)==0)cycle
	  p=p+1
	  if(i==p)cycle
	  bgeop(p)=bgeop(i)
	  pstep(p)=pstep(i)
	  eop(p)=eop(i)
	  px(p)=px(i)
	  py(p)=py(i)
	enddo
	newnp=p
	nsearch=0
	do i=1,ne,1
	  epoch(i)=0
	enddo
	end	subroutine advance_front
!--------------------------------------------------------------
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
!--------------------------------------------------------------
      subroutine search_in_bg(xx,yy,start,result,step)
	include 't0.inc'
	double precision xx,yy,step
	integer start,result,i,j,nei,k,last,next,current,e
	double precision arr(3,3)
	double precision dcp,dep,bulk(3),v1,v2,v3,v4
	
	bgsearch=bgsearch+1
      current=start
	result=0
	next=0
	last=0
	do i=1,bgne,1
	  if(bgepoch(current)==bgsearch)exit
	  gcx=bgcx(current)
	  gcy=bgcy(current)
	  gcr=bgcr(current)
	  dcp=sqrt((xx-gcx)**2+(yy-gcy)**2)
	  if((dcp-gcr)<gcr*fm2)then
	    result=current
		exit
	  endif
	  bgepoch(current)=bgsearch
	  dep=-1
	  next=0
	  do j=1,3,1
	    nei=bgnoe(j,current)
	    if(nei<1.or.nei==last.or.bgepoch(nei)==bgsearch)cycle
		gcx=bgcx(nei)
		gcy=bgcy(nei)
		gcr=bgcr(nei)
     	 	dcp=sqrt((xx-gcx)**2+(yy-gcy)**2) 
	    if((dcp-gcr)<gcr*fm2)then
	      result=current
		  exit
	    endif
	    xce=0
	    yce=0
	    do k=1,3,1
	      xce=xce+bgpx(bgnod(k,nei))/3.0
	      yce=yce+bgpy(bgnod(k,nei))/3.0
	    enddo
	    dcp=sqrt((xx-xce)**2+(yy-yce)**2) 
	    if(dep<0.or.dep>dcp)then
	      dep=dcp
	      next=nei
	    endif
	  enddo
	  if(next==0)exit
	  last=current
	  current=next
	enddo
	 
	current=result
	result=0
      if(current>0)then
	  do i=1,3,1
	    do j=1,3,1
	      arr(1,j)=1
	      arr(2,j)=bgpx(bgnod(j,current))
	      arr(3,j)=bgpy(bgnod(j,current))
	    enddo
	    arr(2,i)=xx
	    arr(3,i)=yy
	    call det(arr,3,bulk(i))
	  enddo
	  v1=bulk(1)+bulk(2)+bulk(3) 
	  v2=abs(bulk(1))+abs(bulk(2))+abs(bulk(3))
	  v3=(v2-v1)/v1
	  if(v3<1e-8)then
	    result=current
	    step=0
	    do i=1,3,1
	      bulk(i)=bulk(i)/v1
	      step=step+bgstep(bgnod(i,result))*bulk(i)
	    enddo
	  else
	    call search_bg_all(xx,yy,current)
	    do i=1,ndel,1
	      e=edel(i)
	      do k=1,3,1
	        do j=1,3,1
	          arr(1,j)=1
	          arr(2,j)=bgpx(bgnod(j,e))
	          arr(3,j)=bgpy(bgnod(j,e))
	        enddo
	        arr(2,k)=xx
	        arr(3,k)=yy
	        call det(arr,3,bulk(k))
	      enddo
	      v1=bulk(1)+bulk(2)+bulk(3)  
	      v2=abs(bulk(1))+abs(bulk(2))+abs(bulk(3))
	      v4=(v2-v1)/v1      
		  if(v4<1e-8)then
	        result=e
	        step=0
	        do k=1,3,1
	          bulk(k)=bulk(k)/v1
	          step=step+bgstep(bgnod(k,result))*bulk(k)
			enddo
	        exit
	      endif
	    enddo
	  endif
	endif

	if(result==0)then
	 do i=bgne,1,-1
	   gcx=bgcx(i)
	   gcy=bgcy(i)
	   gcr=bgcr(i)
	   dcp=sqrt((xx-gcx)**2+(yy-gcy)**2)
         if((dcp-gcr)<gcr*fm2)then
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
		 if(v4<1e-8)then
	        result=i
	        step=0
	        do k=1,3,1
	          bulk(k)=bulk(k)/v1
			  step=step+bgstep(bgnod(k,result))*bulk(k)
			enddo
	        exit
	     endif
	   endif
	  enddo
	endif 
	end subroutine search_in_bg   
!--------------------------------------------------------------
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
	   if(nbg(1)/=0.and.nbg(2)/=0.and.nbg(3)/=0)cycle
	   gcx=bgcx(i)
	   gcy=bgcy(i)
	   gcr=bgcr(i)
	   dcp=sqrt((xx-gcx)**2+(yy-gcy)**2)
!        if((dcp-gcr)<gcr*fm2*1000)then
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
!		 if(v4<1e-4)then
	        if(vmin>v4)then
			  vmin=v4
	          emin=i
	          vol(1)=bulk(1)
	          vol(2)=bulk(2)
	          vol(3)=bulk(3)
	        endif
!	     endif
!	   endif
	 enddo
	 if(emin==0)return
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
!--------------------------------------------------------------
      subroutine search_bg_all(xx,yy,ending)
	include 't0.inc'
	integer ending,ele,nei,num,k
	double precision xx,yy,dcp

	ndel=1
      edel(ndel)=ending
	num=1
	bgepoch(ending)=-bgsearch
      do 
	  ele=edel(num)
	  do k=1,3,1
	    nei=bgnoe(k,ele)
	    if(nei<=0.or.bgepoch(nei)==-bgsearch)cycle
	    gcx=bgcx(nei)
          gcy=bgcy(nei)
	    gcr=bgcr(nei)
	    dcp=sqrt((xx-gcx)**2+(yy-gcy)**2)
	    if((dcp-gcr)<(fm2*gcr))then
	      ndel=ndel+1
	      edel(ndel)=nei
	      bgepoch(nei)=-bgsearch
	    endif
	  enddo
	  num=num+1
	  if(num>ndel)exit
      enddo
	end subroutine search_bg_all
!---------------------------------------------------------------
	subroutine smooth
	include 't0.inc'
	integer m,n,pp(3),i,j,k,p,nei,bgep,ep,p2(3)
      double precision sl(3),slmin,stmin,f4,f1,f2,f3,f5
	double precision xx,yy,step,lcs,x1,x2,y1,y2,vx,vy

	f1=0.5
	f2=1.2
	f5=1.0
	f3=0.3
	f4=0.5
      newnp=np
	bgsearch=0
	do i=1,bgne,1
	  bgepoch(i)=0
	enddo
      do i=1,ne,1
	  pp(1)=nod(1,i)
	  pp(2)=nod(2,i)
	  pp(3)=nod(3,i)
	  do j=1,3,1
	    m=j+1
	    n=j+2
	    if(m>3)m=m-3
	    if(n>3)n=n-3
	    sl(j)=sqrt((px(pp(m))-px(pp(n)))**2
     &		+(py(pp(m))-py(pp(n)))**2)
	  enddo
	  slmin=min(sl(1),sl(2),sl(3))
	  stmin=min(pstep(pp(1)),pstep(pp(2)),pstep(pp(3)))
	  if(cr(i)>slmin*f2.or.cr(i)>stmin*f5)then
	    xx=cx(i)
	    yy=cy(i)
	  	call search_in_bg(xx,yy,bgeop(pp(1)),bgep,step)
		if(bgep==0)cycle
	    p=newnp+1
		px(p)=xx
	    py(p)=yy
	    j=min(i,ne)
	    call search_1st(p,j,ep)
		if(ep==0)cycle
	    p2(1)=nod(1,ep)
	    p2(2)=nod(2,ep)
	    p2(3)=nod(3,ep)
		do k=1,3,1
	      nei=noe(k,ep)
	      if(nei/=0)cycle
	      m=k+1
	      n=m+1
	      if(m>3)m=m-3
	      if(n>3)n=n-3
	      vx=py(p2(m))-py(p2(n))
	      vy=px(p2(n))-px(p2(m))
	      lcs=sqrt(vx*vx+vy*vy)
	      vx=vx/lcs
	      vy=vy/lcs
		  lcs=abs(vx*(xx-px(p2(n)))+vy*(yy-py(p2(n))))
	      if(lcs<f3*step)then
		    p=0
	        exit
	      endif
		enddo
	    if(p==0)cycle
	    call search_all(p,ep)
	    do j=1,ndel,1
	      do k=1,3,1
	        vx=px(nod(k,edel(j)))
	        vy=py(nod(k,edel(j)))
	        lcs=sqrt((vx-xx)**2+(vy-yy)**2)
	        if(lcs<f4*step.and.lcs<slmin*f4)then
	          p=0
	          exit
	        endif
	   	  enddo
	      if(p==0)exit
	    enddo
	    if(p==0)cycle
	    newnp=newnp+1
	    px(newnp)=xx
	    py(newnp)=yy
	    bgeop(newnp)=bgep
	    pstep(newnp)=step
	    eop(newnp)=ep
	  endif  
	enddo      

	do i=np+1,newnp
	  if(eop(i)==0)cycle
	  do j=i+1,newnp,1
	    if(eop(j)==0)cycle
	    if(pstep(j)<pstep(i))then
	      xx=px(i);px(i)=px(j);px(j)=xx
	      yy=py(i);py(i)=py(j);py(j)=yy
	      step=pstep(i);pstep(i)=pstep(j);pstep(j)=step
	      ep=eop(i);eop(i)=eop(j);eop(j)=ep
	      bgep=bgeop(i);bgeop(i)=bgeop(j);bgeop(j)=bgep
	    endif
	  enddo
	  xx=px(i)
	  yy=py(i)
	  step=pstep(i)*f1
	  x1=xx-step
	  x2=xx+step
	  y1=yy-step
	  y2=yy+step
	  do j=i+1,newnp,1
	    if(px(j)<x1.or.px(j)>x2)cycle
	    if(py(j)<y1.or.py(J)>y2)cycle
	    lcs=sqrt((px(j)-xx)**2+(py(j)-yy)**2)
	    if(lcs<step)then
	      eop(j)=0
	    endif
	  enddo
	enddo

	p=np 
      do i=np+1,newnp,1
	  if(eop(i)==0)cycle
	  p=p+1
	  if(i==p)cycle
	  bgeop(p)=bgeop(i)
	  pstep(p)=pstep(i)
	  eop(p)=eop(i)
	  px(p)=px(i)
	  py(p)=py(i)
	enddo
	newnp=p
	
	do i=np+1,newnp,1
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
	print*,'smooth over,  np=',np,',  ne=',ne
	end subroutine smooth
!---------------------------------------------------------------
      subroutine optimize
	include 't0.inc'
	real newx(pmax),newy(pmax)
	integer i,times,newn(pmax),p1,p2,p3
      
	print*,'optimizing start!'
      do times=1,100
        do i=1,np
	    newx(i)=0
	    newy(i)=0
	    newn(i)=0
	  enddo
	  do i=1,ne
	    p1=nod(1,i)
	    p2=nod(2,i)
	    p3=nod(3,i)
	    x=px(p1)+px(p2)+px(p3)
	    y=py(p1)+py(p2)+py(p3)
	    newn(p1)=newn(p1)+3
	    newx(p1)=newx(p1)+x
	    newy(p1)=newy(p1)+y
	    newn(p2)=newn(p2)+3
	    newx(p2)=newx(p2)+x
	    newy(p2)=newy(p2)+y
		newn(p3)=newn(p3)+3
	    newx(p3)=newx(p3)+x
	    newy(p3)=newy(p3)+y
	  enddo
	  do i=nps+1,np
	    px(i)=newx(i)/newn(i)
	    py(i)=newy(i)/newn(i)
	  enddo
	enddo
	    
	end subroutine optimize