	program main
	include 't0.inc'
	integer iloop, npoch(2000)
	double precision arc,x0,y0,x1,y1,pai,arcboun
	double precision xbo(2000), ybo(2000), arcbo(2000)
	character filename*32
	common/mainp/arcboun
      
      sidetype = 3
	open(11,file='boungridt.dat')
	open(12,file='bounpoch.dat')
!sidetype=1 离散形式的边界
	if(sidetype==1)then
        filename='boun3.dat'
        call loadboun(filename)
	  x0=bx(1)
	  y0=by(1)
	  x1=bx(nboun)
	  y1=by(nboun)
	  arc=bs(nboun)
        call arcscatter(x0,y0,x1,y1,arc)
        call savedata(1)
!sidetype=2 单段连续边界，函数形式，须改动
	elseif(sidetype==2)then
	  pai=2*asin(1.0)
        x0=1.0
	  y0=0
	  x1=1.0
	  y1=0
	  arc=2*pai
        call arcscatter(x0,y0,x1,y1,arc)
        call savedata(4)
!sidetype=3 分段连续边界，函数形式，须改动
	elseif(sidetype==3)then
	  call xyarcbo(xbo, ybo, arcbo,npoch)
	  arcboun = 0
	  do iloop = 1, 17
	    x0 = xbo(iloop)
	    y0 = ybo(iloop)
	    x1 = xbo(iloop+1)
	    y1 = ybo(iloop+1)
	    arc = arcbo(iloop)
          call arcscatter(x0,y0,x1,y1,arc)
	    np = np - 1
          call savendata(npoch(iloop))
	    arcboun = arcboun + arcbo(iloop)
	  enddo
	endif
	close(11)
	close(12)
	end program
!----------------------------------------------------------
!单段连续边界，弧长arc和直角坐标x，y的关系，需改动
	subroutine arc_to_xy(arc,x,y)
	double precision arc,x,y,pai
	pai=2*asin(1.0)
	x=cos(-arc)
	y=sin(-arc)
	end subroutine arc_to_xy
!----------------------------------------------------------
!分段连续边界，各段起止点、弧长
      subroutine xyarcbo(xbo, ybo, arcbo, npoch)
	double precision xbo(100), ybo(100), arcbo(100)
	integer i, npoch(100)
	double precision pai, rad, arc
	pai=2*asin(1.0)
	rad=300
	arc = 0

	arcbo(1)  = 15
	arcbo(2)  = rad*(pai/4)-30
	arcbo(3)  = 30
	arcbo(4)  = rad*(pai/4)-30
	arcbo(5)  = 30
	arcbo(6)  = rad*(pai/4)-30
	arcbo(7)  = 30
	arcbo(8)  = rad*(pai/4)-30
	arcbo(9)  = 30
	arcbo(10) = rad*(pai/4)-30
	arcbo(11) = 30
	arcbo(12) = rad*(pai/4)-30
	arcbo(13) = 30
	arcbo(14) = rad*(pai/4)-30
	arcbo(15) = 30
	arcbo(16) = rad*(pai/4)-30
	arcbo(17) = 15

	npoch(1)  = 2
	npoch(2)  = 1
	npoch(3)  = 2
	npoch(4)  = 1
	npoch(5)  = 2
	npoch(6)  = 1
	npoch(7)  = 2
	npoch(8)  = 1
	npoch(9)  = 2
	npoch(10) = 1
	npoch(11) = 2
	npoch(12) = 1
	npoch(13) = 2
	npoch(14) = 1
	npoch(15) = 2
	npoch(16) = 1
	npoch(17) = 2

      do i = 1, 18
	  xbo(i)=rad*cos(arc/rad)
	  ybo(i)=rad*sin(arc/rad)
	  arc = arc+arcbo(i)
	enddo

	end subroutine xyarcbo
!----------------------------------------------------------
!分段连续边界，弧长arc和直角坐标x，y的关系，需改动
	subroutine arcn_to_xy(arc,x,y)
	double precision arc,x,y,rad
	double precision arcboun
	common/mainp/arcboun
	rad = 300
	x=rad*cos((arc+arcboun)/rad)
	y=rad*sin((arc+arcboun)/rad)
	end subroutine arcn_to_xy
!----------------------------------------------------------
      subroutine savedata(npoch)
	include 't0.inc'
	integer i, npoch
	do i=1,np-1
	  write(11,*)px(i),py(i)
	  write(12,*)npoch
	enddo
	end subroutine savedata
!----------------------------------------------------------
      subroutine savendata(npoch)
	include 't0.inc'
	integer i, npoch
	do i=1,np
	  write(11,*)px(i),py(i)
	enddo
	if (npoch .eq. 1) then
	  np = np - 1
	else
	  np = np + 1
	endif

	do i = 1,np
	  write(12,*)npoch
	enddo

	end subroutine savendata