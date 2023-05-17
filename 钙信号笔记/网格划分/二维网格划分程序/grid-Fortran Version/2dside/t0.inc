      implicit none
   
      integer pmax,emax,bmax	
      real fmin,fm2
      parameter(pmax=100000,emax=200000,bmax=50000)
      parameter(fmin=1e-8,fm2=1e-8)                                                  
	  												
      integer np,nboun,sidetype,nbs, bounarc						     
      real px(bmax),py(bmax),pstep(bmax)
      double precision bx(bmax),by(bmax),bs(bmax)			 
      common/grid1/np,nboun,sidetype,nbs,bounarc
      common/grid2/px,py,pstep,bx,by,bs

	  real gx(3),gy(3)
	  double precision gcx,gcy,gcr,x,y,xce,yce
	  common/dat1/gx,gy,gcx,gcy,gcr,x,y,xce,yce

	  integer bgnp,bgne
      real bgpx(pmax),bgpy(pmax),bgstep(pmax)
	  integer bgnod(3,emax),bgnoe(3,emax)
      double precision bgcx(emax),bgcy(emax),bgcr(emax)			 
      common/bg1/bgnp,bgne,bgstep
      common/bg2/bgpx,bgpy,bgcx,bgcy,bgcr,bgnod,bgnoe
