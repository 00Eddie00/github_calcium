      implicit none
   
      integer pmax,emax,delmax,bmax	
      real fmin,fm2
      parameter(pmax=100000,emax=200000,delmax=1000,bmax=2000)
      parameter(fmin=1e-8,fm2=1e-8)                                                  
	  												
      integer np,ne,nsearch						     
      real px(pmax),py(pmax),ipx(4),ipy(4)
	  real centx,centy,length         
      integer nod(3,emax),noe(3,emax),epoch(emax)           
      double precision cx(emax),cy(emax),cr(emax)			 
      common/grid1/np,ne,nsearch
      common/grid2/px,py,cx,cy,cr,nod,noe,epoch 
	  common/convex/ipx,ipy,centx,centy,length

	  real gx(3),gy(3)
	  double precision gcx,gcy,gcr,x,y,xce,yce
	  common/dat1/gx,gy,gcx,gcy,gcr,x,y,xce,yce

	  integer ndel,edel(delmax),neibn(delmax)
	  common/deldat/ndel,edel,neibn

	  integer eos(2,bmax),side(3,bmax),nb,nps
	  common/bound/eos,side,nb,nps

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer newnp,ecrit(emax),eop(pmax),bgeop(pmax)
	  real pstep(pmax)
	  common/adf/newnp,ecrit,eop,bgeop,pstep

	  integer bgnp,bgne,bgfail,bgsearch						     
      real bgpx(pmax),bgpy(pmax),bgstep(pmax)
	  integer bgnod(3,emax),bgnoe(3,emax),bgepoch(emax)          
      double precision bgcx(emax),bgcy(emax),bgcr(emax)			 
      common/bg1/bgnp,bgne,bgstep,bgfail,bgsearch,bgepoch
      common/bg2/bgpx,bgpy,bgcx,bgcy,bgcr,bgnod,bgnoe
