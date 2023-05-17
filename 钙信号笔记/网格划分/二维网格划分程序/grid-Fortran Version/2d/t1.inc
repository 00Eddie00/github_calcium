	  implicit none
      integer pmax,emax
      parameter(pmax=50000,emax=100000)

	  real px(pmax),py(pmax),alv(pmax),pstep(pmax)
	  real be(3,emax),ce(3,emax),vol(emax)
	  real u(pmax),v(pmax),p(pmax),r(pmax),dr(pmax)
	  real ex(emax),ey(emax)
	  real norm
	  integer nod(3,emax),np,ne
	  common/grid/px,py,np,ne,nod,pstep
	  common/coeffi/be,ce
	  common/geo/vol,alv,ex,ey 
	  common/mech/p,r,u,v,dr,norm