		implicit none
        integer pmax,emax
		parameter(pmax=100000,emax=200000)

		integer nod(3,emax),bgeop(pmax),nnp,np,ne
		real px(pmax),py(pmax),npx(pmax),npy(pmax)
		real u(pmax),v(pmax),p(pmax),r(pmax)
		real unew(pmax),vnew(pmax),pnew(pmax),rnew(pmax)
		common/grid/nod,px,py,npx,npy,bgeop,ne,nnp,np
		common/value/u,v,p,r,unew,vnew,pnew,rnew