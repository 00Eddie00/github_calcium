program aa
integer npoch(50000),sideface(3)

OPEN(1,FILE='SIDEGRIDT.dat')
READ(1,*)NBP
CLOSE(1)
OPEN(1,FILE='GRIDT.dat')
READ(1,*)NP
CLOSE(1)
do i = 1,nbp,1
   npoch(i)=1;
enddo
do i = nbp+1,np,1
npoch(i)=0;
enddo
open(1,file='SIDEFACE.dat')
read(1,*)nS
do i=1,nS
  read(1,*)(SIDEFACE(j),j=1,3)
  if (sideface(3) .ne. 1)then
    npoch(sideface(1))=sideface(3);
    npoch(sideface(2))=sideface(3);
  endif
enddo
close(1)

open(1,file='npoch.dat')
do i=1,np
write(1,*)npoch(i)
enddo
close(1)
end program

