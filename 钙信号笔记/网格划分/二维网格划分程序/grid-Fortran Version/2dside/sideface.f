	program main

	integer boun(100), nboun, npoch(10000)

      nboun = 2
	boun(1)=70
	boun(2)=846

	open(1,file='sideface.dat')
	open(2,file='sidegridt.dat')
	open(3,file='sidenpoch.dat')
	read(2,*)nb
      write(1,*)nb
	do i = 1,nb
	  read(3,*)npoch(i)
	enddo

      istart = 1
	do i=1,nboun
	  iend = boun(i)
	  do j = istart,iend-1
	    k = npoch(j)
	    if (k .ne. npoch(j+1)) then
	      k = 1
	    endif
	    write(1,*)j,j+1,k
	  enddo
	  k = npoch(iend)
	  if (k .ne. npoch(istart)) then
	    k = 1
	  endif
	  write(1,*)iend,istart,k
	  istart = iend+1
	enddo
	close(1)
	close(2)
	close(3)
      end