       pi=4*atan(1.)
       open(1,file='tfactor.in',status='old')
       open(2,file='tfactor.out',status='unknown')
       do it=1,15
       read(1,'(t23,f6.4)') eps
       read(1,*)
       read(1,*)
       read(1,*)
       tf=0.
       do i=1,500
       read(1,*) dum,th
       tf=tf+(cos(pi*th/180))**2/500
       enddo
       write(2,*) eps,tf
       enddo
       end
