      subroutine knnimp(x,ximp,p,n,imiss,irmiss,kn,workp,workn,iworkp,
&     iworkn)
      integer p,n,i,j,k,kn,m
      integer imiss(p,n),irmiss(p)
      double precision x(p,n),ximp(p,n)
      integer iworkp(p),iworkn(n), jj
      double precision workp(p),workn(n)
      m=kn+1
      do 23000 i=1,p
      if(.not.(irmiss(i)  .eq.  0))goto 23002
      goto 23000
23002 continue
      do 23004 j=1,n
      workn(j)=x(i,j)
      iworkn(j)=imiss(i,j)
23004 continue
      call misdis(workn,x,p,n,iworkn,imiss,workp,iworkp)
      call porder(m,workp,p,iworkp,workn)
      call misave(x,workn,p,n,iworkn,imiss,iworkp(2),kn)
      do 23006 k=1,n
      if(.not.(iworkn(k) .eq. 0))goto 23008
      goto 23006
23008 continue
      ximp(i,k)=workn(k)
      if(.not.(iworkn(k) .eq. 2))goto 23010
      imiss(i,k)=2
23010 continue
23006 continue
23000 continue
      return
      end
      subroutine misdis(x0,x,p,n,imiss0,imiss,dis,iworkp)
      integer p,n,imiss(p,n),iworkp(p),imiss0(n)
      double precision x0(n),x(p,n),dis(p),dismax
      dismax=1d10
      do 23012 j=1,p 
      dis(j)=0d0
      iworkp(j)=0
23012 continue
      do 23014 k=1,n
      if(.not.(imiss0(k) .gt. 0))goto 23016
      goto 23014
23016 continue
      do 23018 j=1,p 
      if(.not.(imiss(j,k) .gt. 0))goto 23020
      goto 23018
23020 continue
      dis(j)= dis(j)+(x0(k)-x(j,k))**2
      iworkp(j)=iworkp(j)+1
23018 continue
23014 continue
      do 23022 j=1,p
      if(.not.(iworkp(j) .gt. 0))goto 23024
      dis(j)=dis(j)/iworkp(j)
      goto 23025
23024 continue
      dis(j)=dismax
23025 continue
23022 continue
      return
      end
      subroutine misave(x,x0,p,n,imiss0,imiss,index,m)
      integer p,n,imiss0(n),imiss(p,n),index(m),ktot
      double precision x(p,n),x0(n)

      do 23026 k=1,n
      x0(k)=0d0
      if(.not.(imiss0(k) .eq. 0))goto 23028
      goto 23026
23028 continue
      ktot=0
      do 23030 j=1,m
      jj=index(j)
      if(.not.(imiss(jj,k) .gt. 0))goto 23032
      goto 23030
23032 continue
      x0(k)=x0(k)+x(jj,k)
      ktot=ktot+1
23030 continue
      if(.not.(ktot .gt. 0))goto 23034
      x0(k)=x0(k)/ktot
      goto 23035
23034 continue
      imiss0(k)=2
23035 continue
23026 continue
      return
      end
	subroutine porder(kn,dist, ntr,pos,nndist)
        integer kn, ntr
	double precision  dist(ntr), nndist(kn)
	integer pos(kn)
	   do 50 j = 1, ntr
	      if(j .le. kn) then
	         do 30 k = 1, j-1
	            if(dist(j) .lt. nndist(k)) then
	               do 20 k1 = j-1, k, -1
	                  nndist(k1+1) = nndist(k1)
20	                  pos(k1+1) = pos(k1)
	               nndist(k) = dist(j)
	               pos(k) = j
	               goto 50
	            endif
30	         continue
	         nndist(j) = dist(j)
	         pos(j) = j
	      else
	         if(dist(j) .ge. nndist(kn)) go to 50
	         do 40 k = 1, kn
	            if(dist(j) .lt. nndist(k)) then
	               do 35 k1 = kn-1, k, -1
	                  nndist(k1+1) = nndist(k1)
35	                  pos(k1+1) = pos(k1)
	               nndist(k) = dist(j)
	               pos(k) = j
	               goto 50
	            endif
40	         continue
  	      endif
50	   continue
	   return
	   end


      subroutine twomis(x,p,n,imiss,x0,imiss0,maxit,eps,istart,clust, 
&     nsize,dist,ratio,iter,iworkp,iworkn)
      integer p,n,imiss(p,n),imiss0(n,2),maxit,istart(2),clust(p,2),
&     iworkp(p)
      double precision x(p,n),x0(n,2),eps,dist(p,2)
      integer nsize(2),iworkn(n)
      integer iter,imax
      double precision ratio,dold,dnew
      if(.not.(maxit .lt. 1))goto 23000
      maxit=5
23000 continue
      do 23002 i=1,2
      do 23004 j=1,n
      x0(j,i)=x(istart(i),j)
      imiss0(j,i)=imiss(istart(i),j)
23004 continue
23002 continue
      iter =0
      ratio = 1d1
23006 if(.not.((iter .lt. maxit) .and. (ratio .gt. eps)))goto 23007
      iter=iter+1
      do 23008 i=1,2
      call misdis(x0(1,i),x,p,n,imiss0(1,i),imiss,dist(1,i),iworkp)
      nsize(i)=0
23008 continue
      dnew=0d0
      do 23010 j=1,p 
      if(.not.(dist(j,1) .lt. dist(j,2)))goto 23012
      imax=1
      goto 23013
23012 continue
      imax=2
23013 continue
      nsize(imax)=nsize(imax)+1
      clust(nsize(imax),imax)=j
      dnew=dnew+dist(j,imax)
23010 continue
      if(.not.(dnew .eq. 0d0))goto 23014
      goto 23007
23014 continue
      if(.not.(iter .eq. 1))goto 23016
      dold=dnew*1d1
23016 continue
      ratio=dabs(dnew-dold)/dold
      dold=dnew
      do 23018 i=1,2
      do 23020 j=1,n 
      iworkn(j)=1
23020 continue
      call misave(x,x0(1,i),p,n,iworkn,imiss,clust(1,i),nsize(i))
      do 23022 j=1,n
      if(.not.(iworkn(j) .eq. 2))goto 23024
      imiss0(j,i)=1
      goto 23025
23024 continue
      imiss0(j,i)=0
23025 continue
23022 continue
23018 continue
      goto 23006
23007 continue
      return
      end
