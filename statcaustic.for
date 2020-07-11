c This program computes an ensemble of trajectories 
c for the Hamiltonian H=(px^2+py^2+pz^2)/2m-Fz=E
c store the data in multiple files
c also compute the locus of the cuastic
c final version: July 11 2020
c Author: Harindranath Ambalampitiya
c Email: harindranath@huskers.unl.edu
c Free to use with proper Acknowledgments      
	  implicit double precision(a-h,o-z)
      character string*3
      
	  open(14,file='caustic.out')
      open(18,file='alltrajs.out')
	  pi=2.*asin(1.d0)
c energy and field
	  en=0.05
	  fld=0.02
c initial points
	  x10=2.
	  z10=5.
	  p=sqrt(2.*(en+fld*z10))
	  
	  k=0
c loop over ejection angle
      do 1 theta=-pi,pi,pi/40.
	  z1dot=p*cos(theta)
	  x1dot=p*sin(theta)
	  k=k+1
	  kk=100+k
	  write(string,fmt='(I3)')kk
      print*,kk
      open(unit=13,file="stat_cstcs"//string//".out",status='replace')
      do 2 t2nano=0.,1000.,0.5
	  time=t2nano
	  xt=x1dot*time+x10
	  zt=.5*fld*time*time+z1dot*time+z10
      write(13,*),time,xt,zt
c      write(18,*),time,xt,zt
 2    continue
c      write(18,*)
 1    continue
c plot the cuastic surface 
      do 21 xc=-50.,50.,0.5
	  a1=fld*fld*(xc-x10)**2/(4.*en**2)
	  a2=1.+z10*fld/en
	  zc=a1*en/(a2*fld)-en/fld
      write(14,*),xc,zc,-fld*zc
 21   continue
      end