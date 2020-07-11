c This program computes the exact and semiclassical Green's function
c for the Hamiltonian H=(px^2+py^2+pz^2)/2m-Fz=E
c calculations are done using cylindrical coordinates and atomic units
c all subroutines are given
c final version: July 11 2020
c Author: Harindranath Ambalampitiya
c Email: harindranath@huskers.unl.edu
c Free to use with proper Acknowledgments    
      implicit double precision(a-h,o-z)
      complex*16 gree,grescl,gre0
      common /hh/time(2),sj(2),help1,help2,help3
      open(unit=4,file='green_el.out')
c set the initial coordinates
      rho0=2.
      z0=5.
      phi0=0.
c set the final coordinates
      rhofi=5.
      zfi=10.
      phifi=0.
c set the energy and e field
      e0in=0.05
	  field=0.02
      ak=sqrt(2.*abs(e0in))
      pi=2.*asin(1.d0)
c      do 5 phifi=-3.,3.,.5
      do 5 rhofi=-30.,30.,0.05
           
      r=sqrt(zfi*zfi+rhofi*rhofi)
      site=rhofi/r
      cote=zfi/r
      r1=1.-cote*cote
      if(abs(r1).lt.1.d-5)r1=0.d0
      site=sqrt(r1)
c external subroutine for exact Green's functioin 
      call greenel(rho0,z0,phi0,rhofi,zfi,phifi,e0in,field,
     1 gree)
c external subroutine for semiclassical Green's function
      call greenscl(rho0,z0,phi0,rhofi,zfi,phifi,e0in,field,
     1 vrhoi,vzi,grescl)
c free_ particle green's function
      r2=rho0**2+rhofi**2-2.*rho0*rhofi*cos(phifi-phi0)+(zfi-z0)**2
      r=sqrt(r2) 
      t=r/ak
       if(e0in.gt.0.d0)then
      gre0=-exp((0.,1.)*ak*r)/(2.*pi*r)
        else
      gre0=-exp(-ak*r)/(2.*pi*r)
       endif

      greenr=dble(gree)
	  greeni=dimag(gree)
	  grscr=dble(grescl)
	  grsci=dimag(grescl)
      write(4,7) rhofi,greenr,greeni,grscr,grsci
      print 7,rhofi,gree,grescl
    5 continue
    7 format(5e14.5)
    2 format(5x,'ak and gre0',3e14.5)
    1 format(5x,'gree',5e14.5)
      stop
      end
c
      subroutine greenel(rho0,z0,phi0,rhofi,zfi,phifi,e0in,field,
     1 gree)
c for the initial point rho0,z0,phi0,final point rhofi,zfi,phifi, field 
c phifi should be between -pi and pi
c and energy e0in it produces Green function of Slonim and Dalidchik  
c F. I. Dalidchik and V. Z. Slonim,``Strong exchange interaction effects 
c in a homogeneous electric field", JETP \bf{43}, 25 (1976)     
c the sign is opposite to that of S.D.
c force is along positive z axis
      implicit double precision(a-h,o-z)
      complex*16 gree,f2,f2p
      f13=(2.*field)**.333333
      zeta0=-2.*e0in/f13**2
      dr2=rho0*rho0+rhofi*rhofi-2.*rho0*rhofi*cos(phifi-phi0)
     1 +(z0-zfi)**2
      dr=sqrt(dr2)
      zpl=z0+zfi
      zeta1=f13*(zpl-dr)/2.
      zeta2=f13*(zpl+dr)/2.
      call aimat(zeta0-zeta1,v,v1,u,u1,3.d0)
      f1=v
      f1p=v1
      call aimat(zeta0-zeta2,v,v1,u,u1,3.d0)
      f2=u+(0.,1.)*v
      f2p=u1+(0.,1.)*v1
c         print 1,field,zeta0-zeta1,zeta0-zeta2,f1p*f2-f1*f2p,dr
c f1,v,
    1 format(6e13.5)
      gree=(f1p*f2-f1*f2p)/dr/2.
      return
      end
c
      subroutine greenscl(rho0,z0,phi0,rhofi,zfi,phifi,e0in,field,
     1 vrhoi,vzi,gree)
c for the initial point rho0,z0,phi0,final point rhofi,zfi,phifi,field 
c phifi should be between -pi and pi
c and energy e0in it produces initial velocity vrhoi,vzi and the 
      implicit double precision(a-h,o-z)
      complex*16 gree,c1,cs
      dimension vzi(2),vrhoi(2)
      common /hh/time(2),sj(2),help1,help2,help3
      pi=2.*asin(1.d0)
      rhoi=rho0
      zi=z0
      rinitial2=rho0*rho0+z0*z0
      rinitial=sqrt(rinitial2)
      rfinal2=rhofi*rhofi+zfi*zfi
c distance^2
      rd2=rinitial2+rfinal2-2.*(z0*zfi+rho0*rhofi*cos(phifi-phi0)) 
      phii=phi0
      field1=field
      e0in1=e0in
c pure electric field case
c this is k^2/F
        alength=2.*e0in/field
      r1=(1.+2.*z0/alength)*(1.+2.*zfi/alength)
      r4=(rho0**2+rhofi**2-2.*rho0*rhofi*cos(phifi-phi0))/alength**2
c this is s^2 
      r1=r1-r4
      ssm=sqrt(r1)
      q2=1.+(z0+zfi)/alength
c times for two trajectories
      t12=2.*alength/field*(q2-ssm) 
      t22=2.*alength/field*(q2+ssm) 
      time(1)=sqrt(t12)
      time(2)=sqrt(t22)
c 
      q1=rhofi*cos(phifi-phi0)-rho0
      cs=0.
      do 2 i=1,2
c initial velocities
      vzi(i)=(zfi-z0)/time(i)-field*time(i)/2.
      vrhoi(i)=q1/time(i)
c      print 61
      
   
c reduced action
      sj(i)=rd2/(2.*time(i))+time(i)*(field*(z0+zfi)/2.+e0in)
     1 -field*field*time(i)**3/24.        
      
c Green's function
      r1=q2
      if(i.eq.2)r1=-r1  
      r1=r1-ssm
      c1=exp((0.,1.)*sj(i))/sqrt(abs(r1))
      if(i.eq.2)c1=c1/(0.,1.)
    2 cs=cs+c1 
c sign opposite to that in JETP
      gree=-field/(4.*sqrt(2.*ssm)*pi*e0in)*cs
      return
      end
	  
c Airy function
      SUBROUTINE AIMAT(Z,V,V1,U,U1,ZB)
      implicit double precision(a-h,o-z)
c Airy functions and their derivatives, formulae 10.4.2, 10.4.3;
c 10.4.59 - 10.4.64, and 10.4.66, 10.4.67 of Abramowitz
      C1=.3550281
      C2=.2588194
      ZA=ABS(Z)  
      IF(ZA-ZB)1,1,2
    1 A=1.          
      F=1.          
      Z3=Z**3       
      A1=Z*Z/2.     
      F1=A1         
      B=Z           
      G=Z           
      B1=1.         
      G1=1.         
      DO 3 K=1,200  
      R1=3*K        
      A=A*Z3/R1/(R1-1.)
      F=F+A            
      B=B*Z3/R1/(R1+1.)
      G=G+B            
      B1=B1*Z3/R1/(R1-2.)
      G1=G1+B1           
      IF(K.EQ.1)GOTO 3 
      A1=A1*Z3/(R1-1.)/(R1-3.)
      F1=F1+A1                
      IF(ABS(A).LT.1.E-6)GOTO 4
    3 CONTINUE                 
    4 V=C1*F-C2*G              
      U=1.7320508*(C1*F+C2*G)  
      V1=C1*F1-C2*G1           
      U1=1.7320508*(C1*F1+C2*G1)
      GOTO 5                    
    2 DZ=ZA*SQRT(ZA)*2./3.      
      A=1.                      
      S=1.                      
      A1=5/72./DZ               
      S1=A1                     
      B=1.                      
      B1=-1.4*A1                
      P=1.                      
      P1=B1                     
      IND=1                     
      DO 6 K=1,10               
      IF(Z.LT.0.)IND=-IND       
      R=12*K-5                  
      A=A1*R*(R+4.)/144./K/DZ   
      A1=A*(R+6.)*(R+10.)/72./(2.*K+1.)/DZ
      B=-(R+6.)/(R+4.)*A                  
      B1=-(R+12.)/(R+10.)*A1              
      IF(ABS(A1).GT.ABS(A).OR.ABS(B1).GT.ABS(B))GOTO 7
      S=S+IND*A                                       
      S1=S1+IND*A1                                    
      P=P+IND*B                                       
    6 P1=P1+IND*B1                                    
    7 Q=ZA**(-.25)/1.7724539                          
      SZ=SQRT(ZA)                                     
      IF(Z)8,9,9                                      
    9 F=EXP(-DZ)                                      
      V=.5*Q*F*(S-S1)                                 
      V1=-.5*SZ*Q*F*(P-P1)                            
      U=Q/F*(S+S1)                                    
      U1=Q/F*SZ*(P+P1)                                
      GOTO 5                                          
    8 F=SIN(DZ+.7853982)                              
      G=COS(DZ+.7853982)                              
      V=Q*(F*S-G*S1)                                  
      V1=-SZ*Q*(G*P+F*P1)                             
      U=Q*(G*S+F*S1)                                  
      U1=Q*SZ*(F*P-G*P1)                              
    5 RETURN                                          
      END        