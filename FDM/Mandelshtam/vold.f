
      program filter_diagonalization
c!
c! Vladimir Mandelshtam, last modified Jan. 1998
c! Please don't distribute without permission.
c! The author is not responsible for any damage caused by this code.
c! Report any bugs found (mandelsh@uci.edu).
c! Any suggestions are welcome.
c! Please let me know how it works.
c!
c!
c!
c!  This program for a signal c_n=C(t0+tau*n) n=0,1,...,Nsig
c!  obtains the complex frequencies w_k and amplitudes d_k
c!  by fitting to the form
c!                          c_n= sum_k d_k exp(-inw_k)         (1)
c!
c!  The coefficients are phase corrected by d_k=d_k*exp(i*theta)
c!
c!  The spectrum is computed using the formula:
c!
c!            F(w)=sum_k d_k/(w-w_k+i*Gamm)       (2)
c!
c!     Cheating:    replace w_k in (2) by Re{w_k}+i*cheat*Im{w_k}
c!c
c!     If Im{w_k}>0 it is replaced by -Im{w_k} in (2).
c!
c!
c! Remark: for diagonalization of a complex general matrix I use cg(...
c!         which is an old routine and might cause some problems at a high
c!         optimization level on some machines
c!
c! Remark: Try first to use Nb0=-1, error1=error2=1000, 
c!         Gcut=10000, Gamm=0.000001, ispec=1, cheat=1
c!
c! Remark: Try to use such frequency windows for which Nb0 becomes >3 and <300
c!
c! Remark: If there is no noise, i.e. (1) is exact, the choice for Nb0 is crucial
c!             In such a case the best convergence will be achieved when Nb0 is 
c!             just slightly greater than Nb (the actual rank of the U0 matrix).
c!
      implicit real*8(a-h,o-z)
      parameter (Nbmax=300,Nsigmax=270000)
      complex*16 U(Nbmax,Nbmax,-1:1),uk(Nbmax),
     &     w_minus(Nbmax),d_minus(Nbmax),
     &     w_plus(Nbmax),d_plus(Nbmax),
     &     zz(Nbmax,Nbmax),A(Nbmax,Nbmax),uu,phase_corr,
     &     f(Nbmax,-1:1,2),xi,z(Nbmax)
      real*8 err_minus(Nbmax),err_plus(Nbmax)
      complex*16 coef(0:Nsigmax) 
      character*30 par_minus,par_plus,signal,AbsSp,ReSp_minus,ReSp_plus
      character*30 fname,fname1,fname2
!
!       This arrays are needed for cg.f only    
!
      real*8 Zr(Nbmax,Nbmax,2),wr(Nbmax,2)
      open(5,file = 'input.file')
      
!
      read(5,*) signal, idat ! the input data file with the signal given as a complex column vector.
                             ! idat>0 read real signal 
                             ! idat<0 read complex signal 
                             ! |idat|=1 read c_n
                             ! |idat|=2 read t_n,c_n
                             ! idat=-3 read t_n,Re(c_n),Im(c_n)
                             ! idat=-4 read Re(c_n),Im(c_n)


c      read(5,*) Nsig1max,delt1  ! the time step should be given in the units of 1/[frequency units]
                                ! here tau is included in the file with the signal

      read(5,*) t0, theta     ! t0 --- the time delay (we assume that the first time increment is t0
                              ! theta ---- the overall phase correction   exp(i*theta)
      print*,'t0,theta', t0, theta 
      read(5,*) ispec         ! ispec=0  ---- Do spectral analysis with FD and FD/DFT hybrid method
                              ! ispec=1  ---- No DFT
                              ! ispec=-1 ---- DFT only
      print*,'ispec',ispec 

      read(5,*) par_plus      ! spectral parameters from the U_plus diagonalization

      read(5,*) par_minus     ! spectral parameters from the U_minus diagonalization

      read(5,*) AbsSp         ! magnitude spectrum constructed out of U_minus parameters

      read(5,*) ReSp_minus    ! Absorption spectrum constructed out of U_minus parameters
                              ! theta=pi/2 will interchange absorption and dispersion

      read(5,*) ReSp_plus     ! Absorption spectrum constructed out of U_plus parameters

      read(5,*) wmin_old,wmax_old  ! the frequency window. 
                                   ! The units are such that the Nyquist range = 2*pi/tau

      read(5,*) Nsig         ! numner of signal data points to be used in the analysis    

      read(5,*) Nskip        ! number of the signal points to skip before the analysis    

      read(5,*) Nb0          ! The size of the window basis to be used (Nb0<0 -- automatic choice).

      read(5,*) smin         ! the minimum eigenvalue of the overlap matrix used to do the inversion.
                             ! smin=1d-15 for quadruple precision and smin=1d-8 for double precision.
      read(5,*) error1,error2 ! error1 is used to check the errors for the U_plus eigenvalues by
                              ! calculating ||U_minus B_plus_k - u_plus_k B_plus_k||
                              ! error2 is used to check the errors for the U_minus eigenvalues by
                              ! calculating ||w_minus_k - w_plus_k||
                              ! use error1=error2=100000 to avoid missing eigenvalues

      read(5,*) Npower,Gamm,Gcut,cheat  ! Npower --- Number of frequency points to plot the spectrum P(w) 
                                   ! Gamm ----  the smoothing width to construct the spectrum
                              ! Gcut ---- maximum width for a pole to be used to construct the spectrum
                          ! cheat ---- multiply all widths by cheat (be very careful or use cheat=1 !!!)
      print*,'read input'
      if(Nsig.gt.Nsigmax) stop 'increase Nsigmax'
      open(7,file=par_plus)     
      open(9,file=par_minus)    
      open(15,file=signal)
      open(16,file='sizefl')
      if(ispec.le.1) then
         open(18,file=ReSp_plus)
         open(19,file=AbsSp)
         open(22,file=ReSp_minus)
         open(27,file='Spectra')
      endif
      write(6,*) 'read the signal from disk'
      read(16,*) Nsig1max,delt1
      delt=delt1
c      if(Nsig.gt.Nsig1max-1) Nsig=Nsig1max-1
      do n=1,Nskip
         read(15,*)
      enddo
      c2=0d0
      do n=0,Nsig
         print*,'n',n
         if(idat.eq.-3) read(15,*,end=100) tt,c1,c2
         if(idat.eq.-4) read(15,*,end=100) c1,c2
         if(idat.eq.1)  read(15,*,end=100) c1
         if(idat.eq.2)  read(15,*,end=100) tt,c1
c         if(idat.eq.2)  read(15,*) tt,c1
         if(idat.eq.-1) read(15,*,end=100) coef(n)
         if(idat.eq.-2) read(15,*,end=100) tt,coef(n)
         if(idat.ne.-1.and.idat.ne.-2) coef(n)=dcmplx(c1,c2)
c         print*,'Time',tt,'Sig',c1
      enddo
      goto 101
 100  Nsig=n-1
 101  M=(Nsig-2)/2
      Nsig=2*M+2
      write(6,*) 'M=',M,' which means that', 2*M+3, 'c_n are used'
      pi=dacos(-1d0)
      xi=(0d0,1d0)
      phase_corr=cdexp(dcmplx(0d0,theta))
      write(6,*) 'test',exp(dcmplx(0d0,theta)),phase_corr
      wmin=wmin_old*delt
      if(wmin.lt.-pi) wmin=-pi
      wmax=wmax_old*delt
      if(wmax.gt.pi) wmax=pi
      wmin_old=wmin/delt
      wmax_old=wmax/delt
      If (Nb0.lt.0) then
         write(6,*) 'Use the FFT grid'
         Nb0=((wmax-wmin)/2/pi)*(M+1)
         if(Nb0.lt.1) Nb0=1
      endif
      if(Nb0.gt.Nbmax) stop 'Nb0>Nbmax'
      do i=1,Nb0
         z(i)=wmin+(i-0.5d0)*(wmax-wmin)/Nb0
         z(i)=cdexp(xi*z(i)) ! this is actually 1/z_j as defined in the paper
      enddo
      write(6,*) 'Nb0=',Nb0, ' wmin,wmax=',wmin_old,wmax_old
!
!___________________________________________________________________________________
!
      call       FD(Nsig,coef,delt,wmin,wmax,Nb0,error1,error2,
     &     smin,z,Nb,t0,phase_corr,ispec,Npower,cheat,Gamm,Gcut,
     &     w_minus,d_minus,err_minus,w_plus,d_plus,err_plus,
     &     U,uk,f,zz,A,A,Zr,wr)
      stop
      end
!
!___________________________________________________________________________________

      subroutine FD(Nsig,coef,delt,wmin,wmax,Nb0,error1,error2,
     &     smin,z,Nb,t0,phase_corr,ispec,Npower,cheat,Gamm,Gcut,
     &     w_minus,d_minus,err_minus,w_plus,d_plus,err_plus,
     &     U,uk,f,zz,A,Ar,Zr,wr)
!
      implicit real*8(a-h,o-z)
c      implicit none
      complex*16 U(Nb0,Nb0,-1:1),uk(Nb0),coef(0:Nsig),
     &     w_minus(Nb0),d_minus(Nb0),
     &     w_plus(Nb0),d_plus(Nb0),
     &     zz(Nb0,Nb0),A(Nb0,Nb0),
     &     f(Nb0,2,-1:1),z(Nb0),
     &     Z1,Z2,phase_corr,Power,xi,ss,uu,cc
      real*8 err_minus(Nb0),err_plus(Nb0),Gcut,cheat,omr,omi
      integer Nb,Nb0,i,j,n,k,k1,M,Nsig,ispec,ierr,N0,Npower
      real*8 rho,pi,Omega,waver,smin,Gamm,t0,delt,yiter,error1,
     &     error2,wmin,wmax
!
!       This arrays are needed for cg.f only    
!
      real*8 Ar(Nb0,Nb0,2),Zr(Nb0,Nb0,2),wr(Nb0,2)
c      equivalence (A,Ar)
!

      M=(Nsig-2)/2
      pi=dacos(-1d0)
      xi=(0d0,1d0)
      if(ispec.lt.0) goto 399
      write(6,*) 'construction of small U0,U1 and U2 matrices'
      do i=1,Nb0
         f(i,1,0)=(0d0,0)
         f(i,2,0)=(0d0,0)
         Z1=(1d0,0)
         Z2=z(i)**M
         ss=coef(1)
         do k=1,M
            Z1=Z1*z(i)
            Z2=Z2*z(i)
            f(i,1,0)=f(i,1,0)+coef(k+1)*Z1
            f(i,2,0)=f(i,2,0)+coef(k+M+1)*Z2
            ss=ss+(k+1)*(coef(k+1)*Z1)+(M-k+1)*(coef(k+M+1)*Z2)
         enddo
         f(i,1,1)=f(i,1,0)/z(i)-coef(2)+coef(M+2)*Z1
         f(i,2,1)=f(i,2,0)/z(i)-coef(M+2)*Z1+coef(2*M+2)*Z2
         f(i,1,-1)=(f(i,1,0)+coef(1)-coef(M+1)*Z1)*z(i)
         f(i,2,-1)=(f(i,2,0)-coef(2*M+1)*Z2+coef(M+1)*Z1)*z(i)
         U(i,i,0)=ss
         U(i,i,1)=(ss-coef(1)-f(i,1,0)+f(i,2,0))/z(i)
     &        +coef(2*M+2)*Z2
         U(i,i,-1)=(ss+coef(1)-2*coef(M+1)*Z1+f(i,1,0)
     &        -f(i,2,0))*z(i)+coef(0)
      enddo
      do i=1,Nb0
         do j=1,i-1
            cc=z(i)/z(j)
            U(i,j,0)=coef(1)+(f(j,1,0)-cc**(M+1)*f(j,2,0)-
     &           cc*f(i,1,0)+cc**(-M)*f(i,2,0))/(1-cc)
            U(j,i,0)=U(i,j,0)
            U(i,j,-1)=coef(0)+(f(j,1,-1)-cc**(M+1)*f(j,2,-1)-    ! matrix of 1/U
     &           cc*f(i,1,-1)+cc**(-M)*f(i,2,-1))/(1-cc)
            U(j,i,-1)=U(i,j,-1)
            U(i,j,1)=coef(2)+(f(j,1,1)-cc**(M+1)*f(j,2,1)-    ! matrix of U
     &           cc*f(i,1,1)+cc**(-M)*f(i,2,1))/(1-cc)
            U(j,i,1)=U(i,j,1)
         enddo
      enddo
      write(6,*) 'Diagonalization of the small U0 matrix'
c_________________________________________________________________
c
c  This may be replaced by another diagonalization routine
c
c
      do i=1,Nb0
         do j=1,Nb0
            Ar(i,j,1)=real(U(i,j,0))
            Ar(i,j,2)=dimag(U(i,j,0))
         enddo
      enddo
      call CG(Nb0,Nb0,Ar(1,1,1),Ar(1,1,2),
     &     wr(1,1),wr(1,2),1,Zr(1,1,1),Zr(1,1,2),
     &     f(1,1,-1),IERR)
      do i=1,Nb0
         w_minus(i)=dcmplx(wr(i,1),wr(i,2))
         do j=1,Nb0
            zz(i,j)=dcmplx(Zr(i,j,1),Zr(i,j,2))
         enddo
      enddo
c_________________________________________________________________

      waver=0d0
      do n=1,Nb0
         waver=waver+cdabs(w_minus(n))
      enddo
      waver=waver/Nb0
c
      write(6,*) 'Construct U0^{-1/2} U^p U0^{-1/2} using SVD of S'
      Nb=0
      do i=1,Nb0
         if(cdabs(w_minus(i)).gt.smin*waver) then
            Nb=Nb+1
            cc=(0,0d0)
            do k=1,Nb0
               cc=cc+zz(k,i)**2
            enddo
            cc=(1d0,0)/cdsqrt(cc*w_minus(i))
            do k=1,Nb0
               zz(k,Nb)=cc*zz(k,i)
            enddo
         endif
      enddo
      write(6,*) 'Nb0=', Nb0,' Nb=',Nb,
     &     'are left after analysing the eigenvalues of U0',
     &     ' using the smin<u0 criterion.'
      do j=1,Nb
         do n=1,Nb0
            A(j,n)=(0d0,0)
            do k=1,Nb0
               A(j,n)=A(j,n)+U(k,n,-1)*zz(k,j)
            enddo
         enddo
      enddo
      do j=1,Nb
         do i=1,Nb
            U(j,i,-1)=(0d0,0)
            do k=1,Nb0
               U(j,i,-1)=U(j,i,-1)+A(j,k)*zz(k,i)
            enddo
         enddo
      enddo
      do j=1,Nb
         do n=1,Nb0
            A(j,n)=(0d0,0)
            do k=1,Nb0
               A(j,n)=A(j,n)+U(k,n,1)*zz(k,j)
            enddo
         enddo
      enddo
      do j=1,Nb
         do i=1,Nb
            U(j,i,1)=(0d0,0)
            do k=1,Nb0
               U(j,i,1)=U(j,i,1)+A(j,k)*zz(k,i)
            enddo
         enddo
      enddo
c
      write(6,*) 'Diagonalization of the small U_plus-matrix'
c_________________________________________________________________
c
c  This may be replaced by another diagonalization routine
c
c
      do i=1,Nb
         do j=1,Nb
            Ar(i,j,1)=real(U(i,j,1))
            Ar(i,j,2)=dimag(U(i,j,1))
         enddo
      enddo
      call CG(Nb0,Nb,Ar(1,1,1),Ar(1,1,2),
     &     wr(1,1),wr(1,2),1,Zr(1,1,1),Zr(1,1,2),
     &     f(1,1,-1),IERR)
      do k=1,Nb
         uk(k)=dcmplx(wr(k,1),wr(k,2))
         w_plus(k)=(0,1d0)*cdlog(uk(k))/delt    ! get the frequency in the correct units
         do j=1,Nb
            A(k,j)=dcmplx(Zr(k,j,1),Zr(k,j,2))
         enddo
      enddo
c_________________________________________________________________

      do k=1,Nb
         cc=(0,0d0)
         do j=1,Nb
            cc=cc+A(j,k)**2
         enddo
         cc=(1d0,0)/cdsqrt(cc)
         do j=1,Nb
            A(j,k)=cc*A(j,k)
         enddo
c Compute the errors
         err_plus(k)=0d0
         do j=1,Nb
            cc=A(j,k)/uk(k)
            do i=1,Nb
               cc=cc-U(i,j,-1)*A(i,k)
            enddo
            err_plus(k)=err_plus(k)+cc*dconjg(cc)
         enddo
         err_plus(k)=dsqrt(err_plus(k)*dsqrt(dfloat(Nb)))
      enddo
      do i=1,Nb0
         do k=1,Nb
            U(i,k,0)=(0d0,0)
            do j=1,Nb
               U(i,k,0)=U(i,k,0)+A(j,k)*zz(i,j)
            enddo
         enddo
      enddo
      write(6,*) 'compute the coefficients'
      do 87 k=1,Nb
         d_plus(k)=(0,0d0)
         if(err_plus(k).gt.error1) goto 87
         rho=cdabs(uk(k))
         if(rho.gt.0.999999999d0) then
            uu=(1d0,0)/uk(k)
         else
            uu=rho/uk(k)
         endif
         Z1 = coef(M+1)*uu
         Z2 = coef(2*M+1)*uu
         do n=M,2,-1
            Z1 = (Z1+coef(n))*uu
            Z2 = (Z2+coef(M+n))*uu
         enddo
         do j=1,Nb0
            cc=uu/z(j)
            d_plus(k)=d_plus(k)+U(j,k,0)*
     &           (coef(1)+(f(j,1,0)-cc**(M+1)*f(j,2,0)-
     &           cc*Z1+z(j)**M*Z2)/(1-cc))
         enddo
         if(rho.gt.0.999999999d0) then
            d_plus(k)=d_plus(k)/(M+1)
         else
            d_plus(k)=d_plus(k)/((1d0,0)-rho**(M+1))*((1d0,0)-rho)
         endif
         d_plus(k)=d_plus(k)**2*phase_corr
     &        *cdexp((0,1d0)*w_plus(k)*(t0+delt)) !to adjust the t=0
 87   continue
      call cpiksrt(Nb,w_plus,d_plus,err_plus)
      write(9,79) 
      do k=1,Nb
         if(err_plus(k).lt.error1)
     &        write(9,13) real(w_plus(k)),dimag(w_plus(k)),
     &        real(d_plus(k)),dimag(d_plus(k)),err_plus(k)
      enddo
      close(9)
c
      if(ispec.le.1) then
         write(6,*) 'Compute the spectrum'
         do i=0,Npower
            Power=(0d0,0)
            Omega=wmin+i*((wmax-wmin)/Npower)
            Omega=Omega/delt
            do k=1,Nb 
               cc=dcmplx(dreal(w_plus(k)),cheat*dimag(w_plus(k)))
               if(err_plus(k).lt.error1.and.
     &        dimag(w_plus(k)).gt.-Gcut) then
c     &              .and.dreal(w_plus(k)).lt.wmax/delt.and.
c     &              dreal(w_plus(k)).gt.wmin/delt) then
                  if(dimag(cc).gt.0d0) then
                     Power=Power+d_plus(k)/
     &                    (Omega-dconjg(cc)+xi*Gamm)
                  else
                     Power=Power+d_plus(k)/
     &                    (Omega-cc+xi*Gamm)
                  endif
               endif
            enddo
            write(18,*) Omega,real(Power)
         enddo
         close(18)
      endif
      write(6,*) 'Diagonalization of the small U_minus-matrix'
c_________________________________________________________________
c
c  This may be replaced by another diagonalization routine
c
c
      do i=1,Nb
         do j=1,Nb
            Ar(i,j,1)=real(U(i,j,-1))
            Ar(i,j,2)=dimag(U(i,j,-1))
         enddo
      enddo
      call CG(Nb0,Nb,Ar(1,1,1),Ar(1,1,2),
     &     wr(1,1),wr(1,2),1,Zr(1,1,1),Zr(1,1,2),
     &     f(1,1,-1),IERR)
      do k=1,Nb
         uk(k)=1/dcmplx(wr(k,1),wr(k,2))
         w_minus(k)=(0,1d0)*cdlog(uk(k))/delt    ! get the frequency in the correct units
         do j=1,Nb
            A(k,j)=dcmplx(Zr(k,j,1),Zr(k,j,2))
         enddo
      enddo
c_________________________________________________________________
      do k=1,Nb
         cc=(0,0d0)
         do j=1,Nb
            cc=cc+A(j,k)**2
         enddo
         cc=(1d0,0)/cdsqrt(cc)
         do j=1,Nb
            A(j,k)=cc*A(j,k)
         enddo
c Compute the errors
         err_minus(k)=99999999999d0
         do k1=1,Nb
            if(cdabs(w_plus(k1)-w_minus(k)).lt.err_minus(k)) 
     &           err_minus(k)=cdabs(w_plus(k1)-w_minus(k))
         enddo
c         err_minus(k)=0d0
c         do j=1,Nb
c            cc=A(j,k)*uk(k)
c            do i=1,Nb
c               cc=cc-U(i,j,1)*A(i,k)
c            enddo
c            err_minus(k)=err_minus(k)+cc*dconjg(cc)
c         enddo
c         err_minus(k)=dsqrt(err_minus(k)*dsqrt(dfloat(Nb)))
      enddo
      do i=1,Nb0
         do k=1,Nb
            U(i,k,0)=(0d0,0)
            do j=1,Nb
               U(i,k,0)=U(i,k,0)+A(j,k)*zz(i,j)
            enddo
         enddo
      enddo
      write(6,*) 'compute the coefficients'
      do 88 k=1,Nb
         d_minus(k)=(0,0d0)
         if(err_minus(k).gt.error2) goto 88
         rho=cdabs(uk(k))
         if(rho.gt.0.999999999d0) then
            uu=(1d0,0)/uk(k)
         else
            uu=rho/uk(k)
         endif
         Z1 = coef(M+1)*uu
         Z2 = coef(2*M+1)*uu
         do n=M,2,-1
            Z1 = (Z1+coef(n))*uu
            Z2 = (Z2+coef(M+n))*uu
         enddo
         do j=1,Nb0
            cc=uu/z(j)
            d_minus(k)=d_minus(k)+U(j,k,0)*
     &           (coef(1)+(f(j,1,0)-cc**(M+1)*f(j,2,0)-
     &           cc*Z1+z(j)**M*Z2)/(1-cc))
         enddo
         if(rho.gt.0.999999999d0) then
            d_minus(k)=d_minus(k)/(M+1)
         else
            d_minus(k)=d_minus(k)/
     &           ((1d0,0)-rho**(M+1))*((1d0,0)-rho)
         endif
         d_minus(k)=d_minus(k)**2*phase_corr
     &        *cdexp((0,1d0)*w_minus(k)*(t0+delt)) !to adjust the t=0
 88      continue
      call cpiksrt(Nb,w_minus,d_minus,err_minus)
      write(7,79) 
 79   format(8x,'Re w',17x,'Im w',14x,'Re d',13x,'Im d',8x,'error')
      if(ispec.eq.0) write(6,*) 'Split the signal'
      do k=1,Nb
         if(err_minus(k).lt.error2) then
		omr=5305.16d0*real(w_minus(k))
		omi=5305.16d0*dimag(w_minus(k))
	write(27,23)omr,abs(omi),cdabs(d_minus(k))
            write(7,13) omr,omi,
     &        real(d_minus(k)),dimag(d_minus(k)),err_minus(k)
            if(d_minus(k).ne.(0,0d0).and.ispec.eq.0) then    !subtruct the signal part
               if(dimag(w_minus(k)).gt.0d0) then
                  uu=cdexp(-(0,1d0)*dconjg(w_minus(k))*delt)
               else
                  uu=cdexp(-(0,1d0)*w_minus(k)*delt)
               endif
               Z1=d_minus(k)/uu/phase_corr
               do n=0,2*M+2
                  Z1=Z1*uu
                  coef(n)=coef(n)-Z1
               enddo
            endif
         endif
 13      format(E22.16,3E17.10,E9.2)
 23      format(2x,F10.2,3x,F10.2,3x,E10.4)
      enddo
      close(7)
c     
 399  if(ispec.le.1) then
         write(6,*) 'Compute the spectrum'
         if(ispec.lt.0) write(6,*) 'Only DFT is used'
         N0=M/5
         if(ispec.lt.1) then    !before DFT use a window
            do n=0,2*M+2
               if (n.gt.N0) then
                  yiter = dfloat(n-N0)/dfloat(2*M+3-N0)
                  coef(n)=coef(n)*dexp(yiter)*(1D0-yiter)
               endif
            enddo
         endif
         do i=0,Npower
            Omega=wmin+i*((wmax-wmin)/Npower)
            Power=0d0
            if(ispec.lt.1) then !do DFT
               uu=cdexp((0,1d0)*Omega)
               Z1=-delt*(0,1d0)*cdexp(dcmplx(0d0,t0*Omega))
               if(t0.gt.0d0) then
                  Power=-delt*(0,0.5d0)*coef(0) 
               else
                  Power=-delt*coef(0)*cdexp(dcmplx(0d0,t0*Omega)) 
               endif
               do n=1,2*M+2
                  Z1=Z1*uu
                  Power=Power+coef(n)*Z1
               enddo
            endif
            Power=Power*phase_corr ! to adjust the constant phase correction
            Omega=Omega/delt
            write(97,*) 5305.16*Omega,dreal(Power)
            do k=1,Nb 
               cc=dcmplx(dreal(w_minus(k)),cheat*dimag(w_minus(k)))
               if(err_minus(k).lt.error2.and.ispec.ge.0.and.
     &              dimag(w_minus(k)).gt.-Gcut) then
c     &              .and.dreal(w_minus(k)).lt.wmax/delt.and.
c     &              dreal(w_minus(k)).gt.wmin/delt) 
                  if(dimag(cc).gt.0) then
                     Power=Power+d_minus(k)/
     &                    (Omega-dconjg(cc)+xi*Gamm)
                  else
                     Power=Power+d_minus(k)/
     &                    (Omega-cc+xi*Gamm)
                  endif
               endif
            enddo
            write(19,*) 5305.16*Omega,cdabs(Power)
            write(22,*) 5305.16*Omega,real(Power)
         enddo
      endif
c
      return
      end
c
      subroutine cpiksrt(Nb,wr,wi,err)
      implicit real*8(a-h,o-z)
      complex*16 wr(Nb),wi(Nb),war,wai
      real*8 err(Nb),wae
      do j=2,Nb
         war=wr(j)
         wai=wi(j)
         wae=err(j)
         do i=j-1,1,-1
            if(real(wr(i)).ge.real(war))go to 10
            wr(i+1)=wr(i)
            wi(i+1)=wi(i)
            err(i+1)=err(i)
         enddo
         i=0
 10      continue
         wr(i+1)=war
         wi(i+1)=wai
         err(i+1)=wae
      enddo
      return
      end
c  
c      
      SUBROUTINE CG(NM,N,AR,AI,WR,WI,MATZ,ZR,ZI,FV,IERR)
C
      INTEGER N,NM,IS1,IS2,IERR,MATZ
      real*8 AR(NM,N),AI(NM,N),WR(N),WI(N),ZR(NM,N),ZI(NM,N),
     X       FV(N,3)
C
C     THIS SUBROUTINE CALLS THE RECOMMENDED SEQUENCE OF
C     SUBROUTINES FROM THE EIGENSYSTEM SUBROUTINE PACKAGE (EISPACK)
C     TO FIND THE EIGENVALUES AND EIGENVECTORS (IF DESIRED)
C     OF A COMPLEX GENERAL MATRIX.
C
C     ON INPUT
C
C        NM  MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL
C        ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C        DIMENSION STATEMENT.
C
C        N  IS THE ORDER OF THE MATRIX  A=(AR,AI).
C
C        AR  AND  AI  CONTAIN THE REAL AND IMAGINARY PARTS,
C        RESPECTIVELY, OF THE COMPLEX GENERAL MATRIX.
C
C        MATZ  IS AN INTEGER VARIABLE SET EQUAL TO ZERO IF
C        ONLY EIGENVALUES ARE DESIRED.  OTHERWISE IT IS SET TO
C        ANY NON-ZERO INTEGER FOR BOTH EIGENVALUES AND EIGENVECTORS.
C
C     ON OUTPUT
C
C        WR  AND  WI  CONTAIN THE REAL AND IMAGINARY PARTS,
C        RESPECTIVELY, OF THE EIGENVALUES.
C
C        ZR  AND  ZI  CONTAIN THE REAL AND IMAGINARY PARTS,
C        RESPECTIVELY, OF THE EIGENVECTORS IF MATZ IS NOT ZERO.
C
C        IERR  IS AN INTEGER OUTPUT VARIABLE SET EQUAL TO AN ERROR
C           COMPLETION CODE DESCRIBED IN THE DOCUMENTATION FOR COMQR
C           AND COMQR2.  THE NORMAL COMPLETION CODE IS ZERO.
C
C        FV is a TEMPORARY STORAGE ARRAY.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IF (N .LE. NM) GO TO 10
      IERR = 10 * N
      GO TO 50
C
   10 CALL  CBAL(NM,N,AR,AI,IS1,IS2,FV(1,1))
      CALL  CORTH(NM,N,IS1,IS2,AR,AI,FV(1,2),FV(1,3))
      IF (MATZ .NE. 0) GO TO 20
C     .......... FIND EIGENVALUES ONLY ..........
      CALL  COMQR(NM,N,IS1,IS2,AR,AI,WR,WI,IERR)
      GO TO 50
C     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........
   20 CALL  COMQR2(NM,N,IS1,IS2,FV(1,2),FV(1,3),AR,AI,WR,WI,ZR,ZI,IERR)
      IF (IERR .NE. 0) GO TO 50
      CALL  CBABK2(NM,N,IS1,IS2,FV(1,1),N,ZR,ZI)
   50 RETURN
      END
      SUBROUTINE CBABK2(NM,N,LOW,IGH,SCALE,M,ZR,ZI)
C
      INTEGER I,J,K,M,N,II,NM,IGH,LOW
      real*8 SCALE(N),ZR(NM,M),ZI(NM,M)
      real*8 S
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE
C     CBABK2, WHICH IS A COMPLEX VERSION OF BALBAK,
C     NUM. MATH. 13, 293-304(1969) BY PARLETT AND REINSCH.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 315-326(1971).
C
C     THIS SUBROUTINE FORMS THE EIGENVECTORS OF A COMPLEX GENERAL
C     MATRIX BY BACK TRANSFORMING THOSE OF THE CORRESPONDING
C     BALANCED MATRIX DETERMINED BY  CBAL.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY  CBAL.
C
C        SCALE CONTAINS INFORMATION DETERMINING THE PERMUTATIONS
C          AND SCALING FACTORS USED BY  CBAL.
C
C        M IS THE NUMBER OF EIGENVECTORS TO BE BACK TRANSFORMED.
C
C        ZR AND ZI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE EIGENVECTORS TO BE
C          BACK TRANSFORMED IN THEIR FIRST M COLUMNS.
C
C     ON OUTPUT
C
C        ZR AND ZI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE TRANSFORMED EIGENVECTORS
C          IN THEIR FIRST M COLUMNS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IF (M .EQ. 0) GO TO 200
      IF (IGH .EQ. LOW) GO TO 120
C
      DO 110 I = LOW, IGH
         S = SCALE(I)
C     .......... LEFT HAND EIGENVECTORS ARE BACK TRANSFORMED
C                IF THE FOREGOING STATEMENT IS REPLACED BY
C                S=1.0D0/SCALE(I). ..........
         DO 100 J = 1, M
            ZR(I,J) = ZR(I,J) * S
            ZI(I,J) = ZI(I,J) * S
  100    CONTINUE
C
  110 CONTINUE
C     .......... FOR I=LOW-1 STEP -1 UNTIL 1,
C                IGH+1 STEP 1 UNTIL N DO -- ..........
  120 DO 140 II = 1, N
         I = II
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 140
         IF (I .LT. LOW) I = LOW - II
         K = SCALE(I)
         IF (K .EQ. I) GO TO 140
C
         DO 130 J = 1, M
            S = ZR(I,J)
            ZR(I,J) = ZR(K,J)
            ZR(K,J) = S
            S = ZI(I,J)
            ZI(I,J) = ZI(K,J)
            ZI(K,J) = S
  130    CONTINUE
C
  140 CONTINUE
C
  200 RETURN
      END
      SUBROUTINE COMQR2(NM,N,LOW,IGH,ORTR,ORTI,HR,HI,WR,WI,ZR,ZI,IERR)
C
      INTEGER I,J,K,L,M,N,EN,II,JJ,LL,NM,NN,IGH,IP1,
     X        ITN,ITS,LOW,LP1,ENM1,IEND,IERR
      real*8 HR(NM,N),HI(NM,N),WR(N),WI(N),ZR(NM,N),ZI(NM,N),
     X       ORTR(IGH),ORTI(IGH)
      real*8 SI,SR,TI,TR,XI,XR,YI,YR,ZZI,ZZR,NORM,TST1,TST2,
     X       PYTHAG
C
C     THIS SUBROUTINE IS A TRANSLATION OF A UNITARY ANALOGUE OF THE
C     ALGOL PROCEDURE  COMLR2, NUM. MATH. 16, 181-204(1970) BY PETERS
C     AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).
C     THE UNITARY ANALOGUE SUBSTITUTES THE QR ALGORITHM OF FRANCIS
C     (COMP. JOUR. 4, 332-345(1962)) FOR THE LR ALGORITHM.
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS
C     OF A COMPLEX UPPER HESSENBERG MATRIX BY THE QR
C     METHOD.  THE EIGENVECTORS OF A COMPLEX GENERAL MATRIX
C     CAN ALSO BE FOUND IF  CORTH  HAS BEEN USED TO REDUCE
C     THIS GENERAL MATRIX TO HESSENBERG FORM.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
C          SUBROUTINE  CBAL.  IF  CBAL  HAS NOT BEEN USED,
C          SET LOW=1, IGH=N.
C
C        ORTR AND ORTI CONTAIN INFORMATION ABOUT THE UNITARY TRANS-
C          FORMATIONS USED IN THE REDUCTION BY  CORTH, IF PERFORMED.
C          ONLY ELEMENTS LOW THROUGH IGH ARE USED.  IF THE EIGENVECTORS
C          OF THE HESSENBERG MATRIX ARE DESIRED, SET ORTR(J) AND
C          ORTI(J) TO 0.0D0 FOR THESE ELEMENTS.
C
C        HR AND HI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE COMPLEX UPPER HESSENBERG MATRIX.
C          THEIR LOWER TRIANGLES BELOW THE SUBDIAGONAL CONTAIN FURTHER
C          INFORMATION ABOUT THE TRANSFORMATIONS WHICH WERE USED IN THE
C          REDUCTION BY  CORTH, IF PERFORMED.  IF THE EIGENVECTORS OF
C          THE HESSENBERG MATRIX ARE DESIRED, THESE ELEMENTS MAY BE
C          ARBITRARY.
C
C     ON OUTPUT
C
C        ORTR, ORTI, AND THE UPPER HESSENBERG PORTIONS OF HR AND HI
C          HAVE BEEN DESTROYED.
C
C        WR AND WI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE EIGENVALUES.  IF AN ERROR
C          EXIT IS MADE, THE EIGENVALUES SHOULD BE CORRECT
C          FOR INDICES IERR+1,...,N.
C
C        ZR AND ZI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE EIGENVECTORS.  THE EIGENVECTORS
C          ARE UNNORMALIZED.  IF AN ERROR EXIT IS MADE, NONE OF
C          THE EIGENVECTORS HAS BEEN FOUND.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE LIMIT OF 30*N ITERATIONS IS EXHAUSTED
C                     WHILE THE J-TH EIGENVALUE IS BEING SOUGHT.
C
C     CALLS CDIV FOR COMPLEX DIVISION.
C     CALLS CSROOT FOR COMPLEX SQUARE ROOT.
C     CALLS PYTHAG FOR  DSQRT(A*A + B*B) .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IERR = 0
C     .......... INITIALIZE EIGENVECTOR MATRIX ..........
      DO 101 J = 1, N
C
         DO 100 I = 1, N
            ZR(I,J) = 0.0D0
            ZI(I,J) = 0.0D0
  100    CONTINUE
         ZR(J,J) = 1.0D0
  101 CONTINUE
C     .......... FORM THE MATRIX OF ACCUMULATED TRANSFORMATIONS
C                FROM THE INFORMATION LEFT BY CORTH ..........
      IEND = IGH - LOW - 1
      IF (IEND) 180, 150, 105
C     .......... FOR I=IGH-1 STEP -1 UNTIL LOW+1 DO -- ..........
  105 DO 140 II = 1, IEND
         I = IGH - II
         IF (ORTR(I) .EQ. 0.0D0 .AND. ORTI(I) .EQ. 0.0D0) GO TO 140
         IF (HR(I,I-1) .EQ. 0.0D0 .AND. HI(I,I-1) .EQ. 0.0D0) GO TO 140
C     .......... NORM BELOW IS NEGATIVE OF H FORMED IN CORTH ..........
         NORM = HR(I,I-1) * ORTR(I) + HI(I,I-1) * ORTI(I)
         IP1 = I + 1
C
         DO 110 K = IP1, IGH
            ORTR(K) = HR(K,I-1)
            ORTI(K) = HI(K,I-1)
  110    CONTINUE
C
         DO 130 J = I, IGH
            SR = 0.0D0
            SI = 0.0D0
C
            DO 115 K = I, IGH
               SR = SR + ORTR(K) * ZR(K,J) + ORTI(K) * ZI(K,J)
               SI = SI + ORTR(K) * ZI(K,J) - ORTI(K) * ZR(K,J)
  115       CONTINUE
C
            SR = SR / NORM
            SI = SI / NORM
C
            DO 120 K = I, IGH
               ZR(K,J) = ZR(K,J) + SR * ORTR(K) - SI * ORTI(K)
               ZI(K,J) = ZI(K,J) + SR * ORTI(K) + SI * ORTR(K)
  120       CONTINUE
C
  130    CONTINUE
C
  140 CONTINUE
C     .......... CREATE REAL SUBDIAGONAL ELEMENTS ..........
  150 L = LOW + 1
C
      DO 170 I = L, IGH
         LL = MIN0(I+1,IGH)
         IF (HI(I,I-1) .EQ. 0.0D0) GO TO 170
         NORM = PYTHAG(HR(I,I-1),HI(I,I-1))
         YR = HR(I,I-1) / NORM
         YI = HI(I,I-1) / NORM
         HR(I,I-1) = NORM
         HI(I,I-1) = 0.0D0
C
         DO 155 J = I, N
            SI = YR * HI(I,J) - YI * HR(I,J)
            HR(I,J) = YR * HR(I,J) + YI * HI(I,J)
            HI(I,J) = SI
  155    CONTINUE
C
         DO 160 J = 1, LL
            SI = YR * HI(J,I) + YI * HR(J,I)
            HR(J,I) = YR * HR(J,I) - YI * HI(J,I)
            HI(J,I) = SI
  160    CONTINUE
C
         DO 165 J = LOW, IGH
            SI = YR * ZI(J,I) + YI * ZR(J,I)
            ZR(J,I) = YR * ZR(J,I) - YI * ZI(J,I)
            ZI(J,I) = SI
  165    CONTINUE
C
  170 CONTINUE
C     .......... STORE ROOTS ISOLATED BY CBAL ..........
  180 DO 200 I = 1, N
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 200
         WR(I) = HR(I,I)
         WI(I) = HI(I,I)
  200 CONTINUE
C
      EN = IGH
      TR = 0.0D0
      TI = 0.0D0
      ITN = 30*N
C     .......... SEARCH FOR NEXT EIGENVALUE ..........
  220 IF (EN .LT. LOW) GO TO 680
      ITS = 0
      ENM1 = EN - 1
C     .......... LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
C                FOR L=EN STEP -1 UNTIL LOW DO -- ..........
  240 DO 260 LL = LOW, EN
         L = EN + LOW - LL
         IF (L .EQ. LOW) GO TO 300
         TST1 = DABS(HR(L-1,L-1)) + DABS(HI(L-1,L-1))
     X            + DABS(HR(L,L)) + DABS(HI(L,L))
         TST2 = TST1 + DABS(HR(L,L-1))
         IF (TST2 .EQ. TST1) GO TO 300
  260 CONTINUE
C     .......... FORM SHIFT ..........
  300 IF (L .EQ. EN) GO TO 660
      IF (ITN .EQ. 0) GO TO 1000
      IF (ITS .EQ. 10 .OR. ITS .EQ. 20) GO TO 320
      SR = HR(EN,EN)
      SI = HI(EN,EN)
      XR = HR(ENM1,EN) * HR(EN,ENM1)
      XI = HI(ENM1,EN) * HR(EN,ENM1)
      IF (XR .EQ. 0.0D0 .AND. XI .EQ. 0.0D0) GO TO 340
      YR = (HR(ENM1,ENM1) - SR) / 2.0D0
      YI = (HI(ENM1,ENM1) - SI) / 2.0D0
      CALL CSROOT(YR**2-YI**2+XR,2.0D0*YR*YI+XI,ZZR,ZZI)
      IF (YR * ZZR + YI * ZZI .GE. 0.0D0) GO TO 310
      ZZR = -ZZR
      ZZI = -ZZI
  310 CALL CDIV(XR,XI,YR+ZZR,YI+ZZI,XR,XI)
      SR = SR - XR
      SI = SI - XI
      GO TO 340
C     .......... FORM EXCEPTIONAL SHIFT ..........
  320 SR = DABS(HR(EN,ENM1)) + DABS(HR(ENM1,EN-2))
      SI = 0.0D0
C
  340 DO 360 I = LOW, EN
         HR(I,I) = HR(I,I) - SR
         HI(I,I) = HI(I,I) - SI
  360 CONTINUE
C
      TR = TR + SR
      TI = TI + SI
      ITS = ITS + 1
      ITN = ITN - 1
C     .......... REDUCE TO TRIANGLE (ROWS) ..........
      LP1 = L + 1
C
      DO 500 I = LP1, EN
         SR = HR(I,I-1)
         HR(I,I-1) = 0.0D0
         NORM = PYTHAG(PYTHAG(HR(I-1,I-1),HI(I-1,I-1)),SR)
         XR = HR(I-1,I-1) / NORM
         WR(I-1) = XR
         XI = HI(I-1,I-1) / NORM
         WI(I-1) = XI
         HR(I-1,I-1) = NORM
         HI(I-1,I-1) = 0.0D0
         HI(I,I-1) = SR / NORM
C

         DO 490 J = I, N
            YR = HR(I-1,J)
            YI = HI(I-1,J)
            ZZR = HR(I,J)
            ZZI = HI(I,J)
            HR(I-1,J) = XR * YR + XI * YI + HI(I,I-1) * ZZR
            HI(I-1,J) = XR * YI - XI * YR + HI(I,I-1) * ZZI
            HR(I,J) = XR * ZZR - XI * ZZI - HI(I,I-1) * YR
            HI(I,J) = XR * ZZI + XI * ZZR - HI(I,I-1) * YI
  490    CONTINUE
C
  500 CONTINUE
C
      SI = HI(EN,EN)
      IF (SI .EQ. 0.0D0) GO TO 540
      NORM = PYTHAG(HR(EN,EN),SI)
      SR = HR(EN,EN) / NORM
      SI = SI / NORM
      HR(EN,EN) = NORM
      HI(EN,EN) = 0.0D0
      IF (EN .EQ. N) GO TO 540
      IP1 = EN + 1
C
      DO 520 J = IP1, N
         YR = HR(EN,J)
         YI = HI(EN,J)
         HR(EN,J) = SR * YR + SI * YI
         HI(EN,J) = SR * YI - SI * YR
  520 CONTINUE
C     .......... INVERSE OPERATION (COLUMNS) ..........
  540 DO 600 J = LP1, EN
         XR = WR(J-1)
         XI = WI(J-1)
C
         DO 580 I = 1, J
            YR = HR(I,J-1)
            YI = 0.0D0
            ZZR = HR(I,J)
            ZZI = HI(I,J)
            IF (I .EQ. J) GO TO 560
            YI = HI(I,J-1)
            HI(I,J-1) = XR * YI + XI * YR + HI(J,J-1) * ZZI
  560       HR(I,J-1) = XR * YR - XI * YI + HI(J,J-1) * ZZR
            HR(I,J) = XR * ZZR + XI * ZZI - HI(J,J-1) * YR
            HI(I,J) = XR * ZZI - XI * ZZR - HI(J,J-1) * YI
  580    CONTINUE
C
         DO 590 I = LOW, IGH
            YR = ZR(I,J-1)
            YI = ZI(I,J-1)
            ZZR = ZR(I,J)
            ZZI = ZI(I,J)
            ZR(I,J-1) = XR * YR - XI * YI + HI(J,J-1) * ZZR
            ZI(I,J-1) = XR * YI + XI * YR + HI(J,J-1) * ZZI
            ZR(I,J) = XR * ZZR + XI * ZZI - HI(J,J-1) * YR
            ZI(I,J) = XR * ZZI - XI * ZZR - HI(J,J-1) * YI
  590    CONTINUE
C
  600 CONTINUE
C
      IF (SI .EQ. 0.0D0) GO TO 240
C
      DO 630 I = 1, EN
         YR = HR(I,EN)
         YI = HI(I,EN)
         HR(I,EN) = SR * YR - SI * YI
         HI(I,EN) = SR * YI + SI * YR
  630 CONTINUE
C
      DO 640 I = LOW, IGH
         YR = ZR(I,EN)
         YI = ZI(I,EN)
         ZR(I,EN) = SR * YR - SI * YI
         ZI(I,EN) = SR * YI + SI * YR
  640 CONTINUE
C
      GO TO 240
C     .......... A ROOT FOUND ..........
  660 HR(EN,EN) = HR(EN,EN) + TR
      WR(EN) = HR(EN,EN)
      HI(EN,EN) = HI(EN,EN) + TI
      WI(EN) = HI(EN,EN)
      EN = ENM1
      GO TO 220
C     .......... ALL ROOTS FOUND.  BACKSUBSTITUTE TO FIND
C                VECTORS OF UPPER TRIANGULAR FORM ..........
  680 NORM = 0.0D0
C
      DO 720 I = 1, N
C
         DO 720 J = I, N
            TR = DABS(HR(I,J)) + DABS(HI(I,J))
            IF (TR .GT. NORM) NORM = TR
  720 CONTINUE
C
      IF (N .EQ. 1 .OR. NORM .EQ. 0.0D0) GO TO 1001
C     .......... FOR EN=N STEP -1 UNTIL 2 DO -- ..........
      DO 800 NN = 2, N
         EN = N + 2 - NN
         XR = WR(EN)
         XI = WI(EN)
         HR(EN,EN) = 1.0D0
         HI(EN,EN) = 0.0D0
         ENM1 = EN - 1
C     .......... FOR I=EN-1 STEP -1 UNTIL 1 DO -- ..........
         DO 780 II = 1, ENM1
            I = EN - II
            ZZR = 0.0D0
            ZZI = 0.0D0
            IP1 = I + 1
C
            DO 740 J = IP1, EN
               ZZR = ZZR + HR(I,J) * HR(J,EN) - HI(I,J) * HI(J,EN)
               ZZI = ZZI + HR(I,J) * HI(J,EN) + HI(I,J) * HR(J,EN)
  740       CONTINUE
C
            YR = XR - WR(I)
            YI = XI - WI(I)
            IF (YR .NE. 0.0D0 .OR. YI .NE. 0.0D0) GO TO 765
               TST1 = NORM
               YR = TST1
  760          YR = 0.01D0 * YR
               TST2 = NORM + YR
               IF (TST2 .GT. TST1) GO TO 760
  765       CONTINUE
            CALL CDIV(ZZR,ZZI,YR,YI,HR(I,EN),HI(I,EN))
C     .......... OVERFLOW CONTROL ..........
            TR = DABS(HR(I,EN)) + DABS(HI(I,EN))
            IF (TR .EQ. 0.0D0) GO TO 780
            TST1 = TR
            TST2 = TST1 + 1.0D0/TST1
            IF (TST2 .GT. TST1) GO TO 780
            DO 770 J = I, EN
               HR(J,EN) = HR(J,EN)/TR
               HI(J,EN) = HI(J,EN)/TR
  770       CONTINUE
C
  780    CONTINUE
C
  800 CONTINUE
C     .......... END BACKSUBSTITUTION ..........
      ENM1 = N - 1
C     .......... VECTORS OF ISOLATED ROOTS ..........
      DO  840 I = 1, ENM1
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 840
         IP1 = I + 1
C
         DO 820 J = IP1, N
            ZR(I,J) = HR(I,J)
            ZI(I,J) = HI(I,J)
  820    CONTINUE
C
  840 CONTINUE
C     .......... MULTIPLY BY TRANSFORMATION MATRIX TO GIVE
C                VECTORS OF ORIGINAL FULL MATRIX.
C                FOR J=N STEP -1 UNTIL LOW+1 DO -- ..........
      DO 880 JJ = LOW, ENM1
         J = N + LOW - JJ
         M = MIN0(J,IGH)
C
         DO 880 I = LOW, IGH
            ZZR = 0.0D0
            ZZI = 0.0D0
C
            DO 860 K = LOW, M
               ZZR = ZZR + ZR(I,K) * HR(K,J) - ZI(I,K) * HI(K,J)
               ZZI = ZZI + ZR(I,K) * HI(K,J) + ZI(I,K) * HR(K,J)
  860       CONTINUE
C
            ZR(I,J) = ZZR
            ZI(I,J) = ZZI
  880 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- ALL EIGENVALUES HAVE NOT
C                CONVERGED AFTER 30*N ITERATIONS ..........
 1000 IERR = EN
 1001 RETURN
      END

      SUBROUTINE COMQR(NM,N,LOW,IGH,HR,HI,WR,WI,IERR)
C
      INTEGER I,J,L,N,EN,LL,NM,IGH,ITN,ITS,LOW,LP1,ENM1,IERR
      real*8 HR(NM,N),HI(NM,N),WR(N),WI(N)
      real*8 SI,SR,TI,TR,XI,XR,YI,YR,ZZI,ZZR,NORM,TST1,TST2,
     X       PYTHAG
C
C     THIS SUBROUTINE IS A TRANSLATION OF A UNITARY ANALOGUE OF THE
C     ALGOL PROCEDURE  COMLR, NUM. MATH. 12, 369-376(1968) BY MARTIN
C     AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 396-403(1971).
C     THE UNITARY ANALOGUE SUBSTITUTES THE QR ALGORITHM OF FRANCIS
C     (COMP. JOUR. 4, 332-345(1962)) FOR THE LR ALGORITHM.
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A COMPLEX
C     UPPER HESSENBERG MATRIX BY THE QR METHOD.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
C          SUBROUTINE  CBAL.  IF  CBAL  HAS NOT BEEN USED,
C          SET LOW=1, IGH=N.
C
C        HR AND HI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE COMPLEX UPPER HESSENBERG MATRIX.
C          THEIR LOWER TRIANGLES BELOW THE SUBDIAGONAL CONTAIN
C          INFORMATION ABOUT THE UNITARY TRANSFORMATIONS USED IN
C          THE REDUCTION BY  CORTH, IF PERFORMED.
C
C     ON OUTPUT
C
C        THE UPPER HESSENBERG PORTIONS OF HR AND HI HAVE BEEN
C          DESTROYED.  THEREFORE, THEY MUST BE SAVED BEFORE
C          CALLING  COMQR  IF SUBSEQUENT CALCULATION OF
C          EIGENVECTORS IS TO BE PERFORMED.
C
C        WR AND WI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE EIGENVALUES.  IF AN ERROR
C          EXIT IS MADE, THE EIGENVALUES SHOULD BE CORRECT
C          FOR INDICES IERR+1,...,N.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE LIMIT OF 30*N ITERATIONS IS EXHAUSTED
C                     WHILE THE J-TH EIGENVALUE IS BEING SOUGHT.
C
C     CALLS CDIV FOR COMPLEX DIVISION.
C     CALLS CSROOT FOR COMPLEX SQUARE ROOT.
C     CALLS PYTHAG FOR  DSQRT(A*A + B*B) .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IERR = 0
      IF (LOW .EQ. IGH) GO TO 180
C     .......... CREATE REAL SUBDIAGONAL ELEMENTS ..........
      L = LOW + 1
C
      DO 170 I = L, IGH
         LL = MIN0(I+1,IGH)
         IF (HI(I,I-1) .EQ. 0.0D0) GO TO 170
         NORM = PYTHAG(HR(I,I-1),HI(I,I-1))
         YR = HR(I,I-1) / NORM
         YI = HI(I,I-1) / NORM
         HR(I,I-1) = NORM
         HI(I,I-1) = 0.0D0
C
         DO 155 J = I, IGH
            SI = YR * HI(I,J) - YI * HR(I,J)
            HR(I,J) = YR * HR(I,J) + YI * HI(I,J)
            HI(I,J) = SI
  155    CONTINUE
C
         DO 160 J = LOW, LL
            SI = YR * HI(J,I) + YI * HR(J,I)
            HR(J,I) = YR * HR(J,I) - YI * HI(J,I)
            HI(J,I) = SI
  160    CONTINUE
C
  170 CONTINUE
C     .......... STORE ROOTS ISOLATED BY CBAL ..........
  180 DO 200 I = 1, N
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 200
         WR(I) = HR(I,I)
         WI(I) = HI(I,I)
  200 CONTINUE
C
      EN = IGH
      TR = 0.0D0
      TI = 0.0D0
      ITN = 30*N
C     .......... SEARCH FOR NEXT EIGENVALUE ..........
  220 IF (EN .LT. LOW) GO TO 1001
      ITS = 0
      ENM1 = EN - 1
C     .......... LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
C                FOR L=EN STEP -1 UNTIL LOW D0 -- ..........
  240 DO 260 LL = LOW, EN
         L = EN + LOW - LL
         IF (L .EQ. LOW) GO TO 300
         TST1 = DABS(HR(L-1,L-1)) + DABS(HI(L-1,L-1))
     X            + DABS(HR(L,L)) + DABS(HI(L,L))
         TST2 = TST1 + DABS(HR(L,L-1))
         IF (TST2 .EQ. TST1) GO TO 300
  260 CONTINUE
C     .......... FORM SHIFT ..........
  300 IF (L .EQ. EN) GO TO 660
      IF (ITN .EQ. 0) GO TO 1000
      IF (ITS .EQ. 10 .OR. ITS .EQ. 20) GO TO 320
      SR = HR(EN,EN)
      SI = HI(EN,EN)
      XR = HR(ENM1,EN) * HR(EN,ENM1)
      XI = HI(ENM1,EN) * HR(EN,ENM1)
      IF (XR .EQ. 0.0D0 .AND. XI .EQ. 0.0D0) GO TO 340
      YR = (HR(ENM1,ENM1) - SR) / 2.0D0
      YI = (HI(ENM1,ENM1) - SI) / 2.0D0
      CALL CSROOT(YR**2-YI**2+XR,2.0D0*YR*YI+XI,ZZR,ZZI)
      IF (YR * ZZR + YI * ZZI .GE. 0.0D0) GO TO 310
      ZZR = -ZZR
      ZZI = -ZZI
  310 CALL CDIV(XR,XI,YR+ZZR,YI+ZZI,XR,XI)
      SR = SR - XR
      SI = SI - XI
      GO TO 340
C     .......... FORM EXCEPTIONAL SHIFT ..........
  320 SR = DABS(HR(EN,ENM1)) + DABS(HR(ENM1,EN-2))
      SI = 0.0D0
C
  340 DO 360 I = LOW, EN
         HR(I,I) = HR(I,I) - SR
         HI(I,I) = HI(I,I) - SI
  360 CONTINUE
C
      TR = TR + SR
      TI = TI + SI
      ITS = ITS + 1
      ITN = ITN - 1
C     .......... REDUCE TO TRIANGLE (ROWS) ..........
      LP1 = L + 1
C
      DO 500 I = LP1, EN
         SR = HR(I,I-1)
         HR(I,I-1) = 0.0D0
         NORM = PYTHAG(PYTHAG(HR(I-1,I-1),HI(I-1,I-1)),SR)
         XR = HR(I-1,I-1) / NORM
         WR(I-1) = XR
         XI = HI(I-1,I-1) / NORM
         WI(I-1) = XI
         HR(I-1,I-1) = NORM
         HI(I-1,I-1) = 0.0D0
         HI(I,I-1) = SR / NORM
C
         DO 490 J = I, EN
            YR = HR(I-1,J)
            YI = HI(I-1,J)
            ZZR = HR(I,J)
            ZZI = HI(I,J)
            HR(I-1,J) = XR * YR + XI * YI + HI(I,I-1) * ZZR
            HI(I-1,J) = XR * YI - XI * YR + HI(I,I-1) * ZZI
            HR(I,J) = XR * ZZR - XI * ZZI - HI(I,I-1) * YR
            HI(I,J) = XR * ZZI + XI * ZZR - HI(I,I-1) * YI
  490    CONTINUE
C
  500 CONTINUE
C
      SI = HI(EN,EN)
      IF (SI .EQ. 0.0D0) GO TO 540
      NORM = PYTHAG(HR(EN,EN),SI)
      SR = HR(EN,EN) / NORM
      SI = SI / NORM
      HR(EN,EN) = NORM
      HI(EN,EN) = 0.0D0
C     .......... INVERSE OPERATION (COLUMNS) ..........
  540 DO 600 J = LP1, EN
         XR = WR(J-1)
         XI = WI(J-1)
C
         DO 580 I = L, J
            YR = HR(I,J-1)
            YI = 0.0D0
            ZZR = HR(I,J)
            ZZI = HI(I,J)
            IF (I .EQ. J) GO TO 560
            YI = HI(I,J-1)
            HI(I,J-1) = XR * YI + XI * YR + HI(J,J-1) * ZZI
  560       HR(I,J-1) = XR * YR - XI * YI + HI(J,J-1) * ZZR
            HR(I,J) = XR * ZZR + XI * ZZI - HI(J,J-1) * YR
            HI(I,J) = XR * ZZI - XI * ZZR - HI(J,J-1) * YI
  580    CONTINUE
C
  600 CONTINUE
C
      IF (SI .EQ. 0.0D0) GO TO 240
C
      DO 630 I = L, EN
         YR = HR(I,EN)
         YI = HI(I,EN)
         HR(I,EN) = SR * YR - SI * YI
         HI(I,EN) = SR * YI + SI * YR
  630 CONTINUE
C
      GO TO 240
C     .......... A ROOT FOUND ..........
  660 WR(EN) = HR(EN,EN) + TR
      WI(EN) = HI(EN,EN) + TI
      EN = ENM1
      GO TO 220
C     .......... SET ERROR -- ALL EIGENVALUES HAVE NOT
C                CONVERGED AFTER 30*N ITERATIONS ..........
 1000 IERR = EN
 1001 RETURN
      END
      SUBROUTINE CORTH(NM,N,LOW,IGH,AR,AI,ORTR,ORTI)
C
      INTEGER I,J,M,N,II,JJ,LA,MP,NM,IGH,KP1,LOW
      real*8 AR(NM,N),AI(NM,N),ORTR(IGH),ORTI(IGH)
      real*8 F,G,H,FI,FR,SCALE,PYTHAG
C
C     THIS SUBROUTINE IS A TRANSLATION OF A COMPLEX ANALOGUE OF
C     THE ALGOL PROCEDURE ORTHES, NUM. MATH. 12, 349-368(1968)
C     BY MARTIN AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
C
C     GIVEN A COMPLEX GENERAL MATRIX, THIS SUBROUTINE
C     REDUCES A SUBMATRIX SITUATED IN ROWS AND COLUMNS
C     LOW THROUGH IGH TO UPPER HESSENBERG FORM BY
C     UNITARY SIMILARITY TRANSFORMATIONS.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
C          SUBROUTINE  CBAL.  IF  CBAL  HAS NOT BEEN USED,
C          SET LOW=1, IGH=N.
C
C        AR AND AI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE COMPLEX INPUT MATRIX.
C
C     ON OUTPUT
C
C        AR AND AI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE HESSENBERG MATRIX.  INFORMATION
C          ABOUT THE UNITARY TRANSFORMATIONS USED IN THE REDUCTION
C          IS STORED IN THE REMAINING TRIANGLES UNDER THE
C          HESSENBERG MATRIX.
C
C        ORTR AND ORTI CONTAIN FURTHER INFORMATION ABOUT THE
C          TRANSFORMATIONS.  ONLY ELEMENTS LOW THROUGH IGH ARE USED.
C
C     CALLS PYTHAG FOR  DSQRT(A*A + B*B) .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      LA = IGH - 1
      KP1 = LOW + 1
      IF (LA .LT. KP1) GO TO 200
C
      DO 180 M = KP1, LA
         H = 0.0D0
         ORTR(M) = 0.0D0
         ORTI(M) = 0.0D0
         SCALE = 0.0D0
C     .......... SCALE COLUMN (ALGOL TOL THEN NOT NEEDED) ..........
         DO 90 I = M, IGH
   90    SCALE = SCALE + DABS(AR(I,M-1)) + DABS(AI(I,M-1))
C
         IF (SCALE .EQ. 0.0D0) GO TO 180
         MP = M + IGH
C     .......... FOR I=IGH STEP -1 UNTIL M DO -- ..........
         DO 100 II = M, IGH
            I = MP - II
            ORTR(I) = AR(I,M-1) / SCALE
            ORTI(I) = AI(I,M-1) / SCALE
            H = H + ORTR(I) * ORTR(I) + ORTI(I) * ORTI(I)
  100    CONTINUE
C
         G = DSQRT(H)
         F = PYTHAG(ORTR(M),ORTI(M))
         IF (F .EQ. 0.0D0) GO TO 103
         H = H + F * G
         G = G / F
         ORTR(M) = (1.0D0 + G) * ORTR(M)
         ORTI(M) = (1.0D0 + G) * ORTI(M)
         GO TO 105
C
  103    ORTR(M) = G
         AR(M,M-1) = SCALE
C     .......... FORM (I-(U*UT)/H) * A ..........
  105    DO 130 J = M, N
            FR = 0.0D0
            FI = 0.0D0
C     .......... FOR I=IGH STEP -1 UNTIL M DO -- ..........
            DO 110 II = M, IGH
               I = MP - II
               FR = FR + ORTR(I) * AR(I,J) + ORTI(I) * AI(I,J)
               FI = FI + ORTR(I) * AI(I,J) - ORTI(I) * AR(I,J)
  110       CONTINUE
C
            FR = FR / H
            FI = FI / H
C
            DO 120 I = M, IGH
               AR(I,J) = AR(I,J) - FR * ORTR(I) + FI * ORTI(I)
               AI(I,J) = AI(I,J) - FR * ORTI(I) - FI * ORTR(I)
  120       CONTINUE
C
  130    CONTINUE
C     .......... FORM (I-(U*UT)/H)*A*(I-(U*UT)/H) ..........
         DO 160 I = 1, IGH
            FR = 0.0D0
            FI = 0.0D0
C     .......... FOR J=IGH STEP -1 UNTIL M DO -- ..........
            DO 140 JJ = M, IGH
               J = MP - JJ
               FR = FR + ORTR(J) * AR(I,J) - ORTI(J) * AI(I,J)
               FI = FI + ORTR(J) * AI(I,J) + ORTI(J) * AR(I,J)
  140       CONTINUE
C
            FR = FR / H
            FI = FI / H
C
            DO 150 J = M, IGH
               AR(I,J) = AR(I,J) - FR * ORTR(J) - FI * ORTI(J)
               AI(I,J) = AI(I,J) + FR * ORTI(J) - FI * ORTR(J)
  150       CONTINUE
C
  160    CONTINUE
C
         ORTR(M) = SCALE * ORTR(M)
         ORTI(M) = SCALE * ORTI(M)
         AR(M,M-1) = -G * AR(M,M-1)
         AI(M,M-1) = -G * AI(M,M-1)
  180 CONTINUE
C
  200 RETURN
      END
      SUBROUTINE CBAL(NM,N,AR,AI,LOW,IGH,SCALE)
C
      INTEGER I,J,K,L,M,N,JJ,NM,IGH,LOW,IEXC
      real*8 AR(NM,N),AI(NM,N),SCALE(N)
      real*8 C,F,G,R,S,B2,RADIX
      LOGICAL NOCONV
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE
C     CBALANCE, WHICH IS A COMPLEX VERSION OF BALANCE,
C     NUM. MATH. 13, 293-304(1969) BY PARLETT AND REINSCH.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 315-326(1971).
C
C     THIS SUBROUTINE BALANCES A COMPLEX MATRIX AND ISOLATES
C     EIGENVALUES WHENEVER POSSIBLE.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        AR AND AI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE COMPLEX MATRIX TO BE BALANCED.
C
C     ON OUTPUT
C
C        AR AND AI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE BALANCED MATRIX.
C
C        LOW AND IGH ARE TWO INTEGERS SUCH THAT AR(I,J) AND AI(I,J)
C          ARE EQUAL TO ZERO IF
C           (1) I IS GREATER THAN J AND
C           (2) J=1,...,LOW-1 OR I=IGH+1,...,N.
C
C        SCALE CONTAINS INFORMATION DETERMINING THE
C           PERMUTATIONS AND SCALING FACTORS USED.
C
C     SUPPOSE THAT THE PRINCIPAL SUBMATRIX IN ROWS LOW THROUGH IGH
C     HAS BEEN BALANCED, THAT P(J) DENOTES THE INDEX INTERCHANGED
C     WITH J DURING THE PERMUTATION STEP, AND THAT THE ELEMENTS
C     OF THE DIAGONAL MATRIX USED ARE DENOTED BY D(I,J).  THEN
C        SCALE(J) = P(J),    FOR J = 1,...,LOW-1
C                 = D(J,J)       J = LOW,...,IGH
C                 = P(J)         J = IGH+1,...,N.
C     THE ORDER IN WHICH THE INTERCHANGES ARE MADE IS N TO IGH+1,
C     THEN 1 TO LOW-1.
C
C     NOTE THAT 1 IS RETURNED FOR IGH IF IGH IS ZERO FORMALLY.
C
C     THE ALGOL PROCEDURE EXC CONTAINED IN CBALANCE APPEARS IN
C     CBAL  IN LINE.  (NOTE THAT THE ALGOL ROLES OF IDENTIFIERS
C     K,L HAVE BEEN REVERSED.)
C
C     ARITHMETIC IS REAL THROUGHOUT.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      RADIX = 16.0D0
C
      B2 = RADIX * RADIX
      K = 1
      L = N
      GO TO 100
C     .......... IN-LINE PROCEDURE FOR ROW AND
C                COLUMN EXCHANGE ..........
   20 SCALE(M) = J
      IF (J .EQ. M) GO TO 50
C
      DO 30 I = 1, L
         F = AR(I,J)
         AR(I,J) = AR(I,M)
         AR(I,M) = F
         F = AI(I,J)
         AI(I,J) = AI(I,M)
         AI(I,M) = F
   30 CONTINUE
C
      DO 40 I = K, N
         F = AR(J,I)
         AR(J,I) = AR(M,I)
         AR(M,I) = F
         F = AI(J,I)
         AI(J,I) = AI(M,I)
         AI(M,I) = F
   40 CONTINUE
C
   50 GO TO (80,130), IEXC
C     .......... SEARCH FOR ROWS ISOLATING AN EIGENVALUE
C                AND PUSH THEM DOWN ..........
   80 IF (L .EQ. 1) GO TO 280
      L = L - 1
C     .......... FOR J=L STEP -1 UNTIL 1 DO -- ..........
  100 DO 120 JJ = 1, L
         J = L + 1 - JJ
C
         DO 110 I = 1, L
            IF (I .EQ. J) GO TO 110
            IF (AR(J,I) .NE. 0.0D0 .OR. AI(J,I) .NE. 0.0D0) GO TO 120
  110    CONTINUE
C
         M = L
         IEXC = 1
         GO TO 20
  120 CONTINUE
C
      GO TO 140
C     .......... SEARCH FOR COLUMNS ISOLATING AN EIGENVALUE
C                AND PUSH THEM LEFT ..........
  130 K = K + 1
C
  140 DO 170 J = K, L
C
         DO 150 I = K, L
            IF (I .EQ. J) GO TO 150
            IF (AR(I,J) .NE. 0.0D0 .OR. AI(I,J) .NE. 0.0D0) GO TO 170
  150    CONTINUE
C
         M = K
         IEXC = 2
         GO TO 20
  170 CONTINUE
C     .......... NOW BALANCE THE SUBMATRIX IN ROWS K TO L ..........
      DO 180 I = K, L
  180 SCALE(I) = 1.0D0
C     .......... ITERATIVE LOOP FOR NORM REDUCTION ..........
  190 NOCONV = .FALSE.
C
      DO 270 I = K, L
         C = 0.0D0
         R = 0.0D0
C
         DO 200 J = K, L
            IF (J .EQ. I) GO TO 200
            C = C + DABS(AR(J,I)) + DABS(AI(J,I))
            R = R + DABS(AR(I,J)) + DABS(AI(I,J))
  200    CONTINUE
C     .......... GUARD AGAINST ZERO C OR R DUE TO UNDERFLOW ..........
         IF (C .EQ. 0.0D0 .OR. R .EQ. 0.0D0) GO TO 270
         G = R / RADIX
         F = 1.0D0
         S = C + R
  210    IF (C .GE. G) GO TO 220
         F = F * RADIX
         C = C * B2
         GO TO 210
  220    G = R * RADIX
  230    IF (C .LT. G) GO TO 240
         F = F / RADIX
         C = C / B2
         GO TO 230
C     .......... NOW BALANCE ..........
  240    IF ((C + R) / F .GE. 0.95D0 * S) GO TO 270
         G = 1.0D0 / F
         SCALE(I) = SCALE(I) * F
         NOCONV = .TRUE.
C
         DO 250 J = K, N
            AR(I,J) = AR(I,J) * G
            AI(I,J) = AI(I,J) * G
  250    CONTINUE
C
         DO 260 J = 1, L
            AR(J,I) = AR(J,I) * F
            AI(J,I) = AI(J,I) * F
  260    CONTINUE
C
  270 CONTINUE
C
      IF (NOCONV) GO TO 190
C
  280 LOW = K
      IGH = L
      RETURN
      END
      SUBROUTINE CDIV(AR,AI,BR,BI,CR,CI)
      real*8 AR,AI,BR,BI,CR,CI
C
C     COMPLEX DIVISION, (CR,CI) = (AR,AI)/(BR,BI)
C
      real*8 S,ARS,AIS,BRS,BIS
      S = DABS(BR) + DABS(BI)
      ARS = AR/S
      AIS = AI/S
      BRS = BR/S
      BIS = BI/S
      S = BRS**2 + BIS**2
      CR = (ARS*BRS + AIS*BIS)/S
      CI = (AIS*BRS - ARS*BIS)/S
      RETURN
      END
      SUBROUTINE CSROOT(XR,XI,YR,YI)
      real*8 XR,XI,YR,YI
C
C     (YR,YI) = COMPLEX DSQRT(XR,XI)
C     BRANCH CHOSEN SO THAT YR .GE. 0.0 AND SIGN(YI) .EQ. SIGN(XI)
C
      real*8 S,TR,TI,PYTHAG
      TR = XR
      TI = XI
      S = DSQRT(0.5D0*(PYTHAG(TR,TI) + DABS(TR)))
      IF (TR .GE. 0.0D0) YR = S
      IF (TI .LT. 0.0D0) S = -S
      IF (TR .LE. 0.0D0) YI = S
      IF (TR .LT. 0.0D0) YR = 0.5D0*(TI/YI)
      IF (TR .GT. 0.0D0) YI = 0.5D0*(TI/YR)
      RETURN
      END
      real*8 FUNCTION EPSLON (X)
      real*8 X
C
C     ESTIMATE UNIT ROUNDOFF IN QUANTITIES OF SIZE X.
C
      real*8 A,B,C,EPS
C
C     THIS PROGRAM SHOULD FUNCTION PROPERLY ON ALL SYSTEMS
C     SATISFYING THE FOLLOWING TWO ASSUMPTIONS,
C        1.  THE BASE USED IN REPRESENTING FLOATING POINT
C            NUMBERS IS NOT A POWER OF THREE.
C        2.  THE QUANTITY  A  IN STATEMENT 10 IS REPRESENTED TO
C            THE ACCURACY USED IN FLOATING POINT VARIABLES
C            THAT ARE STORED IN MEMORY.
C     THE STATEMENT NUMBER 10 AND THE GO TO 10 ARE INTENDED TO
C     FORCE OPTIMIZING COMPILERS TO GENERATE CODE SATISFYING
C     ASSUMPTION 2.
C     UNDER THESE ASSUMPTIONS, IT SHOULD BE TRUE THAT,
C            A  IS NOT EXACTLY EQUAL TO FOUR-THIRDS,
C            B  HAS A ZERO FOR ITS LAST BIT OR DIGIT,
C            C  IS NOT EXACTLY EQUAL TO ONE,
C            EPS  MEASURES THE SEPARATION OF 1.0 FROM
C                 THE NEXT LARGER FLOATING POINT NUMBER.
C     THE DEVELOPERS OF EISPACK WOULD APPRECIATE BEING INFORMED
C     ABOUT ANY SYSTEMS WHERE THESE ASSUMPTIONS DO NOT HOLD.
C
C     THIS VERSION DATED 4/6/83.
C
      A = 4.0D0/3.0D0
   10 B = A - 1.0D0
      C = B + B + B
      EPS = DABS(C-1.0D0)
      IF (EPS .EQ. 0.0D0) GO TO 10
      EPSLON = EPS*DABS(X)
      RETURN
      END
      real*8 FUNCTION PYTHAG(A,B)
      real*8 A,B
C
C     FINDS DSQRT(A**2+B**2) WITHOUT OVERFLOW OR DESTRUCTIVE UNDERFLOW
C
      real*8 P,R,S,T,U
      P = DMAX1(DABS(A),DABS(B))
      IF (P .EQ. 0.0D0) GO TO 20
      R = (DMIN1(DABS(A),DABS(B))/P)**2
   10 CONTINUE
         T = 4.0D0 + R
         IF (T .EQ. 4.0D0) GO TO 20
         S = R/T
         U = 1.0D0 + 2.0D0*S
         P = U*P
         R = (S/U)**2 * R
      GO TO 10
   20 PYTHAG = P
      RETURN
      END
