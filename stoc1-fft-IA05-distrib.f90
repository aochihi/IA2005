! program stoc1-fft-IA05-distrib.f90 (26th October 2018).
! under GNU Lesse General Public License v3.0. on @GitHub.com
!
! Refer Ide and Aochi (JGR, 2005). doi:1O.1029/2004JB003591.
! See README file or related file.
! >ifort stoc1-fft-IA05-distrib.f90 ran1.f fourn-d.f kernel31s_05Avril.f (verified on Unix and Windows)
IMPLICIT NONE
! 
INTEGER nmax, ndata1, ndata2
INTEGER npower, nscale, ixmax, itmx
PARAMETER ( nmax = 64, nscale = 4, npower = 3, ixmax = nmax*nscale**npower )
PARAMETER ( itmx = 500 )
PARAMETER ( ndata1 = 4*nmax, ndata2 = 4*nmax )
REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: vel, vel2
REAL(8), DIMENSION(:, :), ALLOCATABLE :: tau0, tp, tr, &
	stress, sigma, w, a, tau, dc
REAL(8), DIMENSION(:), ALLOCATABLE :: x0, y0, smrate, smoment
REAL(8) :: pi, mu, const, facbiem, facfft, &
		tp0, tr0, dc0, t0, dsreal, dtreal, coef, &
		p000, ker31s, dtau, dsigma, alpha, &
		ds, dt, rad, r0, rini, xhypo, yhypo, &
		xo, yo, r0dum, dcdum, dcmax, dim, smo, mw, &
		t, piece1, ans
REAL, DIMENSION(:, :), ALLOCATABLE :: dcorg
REAL :: ran1
INTEGER, DIMENSION(:, :), ALLOCATABLE :: iv, irup
INTEGER :: i, j, k, l, m, n, ndir, idata, ix, iy, &
		kmin, iter, kmax, itmx1, icheck, icheck2, it, &
		ndense, nasp, idum, iscale, ihypo, isim, isim0, nhypo, &
		nmax2, ns, i0, j0, i1, j1, k1, nscale2, npower2
INTEGER :: ndata(2)
DOUBLE COMPLEX zdata(ndata1*ndata2), zresp(ndata1*ndata2)
DOUBLE COMPLEX zvel(ndata1*ndata2, itmx)
DOUBLE COMPLEX zker(ndata1*ndata2, itmx)
DOUBLE COMPLEX zans(ndata1*ndata2)
EXTERNAL ker31s, ran1
CHARACTER*40 name2, name3, name4, name5, name6, name7, dir, param_file
CHARACTER*5  num, num2

! PARAMETER FILE
        param_file = "IA05.prm"
        open(11, file=param_file, status="old", err=99)
        read(11,*) isim0
        close(11)
        
!
! FIXED PARAMETER
!
	pi = acos(-1.0d0)
        facbiem = 2.0d0
	ds = 1.0d0
	dt = ds/facbiem
	kmin = itmx/nscale
! for FFT
	facfft = 1.0d0/(ndata1*ndata2)
	ndata(1) = ndata1
	ndata(2) = ndata2
!
! OUTPUT DIRECTORY
!
	dir = '.'
        ndir = index(dir, ' ')-1
!
! SCALE-INDEPENDENT PARAMETER
!
	mu = 32.40d0
	alpha = 6.0d0
	tp0 = 5.0d0
	tr0 = 0.0d0
	t0 = 3.0d0
	const = sqrt(3.0d0)/(4.0d0*pi)*mu
!
! SCALING PARAMETER
! unit (mu) = mu [GPa]/tb[MPa] dc0[m]/ds [km] = 1
! dc0 should be normalized by 0.001 ds.
!
! INITIAL SETTING
	dc0 = 0.250d0*ds
	r0 = 5.6250d0*ds
	rini = 3.75d0*ds
	ndense = 4
	dim = -2.0
	xhypo = (1+nmax)/2.
	yhypo = (1+nmax)/2.

ALLOCATE( vel(nmax, nmax, 0:itmx), vel2(nmax, nmax, 0:itmx) )
ALLOCATE(tau0(nmax, nmax),     tp(nmax, nmax),   dc(nmax, nmax), &
       stress(nmax, nmax),     tr(nmax, nmax),    a(nmax, nmax), &
     	sigma(nmax, nmax),      w(nmax, nmax), &
     	   iv(nmax, nmax),   irup(nmax, nmax),  tau(nmax, nmax) )
ALLOCATE( smrate(0:itmx), smoment(0:itmx) )
ALLOCATE( dcorg(ixmax, ixmax) )

!
! CREATING DC
!
! to be fixed (-411) for obtaining the same result as case1-1
        idum = -411
        dcmax = dc0*nscale**(npower + 1)
	dcorg = real(dcmax)

        nscale2 = nscale/2
        npower2 = npower*2
	nhypo = ndense*(nscale2*nscale2)**npower2
ALLOCATE( x0(nhypo), y0(nhypo) )

        do iscale = 0,  npower2
          nasp = ndense*(nscale2*nscale2)**(npower2 - iscale)
          r0dum = r0*nscale2**iscale
          dcdum = dc0*nscale2**iscale
          do ihypo = 1, nasp
            xo = ran1(idum)*ixmax
            yo = ran1(idum)*ixmax
            if(iscale.eq.0 ) then
                x0(ihypo) = xo
                y0(ihypo) = yo
            endif

            do i = int(xo - r0dum)-1, int(xo + r0dum)+1
              i1 = i
              if(i1.lt.1) i1 = ixmax + i
              if(i1.gt.ixmax) i1 =  i - ixmax
              do j = int(yo - r0dum)-1, int(yo + r0dum)+1
                j1 = j
                if(j1.lt.1) j1 = ixmax + j
                if(j1.gt.ixmax) j1 =  j - ixmax

                rad = sqrt((i-xo)**2 + (j-yo)**2)
                if( rad.le.r0dum ) then
                  if( dcorg(i1,j1).gt.dcdum ) dcorg(i1,j1) = real(dcdum)
                endif
	      enddo
	    enddo

          enddo
        enddo
	name7 = dir(1:ndir)//'/hetero.org'
	ns = ixmax/256
	if(ns.lt.1) ns = 1
	open(71, file=name7)
        do i=1, ixmax, ns
          do j=1, ixmax, ns
            write(71,'(2i7, f10.3)') i, j, dcorg(i,j)
          enddo
        enddo
        close(71)

!
! RESPONSE FUNCTION
!
p000 = ker31s(0.0d0, 0.0d0, 0.0d0, 0.0d0, facbiem)
do k = 1, itmx
  do idata=1, ndata1*ndata2
    zresp(idata) = cmplx(0.0d0, 0.0d0)
  enddo
  do i = 1-nmax, nmax-1
    ix = i + 1
    if(i.lt.0)  ix = ix + ndata1
    do j = 1-nmax, nmax-1
      piece1 = ker31s(dble(i), dble(j), 0.d0, dble(k), facbiem)
      iy = j + 1
      if(j.lt.0) iy = iy + ndata2
      idata = ix + (iy-1)*ndata1
      zresp(idata) = cmplx(piece1, 0.0d0)
    enddo
  enddo
  call fourn(zresp, ndata, 2, 1)
  do idata=1, ndata1*ndata2
    zker(idata, k) = zresp(idata)
  enddo
enddo

name2 = dir(1:ndir)//'/hoge2.dat'

open(12, file=name2)
write(12, '(5i10)') nmax, nscale, npower, ixmax, itmx
write(12, '(5f10.3)') mu, alpha, tp0, tr0, t0
write(12, '(4f10.3)') dc0, dcmax, r0, rini
write(12, '(4i10)') ndense, nscale2, npower2, nhypo
close(12)

!!
!! ITERATION OF HYPOCENTER LOCATION
!!
!! isim0 red by a parameter file
do isim = isim0, isim0
  ihypo = isim
!! test 
  if(isim.eq.0) ihypo = 1
!! two scenarios of Mw3.8 for Aochi & Burnol (2018)
  if(isim.eq.1) ihypo = 537
  if(isim.eq.2) ihypo = 546
! for checking large events
  if(isim.eq.1) ihypo = 806
  if(isim.eq.2) ihypo = 7987 
  if(isim.eq.3) ihypo = 9141
  if(isim.eq.4) ihypo = 10426
  if(isim.eq.5) ihypo = 12746 
  if(isim.eq.6) ihypo = 13375

  write(num, '(i5.5)') ihypo 
  name6 = dir(1:ndir)//'/output'//num(1:5)//'i.dat'
  open(16, file=name6)
  write(16, '(i10, 2f10.3)') ihypo, x0(ihypo), y0(ihypo)
  write(16, '(f10.3)') rini
  close(16)

!!	
!! ITERATION OF STAGE
!!
  STAGE: do iter = 0, npower
    kmax = itmx
    ns = nscale**iter

! RENORMALIZATION
    nmax2 = nmax*nscale**iter
    do i = 1, nmax
      i0 = ixmax/2 - nmax2/2 + ns*(i-1) + 1
      do j = 1, nmax
        j0 = ixmax/2 - nmax2/2 + ns*(j-1) + 1
        dc(i,j) = 0.0d0
        do i1 = 1, ns
          do j1 = 1, ns
            ix = i0 + i1 - 1 + int(x0(ihypo) - (ixmax+1)/2.)
            if(ix.lt.1) ix = ixmax + ix
            if(ix.gt.ixmax) ix =  ix - ixmax
            iy = j0 + j1 - 1 + int(y0(ihypo) - (ixmax+1)/2.)
            if(iy.lt.1) iy = ixmax + iy
            if(iy.gt.ixmax) iy =  iy - ixmax

            dc(i,j) = dc(i,j) + dble(dcorg(ix, iy))
          enddo
        enddo
        dc(i,j) = dc(i,j)/(ns**2)
      enddo
    enddo

! INITIALIZATION OF PARAMETERS
    vel = 0.0d0
    smrate = 0.0d0

    do i = 1, nmax
      do j = 1, nmax
        w(i,j) = 0.0d0
        tau0(i,j) = t0
        tp(i,j) = tp0
        tr(i,j) = tr0
        dc(i,j) = dc(i,j)/ns
        sigma(i,j) = tp(i,j)
        a(i,j) = (tp(i,j) - tr(i,j))/dc(i, j)
        iv(i,j) = 0
	irup(i,j) = -1
        
        if (iter.eq.0) then
          rad = sqrt((i-xhypo)**2 + (j-yhypo)**2)*ds
          if ( rad.le.rini ) then
            tp(i,j) = 0.0d0
            dc(i,j) = 0.0d0
            a(i,j) = 0.0d0
            sigma(i,j) = tp(i,j)
            iv(i,j) = 2
	  endif
	endif
      enddo
    enddo

    if( iter.ne.0 ) then
      kmax = itmx
      if( iter.eq.npower) kmax = itmx
      kmin = itmx1/nscale

!    do k = 1, itmx
      do k = 1, nscale*(itmx1/nscale)
        n = (k-1)/nscale + 1
        do i = 1, nmax
          l = (i-1)/nscale + 1 + (nscale-1)*nmax/(2*nscale) 
          do j = 1, nmax
            m = (j-1)/nscale + 1 + (nscale-1)*nmax/(2*nscale)
            if( vel2(i,j,k).ne.0. )then
              vel(l,m,n) = vel(l,m,n) + &
		vel2(i,j,k)/(nscale*nscale*nscale)
            endif
          enddo
        enddo
      enddo
    endif

    icheck = 0
    itmx1 = kmax
    smoment = 0.0d0

! ITERATION OF TIME
    TIME: do k = 1, kmax
      if( mod(k, 20).eq.0 ) write(6,'(a20, 3i5)') "SIMULATION", ihypo, iter, k
      icheck2 = 0

      if ( k.ne.1 ) then
        do idata=1, ndata1*ndata2
          zdata(idata) = cmplx(0.0d0, 0.0d0)
        enddo
        do j = 1, nmax
          do i = 1, nmax
            idata = i + (j-1)*ndata1
            zdata(idata) = cmplx(dble(vel(i,j,k-1)), 0.0d0)
          enddo
        enddo
        call fourn(zdata, ndata, 2, 1)
        do idata=1, ndata1*ndata2
          zvel(idata, k-1) = zdata(idata)
        enddo

        do idata=1, ndata1*ndata2
          zans(idata) = (0.0d0, 0.0d0)
          do n=1, k-1
            zans(idata) = zans(idata) +  &
		zker(idata, k-n)*zvel(idata, n)
          enddo
        enddo
        call fourn(zans, ndata, 2, -1)
      endif

      do i=1, nmax
        do j=1, nmax
          idata = i + (j-1)*ndata1

          if( k.ne.1 ) then
            ans = facfft*dble(zans(idata)) 
            dtau = tau0(i,j) + const*ans
          else
  	    dtau = tau0(i,j)
          endif
          dsigma = sigma(i,j)

	  if( iter.ne.0.and.k.le.kmin ) then
	    if( vel(i,j,k).ne.0. ) then
	      iv(i,j) = 1
	    else
	      iv(i,j) = 0
	      if( dtau.gt.dsigma ) iv(i,j) = 2
	    endif
	  else
            if( iv(i,j).eq.0 ) then
              vel(i,j,k) = 0.0d0
              if( dtau.gt.dsigma ) iv(i,j) = 2
            else
              iv(i,j) = 1
	      if( (const*p000 + a(i,j)*dt).ge.0. ) then
	        vel(i,j,k) = (tr(i,j) - dtau)/(const*p000)
	      else
                vel(i,j,k)  = (dsigma - dtau)/(const*p000 + a(i,j)*dt)
                if((w(i,j)+vel(i,j,k)*dt).gt.dc(i,j)) then
                  vel(i,j,k) = ( tr(i,j) - dtau)/(const*p000)
                endif
	      endif
              if(vel(i,j,k).lt.0.) then
                vel(i,j,k)  = 0.0d0
                iv(i, j) = 0
                if( dtau.gt.dsigma ) iv(i,j) = 2
              endif
            endif
	  endif

          stress(i,j) = dtau + const*p000*vel(i,j,k)
          tau(i,j) = dtau
          w(i,j) = w(i,j) + vel(i,j,k)*dt
	  if( iter.ne.npower ) vel2(i,j,k) = vel(i,j,k)
          if( w(i,j).gt.dc(i,j) ) then
            sigma(i,j) = tr(i,j)
          else
            sigma(i,j) = tp(i,j) - a(i,j)*w(i,j)
          endif
	  if( irup(i,j).lt.0.and.vel(i,j,k).ne.0.) irup(i,j) = k

          if((i.eq.1.or.i.eq.nmax).and.vel(i,j,k).ne.0.) icheck = 1
          if((j.eq.1.or.j.eq.nmax).and.vel(i,j,k).ne.0.) icheck = 1
	  if( vel(i,j,k).ne.0.) icheck2 = 1

	  smrate(k) = smrate(k) + vel(i,j,k)*alpha
          smoment(k) = smoment(k) + w(i,j)*ns
        enddo
      enddo

      if( iter.le.npower ) then
	write(num, '(i5.5)') ihypo
        write(num2,'(i4.4)') iter*1000 + k
        name3 = dir(1:ndir)//'/'//num(1:5)//'_step'//num2(1:4)//'.dat'
        open(13, file=name3)
	   
        do i = 1, nmax
          do j = 1, nmax
            write(13, '(f7.1, 1x, f7.1, 1x, 3f13.8)') &
     		real(i), real(j), vel(i,j,k)*alpha, w(i,j)*ns, stress(i,j)
          enddo
        enddo
        close(13)
      endif

      if( k.eq.itmx.or.(iter.ne.npower.and.icheck.eq.1).or.icheck2.eq.0) then
        write(num,'(i5.5)') ihypo
        write(num2,'(i1.1)') iter
        name3 = dir(1:ndir)//'/output'//num(1:5)//num2(1:1)//'.dat'
        open(13, file=name3)

	smo = 0.0d0
        do i = 1, nmax
          do j = 1, nmax
            write(13, 105) i, j, dc(i,j)*ns, irup(i,j), w(i,j)*ns
	    if(w(i,j).ne.0.) smo = smo + w(i,j)*ns
 105	  format(2i5, 1x, f9.3, i10, 1x, f9.3)
          enddo
        enddo
        close(13)

        coef = (0.4)**3*(ds*ns)**2*mu*10.0**9
	dsreal = 4.d0*ns*ds
	dtreal = dsreal/(2.*alpha*1000.)

	smo = 4.*(smo/1000.)*dsreal**2*mu*10.0**9
	mw = (log10(smo)-9.1)/1.5
        name4 = dir(1:ndir)//'/output'//num(1:5)//'f.dat'
        open(14, file=name4)
	write(14, '(3i5)') ihypo, iter, k
	write(14, '(e10.4, 2f10.3, f10.4)') smo, mw, dsreal, dtreal
        do it=1, k
          write(14, '(i5,f12.4,2e12.4)') it, it*dtreal, smoment(it)*coef, &
                  (smoment(it)-smoment(it-1))*coef/dtreal
        enddo
	close(14)

	name5 = dir(1:ndir)//'/moment'//num(1:5)//num2(1:1)//'.dat'
	open(15, file=name5)
	smo = 0.0d0
	do k1 = 1, k
	  smo = smo + smrate(k1)*dtreal*dsreal**2*mu*10.0**9
	  mw =  (log10(smo)-9.1)/1.5
	  write(15, '(f15.6, e15.5, e15.5, f15.5)') &
		k1*dtreal, smrate(k1)*dsreal**2*mu*10.0**9, smo, mw
	enddo
	close(15)

      endif

      if(iter.ne.npower.and.icheck.eq.1) then
        itmx1 = k - 1
        write(6,*) "returning at ...", itmx1
        exit TIME
      endif

      if(icheck2.eq.0) then
	write(6,*) "rupture terminated", iter, k
	exit STAGE
      endif

    enddo TIME
  enddo STAGE
enddo

 99     continue
        write(*,*) "END OF SIMULATION"

END

