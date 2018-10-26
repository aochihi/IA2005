c SUBROUTINE ker31s_05Avril.f
c DATE: 05 April 2000
c Refer Aochi et al.(Pageoph, 2000) doi:10.1007/PL00001072
	REAL(8) FUNCTION ker31s(c1, c2, c3, t, w)
	IMPLICIT NONE
	REAL(8) c1, c2, c3, t, w, t1
	REAL(8) kernel1, kernel2
	REAL(8) comp1, comp2, comp3
	REAL(8) p311, q311, e311
	EXTERNAL p311, q311, e311

	t1 = t + 1.0d0
	      kernel1 = p311(c1+0.50d0, c2, c3, t1, w)
     &		      - p311(c1-0.50d0, c2, c3, t1, w)	
              kernel2 = p311(c1+0.50d0, c2, c3, t, w)     
     &                - p311(c1-0.50d0, c2, c3, t, w)   
	      comp1 = kernel1 - kernel2

              kernel1 = q311(c1+0.50d0, c2, c3, t1, w)
     &                - q311(c1-0.50d0, c2, c3, t1, w)
              kernel2 = q311(c1+0.50d0, c2, c3, t, w)
     &                - q311(c1-0.50d0, c2, c3, t, w)
              comp2 = kernel1 - kernel2

              kernel1 = e311(c1+0.50d0, c2+0.50d0, c3, t1, w)
     &                - e311(c1-0.50d0, c2+0.50d0, c3, t1, w)
     &                - e311(c1+0.50d0, c2-0.50d0, c3, t1, w)
     &                + e311(c1-0.50d0, c2-0.50d0, c3, t1, w)
              kernel2 = e311(c1+0.50d0, c2+0.50d0, c3, t, w)
     &                - e311(c1-0.50d0, c2+0.50d0, c3, t, w)
     &                - e311(c1+0.50d0, c2-0.50d0, c3, t, w)
     &                + e311(c1-0.50d0, c2-0.50d0, c3, t, w)
              comp3 = kernel1 - kernel2

	      ker31s = comp1 + comp2 - comp3 

	return
	END
c
        REAL(8) FUNCTION p311(chi1, chi2, chi3, dt, w)
	IMPLICIT NONE
        REAL(8) chi1, chi2, chi3, dt, w
	REAL(8) c2a, c2b, r311b
        REAL(8) c, rc2, rc, x21, x1, r2, r
	REAL(8) chi1d2, fb1, fb2, fb3, fb4
	REAL(8) i11, j1, k1, l2, m

        c = 1.0d0/(w*sqrt(3.0d0))
        rc  = c*dt
        rc2 = rc**2
	chi1d2 = chi1**2 + chi3**2
        x21 = rc2 - chi1d2

	if(x21.le.0.) then
	  p311  = 0.0d0
	else
	  c2a = chi2 + 0.50d0
	  c2b = chi2 - 0.50d0
	  x1 = sqrt(x21)
	  if(abs(c2a).ge.x1.and.abs(c2b).ge.x1) then
	    if((c2a*c2b).lt.0.) then
              i11 = chi1/chi1d2
              k1  = rc2/chi1d2
              fb1 = i11*x1*(1.0d0 - 4.0d0*k1)/3.0d0
              p311 = -2.0*fb1
	      r311b = 8.0d0*chi3**2*fb1/chi1d2
	    else	      
	      p311  = 0.0d0 
	      r311b = 0.0d0
	    endif
	  else if(abs(c2a).le.x1.and.abs(c2b).le.x1) then
            i11 = chi1/chi1d2
	    k1  = rc2/chi1d2

            r2 = chi1**2 + c2a**2 + chi3**2
            r  = sqrt(r2)
	    j1 = 1.0d0/r2 + 2.0d0/chi1d2
	    l2 = rc*c2a/r
	    m  = rc/r
	    fb1 = l2*(1.0d0 - 2.0d0*rc2*j1/3.0d0) 
            fb3 = l2*(j1*(1.0d0 - 4.0d0/3.0d0*k1) - (m/r)**2)

            r2 = chi1**2 + c2b**2 + chi3**2
            r  = sqrt(r2)
            j1 = 1.0d0/r2 + 2.0d0/chi1d2
            l2 = rc*c2b/r
	    m  = rc/r
            fb2 = l2*(1.0d0 - 2.0d0*rc2*j1/3.0d0)
            fb4 = l2*(j1*(1.0d0 - 4.0d0/3.0d0*k1) - (m/r)**2)

	    p311 = -(fb1 - fb2)*i11
            r311b = 2.0d0*i11*chi3**2*(fb3 - fb4)
	  else if(c2a.gt.x1.and.c2b.lt.x1) then
            i11 = chi1/chi1d2
            k1  = rc2/chi1d2

	    fb1 = x1*(1.0d0 - 4.0d0*k1)/3.0d0
            fb3 = 2.0d0/chi1d2*fb1
            r2 = chi1**2 + c2b**2 + chi3**2
            r  = sqrt(r2)
            j1 = 1.0d0/r2 + 2.0d0/chi1d2
            l2 = rc*c2b/r
	    m  = rc/r
            fb2 = l2*(1.0d0 - 2.0d0*rc2*j1/3.0d0) 
            fb4 = l2*(j1*(1.0d0 - 4.0d0/3.0d0*k1) - (m/r)**2)

            p311 = -(fb1 - fb2)*i11
            r311b = 2.0d0*i11*chi3**2*(fb3 - fb4)
	  else if(c2a.gt.(-x1).and.c2b.lt.(-x1)) then
            i11 = chi1/chi1d2   
            k1  = rc2/chi1d2    

	    fb1 = x1*(1.0d0 - 4.0d0*k1)/3.0d0
            fb3 = 2.0d0/chi1d2*fb1
            r2 = chi1**2 + c2a**2 + chi3**2
            r  = sqrt(r2)       
            j1 = 1.0d0/r2 + 2.0d0/chi1d2        
            l2 = rc*c2a/r       
	    m  = rc/r
            fb2 = l2*(1.0d0 - 2.0d0*rc2*j1/3.0d0)       
            fb4 = l2*(j1*(1.0d0 - 4.0d0/3.0d0*k1) - (m/r)**2)

            p311  = -(fb1 + fb2)*i11
            r311b = 2.0d0*i11*chi3**2*(fb3 + fb4)
	  else
            pause 'wrong condition'
          endif
	  if(chi3.ne.0.) p311 = p311 + r311b
        endif

	return
	END
c
        REAL(8) FUNCTION q311(chi1, chi2, chi3, dt, w)
	IMPLICIT NONE
        REAL(8) chi1, chi2, chi3, dt, w
	REAL(8) c2a, c2b, r311a
        REAL(8) c, p, p3, rc2, rc, x21, x1, r2, r
	REAL(8) chi1d2, fb1, fb2, fb3, fb4, const
	REAL(8) i11, j1, k1, l2, m

        c = 1.0d0/w
        p = 1.0d0/sqrt(3.0d0)
	p3 = p**3
        rc  = c*dt
        rc2 = rc**2
	chi1d2 = chi1**2 + chi3**2
        x21 = rc2 - chi1d2
	if(x21.le.0.) then
	  q311 = 0.0d0
	  r311a = 0.0d0
	else
	  c2a = chi2 + 0.50d0
	  c2b = chi2 - 0.50d0
	  x1 = sqrt(x21)
	  if(abs(c2a).ge.x1.and.abs(c2b).ge.x1) then
	    if((c2a*c2b).lt.0.) then
              i11 = chi1/chi1d2
              k1  = rc2/chi1d2
              fb1 = 8.0d0/3.0d0*x1*p3*i11
              q311 = fb1*(1.0d0 - k1)
              r311a = chi3**2*fb1*(1.0d0 - 4.0d0*k1)/chi1d2
	    else
  	      q311 = 0.0d0
	      r311a = 0.0d0
	    endif 
	  else if(abs(c2a).le.x1.and.abs(c2b).le.x1) then
            i11 = chi1/chi1d2
	    k1  = rc2/chi1d2

            r2 = chi1**2 + c2a**2 + chi3**2
            r  = sqrt(r2)
	    j1 = 1.0d0/r2 + 2.0d0/chi1d2
	    l2 = rc*c2a/r
	    m  = rc/r
            fb1 = l2*(1.0d0 - rc2*j1/3.0d0)
            fb3 = l2*(j1*(1.0d0 - 4.0d0/3.0d0*k1) - (m/r)**2)
            r2 = chi1**2 + c2b**2 + chi3**2
            r  = sqrt(r2)
            j1  = 1.0d0/r2 + 2.0d0/chi1d2
            l2  = rc*c2b/r
	    m   = rc/r
            fb2 = l2*(1.0d0 - rc2*j1/3.0d0)
            fb4 = l2*(j1*(1.0d0 - 4.0d0/3.0d0*k1) - (m/r)**2)

	    const = 2.0d0*p3*i11
            q311 = const*(fb1 - fb2)
            r311a = const*chi3**2*(fb3 - fb4)
	  else if(c2a.gt.x1.and.c2b.lt.x1) then
            i11 = chi1/chi1d2
            k1  = rc2/chi1d2

	    const = 2.0d0/3.0d0*x1
            fb1 = const*(1.0d0 - k1)
            fb3 = const*(1.0d0 - 4.0d0*k1)/chi1d2

            r2 = chi1**2 + c2b**2 + chi3**2
            r  = sqrt(r2)
            j1  = 1.0d0/r2 + 2.0d0/chi1d2
            l2  = rc*c2b/r
	    m   = rc/r
            fb2 = l2*(1.0d0 - rc2*j1/3.0d0)
            fb4 = l2*(j1*(1.0d0 - 4.0d0/3.0d0*k1) - (m/r)**2)

	    const = 2.0d0*p3*i11
            q311 = const*(fb1 - fb2)
            r311a = const*chi3**2*(fb3 - fb4)
	  else if(c2a.gt.(-x1).and.c2b.lt.(-x1)) then
	    i11 = chi1/chi1d2
	    k1  = rc2/chi1d2

	    const = 2.0d0/3.0d0*x1
            fb1 = const*(1.0d0 - k1)
            fb3 = const*(1.0d0 - 4.0d0*k1)/chi1d2
            r2 = chi1**2 + c2a**2 + chi3**2
            r  = sqrt(r2)
            j1  = 1.0d0/r2 + 2.0d0/chi1d2
            l2  = rc*c2a/r
	    m  = rc/r
            fb2 = l2*(1.0d0 - rc2*j1/3.0d0)
            fb4 = l2*(j1*(1.0d0 - 4.0d0/3.0d0*k1) - (m/r)**2)

	    const = 2.0d0*p3*i11
	    q311 = const*(fb1 + fb2)
            r311a = const*chi3**2*(fb3 + fb4)
          else
            pause 'bad if'
          endif
	  if(chi3.ne.0.) q311 = q311 - r311a
        endif

        return
        END
c
        REAL(8) FUNCTION e311(chi1, chi2, chi3, dt, w)
	IMPLICIT NONE
        REAL(8) chi1, chi2, chi3, dt, w
        REAL(8) c, rc2, rc, x22, x21, x2, x1, r2, r
	REAL(8) pi, f1, f2
	REAL(8) chi1d2, chi1d, chi2d2, chi2d
	REAL(8) i11, j1, k1, l2, i22, j2, k2, l1

        c = 1.0d0/(w*sqrt(3.0d0))
	pi = acos(-1.0d0)
        rc  = c*dt
        rc2 = rc**2

	chi1d2 = chi1**2 + chi3**2
	chi1d  = sqrt(chi1d2)
	chi2d2 = chi2**2 + chi3**2
	chi2d  = sqrt(chi2d2)
        x21 = rc2 - chi1d2
        x22 = rc2 - chi2d2
        r2 = chi1**2 + chi2**2 + chi3**2
        r = sqrt(r2)

	i11 = chi1/chi1d2
	i22 = chi2/chi2d2
	j1  = 1.0d0/r2 + 2.0d0/chi1d2
	j2  = 1.0d0/r2 + 2.0d0/chi2d2
	k1  = rc2/chi1d2
	k2  = rc2/chi2d2
	l2  = rc*chi2/r
	l1  = rc*chi1/r

	if((chi1.ge.0.).and.(chi2.ge.0.).and.
     &		(rc2.gt.(chi3**2))) then
	  e311 = 2.0d0*pi
	else
	  e311 = 0.0d0
	endif

	if(x21.gt.0.) then
          x1 = sqrt(x21)
          if(abs(chi2).le.x1) then
	    f1 = i11*(x1 + l2)
	    f2 = atan(x1/chi1) + atan(chi2/chi1)
	    e311 = e311 + f1 - f2
          else if(chi2.gt.x1) then
            e311 = e311 + 2.0d0*(i11*x1 - atan(x1/chi1))
          endif
        endif

        if(x22.ge.0.) then
          x2 = sqrt(x22)
          if(abs(chi1).le.x2) then
	    f1 = i22*(x2 + l1)
	    f2 = atan(x2/chi2) + atan(chi1/chi2)
	    e311 = e311 + f1 - f2
          else if(chi1.gt.x2) then
	    e311 = e311 + 2.0d0*(i22*x2 - atan(x2/chi2))
          endif
        endif

        return
        END
c
