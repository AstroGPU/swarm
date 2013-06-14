from math import *
from collections import namedtuple

M_PI = pi

def improve_mean_to_eccentric_annomaly_guess(e, M, x):
    """
    Adapted from C code
    {
    	double sx, cx
    	sincos(x,&sx,&cx)
    	double es = e*sx
    	double ec = e*cx
    	double f = x-es-M
    	double fp  = 1.-ec
    	double fpp = es
    	double fppp = ec
    	double dx = -f/fp
    	dx = -f/(fp+dx*fpp/2.)
    	dx = -f/(fp+dx*fpp/2.+dx*dx*fppp/6.)
    	return x+dx
    }
    """
    sx, cx = sincos(x)
    es = e * sx
    ec = e * cx
    f = x - es - M
    fp = 1.0 - ec
    fpp = es
    fppp = ec
    dx = -f / fp
    dx = -f / (fp+dx*fpp/2.0)
    dx = -f / (fp+dx*fpp/2.0+dx*dx*fppp/6.0)
    return x + dx



def mean_to_eccentric_annomaly(e, M):
    """
    {
    	const int ORBEL_EHIE_NMAX = 3

    	int i
    	int nper = (int)(M/(2.*M_PI))
    	M = M - nper*2.*M_PI
    	if(M<0.)  M = M + 2.*M_PI
    	assert(M>=0.)
    	if(M>M_PI)
    	{
    		M = 2.*M_PI - M
    		double x = pow(6.*M,1./3.)
    		for(i=1i<=ORBEL_EHIE_NMAX++i)
    			x = improve_mean_to_eccentric_annomaly_guess(e,M,x)
    		x = 2.*M_PI-x
    		return x
    	}
    	else
    	{
    		double x = pow(6.*M,1./3.)
    		for(i=1i<=ORBEL_EHIE_NMAX++i)
    			x = improve_mean_to_eccentric_annomaly_guess(e,M,x)
    		return x
    	}
    }"""
    ORBEL_EHIE_NMAX = 3
    nper = (M/(2.*M_PI))
    M = M - nper*2.*M_PI
    if (M<0.0):
        M = M + 2.*M_PI
    assert(M>=0.0)
    x = (6.0 * M) ** (1.0/3.0)
    if (M > M_PI):
        M = 2.0*M_PI - M
        for i in range(1,ORBEL_EHIE_NMAX):
            x = improve_mean_to_eccentric_annomaly_guess(e,M,x)
        x = 2.0*M_PI-x
        return x
    else:
        for i in range(1,ORBEL_EHIE_NMAX):
            x = improve_mean_to_eccentric_annomaly_guess(e,M,x)
        return x

KeplerianCoordinates = namedtuple('KeplerianCoordinates',
    ['a', 'e', 'i', 'O', 'w', 'M' ])



def calc_cartesian_for_ellipse(kepler_coordinates, GM):
    """{
    	double cape = mean_to_eccentric_annomaly(e,M)
    	double scap, ccap
    	sincos(cape,&scap,&ccap)
    	double sqe = sqrt(1.-e*e)
    	double sqgma = sqrt(GM*a)
    	double xfac1 = a*(ccap-e)
    	double xfac2 = a*sqe*scap
    	double ri = 1./(a*(1.-e*ccap))
    	double vfac1 = -ri*sqgma * scap
    	double vfac2 =  ri*sqgma*sqe*ccap

    	double sw, cw, so, co, si, ci
    	sincos(w,&sw,&cw)
    	sincos(O,&so,&co)
    	sincos(i,&si,&ci)
    	double d1[] = { cw*co-sw*so*ci, cw*so+sw*co*ci, sw*si}
    	double d2[] = {-sw*co-cw*so*ci,-sw*so+cw*co*ci, cw*si}
    	*x  = d1[0]*xfac1+d2[0]*xfac2
    	*y  = d1[1]*xfac1+d2[1]*xfac2
    	*z  = d1[2]*xfac1+d2[2]*xfac2
    	*vx = d1[0]*vfac1+d2[0]*vfac2
    	*vy = d1[1]*vfac1+d2[1]*vfac2
    	*vz = d1[2]*vfac1+d2[2]*vfac2
    }"""
    a, e, i, O, w, M = kepler_coordinates
    cape = mean_to_eccentric_annomaly(e,M)
    scap, ccap = sincos(cape)
    sqe = sqrt(1.0-e*e)
    sqgma = sqrt(GM*a)
    xfac1 = a*(ccap-e)
    xfac2 = a*sqe*scap
    ri = 1.0/(a*(1.0-e*ccap))
    vfac1 = -ri*sqgma * scap
    vfac2 =  ri*sqgma*sqe*ccap
    sw, cw = sincos(w)
    so, co = sincos(O)
    si, ci = sincos(i)
    d1 = [ cw*co-sw*so*ci, cw*so+sw*co*ci, sw*si ]
    d2 = [-sw*co-cw*so*ci,-sw*so+cw*co*ci, cw*si ]
    x  = d1[0]*xfac1+d2[0]*xfac2
    y  = d1[1]*xfac1+d2[1]*xfac2
    z  = d1[2]*xfac1+d2[2]*xfac2
    vx = d1[0]*vfac1+d2[0]*vfac2
    vy = d1[1]*vfac1+d2[1]*vfac2
    vz = d1[2]*vfac1+d2[2]*vfac2
    return ((x,y,z),(vx,vy,vz))

def calc_keplerian_for_cartesian(cartesian,GM):
    """{
    	double a,e,i,O,w,M
    	const double EPSILON = 1.e-8

    	double h[] = {y*vz-z*vy, z*vx-x*vz, x*vy-y*vx}
    	double h2 = h[0]*h[0]+h[1]*h[1]+h[2]*h[2]
    	double hh = sqrt(h2)
    	i = acos(h[2]/hh)
    	double fac = sqrt(h[0]*h[0]+h[1]*h[1])/hh
    	double u

    	if(fac<EPSILON)
    	{
    		O = 0.
    		u = atan2(y,x)
    		if(fabs(i-M_PI)<10.*EPSILON) u = -u
    	}
    	else
    	{
    		O = atan2(h[0],-h[1])
    		u = atan2(z/sin(i),x*cos(O)+y*sin(O))
    	}
    	if(O<0.) O += 2.*M_PI
    	if(u<0.) u += 2.*M_PI
    	double r = sqrt(x*x+y*y+z*z)
    	double energy = (vx*vx+vy*vy+vz*vz)*0.5-GM/r

    	if(fabs(energy*r/GM)<sqrt(EPSILON))
    	{ // Parabola
    		a = 0.5*h2/GM
    		e = 1.
    		double ww = acos(2.*a/r-1.)
    		if(vx*x+vy*y+vz*z<0.) w = 2.*M_PI-w
    		double tmpf = tan(0.5*w)
    		M = tmpf*(1.+tmpf*tmpf/3.)
    		w = u-ww
    		if(w<0.) w+= 2.*M_PI
    		w -= round(w/(2.*M_PI))*2.*M_PI
    	}
    	else if (energy<0)
    	{ // Elipse
    		a = -0.5*GM/energy
    		fac = 1.-h2/(GM*a)
    		double ww, cape
    		if(fac>EPSILON)
    		{
    			e = sqrt(fac)
    			double face = (a-r)/(a*e)
    			if(face>1.) cape = 0.
    			else if (face>-1.) cape = acos(face)
    			else cape = M_PI

    			if(vx*x+vy*y+vz*z<0.) cape = 2.*M_PI-cape
    			double cw = (cos(cape)-e)/(1.-e*cos(cape))
    			double sw = sqrt(1.-e*e)*sin(cape)/(1.-e*cos(cape))
    			ww = atan2(sw,cw)
    			if(ww<0.) ww += 2.*M_PI
    		}
    		else
    		{
    			e = 0.
    			ww = u
    			cape = u
    		}
    		M = cape - e*sin(cape)
    		w = u - ww
    		if(w<0.) w += 2.*M_PI
    		w -= round(w/(2.*M_PI))*2.*M_PI
    	}
    	else if (energy>0)
    	{ // Hyperbola
    		a = 0.5*GM/energy
    		fac = h2/(GM*a)
    		double ww, capf
    		if(fac>EPSILON)
    		{
    			e = sqrt(1.+fac)
    			double tmpf = (a+r)/(a*e)
    			capf = log(tmpf+sqrt(tmpf*tmpf-1.))
    			if(vx*x+vy*y+vz*z<0.) capf = -capf
    			double cw = (e-cosh(capf))/(e*cosh(capf)-1.)
    			double sw = sqrt(e*e-1.)*sinh(capf)/(e*cosh(capf)-1.)
    			ww = atan2(sw,cw)
    			if(ww<0.) ww += 2.*M_PI
    		}
    		else
    		{
    			e = 1.
    			double tmpf = 0.5*h2/GM
    			ww = acos(2.*tmpf/r-1.)
    			if(vx*x+vy*y+vz*z<0.) ww = 2.*M_PI-ww
    			tmpf = (a+r)/(a*e)
    			capf = log(tmpf+sqrt(tmpf*tmpf-1.))
    		}
    		M = e*sinh(capf)-capf
    		w = u - ww
    		if(w<0.) w+=2.*M_PI
    		w -= round(w/(2.*M_PI))*2.*M_PI
    	}

    	*pa = a
    	*pe = e
    	*pi = i
    	*pO = O
    	*pw = w
    	*pM = M
    }"""
    (x,y,z), (vx,vy,vz) = cartesian
    EPSILON = 1.e-8
    h = [y*vz-z*vy, z*vx-x*vz, x*vy-y*vx]
    h2 = h[0]*h[0]+h[1]*h[1]+h[2]*h[2]
    hh = sqrt(h2)
    i = acos(h[2]/hh)
    fac = sqrt(h[0]*h[0]+h[1]*h[1])/hh
    if(fac<EPSILON):
        O = 0.0
        u = atan2(y,x)
        if(abs(i-M_PI)<10.0*EPSILON):
            u = -u
    else:
        O = atan2(h[0],-h[1])
        u = atan2(z/sin(i),x*cos(O)+y*sin(O))
    if(O<0.):
        O += 2.*M_PI
    if(u<0.):
        u += 2.*M_PI
    r = sqrt(x*x+y*y+z*z)
    energy = (vx*vx+vy*vy+vz*vz)*0.5-GM/r

    if(fabs(energy*r/GM)<sqrt(EPSILON)):
        # Parabola
        a = 0.5*h2/GM
        e = 1.
        ww = acos(2.*a/r-1.)
        if(vx*x+vy*y+vz*z<0.):
            ww = 2.*M_PI-ww
        tmpf = tan(0.5*ww)
        M = tmpf*(1.+tmpf*tmpf/3.)
        w = u-ww
        if(w<0.):
            w+= 2.*M_PI
        w -= round(w/(2.*M_PI))*2.*M_PI
    elif (energy<0):
        # Elipse
        a = -0.5*GM/energy
        fac = 1.-h2/(GM*a)
        if(fac>EPSILON):
            e = sqrt(fac)
            face = (a-r)/(a*e)
            if(face>1.):
                cape = 0.
            elif (face>-1.):
                cape = acos(face)
            else:
                cape = M_PI

            if(vx*x+vy*y+vz*z<0.):
                cape = 2.*M_PI-cape
            cw = (cos(cape)-e)/(1.-e*cos(cape))
            sw = sqrt(1.-e*e)*sin(cape)/(1.-e*cos(cape))
            ww = atan2(sw,cw)
            if(ww<0.):
                ww += 2.*M_PI
        else:
            e = 0.
            ww = u
            cape = u

        M = cape - e*sin(cape)
        w = u - ww
        if(w<0.):
            w += 2.*M_PI
        w -= round(w/(2.*M_PI))*2.*M_PI
    elif (energy>0):
        # Hyperbola
        a = 0.5*GM/energy
        fac = h2/(GM*a)
        if(fac>EPSILON):
            e = sqrt(1.+fac)
            tmpf = (a+r)/(a*e)
            capf = log(tmpf+sqrt(tmpf*tmpf-1.))
            if(vx*x+vy*y+vz*z<0.):
                capf = -capf
            cw = (e-cosh(capf))/(e*cosh(capf)-1.)
            sw = sqrt(e*e-1.)*sinh(capf)/(e*cosh(capf)-1.)
            ww = atan2(sw,cw)
            if(ww<0.):
                ww += 2.*M_PI
        else:
            e = 1.
            tmpf = 0.5*h2/GM
            ww = acos(2.*tmpf/r-1.)
            if(vx*x+vy*y+vz*z<0.):
                ww = 2.*M_PI-ww
            tmpf = (a+r)/(a*e)
            capf = log(tmpf+sqrt(tmpf*tmpf-1.))
        M = e*sinh(capf)-capf
        w = u - ww
        if(w<0.):
            w+=2.*M_PI
        w -= round(w/(2.*M_PI))*2.*M_PI

    return KeplerianCoordinates(a,e,i,O,w,M)

