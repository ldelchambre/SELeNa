# -*- coding: utf-8 -*-
"""
@author: Olivier W.

NSIE+shear lens model

Parameters b', e, theta_e and s' corresponds to the ones defined in
"gravlens 1.06: Software for Gravitational Lensing" by C. Keeton, Eq.(3.25)
p23.
URL: http://www.physics.rutgers.edu/~keeton/gravlens/manual.pdf

The parameters of the gravlens code are:
[p[1], p[4], p[5],p[8]] = [b', e, theta_e, s']

The link to the previous parameters f, bc and varpi is defined by:

f --> 1 - e
bc --> (1 - e) * s
varpi --> theta_e

where

b' = b * q * sqrt(2 / (1+q**2))
s' = s * q * sqrt(2 / (1+q**2))
e = 1 - q ==> f = q

In addition, we define epsilon such as q**2 = (1-epsilon)/(1+epsilon) which
leads to :
epsilon = 2/(1+(1-e)**2) - 1

As a result, the deflection angle components defined by Keeton are equivalent
to our previous definition (up to the sign convention and the tE
normalization).
Sign convention:
us     --> source = image + alpha
Keeton --> source = image - alpha

The coordinate system is now fixed to the (North, East) orientation. Therefore,
we can link the parameter varpi to theta_e, and w_arb to theta_g.

"""

from numpy import sqrt, arctanh, arctan, arctan2, cos, sin, log, pi
from numpy import array
from scipy.misc import derivative


def qprime(q):
    return (1-q**2)**(0.5)

def psi(x,y,q,s):
    return (q**2*(s**2+y**2) + x**2)**(0.5)


def alpha_x(x, y, b, q, s):
    """
    First component of the NSIE deflection angle at position (x,y) in the lens
    plane without any rotation of the lens (theta_e = 0).
    """
    if q == 1:
        return b*x/(s+(s**2+x**2+y**2)**(0.5))
    else:
        return b*q/qprime(q)*arctanh((qprime(q)*x)/(psi(x,y,q,s)+q**2*s))

def alpha_y(x, y, b, q, s):
    """
    Second component of the NSIE deflection angle at position (x,y) in the lens
    plane without any rotation of the lens (theta_e = 0).
    """
    if q == 1:
        return b*y/(s+(s**2+x**2+y**2)**(0.5))
    else:
        return b*q/qprime(q)*arctan((qprime(q)*y)/(psi(x,y,q,s)+s))


def alpha_x_rot_x(x, y, b, q, s, theta_e):
    """
    First component of the NSIE deflection angle at position (x,y) in the lens
    plane.
    """

    r = (x**2+y**2)**(0.5)
    theta = arctan2(y,x)

    xnew = r*cos(theta-theta_e)
    ynew = r*sin(theta-theta_e)

    return alpha_x(xnew,ynew,b,q,s)*cos(theta_e) - alpha_y(xnew,ynew,b,q,s)*sin(theta_e)

def alpha_x_rot_y(y, x, b, q, s, theta_e):
    """
    First component of the NSIE deflection angle at position (x,y) in the lens
    plane. Note that alpha_x_rot_y(y, x, ...) = alpha_x_rot_x(x, y, ...).
    """
    r = (x**2+y**2)**(0.5)
    theta = arctan2(y,x)

    xnew = r*cos(theta-theta_e)
    ynew = r*sin(theta-theta_e)

    return alpha_x(xnew,ynew,b,q,s)*cos(theta_e) - alpha_y(xnew,ynew,b,q,s)*sin(theta_e)

def alpha_y_rot_y(y, x, b, q, s, theta_e):
    """
    Second component of the NSIE deflection angle at position (x,y) in the lens
    plane.
    """
    r = (x**2+y**2)**(0.5)
    theta = arctan2(y,x)

    xnew = r*cos(theta-theta_e)
    ynew = r*sin(theta-theta_e)

    return alpha_x(xnew,ynew,b,q,s)*sin(theta_e) + alpha_y(xnew,ynew,b,q,s)*cos(theta_e)

def alpha_y_rot_x(x, y, b, q, s, theta_e):
    """
    Second component of the NSIE deflection angle at position (x,y) in the lens
    plane. Note that alpha_y_rot_x(x, y, ...) = alpha_y_rot_y(y, x, ...).
    """
    r = (x**2+y**2)**(0.5)
    theta = arctan2(y,x)

    xnew = r*cos(theta-theta_e)
    ynew = r*sin(theta-theta_e)

    return alpha_x(xnew,ynew,b,q,s)*sin(theta_e) + alpha_y(xnew,ynew,b,q,s)*cos(theta_e)


def kappa_NSIE(x, y, b, q, s, theta_e):
    """
    Surface mass density for the NSIE at the position (x,y) located in the lens
    plane.
    """
    r = sqrt(x**2+y**2)
    theta = arctan2(y,x)

    x_new = r*cos(theta-theta_e)
    y_new = r*sin(theta-theta_e)

    return b/(2*sqrt((x_new)**2/q**2 + (y_new)**2 + s**2))


def alpha_xy_shear(x, y, g, theta_g):
    """
    First and second components of the deflection angle associated with the
    external shear.
    """
    r = sqrt(x**2+y**2)
    theta = arctan2(y,x)

    Gx = g * r * cos(theta-2*theta_g)
    Gy =-g * r * sin(theta-2*theta_g)

    return Gx, Gy

def potential_NSIE(x, y, b, q, s, theta_e):
    """
    The NSIE deflection potential evaluated at the position (x,y) in the lens
    plane.
    """
    r = sqrt(x**2+y**2)
    theta = arctan2(y,x)

    x_new = r*cos(theta-theta_e)
    y_new = r*sin(theta-theta_e)

    alpha_x_ = alpha_x(x_new, y_new, b, q, s)
    alpha_y_ = alpha_y(x_new, y_new, b, q, s)

    if s == 0:
        psi = y_new * alpha_y_ + x_new * alpha_x_
    else:
        psi = y_new * alpha_y_ + \
              x_new * alpha_x_ - \
              b*q*s*log(sqrt((sqrt(q**2*(s**2+y_new**2)+x_new**2) + s)**2 + \
              (1-q**2)*y_new**2)) + b*q*s*log((1+q)*s)

    return psi

def potential_shear(x, y, g, theta_g):
    """
    The deflection potential, associated with the external shear, evaluated at
    the position (x,y) in the lens plane.
    """
    return 0.5*g*cos(2*theta_g)*(x**2-y**2) + g*x*y*sin(2*theta_g)

def potential_lens_2ndDerivative(x, y, b, q, s, theta_e):
    """
    Evaluation of the second derivative of the NSIE deflection potential at the
    (x,y) position in the lens plane.
    """
    xx = derivative(alpha_x_rot_x, x, args=[y,b,q,s,theta_e], dx=1e-8)
    yy = derivative(alpha_y_rot_y, y, args=[x,b,q,s,theta_e], dx=1e-8)
    xy = derivative(alpha_x_rot_y, y, args=[x,b,q,s,theta_e], dx=1e-8)
    #yx = derivative(alpha_y_rot_x, x, args=[y,b,q,s,theta_e], dx=1e-8)
    yx = xy

    return array([[xx,xy],[yx,yy]])


def potiential_shear_2ndDerivative(x, y, g, theta_g):
    """
    The second derivative of the shear deflection potential at the (x,y)
    position in the lens plane.
    """
    return array([[g*cos(2*theta_g),g*sin(2*theta_g)],\
                  [g*sin(2*theta_g),-g*cos(2*theta_g)]])


def amplification_NSIEG(x, y, b, q, s, theta_e, g, theta_g):
    """
    Amplification of the lensed image (x,y) of a point-like source for a NSIE
    lens model + external shear.

    Procedure to obtain the amplification of a lensed image:
    1) Define a NSIE lens model (b, q, s, theta_e) and the external shear (g, theta_g)
    2) Define the position (in the lens plane) of a point-like source (x_s, y_s)
    3) Solve the lens equation (x_s = x_im - alpha_x, y_s = y_im - alpha_y) to
       obtain the position (x_im, y_im) of the lensed image(s)
    4) Determine the amplification of a lensed image:
       mu = amplification_NSIEG(x_im, y_im, b, q, s, theta_e, g, theta_g)
    """
    kappa = kappa_NSIE(x, y, b, q, s, theta_e)

    pot_lens_2ndD = potential_lens_2ndDerivative(x, y, b, q, s, theta_e)
    pot_shear_2ndD = potiential_shear_2ndDerivative(x, y, g, theta_g)

    mu_inv = (1-kappa)**2 - 0.25*((pot_lens_2ndD[0,0]+pot_shear_2ndD[0,0]) - \
             (pot_lens_2ndD[1,1]+pot_shear_2ndD[1,1]))**2 - \
             (pot_lens_2ndD[0,1]+pot_shear_2ndD[0,1])**2

    return 1/mu_inv

# -----------------------------------------------------------------------------
def NSIEG_deflectionAngle(modelParameters, y, x):
    """
    Backwards compability for NSIEG_deflectionAngle
    """
    b, theta_e, q, s, g, theta_g = modelParameters
    ax_lens = alpha_x_rot_x(x, y, b, q, s, theta_e)
    ay_lens = alpha_y_rot_x(x, y, b, q, s, theta_e)
    ax_shear, ay_shear = alpha_xy_shear(x, y ,g, theta_g)

    return ax_lens+ax_shear, ay_lens+ay_shear


def NSIEG_backwards_source(modelParameters, y, x, coord='polar'):
    """
    Model parameters : [b, theta_e, q, s, g, theta_g]
    """
    b, theta_e, q, s, g, theta_g = modelParameters

    ax_lens = alpha_x_rot_x(x, y, b, q, s, theta_e)
    ay_lens = alpha_y_rot_x(x, y, b, q, s, theta_e)

    ax_shear, ay_shear = alpha_xy_shear(x, y, g, theta_g)

    source_x = x - (ax_lens+ax_shear)
    source_y = y - (ay_lens+ay_shear)

    if coord == 'cart':
        return source_x, source_y
    elif coord == 'polar':
        return (source_x**2+source_y**2)**(0.5), arctan2(source_y, source_x)

def shear(y, x, g, theta_g, varpi=0):
    """
    """
    return alpha_xy_shear(x, y, g, theta_g)


def NSIEG_SDF(modelParameters, y, x, source):
    """
    The so-called 'Square Deviation Function' defined by of Schramm & Kaiser
    (1987).

    modelParameters = [b, theta_e, q, s, g, theta_g]
    """
    b, theta_e, q, s, g, theta_g = modelParameters

    ax_lens = alpha_x_rot_x(x, y, b, q, s, theta_e)
    ay_lens = alpha_y_rot_x(x, y, b, q, s, theta_e)

    ax_shear, ay_shear = alpha_xy_shear(x, y, g, theta_g)

    sk1 = source[0] - x + (ax_lens + ax_shear)
    sk2 = source[1] - y + (ay_lens + ay_shear)

    return sqrt(sk1**2 + sk2**2)

if __name__ == '__main__':
    # -*- coding: utf-8 -*-
    """
    Test of GraL_NSIEg.py
    
    @author: Olivier W.
    
    Parameters b', e, theta_e and s' corresponds to the ones defined in 
    "gravlens 1.06: Software for Gravitational Lensing" by C. Keeton, Eq.(3.25) 
    p23. 
    URL: http://www.physics.rutgers.edu/~keeton/gravlens/manual.pdf
    
    The parameters of the gravlens code are:
    [p[1], p[4], p[5],p[8]] = [b', e, theta_e, s']
    
    The link to the previous parameters f, bc and varpi is defined by:
    
    f --> 1 - e
    bc --> (1 - e) * s
    varpi --> theta_e
    
    where 
    
    b' = b * q * sqrt(2 / (1+q**2))
    s' = s * q * sqrt(2 / (1+q**2))
    e = 1 - q ==> f = q
    
    In addition, we define epsilon such as q**2 = (1-epsilon)/(1+epsilon) which 
    leads to :
    epsilon = 2/(1+(1-e)**2) - 1
    
    As a result, the deflection angle components defined by Keeton are equivalent
    to our previous definition (up to the sign convention and the tE 
    normalization).
    Sign convention:
    us     --> source = image + alpha
    Keeton --> source = image - alpha
      
    The coordinate system is now fixed to the (North, East) orientation. Therefore,
    we can link the parameter varpi to theta_e, and w_arb to theta_g.
    
    """ 
    # Parameters
    # Position of a point in the image plane
    # (actually the position (-0.2773999, 1.982248) corresponds to one of the 
    # lensed image position of the source located in (x_s, y_s)=(0.1, 0.08)) for 
    # the NSIEG model defined here below.
    
    # Obtained from gravlens with input
    #   maingrid 1 -2.5 2.5 100 -2.5 2.5 100
    #   setlens 1 1
    #   alpha 1.31 0 0 0.77 121. 0.17 62 0.13 0 1
    #   0 0 0 0 0 0 0 0 0 0
    #   findimg 0.1 0.08
    # x, y = -0.2773999, 1.982248
    
    # x_s, y_s = (0.1, 0.08)
    
    # Model (Keeton parameters)
    # bp = 1.31
    # e = 0.77
    # theta_e = 121.*(pi/180.)
    # sp = 0.13
    # g = 0.17
    # theta_g = 62.*(pi/180.)
    
    x, y = 1.206793e+00, 4.915795e-01
    x_s, y_s = (0.1, 0.1)
    bp = 1.
    e = 0.5
    theta_e = 0
    sp = 0.1
    g = 0.1
    theta_g = 0.1*pi
    
    # The corresponding b, q, s parameters
    q = 1-e; qp = sqrt(1-q**2)
    ep = 2/(1+(1-e)**2) - 1
    b = bp/sqrt(1-ep)
    s = sp/sqrt(1-ep)
    
    print(q, ep, b, s
    
    # Let's compute some quantities
    kappa = kappa_NSIE(x, y, b, q, s, theta_e)
    
    ax_lens = alpha_x_rot_x(x, y, b, q, s, theta_e)
    ay_lens = alpha_y_rot_x(x, y, b, q, s, theta_e)
    
    ax_shear, ay_shear = alpha_xy_shear(x, y, g, theta_g)
    
    pot_lens = potential_NSIE(x, y, b, q, s, theta_e)
    pot_shear = potential_shear(x, y, g, theta_g)
    
    pot_lens_2ndD = potential_lens_2ndDerivative(x, y, b, q, s, theta_e)
    pot_shear_2ndD = potiential_shear_2ndDerivative(x, y, g, theta_g)
    
    mu = amplification_NSIEG(x, y, b, q, s, theta_e, g, theta_g)
    
    source_x, source_y = NSIEG_backwards_source([b, theta_e, q, s, g, theta_g], y, x, coord='cart')
    sdf = NSIEG_SDF([b,theta_e,q,s,g,theta_g], y, x, (x_s,y_s))
    
    # print the results
    print('')
    print('---------------------------------')
    print('Parameters')
    print('')
    print('x = {}, y = {}'.format(x,y))
    print('')
    print('---------------------------------')
    print('kappa')
    print('')
    print('kappa = {}'.format(kappa))
    print('')
    print('---------------------------------')
    print('Deflection')
    print('')
    print('alpha_x_lens = {}'.format(ax_lens))
    print('alpha_y_lens = {}'.format(ay_lens))
    print('')
    print('alpha_x_shear = {}'.format(ax_shear))
    print('alpha_y_shear = {}'.format(ay_shear))
    print('')
    print('alpha_x_total = {}'.format(ax_lens+ax_shear))
    print('alpha_y_total = {}'.format(ay_lens+ay_shear))
    print('')
    print('---------------------------------')
    print('Potential')
    print('pot lens: {}'.format(pot_lens))
    print('pot shear: {}'.format(pot_shear))
    print('pot total: {}'.format(pot_lens+pot_shear))
    print('')
    print('---------------------------------')
    print('Second derivative of the potential')
    print('')
    print('Lens')
    print('pot_xx : {}'.format(pot_lens_2ndD[0,0]))
    print('pot_yy : {}'.format(pot_lens_2ndD[1,1]))
    print('pot_xy : {}'.format(pot_lens_2ndD[0,1]))
    print('pot_yx : {}'.format(pot_lens_2ndD[1,0]))
    print('')
    print('Shear')
    print('pot_xx : {}'.format(pot_shear_2ndD[0,0]))
    print('pot_yy : {}'.format(pot_shear_2ndD[1,1]))
    print('pot_xy : {}'.format(pot_shear_2ndD[0,1]))
    print('pot_yx : {}'.format(pot_shear_2ndD[1,0]))
    print('')
    print('Total')
    print('pot_xx : {}'.format(pot_lens_2ndD[0,0] + pot_shear_2ndD[0,0]))
    print('pot_yy : {}'.format(pot_lens_2ndD[1,1] + pot_shear_2ndD[1,1]))
    print('pot_xy : {}'.format(pot_lens_2ndD[0,1] + pot_shear_2ndD[0,1]))
    print('pot_yx : {}'.format(pot_lens_2ndD[1,0] + pot_shear_2ndD[1,0]))
    
    print('---------------------------------')
    print('Amplification')
    print('')
    print('(x,y)=({},{}) \t mu: {}'.format(x,y,mu))
    print('')
    print('---------------------------------')
    print('Check that (x,y) is a lensed image')
    print('')
    print('True position of the source: ({:.3f},{:.3f})'.format(x_s, y_s))
    print('Position of the source from the lens equation: ({:.3f},{:.3f})'.format(source_x, source_y))
    print('SDF function at position (x,y): {}'.format(sdf)
    print(''
    print('---------------------------------'
