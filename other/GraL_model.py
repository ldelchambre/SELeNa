# -*- coding: utf-8 -*-
"""
Classes dedicated to gravitational lens modeling
"""

from __future__ import division

__author__ = 'O. Wertz'

from numpy import pi
from numpy import array, shape 
from numpy import mod, sqrt, arctanh, arctan, arctan2, cos, sin


class GraL_LensModel(object):
    """
    Main class for lens models where we define global methods:
    + shear
    + SDF
    + backwards_source
    + rotation
    
    ALl models inherit from this class.
    """
    def __init__(self, modelParameters=array([1,0])):
        """
        Mandatory parameters:
        0: tE
        1: varpi
        """
        self._modelParameters = modelParameters
        
    @property    
    def modelParameters(self):
        """
        Model parameters getter.
        """
        return self._modelParameters 

    @modelParameters.setter
    def modelParameters(self, value):
        """
        Model parameters setter.
        """
        if not shape(value) == shape(self._modelParameters):
            text_error_0 = """The new and old model parameter vectors must have
            the same shape: {}"""
            raise IndexError(text_error_0.format(shape(self._modelParameters)))
            
        self._modelParameters = value
             
    def shear(self, t2_, t1_, g, w_arb, varpi):
        """
        Define the constant external shear. It is expressed in the coordinates
        attached to the lens. For the (N)SIE, it corresponds to the semi-major 
        and semi-minor axis of the elliptical iso-contour of the surface mass 
        density (\kappa). 
        
        Parameters
        ----------
        t1_, t2_ : float or numpy array
            Position in the image plane.
        g, w_arb : float
            Magnitude and orientation (rad) of the external shear. The 
            orientation is measured with respect to the horizontal axis 
            (right).
        varpi : float
            Orientation (rad) of the model. If varpi = 0, the coordinate axis
            are (x,y).
        """
        t_ = sqrt(t2_**2+t1_**2)
        phi_ = arctan2(t2_,t1_)
        
        gamma1 = g * t_ * cos(2*w_arb-phi_ + varpi)
        gamma2 = g * t_ * sin(2*w_arb-phi_ + varpi)  

        return (gamma1, gamma2)        

    def rotation(self, t2_, t1_, varpi=0, tE=1, normed=True):
        """
        Change of basis from (x,y) into varpi-rotated basis. 
        
        Parameters
        ----------
        t1_, t2_ : float or numpy array
            Coordinates to rotate.
        varpi : float
            Orientation (rad) of the new basis. If varpi = 0, the coordinate 
            axis are unchanged.
        tE : float > 0
            Angular radius of the Einstein ring (typical normalising factor).
        """
        x1_, x2_ = t1_/tE, t2_/tE        
        
        x1 = x1_*cos(varpi)-x2_*sin(varpi)
        x2 = x2_*cos(varpi)+x1_*sin(varpi) 
        
        x = sqrt(x1**2+x2**2)
        phi = arctan2(x2,x1)

        if normed:
            return x1, x2, x, phi 
        else:
            return x1*tE, x2*tE, x*tE, phi                       

    def SDF(self, position, source, deflectionAngle):
        """
        The so-called 'Square Deviation Function' defined by of Schramm & Kaiser
        (1987), basically given by:
        
        SDF**2 = (y1 - x1 - alpha1)**2 + (y2 - x2 - alpha2)**2 ,
        
        where (y1,y2) represents the position of the source projected in the 
        image plane, (x1,x2) a position in the image plane, and (alpha1, alpha2) 
        the angular deflection angle vector projected in the image plane, all
        in cartesian coordinates. For a given model and position of the source,
        the lensed image positions are the ones for which SDF = 0.

        Parameters
        ----------
        position : list or tuple
            Position(s) in the image plane.
        source : list or tuple
            Position of the source in polar coordinate (angular separation,
            azimuth with respect to the horizontal axis).
        deflectionAngle : callable
            Function which returns the value of the deflection angle at a given
            position in the image plane.
            
        NB: This method is overrided in the sub-classes which inherit from
        GraL_LensModel.
        """
        t1_, t2_ = position[0], position[1]
        
        tE = self._modelParameters[0]
        varpi = self._modelParameters[1] 
        
        alpha1, alpha2 = deflectionAngle(self, t2_, t1_)        
        t1, t2, _, _ = GraL_LensModel.rotation(self, t2_, t1_, varpi, tE, 
                                               normed=False)
        
        sk1 = source[0]*cos(source[1]+varpi) - t1 - alpha1   
        sk2 = source[0]*sin(source[1]+varpi) - t2 - alpha2        
                       
        return sqrt(sk1**2 + sk2**2)
        
        
    def backwards_source(self, t2_, t1_, deflectionAngle):
        """
        Read lens equation 'backwards': mapping from the lensed image position 
        to the source position. Basically, for a given position in the image
        plane and a deflection law, it returns the source position at the 
        origine of the image located at the given position.
        
        Parameters
        ----------
        t1_, t2_ : float or numpy array
            Position in the image plane.  
        deflectionAngle : callable
            Function which returns the value of the deflection angle at a given
            position in the image plane.            
        """
        tE = self._modelParameters[0]
        varpi = self._modelParameters[1]
        
        alpha1, alpha2 = deflectionAngle(self, t2_, t1_)
        
        t1, t2, _, _ = GraL_LensModel.rotation(self, t2_, t1_, varpi, tE, 
                                               normed=False) 
        
        s1, s2 = t1 + alpha1, t2 + alpha2
        
        tS = sqrt(s1**2 + s2**2)
        t_arb = arctan2(s2,s1)-varpi
        
        return tS, mod(t_arb,2*pi)
        
                    
class GraL_NSIE(GraL_LensModel):
    """
    Define the Non-Singular Isothermal Ellipsoid model.
    
    Physical parameters:
    0: tE = Angular radius of the Einstein ring (typical normalising factor).
    1: varpi = Orientation (rad) of the model.
    2: f = elliptical axis ratio. Axi-symmetric model corresponds to f=1.
    3: bc = angular radius of the core.   
    """
    def __init__(self, modelParameters):
        """
        Parameters: tE, varpi, f, bc. 
        """  
        super(GraL_NSIE, self).__init__(modelParameters)
        
                                                                    
    def deflectionAngle(self, t2_, t1_):
        """
        Define the deflection angle which includes the physical information of 
        the lens. 
        
        Parameters:
        -----------
        t1_, t2_ : float or numpy array
            Position in the image plane.        
        """
        f = self._modelParameters[2]
        bc = self._modelParameters[3]
        tE = self._modelParameters[0]        
        varpi = self._modelParameters[1]
        
        fp = sqrt(1-f**2)        
        rc = bc/tE        
        x1, x2, _, _ = GraL_NSIE.rotation(self, t2_, t1_, varpi, tE, normed=True)        
        
        alpha1 = tE*(-sqrt(f)/fp*arctanh(fp*x1/(sqrt(x1**2+f**2*x2**2+rc**2)
                                         +f*rc)))
        alpha2 = tE*(-sqrt(f)/fp*arctan(fp*x2/(sqrt(x1**2+f**2*x2**2+rc**2)
                                        +1/f*rc)))

        return alpha1, alpha2           

    def SDF(self, position, source):
        """
        Override the SDF method from GraL_LensModel to include the deflection
        angle.
        """
        return super(GraL_NSIE, self).SDF(position, source, 
                                          GraL_NSIE.deflectionAngle)

    def backwards_source(self, t2_, t1_):
        """
        Override the source method from GraL_LensModel to include the deflection
        angle.
        """
        return super(GraL_NSIE, self).backwards_source(t2_, t1_, 
                                                     GraL_NSIE.deflectionAngle)                                          
                                                

class GraL_NSIEG(GraL_NSIE):
    """
    Define the Non-Singular Isothermal Ellipsoid model in presence of external
    shear.
    
    Physical parameters:
    0: tE = Angular radius of the Einstein ring (typical normalising factor).
    1: varpi = Orientation (rad) of the model.
    2: f = elliptical axis ratio. Axi-symmetric model corresponds to f=1.
    3: bc = angular radius of the core.
    4: g = intensity of the external shear.
    5: w_arb = orientation of the shear.       
    """
    def __init__(self, modelParameters):
        """
        Parameters: tE, varpi, f, bc, g, w_arb. 
        """  
        super(GraL_NSIE, self).__init__(modelParameters)
                                                                        
    def deflectionAngle(self, t2_, t1_):
        """
        Define the deflection angle which includes the physical information of 
        the lens. 
        
        Parameters:
        -----------
        t1_, t2_ : float or numpy array
            Position in the image plane. 
        """
        varpi = self._modelParameters[1]
        g = self._modelParameters[4]
        w_arb = self._modelParameters[5]
       
        alpha1, alpha2 = super(GraL_NSIEG, self).deflectionAngle(t2_, t1_)
        gamma1, gamma2 = GraL_NSIEG.shear(self, t2_, t1_, g, w_arb, varpi)
        
        return alpha1+gamma1, alpha2+gamma2           

    def SDF(self, position, source):
        """
        Override the SDF method from GraL_LensModel to include the deflection
        angle.        
        """
        return super(GraL_NSIE, self).SDF(position, source, 
                                          GraL_NSIEG.deflectionAngle) 

    def backwards_source(self, t2_, t1_):
        """
        Override the source method from GraL_LensModel to include the deflection
        angle.
        """
        return super(GraL_NSIE, self).backwards_source(t2_, t1_, 
                                                    GraL_NSIEG.deflectionAngle)                                               

             
############################################################################### 
# PROCDURAL VERSION
###############################################################################        
def rotation(t2_, t1_, varpi, tE, normed=True):
    """
    """
    x1_, x2_ = t1_/tE, t2_/tE        
    
    x1 = x1_*cos(varpi)-x2_*sin(varpi)
    x2 = x2_*cos(varpi)+x1_*sin(varpi) 
    
    x = sqrt(x1**2+x2**2)
    phi = arctan2(x2,x1)

    if normed:
        return x1, x2, x, phi 
    else:
        return x1*tE, x2*tE, x*tE, phi

def shear(t2_, t1_, g, w_arb, varpi):
    """
    """
    t_ = sqrt(t2_**2+t1_**2)
    phi_ = arctan2(t2_,t1_)
    
    gamma1 = g * t_ * cos(2*w_arb-phi_ + varpi)
    gamma2 = g * t_ * sin(2*w_arb-phi_ + varpi)  

    return (gamma1, gamma2) 

def NSIEG_deflectionAngle(modelParameters, t2_, t1_):
    """
    """  
    f = modelParameters[2]
    bc = modelParameters[3]
    tE = modelParameters[0]        
    varpi = modelParameters[1]
    g = modelParameters[4]
    w_arb = modelParameters[5]
    
    fp = sqrt(1-f**2)        
    rc = bc/tE        
    x1, x2, _, _ = rotation(t2_, t1_, varpi, tE, normed=True)    

    gamma1, gamma2 = shear(t2_, t1_, g, w_arb, varpi)    
    
    alpha1 = tE*(-sqrt(f)/fp*arctanh(fp*x1/(sqrt(x1**2+f**2*x2**2+rc**2)
                                     +f*rc))) + gamma1
    alpha2 = tE*(-sqrt(f)/fp*arctan(fp*x2/(sqrt(x1**2+f**2*x2**2+rc**2)
                                    +1/f*rc))) + gamma2 

    return alpha1, alpha2    
            
def NSIEG_backwards_source(modelParameters, t2_, t1_):
    """
    Parameters (modelParameters)
    0: tE
    1: varpi
    2: f
    3: bc
    4: g
    5: w_arb    
    """
    tE = modelParameters[0]        
    varpi = modelParameters[1]
    
    x1, x2, _, _ = rotation(t2_, t1_, varpi, tE, normed=True) 
    t1, t2 = x1*tE, x2*tE
       
    alpha1, alpha2 = NSIEG_deflectionAngle(modelParameters, t2_, t1_)
   
    s1, s2 = t1 + alpha1, t2 + alpha2        
    tS = sqrt(s1**2 + s2**2)
    t_arb = arctan2(s2,s1)-varpi
        
    return tS, mod(t_arb,2*pi) 

def NSIEG_SDF(modelParameters, t2_, t1_, source):
    """
    The so-called 'Square Deviation Function' defined by of Schramm & Kaiser
    (1987).
    
    Parameters (modelParameters)
    0: tE
    1: varpi
    2: f
    3: bc
    4: g
    5: w_arb
    """
    tE = modelParameters[0]
    varpi = modelParameters[1] 
    
    alpha1, alpha2 = NSIEG_deflectionAngle(modelParameters, t2_, t1_)        
    t1, t2, _, _ = rotation(t2_, t1_, varpi, tE, normed=False)
    
    sk1 = source[0]*cos(source[1]+varpi) - t1 - alpha1   
    sk2 = source[0]*sin(source[1]+varpi) - t2 - alpha2        
                   
    return sqrt(sk1**2 + sk2**2)    
    
###############################################################################    
