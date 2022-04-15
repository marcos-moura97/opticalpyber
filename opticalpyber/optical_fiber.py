"""
Optical Fiber model

Author: Marcos Moura
Email: gabbyru2@gmail.com

Suported Models:

         , - - - - - ,
     , '    cladding   ' ,
   ,          (n2)         ,
  ,          _ _ _          ,
 ,         ,       ,         ,
 ,        ,  core   ,         ,
 ,         , (n1)  ,          ,
  ,          - - -           ,
   ,                       ,
     ,                  ,
        ' - - - - - - '

Description:
    Model of an optical fiber waveguide. This model calculate the LP Modes that can pass through an step-index fiber.
    The code  find all the L modes for an given M mode and can calculate a meshgrid with the dispersion of the modes.
    The dispersion of the fields are solved in cilindrical coordinates.

References:
    [1] Okamoto, K. (2021). Fundamentals of optical waveguides. Elsevier.
"""

import numpy as np
import math

from scipy.special import jv

from utils.constants import PI
from waveguide import Waveguide

class OpticalFiber(Waveguide):
    def __init__(self, n1: float, n2: float, d: float):
        super().__init__(n1,n2,d)

        self._lbd: float       = None # wavelength of simulation
        self._k0: float        = None # wavenumber
        self._V: float         = None # normalized frequency

        # cilyndrical coordinates (r, theta and z)
        self._r: np.ndarray    = None # dots array on axis R
        self._theta: np.ndarray= None # dots array on angle theta
        self._z: np.ndarray    = None # dots array on axis Z

        # fiber calculation mode
        self._m: int          = None
        self._l:  int         = None

    def updateCoordinates(self):
        # cilyndrical coordinates
        self.r=np.arange(0,self.a,self.dx)
        self.theta=np.linspace(0,2*PI,len(self.r))
        [self.r,self.theta]=np.meshgrid(self.r,self.theta)
        print(len(self.r))
    # find m modes
    def findMModes(self, lbd: float=None, mode: int=0) -> np.array:
        if(type(mode) is not int):
            raise TypeError("Mode Type is always a integer, because is the subindex of M mode")
        self._m = mode

        if not self.lbd and not lbd:
            raise RuntimeError("Wavelength not seted yet")
        elif lbd:
            self.lbd = lbd
            self.k0 = 2*PI/self.lbd
            self.V = self.k0*(2*self.a)*math.sqrt(self.n1**2 - self.n2**2)

        m = []
        for x in range(1,11):
            for t in np.linspace((x-1)*PI,(x+1)*PI,500):
                s = jv(self._m,t)

                if s==0 or abs(s)<=0.01:
                    m.append(t)
            if x*PI>self.V:
                break
        return m

    # calculate coordinates of dispersion field on each mode
    def getDispersionCartesian(self, m0: float)-> tuple[np.array]:
        z= (jv(self._m,m0*self.r/self.a)*np.cos(self._m*self.theta))
        z = list(map(lambda x: x**2, z))

        #convert from cilindric to cartesian
        x=self.r*np.cos(self.theta)
        y=self.r*np.sin(self.theta)

        return (x, y, z)



    # calculate Dispersion Fields
    def calculateDispersionField(self, mode: float, phi: float=0, dx: float=1E-7) -> tuple[np.array]:

        self.phi = phi
        self.dx = dx

        # make an meshgrid of r and phi
        self.updateCoordinates()

        cartesian_dispersion = self.getDispersionCartesian(mode)

        return cartesian_dispersion


