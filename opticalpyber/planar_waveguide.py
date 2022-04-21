"""
Planar Waveguide model

Author: Marcos Moura
Email: gabbyru2@gmail.com

Suported Model:

              cladding(n3)
---------------------------------------
               core (n1)
_______________________________________
             substrate (n2)

Description:
    Model of slab waveguide. This model calculates the effective index of the waveguide and retrieve the dispersion and
    electric fields distribution for the polarizations TE and TM

References:
    [1] Okamoto, K. (2021). Fundamentals of optical waveguides. Elsevier.
"""

import numpy as np
import math

from opticalpyber.utils.bisection_method import Bisection
from opticalpyber.utils.constants import PI

from opticalpyber.waveguide import Waveguide

class PlanarWaveguide(Waveguide):
    def __init__(self, n1: float, n2: float, n3: float, d: float):
        super().__init__(n1,n2,d,n3=n3)

        # dispersion field properties
        self._amplitude: float = None # amplitude of the wave
        self._phi: float       = None # phase shift of wave

    # dispersion equations
    # ref: eq. 2.27 of Okamoto [1]
    def dispersionTE(self, m: int, b: float) -> float:
        if(type(m) is not int):
            raise TypeError("Mode is always a integer, because is the subindex of TE mode")
        return (self.V/2)*math.sqrt(1-b)-m*PI/2 - math.atan(math.sqrt(b/(1-b)))

    # ref: eq. 2.28 of Okamoto [1]
    def dispersionTM(self, m: int, b: float) -> float:
        if(type(m) is not int):
            raise TypeError("Mode is always a integer, because is the subindex of TM mode")
        return (self.V/2)*math.sqrt(1-b)-m*PI/2 - math.atan((self.n1/self.n2)**2*(math.sqrt(b/(1-b))))

    # get Y-axis dispersion equations
    def getDispersionSubstrate(self) -> np.array:
        x = np.arange(-2*self.lbd,-self.a,self.dx)
        k = math.sqrt(self.k0**2 * self.n1**2 - self.beta**2)
        ksi = math.sqrt(self.beta**2 - self.k0**2 * self.n2**2)

        Ey1 = self.amplitude*math.cos((k)*self.a - self.phi)*np.exp(ksi*(x+self.a)) #x<-a
        return Ey1

    def getDispersionCore(self) -> np.array:
        x = np.arange(-self.a,self.a,self.dx)
        k = math.sqrt(self.k0**2* self.n1**2 - self.beta**2)

        Ey2= self.amplitude*np.cos(k*x + self.phi) #-a<x<a
        return Ey2

    def getDispersionCladding(self) -> np.array:
        x = np.arange(self.a,2*self.lbd,self.dx)

        k = math.sqrt(self.k0**2 * self.n1**2 - self.beta**2)
        sigma = math.sqrt(self.beta**2 - self.k0**2 * self.n3**2)

        Ey3 = self.amplitude*math.cos(k*self.a + self.phi)*np.exp(-sigma*(x-self.a)) #x>a
        return Ey3

    # calculate neff
    def calculateNeff(self, mode: str = "TE", mode_number: int = 0) -> None:
        if(type(mode_number) is not int):
            raise TypeError("Mode Number is always a integer, because is the subindex of the polarization mode")
        self.k0 = 2*PI/self.lbd
        self.V = self.k0*(2*self.a)*math.sqrt(self.n1**2 - self.n2**2)

        # to solve dispersion equations
        bisection = Bisection(mode=mode_number)
        if(mode=="TE"):
            bisection.callback = self.dispersionTE # eigenvalue solver
        elif(mode=="TM"):
            bisection.callback = self.dispersionTM # eigenvalue solver
        else:
            raise TypeError("Invalid mode! You can choose between transversal electric(TE) or transversal Magnetic")

        b = bisection.findEigenvalue()

        # updating neff and gamma
        self.neff = math.sqrt(b*(self.n1**2-self.n2**2)+self.n2**2)
        self.beta = self.neff*self.k0

    # get dispersion field
    def calculateDispersionField(self, lbd: float, mode: str = "TE", mode_number: int = 0, amplitude: float=1.0, phi: float=0, dx: float=1E-8)->np.array:
        self.lbd = lbd
        self.amplitude = amplitude
        self.phi = phi
        self.dx = dx

        self.calculateNeff(mode=mode, mode_number=mode_number)

        Ey1 = self.getDispersionSubstrate()
        Ey2 = self.getDispersionCore()
        Ey3 = self.getDispersionCladding()

        return np.concatenate((Ey1,Ey2,Ey3))
