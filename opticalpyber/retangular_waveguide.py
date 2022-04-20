"""
Retangular Waveguide model

Author: Marcos Moura
Email: gabbyru2@gmail.com

Suported Models:

            1) Rib                                           2) Ridge
               __________
 cladding(n3) |   rib   | cladding(n3)                           _____________
______________|   (n2)  |______________            cladding(n3) | ridge (n2) | cladding(n3)
---------------------------------------          _______________|____________|___________
               core (n1)                                           core (n1)
_______________________________________          _______________________________________
               substrate (n2)                                    substrate (n2)


Description:
    Model of retangular waveguide. This model uses the effective index model [2] to analize a mode E^x_pq

References:
    [1] Okamoto, K. (2021). Fundamentals of optical waveguides. Elsevier.
    [2] Knox, R. M. and P. P. Toulios. 1970. Integrated circuits for the millimeter through optical frequency range.
        Symposium on Submillimeter Waves, Polytechnic Institute of Brooklyn, pp. 497–516.

"""

import numpy as np

from utils.bisection_method import Bisection
from utils.constants import PI

from waveguide import Waveguide

class RetangularWaveguide(Waveguide):
    def __init__(self, n1: float, n2: float, n3: float, d: float, h: float, t: float, hc: float):
        super().__init__(n1,n2,d,n3=n3)

        # retangular waveguide dimensions
        self._h  = h     # rib bigger height
        self._t  = t     # rib smaller height, if =0, than is a ridge waveguide
        self._hc = hc    # core height

        self._rib_height = None

        # to solve dispersion equations
        self._bisection = Bisection(divsn=(self.n1-self.n2)/100, eps=1e-2, top_value=self.n1, bottom_value=self.n2)

        # solve dispersion equations
        self._bisection.callback = self.dispersionFunction

        # field dispertion



    @property
    def h(self) -> float:
        return self._h

    @property
    def t(self) -> float:
        return self._t

    @property
    def hc(self) -> float:
        return self._hc

    @Waveguide.neff.getter
    def neff(self) -> np.ndarray:
        return self._neff

    # setters
    @h.setter
    def h(self, new_h: float) -> None:
        if new_h < 0:
            raise TypeError("Only values bigger than 0 are allowed")
        elif new_h < self.t:
            raise TypeError("To be caracterized as rectangular waveguide tha value of h must be bigger than t")

        self._h = new_h

    @t.setter
    def t(self, new_t: float) -> None:
        if new_t < 0:
            raise TypeError("Only values bigger than 0 are allowed")
        elif new_t > self.h:
            raise TypeError("To be caracterized as rectangular waveguide tha value of t must be smaller than h")

        self._t = new_t

    @hc.setter
    def hc(self, new_hc: float) -> None:
        if new_hc < 0:
            raise TypeError("Only values bigger than 0 are allowed")

        self._hc = new_hc

    @Waveguide.neff.setter
    def neff(self, new_neff: np.ndarray) -> None:
        if(type(new_neff) is not np.ndarray and not isinstance(new_neff, list)):
            raise TypeError("The effective index of a retangular waveguide must be an array of shape (2,)")
        if(isinstance(new_neff, list)):
            new_neff = np.ndarray(new_neff)
        if(new_neff.shape==(2,)):
            self._neff = new_neff
        else:
            raise AttributeError("The effective index of a retangular waveguide must be an array of shape (2,), \
                because one value is for the region -a<x<a and another for the regions x>a and x<-a")

    # parameters calculation
    # ref: eq. 2.89a-2.89e of Okamoto [1]
    def kappa(self, ne: float) -> np.complex128:
        return self.k0*np.sqrt(self.n1**2 - ne**2)

    def sigma(self, ne: float) -> np.complex128:
        return self.k0*np.sqrt(ne**2 - self.n3**2)

    def gamma(self, ne: float) -> np.complex128:
        return self.k0*np.sqrt(ne**2 - self.n2**2)

    def alpha(self, ne: float) -> np.complex128:
        return np.arctan((self.sigma(ne)/self.kappa(ne)))

    def psi(self, ne: float) -> np.complex128:
        return np.arctanh(self.sigma(ne)/self.gamma(ne))


    # dispersion equation of four-layer slab waveguide
    # ref: eq. 2.88 of Okamoto [1]
    def dispersionFunction(self, mode: float, x: float) -> float:
        return np.sin(self.kappa(x)*self.hc - 2*self.alpha(x)) - \
           np.sin(self.kappa(x)*self.hc)*np.exp(-2*(self.sigma(x)*self.rib_height+self.psi(x)))

    def calculateNeff(self, lbd: float) -> np.array:
        self.lbd = lbd
        self.k0 = 2*PI/self.lbd

        s = [self.h,self.t] # heights of rib inside and outside the core region

        n_e = []
        # run it for both heights
        for rib_index, rib_height in enumerate(s):
            self.rib_height = rib_height

            n_e.append(self._bisection.findEigenvalue())
        self.neff = np.array(n_e)
        return n_e


