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
import math

from opticalpyber.utils.bisection_method import Bisection
from opticalpyber.utils.constants import PI, C_0, MU_0, EPSILON_0

from opticalpyber.waveguide import Waveguide


class RetangularWaveguide(Waveguide):
    def __init__(self, n1: float, n2: float, n3: float, d: float, h: float, t: float, hc: float):
        super().__init__(n1, n2, d, n3=n3)

        # retangular waveguide dimensions
        self._h = h  # rib bigger height
        self._t = t  # rib smaller height, if =0, than is a ridge waveguide
        self._hc = hc  # core height

        self._rib_height = None

        # to solve dispersion equations
        self._bisection = Bisection(divsn=(self.n1 - self.n2) / 100, eps=1e-6, top_value=self.n1, bottom_value=self.n2)

        # solve dispersion equations
        self._bisection.callback = self.dispersionFunction

        # mode
        self._mode: int = None

    @property
    def h(self) -> float:
        return self._h

    @property
    def t(self) -> float:
        return self._t

    @property
    def hc(self) -> float:
        return self._hc

    @property
    def mode(self) -> int:
        return self._mode

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

    @mode.setter
    def mode(self, new_mode: int) -> None:
        if (type(new_mode) is not int):
            raise TypeError("Mode Type is always a integer, because is the subindex of M mode")
        self._mode = new_mode

    @Waveguide.neff.setter
    def neff(self, new_neff: np.ndarray) -> None:
        if (type(new_neff) is not np.ndarray and not isinstance(new_neff, list)):
            raise TypeError("The effective index of a retangular waveguide must be an array of shape (2,)")
        if (isinstance(new_neff, list)):
            new_neff = np.ndarray(new_neff)
        if (new_neff.shape == (2,)):
            self._neff = new_neff
        else:
            raise AttributeError("The effective index of a retangular waveguide must be an array of shape (2,), \
                because one value is for the region -a<x<a and another for the regions x>a and x<-a")

    # X field functions
    def X_1(self, u1: np.complex128, x: float) -> np.complex128:
        return np.cos(u1 / self.a * x - (self.mode * PI) / 2)

    def dX_1(self, u1: np.complex128, x: float) -> np.complex128:
        return -u1 / self.a * np.sin(u1 / self.a * x - (self.mode * PI) / 2)

    def dX2_1(self, u1: np.complex128, x: float) -> np.complex128:
        return -((u1 / self.a) ** 2) * np.cos(u1 / self.a * x - (self.mode * PI) / 2)

    def X_2(self, u1: np.complex128, w1: np.complex128, x: float) -> np.complex128:
        return np.cos(u1 + (self.mode * PI) / 2) * np.exp(w1 / self.a * (x + self.a))

    def dX_2(self, u1: np.complex128, w1: np.complex128, x: float) -> np.complex128:
        return w1 / self.a * np.cos(u1 + (self.mode * PI) / 2) * np.exp(w1 / self.a * (x + self.a))

    def dX2_2(self, u1: np.complex128, w1: np.complex128, x: float) -> np.complex128:
        return ((w1 / self.a) ** 2) * np.cos(u1 + (self.mode * PI) / 2) * np.exp(w1 / self.a * (x + self.a))

    def X_3(self, u1: np.complex128, w1: np.complex128, x: float) -> np.complex128:
        return np.cos(u1 - (self.mode * PI) / 2) * np.exp(-w1 / self.a * (x - self.a))

    def dX_3(self, u1: np.complex128, w1: np.complex128, x: float) -> np.complex128:
        return -w1 / self.a * np.cos(u1 - (self.mode * PI) / 2) * np.exp(-w1 / self.a * (x - self.a))

    def dX2_3(self, u1: np.complex128, w1: np.complex128, x: float) -> np.complex128:
        return ((w1 / self.a) ** 2) * np.cos(u1 - (self.mode * PI) / 2) * np.exp(-w1 / self.a * (x - self.a))

    def XfieldFunc(self, xx: float, u1: np.complex128, w1: np.complex128, regx: int) -> np.array:

        if regx == 1:
            X = self.X_1(u1, xx)
            dXdx = self.dX_1(u1, xx)
            dXdx2 = self.dX2_1(u1, xx)

        elif regx == 2:
            X = self.X_2(u1, w1, xx)
            dXdx = self.dX_2(u1, w1, xx)
            dXdx2 = self.dX2_2(u1, w1, xx)

        elif regx == 3:
            X = self.X_3(u1, w1, xx)
            dXdx = self.dX_3(u1, w1, xx)
            dXdx2 = self.dX2_3(u1, w1, xx)

        return [X, dXdx, dXdx2]

    # Y field functions
    def K1(self, ne: float) -> float:
        return np.sqrt(1 + (self.kappa(ne) / self.sigma(ne)) ** 2) / 2

    def Y_1(self, ne: float, y: float) -> float:
        return np.cos(self.alpha(ne)) * np.exp(self.sigma(ne) * y)

    def dY_1(self, ne: float, y: float) -> float:
        return self.sigma(ne) * np.cos(self.alpha(ne)) * np.exp(self.sigma(ne) * y)

    def Y_2(self, ne: float, y: float) -> float:
        return np.cos(self.kappa(ne) * y - self.alpha(ne))

    def dY_2(self, ne: float, y: float) -> float:
        return -self.kappa(ne) * np.sin(self.kappa(ne) * y - self.alpha(ne))

    def Y_3(self, ne: float, y: float) -> float:
        return self.K1(ne) * (np.sin(self.kappa(ne) * 2 * self.a) * np.exp(-self.sigma(ne) * (y - 2 * self.a)) - \
                              np.sin(self.kappa(ne) * 2 * self.a - 2 * self.alpha(ne)) * np.exp(
                    self.sigma(ne) * (y - 2 * self.a)))

    def dY_3(self, ne: float, y: float) -> float:
        return -self.K1(ne) * self.sigma(ne) * (
                    np.sin(self.kappa(ne) * 2 * self.a) * np.exp(-self.sigma(ne) * (y - 2 * self.a)) + \
                    np.sin(self.kappa(ne) * 2 * self.a - 2 * self.alpha(ne)) * \
                    np.exp(self.sigma(ne) * (y - 2 * self.a)))

    def Y_4(self, ne: float, y: float, s: float) -> float:
        return self.K1(ne) * (np.sin(self.kappa(ne) * 2 * self.a) * np.exp(-self.sigma(ne) * s) - np.sin(
            self.kappa(ne) * 2 * self.a - \
            2 * self.alpha(ne)) * np.exp(self.sigma(ne) * s)) * np.exp(-self.gamma(ne) * (y - 2 * self.a - s))

    def dY_4(self, ne: float, y: float, s: float) -> float:
        return -self.gamma(ne) * self.K1(ne) * (np.sin(self.kappa(ne) * 2 * self.a) * np.exp(-self.sigma(ne) * s) - \
                                                np.sin(self.kappa(ne) * 2 * self.a - 2 * self.alpha(ne)) * np.exp(
                    self.sigma(ne) * s)) * \
               np.exp(-self.gamma(ne) * (y - 2 * self.a - s))

    def YfieldFunc(self, yy: float, ne: float, s: float, regy: int) -> np.array:

        if regy == 1:
            Y = self.Y_1(ne, yy)
            dYdy = self.dY_1(ne, yy)
        elif regy == 2:
            Y = self.Y_2(ne, yy)
            dYdy = self.dY_2(ne, yy)
        elif regy == 3:
            Y = self.Y_3(ne, yy)
            dYdy = self.dY_3(ne, yy)
        elif regy == 4:
            Y = self.Y_4(ne, yy, s)
            dYdy = self.dY_4(ne, yy, s)

        return [Y, dYdy]

    # parameters calculation
    # ref: eq. 2.89a-2.89e of Okamoto [1]
    def kappa(self, ne: float) -> np.complex128:
        return self.k0 * np.sqrt(self.n1 ** 2 - ne ** 2)

    def sigma(self, ne: float) -> np.complex128:
        return self.k0 * np.sqrt(ne ** 2 - self.n2 ** 2)

    def gamma(self, ne: float) -> np.complex128:
        return self.k0 * np.sqrt(ne ** 2 - self.n3 ** 2)

    def alpha(self, ne: float) -> np.complex128:
        return np.arctan((self.sigma(ne) / self.kappa(ne)))

    def psi(self, ne: float) -> np.complex128:
        return np.arctanh(self.sigma(ne) / self.gamma(ne))

    # dispersion equation of four-layer slab waveguide
    # ref: eq. 2.88 of Okamoto [1]
    def dispersionFunction(self, mode: int, x: float) -> np.complex128:
        return np.sin(self.kappa(x) * self.hc - 2 * self.alpha(x)) - \
               np.sin(self.kappa(x) * self.hc) * np.exp(-2 * (self.sigma(x) * self.rib_height + self.psi(x)))

    # dispersion equation E^x_PQ
    # ref: section 2.2.4 of Okamoto [1]
    def dispersionX(self, m: int, b: float) -> np.complex128:
        return self.V * np.sqrt(1 - b) - m * PI / 2 - np.arctan(
            self.neff[0] ** 2 / self.neff[1] ** 2 * np.sqrt(b / (1 - b)))

    def calculateNeff(self, lbd: float) -> np.array:
        self.lbd = lbd
        self.k0 = 2 * PI / self.lbd

        s = [self.h, self.t]  # heights of rib inside and outside the core region

        n_e = []
        # run it for both heights
        for rib_index, rib_height in enumerate(s):
            self.rib_height = rib_height

            n_e.append(self._bisection.findEigenvalue())
            self._bisection = Bisection(divsn=(self.n1 - self.n2) / 100, eps=1e-6, top_value=self.n1, bottom_value=self.n2)
            self._bisection.callback = self.dispersionFunction

        self.neff = np.array(n_e)
        self.V = self.k0 * (2 * self.a) * np.sqrt(self.neff[0] ** 2 - self.neff[1] ** 2)
        return n_e

    def calculateDispersionFieldX(self, mode_number: int = 0, amplitude: float = 1.0, phi: float = 0,
                                  dx: float = 1E-8) -> tuple[np.array]:

        self.amplitude = amplitude
        self.phi = phi
        self.dx = dx
        self.mode = mode_number

        # updae bisection parameters
        self._bisection = Bisection(divsn=1E-3,eps=1e-2,mode=mode_number)

        self._bisection.callback = self.dispersionX

        b = self._bisection.findEigenvalue()

        # calculate beta througth unitary variables
        u = self.V * np.sqrt(1 - b)
        w = self.V * np.sqrt(b)

        self.beta = self.k0 * np.sqrt((w / self.k0 / self.a) ** 2 + self.neff[1] ** 2)

        return self.getDispersionField(u, w)

    def getDispersionField(self, u: float, w: float) -> tuple[np.array]:
        # create the output arrays
        XX = np.linspace(-10*self.a, 10*self.a, 50,dtype="complex_")
        YY = np.linspace(-3*self.a, 2*self.a+self.h+3e-6, 50,dtype="complex_")

        l = [len(XX), len(YY)]
        X = np.zeros(l[0],dtype="complex_")
        dXdx = X
        dXdx2 = X
        Y = np.zeros(l[1],dtype="complex_")
        dYdy = Y
        Zx = np.zeros([l[1], l[0]],dtype="complex_")

        for i in range(l[0]):
            if abs(XX[i]) <= self.a:
                n_ = self.neff[0]
                RegX = 1
                s = h
            elif XX[i] < -self.a:
                n_ = self.neff[1]
                RegX = 2
                s = t
            elif XX[i] > self.a:
                n_ = self.neff[1]
                RegX = 3
                s = t

            X_ = self.XfieldFunc(XX[i], u, w, RegX)
            X[i] = X_[0]
            dXdx[i] = X_[1]
            dXdx2[i] = X_[2]

            for k in range(l[1]):
                if YY[k] < 0:
                    RegY = 1
                    n = self.n2
                elif (YY[k] >= 0) and (YY[k] < 2*self.a):
                    RegY = 2
                    n = self.n1
                elif (YY[k] >= d) and (YY[k] < 2*self.a+s):
                    RegY = 3
                    n = self.n2
                else:
                    RegY = 4
                    n = self.n3

                K = (self.k0*C_0)*EPSILON_0*n**2*self.beta

                Y_ = self.YfieldFunc(YY[k],n_,s, RegY)

                Y[k] = Y_[0]
                dYdy[k] = Y_[1]
                Hy = X_[0]*Y_[0]
                dHydx2 = X_[2]*Y_[0]
                Ex = (self.k0*C_0)*MU_0/self.beta*Hy + dHydx2/K
                Zx[k,i] = Ex

        return (XX, YY, Zx)
