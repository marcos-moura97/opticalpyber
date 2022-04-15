"""

Basic properties of a waveguide

"""
from typing import Optional


class Waveguide:
    def __init__(self, n1: float, n2: float, d: float, n3: Optional[float]=None):
        self._n1 = n1       # core
        self._n2 = n2       # substrate
        self._n3 = n3       # cladding

        self._a = d/2       # lenght of waveguide

        # initialize other properties
        self._lbd: float       = None # wavelength of simulation
        self._k0: float        = None # wavenumber
        self._V: float         = None # normalized frequency

        self._neff: float      = None # effective refractive index
        self._beta: float      = None
        self._gamma: float     = (n2**2-n3**2)/(n1**2 - n2**2) if n3 else n2**2/(n1**2 - n2**2)

        self._dx: float        = None # discretized X axis

    @property
    def n1(self) -> float:
        return self._n1

    @property
    def n2(self) -> float:
        return self._n2

    @property
    def n3(self) -> float:
        return self._n3

    @property
    def neff(self) -> float:
        return self._neff

    @property
    def a(self) -> float:
        return self._a

    @property
    def lbd(self) -> float:
        return self._lbd

    @property
    def k0(self) -> float:
        return self._k0

    @property
    def V(self) -> float:
        return self._V

    @property
    def dx(self) -> float:
        return self._dx


    # setters
    @n1.setter
    def n1(self, new_n1: float) -> None:
        if new_n1 < self.n2 or new_n1 < self.n3:
            raise TypeError("Core refractive index must have the biggest refractive index")
        elif new_n1 < 0:
            raise TypeError("Only values bigger or equal 0 are allowed")
        else:
            self._n1 = new_n1

    @n2.setter
    def n2(self, new_n2: float) -> None:
        if new_n2 > 0:
            self._n2 = new_n2
        else:
            raise TypeError("Only values bigger or equal 0 are allowed")

    @n3.setter
    def n3(self, new_n3: float) -> None:
        if new_n3 > 0:
            self._n3 = new_n3
        else:
            raise TypeError("Only values bigger or equal 0 are allowed")

    @neff.setter
    def neff(self, new_neff: float) -> None:
        if new_neff > 0:
            self._neff = new_neff
        else:
            raise TypeError("Only values bigger or equal 0 are allowed")

    @a.setter
    def a(self, new_d: float) -> None:
        if new_d > 0:
            self._a = new_d/2
        else:
            raise TypeError("Dimensions can not be negatives")

    @k0.setter
    def k0(self, new_k0: float) -> None:
        self._k0 = new_k0

    @V.setter
    def V(self, new_V: float) -> None:
        self._V = new_V

    @lbd.setter
    def lbd(self, new_lbd: float) -> None:
        if new_lbd > 0:
            self._lbd = new_lbd
        else:
            raise TypeError("Dimensions can not be negatives")

    @dx.setter
    def dx(self, new_dx: float) -> None:
        if new_dx <= self.lbd/1E2:
            self._dx = new_dx
        else:
            raise TypeError("The minimum step must be at least 100 smaller than the wavelength")
