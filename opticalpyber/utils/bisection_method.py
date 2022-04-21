"""

Bisection Method to Calculate EigenValue

"""
from typing import Callable

import numpy as np


class Bisection:
    def __init__(self,divsn: float = .01, eps: float = 1E-6, mode: int = 0, top_value: float = 1, bottom_value: float = 0):
        self.divsn = divsn
        self.eps = eps

        self.b0: np.complex128 = np.complex128(top_value-self.eps) # initial guess
        self.b1: np.complex128 = self.b0
        self.b2: np.complex128 = self.b1-self.divsn
        self.bottom_value = bottom_value

        self.mode = mode
        self._callback: Callable = None


    @property
    def callback(self):
        return self._callback

    @callback.setter
    def callback(self,new_callback: Callable) -> None:
        self._callback = new_callback

    def findEigenvalue(self) -> float:
        F1 = self.callback(self.mode, self.b1)
        while self.divsn>self.eps:
            F2 = self.callback(self.mode, self.b2)
            if F1*F2<0:
                self.b2 = (self.b1+self.b2)/2
                self.divsn = self.divsn/2
            else:
                self.b1=self.b2
                self.b2 = self.b1-self.divsn
                F1=F2
            if self.b2<self.bottom_value and self.b1<self.bottom_value:
                # raise warning
                return self.bottom_value
        b = (self.b1+self.b2)/2
        if b<0:
            b=0
        return b
