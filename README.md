# Optical Pyber

Package with tools to analize waveguides and optical fibers. This package uses numerical methods to find the eigenvalues related to the propagation modes and other characteristics of the waveguide. Currently, this package supports three types of waveguide: Planar, Retangular and Optical Fibers.

## Requirements
- numpy
- scipy
- matplotlib (if you want to plot something)


## Examples

### Planar Waveguide

(todo)

### Retangular Waveguide

(todo)

### Optical Fiber

- Find and plot the modes LP01, LP11 and LP21 of an Step-index fiber[4].

```
""" Find and plot LP Modes of an Optical Fiber """
from opticalpyber import OpticalFiber
import matplotlib.pyplot as plt

# materials
n1 = 1.45 # core
n2 = 1.4  # cladding

# diameter of the core
d = 16E-6 # um

# wavelength of analysis
lbd = 8.11E-6 # um

fiber = OpticalFiber(n1, n2, d)

for mode in range(3):
    modes = fiber.findMModes(lbd, mode)
    print(f'For mode M = {mode} we found {len(modes)} values.')
    print(modes)

    # plot the modes
    (x, y, z) = a.calculateDispersionField(modes[-1], dx=1e-8)

    figure(i)
    plt.contourf(x, y, z)
    plt.colorbar()
    plt.show()
```

## References

[1] Okamoto, K. (2021). Fundamentals of optical waveguides. Elsevier.

[2] GHATAK, A. A., Ghatak, A., Thyagarajan, K., & Thyagarajan, K. (1998). An introduction to fiber optics. Cambridge university press.

[3] Orfanidis, S. J. (2002). Electromagnetic waves and antennas.

[4] Paschotta, D., 2022. LP Modes. [online] Rp-photonics.com. Available at: <https://www.rp-photonics.com/lp_modes.html> [Accessed 15 April 2022].
