<h1 align="center">
    <img alt="OpticalPyber" ttle="OpticalPyber" src="https://github.com/marcos-moura97/opticalpyber/blob/main/docs/logo.png" width="30%" height="auto"/>
</h1>

<p align="center">
  <a href="#about-the-project">About the Project</a>&nbsp;&nbsp;&nbsp;|&nbsp;&nbsp;&nbsp;
  <a href="#requirements">Requirements</a>&nbsp;&nbsp;&nbsp;|&nbsp;&nbsp;&nbsp;
  <a href="#examples">Examples</a>&nbsp;&nbsp;&nbsp;|&nbsp;&nbsp;&nbsp;
  <a href="#next-steps">Next Steps</a>&nbsp;&nbsp;&nbsp;|&nbsp;&nbsp;&nbsp;
  <a href="#references">References</a>
</p>

<br />

## About the Project

Package with tools to analize waveguides and optical fibers. This package uses numerical methods to find the eigenvalues related to the propagation modes and other characteristics of the waveguide. Currently, this package supports three types of waveguide: Planar, Retangular and Optical Fibers.

This package is made based on python codes used during my graduation, this codes can be found in [this repository](https://github.com/marcos-moura97/eletromagnetism_python.git).

## Requirements
- numpy
- scipy
- matplotlib (if you want to plot something)

## Examples

### Planar Waveguide

- Find the effective index and plot the profile dispersion curve for modes TM0, TM1 and TM2 in an symetric guide [1-2].
```
from opticalpyber import PlanarWaveguide
import matplotlib.pyplot as plt

# materials
n1 = 3.38       #core
n2 = 3          #cladding
n3 = n2         #substrate

d = 2.5E-6      #core length
lbd = 1E-6      # wavelength

guide = PlanarWaveguide(n1, n2, n3, d)

for mode in range(3):
    Ey = guide.calculateDispersionField(lbd, mode="TM", mode_number=mode)

    # make the x arrays
    x = np.arange(-2*guide.lbd,2*guide.lbd,guide.dx)

    # ploting
    plt.figure(mode)
    plt.plot(x,Ey[:-1],'b')
    # plot the bars that divide the materials
    plt.plot([-guide.a,-guide.a],[-guide.amplitude,guide.amplitude],'k--') # substracte / core
    plt.plot([guide.a,guide.a],[-guide.amplitude,guide.amplitude],'k--')   # core / cladding
    # some info
    plt.title(f'Dispersion of Ey for mode TM{mode}')
    plt.xlabel('y')
    plt.ylabel('Ey')
    plt.xticks([-guide.a,guide.a], ['-a','a'])
    plt.show()

```

### Retangular Waveguide

- Find the effective index and plot the profile dispersion curve for modes TM0 in an symetric rib waveguide[1].

```
from opticalpyber import RetangularWaveguide
import matplotlib.pyplot as plt

## materials

n1 = 3.50125  # core
n2 = 3.5  # substrate/rib
n3 = 3  # cladding

d = 2 * 8e-6  # core length

h = 3e-6  # rib bigger height
t = 1e-6  # rib smaller height, if =0, than is a ridge waveguide
hc = 6.5e-6  # core height

lbd = 1.32e-6  # wavelength

guide = RetangularWaveguide(n1, n2, n3, d, h, t, hc)
all_modes = guide.calculateNeff(lbd)
print(guide.neff)

(XX, YY, Zx) = guide.calculateDispersionFieldX()

# Ploting
plt.figure()
# plot
plt.contourf(XX, YY, Zx)
# structure plot
plt.plot([-(3*guide.a), -guide.a, -guide.a, guide.a, guide.a, 3*guide.a],
         [guide.hc+guide.t, guide.hc+guide.t, guide.hc+guide.h, guide.hc+guide.h, guide.hc+guide.t, guide.hc+guide.t], 'k')
plt.plot([-(3*guide.a), 3*guide.a], [guide.hc, guide.hc], 'k')
plt.plot([-(3*guide.a), 3*guide.a], [0, 0], 'k')
# labels
plt.xticks([-guide.a,guide.a], ['-a','a'])
plt.yticks([0,guide.hc,guide.hc+guide.t, guide.hc+guide.h], ['0','hc','hc+t','hc+h'])
# limits
plt.xlim([-(3*guide.a), 3*guide.a])
plt.ylim([-guide.t, guide.hc+guide.h+guide.t])

plt.show()

```

### Optical Fiber

- Find and plot the modes LP01, LP11 and LP21 of an Step-index fiber[4].

```
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
    (x, y, z) = fiber.calculateDispersionField(modes[-1], dx=1e-8)

    plt.figure(mode)
    plt.contourf(x, y, z)
    plt.colorbar()
    plt.show()
```

## Next Steps

- Unit Tests
- Add Slot Waveguides
- More Robust Optical Fiber Model

## References

[1] Okamoto, K. (2021). Fundamentals of optical waveguides. Elsevier.

[2] GHATAK, A. A., Ghatak, A., Thyagarajan, K., & Thyagarajan, K. (1998). An introduction to fiber optics. Cambridge university press.

[3] Orfanidis, S. J. (2002). Electromagnetic waves and antennas.

[4] Paschotta, D., 2022. LP Modes. [online] Rp-photonics.com. Available at: <https://www.rp-photonics.com/lp_modes.html> [Accessed 15 April 2022].
