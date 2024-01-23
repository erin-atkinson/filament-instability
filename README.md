# filament-instability
Detangling the interactions between turbulence and submesoscale flow

Simulations of an oceanic mixed layer filament.
Inputs and default (nominal) values
Ro=1: Maximum magnitude of Rossby number
Ri=0.6: Minimum Richardson number
Ek=nothing: (Turbulent) Ekman number, for closure. If nothing then a Smagorinsky-Lilly default is used
Pr=1: Prandtl number for buoyancy diffusion
α=0.25: Front width - separation ratio
λ=0.05: Fractional width of transition to deep water
δ=-0.25: Fractional height change of transition to deep water across filament

The boundary layer height is defined by a shape function
$$\gamma(s; \alpha, \delta) = -1 + \frac{\delta}{2} \left(\text{erf}\left(s + \frac{1}{2\alpha}\right) - \text{erf}\left(s - \frac{1}{2\alpha}\right)\right)$$

The simulation initially develops some boundary layer turbulence. The fluctuation fields $\vec {u}'$ and $b'$ are then added to the nominal filament state.
