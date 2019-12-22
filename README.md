# Least square IIR filter design
## Design IIR filter with least square based method

Compute L2 solution by iteratively solving overdetermined linear equations.

User must specify desired complex frequency response and inital weights to weight the error on frequency grid.

# Examples:
## Passband linear phase low pass IIR filter
![Diagram1](./graph/LPF.svg)

## Passband linear phase band pass IIR filter
![Diagram1](./graph/BPF.svg)

## Passband linear phase arbitrary bands IIR filter
![Diagram1](./graph/BPF_HPF.svg)

# Discussion

## FDLS vs Eqnerror
All 3 methods design digital filter base on arbitrary frequency grid and gain vector as input.

From filter designer perspective, this property is very desirable.

Designer once have to craft frequency response equations on s-plane and convert them to z-plane using bilinear transform.

Using these 3 methods, designers can convert their analog frequency requirements to digital IIR filter directly.

Frequency domain least square(FDLS) method did pretty good job at preserving high frequency of the analog filter.

In my opinion, FDLS is easy version of Matched Z-transform, they both has potentials on designing filter that need to preserve frequency response around Nyquist. Although they work in completely different way.

However, FDLS is shouldn't be used if you need linear phase passband requirement, FDLS simply cannot handle sudden changed peaks in desire response.
