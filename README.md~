# Contributions to Jupiter’s gravity field from dynamics in the dynamo region

## Requirements 
- `Python 3+` 
- `NumPy`
- 'SciPy'

## 1 Overview 

<p align="center">
  <img src="figures/jupiter_flows.jpg" width="1000">
    <br>
 <em> <font size = "4"> Zonal flows occur in Jupiter's deep atmosphere and dynamo region (i.e., where the magnetic
 field is <br> generated). These flows are associated with density perturbations that alter the planet's gravity field.</font> </em>  
</p>


Zonal (east-west) flows occur throughout Jupiter's interior. Near the "surface" of the planet, we observe alternating east-west jets with peak speeds of ~100 m/s. The surface winds may extend thousands of kilometers deep into the planet, or be confined to a thin weather layer on top of weaker convective flows. Deeper down, in the dynamo region where the magnetic field is generated, there are slow zonal flows with speeds on the order of 0.1-20 cm/s. The zonal flows within the deep atmosphere and dynamo region are associated with high and low pressures, and thus with density perturbations relative to the mean background density. The density perturbations induced by zonal flows alter Jupiter's gravity field.

Zonal flows in the upper 3000 km of Jupiter's non-conducting atmosphere with amplitudes greater than a few meters per second can produce gravity signals that are consistent with the recent Juno gravity observations. The gravity contribution from flows deeper down in the dynamo region, however, have not yet been investigated. Zonal flow in the dynamo region could have a non-negligible contribution to Jupiter's gravity field. Although the dynamo zonal flow is expected to be only on the order of 0.1-20 cm/s the flow occupies most of the volume of the planet and resides in a region with a large background density. The slow motion of a large volume of dense fluid has the potential to produce a large gravity signal. 

In this project, we calculate the gravity signal produced physically motivated dynamo zonal flow profiles. We begin by constructing zonal flow profiles for the dynamo region in [Section 2](#2-determining-the-zonal-flow-in-the-dynamo-region). We then develop the mathematical framework to calculate the density perturbation and gravity signal associated with the zonal flow in [Section 3](#3-gravity-calculation). Then, we calculate the gravity signal for a few example flows in [Section 4](#4-Results). A more detailed discussion of the methods and results be found in Kulowski et al. (2020) (in review). 
  
## 2 Determining the zonal flow in the dynamo region 

The first step in this project is to determine the zonal flow in the dynamo region. We simplify the dyamics so that the flow obeys Ferraro's law of isorotation (i.e., fluid parcels connected by the same magnetic field line rotate at the same angular rate). Given the internal magnetic field and the angular velocity at the dynamo surface, we can specify the zonal flow throughout the dynamo region by extending the angular velocity at the dynamo surface into the interior along the magnetic field lines. 

A prior, we do not know the structure of Jupiter's internal magnetic field or the angular velocity at the dynamo surface. We therefore assume that they have simple forms. We consider two different models for the internal magnetic field: a simple polynomial model (i.e., the spectral poloidal scalar is has a polynomial form) and a minimal Ohmic dissipation model (i.e., where the magnetic field has a configuration that minimizes Ohmic dissipation). At the dynamo surface, we assume that the zonal flow is proportional to the zonally averaged, theta-component of Jupiter's magnetic field. We allow for some small adjustments to this flow so that we can satisfy Ferraro's law. We extend the angular velocity associated with the zonal flow at the dynamo surface into the interior along magnetic field lines to the dynamo zonal flow profile. The constructed zonal flow profiles for each magnetic field model are shown below. 
<p align="center">
  <img src="figures/ferraro.jpg" width="800">
    <br>
</p>

## 3 Gravity calculation

There are two steps to calculate the gravity signal produced by a zonal flow profile. First, we calculate the density perturbation associated with the zonal flow using the vorticity equation, which is given by 
<p align="center">
<img src="https://latex.codecogs.com/gif.latex?2(\boldsymbol{\Omega}&space;\cdot&space;\nabla)(\rho_{0}&space;\boldsymbol{u})&space;=&space;-&space;\nabla&space;\rho&space;'&space;\times&space;\boldsymbol{g_{\rm{eff}_{0}}}" title="2(\boldsymbol{\Omega} \cdot \nabla)(\rho_{0} \boldsymbol{u}) = - \nabla \rho ' \times \boldsymbol{g_{\rm{eff}_{0}}}" />
</p>

where <img src="https://latex.codecogs.com/gif.latex?\boldsymbol{\Omega}" title="\boldsymbol{\Omega}" /> is the planet's rotational angular velocity, <img src="https://latex.codecogs.com/gif.latex?\rho_{0}" title="\rho_{0}" /> is the hydrostatic background density, <img src="https://latex.codecogs.com/gif.latex?\rho'" title="\rho'" /> is the density perturbation arising from the flow, and <img src="https://latex.codecogs.com/gif.latex?\boldsymbol{g_{\rm{eff}_{0}}}" title="\boldsymbol{g_{\rm{eff}_{0}}}" /> is the backgound effective gravity. Given a zonal flow profile, we integrate vorticity equation to obtain the density perturbation.

Having obtained the density perturbation, we can now calculate the gravity signal associated with it. We compute the gravity signal in spectral space, so that the gravity field is represented by zonal gravity harmonics. The zonal gravity harmonics are given by 

<p align="center">
<img src="https://latex.codecogs.com/gif.latex?\Delta&space;J_{n}&space;=&space;-&space;\frac{1}{Ma^{n}}&space;\int_{V}\rho'r^{n}P_{n}(\cos\theta)dV" title="\Delta J_{n} = - \frac{1}{Ma^{n}} \int_{V}\rho'r^{n}P_{n}(\cos\theta)dV" />
</p>

where <img src="https://latex.codecogs.com/gif.latex?n" title="n" /> is the spherical harmonic degree, <img src="https://latex.codecogs.com/gif.latex?M" title="M" /> is Jupiter's mass, <img src="https://latex.codecogs.com/gif.latex?a" title="a" /> is the equatorial radius, <img src="https://latex.codecogs.com/gif.latex?r" title="r" /> is the radial position, <img src="https://latex.codecogs.com/gif.latex?P_{n}(\cos\theta)" title="P_{n}(\cos\theta)" /> is the Legendre polynomial of degree <img src="https://latex.codecogs.com/gif.latex?n" title="n" />, and <img src="https://latex.codecogs.com/gif.latex?V" title="V" /> is Jupiter's volume.

## 4 Results 

We calculate the zonal gravity harmonics (J2-J10) for our physically motivated dynamo zonal flow profiles. The zonal flow profiles with RMS velocities of 10 cm/s for the simple polynomial and minimal Ohmic dissipation magnetic field models are shown below. For each model, we allow Jupiter's core to be compact or dilute.

<p align="center">
  <img src="figures/jns.png" width="900">
</p>

These results can be produced by the files `simple_polynomial.py` and `minimal_ohmic_dissipation.py`. For all cases, the dynamo zonal flow would produce J3 values on the same order of magnitude as the Juno inferred value and J2 and J4 values on the same order as 3000 km deep atmospheric zonal flow, but would not contribute much to higher order gravity harmonics.
