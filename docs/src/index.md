

```@meta
Author = "Andrea Zunino"
Author = "Alessandro Ghirotto"
```

# Contents

```@contents
Pages = ["index.md"]
Depth = 4
```

# User guide

**MagGrav2Dpoly** is a Julia package to perform magnetic and gravity anomaly calculations using a 2D or 2.75D parameterization in terms of polygons with uniform arbitrary magnetizations and density contrasts. It provides functions to both solve the forward problem and to calculate the gradient of a given misfit function. Such functions can be used to solve inverse problems both in the deterministic and probabilistic approach. In particular, this package provides some functions to solve inverse problems using the Hamiltonian Monte Carlo (HMC) method, as part of the `HMCLab` project (see the  `HMCSampler.jl` package). Gradients are calculated using the technique of automatic differentiation.
With this package it is also possible to perform joint magnetic and gravity forward and gradient calculations and hence solve joint inverse problems, see the tutorials below.

The forward problem formulations for the magnetic case implemented in this package are the following:
* 2D case: Talwani & Heirtzler (1962, 1964), Won & Bevis (1987) and revised Kravchinsky, Hnatyshin, Lysak, & Alemie (2019);
* 2.75D case: revised Rasmussen & Pedersen (1979) and Campbell (1983).

The forward problem formulations for the gravity case implemented in this package are the following:
* 2D case: Talwani, Worzel, & Landisman (1959), with background theory derived from the paper of Hubbert (1948);
* 2.75D case:  Rasmussen & Pedersen (1979).

  
If you use this code for research or else, please cite the related papers:

* Ghirotto, Zunino, Armadillo &, Mosegaard (2021). **Magnetic Anomalies Caused by 2D Polygonal Structures with Uniform Arbitrary Polarization: new insights from analytical/numerical comparison among available algorithm formulations**. *Geophysical Research Letters, 48*(7), e2020GL091732, https://doi.org/10.1029/2020GL091732.

* Zunino, Ghirotto, Armadillo, & Fichtner (2022). **Hamiltonian Monte Carlo probabilistic joint inversion of 2D (2.75D) gravity and magnetic data**. *Geophysical Research Letters*,  49, e2022GL099789. https://doi.org/10.1029/2022GL099789.

Regarding solving the inverse problem with the HMC method, please see the following paper and check out the package `HMCSampler.jl`:

* Zunino, Gebraad, Ghirotto, & Fichtner (2023). **HMCLab: a framework for solving diverse geophysical inverse problems using the Hamiltonian Monte Carlo method**. in prep.

In addition, a tutorial about the use of forward formulations and the basic tuning strategies for HMC inversion is presented in detail in the below section.


## Installation

To install the package first enter into the package manager mode in Julia by typing "`]`" at the 
REPL prompt and add the "JuliaGeoph" registry as
```
(@v1.8) pkg> registry add https://gitlab.com/JuliaGeoph/JuliaGeophRegistry
```
Then add the package by simply issuing
```
(@v1.8) pkg> add MagGrav2Dpoly
```
The package will be automatically downloaded from the web and installed.


## Theoretical Background

### Forward calculation

For a theoretical explanation, let us consider a three-dimensional non-magnetic and zero-density
space in which a body infinitely extended in the ``y`` direction is immersed. 

The common aim of all formulations is the calculation of the total-field magnetic intensity response
and vertical attraction of this body upon an observation point ``(x_0,z_0)`` located along a profile
aligned to the ``x`` direction (the positive ``z`` axis is assumed pointing downward).

The starting assumption is that our body can be considered as discretized by an 
infinite number of elementary volumes with uniform magnetizazion and density contrast
and infinitesimal dimensions ``dx``, ``dy``, ``dz``.

Within this assumption, the magnetic and gravity fields associated to the body can be mathematically 
expressed in terms of a line integral around its periphery, represented in two dimensions 
as its polygonal cross-section (in red).
![](images/intro.svg)


### Inverse calculations

See the related papers and examples for inverse calculations using the HMC strategy.


## Tutorial for magnetic calculations

First load the module,
```@example ex1
using Mag2Dpoly
nothing # hide
```
then define an array containing the location of the observation points, where the first column represents the ``x`` position and the second the  ``z`` position. Remark: ``z`` points *downward*! So the observations have a negative ``z`` in this case.
```@example ex1
# number of observations
N=101
xzobs = [LinRange(0.0,120.0,N) -1.0*ones(N)]
nothing # hide
```
In order to describe the polygonal bodies, two objects need to be specifies: 1)  an array containing all the positions of the vertices (first column represents the ``x`` position and the second the  ``z`` position) and 2) a mapping relating each polygonal body to its vertices. The position of vertices *must* be specifyied in a **counterclockwise** order.
```@example ex1
# vertices of the poligonal bodies
vertices  = [35.0 50.0;
             65.0 50.0;
             80.0 35.0;
             65.0 20.0;
             35.0 20.0;
             20.0 35.0;
             90.0 60.0;
             95.0 40.0]
			 
# indices of vertices for the first polygon
ind1 = collect(1:6)
# indices of vertices for the second polygon
ind2 = [2,7,8,3]
# define the two bodies in term of indices
bodyindices = [ind1,ind2]
nothing # hide
```
Now, we specify the magnetic properties for each of the polygonal bodies and the angle of the reference system with the north axis:
```@example ex1
# induced magnetization
Jind = MagnetizVector(mod=[4.9,3.5],Ideg=[0.0,0.0],Ddeg=[5.0,5.0])
# remanent magnetization
Jrem = MagnetizVector(mod=[3.1,2.5],Ideg=[45.0,30.0],Ddeg=[0.0,10.0])

# angle with the North axis
northxax = 90.0
nothing # hide
```
Finally, construct the  poligonal body object by instantiating a `MagPolygBodies2D` structure. Here we can determine the type of forward calculation, i.e., 2D, 2.5D or 2.75D by specifying the variable `ylatext`. There are three cases:
1. if `ylatext=nothing` then the polygonal bodies are considered to be extending laterally to infinity, hence a pure 2D forward calculation
2. if `ylatext` is a single real number, then the forward computation is 2.5D, i.e., the polygonal bodies extend laterally on both sides by an amount specified by the value of `ylatext`
3. if `ylatext` is a two-element vector, then the forward computation is 2.75D, i.e., the polygonal bodies extend laterally from `ylatext[1]` to `ylatex[2]` (with the condition `ylatex[2]>ylatex[1]`). 
Here we choose to run a 2.75D forward computation:
```@example ex1
pbody = MagPolygBodies2D(bodyindices,vertices,Jind,Jrem,ylatext=[50.0,90.0])
nothing # hide
```
At this point the magnetic field can be computed:
```@example ex1
tmag = tmagpolybodies2D(xzobs,northxax,pbody)
```
The output vector is the magnetic anomaly at each of the observation points specified above.

Now we can plot the results (e.g., using `CairoMakie`):
```
using CairoMakie

fig = Figure()

ax1 = Axis(fig[1,1],title="Magity anomaly",xlabel="x",ylabel="mag. anom.")
scatter!(ax1,xzobs[:,1],tmag)

ax2 = Axis(fig[2,1],title="Polygonal bodies",xlabel="x",ylabel="z")
x1 = [pbody.geom.bo[1].ver1[:,1]...,pbody.geom.bo[1].ver2[end,1]]
y1 = [pbody.geom.bo[1].ver1[:,2]...,pbody.geom.bo[1].ver2[end,2]]
scatterlines!(ax2,x1,y1)
x2 = [pbody.geom.bo[2].ver1[:,1]...,pbody.geom.bo[2].ver2[end,1]]
y2 = [pbody.geom.bo[2].ver1[:,2]...,pbody.geom.bo[2].ver2[end,2]]
scatterlines!(ax2,x2,y2)
ax2.yreversed=true

linkxaxes!(ax1,ax2)
display(fig)

save("mag.svg",fig) # hide
```

![](images/mag.svg)


## Tutorial for gravity calculations

First load the module,
```@example ex1
using Grav2Dpoly
nothing # hide
```
then define an array containing the location of the observation points, where the first column represents the ``x`` position and the second the  ``z`` position. Remark: ``z`` points *downward*! So the observations have a negative ``z`` in this case.
```@example ex1
# number of observations
N=101
xzobs = [LinRange(0.0,120.0,N) -1.0*ones(N)]
nothing # hide
```
In order to describe the polygonal bodies, two objects need to be specifies: 1)  an array containing all the positions of the vertices (first column represents the ``x`` position and the second the  ``z`` position) and 2) a mapping relating each polygonal body to its vertices. The position of vertices *must* be specifyied in a **counterclockwise** order.
```@example ex1
# vertices of the poligonal bodies
vertices  = [35.0 50.0;
             65.0 50.0;
             80.0 35.0;
             65.0 20.0;
             35.0 20.0;
             20.0 35.0;
             90.0 60.0;
             95.0 40.0]
			 
# indices of vertices for the first polygon
ind1 = collect(1:6)
# indices of vertices for the second polygon
ind2 = [2,7,8,3]
# define the two bodies in term of indices
bodyindices = [ind1,ind2]
nothing # hide
```
Now, we specify the density for each of the polygonal bodies
```@example ex1
# two bodies in this case
rho = [2000.0,3000.0]
nothing # hide
```
Finally, construct the  poligonal body object by instantiating a `GravPolygBodies2D` structure. Here we can determine the type of forward calculation, i.e., 2D, 2.5D or 2.75D by specifying the variable `ylatext`. There are three cases:
1. if `ylatext=nothing` then the polygonal bodies are considered to be extending laterally to infinity, hence a pure 2D forward calculation
2. if `ylatext` is a single real number, then the forward computation is 2.5D, i.e., the polygonal bodies extend laterally on both sides by an amount specified by the value of `ylatext`
3. if `ylatext` is a two-element vector, then the forward computation is 2.75D, i.e., the polygonal bodies extend laterally from `ylatext[1]` to `ylatex[2]` (with the condition `ylatex[2]>ylatex[1]`). 
Here we choose to run a 2.75D forward computation:
```@example ex1
pbody = GravPolygBodies2D(bodyindices,vertices,rho,ylatext=[50.0,90.0])
nothing # hide
```

At this point the gravity field can be computed:
```@example ex1
tgrav = tgravpolybodies2D(xzobs,pbody)
```
The output vector is the gravity anomaly at each of the observation points specified above.

Now we can plot the results (e.g., using `CairoMakie`):
```
using CairoMakie

fig = Figure()

ax1 = Axis(fig[1,1],title="Gravity anomaly",xlabel="x",ylabel="grav. anom.")
scatter!(ax1,xzobs[:,1],tgrav)

ax2 = Axis(fig[2,1],title="Polygonal bodies",xlabel="x",ylabel="z")
x1 = [pbody.geom.bo[1].ver1[:,1]...,pbody.geom.bo[1].ver2[end,1]]
y1 = [pbody.geom.bo[1].ver1[:,2]...,pbody.geom.bo[1].ver2[end,2]]
scatterlines!(ax2,x1,y1)
x2 = [pbody.geom.bo[2].ver1[:,1]...,pbody.geom.bo[2].ver2[end,1]]
y2 = [pbody.geom.bo[2].ver1[:,2]...,pbody.geom.bo[2].ver2[end,2]]
scatterlines!(ax2,x2,y2)
ax2.yreversed=true

linkxaxes!(ax1,ax2)
display(fig)

save("grav.svg",fig) # hide
```

![](images/grav.svg)


## Tutorial for joint mag and grav
First load the modules,
```@example ex1
using MagGrav2Dpoly
```
then define 1) the angle between the Magnetic Field's North and the model profile, 2) an array containing the location of the observation points along it, where the first column represents the ``x`` position and the second the  ``z`` position. Remark: ``z`` points *downward*! So the observations have a negative ``z`` in this case.
```@example ex1
# angle with the North axis
northxax = 90.0

# number of observations
N=101
xzobs_mag = [LinRange(0.0,120.0,N) -1.0*ones(N)]
xzobs_grav = copy(xzobs_mag)
nothing # hide
```
In order to describe the polygonal bodies, two objects need to be specifies: 1)  an array containing all the positions of the vertices (first column represents the ``x`` position and the second the  ``z`` position) and 2) a mapping relating each polygonal body to its vertices. The position of vertices *must* be specifyied in a **counterclockwise** order.
```@example ex1
# vertices of the poligonal bodies
vertices  = [35.0 50.0;
             65.0 50.0;
             80.0 35.0;
             65.0 20.0;
             35.0 20.0;
             20.0 35.0;
             90.0 60.0;
             95.0 40.0]
			 
# indices of vertices for the first polygon
ind1 = collect(1:6)
# indices of vertices for the second polygon
ind2 = [2,7,8,3]
# define the two bodies in term of indices
bodyindices = [ind1,ind2]
nothing # hide
```
Now, we specify the density and magnetization for each of the polygonal bodies
```@example ex1
# two bodies in this case
# densities
rho = [2000.0,3000.0]
# induced magnetizations
Jind = Mag2Dpoly.MagnetizVector(mod=[4.9,1.0],Ideg=[0.0,0.0],Ddeg=[5.0,5.0])
# remanent magnetizations
Jrem = Mag2Dpoly.MagnetizVector(mod=[3.1,1.5],Ideg=[45.0,-45.0],Ddeg=[0.0,0.0])
nothing # hide
```
Finally, construct the poligonal body object by instantiating a `JointPolygBodies2D` structure. Here we can determine the type of forward calculation, i.e., 2D, 2.5D or 2.75D by specifying the variable `ylatext`. There are three cases:
1. if `ylatext=nothing` then the polygonal bodies are considered to be extending laterally to infinity, hence a pure 2D forward calculation
2. if `ylatext` is a single real number, then the forward computation is 2.5D, i.e., the polygonal bodies extend laterally on both sides by an amount specified by the value of `ylatext`
3. if `ylatext` is a two-element vector, then the forward computation is 2.75D, i.e., the polygonal bodies extend laterally from `ylatext[1]` to `ylatex[2]` (with the condition `ylatex[2]>ylatex[1]`). 
Here we choose to run a 2.75D forward computation:
```@example ex1
pbody = JointPolygBodies2D(bodyindices,vertices,Jind,Jrem,rho,ylatext=[50.0,90.0])
nothing # hide
```

At this point the gravity and magnetic fields can be computed:
```@example ex1
# compute the gravity anomaly and total field magnetic anomaly
tgrav,tmag = tjointpolybodies2D(xzobs_grav,xzobs_mag,northxax,pbody)
nothing # hide
```
The output vectors are the gravity and magnetic anomalies at each of the observation points specified above.

Alternatively, in the 2D case we could choose among other forward formulations implemented, specified as strings. Now we choose as forward types for the gravity and magnetic case `"wonbev"` and `"talwani_red"` respectively (see docs of `Grav2Dpoly` and `Mag2Dpoly` for details):
```@example ex1
pbody = JointPolygBodies2D(bodyindices,vertices,Jind,Jrem,rho,ylatext=nothing)
nothing # hide
# type of forward algorithms of the gravity and magnetic case
forwtype_grav = "wonbev"
forwtype_mag = "talwani_red"
# compute the gravity anomaly and total field magnetic anomaly 
tgrav,tmag = tjointpolybodies2Dgen(xzobs_grav,xzobs_mag,northxax,pbody,forwtype_grav,forwtype_mag)
nothing # hide
```

Now we can plot the results (e.g., using `CairoMakie`):
```
using CairoMakie

fig = Figure()

ax1 = Axis(fig[1,1],title="Gravity anomaly",xlabel="x",ylabel="grav. anom.")
scatter!(ax1,xzobs_grav[:,1],tgrav)

ax2 = Axis(fig[2,1],title="Total-field magnetic anomaly",xlabel="x",ylabel="mag. anom.")
scatter!(ax2,xzobs_mag[:,1],tmag)

ax3 = Axis(fig[3,1],title="Polygonal bodies",xlabel="x",ylabel="z")
x1 = [pbody.geom.bo[1].ver1[:,1]...,pbody.geom.bo[1].ver2[end,1]]
y1 = [pbody.geom.bo[1].ver1[:,2]...,pbody.geom.bo[1].ver2[end,2]]
scatterlines!(ax3,x1,y1)
x2 = [pbody.geom.bo[2].ver1[:,1]...,pbody.geom.bo[2].ver2[end,1]]
y2 = [pbody.geom.bo[2].ver1[:,2]...,pbody.geom.bo[2].ver2[end,2]]
scatterlines!(ax3,x2,y2)
ax3.yreversed=true

linkxaxes!(ax1,ax2,ax3)
display(fig)

save("jointgravmag.svg",fig) # hide
```

![](images/jointgravmag.svg)

## Public API

### Module
```@docs
MagGrav2Dpoly
```

### Data structures
```@docs
MagPolygBodies2D
GravPolygBodies2D
JointPolygBodies2D
```

!!! warning 
    Vertices of the polygonal bodies must be provided 
    in **counterclockwise order** to the function `JointPolygBodies2D`
	in order to perform gravity and magnetic anomaly calculations.



### Forward functions for magnetics
```@docs
tmagpolybodies2D
tmagpolybodies2Dgen
```

### HMC helper functions for magnetics
```@docs
Mag2Dpoly.HMCMag2Dpoly
Mag2Dpoly.Mag2DpolyProb
```

#### Misfit structure & functions for magnetics
```@docs
Mag2DPolyMisf
precalcADstuffmag
calcmisfmag
calc∇misfmag
```

### Useful functions for magnetics
!!! note
    These functions are not exported. To call them
    type `Mag2Dpoly.` before the name of the functions.
	
```@docs
Mag2Dpoly.convert_H_to_B_nT
Mag2Dpoly.convert_B_nT_to_H
Mag2Dpoly.magcomp
```


### Forward functions for gravity
```@docs
tgravpolybodies2D
tgravpolybodies2Dgen
```

### HMC helper functions for gravity
```@docs
Grav2Dpoly.HMCGrav2Dpoly
Grav2Dpoly.Grav2DpolyProb
```

#### Misfit structure & functions for gravity
```@docs
Grav2DPolyMisf
precalcADstuffgrav
calcmisfgrav
calc∇misfgrav
```

### Forward functions for joint mag and grav calculations
```@docs
tjointpolybodies2D
tjointpolybodies2Dgen
```

### HMC helper functions for joint mag and grav
```@docs
JointMagGrav2Dpoly.HMCJointMagGrav2Dpoly
JointMagGrav2Dpoly.Joint2DpolyProb
```

