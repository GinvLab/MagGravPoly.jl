
# GeoPolygons


## User guide

**GeoPolygons** *is a Julia package utility developed for handling polygonal shapes in the framework of 2D to 2.75D potential-fields forward and inverse modeling*. 


## Documentation


### Installation

To install the package first enter into the package manager mode in Julia by typing "`]`" at the 
REPL prompt and add the "JuliaGeoph" registry as
```
(@v1.9) pkg> registry add https://github.com/GinvLab/GinvLabRegistry
```
Then add the package by simply issuing
```
(@v1.9) pkg> add GeoPolygons
```
The package will be automatically downloaded from the web and installed.


### Tutorial

#### Check the orientation of a polygon
First load the module and define a list of vertices of the poligonal bodies 
and the relative indices mapping each body to its vertices:
```@example ex1
using MagGravPoly.GeoPolygons
# vertices of the poligonal bodies
vertices  = [35.0 50.0;
             65.0 50.0;
             80.0 35.0;
             65.0 20.0;
             35.0 20.0;
             20.0 35.0]
			 
# indices of vertices for the body
ind1 = collect(1:6)
bodyindices = [ind1]
# construct the poligonal body object
pbody = PolygBodies2D(bodyindices,vertices)
nothing # hide
```
and then checks if the polygon is defined in a clockwise or counter-clockwise order: 
```@example ex1
# check about the order
for i=1:pbody.nbo
aclk = checkanticlockwiseorder(pbody.bo[i])
@show aclk
end
```

In fact plotting the polygon:

    using PyPlot
    figure()
	for i=1:length(bodyindices)
	x = copy(pbody.bo[i].ver1[:,1])
	y = copy(pbody.bo[i].ver1[:,2])
	append!(x,pbody.bo[i].ver1[1,1])
	append!(y,pbody.bo[i].ver1[1,2])
    p=PyPlot.plot(x,y,"-o")
    for i=1:length(x)-1
	s=string(i)
	if i<=1
	txtalign = "left"
	else
	txtalign = "right"
    end
    text(x[i],y[i],s=s,horizontalalignment=txtalign,
    fontsize=16,verticalalignment="top")
	end
    end
    gca().invert_yaxis()

is easy to see that is counter-clockwise oriented.
![](images/plotex1.svg)

#### Check and fixing for crossing polygon sides
First load the module and define a list of vertices of the poligonal bodies 
and the relative indices mapping each body to its vertices:
```@example ex2
using MagGravPoly.GeoPolygons
# vertices of the poligonal bodies
vertices  = [35.0 50.0;
             80.0 35.0;
             65.0 50.0;
             65.0 20.0;
             35.0 20.0;
             20.0 35.0]
			 
# indices of vertices for the body
ind1 = collect(1:6)
bodyindices = [ind1]
# construct the poligonal body object
pbody = PolygBodies2D(bodyindices,vertices)
nothing # hide
```

Now checks for any crossing polygon side:

```@example ex2
chk = checkpoly(pbody.bo)
@show chk[1]
```

In fact plotting the polygon:

	using PyPlot
    figure()
	for i=1:length(bodyindices)
	x = copy(pbody.bo[i].ver1[:,1])
	y = copy(pbody.bo[i].ver1[:,2])
	append!(x,pbody.bo[i].ver1[1,1])
	append!(y,pbody.bo[i].ver1[1,2])
    p=PyPlot.plot(x,y,"-o")
    for i=1:length(x)-1
	s=string(i)
	if i<=1
	txtalign = "left"
	else
	txtalign = "right"
    end
    text(x[i],y[i],s=s,horizontalalignment=txtalign,
    fontsize=16,verticalalignment="top")
	end
    end
    gca().invert_yaxis()
	
![](images/plotex2.svg)

We can try to fix the polygonal geometries using the `verpolyshift!` function:

```@example ex2
# polygon fixing
GeoPolygons.verpolyshift!(pbody.bo)
```

and the results will be the following:

![](images/plotex2fix.svg)

## Public API
```@docs
MagGravPoly.GeoPolygons
```

### Data structures
```@docs
MagGravPoly.MG2D.BodySegments2D
MagGravPoly.MG2D.PolygBodies2D
```

!!! warning 
    Vertices of the polygonal bodies must be provided 
    counterclockwise to the structure `BodySegments2D`
    to perform magnetic anomaly calculation using the
    functions in the packages `Mag2Dpoly` and `Grav2Dpoly`.
	To assess this use the function `checkanticlockwiseorder`.


```@docs
MagGravPoly.GeoPolygons.TopoEdges
```

### Checking-geometries functions
#### Single polygonal body
```@docs
MagGravPoly.GeoPolygons.intersectpairpoly
MagGravPoly.GeoPolygons.selfintersectpoly
MagGravPoly.GeoPolygons.checkall
MagGravPoly.GeoPolygons.checkanticlockwiseorder
```

#### Multiple polygonal bodies
```@docs
MagGravPoly.GeoPolygons.checktopo
MagGravPoly.GeoPolygons.checkpoly
```

#### Fixing-geometries functions
```@docs
MagGravPoly.GeoPolygons.verpolyshift!
MagGravPoly.GeoPolygons.verpolyallshift!
MagGravPoly.GeoPolygons.vertoposhift!
MagGravPoly.GeoPolygons.fixall!
```

### Useful functions
```@docs
MagGravPoly.GeoPolygons.checkbodyindices
MagGravPoly.GeoPolygons.Inter2Segm
MagGravPoly.GeoPolygons.isInternal
MagGravPoly.GeoPolygons.checkmodelizdim
MagGravPoly.GeoPolygons.calcareapoly
MagGravPoly.GeoPolygons.calcareamanypoly
```
