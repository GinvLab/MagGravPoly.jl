 
using GLMakie


mutable struct InteractionState
    points::Observable{Vector{<:Point}}
end

###############################################################


function drawpoly( ; markersize=20.0)
    
    points = Observable(Point2f[(-2,-2), (2,-2), (2,2),]) #Node

    # initialize struct to hold variable states
    state = InteractionState(points)


    # create the figure
    fig = Figure(size = (1000, 900))
    ax1 = Axis(fig[1,1])

    # plot the initial polygon
    poly!(ax1,state.points, strokewidth=2, strokecolor=:black, color=:skyblue2) 
          #scale_plot = false)
    # plot vertices as a scatter plot
    scatplt = scatter!(ax1,state.points, color=:white, strokewidth=3, markersize=markersize,
             strokecolor=:black)  #, raw=true)

    # resize figure to polygon
    #autolimits!(ax1)
    limits!(ax1,-10,10,-10,10)

    # remove the click and drag zoom action
    deregister_interaction!(ax1, :rectanglezoom)
    
    # # add_move!(scene, points, pplot)
    add_remove_move_points!(ax1, points)

    #mouseevents = addmouseevents!(ax1.scene,scatplt)
    #move_points!(ax1,points,scatplt,mouseevents)

    #move_points!(ax1,points)
    # center!(fig.scene)

    display(fig)

    # ##--------------------------------------------
    ## record a movie for given number of frames
    # nfp = 100
    # fps = 10
    # sleep(5)
    # record(fig.scene, "test.mp4"; framerate = fps) do io
    #     println("Start recording")
    #     for i = 1:nfp
    #         @show i
    #         sleep(1/fps)
    #         recordframe!(io)
    #     end
    #     println("End recording")
    # end
    
    return state
end

#######################################################

"""

Compute shortest distance from point to segment.
a"""
@inline function distPt2Segm(xpt::Real,ypt::Real,
                             xs1::Real,ys1::Real,xs2::Real,ys2::Real)
    # xpt,ypt is the point, xs,ys the segment

    ## http://paulbourke.net/geometry/pointlineplane/
    px = xs2-xs1
    py = ys2-ys1
    pnor = px^2 + py^2

    u =  ((xpt - xs1) * px + (ypt - ys1) * py) / pnor
    if u > 1
        u = 1
    elseif u < 0
        u = 0
    end
    x = xs1 + u * px
    y = ys1 + u * py

    dx = x - xpt
    dy = y - ypt

    dist = sqrt(dx^2 + dy^2)
    return dist
end

#######################################################

"""

Add point to the closest nearest egde, splitting the latter in two.
"""
function addpoint2nearestedge!(ax,points)

    # moupos = mouseposition(ax.scene,.events.mouseposition[]
    # pos = to_world(ax.scene, Point2f0(moupos))
    pos = mouseposition(ax.scene)

    N = length(points[])
    mindist = +Inf
    xpt,ypt = pos
    idx = nothing
    twoedgessamedist = false
    dist = zeros(N)
    for i=1:N
        if i<N
            xs1,ys1 = points[][i]
            xs2,ys2 = points[][i+1]
        else
            xs1,ys1 = points[][i]
            xs2,ys2 = points[][1]
        end
        curdist = distPt2Segm(xpt,ypt,xs1,ys1,xs2,ys2)
        dist[i] = curdist
        if curdist<mindist
            mindist=curdist
        end
    end
    # get the index of for minimum distance and check
    #   whether the minimun distance is unique
    counter = 0
    for i=1:N
        if dist[i]==mindist
            idx=i
            counter+=1
        end
    end
    if counter>1
        ## special case when two (or more) points have the
        ##   same distance
        println("Point has the same distance to two (or more) edges, cannot determine \n to which edge should add it. Please select another point.")
        #@show dist
    elseif counter==0
        error("addpoint2nearesteadge(): Something went wrong, no shortest distance found.")
    else
        # index of the closest edge
        points[] = insert!(points[],idx+1,pos)
    end
    return 
end

#######################################################

function removepoint!(ax,points)

    # get which plot and id the mouse is over
    plot, idx = pick(ax.scene)

    ##########################################
    ##   Check the following!!!             ##   
    ##########################################
    # take action only if the plot is of type Scatter
    targetplot = Scatter{Tuple{Vector{Point{2, Float32}}}}  #<<<<<--------<<<<<<

    if typeof(plot) == targetplot && checkbounds(Bool, points[], idx)
        
        if length(points[]) > 3
            # /\ only if there are already at least 4 points
            deleteat!(points[], idx)
            # next line needed to trigger Observable
            points[] = points[]
        else
            println("Only 3 points left, not deleting. Try moving edges instead.")
        end
    end

    return
end


#######################################################
#=
 function movepoint!(ax,points)

     # get which plot and id the mouse is over
     plot, idx = pick(ax.scene)

     ##########################################
     ##   Check the following!!!             ##   
     ##########################################
     # take action only if the plot is of type Scatter
     targetplot = Scatter{Tuple{Vector{Point{2, Float32}}}}  #<<<<<--------<<<<<<

     if typeof(plot) == targetplot && checkbounds(Bool, points[], idx)
         pos = mouseposition(ax.scene)
         points[][idx] = Point2f(pos)
     end
     points[] = points[]

     return
 end
=#

##########################################################

"""
Add, remove or move points in a polygon interactively.

Controls:
- [A] + Left Click: Add point near closest edge
- [D] + Left Click: Remove point
- [S] + Left Click and Drag: Move point
"""
function add_remove_move_points!(ax, points)

    global idx
    idx = 0  # initialize selected index for dragging
    scene = ax.scene

    # Mouse press handler
    on(events(scene).mousebutton) do mouse_event
        if mouse_event.button == Mouse.left
            pos = mouseposition(scene)

            if mouse_event.action == Mouse.press
                if ispressed(scene, Keyboard.a)
                    # Add a new vertex to polygon
                    addpoint2nearestedge!(ax, points)

                elseif ispressed(scene, Keyboard.d)
                    # Remove points from polygon
                    removepoint!(ax, points)

                elseif ispressed(scene, Keyboard.s)
                    # Select point for dragging
                    plot, i = pick(scene)
                    if plot isa Scatter && checkbounds(Bool, points[], i)
                        idx = i
                    else
                        idx = 0  # No point selected
                    end
                end

            elseif mouse_event.action == Mouse.release
                # Release point after dragging
                idx = 0
                notify(points)
            end
        end
        return Consume(false)
    end

    # Mouse motion handler (for dragging)
    on(events(scene).mouseposition) do _
        if idx > 0 && ispressed(scene, Keyboard.s)
            pos = mouseposition(scene)
            points[][idx] = Point2f(pos)
            notify(points)
        end
    end

    return nothing
end



#######################################################
#=
"""

Move vertices of polygon.
"""

function move_points!(ax, points, scatplt, mouseevents)

    # idx must be passed somewhat from "onmouseleftdragstart" to "onmouseleftdrag"
    idx = nothing

    # left drag starts
    onmouseleftdragstart(mouseevents) do event
        # get which plot and id the mouse is over
        plot, idx = pick(ax.scene) #mouse_selection

        ##########################################
        ##   Check the following!!!             ##   
        ##########################################
        # take action only if the plot is of type Scatter
        targetplot = Scatter{Tuple{Vector{Point{2, Float32}}}}  #<<<<<------

        if typeof(plot) == targetplot && checkbounds(Bool, points[], idx)
            pos = mouseposition(ax.scene)
            points[][idx] = Point2f(pos)
        end
        points[] = points[]
    end
 
    # left drag continues
    onmouseleftdrag(mouseevents) do event
        pos = mouseposition(ax.scene)
        points[][idx] = Point2f(pos) 
        points[] = points[]
    end

    return
end
=#



# function move_points!(ax, points)
#     on(events(ax.scene).mouseposition) do moudrag
#         moueve = events(ax.scene).mousebutton[]
#         if moueve.button==Mouse.left && moueve.action==Mouse.press
#             movepoint!(ax,points)
#         end
#     end
#     return
# end



#########################################################

# function add_move!(scene, points, pplot)
#     idx = Ref(0); dragstart = Ref(false); startpos = Base.RefValue(Point2f0(0))
#     on(events(scene).mousedrag) do drag
#         if ispressed(scene, Mouse.left)
#             if drag == Mouse.down
#                 plot, _idx = mouse_selection(scene)
#                 if plot == pplot
#                     idx[] = _idx; dragstart[] = true
#                     startpos[] = to_world(scene, Point2f0(scene.events.mouseposition[]))
#                 end
#             elseif drag == Mouse.pressed && dragstart[] && checkbounds(Bool, points[], idx[])
#                 pos = to_world(scene, Point2f0(scene.events.mouseposition[]))
#                 points[][idx[]] = pos
#                 points[] = points[]
#             end
#         else
#             dragstart[] = false
#         end
#         return
#     end
# end

#######################################################
