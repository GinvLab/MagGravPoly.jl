
# http://makie.juliaplots.org/v0.15.2/examples/layoutables/axis/index.html#custom_interactions
# http://makie.juliaplots.org/v0.15.2/documentation/events/index.html

# http://juliaplots.org/MakieReferenceImages/gallery//edit_polygon/index.html
# interactions(ax1)


# http://juliaplots.org/MakieReferenceImages/gallery//edit_polygon/index.html

## custom interaction - axis
## https://makie.juliaplots.org/v0.15.2/examples/layoutables/axis/index.html#custom_interactions

using GLMakie

############################################

mutable struct InteractionState
    polygon::Union{Matrix{Float64},Nothing}
    clicks::Observable{Matrix{Float64}}
    points1::Observable{Matrix{Float64}}
    points2::Observable{Matrix{Float64}}
    points3::Observable{Matrix{Float64}}
end

############################################

function plotincrpolyg(ax1,points1,points2,points3,markersize)
    
    # plot edges up to previous point
    lines!(ax1,points1, color=:red,linewidth=3)
    # plot vertices up to previous point
    scatter!(ax1,points1, color=:red, markersize=markersize)
    # plot edges for last two points
    lines!(ax1,points2, color=:orange,linewidth=3)
    # plot vertex for last point
    scatter!(ax1,points3, color=:orange, markersize=markersize)       

    return
end

############################################################
#
# Do some stupid splitting of drawpolyg into
#  several other observables to be able to
#  automatically update plots. This is because
#  the plotting functions require the argument to
#  to be an Observable/Node and not its value in
#  order to be able to actually perform the update.
#

function splitobs1!(drawpolyg)
    if size(drawpolyg,1) >= 1
        # transpose because the plotting routines expect
        #  x and y as ROWs (for dim=2x2, otherwise...)
        points1 = transpose(drawpolyg[1:end-1,:])
    else
        points1 = drawpolyg[:,:]
    end
    return points1
end

##---------------------------------------------------

function splitobs2!(drawpolyg)
    
    if size(drawpolyg,1) >= 2
        # transpose because the plotting routines expect
        #  x and y as ROWs (for dim=2x2, otherwise...)
        points2 = transpose(drawpolyg[end-1:end,:])
    else
        points2 = drawpolyg[:,:]
    end
    return points2
end

##---------------------------------------------------

function splitobs3!(drawpolyg)
    if size(drawpolyg,1) >= 1
        # transpose because the plotting routines expect
        #  x and y as ROWs (for dim=2x2, otherwise...)
        points3 = transpose(drawpolyg[end:end,:])
    else
        points3 = drawpolyg[:,:]
    end
    return points3
end


#########################################################

function findnearestvertex(p::Vector{Float64},pts::Matrix{Float64})
    npt=size(pts,1)
    mindist = Inf
    dist = 0.0
    idx = 1
    for i=1:npt
        dist = sqrt( (p[1] .- pts[i,1]).^2 + (p[2] .- pts[i,2]).^2 )
        if dist<mindist
            mindist = dist
            idx = i
        end
    end
    return idx,mindist
end


#########################################################

function isonmarker(moupos,distmm_wld,markersize_px,ax1)

    testpt_px = Point2f(0.0,markersize_px)
    orig_wld = to_world(ax1.scene,Point2f((0.0,0.0)))
    markersize_wld = Float64( sqrt( sum( (to_world(ax1.scene,testpt_px) .- orig_wld).^2 )) )

    if distmm_wld <= markersize_wld
        onmar = true
    else
        onmar = false
    end
    return onmar
end

#########################################################

function dragvertex!(clicks,mpos,markersize_px,ax1)
    # compute closest vertex to mouse position
    # outputs index of clicks and distance
    idx,dist_wld = findnearestvertex(mpos,clicks)

    # is the mouse within the marker area?
    onmar = isonmarker(mpos,dist_wld,markersize_px,ax1)

    # if mouse inside marker, then relocate the vertex
    if onmar
        #println()
        # @show mpos
        # @show clicks
        #println("on marker")
        clicks[idx,:] = mpos
        #@show clicks
    end

    return clicks
end

#########################################################
    

function MakieLayout.process_interaction(state::InteractionState, event::Union{MouseEvent,KeysEvent}, ax1)
    

    markersize = 30
    markersize_px = Float64(markersize)

    if typeof(event)==MouseEvent

        ##=======================================
        ## mouse events

        if event.type === MouseEventTypes.leftclick

            ##------------------------------
            ## remove last point
            #leftclickcounter[] += 1
            # convert to array of Float64
            pos = convert(Array{Float64},event.data)
            # add a row to clicks with the new point
            state.clicks[] = vcat(state.clicks[], pos')                
            println(" New point: $(pos)") #, # $(leftclickcounter[])")
            ## trigger plotting on first click(s?)
            if size(state.clicks[],1) <= 3
                plotincrpolyg(ax1,state.points1,state.points2,state.points3,
                              markersize)
            end
            
            
        elseif event.type === MouseEventTypes.rightclick

            ##------------------------------
            ## remove last point
            if size(state.clicks[],1) >= 1
                println("Removing last point from polygon.")
                if size(state.clicks[],1) == 1
                    # reset to 0 size
                    state.clicks[] = state.clicks[][1:0,:]
                elseif size(state.clicks[],1) > 1
                    # remove last row
                    state.clicks[] = state.clicks[][1:end-1,:]                        
                end
            else
                println("Attempting to remove last point, but no points registered yet.")
            end

            
        elseif event.type === MouseEventTypes.leftdrag

            
            ##------------------------------
            ## relocate vertex
            if size(state.clicks[],1) >= 1

                mpos = mouseposition(ax1.scene)
                mpos1 = mouseposition(ax1.scene)
                mpos = convert(Vector{Float64},mpos1)

                state.clicks[] = dragvertex!(state.clicks[],mpos,markersize_px,ax1)
            end

        end


    elseif typeof(event)==KeysEvent
        ##=======================================
        # keyboard events
 
        if Keyboard.a in event.keys
            println("a - autolimits")
            autolimits!(ax1)

        elseif Keyboard.f in event.keys
            println("f - close polygon")
            autolimits!(ax1)
            # add fist point to the end to trigger plotting of closed polygon
            state.clicks[] = vcat(state.clicks[],state.clicks[][1:1,:])
            # store a non-closed polygon
            #state.polygon = nothing
            state.polygon = state.clicks.val[1:end-1,:]
            #return drawpolyg
            #deregister_interaction!(ax1, :my_interaction)
            
         elseif Keyboard.space in event.keys
            println("f - close polygon")


        # elseif Keyboard.d in event.keys
        #     # deactivate the interaction if "d" is pressed
        #     deactivate_interaction!(ax1, :my_interaction)

        # elseif Keyboard.r in event.keys
        #     # reactivate the interaction if "r" is pressed
        #     reactivate_interaction!(ax1, :my_interaction)

        # elseif Keyboard.e in event.keys
        #     autolimits!(ax1)
        #     # deregister the interaction if "a" is pressed
        #     deregister_interaction!(ax1, :my_interaction)

        end        
    end
    

    return 
end

#########################################################################

############################################ 

function drawpolygons()
    
    ###################
    ##  Init
    fig = Figure(resolution = (1000, 700))
    ax1 = Axis(fig[1,1])

    xlims!(1,100); ylims!(1,100)



    ########################################
    ## Interactions
    clicks = Node(Array{Float64,2}(undef,0,2))
   
    ## LIFT function splitobs!
    points1 = lift(splitobs1!,clicks)
    points2 = lift(splitobs2!,clicks)
    points3 = lift(splitobs3!,clicks)

    polygon = Array{Float64,2}(undef,0,2)
    #leftclickcounter = Node(0)

    # register_interaction!(ax1, :my_interaction) do event::Union{MouseEvent,KeysEvent}, axis 

    # Removing the rectangle zoom:  deregister_interaction!(ax, :rectanglezoom).
    deregister_interaction!(ax1, :rectanglezoom)
        
    register_interaction!(ax1, :itracepoly, InteractionState(polygon,clicks,points1,points2,points3))

    
    display(fig)
    
    #--------------------------------------------
    ## record stuff

    @show interactions(ax1)
    @show typeof(:itracepoly)

    #--------------------------------------------
    # record events
    #record_events(:itracepoly, fig.scene, "testevents.mp4")
    
    
    #--------------------------------------------
    # record a movie for given number of frames
    # nfp = 300
    # fps = 10
    # record(fig.scene, "test.mp4"; framerate = fps) do io
    #     println("Start recording")
    #     for i = 1:nfp
    #         @show i
    #         sleep(1/fps)
    #         recordframe!(io)
    #     end
    #     println("End recording")
    # end
    

    return 
end


##########################################################
