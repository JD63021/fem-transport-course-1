using LinearAlgebra
include("build_residual.jl")
include("build_jacobian.jl")
include("precomputeShapeFunctions1D.jl")
# =============================================================================
# main1_SS - Stationary solver for the 1D Poisson/heat problem.
#
# Inputs:
#   U              : Initial solution vector (guess).
#   D              : Diffusivity.
#   nodeInfo       : Dictionary with node coordinates (e.g., nodeInfo[:x]).
#   elemInfo       : Dictionary with element connectivity.
#   boundaryInfo   : Dictionary with boundary sets (e.g., boundaryInfo[:flag_...]).
#   time           : Simulation time (used for evaluating BCs).
#   bcFlags        : Dictionary with boundary flag info (e.g., :inlet and :wall).
#   inletProfile   : Either a numeric constant or a function (t, x) → value.
#   order          : Element order (1, 2, or 3).
#   source         : Source term.
#
# Outputs:
#   U              : Final solution vector.
#   x              : Node coordinates.
#   nodeInfo, elemInfo, boundaryInfo : (Possibly updated) mesh data.
# =============================================================================
function main1_SS(U, D, nodeInfo, elemInfo, boundaryInfo, time, bcFlags, inletProfile, order, source)
    # Number of degrees of freedom (assumed one per node in 1D)
    Nxy = length(nodeInfo[:x])
    currentTime = time  # For a stationary problem, time is arbitrary

    # Set up Newton iteration parameters.
    oldResNorm = 1e6
    i = 0
    maxNewtonIters = 200
    tolRes = 1e-6         # Residual tolerance
    tolSolBlowUp = 1e40   # Blow-up threshold
    rnorm1 = 1e9
    rnorm2 = 1e9

    # Update boundary conditions before starting the iteration.
    U = update_bc(U, boundaryInfo, nodeInfo, Nxy, currentTime, bcFlags, inletProfile)

    while (rnorm1 + rnorm2 >= tolRes)
        # Build the Jacobian.
        # (Assumes build_jacobian_s3 is defined elsewhere.)
        J1 = build_jacobian_s3(nodeInfo, elemInfo, boundaryInfo, order, D, U, bcFlags, inletProfile, source)

         

        # Check iteration limits or potential blow-up.
        if i > maxNewtonIters
            println("Max Newton iterations reached.")
            break
        end
        if oldResNorm >= tolSolBlowUp
            println("Solution blew up, stopping.")
            break
        end

        # Build the residual.
        # (Assumes build_residual is defined elsewhere.)
        FF = build_residual(nodeInfo, elemInfo, boundaryInfo, order, D, U, bcFlags, inletProfile, source)

        # Compute the Newton update.
        deltaU = J1 \ FF
        Uprev = copy(U)

        # Use a constant under-relaxation parameter alpha.
        alpha = 1.0  # Adjust as needed (e.g., 0.8 for under-relaxation)
        U = U - alpha * deltaU

        rnorm1 = norm(Uprev - U)
        rnorm2 = norm(FF)

        i += 1
        println("NS iteration ", i, 
                ": relative norm = ", rnorm1, 
                ", residual norm = ", rnorm2, 
                ", alpha = ", alpha)

        if rnorm1 < tolRes
            break
        end
    end

    # Set x as the node coordinates.
    x = nodeInfo[:x]

    return U, x, nodeInfo, elemInfo, boundaryInfo
end

# =============================================================================
# update_bc - Helper function that applies Dirichlet boundary conditions.
#
# This function updates the solution vector uv by imposing BCs according to the 
# boundary information provided. In 1D, each node has one DOF.
#
# Inputs:
#   uv           : Solution vector (e.g., temperature or velocity).
#   boundaryInfo : Dictionary with boundary sets (e.g., boundaryInfo[:flag_...]).
#   nodeInfo     : Dictionary with node coordinates (e.g., nodeInfo[:x]).
#   Nxy          : Number of degrees of freedom (number of nodes).
#   t            : Current time (if the inlet profile depends on time).
#   bcFlags      : Dictionary with keys :inlet and :wall, each mapping to an array of flag IDs.
#   inletProfile : Either a numeric constant or a function (t, x) → value.
#
# Output:
#   uv_new       : Updated solution vector with BCs imposed.
# =============================================================================
function update_bc(uv, boundaryInfo, nodeInfo, Nxy, t, bcFlags, inletProfile)
    uv_new = copy(uv)

    # Helper: For 1D, the global degree-of-freedom index is the node number.
    globalRow_1D(vNode) = vNode

    # -------------------------
    # 1) Apply "inlet" BC (Dirichlet)
    if haskey(bcFlags, :inlet)
        inFlagList = bcFlags[:inlet]
        # Ensure inFlagList is an array.
        if !(inFlagList isa AbstractArray)
            inFlagList = [inFlagList]
        end
        for flagVal in inFlagList
            flagStr = Symbol("flag_" * string(flagVal))
            inletNodes = haskey(boundaryInfo, flagStr) ? boundaryInfo[flagStr] : Int[]
            for nodeID in inletNodes
                rowX = globalRow_1D(nodeID)
                xVal = nodeInfo[:x][nodeID]
                # If inletProfile is a function, call it with (t, xVal); else, use as constant.
                if inletProfile isa Function
                    Uin = inletProfile(t, xVal)
                else
                    Uin = inletProfile
                end
                uv_new[rowX] = Uin
            end
        end
    end

    # -------------------------
    # 2) Apply "wall" BC (e.g., zero value)
    if haskey(bcFlags, :wall)
        wFlagList = bcFlags[:wall]
        if !(wFlagList isa AbstractArray)
            wFlagList = [wFlagList]
        end
        for flagVal in wFlagList
            flagStr = Symbol("flag_" * string(flagVal))
            if !haskey(boundaryInfo, flagStr)
                continue
            end
            sideNodes = unique(boundaryInfo[flagStr])
            for nodeID in sideNodes
                rowX = globalRow_1D(nodeID)
                uv_new[rowX] = 0.0
            end
        end
    end

    return uv_new
end
