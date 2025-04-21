function build_residual(nodeInfo, elemInfo, boundaryInfo, order, kappa, T, bcFlags, bcFunction, f)
    # (1) Basic data
    allNodes    = nodeInfo[:x]               # 1D coordinates of each node
    numNodes    = length(allNodes)
    elements    = elemInfo[:elements]          # (#elements x (order+1))
    numElements = size(elements, 1)
    
    # Initialize global residual vector.
    R = zeros(numNodes)
    
    # (2) Precompute 1D shape functions
    N, dN, gp, w = precomputeShapeFunctions1D(order)
    numGauss = length(w)
    nLoc = order + 1   # number of local nodes per element
    
    # (3) Element loop: assemble interior residual.
    for e in 1:numElements
        Rloc = zeros(nLoc)
        # Local node indices for element e.
        K = elements[e, :]
        # Physical coordinates for these nodes.
        xcoords = allNodes[K]
        # Local T values.
        T_loc = T[K]
        
        for ig in 1:numGauss
            xi = gp[ig]
            wt = w[ig]
            
            # Compute dx/dxi (Jacobian)
            dx_dxi = 0.0
            for iNode in 1:nLoc
                dx_dxi += xcoords[iNode] * dN[iNode, ig]
            end
            detJ = abs(dx_dxi)
            
            # Compute dN/dx = (1/dx_dxi) * dN/dξ
            dNdx = (1 / dx_dxi) .* dN[:, ig]
            
            # Compute dT/dx at this Gauss point.
            dTdx = 0.0
            for iNode in 1:nLoc
                dTdx += dNdx[iNode] * T_loc[iNode]
            end
            
            # Compute physical coordinate at Gauss point: x_gp = Σ (N_i(ξ) * x_i)
            x_gp = dot(N[:, ig], xcoords)
            # Evaluate forcing term. If f is a function, call it; otherwise assume it's a constant.
            force_val = f isa Function ? f(x_gp) : f
            
            # Contribution to local residual.
            for iNode in 1:nLoc
                Rloc[iNode] += (kappa * (dTdx * dNdx[iNode])) * wt * detJ - force_val * N[iNode, ig] * wt * detJ
            end
        end
        
        # Accumulate local residual into global residual.
        for iNode in 1:nLoc
            global_i = K[iNode]
            R[global_i] += Rloc[iNode]
        end
    end
    
    # (4) Impose Dirichlet BCs: for each flagged boundary, set R(node) = T(node) - bcValue.
    # (A) Inlets.
    if haskey(bcFlags, :inlet)
        inFlags = bcFlags[:inlet]
        if !(inFlags isa AbstractVector)
            inFlags = [inFlags]
        end
        allInletNodes = Int[]
        for flagVal in inFlags
            fn = Symbol("flag_" * string(flagVal))
            if haskey(boundaryInfo, fn)
                append!(allInletNodes, boundaryInfo[fn])
            end
        end
        allInletNodes = unique(allInletNodes)
        for nd in allInletNodes
            # Here, bcFunction is used to supply the desired Dirichlet value.
            R[nd] = T[nd] - bcFunction
        end
    end

    # (B) Walls.
    if haskey(bcFlags, :wall)
        wFlags = bcFlags[:wall]
        if !(wFlags isa AbstractVector)
            wFlags = [wFlags]
        end
        allWallNodes = Int[]
        for flagVal in wFlags
            fn = Symbol("flag_" * string(flagVal))
            if haskey(boundaryInfo, fn)
                append!(allWallNodes, boundaryInfo[fn])
            end
        end
        allWallNodes = unique(allWallNodes)
        for nd in allWallNodes
            R[nd] = T[nd]  # Assuming Dirichlet wall BC of zero.
        end
    end

    return R
end
