using SparseArrays

function build_jacobian_s3(nodeInfo, elemInfo, boundaryInfo, order, kappa, T, bcFlags, bcFunction, f)
    # (1) Basic data
    allNodes    = nodeInfo[:x]               # 1D coordinates for all nodes
    numNodes    = length(allNodes)
    elements    = elemInfo[:elements]          # (#elements x (order+1)) connectivity matrix
    numElements = size(elements, 1)
    
    # (2) Precompute shape functions on [0,1]
    N, dN, gp, w = precomputeShapeFunctions1D(order)
    numGauss = length(w)
    nLoc = order + 1  # number of local nodes per element
    
    # Preallocate triplet arrays for assembling the global sparse matrix.
    nTrip = nLoc * nLoc * numElements
    I = zeros(Int, nTrip)
    J = zeros(Int, nTrip)
    V = zeros(nTrip)
    offset = 0

    # (3) Loop over each element to assemble its local stiffness matrix.
    for e in 1:numElements
        # Initialize local stiffness matrix (Kloc)
        Kloc = zeros(nLoc, nLoc)
        # Get local node indices for element e.
        Kidx = elements[e, :]
        # Get physical coordinates of these nodes.
        xcoords = allNodes[Kidx]
        
        for ig in 1:numGauss
            xi = gp[ig]
            wt = w[ig]
            
            # Compute the Jacobian dx/dξ = Σ (x_i * dN_i/dξ)
            dx_dxi = 0.0
            for iNode in 1:nLoc
                dx_dxi += xcoords[iNode] * dN[iNode, ig]
            end
            detJ = abs(dx_dxi)
            
            # Compute dN/dx = (1/dx_dxi) * dN/dξ at this Gauss point.
            dNdx = (1 / dx_dxi) .* dN[:, ig]
            
            # Accumulate contributions into the local stiffness matrix.
            for iNode in 1:nLoc
                for jNode in 1:nLoc
                    Kloc[iNode, jNode] += kappa * (dNdx[iNode] * dNdx[jNode]) * wt * detJ
                end
            end
        end
        
        # Store local contributions using triplet format.
        # Create index arrays for the local block.
        rr = repeat(Kidx, inner = nLoc)   # Each local node index repeated nLoc times.
        cc = repeat(Kidx, outer = nLoc)    # Each local node index repeated for each entry.
        blockSize = nLoc * nLoc
        idxRange = offset + 1 : offset + blockSize
        I[idxRange] = rr
        J[idxRange] = cc
        V[idxRange] = vec(Kloc)
        offset += blockSize
    end
    
    # (4) Assemble the global sparse matrix.
    K = sparse(I, J, V, numNodes, numNodes)
    
    # (5) Impose Dirichlet boundary conditions: for each BC flag,
    #      zero out the row and set the diagonal entry to 1.
    # (A) Inlet BCs
    if haskey(bcFlags, :inlet)
        inFlags = bcFlags[:inlet]
        if !(inFlags isa AbstractVector)
            inFlags = [inFlags]
        end
        allInletNodes = Int[]
        for fVal in inFlags
            fn = Symbol("flag_" * string(fVal))
            if haskey(boundaryInfo, fn)
                append!(allInletNodes, boundaryInfo[fn])
            end
        end
        allInletNodes = unique(allInletNodes)
        for nd in allInletNodes
            K[nd, :] .= 0.0  # Zero out the entire row.
            K[nd, nd] = 1.0  # Set the diagonal to 1.
        end
    end

    # (B) Wall BCs
    if haskey(bcFlags, :wall)
        wFlags = bcFlags[:wall]
        if !(wFlags isa AbstractVector)
            wFlags = [wFlags]
        end
        allWallNodes = Int[]
        for fVal in wFlags
            fn = Symbol("flag_" * string(fVal))
            if haskey(boundaryInfo, fn)
                append!(allWallNodes, boundaryInfo[fn])
            end
        end
        allWallNodes = unique(allWallNodes)
        for nd in allWallNodes
            K[nd, :] .= 0.0
            K[nd, nd] = 1.0
        end
    end

    return K
end
