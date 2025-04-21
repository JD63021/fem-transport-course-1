function precomputeShapeFunctions1D(order::Int, nGP::Union{Int,Nothing}=nothing)
    # 1) Determine number of Gauss points if not provided.
    if nGP === nothing
        if order == 3
            nGP = 5
        else
            nGP = 3
        end
    end

    # 2) Define Gauss integration rule on [0,1]
    gp = Float64[]
    w = Float64[]
    if nGP == 3
        gp = [0.112701665379258, 0.5, 0.887298334620742]
        w  = [0.277777777777778, 0.444444444444444, 0.277777777777778]
    elseif nGP == 5
        gp = [0.04691008, 0.23076534, 0.5, 0.76923466, 0.95308992]
        w  = [0.11846344, 0.23931434, 0.28444444, 0.23931434, 0.11846344]
    else
        error("nGP must be either 3 or 5.")
    end
    numGauss = length(gp)

    # 3) Set number of local nodes based on order.
    numNodes = 0
    if order == 1
        numNodes = 2  # Linear
    elseif order == 2
        numNodes = 3  # Quadratic
    elseif order == 3
        numNodes = 4  # Cubic
    else
        error("Order must be 1, 2, or 3.")
    end

    # 4) Allocate output arrays for shape function values (N) and derivatives (dN)
    N = zeros(numNodes, numGauss)
    dN = zeros(numNodes, numGauss)

    # 5) Compute shape functions and their derivatives based on element order.
    if order == 1
        # Linear elements: N₁ = 1 - ξ, N₂ = ξ with constant derivatives.
        for ig in 1:numGauss
            xi = gp[ig]
            N[1, ig] = 1 - xi
            N[2, ig] = xi

            dN[1, ig] = -1.0
            dN[2, ig] =  1.0
        end

    elseif order == 2
        # Quadratic elements (GMsh ordering: N1, N3, N2)
        # Polynomials are chosen so that:
        #   At ξ = 0: N1 = 1, N2 = 0, N3 = 0;
        #   At ξ = 0.5: N1 = 0, N2 = 0, N3 = 1;
        #   At ξ = 1: N1 = 0, N2 = 1, N3 = 0.
        for ig in 1:numGauss
            xi = gp[ig]
            N[1, ig] = 1 - 3*xi + 2*xi^2   # Node at 0.
            N[3, ig] = 4*xi - 4*xi^2         # Mid node (at 0.5).
            N[2, ig] = -xi + 2*xi^2          # Node at 1.
            
            dN[1, ig] = -3 + 4*xi
            dN[3, ig] = 4 - 8*xi
            dN[2, ig] = -1 + 4*xi
        end

    elseif order == 3
        # Cubic elements.
        # Standard nodal positions for a cubic element: [0, 1/3, 2/3, 1]
        nodes_std = [0.0, 1/3, 2/3, 1.0]
        
        # Precompute denominators for each basis function.
        denom = zeros(4)
        for i in 1:4
            prod_val = 1.0
            for j in 1:4
                if j != i
                    prod_val *= (nodes_std[i] - nodes_std[j])
                end
            end
            denom[i] = prod_val
        end

        # Allocate temporary arrays for standard Lagrange basis functions.
        L = zeros(4, numGauss)
        dL = zeros(4, numGauss)
        
        # Compute the standard Lagrange basis functions and their derivatives.
        for ig in 1:numGauss
            xi_val = gp[ig]
            for i in 1:4
                # Compute L_i(ξ)
                prod_val = 1.0
                for j in 1:4
                    if j != i
                        prod_val *= (xi_val - nodes_std[j])
                    end
                end
                L[i, ig] = prod_val / denom[i]
                
                # Compute derivative dL_i(ξ) using the product rule.
                sum_val = 0.0
                for k in 1:4
                    if k != i
                        prod_der = 1.0
                        for j in 1:4
                            if j != i && j != k
                                prod_der *= (xi_val - nodes_std[j])
                            end
                        end
                        sum_val += prod_der
                    end
                end
                dL[i, ig] = sum_val / denom[i]
            end
        end

        # Reorder the basis functions to match GMsh ordering.
        # Standard ordering: [L1, L2, L3, L4] corresponding to nodes [0, 1/3, 2/3, 1].
        # GMsh ordering: [N1, N4, N2, N3] i.e. first corner, second corner, then interior nodes.
        N[1, :] = L[1, :]  # Node at 0.
        N[2, :] = L[4, :]  # Node at 1.
        N[3, :] = L[2, :]  # Node at 1/3.
        N[4, :] = L[3, :]  # Node at 2/3.
        
        dN[1, :] = dL[1, :]
        dN[2, :] = dL[4, :]
        dN[3, :] = dL[2, :]
        dN[4, :] = dL[3, :]
    end

    return N, dN, gp, w
end
