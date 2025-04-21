using Plots

function plot1D_FEM_BasisAndLegendre(nodeInfo, elemInfo, order, U=nothing)
    # (A) Plot shape functions over the physical domain.
    xNodes   = nodeInfo[:x]
    elements = elemInfo[:elements]
    numElems = size(elements, 1)
    nLoc     = order + 1
    nPlot    = 25
    xiPlot   = range(0, stop=1, length=nPlot)

    p1 = plot(title="Lagrange Shape Functions in Physical Domain",
              xlabel="x", ylabel="φᵢ(x)", legend=:best)
    for e in 1:numElems
        K  = elements[e, :]
        xE = xNodes[K]
        xMin = minimum(xE)
        xMax = maximum(xE)
        for iLoc in 1:nLoc
            yVals = zeros(nPlot)
            xPlotEl = zeros(nPlot)
            for (ip, xiVal) in enumerate(xiPlot)
                (phi_i, _) = shapeFunction_i(order, iLoc, xiVal)
                Xphys = xMin + xiVal * (xMax - xMin)
                yVals[ip] = phi_i
                xPlotEl[ip] = Xphys
            end
            plot!(p1, xPlotEl, yVals, lw=1.2,
                  label="ϕ_{(elem$(e), loc$(iLoc))}")
        end
    end
    display(p1)

    # (B) Plot derivatives of shape functions.
    p2 = plot(title="Derivatives of Lagrange Shape Functions in Physical Domain",
              xlabel="x", ylabel="dφᵢ/dx", legend=:best)
    for e in 1:numElems
        K  = elements[e, :]
        xE = xNodes[K]
        xMin = minimum(xE)
        xMax = maximum(xE)
        for iLoc in 1:nLoc
            yVals = zeros(nPlot)
            xPlotEl = zeros(nPlot)
            for (ip, xiVal) in enumerate(xiPlot)
                (_, dphi_dxi) = shapeFunction_i(order, iLoc, xiVal)
                # Compute dx/dξ for the element:
                dx_dxi = 0.0
                for jLoc in 1:nLoc
                    (_, dNj) = shapeFunction_i(order, jLoc, xiVal)
                    dx_dxi += xE[jLoc] * dNj
                end
                dphi_dx = dphi_dxi / dx_dxi  # Convert derivative to physical coordinate
                Xphys = xMin + xiVal * (xMax - xMin)
                yVals[ip] = dphi_dx
                xPlotEl[ip] = Xphys
            end
            plot!(p2, xPlotEl, yVals, lw=1.2,
                  label="dφ_{(elem$(e), loc$(iLoc))}/dx")
        end
    end
    display(p2)
end


function shapeFunction_i(order::Int, iLoc::Int, xi::Float64)
    if order == 1
        if iLoc == 1
            return (1 - xi, -1.0)
        elseif iLoc == 2
            return (xi, 1.0)
        else
            error("iLoc out of range for order 1.")
        end
    elseif order == 2
        if iLoc == 1
            return (1 - 3*xi + 2*xi^2, -3 + 4*xi)
        elseif iLoc == 2
            return (-xi + 2*xi^2, -1 + 4*xi)
        elseif iLoc == 3
            return (4*xi - 4*xi^2, 4 - 8*xi)
        else
            error("iLoc out of range for order 2.")
        end
    elseif order == 3
        # Cubic: use standard nodal positions and GMsh ordering.
        nodes_std = [0.0, 1/3, 2/3, 1.0]
        idx = 0
        if iLoc == 1
            idx = 1
        elseif iLoc == 2
            idx = 4
        elseif iLoc == 3
            idx = 2
        elseif iLoc == 4
            idx = 3
        else
            error("iLoc out of range for order 3.")
        end
        # Compute φ(ξ) = ∏_{j≠idx} (ξ - nodes_std[j]) / ∏_{j≠idx} (nodes_std[idx] - nodes_std[j])
        numerator = 1.0
        for j in 1:4
            if j != idx
                numerator *= (xi - nodes_std[j])
            end
        end
        denom = 1.0
        for j in 1:4
            if j != idx
                denom *= (nodes_std[idx] - nodes_std[j])
            end
        end
        phiVal = numerator / denom
        # Compute derivative dφ/dξ:
        dphi = 0.0
        for k in 1:4
            if k != idx
                prod_term = 1.0
                for j in 1:4
                    if j != idx && j != k
                        prod_term *= (xi - nodes_std[j])
                    end
                end
                dphi += prod_term
            end
        end
        dphi_dxi = dphi / denom
        return (phiVal, dphi_dxi)
    else
        error("Order must be 1, 2, or 3.")
    end
end

