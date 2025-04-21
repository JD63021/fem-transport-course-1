using Plots
using LinearAlgebra
using Statistics  # for mean()
using Printf      # for @sprintf

# --------------------------
# Define a simple linear regression function.
function linreg(x::Vector{Float64}, y::Vector{Float64})
    n = length(x)
    xmean = mean(x)
    ymean = mean(y)
    slope = sum((x .- xmean) .* (y .- ymean)) / sum((x .- xmean).^2)
    intercept = ymean - slope * xmean
    return slope, intercept
end

# Include required function files.
include("mesh5_gmsh.jl")
include("main1_SS.jl")
include("precomputeShapeFunctions1D.jl")
# (Include any other necessary files.)

# ===========================
# Test case configuration function for MMS
# ===========================
function testCase_config()
    testCase = Dict{Symbol, Any}()
    # The mesh file will be updated in the loop.
    testCase[:gmshVelFile] = ""
    # For MMS, we use Dirichlet BCs at both endpoints.
    testCase[:boundaryFlags] = Dict(:wall => [2, 3])
    # Set the element order here (1 for linear, 2 for quadratic, 3 for cubic).
    testCase[:order] = 2   # Modify as needed.
    testCase[:inletProfile] = x -> 0   # not used
    testCase[:inletProfileSS] = 0      # not used (Dirichlet 0)
    testCase[:constantsource] = 1      # not used here
    return testCase
end
# Set mesh file prefix.
# For example, use "oneDla" for a series of linear meshes,
# or "oneDca" for cubic meshes. Adjust testCase[:order] accordingly.
meshPrefix = "oneDqa"   # Modify as needed.

# Manufactured forcing function: f(x) = 100*sin(10*x)
f_func(x) = 100 * sin(10 * x)

# -------------------------
# Start timer
start_time = time()

# User parameter
D = 1.0

# Prepare arrays to store mesh size and L2 error norms.
h_arr = Float64[]
error_arr = Float64[]

# Load base test case configuration.
testCase = testCase_config()



# Loop over mesh files: e.g., oneDla1.m, oneDla2.m, ..., oneDla4.m.
for i in 1:4
    fileName = "$(meshPrefix)$(i).m"
    testCase[:gmshVelFile] = fileName
    println("\nProcessing mesh file: $fileName")
    
    # Generate Mesh & Setup.
    nodeInfo, elemInfo, boundaryInfo = mesh5_gmsh(fileName, testCase[:order])
    Nxy = length(nodeInfo[:x])
    
    # Compute element size h for an equally spaced mesh.
    h = (maximum(nodeInfo[:x]) - minimum(nodeInfo[:x])) / (Nxy - 1)
    println("Element size h = $h")
    
    # Set initial guess.
    U0 = zeros(Nxy)
    
    # Manufactured solution: u(x) = sin(10*x) evaluated at the nodes.
    U_ms = sin.(10 .* nodeInfo[:x])
    
    # Solve the stationary problem using the forcing function.
    U_sol, x_sol, nodeInfo, elemInfo, boundaryInfo = main1_SS(U0, D, nodeInfo, elemInfo, boundaryInfo,
        1, testCase[:boundaryFlags], testCase[:inletProfileSS], testCase[:order], f_func)
    
    # Compute integrated L2 error norm using 5-point Gauss quadrature.
    Nshape, dNshape, gp, wq = precomputeShapeFunctions1D(testCase[:order], 5)
    nLoc = testCase[:order] + 1
    numElems = size(elemInfo[:elements], 1)
    L2_error_sq = 0.0
    
    for e in 1:numElems
        nodes = elemInfo[:elements][e, :]
        xcoords = nodeInfo[:x][nodes]
        Uh_loc = U_sol[nodes]
        for ig in 1:length(wq)
            # Map quadrature point from reference [0,1] to physical element.
            x_gp = dot(Nshape[:, ig], xcoords)
            # Reconstruct FEM solution at quadrature point.
            u_h = dot(Nshape[:, ig], Uh_loc)
            # Exact solution at x_gp.
            u_ex = sin(10 * x_gp)
            # Compute Jacobian: dx/dξ = Σ xcoords[i] * dNshape[i,ig]
            dx_dxi = sum(xcoords .* dNshape[:, ig])
            J = abs(dx_dxi)
            L2_error_sq += wq[ig] * (u_h - u_ex)^2 * J
        end
    end
    L2_error = sqrt(L2_error_sq)
    println("Integrated L2 norm error = $(L2_error)")
    
    push!(h_arr, h)
    push!(error_arr, L2_error)
end

# -------------------------
# Plot L2 error norm vs. element size h on a log-log scatter plot.
p_err = scatter(h_arr, error_arr, xscale=:log10, yscale=:log10,
             xlabel="Element size h", ylabel="L2 Error Norm",
             title="L2 Error vs. Element Size (Manufactured Solution)", legend=false)

# Compute slope via linear regression using only the last 4 values.
logh = log10.(h_arr)
logerr = log10.(error_arr)
subset = (length(h_arr)-3):length(h_arr)
slope, intercept = linreg(logh[subset], logerr[subset])
annot_text = @sprintf("Slope = %.2f", slope)
annotate!(p_err, minimum(h_arr), maximum(error_arr), text(annot_text, :left, 12))
display(p_err)

elapsed_time = time() - start_time
println("Total elapsed time: $(elapsed_time) seconds")
