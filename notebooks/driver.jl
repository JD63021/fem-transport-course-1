# driver.jl
# =========================================================================
# This driver sets up the mesh, initial conditions, and calls the stationary
# solver (main1_SS) for the 1D Poisson/heat problem, then visualizes the 
# final solution and (optionally) the shape functions.
#
# IMPORTANT: This file assumes that the following functions are available:
#   - mesh5_gmsh(gmshFile::String, order::Int)
#   - main1_SS(U, D, nodeInfo, elemInfo, boundaryInfo, flag, boundaryFlags,
#              inletProfileSS, order, constantsource)
#   - plot1D_FEM_BasisAndLegendre(nodeInfo, elemInfo, order, U)
# =========================================================================

using Plots

include("mesh5_gmsh.jl")
include("main1_SS.jl")
include("plot1D_FEM_BasisAndLegendre.jl")


# ===========================
# Test case configuration function
# ===========================
function testCase_config()
    # This function sets up the test case configuration for the 1D Poisson problem.
    # It returns a dictionary with all necessary parameters.
    testCase = Dict{Symbol, Any}()
    testCase[:gmshVelFile]  = "oneD6pt.m"   # gmsh file for velocity grid (Q2)
    # testCase[:gmshPresFile] = "quadp1.m"  # gmsh file for pressure grid (Q1) -- not used here

    testCase[:excludedVelFlags] = []             # velocity boundary flags to exclude
    testCase[:pBoundFlags]      = []              # no pressure BCs in this case
    testCase[:pBoundVals]       = []
    # Define the boundary flags for enforcing Dirichlet conditions:
    testCase[:boundaryFlags] = Dict(:inlet => [2], :wall => [3])
    testCase[:order] = 1
    testCase[:constantsource] = 1

    # Specify the inlet velocity profile as an anonymous function.
    # For this test, it simply returns 1.
    testCase[:inletProfile]   = x -> 1
    testCase[:inletProfileSS] = 1
    return testCase
end

# -------------------------
# Start timer
start_time = time()

# -------------------------
# User parameters
D = 1.0    # Diffusivity
t1 = 0.0   # Not used if purely stationary

# -------------------------
# Load test case configuration
testCase = testCase_config()  # Now the function is defined.

# -------------------------
# Generate Mesh & Setup
nodeInfo, elemInfo, boundaryInfo = mesh5_gmsh(testCase[:gmshVelFile], testCase[:order])
Nxy = length(nodeInfo[:x])  # number of nodes

# -------------------------
# Set initial condition (if U does not already exist, here we initialize to zeros)
U = zeros(Nxy)

# -------------------------
# Stationary solver
# Call the main1_SS function which is assumed to return:
#   [U, x, nodeInfo, elemInfo, boundaryInfo]
U, x, nodeInfo, elemInfo, boundaryInfo = main1_SS(U, D, nodeInfo, elemInfo, boundaryInfo,
    1, testCase[:boundaryFlags], testCase[:inletProfileSS], testCase[:order], testCase[:constantsource])

# -------------------------
# Visualization of final solution
p1 = plot(x, U, seriestype=:scatter, markersize=8, label="Final Temperature", grid=true,
          title="Final Temperature (1D Poisson)", xlabel="x", ylabel="T(x)")
display(p1)


# -------------------------
# OPTIONAL: Plot the shape functions.
plot1D_FEM_BasisAndLegendre(nodeInfo, elemInfo, testCase[:order], U)

# -------------------------
# Print elapsed time
elapsed_time = time() - start_time
println("Elapsed time: $(elapsed_time) seconds")
