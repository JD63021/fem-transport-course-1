# mesh5_gmsh.jl
# This function reads a MATLAB-format Gmsh file for 1D meshes and parses
# either linear (order 1), quadratic (order 2), or cubic (order 3) elements.
#
# It expects the file to define:
#   - msh.POS : node coordinates (Nx2 or Nx3; only the first column is used)
#   - msh.LINES (for order 1), msh.LINES3 (for order 2), or msh.LINES4 (for order 3)
#   - msh.PNT (optional): boundary points in the form [nodeID, boundaryFlag]
#
# Usage:
#   nodeInfo, elemInfo, boundaryInfo = mesh5_gmsh("your_mesh_file.m", order)

function mesh5_gmsh(gmshFile::String, order::Int)
    # Read all lines from the file.
    all_lines = readlines(gmshFile)

    # Helper function: Extract a block of text following a marker (e.g., "msh.POS = [")
    function extract_block(marker::String)
        block = String[]
        start_idx = 0
        for (i, line) in enumerate(all_lines)
            if occursin(marker, line)
                start_idx = i + 1  # Data begins on the next line.
                break
            end
        end
        if start_idx == 0
            return nothing  # Marker not found.
        end
        for j in start_idx:length(all_lines)
            line = strip(all_lines[j])
            # Stop reading when the closing bracket is found.
            if occursin("];", line)
                line_clean = replace(line, "];" => "")
                if !isempty(strip(line_clean))
                    push!(block, strip(line_clean))
                end
                break
            else
                push!(block, line)
            end
        end
        return block
    end

    # -------------------------
    # 1) Extract node coordinates from msh.POS.
    pos_block = extract_block("msh.POS")
    if pos_block === nothing
        error("Could not find the msh.POS block in the file.")
    end

    pos_data = []
    for l in pos_block
        l_clean = replace(l, ";" => "")
        parts = split(l_clean)
        row = parse.(Float64, parts)
        push!(pos_data, row)
    end
    pos_matrix = reduce(vcat, [reshape(row, 1, :) for row in pos_data])
    num_nodes = size(pos_matrix, 1)

    nodeInfo = Dict{Symbol, Any}()
    nodeInfo[:x] = pos_matrix[:, 1]       # Only the x-coordinate is needed.
    nodeInfo[:id] = collect(1:num_nodes)    # Node IDs (1, 2, ..., num_nodes)

    # -------------------------
    # 2) Determine the marker for elements based on 'order'
    element_marker = ""
    if order == 1
        element_marker = "msh.LINES"
    elseif order == 2
        element_marker = "msh.LINES3"
    elseif order == 3
        element_marker = "msh.LINES4"
    else
        error("Order must be 1, 2, or 3 for a 1D Lagrange element.")
    end

    # Extract the element block.
    lines_block = extract_block(element_marker)
    if lines_block === nothing
        error("Expected $element_marker for an order=$order 1D mesh.")
    end

    # -------------------------
    # 3) Process element data.
    lines_data = []
    for l in lines_block
        l_clean = replace(l, ";" => "")
        parts = split(l_clean)
        row = parse.(Int, parts)
        push!(lines_data, row)
    end
    lines_matrix = reduce(vcat, [reshape(row, 1, :) for row in lines_data])
    
    # For order 1, take the first 2 columns; for order 2, first 3 columns; for order 3, first 4 columns.
    elem_conn = nothing
    if order == 1
        elem_conn = lines_matrix[:, 1:2]
    elseif order == 2
        elem_conn = lines_matrix[:, 1:3]
    elseif order == 3
        elem_conn = lines_matrix[:, 1:4]
    end
    elemInfo = Dict{Symbol, Any}()
    elemInfo[:elements] = elem_conn

    # -------------------------
    # 4) Process boundary points from msh.PNT if available.
    boundaryInfo = Dict{Symbol, Any}()
    boundaryInfo[:allBoundaryNodes] = Int[]
    pnt_block = extract_block("msh.PNT")  # May be nothing if not present.
    if pnt_block !== nothing
        pnt_data = []
        for l in pnt_block
            l_clean = replace(l, ";" => "")
            parts = split(l_clean)
            row = parse.(Int, parts)
            push!(pnt_data, row)
        end
        pnt_matrix = reduce(vcat, [reshape(row, 1, :) for row in pnt_data])
        # Each row of pnt_matrix is assumed to be [nodeID, boundaryFlag]
        flags = unique(pnt_matrix[:, 2])
        for flag in flags
            mask = pnt_matrix[:, 2] .== flag
            theseNodes = pnt_matrix[mask, 1]
            fieldName = Symbol("flag_$(flag)")
            boundaryInfo[fieldName] = theseNodes
            append!(boundaryInfo[:allBoundaryNodes], theseNodes)
        end
        boundaryInfo[:allBoundaryNodes] = unique(boundaryInfo[:allBoundaryNodes])
    else
        @warn "No msh.PNT found => no boundary flags read."
    end

    return nodeInfo, elemInfo, boundaryInfo
end
