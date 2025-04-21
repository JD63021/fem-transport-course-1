### A Pluto.jl notebook ###
# ╔═╡ Cell1
begin
    using Pluto
end

# ╔═╡ markdown1
"""
# FEM Transport & SUPG Demo

This notebook will run your full Julia solver—all in the browser.
"""

# ╔═╡ markdown2
"""
**How to execute:**  
Click ▶️ on the next cell. It will `include("driver.jl")` and run everything.
"""

# ╔═╡ run_driver
begin
    # notebook lives in notebooks/, so go one level up to repo root
    repo_root = joinpath(@__DIR__, "..")
    include(joinpath(repo_root, "notebooks", "driver.jl"))
end
