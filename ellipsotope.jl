using LazySets
using Plots
using LinearAlgebra

# Revised struct definition with explicit type parameter usage
struct Ellipsotope{N<:AbstractFloat, VN<:AbstractVector{N}, MN<:AbstractMatrix{N}, I<:Vector{Vector{Int}}} <: AbstractCentrallySymmetric{N}
    center::VN
    generators::MN
    index_set::I
    p::Float64
    A::Union{MN, Nothing}
    b::Union{VN, Nothing}
end

# Constructor without constraints (fixed type parameter binding)
function Ellipsotope(c::VN, G::MN, index_set::I, p::Float64) where {N<:AbstractFloat, VN<:AbstractVector{N}, MN<:AbstractMatrix{N}, I<:Vector{Vector{Int}}}
    return Ellipsotope{N, VN, MN, I}(c, G, index_set, p, nothing, nothing)
end

# Membership check (same as before)
function ∈(x::AbstractVector, E::Ellipsotope)
    β = E.generators \ (x - E.center)
    if E.A !== nothing && E.b !== nothing && !isapprox(E.A * β, E.b)
        return false
    end
    for J in E.index_set
        norm(β[J], E.p) > 1 + eps() && return false
    end
    return true
end

# Linear map (same as before)
function linear_map(M::AbstractMatrix, E::Ellipsotope)
    c_new = M * E.center
    G_new = M * E.generators
    return Ellipsotope(c_new, G_new, E.index_set, E.p, E.A, E.b)
end

# Support vector (same as before)
function σ(d::AbstractVector, E::Ellipsotope)
    c = E.center
    β_opt = similar(d)
    for J in E.index_set
        if E.p == Inf
            β_opt[J] .= sign.(d[J])
        else
            β_opt[J] .= d[J] / norm(d[J], E.p)
        end
    end
    return c + E.generators * β_opt
end

# Example usage (now works)



function plot_ellipsotope(E::Ellipsotope; npoints=100)
    # For 2D Ellipsotopes
    if size(E.generators, 1) == 2
        # Create a unit circle/square of points based on the p-norm
        θ = range(0, 2π, length=npoints)
        
        # Generate boundary points based on p-norm
        boundary_points = zeros(2, npoints)
        
        for i in 1:npoints
            if E.p == Inf
                # For infinity norm (box)
                x = cos(θ[i]) > 0 ? 1.0 : -1.0
                y = sin(θ[i]) > 0 ? 1.0 : -1.0
                point = [x, y]
            elseif E.p == 1
                # For 1-norm (diamond)
                t = θ[i]
                r = 1.0 / (abs(cos(t)) + abs(sin(t)))
                point = r * [cos(t), sin(t)]
            elseif E.p == 2
                # For 2-norm (circle)
                point = [cos(θ[i]), sin(θ[i])]
            else
                # For general p-norm
                t = θ[i]
                # Find r such that (r*cos(t))^p + (r*sin(t))^p = 1
                cosₜ, sinₜ = cos(t), sin(t)
                r = (abs(cosₜ)^E.p + abs(sinₜ)^E.p)^(-1/E.p)
                point = r * [cosₜ, sinₜ]
            end
            
            # Transform point according to the index set structure
            β = zeros(size(E.generators, 2))
            for (idx, J) in enumerate(E.index_set)
                if length(J) == 1
                    β[J[1]] = point[idx]
                else
                    for j in J
                        β[j] = point[idx]  # This is simplified; more complex mappings may be needed
                    end
                end
            end
            
            # Apply generators and add center
            boundary_points[:, i] = E.center + E.generators * β
        end
        
        # Create the plot
        plt = plot(boundary_points[1, :], boundary_points[2, :], 
                  label="Ellipsotope (p=$(E.p))", 
                  aspect_ratio=:equal, 
                  fill=true, fillalpha=0.3,
                  legend=:outertopright,
                  title="Ellipsotope Plot")
        
        # Plot the center point
        scatter!([E.center[1]], [E.center[2]], label="Center", markersize=4)
        
        return plt
    else
        error("Currently only 2D Ellipsotopes can be plotted with this function")
    end
end
E = Ellipsotope([0.0, 0.0], [1.0 0.0; 0.0 1.0], [[1], [2]], 4.0)

plt = plot_ellipsotope(E)
savefig(plt, "ellipsotope.png")
display(plt)