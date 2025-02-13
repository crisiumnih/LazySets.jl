using LazySets, Plots, LinearAlgebra

function minkow_intersection(R::Hyperrectangle{N}, E::Ellipsoid{N}) where {N}
    R_c = R.center
    E_c = E.center
    K = E_c - R_c

    M = (E.shape_matrix)
    L = zeros(N, 2)
    for i in 1:2
        axis = zeros(N, 2)
        axis[i] = 1.0
        L[i] = sqrt(dot(axis', inv(M) * axis))
    end
    
    xi = K
    
    if all(abs.(xi) .<= R.radius .+ L)
        s = [xi[1] >= 0 ? 1 : -1, xi[2] >= 0 ? 1 : -1]
        P = s .* R.radius
        Delta = K - P
        MDelta = M * Delta

        for i in 1:2
            if s[i] * MDelta[i] <= zero(N)
                return true
            end
        end

        quadratic = dot(Delta, M * Delta)
        return quadratic <= one(N)
    else
        return false  
    end
end

H = Hyperrectangle([-2.0, -1.0], [2.0, 1.0])
Q = [1.0 0.5; 0.5 2.0]
F = Ellipsoid([0.0, 1.0], Q)

S = H ⊕ F
neg_E = LinearMap(-I, F)
O = H ⊕ neg_E

intersects = minkow_intersection(H, F)
F_poly = overapproximate(F, 0.1)
println("minkowski: $intersects")
inter2 = !is_intersection_empty(H, F_poly)
println("is intersect inbuilt: $inter2")

plot(H, aspect_ratio=:equal, label="Hyperrectangle (Rectangle)", linecolor=:blue, fillalpha=0.2)
plot!(F, label="Ellipsoid", linecolor=:red, fillalpha=0.2)
plot!(O, label="Minkowski Sum (H ⊕ -E)", linecolor=:green, fillalpha=0.2)
plot!(S, label="Minkowski Sum (H ⊕ E)", linecolor=:blue, fillalpha=0.2)

title!("2D Hyperrectangle and Ellipsoid")
xlabel!("x-axis")
ylabel!("y-axis")