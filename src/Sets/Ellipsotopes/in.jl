function ∈(x::AbstractVector, E::Ellipsotope)
    c, G, p, I, A, b = E.center, E.generators, E.p_norm, E.index_set, E.A, E.b
    n, m = size(G)

    # define new constraints and solve emptiness

end