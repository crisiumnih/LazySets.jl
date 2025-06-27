function is_empty(E::Ellipsotope)
    c, G, p, I, A, b = E.center, E.generators, E.p_norm, E.index_set, E.A, E.b
    m = ngens(E)

    if isempty(A)
        return false, 0.0 # Not empty, cost = 0
    end

    # Use some solver or write helper fucntion
    

    return
end
