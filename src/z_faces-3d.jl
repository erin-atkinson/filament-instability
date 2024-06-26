# Want 3/4 of cells in the boundary layer
g(s) = (s + log(2*cosh(s)))/2

@inline function get_z_faces(simulation_parameters)
    # Want 3/4 of cells in the boundary layer
    sp = simulation_parameters
    α₁ = (sp.Lz - sp.H) / (sp.Nz/4-1)
    α₂ = 4*sp.H / (3sp.Nz)
    z_faces_1(n) = -sp.Lz + (n-1) * α₁ + 20*(α₂ - α₁)* g((n - sp.Nz / 4) / 20)
    # This will usually be off a bit, so we rescale for consistency (129 is 0, 1 is -5)
    return z_faces(n) = -sp.Lz * (z_faces_1(n) - z_faces_1(sp.Nz+1)) / (z_faces_1(1)- z_faces_1(sp.Nz+1))
end