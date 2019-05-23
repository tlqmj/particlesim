const G = 4.302e-3 # G is in pc*M⊙⁻¹ (km/s)²

"""
Returns tuple of forces felt by p₁ and p₂ due to their gravitational interaction
"""
function gravitationalattraction(
    p₁::Particle{N},
    p₂::Particle{N}) where {N}

  d₂₁ = p₂.position - p₁.position
  f₁  = G*p₁.mass*p₂.mass*d₂₁/dot(d₂₁,d₂₁)^(3/2) # G is in pc*M⊙⁻¹ (km/s)²

  return f₁, -f₁
end
