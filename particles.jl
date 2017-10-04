"""
Type for representing particles in N dimensional space.

Can be created specifying all the parameters using SVectors of the same type through the default constructor:
  Particle(mass::T, position::SVector{N,T}, velocity::SVector{N,T}, force::SVector{N,T})
Or skipping the initial force in order to initialize it at rest:
  Particle(mass::T, position::SVector{N,T}, velocity::SVector{N,T})
Or mixing SVectors of different types (with no force):
  Particle(mass::R, position::SVector{N,S}, velocity::SVector{N,T})
Or, generally, specifying the position and velocity as tuples or arrays of any type:
  Particle(mass, position, velocity)
"""
mutable struct Particle{N, T<:AbstractFloat}
    mass::T
    position::SVector{N, T}
    velocity::SVector{N, T}
    force::SVector{N, T}
end

function Particle{N, T<:AbstractFloat}(
    mass::T,
    position::SVector{N,T},
    velocity::SVector{N,T})

  Particle(mass, position, velocity, zeros(SVector{N,T}))
end

function Particle{R, S, T, N}(
    mass::R,
    position::SVector{N,S},
    velocity::SVector{N,T})

  U = promote_type(R, S, T)
  if !(U <: AbstractFloat) U = typeof(1.0) end
  Particle(convert(U, mass), convert(SVector{N,U}, position), convert(SVector{N,U}, velocity))
end

function Particle(
    mass,
    position,
    velocity)

  Particle(mass, SVector(position...), SVector(velocity...))
end

"""
Plot recipe for particle vectors (trayectories or collections of particles)
"""
@recipe function f{N,T}(
    arr::AbstractVector{Particle{N,T}},
    args...)

  tuple( map(p->convert(Tuple, p.position), arr), args...)
end
