#=
  todo:
  Ver qué onda el tema de los recipes de plots
  Ver qué onda el tema de los módulos
=#

using StaticArrays, RecipesBase, Plots
#import Plots.plot

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
type Particle{N, T<:AbstractFloat}
    mass::T
    position::SVector{N, T}
    velocity::SVector{N, T}
    force::SVector{N, T}
end

function Particle{N, T<:AbstractFloat}(mass::T, position::SVector{N,T}, velocity::SVector{N,T})
  Particle(mass, position, velocity, zeros(SVector{N,T}))
end

function Particle{R, S, T, N}(mass::R, position::SVector{N,S}, velocity::SVector{N,T})
  U = promote_type(R, S, T)
  if !(U <: AbstractFloat) U = typeof(1.0) end
  Particle(convert(U, mass), convert(SVector{N,U}, position), convert(SVector{N,U}, velocity))
end

function Particle(mass, position, velocity)
  Particle(mass, SVector(position...), SVector(velocity...))
end

plottrayectory{N,T}(arr::Array{Particle{N,T},1}) = plot(map(p->convert(Tuple, p.position), arr))

"""
Update position, velocity and force of particle
affected by a force 'forcefunc(particle)' after a time interval 'dt'
"""
function step!(particle::Particle, forcefunc, dt::Real)
    particle.position += particle.velocity*dt
    particle.velocity += (particle.force/particle.mass)*dt
    particle.force = forcefunc(particle)
    return particle
end

"""
Run a simulation of a single particle and save it to arr
"""
function sim!{N, T<:AbstractFloat}(particle::Particle{N,T}, forcefunc, arr::Array{Particle{N,T},1}, dt::Real=0.1, nsteps::Integer=1000)
    particle.force = forcefunc(particle)
    push!(arr, deepcopy(particle))
    for i in 1:nsteps
      step!(particle, forcefunc, dt)
      push!(arr, deepcopy(particle))
    end
    return arr
end

"""
Returns tuple of forces felt by p₁ and p₂ because of their gravitational interaction
"""
function gravitationalattraction{N}(p₁::Particle{N}, p₂::Particle{N})
    d₂₁ = p₂.position - p₁.position
    f₁  = 6.67408E-11*p₁.mass*p₂.mass*d₂₁/dot(d₂₁,d₂₁)^(3/2)
    return f₁, -f₁
end


#=
earth = Particle(5.9721986e+24, (0,0,0), (0,0,0))
p = Particle(200, (6371008 + 3000000, 0, 0), (0, 6520.2, 0))
f(particle) = gravitationalattraction(particle, earth)[1]
arr = Array{typeof(p),1}()
sim!(p, f, arr, 1, 10000)
plottrayectory(arr)
plot!(Plots.partialcircle(0, 2π, 100, 6371008))
=#
