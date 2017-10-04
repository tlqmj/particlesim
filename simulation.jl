"""
Update position, velocity and force of particles affected by a
force 'forcefield(particle)' after a time interval 'dt'
"""
function step!(
    particle::Particle,
    forcefield,
    dt::Real)

  particle.position += particle.velocity*dt
  particle.velocity += (particle.force/particle.mass)*dt
  particle.force = forcefield(particle)

  return particle
end

"""
Update position, velocity and force of particles affected by a
force 'forcefield(particle[i], particle[j])' after a time interval 'dt'
"""
function step!{N,T}(
    particles::AbstractVector{Particle{N,T}},
    forcefield,
    dt::Real)

  for particle in particles
    particle.position += particle.velocity*dt
    particle.velocity += (particle.force/particle.mass)*dt
    particle.force     = zeros(SVector{N,T})
  end

  for i=1:size(particles)[1], j=i+1:size(particles)[1]
    forces = forcefield(particles[i], particles[j])
    particles[i].force += forces[1]
    particles[j].force += forces[2]
  end

  return particles
end

"""
Run a simulation of a single particle and save it to arr
"""
function sim!{N, T<:AbstractFloat}(
    particle::Particle{N,T},
    forcefield,
    arr::AbstractVector{Particle{N,T}}=Array{Particle{N,T},1}(),
    dt::Real=0.1,
    nsteps::Integer=1000)

  particle.force = forcefield(particle)
  push!(arr, deepcopy(particle))

  for i in 1:nsteps
    step!(particle, forcefield, dt)
    push!(arr, deepcopy(particle))
  end

  return arr
end

"""
Run a simulation and save it to mat
"""
function sim!{N, T<:AbstractFloat}(
    particles::AbstractVector{Particle{N,T}},
    forcefield,
    mat::Matrix{Particle{N,T}}=Matrix{Particle{N,T}}(),
    dt::Real=0.1,
    nsteps::Integer=1000)

  for particle in particles
    particle.force = zeros(SVector{N,T})
  end

  for i=1:size(particles), j=i+1:size(particles)
    forces = forcefield(particles[i], particles[j])
    particles[i].force += forces[1]
    particles[j].force += forces[2]
  end

  if size(mat) == (0,0)
    mat = deepcopy(particles)
  else
    mat = vcat(mat, deepcopy(particles))
  end

  for i in 1:nsteps
    step!(particles, forcefield, dt)
    mat = vcat(mat, deepcopy(particles))
  end

  return mat
end
