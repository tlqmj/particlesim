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
function step!(
    particles::AbstractVector{Particle{N,T}},
    forcefield,
    dt::Real) where {N,T}

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
function sim!(
    particle::Particle{N,T},
    forcefield,
    arr::AbstractVector{Particle{N,T}}=Array{Particle{N,T},1}(),
    dt::Real=0.1,
    nsteps::Integer=1000) where {N, T<:AbstractFloat}

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
function sim!(
    particles::AbstractVector{Particle{N,T}},
    forcefield,
    mat::Matrix{Particle{N,T}}=Matrix{Particle{N,T}}(),
    dt::Real=0.1,
    nsteps::Integer=1000) where {N, T<:AbstractFloat}

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

function sim(
    particles::AbstractVector{Particle{N,T}},
    forcefield,
    dt::Number,
    fps=24::Number;
    markerscale=10) where {N, T<:AbstractFloat}

  for particle in particles
    particle.force = zeros(SVector{N,T})
  end

  for i=1:size(particles,1), j=(i+1):size(particles,1)
    forces = forcefield(particles[i], particles[j])
    particles[i].force += forces[1]
    particles[j].force += forces[2]
  end

  frame_duration = 1.0/fps
  @info frame_duration
  markersize = markerscale/maximum(p.mass for p in particles)

  scene = Scene(backgroundcolor = :black)
  scatter!(scene, [p.position for p in particles],
    glowwidth = markerscale, glowcolor = (:white, 0.1), color = :white,
    markersize = [p.mass*markersize for p in particles],
    show_axis = false,
    limits = FRect(-10, -10, 20, 20))
  display(scene)


  Δt_sim = 0
  t1 = time()

  while true
    if !scene.events.window_open.val
      return
    end

    step!(particles, forcefield, dt)

    Δt_sim += dt
    if Δt_sim > frame_duration
      push!(scene[end][1], [p.position for p in particles])

      t2 = time()
      diff = frame_duration - (t2 - t1)

      if diff > 0.0
        @info diff
        sleep(diff)
      else
        yield()
      end

      Δt_sim = 0
      t1 = t2
    end

  end
end
