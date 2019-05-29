function update_forces!(particles::AbstractVector{Particle{N,T}}, forcefield) where {N, T<:AbstractFloat}

  for particle in particles
    particle.force = zeros(SVector{N,T})
  end

  @inbounds for i=1:(size(particles, 1)-1), j=(i+1):size(particles, 1)
    forces = forcefield(particles[i], particles[j])
    particles[i].force += forces[1]
    particles[j].force += forces[2]
  end

  return nothing
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
  end

  update_forces!(particles, forcefield)

  for particle in particles
    particle.velocity += (particle.force/particle.mass)*dt
  end

  return particles
end


function first_step!(
    particles::AbstractVector{Particle{N,T}},
    forcefield,
    dt::Real) where {N,T}

  update_forces!(particles, forcefield)

  for particle in particles
    particle.velocity += (particle.force/particle.mass)*dt/2
  end

  step!(particles, forcefield, dt)
end


"""
  function advance_sim!(particles, forcefield, dt, Δt)

Advances the simulation a time t ∈ [Δt, Δt+dt).
Returns t - Δt.
"""
function advance_sim!(
    particles::AbstractVector{Particle{N,T}},
    forcefield,
    dt,
    Δt) where {N, T<:AbstractFloat}

  Δt_partial = zero(dt)

  while Δt_partial <= Δt
    step!(particles, forcefield, dt)
    Δt_partial += dt
  end

  return Δt_partial - Δt
end


function sim(
    particles::AbstractVector{Particle{N,T}},
    forcefield;
    dt,
    speedup=1.0,
    fps=24,
    limit_fps::Bool=true,
    save_video::Bool=false,
    save_path::String="simulation.mp4",
    glowalpha=0.02,
    markersize=10,
    glowwidth=markersize) where {N, T<:AbstractFloat}

  frame_duration = 1.0/fps
  Δt_sim_per_frame = frame_duration*speedup
  t_overshoot = 0.0
  2Δt_sim_per_frame < dt && error("dt is way too high. Decrease it, increase speedup or lower fps")
  Δt_sim_per_frame < dt  && @warn "dt is too high. Decrease it, increase speedup or lower fps"

  scene = Scene(backgroundcolor = :black)
  scatter!(scene, [p.position for p in particles],
    glowwidth = glowwidth, glowcolor = (:white, glowalpha), color = :white,
    markersize = markersize,
    show_axis = true,
    transparency=true)

  cameracontrols(scene).fov[] = 15
  update_cam!(scene, cameracontrols(scene))

  display(scene)

  if save_video
    io = VideoStream(scene; framerate = fps)
  end

  limiter = Limiter(frame_duration)
  first_step!(particles, forcefield, dt)

  while true
    t_overshoot = advance_sim!(particles, forcefield, dt, Δt_sim_per_frame - t_overshoot)

    scene[end][1] = [p.position for p in particles]
    if save_video
      recordframe!(io)
    end

    if limit_fps
      limiter()
    end

    if !scene.events.window_open.val
      if save_video
        save(save_path, io, framerate=fps)
        @info "Video saved at $(pwd())/$(save_path)"
      end
      return
    end

    @info angular_momentum(particles)
  end
end


mutable struct Limiter
  Δt::Float64
  t₀::Float64
end

function Limiter(Δt::Number)
  return Limiter(Float64(Δt), time())
end

function (lim::Limiter)()
  t = time()
  diff = lim.Δt - (t - lim.t₀)

  if diff > 0.0
    sleep(diff)
  else
    yield()
  end

  lim.t₀ = t
end
