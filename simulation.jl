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

function initialize_forces!(particles::AbstractVector{Particle{N,T}}, forcefield) where {N, T<:AbstractFloat}

  for particle in particles
    particle.force = zeros(SVector{N,T})
  end

  for i=1:size(particles,1), j=(i+1):size(particles,1)
    forces = forcefield(particles[i], particles[j])
    particles[i].force += forces[1]
    particles[j].force += forces[2]
  end

  return nothing
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

  initialize_forces!(particles, forcefield)

  scene = Scene(backgroundcolor = :black)
  scatter!(scene, [p.position for p in particles],
    glowwidth = glowwidth, glowcolor = (:white, glowalpha), color = :white,
    markersize = markersize,
    show_axis = true,
    transparency=true)
  display(scene)

  if save_video
    io = VideoStream(scene; framerate = fps)
  end

  frame_duration = 1.0/fps
  Δt_sim_per_frame = frame_duration*speedup
  if Δt_sim_per_frame < dt
    @warn "dt is too high or speedup to low."
  end
  Δt_sim = 0
  t1 = time()

  while true
    if !scene.events.window_open.val
      if save_video
        save(save_path, io, framerate=fps)
        @info "Video saved at $(pwd())/$(save_path)"
      end
      return
    end

    step!(particles, forcefield, dt)

    Δt_sim += dt
    if Δt_sim > Δt_sim_per_frame
      push!(scene[end][1], [p.position for p in particles])
      if save_video
        recordframe!(io)
      end

      if limit_fps
        t2 = time()
        diff = frame_duration - (t2 - t1)
        t1 = t2

        if diff > 0.0
          sleep(diff)
        end

      end

      @info center_of_mass_velocity(particles)
      Δt_sim = 0
    end

  end
end
