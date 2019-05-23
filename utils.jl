function center_of_mass(particles::AbstractVector{Particle{N,T}}) where {N,T}
    sum(p.mass*p.position for p in particles)./sum(p.mass for p in particles)
end

function center_of_mass_velocity(particles::AbstractVector{Particle{N,T}}) where {N,T}
    sum(p.mass*p.velocity for p in particles)./sum(p.mass for p in particles)
end

function avg_distance_from_cm(particles::AbstractVector{Particle{N,T}}) where {N,T}
    cm = center_of_mass(particles)
    mean((norm(p.position - cm) for p in particles))
end

function normalize_cm!(particles::AbstractVector{Particle{N,T}}) where {N,T}
    cm = center_of_mass(particles)
    v_cm = center_of_mass_velocity(particles)

    for p in particles
      p.position = p.position - cm
      p.velocity = p.velocity - v_cm
    end

    return
end

function load_particles(path="particles.bson")
  @load path particles
  return particles = [particles...]
end

function prune!(particles::AbstractVector{Particle{N,T}}, diameter) where {N,T}
    normalize_cm!(particles)
    keep = [norm(p.position) < diameter for p in particles]
    particles = deepcopy(particles[keep])
    normalize_cm!(particles)
    return particles
end

function spherical_to_cartesian(r, θ, ϕ)
    return (r*cos(θ)*sin(ϕ), r*sin(θ)*sin(ϕ), r*cos(ϕ))
end

function plummer_sphere(n, m, a)
    particles = [sample_plummer_sphere(n,m,a) for _ in 1:n]

    normalize_cm!(particles)
    return particles
end

function sample_plummer_sphere(n, m, a)
    # See https://stackoverflow.com/questions/28434765/how-to-randomly-distribute-n-masses-such-that-they-follow-a-plummer-density-dis
    radius = a*(rand()^(-2/3) - 1)^(-1/2)
    position = spherical_to_cartesian(radius, sample_uniform_angular_distribution()...)

    x = 0.0
    y = 0.1
    while y > x^2*(1.0-x^2)^3.5
        x = rand()
        y = 0.1*rand()
    end
    speed = x * sqrt(2*G*n*m) * (radius^2 + a^2)^(-1/4)
    velocity = spherical_to_cartesian(speed, sample_uniform_angular_distribution()...)

    return Particle(m, position, velocity)
end

function sample_uniform_angular_distribution()
    return (acos(2*rand() - 1), 2π*rand())
end
