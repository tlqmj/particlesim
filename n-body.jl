using StaticArrays, LinearAlgebra, Makie, BenchmarkTools, AbstractPlotting, Statistics
using BSON: @save, @load

include("particles.jl")
include("simulation.jl")
include("forcefields.jl")
include("utils.jl")

particles = plummer_sphere(500, 1.0, 10)
particles = prune!(particles, 4avg_distance_from_cm(particles))

particles2 = plummer_sphere(500, 1.0, 10)
particles2 = prune!(particles2, 4avg_distance_from_cm(particles2))

for p in particles2
  p.position = p.position + SVector(80.0, 20.0, 50.0)
  p.velocity = p.velocity - 0.0001*SVector(8.0, 2.0, 0.0)
end

append!(particles, particles2)
normalize_cm!(particles)

sim(particles, gravitationalattraction, dt=0.1,
    fps=10, speedup=5, limit_fps=false, save_video=false,
    markersize=0.5, glowwidth=20)

scene = Scene(backgroundcolor = :black)
scatter!(scene, [p.position for p in particles],
  glowwidth = 1, glowcolor = (:white, 0.1), color = :white,
  markersize=[25*p.mass/10000000 for p in particles],
  show_axis = true,
  limits = FRect3D([-1500, -1500, -1500], [3000, 3000, 3000]))


@save "./particles.bson" particles
