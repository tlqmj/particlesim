using StaticArrays, LinearAlgebra, Makie, BenchmarkTools, AbstractPlotting, Statistics
using BSON: @save, @load

include("particles.jl")
include("simulation.jl")
include("forcefields.jl")
include("utils.jl")

particles = plummer_sphere(2000, 100.0, 1500)
particles = prune!(particles, 3avg_distance_from_cm(particles))

sim(particles, gravitationalattraction, dt=5, fps=10, speedup=10, markerscale=20)

scene = Scene(backgroundcolor = :black)
scatter!(scene, [p.position for p in particles],
  glowwidth = 1, glowcolor = (:white, 0.1), color = :white,
  markersize=[25*p.mass/10000000 for p in particles],
  show_axis = true,
  limits = FRect3D([-1500, -1500, -1500], [3000, 3000, 3000]))


@save "./particles.bson" particles
