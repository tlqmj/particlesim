using StaticArrays, LinearAlgebra, Makie, BenchmarkTools, AbstractPlotting
using BSON: @save, @load

include("particles.jl")
include("simulation.jl")
include("forcefields.jl")

#=
particles = [
  Particle(
    10000000*rand(),
    100*r .* (cos(θ), sin(θ)),
    0.05r .* (-sin(θ) + 0.1*rand(),cos(θ) + 0.1*rand()) .- 0.1*rand().*(cos(θ), sin(θ)))
  for θ = 0:0.1:2π for r = 0.5:0.1:5]
=#
@load "particles.bson" particles

scene = Scene(backgroundcolor = :black, limits = FRect(-5, -5, 5, 5))
scatter!(scene, [p.position for p in particles],
  glowwidth = 1, glowcolor = (:white, 0.1), color = :white,
  markersize=[25*p.mass/10000000 for p in particles],
  show_axis = true,
  limits = FRect(-5, -5, 5, 5))
scene
sim(particles, gravitationalattraction, 0.1, 24)

@save "./particles.bson" particles
