using StaticArrays, RecipesBase, Plots
gr(show = :juno)

include("particles.jl")
include("simulation.jl")
include("forcefields.jl")

particles = [
  Particle(
    6e24,
    r .* (cos(θ), sin(θ)),
    5e-7 * r^0.99 .* (-sin(θ), cos(θ)))
  for θ = 0:0.1:2π for r = 1e6:1e6:1e8]

for i in 1:10
  for j in 1:1
    step!(particles, gravitationalattraction, 100)
  end
  println(i)
  display(scatter(
    particles,
    aspectratio=1,
    xlims=(-2*10^8, 2*10^8),
    ylims=(-2*10^8, 2*10^8),
    legend=false))
end

scatter(
  particles,
  aspectratio=1,
  xlims=(-2*10^8, 2*10^8),
  ylims=(-2*10^8, 2*10^8),
  legend=false)
  
