using Distributions
using CairoMakie
using Makie.Colors

"""
    lognormal_meanone(σ)

LogNormal distribution with standard deviation σ on the log scale, and mean 1
on the raw scale. For this, the mean on the log scale is set to σ²/2.
"""
function lognormal_meanone(sigma)
    mu = -sigma^2 / 2
    return LogNormal(mu, sigma)
end
"""
    expected_percentchange(d::Distribution; n=10_000)

Estimate of  the expected percent change in rate, that is:
E(|R-1|*100) where R ~ d.
This estimate is obtained by simulating `n` rate values from `d`.

Note that E(R)=1 for a `lognormal_meanone`, so E((R-1)*100) = 0,
if we remove the absolute value of R-1.
"""
function expected_percentchange(d::Distribution; n::Int=10_000)
    rates = rand(d, n)
    mean(abs.(rates .- 1)) * 100
end

ratedist15 = lognormal_meanone(0.15);
# find 0.1 and 0.9 quantiles: 0.82, 1.2
percentile = round.(quantile(ratedist15, [0.10,0.90]), digits=2)
expected_percentchange(ratedist15, n=1_000_000) # 12.0%
# called 3 times: 11.961996341224626, 11.965410785526608, 11.96700192051297

ratedist3 = lognormal_meanone(0.3);
# find quantiles, to show on the plot: 0.65, 1.4
percentile = round.(quantile(ratedist3, [0.10,0.90]), digits=2)
expected_percentchange(ratedist3, n=1_000_000) # 23.84649, 23.82303

ratedist5 = lognormal_meanone(0.5);
# find quantiles, to show on the plot: 0.46, 1.67
percentile = round.(quantile(ratedist5, [0.10,0.90]), digits=2)
expected_percentchange(ratedist5, n=1_000_000) # 39.501648, 39.421081

ratedist = lognormal_meanone(0.7);
# find quantiles, to show on the plot: 0.32, 1.92
percentile = round.(quantile(ratedist, [0.10,0.90]), digits=2)
expected_percentchange(ratedist, n=1_000_000) # 54.65084, 54.74559

size_inches = (4.5,3) # width,height, in inches. 72 pts per inch
fig = Figure(resolution = 72 .* size_inches, fontsize = 12)
ax = Axis(fig[1, 1],
  xlabel = "rate R\nmean 1",
  #xlabel = "rate R\nmean 1, percentiles: $(percentile[1]) (10%) - $(percentile[2]) (90%)",
  ylabel = "density",
  xticks = [0,0.5,1,2,3,4,5])
xlims!(ax, 0,4)
ylims!(ax, 0,3)
hideydecorations!(ax, label=false)

dists = [ratedist15, ratedist3, ratedist5, ratedist]
colors = [Colors.HSV(19, 87, 99), Colors.HSV(44, 100, 96), Colors.HSV(217, 77, 100), Colors.HSV(265, 71, 93)]

num = 1
for dist in dists
    plot!(dist, color = colors[num])
    num += 1
end
    
elem_1 = [LineElement(color = Colors.HSV(19, 87, 99), linestyle = nothing)] #FC6722
elem_2 = [LineElement(color = Colors.HSV(44, 100, 96), linestyle = nothing)] #F5B400
elem_3 = [LineElement(color = Colors.HSV(217, 77, 100), linestyle = nothing)] #3A86FF
elem_4 = [LineElement(color = Colors.HSV(265, 71, 93), linestyle = nothing)] #8B44EE

Legend(fig[1, 1],
    [elem_1, elem_2, elem_3, elem_4],
    ["σ = 0.15", "σ = 0.3", "σ = 0.5", "σ = 0.7"],
    tellheight = false,
    tellwidth = false,
    patchsize = (15, 15), rowgap = 2,
    halign = :right, valign = :top, orientation = :vertical)

fig
save("fig_f4lognormal.pdf", fig)