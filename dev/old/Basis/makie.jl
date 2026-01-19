using AtomicKohnSham
using CairoMakie
using LaTeXStrings

m = expmesh(0,10,5;s=1.0)
b = P1IntLegendreBasis(m;ordermax = 3)
X = LinRange(0,10,10000)

fig = Figure(size = (1200, 1000), fontsize = 30)

xtickcustom = [L"r_0", L"r_1", L"r_2", L"r_3",L"r_4"]
ax = Axis(fig[1, 1],
        xticks = m.points,
        xtickformat = values -> xtickcustom,
        ylabel = L"y",
        xticklabelsize = 70,
        yticklabelsize = 70,
        xlabelsize = 70,
        ylabelsize = 60,
        limits=((0, last(m)), nothing)
    )


plots = []
labels = []

for i ∈ 1:b.size
    biX = [b(i,x)[1] for x ∈ X]
    push!(plots, lines!(ax, X, biX; linewidth = 9))
    push!(labels, LaTeXStrings.LaTeXString("Q_{$i}"))
end

Legend(
    fig[1, 2],
    plots,
    labels,
    halign = :left,
    valign = :top,
    margin = (20, 20, 20, 20),
    orientation = :vertical,
    labelsize = 70,
)

fig
save("./dev/Basis/illustration_basis.png",fig)
