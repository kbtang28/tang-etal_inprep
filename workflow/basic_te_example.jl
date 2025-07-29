# activate environment
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using CairoMakie, Random
import ColorSchemes: seaborn_colorblind
import StatsBase: cor

Random.seed!(11)

function gen_X_and_Y(c, nt)
    X = repeat([-1; 1], outer=round(Int, nt/2))
    Y = Float64[]
    for t in 2:length(X)
        if rand() < (1+c)/2
            push!(Y, X[t-1])
        else
            push!(Y, -X[t-1])
        end
    end
    return X[2:end], Y
end

function te(c)
    if c == 1.0
        return 0.0
    else
        return 0.5*((1+c)*log(1+c) + (1-c)*log(1-c)) - 0.5*((1+c^2)*log(1+c^2) + (1-c^2)*log(1-c^2))
    end
end

inch=96
pt = 4/3

fig = Figure(size=(0inch, 2.75inch), fontsize=8pt)
ax1 = Axis(fig[1,3], limits=((0,51), (-1.2, 1.2)))
hidedecorations!(ax1)
X1, Y1 = gen_X_and_Y(0.0, 1e4)
stX = stem!(ax1, 1:50, X1[50:99], color=seaborn_colorblind[2], markersize=5.0, stemcolor=seaborn_colorblind[2], stemwidth=2.0, trunkwidth=0.0)
stY = stem!(ax1, 1:50, 0.75*Y1[50:99], color=seaborn_colorblind[1], markersize=8.0, stemcolor=seaborn_colorblind[1], stemwidth=3.5, trunkwidth=0.0)
hlines!(ax1, [0.0], color=:black, linewidth=0.5)
Box(fig[1,2], color=:gray90)
Label(fig[1,2], "c = 0.0", rotation=pi/2, tellheight=false)
Label(fig[1,1], "(a)", tellheight=false, font=:bold, valign=:top)
gl = fig[1,4] = GridLayout()
te1 = te(0.0)
r1 = round(cor(X1, Y1), digits=2) ≈ 0.0 ?  0.0 : abs(cor(X1,Y1))
Label(gl[2,1], "TE(X → Y) = $(round(te1, digits=2)),\nr = $(r1)", justification=:left, tellheight=false)
Legend(gl[1,1], 
       [MarkerElement(color=seaborn_colorblind[2], marker=:circle), MarkerElement(color=seaborn_colorblind[1], marker=:circle)], ["Xₜ", "Yₜ"], 
       orientation=:horizontal, 
       patchlabelgap=0.0, 
       colgap=10.0, 
       halign=:left, 
       padding=(0.0, 8.0, 0.0, 0.0)
)

ax2 = Axis(fig[2,3], limits=((0,51), (-1.2, 1.2)))
hidedecorations!(ax2)
X2, Y2 = gen_X_and_Y(0.75, 1e4)
stem!(ax2, 1:50, X2[50:99], color=seaborn_colorblind[2], markersize=5.0, stemcolor=seaborn_colorblind[2], stemwidth=2.0, trunkwidth=0.0)
stem!(ax2, 1:50, 0.75*Y2[50:99], color=seaborn_colorblind[1], markersize=8.0, stemcolor=seaborn_colorblind[1], stemwidth=3.5, trunkwidth=0.0)
hlines!(ax2, [0.0], color=:black, linewidth=0.5)
Box(fig[2,2], color=:gray90)
Label(fig[2,2], "c = 0.75", rotation=pi/2, tellheight=false)
Label(fig[2,1], "(b)", tellheight=false, font=:bold, valign=:top)
te2 = te(0.75)
r2 = cor(X2, Y2)
Label(fig[2,4], "TE(X → Y) = $(round(te2, digits=2)),\nr = $(round(r2, digits=2))", justification=:left, tellheight=false)

ax3 = Axis(fig[3,3], limits=((0,51), (-1.2, 1.2)), xlabel="t")
hidedecorations!(ax3, label=false)
X3, Y3 = gen_X_and_Y(1.0, 1e4)
stem!(ax3, 1:50, X3[50:99], color=seaborn_colorblind[2], markersize=5.0, stemcolor=seaborn_colorblind[2], stemwidth=2.0, trunkwidth=0.0)
stem!(ax3, 1:50, 0.75*Y3[50:99], color=seaborn_colorblind[1], markersize=8.0, stemcolor=seaborn_colorblind[1], stemwidth=3.5, trunkwidth=0.0)
hlines!(ax3, [0.0], color=:black, linewidth=0.5)
Box(fig[3,2], color=:gray90)
Label(fig[3,2], "c = 1.0", rotation=pi/2, tellheight=false)
Label(fig[3,1], "(c)", tellheight=false, font=:bold, valign=:top)
te3 = te(1.0)
r3 = cor(X3, Y3)
Label(fig[3,4], "TE(X → Y) = $(round(te3, digits=2)),\nr = $(round(r3, digits=2))", justification=:left, tellheight=false)

colgap!(fig.layout, 1, 4)
colgap!(fig.layout, 2, 0)
colgap!(fig.layout, 3, 8)
colsize!(fig.layout, 3, Aspect(1, 6.0))
resize_to_layout!(fig)

save(joinpath(@__DIR__, "..", "figures", "basic_te_example.png"), fig, dpi=300)