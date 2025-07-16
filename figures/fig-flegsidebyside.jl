using DataFrames
include("input/fleg-net.jl"); # loads packages & defines network "fleg" + "fleg_pruned"
using Format

# figure for v1 of manuscript, done with PhyloPlots v1
figure_filename = "figures/fig_flegsidebyside.pdf"
"""
(a) fleg (unpruned)
(b) fleg (pruned)
"""
# data frame to annotate edges in the full network
edgelabel = ["$(format(Int(e.length),commas=true)) gen.\nN=$(format(flegpopsizes[e.number], commas=true))\n" for e in fleg.edge];
for (i,e) in enumerate(fleg.edge)
    e.ismajor && continue
    edgelabel[i] *= "γ=$(round(e.gamma, digits=2))"
end
edgelabel_df = DataFrame(num=[e.number for e in fleg.edge], lab=edgelabel);

# fuse edges where hybridization events used to be
for i in [22,19,18,12,6]
  PhyloNetworks.fuseedgesat!(i, fleg_pruned)
end

# data frame to annotate edges in the pruned network
edgelabel2 = ["$(format(Int(e.length),commas=true)) gen.\nN=$(format(flegpopsizes[e.number], commas=true))\n" for e in fleg_pruned.edge];
for (i,e) in enumerate(fleg_pruned.edge)
  e.ismajor && continue
  edgelabel2[i] *= "γ=$(round(e.gamma, digits=2))"
end
edgelabel_df2 = DataFrame(num=[e.number for e in fleg_pruned.edge], lab=edgelabel2);

# capitalize Neanderthal, Africa, Denisovan; remove underscores
for net in (fleg_pruned, fleg), n in net.leaf
  n.name = replace(n.name, "non_" => "non-", "_" => " ","chimp" => "chimpanzee",
    "nean" => "Nean", "africa" => "Africa", "deni" => "Deni")
end

# create the figure
# R"capabilities"("cairo") # if TRUE, use cairo_pdf below to get the gamma's
# R"pdf"(figure_filename, height=8, width=10);
R"cairo_pdf"(figure_filename, height=6, width=8);
R"layout"([1 1; 2 2]);
R"par"(mar=[0,0,0,0]);
plot(fleg_pruned, xlim=[0.5,9.3], ylim=[0.8, 9.8], tipoffset=0.2,
     edgelabelcolor="red4", edgelabel = edgelabel_df2, edgecex=0.6);
R"text"(x=0.8, y=2.5, "Ne=$(format(flegpopsizes[28], commas=true))", adj=1, cex=0.6, col="red4");
R"mtext"("a)", line=-2, adj=0.01);
plot(fleg, xlim=[0.2,14], ylim = [0.8, 12.8], tipoffset=0.2,
     edgelabelcolor="red4", edgelabel = edgelabel_df, edgecex=0.6);
R"text"(x=0.8, y=4.5, "Ne=$(format(flegpopsizes[28], commas=true))", adj=1, cex=0.6, col="red4");
R"mtext"("b)", line=-2, adj=0.01);
R"dev.off"();

# post-processing: edge labels are too high. edited in illustrator to move them
# down, all by the same amount. (select 1, then select all with same fill color)

#=
figure for v2 of manuscript

done with PhyloPlots v2.0.1-dev, to use new feature for adjusting the placement
of edge labels, removing the need for post-processing with illustrator.
=#
figure_filename = "figures/fig_g4_annotated.pdf"
"""
fleg network g4, with bottleneck, annotated with Ne and # of generations,
and with edge colors to show which edges are absent in g1.
"""

edgecol = Dict(e.number => "black" for e in fleg.edge)
for i in [4,3] edgecol[i] = "deepskyblue"; end
for i in [7,11,12,26] edgecol[i] = "orangered"; end

R"cairo_pdf"(figure_filename, height=4, width=8);
R"par"(mar=[0,0,0,0]);
plot(fleg, xlim=[0.6,13.5], ylim = [0.8, 12.2], tipoffset=0.2, edgecolor=edgecol,
     edgelabelcolor="red4", edgelabel=edgelabel_df, edgecex=0.65, tipcex=0.8,
     edgelabeladj=[0.5,0.7]);
R"text"(x=0.95, y=4.5, "N=$(format(flegpopsizes[28], commas=true))", adj=1,
        cex=0.65, col="red4");
R"dev.off"();

# fixit: edge annotations are too high. edited in illustrator to move them down,
# all by the same amount. (select 1, then select all with same fill color)
