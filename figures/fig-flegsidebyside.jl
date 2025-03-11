using DataFrames
include("../../input/fleg-net.jl"); # loads packages & defines network "fleg" + "fleg_pruned"
using Format

figure_filename = "figures/fig_flegsidebyside.pdf"
"""
(a) fleg (unpruned)
(b) fleg (pruned)
"""
# data frame to annotate edges in the full network
edgelabel = ["$(format(Int(e.length),commas=true)) gen.\nN=$(format(flegpopsizes[e.number], commas=true))\n" for e in fleg.edge];
for (i,e) in enumerate(fleg.edge)
    e.isMajor && continue
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
  e.isMajor && continue
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

# fixit: edge annotations are too high. edited in illustrator to move them down,
# all by the same amount. (select 1, then select all with same fill color)
