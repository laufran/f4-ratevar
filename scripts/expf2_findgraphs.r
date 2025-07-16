library(admixtools)
library(ggplot2)
library(dplyr)

#' 3d-array of f2 values that can be used by qpgraph and find_graphs
#'
#' The same f2 matrix is repeated twice (3rd dimension) to pretend that there
#' were 2 genomic blocks, and both have the same f2 values.
#' @param csvfile Name of csv file containing the f2 values. There should be one
#' column per population, with population names given in the header. Rows should
#' consider populations in the same order.
make_f2input <- function(csvfile){
  f2_mat <- as.matrix(read.csv(csvfile))
  popnames <- colnames(f2_mat)
  rownames(f2_mat) <- popnames
  array(f2_mat, c(8,8,2), list(popnames, popnames, c("l1","l2")))
}

#' fit f2 data to a fixed graph topology
#'
#' The edge weights are save to a csv file, and a figure of the fitted graph is
#' saved as a pdf.
#' @param f2array f2 data in the form a 3d-array
#' @param graph fixed graph in edgelist format. only the topology is used
#' @param rootname the out files will be named rootname_seed.{csv,pdf}
#' @param seed passed on to qpgraph
fitf2_qpgraph <- function(f2array, graph, rootname, seed=NULL){
  qpres <- qpgraph(f2array, graph, lsqmode=2, seed=seed, return_fstats=TRUE)
  rootname <-  paste0(rootname, "_", seed)
  write.csv(qpres$edges, paste0(rootname, ".csv"), quote=FALSE, row.names=FALSE)
  ggsave(paste0(rootname, ".pdf"), plot_graph(qpres$edges), height=6, width=6)
  return(qpres)
}

#' fit f2 data to a list of candidate network topologies
#'
#' output: data frame with 2 columns: network names & scores.
#' @param f2array f2 data in the form a 3d-array
#' @param graphs list of candidate graphs, each one as an edge list data frame.
#' This list should have names, and these names form a column in the output.
#' @param rootdir the output data frame is saved to: rootdir/netscores_seed.csv
#' @param seed passed on to qpgraph, the same for fitting all networks.
fitf2_candidatenets <- function(f2array, graphs, rootdir, seed=NULL){
  gname <- names(graphs)
  ng <- length(gname)
  gscore <- vector("numeric", ng)
  gdiff <- vector("numeric", ng)
  for (i in 1:ng) {
    g <- graphs[[i]]
    qpres <- qpgraph(f2array, g, lsqmode=2, seed=seed)
    gscore[i] <- qpres$score
    gdiff[i] <- max(abs(qpres$f3$diff))
  }
  dat <- data.frame(net=gname, score=gscore, worstresid=gdiff)
  outcsv <- paste0(rootdir, "netscores_", seed, ".csv")
  write.csv(dat, outcsv, quote=FALSE, row.names=FALSE)
  dat
}

#----------------------------------------------------------#
#     read 2 lists of all candidate networks & graphs      #
# population networks g1-*, g2-*, g3-*, g4-*
graphfile <- grep("^g.*_cu.csv", list.files("input/pop-net_edgelist/"), value=TRUE)
ng <- length(graphfile)
pop_edgelist <- vector("list", ng)
for (i in 1:ng){
  gf <- graphfile[i]
  pop_edgelist[[i]] <- read.csv(paste0("input/pop-net_edgelist/", gf))
  names(pop_edgelist)[i] <- sub("_cu.csv", "", gf)
}
# graphs found as best-scoring on individual simulated reps
d <- "output/g14_topgraphs_byrep_scoretol0.05_correctToB_bestscenario_edgelist/"
graphfile <- grep(".csv", list.files(d), value=TRUE)
ng <- length(graphfile)
g_edgelist <- vector("list", ng+1) # +1 to add g4 at the end
for (i in 1:ng){
  gf <- graphfile[i]
  g_edgelist[[i]] <- read.csv(paste0(d, gf))
  names(g_edgelist)[i] <- sub(".csv", "", gf) # name = hash
}
# add g4 at the end
flegnobott_hash = "b66401d14d13e41dfd92363767147547"
i = ng+1
g_edgelist[[i]] <- read.csv("input/pop-net_edgelist/g4-sg_cu.csv")
names(g_edgelist)[i] <- flegnobott_hash

stopifnot(length(g_edgelist) == 31+1)
rm(i, gf, graphfile, ng, d)

#----------------------------------------------------------#
#          average f2 data from best-case scenario         #
nind <- 10
outdir <- paste0("output/fleg_", nind,
    "-indls_118000-genes_100000-biall_0.0-lindist_1.25e-8-subrate/")

f2csvfile <- paste0(outdir,"avgf2.csv")
qpdir_obsf2 <- paste0(outdir, "avgf2_fit/")
f2_3d <- make_f2input(f2csvfile)

scorelims =  c(1/1.1, 1.1) * c(1.957247e-09, 6.281808e-05) # from range(res_exp$score) later
wrlims = c(1/1.1, 1.1) * c(0.0000185872, 0.004956407)

#----------------------------------------------------------#
# fit average observed f2 to candidate networks & graphs   #

res_n <- fitf2_candidatenets(f2_3d, pop_edgelist, qpdir_obsf2, seed=7630)
# res_n <- read.csv(paste0(qpdir_obsf2, "netscores_7630.csv"))
res_n$isg4 <- res_n$net %in% c("g4-lg", "g4-sg")
res_n$h <- as.numeric(sub("g(\\d)-.*", "\\1", res_n$net))
res_n$row <- NA

#---- fit on graphs found as best-scoring on individual simulated reps
res_g <- fitf2_candidatenets(f2_3d, g_edgelist, qpdir_obsf2, seed=3722)
# res_g <- read.csv(paste0(qpdir_obsf2, "netscores_3722.csv"))
res_g$isg4 <- res_g$net == flegnobott_hash
res_g$h <- 4
res_g$row <- seq_len(nrow(res_g))

#---- merge the two lists and plot
res_avg <- bind_rows(res_n, res_g) %>% arrange(score)
write.csv(res_avg, paste0(qpdir_obsf2, "netscores.csv"), quote=FALSE, row.names=FALSE)
# graphs that fit better than g4 are on rows: 20,28,19,g4-l4-c34,2,30,14
# network with h=3 with worse but similar fit on rows: 8,1,g3-tc
p_avg <- ggplot(res_avg,
    aes(x=score, y=worstresid, color=as.factor(h), shape=isg4)) +
  geom_point() +
  scale_x_log10(minor_breaks=NULL, name="network score", limits=scorelims) +
  scale_y_log10(minor_breaks=NULL, name="worst residual", limits=wrlims) +
  labs(caption = "f2: averaged over simulated replicates") +
  scale_color_manual(
    values = c("sienna1","tomato3","tomato4","black"),
    labels = paste("h =",1:4)) +
  scale_shape_manual(values=c(21,16)) +
  guides(
    color=guide_legend(title="", override.aes=list(shape=21)),
    shape="none") +
  theme_minimal() + theme(
    legend.position = c(1,0), legend.justification=c(1,0),
    plot.margin = margin(.1,.2,.1,.1,"in"))
ggsave(paste0(qpdir_obsf2,"netscores.pdf"), p_avg, width=4, height=4)

#----------------------------------------------------------#
#   read g4, aka g4-sg or fleg: are parameters identifiable?
g4el = read.csv("input/pop-net_edgelist/g4-sg_cu.csv")

plot_graph(g4el, highlight_unidentifiable=TRUE)
ggsave("figures/fleg-unidentifiable.pdf", height=6, width=6)
# all but 1 edge weights are unidentifiable!!
print(unidentifiable_edges(edges_to_igraph(g4el)), n=32)

#   fit on g4 then extract f2's expected from this fitted graph
# unstable estimates of branch lengths and gammas, but excellent fit always:
nameroot <- paste0(qpdir_obsf2, "g4")

res <- fitf2_qpgraph(f2_3d, g4el, nameroot, seed=71)
res$score # 5.631094e-09, but similar edge lengths as with seed 651 below

res <- fitf2_qpgraph(f2_3d, g4el, nameroot, seed=9148)
res$score # 7.707295e-09, very different edge lengths near the root
min(abs(res$f3$diff)) # 8.48807e-07
max(abs(res$f3$diff/res$f3$est)) # 0.0001561513
min(abs(res$f2$diff)) # 8.48807e-07
max(abs(res$f2$diff/res$f2$est)) # 0.001671811
f2expg4_s9148 <- res$f2[,c("pop1","pop2","fit")]
summary(f2expg4_s9148$fit/f2expg4$fit)
# min=0.9985, Q1=median=1.0000, Q3=1.0002, max=1.0013

res <- fitf2_qpgraph(f2_3d, g4el, nameroot, seed=651)
res$score # 4.83908e-09: best (lowest) of all 3 seeds that were tried
print(res$edges, n=33)
# to look at identifiability & instability:
res$opt # i = initial weights, e = estimated weights
min(abs(res$f3$est))  # 0.279371
max(abs(res$f3$diff)) # 3.386108e-05: that's a really good fit
max(abs(res$f3$diff/res$f3$est)) # 0.0001203557
min(abs(res$f2$est))  # 0.00677829
max(abs(res$f2$diff)) # 7.358845e-05
max(abs(res$f2$diff/res$f2$est)) # 0.001801735
f2expg4 <- res$f2[,c("pop1","pop2","fit")]

p_f2s <- ggplot(res$f2, aes(x=est, y=fit)) +
  geom_abline(intercept=0, slope=1, color=grey(0.7)) +
  geom_point(shape=21) +
  scale_x_log10(name="f2: averaged over simulated replicates") +
  scale_y_log10(name="f2: expected from g4") +
  theme_minimal()
ggsave(paste0(outdir,"avgf2_vs_fitted.pdf"), p_f2s, width=4, height=4)

popnames <- rownames(f2_3d)
f2expg4_3d <- array(0, c(8,8,2), list(popnames, popnames, c("l1","l2")))
for (i in 2:8){
  for (j in 1:(i-1)){
    r = which((f2expg4$pop1 == popnames[i] & f2expg4$pop2 == popnames[j]) |
              (f2expg4$pop1 == popnames[j] & f2expg4$pop2 == popnames[i]) )
    v = f2expg4$fit[r]
    f2expg4_3d[i,j,1] <- v; f2expg4_3d[i,j,2] <- v
    f2expg4_3d[j,i,1] <- v; f2expg4_3d[j,i,2] <- v
  }
}
range(f2_3d-f2expg4_3d) # -4.887535e-05  7.358845e-05

write.csv(data.frame(f2expg4_3d[,,1]), paste0(outdir, "fittedf2.csv"),
  quote=FALSE, row.names=FALSE)
qpdir_expf2 <- paste0(outdir, "fittedf2_fit/")
if (!dir.exists(qpdir_expf2)){
  dir.create(qpdir_expf2)
}

#----------------------------------------------------------#
#   fit g4-*expected* f2 to candidate networks & graphs    #

res_n <- fitf2_candidatenets(f2expg4_3d, pop_edgelist, qpdir_expf2, seed=299)
# res_n <- read.csv(paste0(qpdir_expf2, "netscores_299.csv"))
res_n$isg4 <- res_n$net %in% c("g4-lg", "g4-sg")
res_n$h <- as.numeric(sub("g(\\d)-.*", "\\1", res_n$net))
res_n$row <- NA

#---- fit on graphs found as best-scoring on individual simulated reps
res_g <- fitf2_candidatenets(f2expg4_3d, g_edgelist, qpdir_expf2, seed=129)
# res_g <- read.csv(paste0(qpdir_expf2, "netscores_129.csv"))
res_g$isg4 <- res_g$net == flegnobott_hash
res_g$h <- 4
res_g$row <- seq_len(nrow(res_g))

#---- merge the two lists and plot
res_exp <- bind_rows(res_n, res_g) %>% arrange(score)
write.csv(res_exp, paste0(qpdir_expf2, "netscores.csv"), quote=FALSE, row.names=FALSE)
# graphs that fit better than g4 are on rows: g4-l4-c34,20,28,19,30,2,8
# network with h=3 with worse but similar fit on rows: 14,g3-tc
p_exp <- ggplot(res_exp,
    aes(x=score, y=worstresid, color=as.factor(h), shape=isg4)) +
  geom_point() +
  scale_x_log10(minor_breaks=NULL, name="network score", limits=scorelims) +
  scale_y_log10(minor_breaks=NULL, name="worst residual", limits=wrlims) +
  labs(caption = "f2: expected from g4") +
  scale_color_manual(
    values = c("sienna1","tomato3","tomato4","black"),
    labels = paste("h =",1:4)) +
  scale_shape_manual(values=c(21,16)) +
  guides(
    color=guide_legend(title="", override.aes=list(shape=21)),
    shape="none") +
  theme_minimal() + theme(
    legend.position = c(1,0), legend.justification=c(1,0),
    plot.margin = margin(.1,.2,.1,.1,"in"),
    # plot.subtitle = element_text(size = 15)
  )
ggsave(paste0(qpdir_expf2,"netscores.pdf"), p_exp, width=4, height=4)

#----------------------------------------------------------#
# figure: combine the main plots

library(ggpubr)
p <- ggarrange(p_f2s, p_avg, p_exp, ncol=3, nrow=1)
ggsave(paste0(outdir,"SMfig_netscores.pdf"), p, width=10, height=4)
