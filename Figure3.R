# Note: Run sim-study-sec4-sec5.R, and 
# sim-study-sec4-sec5-collate-results.R to 
# do the simulation study and format results first

#### Load libraries & functions ####
library(ggplot2)
library(gridExtra)
library(grid)
library(scales)

# Function to arrange plots in a row
grid_arrange_shared_legend <- function(..., nrow = 1, ncol = length(list(...)), position = c("bottom", "right")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position = "none"))
  gl <- c(gl, nrow = nrow, ncol = ncol)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  grid.draw(combined)
  
}


##### Generate plots comparing performance of Meila to PLRT, p = 10, K = 6####
K=6
p=10

# Figure S5a
sig=2.4
load(file=paste("df-K", K, "-p", p, "-sig", 
                sig, ".Rdata", sep=""))

df_comp_meila <- df[df$Method %in% c("test_indep_per", "meila_mclust", 
                                     "meila_mclust_per"), ]
df_comp_meila$Method <- factor(df_comp_meila$Method, 
                               levels=c("test_indep_per", "meila_mclust", 
                                        "meila_mclust_per"))
lty <- setNames(c(5, 4, 3, 2, 1), levels(factor(df_comp_meila$n)))

p1_comp_meila <- ggplot(aes(x = delta, y = power, 
                            ymin = power - se, ymax = power + se, 
                            color = df_comp_meila$Method, 
                            group = interaction(df_comp_meila$n, df_comp_meila$Method)),
                        linetype=lty,
                        data=df_comp_meila) + 
  geom_hline(yintercept=0.05, linetype="dotted") + 
  annotate(geom="text", label=as.character(expression(alpha==~"0.05")),parse=T,  
           x=0.08, y=0, vjust=0.3) + 
  geom_line(aes(linetype=as.factor(df_comp_meila$n))) +  geom_linerange() + 
  xlab("Dependence between views") +
  ylab("Power at nominal 5% significance level") +
  scale_colour_manual(name = "Methods", 
                        labels=c( "Pseudo LRT", 
                                  "G-test (Chi squared)", 
                                  "G-test (Permutation)"), 
                        values = c("#1f78b4", "#b2df8a", "#33a02c"))   + 
  scale_linetype_manual(name = "  n ", values=lty) +
  ggtitle(expression(sigma==~2.4)) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_gray(base_size = 15.5) + coord_cartesian(ylim = c(0, 1))
p1_comp_meila

# Figure S5b
sig=4.8
load(file=paste("df-K", K, "-p", p, "-sig", 
                sig, ".Rdata", sep=""))

df_comp_meila2 <- df[df$Method %in% c("test_indep_per", "meila_mclust", 
                                      "meila_mclust_per"), ]
df_comp_meila2$Method <- factor(df_comp_meila2$Method, 
                                levels=c("test_indep_per", "meila_mclust", 
                                         "meila_mclust_per"))
lty2 <- setNames(c(5, 4, 3, 2, 1), levels(factor(df_comp_meila2$n)))
p2_comp_meila <- ggplot(aes(x = delta, y = power, 
                            ymin = power - se, ymax = power + se, 
                            color = df_comp_meila2$Method, 
                            group = interaction(df_comp_meila2$n, 
                                                df_comp_meila2$Method)),
                        linetype=lty2, 
                        data=df_comp_meila2) + 
  geom_hline(yintercept=0.05, linetype="dotted") + 
  geom_line(aes(linetype=as.factor(df_comp_meila2$n))) +  geom_linerange() + 
  xlab("Dependence between views") + ylab("") + 
  scale_colour_manual(name = "Methods", 
                      labels=c( "Pseudo LRT", 
                                "G-test (Chi squared)", 
                                "G-test (Permutation)"), 
                      values = c("#1f78b4", "#b2df8a", "#33a02c"))   + 
  scale_linetype_manual(name = "  n ", values=lty2) +
  ggtitle(expression(sigma==~4.8)) +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(), 
        axis.text.y=element_blank()) + 
  theme_gray(base_size = 15.5) + coord_cartesian(ylim = c(0, 1))
p2_comp_meila

# Figure S5c
sig=9.6
load(file=paste("df-K", K, "-p", p, "-sig", 
                sig, ".Rdata", sep=""))
df_comp_meila3 <- df[df$Method %in% c("test_indep_per", "meila_mclust", 
                                      "meila_mclust_per"), ]
df_comp_meila3$Method <- factor(df_comp_meila3$Method, 
                                levels=c("test_indep_per", "meila_mclust", 
                                         "meila_mclust_per"))
lty3 <- setNames(c(5, 4, 3, 2, 1), levels(factor(df_comp_meila3$n)))
p3_comp_meila <- ggplot(aes(x = delta, y = power, 
                            ymin = power - se, ymax = power + se, 
                            color = df_comp_meila3$Method, 
                            group = interaction(df_comp_meila3$n,
                                                df_comp_meila3$Method)), 
                        linetype=lty3,
                        data=df_comp_meila3) + 
  geom_hline(yintercept=0.05, linetype="dotted") + 
  geom_line(aes(linetype=as.factor(df_comp_meila3$n))) +  geom_linerange() + 
  xlab("Dependence between views") +ylab("") + 
  scale_colour_manual(name = "Methods", 
                      labels=c( "Pseudo LRT", 
                                "G-test (Chi squared)", 
                                "G-test (Permutation)"), 
                      values = c("#1f78b4", "#b2df8a", "#33a02c"))   + 
  scale_linetype_manual(name = "  n ", values=lty3) +
  ggtitle(expression(sigma==~9.6)) + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(), 
        axis.text.y=element_blank()) + 
  theme_gray(base_size = 15.5) + coord_cartesian(ylim = c(0, 1))
p3_comp_meila

# Arrange in a row
pdf("Figure3.pdf", width=12, height=6)
grid_arrange_shared_legend(p1_comp_meila, p2_comp_meila, 
                           p3_comp_meila, nrow=1, ncol=3)
dev.off()


