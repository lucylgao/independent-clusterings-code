# run dense-diagonal-gaussian-sim-code.R, dense-gaussian-sim-code.R, 
# dense-t-sim-code.R first  
library(gridExtra)
library(grid)
library(scales)
library(ggplot2)
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

make_power_plot <- function(simId, caseId) { 
  delta_seq <- seq(0, 1, length.out=9)
  n_seq <- c(25, 50, 100, 250, 500)
  
  summary_table <- data.frame(n = numeric(), 
                              delta = numeric(), 
                              power = numeric(), 
                              nSim = numeric(), 
                              methods=character())
  
  for(i in 1:length(n_seq)) { 
    for(j in 1:length(delta_seq)) { 
      results_cat <- c()
      for(seed in 1:12) { 
        load(file=paste(simId, "-results-n", n_seq[i], 
                        "-delta", delta_seq[j], 
                        "-seed", seed, ".Rdata", sep=""))
        results_cat <- cbind(results_cat, results)
      }  
      
      nSim <- (12*167 - apply(results_cat[1:4, ], 1, 
                              function(x) sum(is.na(x))))
      power <- rowSums(results_cat[1:4, ] <= 0.05, 
                       na.rm=T)/nSim
      
      summary_table <- rbind(summary_table, 
                             data.frame(n = n_seq[i], 
                                        delta = delta_seq[j], 
                                        power = power, nSim = nSim, 
                                        methods=c("PLRT-K3", "PLRT-KBIC", 
                                                  "meila-K3", "meila-KBIC")))
    }
  }
  rownames(summary_table) <- NULL 
  
  summary_table$power.se <- 
    sqrt(summary_table$power*(1 - summary_table$power)/
           summary_table$nSim)
  
  summary_table <- summary_table[summary_table$methods %in% c("PLRT-K3", "meila-K3"), ]
  
  
  resultPlot <- ggplot(aes(x = delta, y = power, 
                           ymin = power - power.se, ymax = power + power.se, 
                           color = methods, 
                           group = interaction(as.factor(n), as.factor(methods))),
                       data=summary_table) + 
    geom_hline(yintercept=0.05, linetype="dotted") + 
    annotate(geom="text", label=as.character(expression(alpha==~"0.05")),parse=T,  
             x=0.08, y=0.03) + 
    geom_line(aes(linetype=as.factor(-n)), 
              data=summary_table) +  
    geom_linerange() + 
    xlab("Dependence between views") +
    ylab("Power at nominal 5% significance level") +
    scale_colour_manual(name = "Methods", 
                        labels=c( "G-test",
                                  "Pseudo LRT"), 
                        values = c("#a6611a", "#80cdc1")) + 
    scale_linetype_discrete(name = "n", 
                            labels=c("500", "250", "100", "50", "25"))+ 
    theme_gray(base_size = 15.5)
  
  return(resultPlot)
}

p1 <- make_power_plot("dense-gaussian", 1)
p2 <- make_power_plot("dense-diagonal-gaussian", 2)
p3 <- make_power_plot("dense-t", 3)

pdf("FigureS9.pdf", width=12, height=6)
grid_arrange_shared_legend(p1, p2, p3, nrow=1, ncol=3)
dev.off()


  