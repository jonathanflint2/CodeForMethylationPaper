 
 invnorm <- function (x) {
 y = (rank(x,na.last="keep")-0.5)/sum(!is.na(x))
 return (qnorm(y))
 }
 
 par (mfrow = c (1,2))
 
 PlotTissueWithAndWithoutMutations <- function(x, tissue = "", logp = 1.3, threshold = 40){
 	x = dg[dg$rna.logp > logp & abs(dg$rna.log2fold) <=2 & !is.na(dg$rna.logp) ,]
  x$mut.b6.strict = 0
  x$mut.b6.strict[x$b6.var == "CG"] = 0
  x$mut.b6.strict[x$b6.var == "T"] = 1
  x$mut.b6.strict[x$b6.var == "A"] = 1
   # set density threshold 
  x$high = 0
  x$high[x$density  > threshold] = 1
  x$high[x$density <= threshold] = 0
  #make direction of effect increasing for B6
  x$rna.log2fold  = -1*x$rna.log2fold
  #quantile normalize the RNA change just to be sure...
  x$normalized.log2fold = invnorm(x$rna.log2fold)
  y = subset (x, d2.me > 0 & !is.na(d2.me))
  high  = y[y$high == 1,]
  low  =  y[y$high == 0,]
  
   # plot the RNA abundance for data WITH mutations 
  high.mean.mut <- mean(na.omit(high$normalized.log2fold[high$mut.b6.strict==1]))
  high.mut.lower_ci <- high.mean.mut - 1.96 * sd(na.omit(high$normalized.log2fold[high$mut.b6.strict==1]))/sqrt(length(na.omit(high$normalized.log2fold[high$mut.b6.strict==1])))
  high.mut.upper_ci <- high.mean.mut + 1.96 * sd(na.omit(high$normalized.log2fold[high$mut.b6.strict==1])) / sqrt(length(na.omit(high$normalized.log2fold[high$mut.b6.strict==1])))
  low.mean.mut <- mean(na.omit(low$normalized.log2fold[low$mut.b6.strict==1]))
  low.mut.lower_ci <- low.mean.mut - 1.96 * sd(na.omit(low$normalized.log2fold[low$mut.b6.strict==1])) / sqrt(length(na.omit(low$normalized.log2fold[low$mut.b6.strict==1])))
  low.mut.upper_ci <- low.mean.mut + 1.96 * sd(na.omit(low$normalized.log2fold[low$mut.b6.strict==1])) / sqrt(length(na.omit(low$normalized.log2fold[low$mut.b6.strict==1])))
  greycol = rgb(0.34,0.34,0.34)
  stripchart(list (high$normalized.log2fold[high$mut.b6.strict==1],low$normalized.log2fold[low$mut.b6.strict==1]), at = c(1.2,1.6),  method = "jitter", vertical = TRUE, col = "darkblue", pch = 16, cex = 0.4, group.names = c ("CpG high density", "CpG low density"), cex.axis = 0.7, ylab = "", add = F, xlim = c (1,1.8),las = 1, ylim = c(-2,2), main = tissue)
  segments (1.05, high.mean.mut, 1.35, high.mean.mut,col = "black", lwd = 4, lty = 1) 
  segments (1.05, high.mut.upper_ci, 1.35, high.mut.upper_ci,col = greycol, lwd = 3, lty = 3 )
  segments (1.05, high.mut.lower_ci, 1.35, high.mut.lower_ci,col = greycol, lwd = 3, lty = 3) 
   segments (1.45, low.mean.mut, 1.75, low.mean.mut,col = "black", lwd = 4, lty = 1)
   segments (1.45, low.mut.upper_ci, 1.75, low.mut.upper_ci,col = greycol, lwd = 3, lty = 3)
  segments (1.45, low.mut.lower_ci, 1.75, low.mut.lower_ci,col = greycol, lwd = 3, lty = 3)
  
 # plot the RNA abundance for data WITHOUT mutations 
  high.mean.mut <- mean(na.omit(high$normalized.log2fold[high$mut.b6.strict==0]))
  high.mut.lower_ci <- high.mean.mut - 1.96 * sd(na.omit(high$normalized.log2fold[high$mut.b6.strict==0]))/sqrt(length(na.omit(high$normalized.log2fold[high$mut.b6.strict==0])))
  high.mut.upper_ci <- high.mean.mut + 1.96 * sd(na.omit(high$normalized.log2fold[high$mut.b6.strict==0])) / sqrt(length(na.omit(high$normalized.log2fold[high$mut.b6.strict==0])))
  low.mean.mut <- mean(na.omit(low$normalized.log2fold[low$mut.b6.strict==0]))
  low.mut.lower_ci <- low.mean.mut - 1.96 * sd(na.omit(low$normalized.log2fold[low$mut.b6.strict==0])) / sqrt(length(na.omit(low$normalized.log2fold[low$mut.b6.strict==0])))
  low.mut.upper_ci <- low.mean.mut + 1.96 * sd(na.omit(low$normalized.log2fold[low$mut.b6.strict==0])) / sqrt(length(na.omit(low$normalized.log2fold[low$mut.b6.strict==0])))
  greycol = rgb(0.34,0.34,0.34)
  stripchart(list (high$normalized.log2fold[high$mut.b6.strict==0],low$normalized.log2fold[low$mut.b6.strict==0]), at = c(1.2,1.6),  method = "jitter", vertical = TRUE, col = "darkblue", pch = 16, cex = 0.4, group.names = c ("CpG high density", "CpG low density"), cex.axis = 0.7, ylab = "", add = F, xlim = c (1,1.8),las = 1, ylim = c(-2,2), main = tissue)
  segments (1.05, high.mean.mut, 1.35, high.mean.mut,col = "black", lwd = 4, lty = 1) 
  segments (1.05, high.mut.upper_ci, 1.35, high.mut.upper_ci,col = greycol, lwd = 3, lty = 3 )
  segments (1.05, high.mut.lower_ci, 1.35, high.mut.lower_ci,col = greycol, lwd = 3, lty = 3) 
   segments (1.45, low.mean.mut, 1.75, low.mean.mut,col = "black", lwd = 4, lty = 1)
   segments (1.45, low.mut.upper_ci, 1.75, low.mut.upper_ci,col = greycol, lwd = 3, lty = 3)
  segments (1.45, low.mut.lower_ci, 1.75, low.mut.lower_ci,col = greycol, lwd = 3, lty = 3)
  
  }