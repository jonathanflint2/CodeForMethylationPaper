#############################################
# First create summed density plots using the script SumByDensity.pl on each file that contains 
# pairwise data from each strain comparison
#
#SumByDensity.pl -f $file > $summed
#############################################
# Read in the file (here using Exc-DG cell type)
 x = read.delim("Exc-DG.summed.txt")
 # for other cell types....  x = read.delim ("Exc-CA3.summed.txt") etc
 ###############
# Column headings in the summed file mean the following:
###############

  # how many sites are CpG?     CG.all  
  # how many of these sites are methylated? CG.me.1 blue
  # how many of these sites are unmethylated? CG.me.0  grey 
  # how many of the methylated CpG sites are mutated? mutant.me.1           
  # how many of the unmethylated CpG sites are mutated? mutant.me.0   
  
  ###############  
  # set the values to be plotted as percentages
  ###############
 x$pct.mutant.me.0  = 100 * x$mutant.me.0 / x$CG.all
 x$pct.mutant.me.1  = 100 * x$mutant.me.1 / x$CG.all
 x$pct.CG.me.1  = 100 * x$CG.me.1 / x$CG.all
 x$pct.CG.me.0  = 100 * x$CG.me.0 / x$CG.all
 
 ###############
  # Each plot total uses CAST as an outgroup so first extract those where the first strain (s1) is cast
  ###############
a = subset (x, s1 == "cast")
strain2 = unique(a$s2)

# color palette
me0colors = c ("gray1", "gray10", "gray20", "gray30", "gray40", "gray50", "gray60", "gray70")
me1colors = c ("dodgerblue1", "deepskyblue1", "dodgerblue2", "deepskyblue2", "dodgerblue3", "deepskyblue3" ,"dodgerblue4", "deepskyblue4")

#dev.new()
#par (mfrow = c (3,1))
############################################################
# First plot: decay of methylation with sequence density
############################################################
###############
# set the plots to show a density of CpGs per 1Kb - 
# most files show density per 2kb so divide by 2 for the plotting axis
#
# For other segment densities correct by the relevant factor
#factor = 2
###############
plot (x$seg/factor, x$pct.mutant.me.0, pch = 19, xlim = c (0,100) , ylim = c (0, 100), xlab = "CpG density (CpGs/Kb)", ylab = "Percent methylated" , type =  "n")

# plot for each strain

for (i in 1:length(strain2)) {
	y = subset (a , s2 == strain2[i])
	lines (y$seg/factor, y$pct.CG.me.1, pch = 19, col = me1colors[i])	
	}
	legend ("topright", bty = "n", pch = c(16,16,16,16,16,16,16), col = me1colors[1:7], legend = strain2)
		
 ############################################################
 # Second plot: relationship between mutation and methylation, with respect to CpG density
 ############################################################

	  
  plot (x$seg/factor, x$pct.mutant.me.0, pch = 19, xlim = c (0,100) , ylim = c (0, 0.75), xlab = "CpG density (CpGs/Kb)", ylab = "Percent mutant" , type =  "n")
 for (i in 1:length(strain2)) {
 	y = subset (a, s2 == strain2[i])
 	lines (y$seg/factor, y$pct.mutant.me.0, pch = 19, col = me0colors[i])
	lines (y$seg/factor, y$pct.mutant.me.1, pch = 19, col = me1colors[i])		

 }
 legend ("topright", bty = "n", pch = c(16,16), col = c(me1colors[1], me0colors[1]), legend = c("Methylated",  "Unmethylated"))
  
  
 ############################################################
 # Third plot: relationship between probability of methylation state and ratio of the number of methylated sites
 # in a pairwise comparison
 ############################################################	
 
 plot (x$seg/factor, x$ratio, pch = 19, xlim = c (0,100) , ylim = c (0.6, 1), xlab = "CpG density (CpGs/Kb)", ylab = "Ratio or probability" , type =  "n")

  for (i in 1:length(strain2)) {
 	y = subset (a, s2 == strain2[i])
 lines (y$seg/factor, y$prob,  col = me0colors[i])
 lines (y$seg/factor, y$ratio, col = me1colors[i])	
 }
 legend ("bottom", bty = "n", pch = c(16,16), col = c(me1colors[1], me0colors[1]), legend = c("ratio of number of CpGs",  "probability of being the same state"))


dev.new()
par (mfrow = c (3,1))

