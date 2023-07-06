
invnorm <- function (x) {
 y = (rank(x,na.last="keep")-0.5)/sum(!is.na(x))
 return (qnorm(y))
 }
dg = read.delim ("Exc-DG.combined.modalities.230623.txt")
threshold = 40
# select out the transcripts we are interested in
x = dg[dg$rna.logp > 1.3 & abs(dg$rna.log2fold) <= 2 & !is.na(dg$rna.logp) ,]
x$mut.b6.strict = 0
x$mut.b6.strict[x$b6.var == "CG"] = 0
x$mut.b6.strict[x$b6.var == "T"] = 1
x$mut.b6.strict[x$b6.var == "A"] = 1
x$high[x$density  > threshold] = 1
x$high[x$density <= threshold] = 0
#make direction of effect increasing for B6
x$rna.log2fold  = -1*x$rna.log2fold
x$normalized.log2fold = invnorm(x$rna.log2fold)
#make states specific to the ATAC data - hard coded for count number threshold
x$specific.states  = 0
x$specific.states[!is.na(x$atac.counts) & x$atac.counts > 10]  = x$states[!is.na(x$atac.counts) & x$atac.counts > 10]
#we are only interested in sites that are methylated in D2
y = subset (x, d2.me > 0 & !is.na(d2.me))
high  = y[y$high == 1,]
low  =  y[y$high == 0,]

# compare effect of mutations in high and low density 

t.test (high$normalized.log2fold[high$mut.b6.strict==1], low$normalized.log2fold[low$mut.b6.strict==1])

# compare effect of mutations in high with effect of high density alone (no mutations)

t.test (high$normalized.log2fold[high$mut.b6.strict==1], high$normalized.log2fold[high$mut.b6.strict==0])

# comparing the effect attributable to other confounds

fit0 = lm (normalized.log2fold ~ specific.states +  high + mut.b6.strict, data = y)
fit1 = lm (normalized.log2fold ~ specific.states +  high + mut.b6.strict + high:mut.b6.strict, data = y)

anova (fit0, fit1)

# are the p-values well calibrated? Test the above by sampling:

permute_and_fit <- function(df) {	
# Extract unique gene names from the data frame
gene_names <- unique(df$gene)
# Extract corresponding values from the 'normalized.log2fold' column
values <- df$normalized.log2fold
# Create the lookup table
lookup_table <- setNames(values, gene_names)
gene_names <- names(lookup_table)
# Randomly shuffle the values in the lookup table
shuffled_values <- sample(lookup_table)
# Create a new lookup table with shuffled values
shuffled_lookup_table <- setNames(shuffled_values, gene_names)
df$values <- shuffled_lookup_table[df$gene]
# run the two linear models and report back the P-value
 fit0 <- lm(values ~ specific.states +  high + mut.b6.strict, data = df)
  fit1 <- lm(values~ specific.states +  high + mut.b6.strict + high:mut.b6.strict, data = df)
  a = anova(fit0,fit1)
  p_value <- a$"Pr(>F)"[2]
  return(p_value)
}



