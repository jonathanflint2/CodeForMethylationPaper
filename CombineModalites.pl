#!/usr/bin/perl
use strict;
use Getopt::Std;

our($opt_h, $opt_p,$opt_f, $opt_z,$opt_a,$opt_r, $opt_c, $opt_d, $opt_F, $opt_S, $opt_i, $opt_A, $opt_R);

getopts("hf:a:p:r:c:vd:F:S:i:A:R:");
&usage if $opt_h;
sub usage {
  die "CombineModalities.pl -d density file (e.g. pairwise) -r rna DESEQ file -R RNA count file -a atac DESEQ file -A atac count file -f counted file -c chromosome
counted file should be for a single chromosome, same as give by -c option
counted file must have a column headed \"gene\"
strains are assumed to be b6 and d2
change this with -F and -S
-i interval size of the pairwise density file (default is 2000 - assumes density is reported in 2kb chuncks )
\n";
}
sub mean {
    my ($array_ref) = @_;
    my $sum = 0;
    foreach my $value (@$array_ref) {
        $sum += $value;
    }
    my $mean = $sum / scalar(@$array_ref);
    return $mean;
  }
sub round {
    my ($number) = @_;

    my $rounded = int($number);
    $rounded++ if ($number - int($number)) >= 0.5;

    return $rounded;
  }




my $AllowedStrains = {};
$AllowedStrains->{'aj'} = 0;
$AllowedStrains->{'balb'} = 1;
$AllowedStrains->{'b6'} = 2;
$AllowedStrains->{'cast'} = 3;
$AllowedStrains->{'d2'} = 4;
$AllowedStrains->{'fvb'} = 5;
$AllowedStrains->{'pwk'} = 6;
$AllowedStrains->{'wsb'} = 7 ;

# Define interval size

my $interval = defined ($opt_i)? $opt_i : 2000;

# Define strains to get
my $firststrain = defined ($opt_F) ? $opt_F : "b6";
my $secondstrain = defined ($opt_S) ? $opt_S : "d2";

unless (exists ($AllowedStrains->{$firststrain}) && exists ($AllowedStrains->{$secondstrain}) ) {
  die "can't find $firststrain or $secondstrain in the list of allowed strains\n";
}
my $chr_to_use;
if (defined($opt_c)) {
  $chr_to_use = $opt_c;
}
else {
  &usage;
}



my $densityfile = $opt_d;
my $rnafile = $opt_r;
my $rnacountfile = $opt_R;
my $atacfile = $opt_a;
my $ataccountfile = $opt_A;
my $countedfile = $opt_f;

my $Density = {};
my $RNA = {};
my $ATAC = {};
my $RNACount = {};
my $ATACCount = {};

open (DENSITYFILE, $densityfile) || die "Cannot open density file $densityfile\n";
my $Header = {};

my $head = <DENSITYFILE>;
chomp ($head);
my @data = split;

my @data = split (/\t/, $head);
for (my $n= 0; $n <= $#data; $n++) {
  $Header->{$data[$n]} = $n;
}

my $countline = 0;
my $lastseg = 0;
my $lastchr;
my $lastdensity;
my  $lastcgall;
my $donechr = 0;

while (<DENSITYFILE>) {
  chomp;
  my @data = split;
  my $f1strain = $data[$Header->{'s1'}];
  my $s2strain = $data[$Header->{'s2'}];
  unless ($f1strain eq $firststrain && $s2strain eq $secondstrain) {
    next;
  }
  my $chr = $data[$Header->{'chr'}];
  
##############
# match the position to the intervals defined in the pairwise file
##############

  if ($chr eq $chr_to_use) {
    $donechr++;
    my $seg       =  $data[$Header->{'seg'}]; 
    my $density   =  $data[$Header->{'n'}];
    my $CGall     =  $data[$Header->{'CG.all'}];

    if ($lastseg == 0) {
      $Density->{$chr}{0}{nextseg} = $seg;
      #print "seg $seg\n";
      $Density->{$chr}{0}{density} = "NA";
      $lastseg = $seg;
    }
    if ($seg ne $lastseg) {   
      $Density->{$chr}{$lastseg}{nextseg} = $seg;
      #print "$chr $seg $lastseg \t$density\t$CGall\n";
      $Density->{$chr}{$lastseg}{'density'} = $density;
      $Density->{$chr}{$lastseg}{'CG.all'} = $CGall; 
    }
    $lastseg = $seg;
  }
  else {
    last unless ($donechr == 0);
  }
}
close (DENSITYFILE);

foreach my $chr (sort keys %{$Density}) {
  foreach my $seg (sort{$a<=> $b} keys %{$Density->{$chr}}) {
    my $nextseg =  $Density->{$chr}{$seg}{'nextseg'};
    my $density =  $Density->{$chr}{$seg}{'density'};
    my $CGall = $Density->{$chr}{$seg}{'CG.all'};
#print "$seg $nextseg, $density , $CGall\n";
    for (my $n = $seg; $n < $nextseg; $n+= $interval) {      
      $Density->{$chr}{$n}{'nextseg'} = $seg;
      $Density->{$chr}{$n}{'density'} = $density;
      $Density->{$chr}{$n}{'CG.all'} = $CGall;
 # print "$chr\t$n\t$seg\t$nextseg\t$density\n";
    }
  }
}


open (RNA, $rnafile) || die "Cannot open the rna file $rnafile\n";
#baseMean	log2FoldChange	lfcSE	stat	pvalue	padj	neglog10P	neglog10Padj	gene
my $Header = {};
my $head = <RNA>;
chomp ($head);
my @data = split;

my @data = split (/\t/, $head);
for (my $n= 0; $n <= $#data; $n++) {
  $Header->{$data[$n]} = $n;
}
while (<RNA>) {
  chomp;
  my @data = split;
  my $logp =  $data[$Header->{'neglog10Padj'}];
  my $log2fold =  $data[$Header->{'log2FoldChange'}];
  $RNA->{$data[$Header->{'gene'}]}{'logp'} = $logp;;
  $RNA->{$data[$Header->{'gene'}]}{'log2fold'} = $log2fold; 
}

close (RNA);

open (RNACOUNT, $rnacountfile) || die "Cannot open the rna file $rnacountfile\n";
#baseMean	log2FoldChange	lfcSE	stat	pvalue	padj	neglog10P	neglog10Padj	gene
my $Header = {};
my $head = <RNACOUNT>;
chomp ($head);
my @data = split;

my @data = split (/\t/, $head);
for (my $n= 0; $n <= $#data; $n++) {
  $Header->{$data[$n]} = $n;
}
while (<RNACOUNT>) {
  chomp;
  my @data = split;
  my $b6 = $data[$Header->{'b6.rna.counts'}];
  my $d2 = $data[$Header->{'d2.rna.counts'}];
  $RNACount->{$data[$Header->{'gene'}]}{'b6'} = $b6;
  $RNACount->{$data[$Header->{'gene'}]}{'d2'} = $d2; 
}

close (RNACOUNT);



open (ATAC, $atacfile) || die "Cannot open the atac file $atacfile\n";
#baseMean	log2FoldChange	lfcSE	stat	pvalue	padj	neglog10P	neglog10Padj	gene
my $Header = {};
my $head = <ATAC>;
chomp ($head);
my @data = split (/\t/, $head);
for (my $n= 0; $n <= $#data; $n++) {
  $Header->{$data[$n]} = $n;
}




#chr	start	end	peak	baseMean	log2FoldChange	lfcSE	stat	pvalue	padj	neglog10P	neglog10Padj
while (<ATAC>) {
  chomp;
  my @data = split;
  my $chr = $data[$Header->{'chr'}];
  next unless ($chr eq $chr_to_use);
  
  my $start =  $data[$Header->{'start'}];
  my $end =  $data[$Header->{'end'}];
  my $logp =  $data[$Header->{'neglog10Padj'}];
  my $log2fold =  $data[$Header->{'log2FoldChange'}];
  
  for (my $n = $start; $n <= $end; $n++) {
    $ATAC->{$chr}{$n}{logp} = $logp; 
    $ATAC->{$chr}{$n}{log2fold} =$log2fold;
  }
}
close (ATAC);

open (ATACCOUNT, $ataccountfile) || die "Cannot open the ATAC count file $ataccountfile\n";
#baseMean	log2FoldChange	lfcSE	stat	pvalue	padj	neglog10P	neglog10Padj	gene
my $Header = {};
my $head = <ATACCOUNT>;
chomp ($head);
my @data = split (/\t/, $head);
for (my $n= 0; $n <= $#data; $n++) {
  $Header->{$data[$n]} = $n;
}

#chr	start	end	peak	baseMean	log2FoldChange	lfcSE	stat	pvalue	padj	neglog10P	neglog10Padj
while (<ATACCOUNT>) {
  chomp;
  my @data = split;
  my $chr = $data[$Header->{'chr'}];
  next unless ($chr eq $chr_to_use);
  
  my $start =  $data[$Header->{'start'}];
  my $end =  $data[$Header->{'end'}];

  my $b6 = $data[$Header->{'b6.atac.counts'}];
  my $d2 = $data[$Header->{'d2.atac.counts'}];
  
  for (my $n = $start; $n <= $end; $n++) {
    $ATACCount->{$chr}{$n}{b6} = $b6; 
    $ATACCount->{$chr}{$n}{d2} = $d2;
  }
}
close (ATACCOUNT);



open (COUNTED, $countedfile ) || die "Cannot open counted file $countedfile file\n";

my $Header = {};
my $head = <COUNTED>;
chomp ($head);
my @data = split (/\t/, $head);
for (my $n= 0; $n <= $#data; $n++) {
  $Header->{$data[$n]} = $n;
}
print "$head\trna.logp\trna.log2fold\trna.b6.counts\trna.d2.count\tatac.logp\tatac.log2fold\tatac.b6.counts\tatac.d2.counts\tdensity\tCG.all\n";

#chr     pos     context b6.var  d2.var  b6.me   b6.cov  d2.me   d2.cov  annotations     cgi     tss     gene
while (<COUNTED>) {
  chomp;
  my @data = split;
  
  my $chr =  $data[$Header->{'chr'}];
  my $gene = $data[$Header->{'gene'}];
  my $pos = $data[$Header->{'pos'}];

   
  print "$_";
  if (exists ($RNA->{$gene})) {
    printf "\t%3.2f\t%3.2f", $RNA->{$gene}{'logp'}, $RNA->{$gene}{'log2fold'};
    print "\t$RNACount->{$gene}{'b6'}\t$RNACount->{$gene}{'d2'}";
  }
  else {
    print "\tNA\tNA\tNA\tNA";
  }
  if (exists ($ATAC->{$chr}{$pos})) {
    printf "\t%3.2f\t%3.2f", $ATAC->{$chr}{$pos}{'logp'}, $ATAC->{$chr}{$pos}{'log2fold'};
     print "\t$ATACCount->{$chr}{$pos}{'b6'}\t$ATACCount->{$chr}{$pos}{'d2'}";
  }
  else {
    print "\tNA\tNA\tNA\tNA";
  }
  
  ####################
  #
  ####################

   my $x = $pos % $interval;
    my $seg = $pos - $x;
 # print "\nseg $seg pos $pos\n";
  if (exists ($Density->{$chr}{$seg})) {
    print "\t$Density->{$chr}{$seg}{'density'}\t$Density->{$chr}{$seg}{'CG.all'}";
  }
  # if not could because there was a missing segment - try the one below
  else {
    my $nextseg = $seg - $interval;
      if (exists ($Density->{$chr}{$nextseg})) {
	print "\t$Density->{$chr}{$nextseg}{'density'}\t$Density->{$chr}{$nextseg}{'CG.all'}";
      }
    else {
      print "\tNA\tNA";
    }
  }
  print "\n";
}
