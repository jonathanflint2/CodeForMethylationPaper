#!/usr/bin/perl

use strict;
use Getopt::Std;


our($opt_h, $opt_t,$opt_f, $opt_z, $opt_b, $opt_F, $opt_s, $opt_v, $opt_c, $opt_T, $opt_C, $opt_V, $opt_H, $opt_M, $opt_m, $opt_S, $opt_i, $opt_A, $opt_a, $opt_H, $opt_Y);

getopts("hf:b:T:c:C:VvHS:s:Fi:Aa:Hm:M:Y");
&usage if $opt_h;
sub usage {
  die "PairwiseComparison.pl -f <input file>
and calculate pairwise correlations
-i segments (e.g.  <1000> )
-T <type to get > CG AT 

-A all pairs
-H add a header
-m minimum allowed fraction of methylation to count as methylated (> 0.1)
-M minimum read coverage (4)

-F use fractional methylation for correlational analysis rather than  binary (default)

input format  (note that there are columns for strain methylation  and strain coverage)

chr     pos     context aj.var  b6.var  balb.var cast.var d2.var  fvb.var pwk.var wsb.var 
aj.me   aj.cov  b6.me   b6.cov  balb.me balb.cov cast.me cast.cov d2.me   d2.cov  fvb.me  fvb.cov pwk.me  pwk.cov wsb.me  wsb.cov 
annotations


\n";
}

#################################################################
my $minimum = defined ($opt_m) ? $opt_m : 0.1;
my $minimumcoverage = defined ($opt_M)? $opt_M : 5;
#################################################################

sub mean {
    my ($array_ref) = @_;
    my $sum = 0;
    foreach my $value (@$array_ref) {
        $sum += $value;
    }
    my $mean = $sum / scalar(@$array_ref);
    return $mean;
}

sub sd {
    my ($array_ref, $mean) = @_;
    my $sum_sq_diff = 0;
    foreach my $value (@$array_ref) {
        $sum_sq_diff += ($value - $mean) ** 2;
    }
    my $sd = sqrt($sum_sq_diff / (scalar(@$array_ref) - 1));
    return $sd;
}

sub covariance {
    my ($array1_ref, $array2_ref, $mean1, $mean2) = @_;
    my $sum_product = 0;
    foreach my $i (0..$#{$array1_ref}) {
        $sum_product += ($array1_ref->[$i] - $mean1) * ($array2_ref->[$i] - $mean2);
    }
    my $covariance = $sum_product / (scalar(@$array1_ref) - 1);
    return $covariance;
}



# Define the sum function
sub sum {
  my ($array) = @_;
  my $sum = 0;
  foreach my $element (@$array) {
    $sum += $element;
  }
  return $sum;
}


my $Strains = {};
$Strains->{'aj'} = 1;
$Strains->{'b6'} = 2;
$Strains->{'balb'} = 3;
$Strains->{'cast'} = 4;
$Strains->{'d2'} = 5;
$Strains->{'fvb'} = 6;
$Strains->{'pwk'} = 7;
$Strains->{'wsb'} = 8;

my $PairsOfStrains = {};
foreach my $s1 (sort keys %{$Strains} ) {
  foreach my $s2 (sort keys %{$Strains} ) {
    ############
    # omit the next two lines to get all strain contributions
    ############
    next if ($s1 eq $s2);
    unless ($opt_A) {    
    next if (exists($PairsOfStrains->{$s2}{$s1}));
  }
    $PairsOfStrains->{$s1}{$s2}++;
  }
}
sub sum {
  my $total = 0;
  foreach my $number (@_) {
    $total += $number;
  }
  return $total;
}




sub prop_occurrence {
 my  ($array1_ref, $array2_ref) = @_;
 my $p = 1;
 my $q = 0;
  # Convert binary data to numeric data
  my @array1 = map { $_ ? 1 : 0 } @$array1_ref;
  my @array2 = map { $_ ? 1 : 0 } @$array2_ref;

  # Count the occurrences of p and q in array1
  my $count_p_array1 = grep { $_ == $p } @array1;
  my $count_q_array1 = grep { $_ == $q } @array1;

  # Count the occurrences of p and q in array2
  my $count_p_array2 = grep { $_ == $p } @array2;
  my $count_q_array2 = grep { $_ == $q } @array2;

  # Count the occurrences of p and q at the same site in both arrays
 my $count_pp = 0;
 my $count_qq = 0;
  for my $i (0..$#array1) {
    $count_pp++ if ($array1[$i] == $p && $array2[$i] == $p);
   $count_qq++ if ($array1[$i] == $q && $array2[$i] == $q);
  }

  # Calculate the proportions
  my $prop_p_array1 = $count_p_array1 / scalar(@array1); 
 my $prop_p_array2 = $count_p_array2 / scalar(@array2);
  my $prop_q_array1 = $count_q_array1 / scalar(@array1); 
  my $prop_q_array2 = $count_q_array2 / scalar(@array2);  
 my $prop_pp = $count_pp / scalar(@array1);
 my $prop_qq = $count_qq / scalar(@array1);
  return ($prop_p_array1,  $prop_p_array2,$prop_q_array1,$prop_q_array2, $prop_pp, $prop_qq);
}






#################################################################


my $typetoget = defined($opt_T)? $opt_T : "CG";
my $AllowedTypes = {};
$AllowedTypes->{'C'} = 'G';
$AllowedTypes->{'CG'} = 'G';
$AllowedTypes->{'A'} = 'T';
$AllowedTypes->{'AT'} = 'T';

unless (exists($AllowedTypes->{$typetoget})) {
  die "Allowed types to get are \"C\", \"CG\" Default is \"CG\"\n";
}

my $Var = {};
$Var->{"aj.var"}++;
$Var->{"b6.var"}++;
$Var->{"balb.var"}++;
$Var->{"cast.var"}++;
$Var->{"d2.var"}++;
$Var->{"fvb.var"}++;
$Var->{"pwk.var"}++;
$Var->{"wsb.var"}++;

my $Me = {};
$Me->{"aj.me"}++;
$Me->{"b6.me"}++;
$Me->{"balb.me"}++;
$Me->{"cast.me"}++;
$Me->{"d2.me"}++;
$Me->{"fvb.me"}++;
$Me->{"pwk.me"}++;
$Me->{"wsb.me"}++;

my $Sequences ={};
$Sequences->{'CG'} = 1;
$Sequences->{'C'} = 2;
$Sequences->{'G'} = 3;
$Sequences->{'A'} = 4;
$Sequences->{'T'} = 5;
$Sequences->{'NA'} = 6;

my $Methylation = {};
$Methylation->{'1'}= 1;
$Methylation->{'0'}= 2;
$Methylation->{'NA'}= 3;

my $Repeats = {};
$Repeats->{"SINE"}++;
$Repeats->{"LTR"}++; 
$Repeats->{"LINE"}++; 
$Repeats->{"LINE_L2"}++; 
$Repeats->{"LINE_L1"}++;
$Repeats->{"SINE_B2"}++; 
$Repeats->{"SINE_B4"}++; 
$Repeats->{"Simple_repeat"}++;

my $RegContext = {};
$RegContext->{"Enhancer"}++;
$RegContext->{"Promoter"}++;;
$RegContext->{"Het"}++;;        
$RegContext->{"Transcription"}++;;      
$RegContext->{'enhD'} = "enhD";
$RegContext->{'enhP'} = "enhP";
$RegContext->{'prom'} = "prom";
$RegContext->{'K4m3'} = "K4m3";
$RegContext->{'CTCF'} = "CTCF";
my $Genes = {};
$Genes->{'exon'}++;;
$Genes->{'intron'}++;
$Genes->{'fiveprime'}++;
$Genes->{'threeprime'}++;


#keep  arrays for storing h2, cgi and tss if a column heading is present in the file

my $Heritability  = {};
my $TSS = {};
my $CpGIslands = {};
my $Annotations = {};

my $Types = {};
$Types->{$typetoget}{'me.0'}++;
$Types->{$typetoget}{'me.1'}++;
$Types->{$typetoget}{'all'}++;
$Types->{$typetoget}{'all'}++;
$Types->{'mutant'}{'me.0'}++;
$Types->{'mutant'}{'me.1'}++;


$Types->{'s1'}{'me.0'}++;
$Types->{'s1'}{'me.1'}++;

$Types->{'s2'}{'me.0'}++;
$Types->{'s2'}{'me.1'}++;
$Types->{'both'}{'me.0'}++;
$Types->{'both'}{'me.1'}++;
$Types->{'one'}{'me.1'}++;
;
my $Correlations = {};
my $SeqCorrelations = {};

my $Results = {};
my $Queries = {};

my $Header = {};
open (FILE, $opt_f) || die "Cannot open file $opt_f\n";
my $head = <FILE>;
chomp($head);
my @data = split (/\t/, $head);
for (my $n = 0; $n <= $#data; $n++) {
  $Header->{$data[$n]} = $n;
}

while (<FILE>) {
  chomp;
  my @data = split;
  my $good = 0;
  my $varcount = 0;
  my $STypes = {};
  my $MTypes = {};
  my $pos = $data[$Header->{'pos'}];
  my $seg = 1;
  my $chr = $data[$Header->{'chr'}];
  
  if ($opt_i) {
    my $x = $pos % $opt_i;
    $seg = $pos - $x;
  }
  
  #$seg = $chr . "\t" . $seg;
  # print "seg $seg\n";
  
  my $Done = {};
  
  #####################    
  # Process annotations
  #####################
 
  #print "$pos $s1 $s2 $s1var  $s2var\n";
  
  foreach my $s1 (sort keys %{$PairsOfStrains}) {
    foreach my $s2 (sort keys %{$PairsOfStrains->{$s1}}) {    
      my $s1var = $s1 . ".var";
      my $s2var = $s2 . ".var";
      
      $s1var = $data[$Header->{$s1var}];
      $s2var = $data[$Header->{$s2var}];
      
      my $s1me = $s1 . ".me";
      my $s2me = $s2 . ".me";
      
      $s1me = $data[$Header->{$s1me}];
      $s2me = $data[$Header->{$s2me}];
      
      my $s1cov = $s1 . ".cov";
      my $s2cov = $s2 . ".cov";
      $s1cov = $data[$Header->{$s1cov}];
      $s2cov = $data[$Header->{$s2cov}];
      
      #   print "$s1 $s1me $s1cov $s2 $s2me $s2cov\n";
      #############################
      #skip if neither site is type to get - ie neither site is a CpG
      #############################
      next if ( $s1var ne $typetoget && $s2var ne  $typetoget);   
      # count the total number of $typtoget in each segment 
      if ($s1var eq $typetoget || $s2var eq $typetoget) {
        $Queries->{$s1}{$s2}{$chr}{$seg}++;
      }
      #############################
      #skip if either site has two few reads to call
      #############################
      
      next if ( $s1cov  < $minimumcoverage || $s2cov < $minimumcoverage);   
      
      #deal with annotations
      my   @annot = split (/,/, $data[$Header->{'annotations'}]);
      # check:
      #somehow this is missing simple repeats
      
      
      for (my $n = 0; $n<= $#annot; $n++) {
        next if ($annot[$n] eq "NA") ;
        $Annotations->{$chr}{$seg}{$s1}{$s2}{$annot[$n]}++
      }

       ####################
      # collect the heritability values, if present, and keep the mean and n for each segment -
      # looks like there are not more than than 6 for each 2Kb segment so keep them in an array
      # later can calculate a weighted average in the SumByDensity function
      ####################

      if (exists( $Header->{'h2'})) {
	if ( $data[$Header->{'h2'}] ne "NA") {
	  push (@{$Heritability->{$chr}{$seg}{$s1}{$s2}},  $data[$Header->{'h2'}]);
	  #	  print "@{$Heritability->{$chr}{$seg}{$s1}{$s2}}\n";
	}	
      }
      
      if (exists($Header->{'tss'})) {
	if ( $data[$Header->{'tss'}] ne "NA") {
	 # print "tss $data[$Header->{'tss'}]\n";
	  $TSS->{$chr}{$seg}{$s1}{$s2}++;
	}
      }
      if (exists($Header->{'cgi'})) {
	if ( $data[$Header->{'cgi'}] ne "NA") {
	 # print "cgi $data[$Header->{'cgi'}]\n";
	  $CpGIslands->{$chr}{$seg}{$s1}{$s2}++;
	}
      }
      
      
      ####################
      # Equalize the coverage
      ####################
      my $skipreduction = 0;
      ####################################
      # deal with missing data - could be due to a deletion , in which case skip the reduction,
      # or could be too low coverage - in which case delete this line
      # Note that the first possibility means there will be unequal coverage, due to the presence of the non-deleted strain
      ####################################      
      if ($s1me eq "NA") {
        if ($s1var eq "NA" || $s1var eq "N") {
          $skipreduction++;
          # try removing line...
        }
      }
      
      if ($s2me eq "NA") {
        if ($s2var eq "NA" || $s2var eq "N") {
          $skipreduction++;
          # try removing line...
        }
      }
      ######
      # equalization here
      ######
      unless ($skipreduction > 0) { 
        if ($s1cov > $s2cov) {
	  #################################################################
          #reduce methylated reads by probability that a read is methylated and re-calculate methylated reads
          #################################################################
          my $methylatedreads = $s1me;
	  my $s1fraction = $methylatedreads/$s1cov;
          my $diff = $s1cov - $s2cov;
          $s1cov = $s1cov - $diff;
	  # print "diff $diff and s1cov is now $s1cov\n";
          for ( my $i = 1; $i <= $diff; $i++) {
            my $removep = rand(1);
            if ($removep <= $s1fraction) {
              $methylatedreads--;            
            }
	    #    print "i = $i prob = $removep $methylatedreads\n";
          }
          if ($s1cov > 0) {
            $s1me =  $methylatedreads;
          }
          else {
            $s1cov = "NA";
            $s1me = "NA";
          }      
        }
        
        if ($s2cov > $s1cov) {
          #################################################################
          #reduce methylated reads by probability that a read is methylated and re-calculate methylated reads
          #################################################################
          my $methylatedreads = $s2me;
	  my $s2fraction = $methylatedreads/$s2cov;
          my $diff = $s2cov - $s1cov;
          $s2cov = $s2cov - $diff;
          # print "diff $diff\n";
          for ( my $i = 1; $i <= $diff; $i++) {
            my $removep = rand(1);
            if ($removep <= $s2fraction) {
              $methylatedreads--;            
            }
            # print "i = $i prob = $removep $methylatedreads\n";
          }
          if ($s2cov > 0) {
            $s2me =  $methylatedreads;
          }
          else {
            $s2cov = "NA";
            $s2me = "NA";
          }      
        }           
      }
      ############
      # end equalization
      ############
      
      #####################    
      # Process methylation
      #####################
      
      ###################
      #Decide whether a site is methylated or not, set it to 0 or 1,
      # set to NA if coverage is too low or if the seqence is N or NA
      # set to 0 if the coverage is good enough to call
      # set to 0 if the sequence is not C or CG
      ###################
      my $m1bin = 0;
      my $m2bin = 0;      
      my $s1fraction = 0;
      if ($s1cov ne "NA" && $s1cov > 0) {
        $s1fraction = $s1me/$s1cov;
      }
      my $s2fraction = 0;
      if ($s2cov ne "NA" && $s2cov > 0) {
        $s2fraction = $s2me/$s2cov;
      }
      
      if ($s1me eq "NA" ) {
	if ($s1var eq  "A" || $s2var eq "T") {
	  $m1bin = 0;
	}
	else {
	  $m1bin = 0;
	}
      }
      elsif ( $s1me >= $minimumcoverage && $s1fraction > $minimum ) {
        $m1bin = 1;
      }
      else {
	$m1bin = 0;
      }
      
      if ($s2me eq "NA") {
	if ($s2var eq  "A" || $s2var eq "T") {
	  $m2bin = 0;
	}
	else {
	  $m2bin = "NA";
	}
      }
      
      elsif ($s2me >= $minimumcoverage && $s2fraction > $minimum ) {
        $m2bin = 1;
      }
      
      else {
	$m2bin = 0 ;
      }
      #chr1	3020945	CCCGG	CG	CG	CG	CG	T	T	CG	CG	1	1	NA	NA	2	4	NA	NA
      # Keep fractional methylation for correlation analysis
      
      # print "$s1 $s1me $s1cov $m1bin $s2 $s2me $s2cov  $m2bin\n";
      
      if ($opt_F) {
	if ($m1bin ne "NA" && $m2bin ne "NA") {
	  push (@{$Correlations->{$chr}{$seg}{$s1}{$s2}{'s1'}}, $s1fraction);
	  push (@{$Correlations->{$chr}{$seg}{$s1}{$s2}{'s2'}}, $s2fraction);        
	}
      }
      
      # OR  keep binary methylation  for correlation analysis
      else {
	#deal with NAs - if sequence is A or T or G or NA  then set the mbin value to 0
	#only include if neither are "NA"
	if ($m1bin ne "NA" && $m2bin ne "NA") {
	  
	  push (@{$Correlations->{$chr}{$seg}{$s1}{$s2}{'s1'}}, $m1bin);	  
	  push (@{$Correlations->{$chr}{$seg}{$s1}{$s2}{'s2'}}, $m2bin);
	}
      }     
      ###################
      # Generate sequence comparison for the region,
      ###################
      
      if ($opt_Y) {
	#simple equivalence
	if ($s1var eq  $s2var ) {
	  push (@{$SeqCorrelations->{$chr}{$seg}{$s1}{$s2}{'s1'}}, 1);
	  push (@{$SeqCorrelations->{$chr}{$seg}{$s1}{$s2}{'s2'}}, 1);
	}
	else {
	  
	  push (@{$SeqCorrelations->{$chr}{$seg}{$s1}{$s2}{'s1'}}, 1);
	  push (@{$SeqCorrelations->{$chr}{$seg}{$s1}{$s2}{'s2'}}, 0);
	}
	
      }
      else {	
	if ($s1var eq "N" || $s2var eq "N" ) {
	   #Do nothing?
	}
	elsif ($s1var eq "C" && $s2var eq "C") {
	  push (@{$SeqCorrelations->{$chr}{$seg}{$s1}{$s2}{'s1'}}, 1);
	  push (@{$SeqCorrelations->{$chr}{$seg}{$s1}{$s2}{'s2'}}, 1);
	}
	elsif ($s1var eq "CG" && $s2var eq "CG") {
	  push (@{$SeqCorrelations->{$chr}{$seg}{$s1}{$s2}{'s1'}}, 1);
	  push (@{$SeqCorrelations->{$chr}{$seg}{$s1}{$s2}{'s2'}}, 1);
	}
	elsif ($s1var eq "C" && $s2var ne "C") {
	  push (@{$SeqCorrelations->{$chr}{$seg}{$s1}{$s2}{'s1'}}, 1);
	  push (@{$SeqCorrelations->{$chr}{$seg}{$s1}{$s2}{'s2'}}, 0);
	}
	elsif ($s2var eq "C" && $s1var ne "C") {
	  push (@{$SeqCorrelations->{$chr}{$seg}{$s1}{$s2}{'s1'}}, 0);
	  push (@{$SeqCorrelations->{$chr}{$seg}{$s1}{$s2}{'s2'}}, 1);
	}
	
	elsif ($s1var eq "CG" && $s2var ne "CG") {
	  push (@{$SeqCorrelations->{$chr}{$seg}{$s1}{$s2}{'s1'}}, 1);
	  push (@{$SeqCorrelations->{$chr}{$seg}{$s1}{$s2}{'s2'}}, 0);
	}
	elsif ($s2var eq "CG" && $s1var ne "CG") {
	  push (@{$SeqCorrelations->{$chr}{$seg}{$s1}{$s2}{'s1'}}, 0);
	  push (@{$SeqCorrelations->{$chr}{$seg}{$s1}{$s2}{'s2'}}, 1);
	}
	elsif ($s2var eq  $s1var ) {
	  push (@{$SeqCorrelations->{$chr}{$seg}{$s1}{$s2}{'s1'}}, 0);
	  push (@{$SeqCorrelations->{$chr}{$seg}{$s1}{$s2}{'s2'}}, 0);
	}
      }
      
      
      #how many sites are CpG?
      if ($s1var eq $typetoget) {
        $Results->{$chr}{$seg}{$s1}{$s2}{$typetoget}{'all'}++;
      }
      
      # how many of these sites are methylated?
      if ($s1var eq $typetoget && $m1bin == 1) {        
        $Results->{$chr}{$seg}{$s1}{$s2}{$typetoget}{'me.1'}++;
      }
      
      # how many of these sites are unmethylated?
      if ($s1var eq $typetoget && $m1bin == 0) {
        $Results->{$chr}{$seg}{$s1}{$s2}{$typetoget}{'me.0'}++;
      }
      # how many of the methylated CpG sites are mutated?
      if ($s1var eq $typetoget && $m1bin == 1 && $s2var ne $typetoget ) {
        $Results->{$chr}{$seg}{$s1}{$s2}{'mutant'}{'me.1'}++;
      }
      
      # how many of the unmethylated CpG sites are mutated?
      
      if ($s1var eq $typetoget && $m1bin == 0 && $s2var ne $typetoget ) {
        $Results->{$chr}{$seg}{$s1}{$s2}{'mutant'}{'me.0'}++;
      }
      
      #Look at the site that is mutant - 
      
      # if ($s1var eq $typetoget && $m1bin == 1 && $s2var ne $typetoget ) {
      #  $Results->{$chr}{$seg}{$s1}{$s2}{'mutant'}{'me.1'}++;
      #}
      
      #what is the proportion of sites in the segment that are methylated in the two strains?
      # keep a record for each strain and calculate the ratio at the end
      if ($s1var eq $typetoget && $m1bin == 1) {
        $Results->{$chr}{$seg}{$s1}{$s2}{'s1'}{'me.1'}++;
      }
      
      if ($s2var eq $typetoget && $m2bin == 1) {
        $Results->{$chr}{$seg}{$s1}{$s2}{'s2'}{'me.1'}++;
      }
      
      if ($s1var eq $typetoget && $m1bin == 0) {
        $Results->{$chr}{$seg}{$s1}{$s2}{'s1'}{'me.0'}++;
      }
      
      if ($s2var eq $typetoget && $m2bin == 0) {
        $Results->{$chr}{$seg}{$s1}{$s2}{'s2'}{'me.0'}++;
      }
      
      if ($s1var eq $typetoget && $m1bin == 1 && $m2bin == 1) {
        $Results->{$chr}{$seg}{$s1}{$s2}{'both'}{'me.1'}++;
      }
      
      if ($s1var eq $typetoget && $m1bin == 0 && $m2bin == 0) {
        $Results->{$chr}{$seg}{$s1}{$s2}{'both'}{'me.0'}++;
      }
      if ($s1var eq $typetoget) {
	my $onesite = 0;
        if ( $m1bin == 1) {
	  $onesite++;
	}
	if ( $m2bin == 1) {
	  $onesite++;
	}
	if ($onesite == 1) {
          $Results->{$chr}{$seg}{$s1}{$s2}{'one'}{'me.1'}++;
        }
      }
      $Done->{$s1}{$s2}++;
    }
  }
}


if ($opt_H) {
  print "chr\tseg\ts1\ts2\tn\tm1.m2\tme.ratio\tme.same\tseq1.seq2\tseq.ratio\tseq.same";
  foreach my $var (sort keys %{$Types}) {
    foreach my $me (sort keys %{$Types->{$var}}) {
      print "\t$var.$me";
    }
  }
  foreach my $gene (sort keys %{$Genes}) {
    print "\t$gene";
  }
  foreach my $reg (sort keys %{$RegContext}) {
    print "\t$reg";
  }
  foreach my $rep (sort keys %{$Repeats}) {
    print "\t$rep";
  }
  if (exists($Header->{'tss'})) {
    print "\ttss";
  }
  if (exists($Header->{'cgi'})) {
    print "\tcgi";
  }
  if (exists($Header->{'h2'})) {
    print "\tmean.h2\tn.h2";
  }
  
  print "\n";
}
my $Totals = {};

foreach my $chr (sort keys %{$Results}) {
  foreach my $seg (sort {$a <=> $b} keys %{$Results->{$chr}}) {
    foreach my $s1 (sort keys %{$Results->{$chr}{$seg}}) {
      foreach my $s2 (sort keys %{$Results->{$chr}{$seg}{$s1}}) {
	print "$chr\t$seg\t$s1\t$s2\t$Queries->{$s1}{$s2}{$chr}{$seg}";
	my @array1 = @{$Correlations->{$chr}{$seg}{$s1}{$s2}{'s1'}};
	my @array2 = @{$Correlations->{$chr}{$seg}{$s1}{$s2}{'s2'}};
	
	
	#####################
	# Estimate ratio of the two sites, where ratio is the number of sites that are methylated in one region
	# over the the number of sites that are methylated over a second region, and n1 is < n2
	#####################
	my $sum1 = sum(@array1);
	# my $size1 = scalar @array1;
	my $sum2 = sum(@array2);
	# my $size2 = scalar @array2;
	my $binratio = "NA";
	
	if ($sum1 == $sum2) {
	  $binratio = 1;
	}
	elsif ($sum1 >= $sum2 && $sum1 > 0) {	 
	  $binratio = $sum2/$sum1;
	}
	elsif ($sum2 >= $sum1 && $sum2 > 0) {
	  $binratio = $sum1/$sum2;
	}
	
	my $nratio = "$sum1:$sum2";
	# print "\n$sum1\t$size1\t$sum2\t$size2\n";
	if ($binratio eq "NA") {
	  print "\t$nratio\tNA";
	}
	else {
	  printf "\t$nratio\t%3.3f", $binratio;
	}
	#####################
	# Estimate the probability that the sites are equal (either 0, 0 or 1,1)
	#
	#####################
	my ($p_array1,  $p_array2,$q_array1,$q_array2, $prop_p, $prop_q) = prop_occurrence(\@array1, \@array2 );
	my $equal = $prop_p + $prop_q;
	printf "\t%3.4f", $equal;
	
	
	
	#####################
	# Estimate sequence divergence
	#####################
	if (defined($SeqCorrelations->{$chr}{$seg}{$s1}{$s2}{'s1'}) && defined ($SeqCorrelations->{$chr}{$seg}{$s1}{$s2}{'s2'})) {
	my @array1 = @{$SeqCorrelations->{$chr}{$seg}{$s1}{$s2}{'s1'}};
	my @array2 = @{$SeqCorrelations->{$chr}{$seg}{$s1}{$s2}{'s2'}};

	my $sum1 = sum(@array1);
	my $sum2 = sum(@array2);
	my $binratio = "NA";
	
	if ($sum1 == $sum2) {
	  $binratio = 1;
	}
	elsif ($sum1 >= $sum2 && $sum1 > 0) {	 
	  $binratio = $sum2/$sum1;
	}
	elsif ($sum2 >= $sum1 && $sum2 > 0) {
	  $binratio = $sum1/$sum2;
	}
	
	my $nratio = "$sum1:$sum2";
	if ($binratio eq "NA") {
	  print "\t$nratio\tNA";
	}
	else {
	  printf "\t$nratio\t%3.3f", $binratio;
	}
	my ($p_array1,  $p_array2,$q_array1,$q_array2, $prop_p, $prop_q) = prop_occurrence(\@array1, \@array2 );
	my $equal = $prop_p + $prop_q;
	printf "\t%3.4f", $equal;
      }
	else {
	  print "\tNA\tNA\tNA";
	}
	#################################
	
	
	foreach my $var (sort keys %{$Types}) {
	  foreach my $me (sort keys %{$Types->{$var}}) {
	    if (exists($Results->{$chr}{$seg}{$s1}{$s2}{$var}{$me})) {     
	      print "\t$Results->{$chr}{$seg}{$s1}{$s2}{$var}{$me}";
	    }
	    else {
	      print "\t0";
	    }                
	  }
	}
	foreach my $gene (sort keys %{$Genes}) {
	  if (exists ( $Annotations->{$chr}{$seg}{$s1}{$s2}{$gene})) {
	    print "\t$Annotations->{$chr}{$seg}{$s1}{$s2}{$gene}";
	  }
	  else {
	    print "\t0";
	  }
	}
	foreach my $reg (sort keys %{$RegContext}) {
	  
	  if (exists ( $Annotations->{$chr}{$seg}{$s1}{$s2}{$reg})) {
	    print "\t$Annotations->{$chr}{$seg}{$s1}{$s2}{$reg}";
	  }
	  else {
	    print "\t0";
	  }
	}
	foreach my $rep (sort keys %{$Repeats}) {
	  if (exists ( $Annotations->{$chr}{$seg}{$s1}{$s2}{$rep})) {
	    print "\t$Annotations->{$chr}{$seg}{$s1}{$s2}{$rep}";
	  }
	  else {
	    print "\t0";
	  }
	}
	
	if (exists( $Header->{'tss'})) {
	  if (defined ( $TSS->{$chr}{$seg}{$s1}{$s2})) {
	    print "\t$TSS->{$chr}{$seg}{$s1}{$s2}";
	  }
	  else {
	    print "\t0";
	  }
	}
	
	if (exists( $Header->{'cgi'})) {
	  if (defined ( $CpGIslands->{$chr}{$seg}{$s1}{$s2})) {
	    print "\t$CpGIslands->{$chr}{$seg}{$s1}{$s2}";
	  }
	  else {
	    print "\t0";
	  }
	}
	
	if (exists( $Header->{'h2'})) {
	  if (defined ( $Heritability->{$chr}{$seg}{$s1}{$s2})) {
	    my @h2array = @{$Heritability->{$chr}{$seg}{$s1}{$s2}};
	    my $nh2 = scalar(@h2array);
	    my $meanh2  = "NA";
	    if ($nh2 > 1) {
	      $meanh2 = mean (\@h2array);
	    }
	    else {
	      $meanh2 = $h2array[0];
	    }
	    #  print "\t";
	    #print join (",", @h2array);  
	    printf "\t%3.4f\t$nh2", $meanh2;
	  }
	  else {
	    print "\tNA\tNA";
	  }
	}
	
	print "\n";
	
      }
    }
  }  
}

