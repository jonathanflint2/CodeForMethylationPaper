#!/usr/bin/perl
use strict;
use Getopt::Std ;
use Compress::Zlib;
use List::Util qw(max);
use Carp qw< croak >;
use List::Util qw< sum >;
use POSIX;

our ($opt_h, $opt_f, $opt_t, $opt_o, $opt_m, $opt_z, $opt_s, $opt_S, $opt_v, $opt_H, $opt_M, $opt_d, $opt_M, $opt_z, $opt_l, $opt_R, $opt_D);
getopts("hf:t:z:sSvHM:m:lRD");
&usage if $opt_h;   
sub usage {

  die "IdentifVariableMethylation.pl -z  <gzipped file> (this is assumed to be the composite methylation file, not the annotated one)
or -f for text file

-M exclude entries with more than a total coverage of M  (default is 500)
-m exclude entries with less than a total coverage of m  (default is 1)
-H file contains a header
-t threshold for number of reads to call for difference analysis (default is 5)
-z open gzipped file
-l just add logp to each line
-R process replicates file

-D only d2 and b6 strains

\n";
}
my $threshold  = defined ($opt_t) ? $opt_t : 5;
&usage unless defined($opt_f) || defined ($opt_z);
my $minimum = defined ($opt_m) ? $opt_m : 1;
my $maximum = defined ($opt_M) ? $opt_M : 500;
sub round {
  my $number = shift || 0;
  my $dec = 10 ** (shift || 0);
  return int( $dec * $number + .5 * ($number <=> 0)) / $dec;
}
my $SeqContext = {};
$SeqContext->{'CGA'} = "CG";
$SeqContext->{'CGC'} = "CG";
$SeqContext->{'CGG'} = "CG";
$SeqContext->{'CGT'} = "CG";
$SeqContext->{'CAA'} = "CH";
$SeqContext->{'CAC'} = "CH";
$SeqContext->{'CAG'} = "CHG";
$SeqContext->{'CAT'} = "CH";
$SeqContext->{'CCA'} = "CH";
$SeqContext->{'CCC'} = "CH";
$SeqContext->{'CCG'} = "CHG";
$SeqContext->{'CCT'} = "CH";
$SeqContext->{'CTA'} = "CH";
$SeqContext->{'CTC'} = "CH";
$SeqContext->{'CTG'} = "CHG";
$SeqContext->{'CTT'} = "CH";
$SeqContext->{'CAN'} = "CH";
$SeqContext->{'CCN'} = "CH";
$SeqContext->{'CGN'} = "CG";
$SeqContext->{'CTN'} = "CH";
$SeqContext->{'CNA'} = "CN";
$SeqContext->{'CNC'} = "CN";
$SeqContext->{'CNG'} = "CHG";
$SeqContext->{'CNT'} = "CN";
$SeqContext->{'CNN'} = "CN";

use constant PI => 3.1415926536;
use constant SIGNIFICANT => 5; # number of significant digits to be returned

my $lastchr = "chr1";
my $lastpos = 0;
my $lastline;
my $count = 0;
my $pairmark = "P";
my $paircount;

my $Strains = {};

my $GoodStrains = {};
my $lastpos = 0;
my $found = 0;
my $file_handle;

#my $StrainLookUp = {};

#$StrainLookUp->{'aj'}= "aj";
#$StrainLookUp->{'balb'}= "balb";
#$StrainLookUp->{'blk6'} = "b6";
#$StrainLookUp->{'cast'} = "cast";
#$StrainLookUp->{'dba'} = "d2";
#$StrainLookUp->{'fvb'} = "fvb";
#$StrainLookUp->{'pwk'} = "pwk";
#$StrainLookUp->{'wsb'} = "wsb";
if ($opt_R) {
  $Strains->{'b61'}{'mc'}++;
  $Strains->{'d21'}{'mc'}++;
  $Strains->{'cast1'}{'mc'}++;
  $Strains->{'fvb1'}{'mc'}++;
  $Strains->{'aj1'}{'mc'}++;
  $Strains->{'balb1'}{'mc'}++;;
  $Strains->{'pwk1'}{'mc'} ++;
  $Strains->{'wsb1'}{'mc'}++;
  $Strains->{'b62'}{'mc'}++;
  $Strains->{'d22'}{'mc'}++;;
  $Strains->{'cast2'}{'mc'}++;
  $Strains->{'fvb2'}{'mc'}++;
  $Strains->{'aj2'}{'mc'}++;
  $Strains->{'balb2'}{'mc'}++;;
  $Strains->{'pwk2'}{'mc'}++;;
  $Strains->{'wsb2'}{'mc'}++;;
  
}
#elsif ($opt_D) {
# $Strains->{'b6'}{'mc'}++;
#  $Strains->{'d2'}{'mc'} ++;
#}
else {
  $Strains->{'b6'}{'mc'}++;
  $Strains->{'d2'}{'mc'} ++;
  $Strains->{'cast'}{'mc'} ++;
  $Strains->{'fvb'}{'mc'}++; 
  $Strains->{'aj'}{'mc'} ++;
  $Strains->{'balb'}{'mc'}++;
  $Strains->{'pwk'}{'mc'} ++;
  $Strains->{'wsb'}{'mc'} ++;
}


if ($opt_z) {
open  $file_handle, "gzip -cd $opt_z |";
}
elsif ($opt_f)  {
  open ( $file_handle, $opt_f) || die "Can't open file $opt_f\n";
}
else {
  &usage;
}

if ($opt_H) {
  my $head = <$file_handle>;
  chomp ($head);
  if ($opt_m) {
    print "mark\t";
  }
  print "$head";
  
  if ($opt_l) {
    print "\tn.logp\n";
  }
  else {
    print "\tstrand\ttype\ttotal.mc\ttotal.cov\tfraction.mc\tn.strains\tlogp\tz.diff\tpair\n";
  }
}
else {
  print "chr\tpos";
  foreach my $strain (sort {$a cmp $b} keys %{$Strains}) {
    print "\t$strain.mc\t$strain.cov";
  }
  print "\tstrand\ttype\ttotal.mc\ttotal.cov\tfraction.mc\tn.strains\tlogp\tz.diff\tpair\n";
}

while (<$file_handle>) {
 #  chr1	3000173	aj	NA	NA	NA	NA	balb	+	CTT	1	4	blk6	+
  chomp;
  my ($chr,$pos, $ajstrain, $ajstrand, $ajseq, $ajmc, $ajcov, $balbstrain, $balbstrand, $balbseq, $balbmc, $balbcov,  $b6strain, $b6strand, $b6seq, $b6mc, $b6cov, $caststrain, $caststrand, $castseq, $castmc,$castcov, $d2strain, $d2strand, $d2seq, $d2mc,  $d2cov, $fvbstrain, $fvbstrand, $fvbseq, $fvbmc,$fvbcov, $pwkstrain, $pwkstrand, $pwkseq, $pwkmc, $pwkcov, $wsbstrain, $wsbstrand, $wsbseq, $wsbmc, $wsbcov, @rest) = split;
   my $KeepLine = $_; 
  
  if ($chr ne $lastchr) {
    $lastchr = $chr;
    next;
  }
  
  $paircount = $chr . "_" . $pos;

  if ($opt_D) {
  
     ($chr,$pos,  $b6strain, $b6strand, $b6seq, $b6mc, $b6cov, $d2strain, $d2strand, $d2seq, $d2mc,  $d2cov,  @rest) = split;
   }
  
  if ($opt_l) {
   # my ($chr,$pos, $ajmc,  $ajcov, $b6mc, $b6cov, $balbmc, $balbcov, $castmc, $castcov, $d2mc,	$d2cov,	$fvbmc,	$fvbcov,$pwkmc,	$pwkcov, $wsbmc, $wsbcov, $strand, $type, $totalmc, $totalcov,	$fractionmc,	$nstrains, $logp, $zdiff, $pair,  $gene,  $feature, $varpos, $varsdp,	$ccre,	$atac) = split;
     ($chr,$pos, $ajmc,  $ajcov, $b6mc, $b6cov, $balbmc, $balbcov, $castmc, $castcov, $d2mc,	$d2cov,	$fvbmc,	$fvbcov,$pwkmc,	$pwkcov, $wsbmc, $wsbcov ,@rest) = split;
   }

  if ($opt_R) {
    
    my  ($chr, $pos, $aj1strain, $aj1strand, $aj1seq, $aj1mc, $aj1cov,  $aj2strain, $aj2strand, $aj2seq, $aj2mc, $aj2cov,$balb1strain, $balb1strand, $balb1seq, $balb1mc, $balb1cov, $balb2strain, $balb2strand, $balb2seq, $balb2mc, $balb2cov, $b61strain, $b61strand, $b61seq, $b61mc, $b61cov, $b62strain, $b62strand, $b62seq, $b62mc, $b62cov, $cast1strain, $cast1strand, $cast1seq, $cast1mc,$cast1cov, $cast2strain, $cast2strand, $cast2seq, $cast2mc,$cast2cov,$d21strain, $d21strand, $d21seq, $d21mc,  $d21cov,$d22strain, $d22strand, $d22seq, $d22mc,  $d22cov,  $fvb1strain, $fvb1strand, $fvb1seq, $fvb1mc,$fvb1cov,$fvb2strain, $fvb2strand, $fvb2seq, $fvb2mc,$fvb2cov, $pwk1strain, $pwk1strand, $pwk1seq, $pwk1mc, $pwk1cov, $pwk2strain, $pwk2strand, $pwk2seq, $pwk2mc, $pwk2cov,  $wsb1strain, $wsb1strand, $wsb1seq, $wsb1mc, $wsb1cov,$wsb2strain, $wsb2strand, $wsb2seq, $wsb2mc, $wsb2cov, @rest) = split;

  $Strains->{'b61'}{'mc'} = $b61mc;
  $Strains->{'d21'}{'mc'} = $d21mc;
  $Strains->{'cast1'}{'mc'} = $cast1mc;
  $Strains->{'fvb1'}{'mc'} = $fvb1mc; 
  $Strains->{'aj1'}{'mc'} = $aj1mc;
  $Strains->{'balb1'}{'mc'} = $balb1mc;
  $Strains->{'pwk1'}{'mc'} = $pwk1mc;
  $Strains->{'wsb1'}{'mc'} = $wsb1mc;


      $Strains->{'b62'}{'mc'} = $b62mc;
  $Strains->{'d22'}{'mc'} = $d22mc;
  $Strains->{'cast2'}{'mc'} = $cast2mc;
  $Strains->{'fvb2'}{'mc'} = $fvb2mc; 
  $Strains->{'aj2'}{'mc'} = $aj2mc;
  $Strains->{'balb2'}{'mc'} = $balb2mc;
  $Strains->{'pwk2'}{'mc'} = $pwk2mc;
    $Strains->{'wsb2'}{'mc'} = $wsb2mc;

    
  $Strains->{'b61'}{'nonmc'}   =  $b61cov - $b61mc;
  $Strains->{'d21'}{'nonmc'}   =  $d21cov - $d21mc;
  $Strains->{'cast1'}{'nonmc'} = $cast1cov- $cast1mc;
  $Strains->{'fvb1'}{'nonmc'}  = $fvb1cov - $fvb1mc; 
  $Strains->{'aj1'}{'nonmc'}   =  $aj1cov - $aj1mc;
  $Strains->{'balb1'}{'nonmc'} = $balb1cov - $balb1mc;
  $Strains->{'pwk1'}{'nonmc'}  = $pwk1cov - $pwk1mc;
    $Strains->{'wsb1'}{'nonmc'}  = $wsb1cov - $wsb1mc;


  $Strains->{'b62'}{'nonmc'}   =  $b62cov - $b62mc;
  $Strains->{'d22'}{'nonmc'}   =  $d22cov - $d22mc;
  $Strains->{'cast2'}{'nonmc'} = $cast2cov- $cast2mc;
  $Strains->{'fvb2'}{'nonmc'}  = $fvb2cov - $fvb2mc; 
  $Strains->{'aj2'}{'nonmc'}   =  $aj2cov - $aj2mc;
  $Strains->{'balb2'}{'nonmc'} = $balb2cov - $balb2mc;
  $Strains->{'pwk2'}{'nonmc'}  = $pwk2cov - $pwk2mc;
    $Strains->{'wsb2'}{'nonmc'}  = $wsb2cov - $wsb2mc;
    
  
  $Strains->{'b61'}{'cov'} = $b61cov;
  $Strains->{'d21'}{'cov'} = $d21cov;
  $Strains->{'cast1'}{'cov'} = $cast1cov;
  $Strains->{'fvb1'}{'cov'} = $fvb1cov; 
  $Strains->{'aj1'}{'cov'} = $aj1cov;
  $Strains->{'balb1'}{'cov'} = $balb1cov;
  $Strains->{'pwk1'}{'cov'} = $pwk1cov;
  $Strains->{'wsb1'}{'cov'} = $wsb1cov;

  $Strains->{'b62'}{'cov'} = $b62cov;
  $Strains->{'d22'}{'cov'} = $d22cov;
  $Strains->{'cast2'}{'cov'} = $cast2cov;
  $Strains->{'fvb2'}{'cov'} = $fvb2cov; 
  $Strains->{'aj2'}{'cov'} = $aj2cov;
  $Strains->{'balb2'}{'cov'} = $balb2cov;
  $Strains->{'pwk2'}{'cov'} = $pwk2cov;
  $Strains->{'wsb2'}{'cov'} = $wsb2cov;

    
  $Strains->{'b61'}{'strand'} = $b61strand;
  $Strains->{'d21'}{'strand'} = $d21strand;
  $Strains->{'cast1'}{'strand'} = $cast1strand;
  $Strains->{'fvb1'}{'strand'} = $fvb1strand; 
  $Strains->{'aj1'}{'strand'} = $aj1strand;
  $Strains->{'balb1'}{'strand'} = $balb1strand;
  $Strains->{'pwk1'}{'strand'} = $pwk1strand;
    $Strains->{'wsb1'}{'strand'} = $wsb1strand;



    $Strains->{'b62'}{'strand'} = $b62strand;
  $Strains->{'d22'}{'strand'} = $d22strand;
  $Strains->{'cast2'}{'strand'} = $cast2strand;
  $Strains->{'fvb2'}{'strand'} = $fvb2strand; 
  $Strains->{'aj2'}{'strand'} = $aj2strand;
  $Strains->{'balb2'}{'strand'} = $balb2strand;
  $Strains->{'pwk2'}{'strand'} = $pwk2strand;
  $Strains->{'wsb2'}{'strand'} = $wsb2strand;

  $Strains->{'b61'}{'seq'} = $b61seq;
  $Strains->{'d21'}{'seq'} = $d21seq;
  $Strains->{'cast1'}{'seq'} = $cast1seq;
  $Strains->{'fvb1'}{'seq'} = $fvb1seq; 
  $Strains->{'aj1'}{'seq'} = $aj1seq;
  $Strains->{'balb1'}{'seq'} = $balb1seq;
  $Strains->{'pwk1'}{'seq'} = $pwk1seq;
  $Strains->{'wsb1'}{'seq'} = $wsb1seq;

 $Strains->{'b62'}{'seq'} = $b62seq;
  $Strains->{'d22'}{'seq'} = $d22seq;
  $Strains->{'cast2'}{'seq'} = $cast2seq;
  $Strains->{'fvb2'}{'seq'} = $fvb2seq; 
  $Strains->{'aj2'}{'seq'} = $aj2seq;
  $Strains->{'balb2'}{'seq'} = $balb2seq;
  $Strains->{'pwk2'}{'seq'} = $pwk2seq;
  $Strains->{'wsb2'}{'seq'} = $wsb2seq;
    
  }
  
  elsif ($opt_D) {
    
    $Strains->{'b6'}{'mc'} = $b6mc;
    $Strains->{'d2'}{'mc'} = $d2mc;
    $Strains->{'b6'}{'nonmc'}   =  $b6cov - $b6mc;
    $Strains->{'d2'}{'nonmc'}   =  $d2cov - $d2mc;
    
    $Strains->{'b6'}{'cov'} = $b6cov;
    $Strains->{'d2'}{'cov'} = $d2cov;
    $Strains->{'b6'}{'strand'} = $b6strand;
    $Strains->{'d2'}{'strand'} = $d2strand;
    $Strains->{'b6'}{'seq'} = $b6seq;
    $Strains->{'d2'}{'seq'} = $d2seq;
  }
  else {

      $Strains->{'b6'}{'mc'} = $b6mc;
  $Strains->{'d2'}{'mc'} = $d2mc;
  $Strains->{'cast'}{'mc'} = $castmc;
  $Strains->{'fvb'}{'mc'} = $fvbmc; 
  $Strains->{'aj'}{'mc'} = $ajmc;
  $Strains->{'balb'}{'mc'} = $balbmc;
  $Strains->{'pwk'}{'mc'} = $pwkmc;
  $Strains->{'wsb'}{'mc'} = $wsbmc;
  
  $Strains->{'b6'}{'nonmc'}   =  $b6cov - $b6mc;
  $Strains->{'d2'}{'nonmc'}   =  $d2cov - $d2mc;
  $Strains->{'cast'}{'nonmc'} = $castcov- $castmc;
  $Strains->{'fvb'}{'nonmc'}  = $fvbcov - $fvbmc; 
  $Strains->{'aj'}{'nonmc'}   =  $ajcov - $ajmc;
  $Strains->{'balb'}{'nonmc'} = $balbcov - $balbmc;
  $Strains->{'pwk'}{'nonmc'}  = $pwkcov - $pwkmc;
  $Strains->{'wsb'}{'nonmc'}  = $wsbcov - $wsbmc;
  
  $Strains->{'b6'}{'cov'} = $b6cov;
  $Strains->{'d2'}{'cov'} = $d2cov;
  $Strains->{'cast'}{'cov'} = $castcov;
  $Strains->{'fvb'}{'cov'} = $fvbcov; 
  $Strains->{'aj'}{'cov'} = $ajcov;
  $Strains->{'balb'}{'cov'} = $balbcov;
  $Strains->{'pwk'}{'cov'} = $pwkcov;
  $Strains->{'wsb'}{'cov'} = $wsbcov;

  $Strains->{'b6'}{'strand'} = $b6strand;
  $Strains->{'d2'}{'strand'} = $d2strand;
  $Strains->{'cast'}{'strand'} = $caststrand;
  $Strains->{'fvb'}{'strand'} = $fvbstrand; 
  $Strains->{'aj'}{'strand'} = $ajstrand;
  $Strains->{'balb'}{'strand'} = $balbstrand;
  $Strains->{'pwk'}{'strand'} = $pwkstrand;
  $Strains->{'wsb'}{'strand'} = $wsbstrand;

  $Strains->{'b6'}{'seq'} = $b6seq;
  $Strains->{'d2'}{'seq'} = $d2seq;
  $Strains->{'cast'}{'seq'} = $castseq;
  $Strains->{'fvb'}{'seq'} = $fvbseq; 
  $Strains->{'aj'}{'seq'} = $ajseq;
  $Strains->{'balb'}{'seq'} = $balbseq;
  $Strains->{'pwk'}{'seq'} = $pwkseq;
      $Strains->{'wsb'}{'seq'} = $wsbseq;
    }
  


 my $Orientation = {};
  foreach my $strain (keys %{$Strains}) {
    next if ($Strains->{$strain}{'strand'} eq "NA");
     next if ($Strains->{$strain}{'strand'} eq '');
    $Orientation->{$Strains->{$strain}{'strand'}}++;
  }
  
  my @strands;
  my $highest = max values %{$Orientation};  
  my $o;
 
  foreach my $strand (keys %{$Orientation}) {
    if ($Orientation->{$strand} == $highest) {
      $o = $strand;
    }
  }

  my $strand = $o;
    
  my $Total = 0;
  my $Meth = 0;
  my $straincount = 0;
  my $Types = {};
 
  foreach my $strain (keys %{$Strains}) {
    if ($Strains->{$strain}{'cov'} eq '') {
      $Strains->{$strain}{'mc'} = "NA";
      $Strains->{$strain}{'cov'} = "NA";
      $Strains->{$strain}{'nonmc'} = "NA";
      next;
    }
    next if ($Strains->{$strain}{'mc'} eq "NA");
    next if ($Strains->{$strain}{'cov'} eq "NA");
  
    $Total+=$Strains->{$strain}{'cov'};
   $Meth += $Strains->{$strain}{'mc'};
    if (exists($SeqContext->{$Strains->{$strain}{'seq'}})) {
      $Types->{$SeqContext->{$Strains->{$strain}{'seq'}}}++;
    }
    $straincount++;
  }
  next if ($Total <= $minimum);
  next if ($Total >= $maximum);

  my $type;
  my $typecount = 0;
  foreach my $seqtype (keys %{$Types}) {
    $typecount++;
    # print "$type\t$Types->{$t}";
    $type = $seqtype;
  }
  if ($typecount > 1) {
    $type = "MIXED";
  }
   my $line=  "$chr\t$pos";
 foreach my $strain (sort {$a cmp $b} keys %{$Strains}) {
    if ($Strains->{$strain}{'cov'} eq '') {
        $line .= "\tNA\tNA";
    }
    else {
    $line .= "\t$Strains->{$strain}{'mc'}\t$Strains->{$strain}{'cov'}";
  }
  }
  my $Diff = {};
   $Diff =  DifferentialSites($Strains, $threshold);
  
  my $logp = round($Diff->{'logp'}, 3);
  my $a = round($Diff->{'a'}, 3);
 # tp.diff\tz.diff\tr.diff\tall.mc\t
  $line .= "\t$strand\t$type\t$Meth\t$Total\t$a\t$straincount\t$logp\t$Diff->{'z'}";

  my $diff = $pos - $lastpos;
 # print "$pos, $strand, $diff, $found\n";
  if ($opt_l) {
 print "$KeepLine\t$logp\n";
  }
  else {  
  if ($found == 1 ) {
    my $diff = $pos - $lastpos; 
    
    if ($strand eq "-") {
      if ($diff == 1) {
	my $pair = $pairmark . $paircount;
	print "$lastline\t$pair\n$line\t$pair\n";
      }
      else {
	print "$lastline\tNA\n$line\tNA\n";
      }
    }
    else {
      print "$lastline\tNA\n";
      $lastline = $line;
      $found = 1;
      $lastpos = $pos;
      next;
    }
    $found = 3;
    next;
  }
  else {
    #	print "$lastline\tNA\n$line\tNA\n";
  }
  
  if ($strand eq "+") {
    $found = 1;
  }
  else {
    print "$line\tNA\n";
    # my $check = length($line);
    # printf "check $check\n";
    # if ($check > 60 ) {
    #   print "$line\tNA\n";
    # }
    
  }
  
}
  $lastpos = $pos;
  $lastline = $line;
  
}


#############################

sub log10 {
  my $n = shift;
  if ($n > 0) {
    return log($n) / log(10);
  }
  else {
    return (0)
  }
}

sub chi_squared_test {
  my %args = @_;
  my $observed = delete $args{observed};
  my $expected = delete $args{expected};
  @$observed == @$expected || die q(Input arrays must have same length);

 # print "O @$observed o $observed->[$_] E  @$expected e $expected->[$_] \n";
  
    my $chi_squared = sum map {
      ($observed->[$_] - $expected->[$_])**2 / $expected->[$_];
    } 0 .. $#$observed;

# CHECK THIS
  my $degrees_of_freedom = @$observed/2 - 1;
  
# print "df $degrees_of_freedom  O $#{$observed},  E $#{$expected}\n";
     # degrees of freedom is equal to (r-1)(c-1),
   
  my $probability = chisqrprob($degrees_of_freedom, $chi_squared);
  # print "$degrees_of_freedom, $chi_squared\t$probability\n";
  return $probability;
  
  
}
sub chisqrprob { # Upper probability   X^2(x^2,n)
	my ($n,$x) = @_;
	if (($n <= 0) || ((abs($n) - (abs(int($n)))) != 0)) {
		die "Invalid n: $n\n"; # degree of freedom
	}
	return precision_string(_subchisqrprob($n, $x));
}

sub precision_string {
	my ($x) = @_;
	if ($x) {
		return sprintf "%." . precision($x) . "f", $x;
	} else {
		return "0";
	}
}

sub _subchisqrprob {
	my ($n,$x) = @_;
	my $p;

	if ($x <= 0) {
		$p = 1;
	} elsif ($n > 100) {
		$p = _subuprob((($x / $n) ** (1/3)
				- (1 - 2/9/$n)) / sqrt(2/9/$n));
	} elsif ($x > 400) {
		$p = 0;
	} else {   
		my ($a, $i, $i1);
		if (($n % 2) != 0) {
			$p = 2 * _subuprob(sqrt($x));
			$a = sqrt(2/PI) * exp(-$x/2) / sqrt($x);
			$i1 = 1;
		} else {
			$p = $a = exp(-$x/2);
			$i1 = 2;
		}

		for ($i = $i1; $i <= ($n-2); $i += 2) {
			$a *= $x / $i;
			$p += $a;
		}
	}
	return $p;
      }

sub _subuprob {
	my ($x) = @_;
	my $p = 0; # if ($absx > 100)
	my $absx = abs($x);

	if ($absx < 1.9) {
		$p = (1 +
			$absx * (.049867347
			  + $absx * (.0211410061
			  	+ $absx * (.0032776263
				  + $absx * (.0000380036
					+ $absx * (.0000488906
					  + $absx * .000005383)))))) ** -16/2;
	} elsif ($absx <= 100) {
		for (my $i = 18; $i >= 1; $i--) {
			$p = $i / ($absx + $p);
		}
		$p = exp(-.5 * $absx * $absx) 
			/ sqrt(2 * PI) / ($absx + $p);
	}

	$p = 1 - $p if ($x<0);
	return $p;
}
sub precision {
	my ($x) = @_;
	return abs int(log10(abs $x) - SIGNIFICANT);
}

sub DifferentialSites {
  my $hash_ref  = shift;
  my $threshold = shift;
 
  my %Strains = %{ $hash_ref }; 
  my $Out = {};
  my $Rows = {};
  my $Columns = {};
  my $Observed = {};
  my $Expected = {};
  my $row =0;
  my $c1;
  my $total = 0;
    my $logp = "NA";
  my $p = 0;
  my $z = 0;
  my $y = 0;
  my $a = 0;
   my $fullymethylated = 0;
  foreach my $strain (sort {$a cmp $b} keys %{$Strains}) {
   
    
    unless ($Strains->{$strain}{'cov'} eq "NA") {
      $row++;
      ############################
      # Code below makes 0 into 0.5 so that the chi square is a bit more balanced but it messes up the zdiff scores! so i have disabled it
      ############################
      
    #  if ($Strains->{$strain}{'mc'} == 0) {
#	$Strains->{$strain}{'mc'} = 0.5;
 #     }
  #    if ($Strains->{$strain}{'nonmc'} == 0) {
#	$Strains->{$strain}{'nonmc'} = 0.5;
 #     }
      $Observed->{$row}{c1} = $Strains->{$strain}{'mc'};
      $Observed->{$row}{c2} = $Strains->{$strain}{'nonmc'};
      $Rows->{$row} =  $Strains->{$strain}{'mc'} + $Strains->{$strain}{'nonmc'};
      $Columns->{c1} += $Strains->{$strain}{'mc'};
      $Columns->{c2} += $Strains->{$strain}{'nonmc'};
      $total+= $Strains->{$strain}{'mc'} + $Strains->{$strain}{'nonmc'};
      if ( $Strains->{$strain}{'mc'}  == $Strains->{$strain}{'cov'}) {
	$fullymethylated++;
      }    
    }
   
   #  print "$strain\t$Strains->{$strain}{'mc'}\t$Strains->{$strain}{'nonmc'}\t$fullymethylated\t$row\n";
  }
  my @obs;
  my @exp;
  foreach my $row (sort keys %{$Observed}) {
    push (@exp,   $Rows->{$row} *  $Columns->{c1} / $total);
    push (@obs,  $Observed->{$row}{c1});  
    push (@exp,   $Rows->{$row} *  $Columns->{c2} / $total);
    push (@obs,  $Observed->{$row}{c2});
    #$Expected->{$row}{c1} =  $Rows->{$row} *  $Columns->{c1} / $total;
    # $Expected->{$row}{c2} =  $Rows->{$row} *  $Columns->{c2} / $total;
    }

 
  if ($fullymethylated == $row) {
    $logp = 0;
  }
  elsif ($row > 1) {
     
    my $pval = chi_squared_test ( observed => \@obs,expected => \@exp);
    if ($pval > 0) {
      $logp = -1 * log10($pval);
    }
    else {
      $logp = "NA";
    }
  }
  
 
  
  ### r differential and partly methylated = at least one strain is completely methylated and at least one is partly methylated
  ### p partly methylated = at least one strain is partly methylated
  ### -> z differential   = at least one strain is completely methylated and at least one is unmethylated
  ### y differential   = at least one strain is completely unmethylated and at least one is partly methylated
  ### -> a  report proportion of methylation to total coverage

  my $ncount = 0;
  my $strainmc = 0;
  my $strainpartial = 0;
  my $strainzero = 0;
  my $straincomplete = 0;
  my $allstrainmethylated = 0;
 my $anystrainpartial = 0;
  my $mctotalcov = 0;
  my $totalcov = 0;
  
  foreach my $strain (sort keys %{$Strains}) {
    my $mc = $Strains->{$strain}{'mc'};
    my $cov = $Strains->{$strain}{'cov'};
    next if ($mc eq "NA" || $cov eq "NA");
    $ncount++;
    if ($mc == $cov) {
      $allstrainmethylated ++;
    }
    
    $mctotalcov += $mc;
    $totalcov+= $cov;
    
   
     if ($mc < $cov) {
	$anystrainpartial++;
      }
    if ($cov >= $threshold) {
      
      if ($mc == $cov) {
	$strainmc++;
      }
     if ($mc == 0 ) {
	$strainzero++;
      }
      if ($mc < $cov) {
	$strainpartial++;
      }
    }
    #print " $strain\t$mc\t$cov ";
  }
  
  if ($strainmc > 0 && 	$strainpartial > 0 && $strainzero == 0) {
    $p = 1;
  }
  if ($anystrainpartial > 0 && $strainzero > 0) {
    $y = 1;
  }
  if ($strainzero > 0 && $strainmc > 0) {
    $z = 1;
  }
  
  if ($totalcov >0) {
    $a =  $mctotalcov/$totalcov;
   
  }
    
 # print " $a, $p, $y, $z, $logp\n";
  
  $Out->{'p'} = $p;
  $Out->{'z'} = $z;
  $Out->{'y'} = $y;
  $Out->{'logp'} = $logp;
  $Out->{'a'} = $a;
  
  return ($Out);
}



