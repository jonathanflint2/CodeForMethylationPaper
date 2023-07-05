#!/usr/bin/perl
use strict;

use Getopt::Std ;  
our ($opt_h, $opt_f,$opt_q, $opt_v,$opt_c, $opt_t, $opt_w, $opt_H, $opt_m, $opt_a, $opt_c);
getopts("hf:v:qcw:Hmac");
&usage if $opt_h;   
sub usage {
  die "AddVariantsToDiffFile8Strains.pl -f diff file  
-v <variants file>  (must have a header with strain names)
-w window to search either side of the methylation site (default is 3)
default variants file is /Users/jf/MyWorkDocuments/MouseAndRatWork/MOUSE.SEQUENCE.VARIANTS/EightStrains.vars.txt
-H add header (hard coded)
-a add two column headings, var.pos and var.sdp, to the header
-m first column is a marker and last column is a diff marker
-c correct orientation - assumes there is a column headed 'strand' for this to work
\n";
}
#my $header = "chr\tpos\ts1\ts1.strand\ts1.seq\ts1.me\ts1.cov\ts2\ts2.strand\ts2.seq\ts2.me\ts2.cov\ts3\ts3.strand\ts3.seq\ts3.me\ts3.cov\ts4\ts4.strand\ts4.seq\ts4.me\ts4.cov\ts5\ts5.strand\ts5.seq\ts5.me\ts5.cov\ts6\ts6.strand\ts6.seq\ts6.me\ts6.cov\ts7\ts7.strand\ts7.seq\ts7.me\ts7.cov\ts8\ts8.strand\ts8.seq\ts8.me\ts8.cov\ttotal.cov\tvar.pos\ts1.s\ts2.s\ts3.s\ts4.s\ts5.s\ts6.s\ts7.s\ts8.s";

my $header = "chr\tpos\ts1\ts1.strand\ts1.seq\ts1.me\ts1.cov\ts2\ts2.strand\ts2.seq\ts2.me\ts2.cov\ts3\ts3.strand\ts3.seq\ts3.me\ts3.cov\ts4\ts4.strand\ts4.seq\ts4.me\ts4.cov\ts5\ts5.strand\ts5.seq\ts5.me\ts5.cov\ts6\ts6.strand\ts6.seq\ts6.me\ts6.cov\ts7\ts7.strand\ts7.seq\ts7.me\ts7.cov\ts8\ts8.strand\ts8.seq\ts8.me\ts8.cov\ttotal.cov\tlogp\tvar.pos\tseq.sdp";

if ($opt_m) {
   $header = "marker\tchr\tpos\ts1\ts1.strand\ts1.seq\ts1.me\ts1.cov\ts2\ts2.strand\ts2.seq\ts2.me\ts2.cov\ts3\ts3.strand\ts3.seq\ts3.me\ts3.cov\ts4\ts4.strand\ts4.seq\ts4.me\ts4.cov\ts5\ts5.strand\ts5.seq\ts5.me\ts5.cov\ts6\ts6.strand\ts6.seq\ts6.me\ts6.cov\ts7\ts7.strand\ts7.seq\ts7.me\ts7.cov\ts8\ts8.strand\ts8.seq\ts8.me\ts8.cov\ttotal.cov\tdiffmarker\tlogp\tvar.pos\tseq.sdp";
}
my $Variants = {};
my $varfile = defined ($opt_v) ? $opt_v : "/Users/jf/MyWorkDocuments/MouseAndRatWork/MOUSE.SEQUENCE.VARIANTS/EightStrains.vars.txt";
my $window = defined ($opt_w)? $opt_w : 3;
my $LookUp = {};

##############
# note that blk6, not b6, is used to be compatible with Heffel's me files
# order is
# aj balb b6 cast dba fvb pwk wsb
#############

$LookUp->{'blk6'} = "ref";
$LookUp->{'aj'} = "A_J";
$LookUp->{'balb'} = "BALB_cJ";
$LookUp->{'cast'} = "CAST_EiJ";
$LookUp->{'dba'} = "DBA_2J";
$LookUp->{'fvb'} = "FVB_NJ";
$LookUp->{'pwk'} = "PWK_PhJ";
$LookUp->{'wsb'} = "WSB_EiJ";

my $Seq = {};
$Seq->{"A"} = "T";
$Seq->{"C"} = "G";
$Seq->{"G"} = "C";
$Seq->{"T"} = "A";
$Seq->{"H"} = "H";
$Seq->{"N"} = "N";

open (VAR, $varfile) || die "Can't open variants file $varfile\n";
my $VHeader = {};
my $line = <VAR>;
chomp ($line);
my @data = split (/\t/, $line);
for (my $n = 0; $n <=$#data; $n++) {
  #print "$data[$n] $n\n";
  $VHeader->{$data[$n]} = $n;
}

# order of variants in the  seqsdp  from: foreach my $strain (keys %{$LookUp}) {
# gives
# aj balb b6 cast dba fvb pwk wsb

while (<VAR>) {
  chomp;
  my @data = split(/\t/, $_);
  #chrom	pos	id	ref	alt	A_J	BALB_cJ	CAST_EiJ	DBA_2J	FVB_NJ	PWK_PhJ	WSB_EiJ
  foreach my $strain (sort keys %{$LookUp}) {
  
    
    my $varstrain = $LookUp->{$strain};
    
    my $col = $VHeader->{$varstrain};
    my $chr = "chr" . $data[0];
my @seq = split (/\s+/, $data[$col]);
   # print "$strain\t$varstrain\t$col\t$data[$col]\t$seq[0]\n";
    my $var = $seq[0];
    if ($#seq == 1) {
    if ($seq[0] ne $seq[1]) {
      $var = "H";
    }
  }
    if ($var eq "NA") {
      $var = "N";
    }
    $Variants->{$chr}{$data[1]}{$strain} = $var;
  }
}
#
#chr1	3003898	aj	NA	NA	NA	NA	balb	+	CGG	3	3	blk6	NA	NA	NA	NA	cast	NA	NA	NA	NA	dba	NA	NA	NA	NA	fvb	+	CGG	pwk	+	CGG	1	1	wsb	NA	NA	NA	NA
if ($opt_H) {
  print "$header\n";
}
open (DIFF, $opt_f) || die "Can't open diff file $opt_f\n";
my $head;
if ($opt_a) {
  $head = <DIFF>;
  chomp($head);
  print "$head\tvar.pos\tvar.sdp\n";
}

my $DHeader = {};
if ($opt_c) {
  my @data = split (/\t/, $head);
  for (my $n = 0; $n <= $#data; $n++) {
    $DHeader->{$data[$n]} = $n;
  }
  unless (exists ($DHeader->{'strand'})) {
    die "Can't find a column labelled strand in the header\n";
  }
}
while (<DIFF>) {
  chomp;
  my @data = split;
 
  my $chr = $data[0];
  my $start = $data[1] - $window;
  my $end =  $data[1] + $window;

if ($opt_m) {
  $chr = $data[1];
  $start = $data[2] - $window;
  $end =  $data[2] + $window;
}
  my $found = 0;
  my $loc;
  my $sdp;
  
  for (my $n = $start; $n <= $end; $n++) {
    my $location;
    my $seqsdp;
    if (exists ($Variants->{$chr}{$n})) {   
      $location = $n - $data[1];
      if ($opt_m) {
	$location = $n - $data[2];
      }    
     # print "$location\n";
      $found++;   
   
      foreach my $strain (sort {$a cmp $b} keys %{$LookUp}) {
	   #############
      # see above for ordering issue
      ###########
      #aj
      #balb
      #blk6
      ##cast
      #dba
      #fvb
      #pwk
      #wsb
      ##############
      # while meth is aj b6 balb   
      #############
	$seqsdp .= $Variants->{$chr}{$n}{$strain};
#	print "$Variants->{$chr}{$n}{$strain}";
      }

      if ($opt_c) {
	my $strand = $data[$DHeader->{'strand'}];
	if ($strand eq "-") {
	  $location = $location * -1;
	}
	my $newsdp;
	my @sdp = split (//, $seqsdp);
	foreach my $var (@sdp) {
	  $newsdp .= $Seq->{$var};
	}
	$seqsdp = $newsdp;
      }
       print "$_\t$location\t$seqsdp\n";
     # $loc = $location;
     # $sdp = $seqsdp;
    }
  }
  if ($found == 0) {
    print "$_\tNA\tNA\n";
  }
  
  #if ($found > 0) {
  #  print "$_\t$loc\t$sdp\n";
  #}
  #else {
  #  print "$_\tNA\tNA\n";
  #}
}

