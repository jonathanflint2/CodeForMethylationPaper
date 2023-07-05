#!/usr/bin/perl
use strict;
use warnings;
#use autodie ':all';

# use Compress::Zlib <- don't use this - too slow
use Time::HiRes 'time';
use Getopt::Std ;  
our ($opt_h, $opt_f,$opt_q, $opt_v, $opt_t);
getopts("hvqf:t:");
&usage if $opt_h;   
sub usage {
  die "FindMethylatedCoordinatesForAllStrainsInOneTissue.pl -f <fofn>

outputs coordinates of methylated sites
file of file names has two columns , strain and file (strain is not used, included for forward compatibility)
-t <int> threshold for coverage. Default is greater than or equal to 1\n";
}

my $threshold = defined ($opt_t) ? $opt_t : 1;
#chr16   3006198 +       CCT     0       6       1
my $Coordinates = {};
my $Files = {};
my $count = 0;
open (FOFN, $opt_f) || die "Can't open file of file names\n";
while (<FOFN>) {
chomp;
my @data = split;
if ($#data == 1) {
  $Files->{$data[1]} = $data[0];
  $count++;
}
}
if ($count ==0) {
  die "No useful files found in $opt_f\n";
}
foreach my $file (sort keys %{$Files}) {
if ($opt_v) {
warn "Processing $file\n";
}
if (-e $file) {
open my $gzip, "gzip -cd $file |";
while (<$gzip>) {

  #for files produced in 09/22
  s/_/ /g;
  s/,/ /g;
  
  my ($chr,$pos, $strand, $context, $methyl, $coverage, $junk) = split;
  if ($methyl > 0 && $coverage >= $threshold) {
#    print "$chr\t$pos\t$strand\t$context\t$methyl\t$coverage\n";
	$Coordinates->{$chr}{$pos}++;
  }
  
}
close ($gzip);
}
else {
  if ($opt_v) {
    warn "Can't find $file\n";
  }
}
}

foreach my $chr (sort {$a cmp $b} keys %{$Coordinates}) {
foreach my $pos (sort {$a <=> $b } keys %{$Coordinates->{$chr}}) {
print "$chr\t$pos\n";
}
}

