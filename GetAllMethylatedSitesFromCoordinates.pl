#!/usr/bin/perl
use strict;
#use warnings;
#use autodie ':all';
use Compress::Zlib;

use Getopt::Std ;  
our ($opt_h, $opt_f,$opt_q, $opt_v,$opt_c, $opt_t, $opt_C);
getopts("hvqf:t:c:C:");
&usage if $opt_h;   
sub usage {
  
  die "GetAllMethylatedSitesFromCoordinates.pl -c <coordinate file> -f <file of gzipped filenames >
output is to stdout
 coordinate file structure is 
chr pos
  file of gzipped file structure is
strain filename
-v verbose
-t <int> threshold for coverage. Default is 1
-C <chromosome as chr1, chr2 etc>  process just one chromosome (assumes the gzipped files are in chromosome order)\n";
}

my $threshold = defined ($opt_t) ? $opt_t : 1;
my $Files = {};
my $Coordinates= {};
open (FOFN, $opt_f) || die "Can't open file of file names\n";
my $count = 0;
while (<FOFN>) {
  chomp;
  my @data = split;
  if ($#data == 1) {
    $Files->{$data[1]}{$data[0]}++;
    $count++;
  }
}
if ($count == 0) {
  die "No useful data in $opt_f\n";
}


open (COORDINATES, $opt_c) || die "Can't open file of coordinates \n";
while (<COORDINATES>) {
  chomp;
  my @data = split;
  $Coordinates->{$data[0]}{$data[1]}++;
}
my $file = $opt_f;
my $Data = {};
my $Strains = {};
#chr16   3006198 +       CCT     0       6       1

foreach my $file (sort keys %{$Files}) {
  foreach my $strain (sort keys %{$Files->{$file}}) {
    $Strains->{$strain}++;
    if ($opt_v) {
      warn "Processing $file for $strain\n";
    }
    if (-e $file) {
      open my $gzip, "gzip -cd $file |";
      my $start = 0;
      while (<$gzip>) {
	chomp;
	# for files generated 09/22
	s/_/ /g;
	s/,/ /g;
	
	my ($chr,$pos, $strand, $context, $methyl, $coverage, $junk) = split;
	
	# where coverage and methylation is set to a float
	$methyl =~ s/\.\d+//;
	$coverage =~ s/\.\d+//;
	# print "$methyl\t$coverage\n";
	
	if ($opt_C) {
	  next if ($chr ne $opt_C && $start == 0) ;
	  if ($opt_v) {
	    if ($chr ne $opt_C && $start == 1) {
	      warn "Finished chromosome $opt_c\n";
	    }
	  }
	  last if ($chr ne $opt_C && $start == 1) ;
	  if ($chr eq $opt_C) {
	    $start = 1;
	  }
	}
	my $line = "$strand\t$context\t$methyl\t$coverage";
	# Add on data if there are more than one file for the strain
	if (exists ($Coordinates->{$chr}{$pos})) {
	  if (exists($Data->{$chr}{$pos}{$strain})) {
	    #my @data = split (/\t/ ,  $Data->{$chr}{$pos}{$strain});
	    my ($strand1, $context1, $methyl1, $coverage1) = split (/\t/ ,  $Data->{$chr}{$pos}{$strain});
	    $methyl +=  $methyl1;
	    $coverage += $coverage1;
	   # print "$chr ,$pos $strand1, $context1, $methyl1, $methyl, $coverage1, $coverage\n";
	    my $line = "$strand1\t$context1\t$methyl\t$coverage";
	    $Data->{$chr}{$pos}{$strain} = $line;
	  }
	    else {
	      $Data->{$chr}{$pos}{$strain} = $line;
	    }
	  }  
      }
      close $gzip;
    }
    else {
      if ($opt_v) {
	warn "Can't open file $file\n";
      }
      }
  }
}

foreach my $chr (sort {$a cmp $b} keys %{$Data}) { 
    foreach my $pos (sort {$a <=> $b} keys %{$Data->{$chr}}) {
      print "$chr\t$pos";
      foreach my $strain (sort {$a cmp $b } keys %{$Strains}) {
	if (exists ($Data->{$chr}{$pos}{$strain})) {
	  print "\t$strain\t$Data->{$chr}{$pos}{$strain}";
	}
	else {
	  print "\t$strain\tNA\tNA\tNA\tNA";
	}
      }
      print "\n";			 
    }
  } 
  
  
