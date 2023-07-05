#!/usr/bin/perl
use strict;
use Getopt::Std ;  
our ($opt_h,$opt_i, $opt_g, $opt_G, $opt_e, $opt_c, $opt_b, $opt_m, $opt_a, $opt_R, $opt_K, $opt_x);
getopts("hig:c:Gb:mae:RKx");
&usage if $opt_h;
sub usage {
  die "MatchToGenes -c <coordinates file> -g <gene file>  (refseq format)

e.g. ./MatchToGenes.pl -c coord 
default gene file is:
-g /Users/jf/Documents/MouseAndRatWork/GENE.INFO/refGene.txt

for bed format use -b

-b /Users/jf/MyWorkDocuments/MouseAndRatWork/HiCAndMethylomes/annotations/knownGeneM23.bed 

coordinates file can be a feature file , must start with  chr pos ... 
other information will be output with features appended at the end (not check on expected line length!)
-x order of first columns is pos then chr

-a include header from the coordinate file and add on a column heading for the gene and feature

-m first column contains a marker not a chromosome
-e <50> minimum expected line length in the input file

-R remove lines that hit a gene feature
-K print out only lines that hit a gene feature 
\n";

}


my $Coordinates = {};
my $Locations = {};
my $Genes = {};
my $genefile = defined ($opt_g) ? $opt_g :  "/Users/jf/Documents/MouseAndRatWork/GENE.INFO/refGene.txt";
my  $expectedlength - defined ($opt_e) ? $opt_e : 50;

  open (COORD, $opt_c) || die "Can't find coordinates file $opt_c\n";
  if ($opt_a) {
    my $head  = <COORD>;
    chomp ($head);
    print "$head\tgene\tfeature\n";
  }
  
  while (<COORD>) {
    chomp;
    my @data = split;
    my $chr;
    my $tmp;
    my $pos;
   
      $tmp = $data[0];
    $pos = $data[1];
    if ($opt_x) {
      $tmp = $data[1];
      $pos = $data[0];
    }

    if ($tmp =~ /chr/) {
      $chr = $tmp ;
    }
    else {
      $chr = "chr" . $tmp;
    }
    $Coordinates->{$chr}{$pos} = $_;
    $Locations->{$chr}{$pos}++;
    
   #print "$chr, $pos\n";
  }

#chr1    3905738 3986215 ENSMUST00000194643.1    3       -       3905738 3905738 8553170 3       396,192,69,     0,79421,80408,  uc287gdp.1      none    none    -1,-1,-1,               Gm37381   
if ($opt_b) {
open (GENES, $opt_b) || die "Can't find gene file $opt_b\n";
while (<GENES>) {
  chomp;
  my @data = split (/\t/, $_);
  my $chr = $data[0];
  my $cstart = $data[1];
  my $cend = $data[2];
  my $transcript = $data[3];
  my $strand = $data[5];
  my $gstart = $data[6];
  my $gend = $data[7];
  my $nexons = $data[9];
  my @exonsize = split (",", $data[10]);
  my @exonstart = split (",", $data[11]);

  my $genename = $data[$#data];
  my $name = $transcript . "_" . $genename ;
  
  next if ($#exonsize + 1 != $nexons);
 
  for (my $n = 0; $n <=  $#exonstart; $n++) {
    my  $startf  = $gstart + $exonstart[$n];
    my $endf = $startf + $exonsize[$n];
 
    my $intstart = $endf + 1;
    #my $intend  = $gstart + $exonstart[$n + 1] -1;
     my $intend  = $exonstart[$n+1] -1;
    my $feature = "NA";


    
 if ($n ==0 && $n == $#exonstart ) {  
  # print "$chr\t$name\t$strand\texon\t$startf\t$endf";
	$Genes->{$chr}{$startf}{end} = $endf;
	$Genes->{$chr}{$startf}{name} = $name;
       	$Genes->{$chr}{$startf}{feature} = "exon";
	$Locations->{$chr}{$startf}++; 
 }
   elsif ($n ==0 ) {  
      if ($strand eq "-") {
	$feature = "threeprime";
	my $start = $startf - 2000;
	my $end = $startf -1;
#	print "$chr\t$name\t$strand\t$feature\t$start\t$end\n";
       	$Genes->{$chr}{$start}{end} = $end;
	$Genes->{$chr}{$start}{name} = $name;
       	$Genes->{$chr}{$start}{feature} = $feature;
	$Locations->{$chr}{$start}++;
#	print "$chr\t$name\t$strand\texon\t$startf\t$endf";
	$Genes->{$chr}{$startf}{end} = $endf;
	$Genes->{$chr}{$startf}{name} = $name;
       	$Genes->{$chr}{$startf}{feature} = "exon";
	$Locations->{$chr}{$startf}++;
#	print "\tintron\t$intstart\t$intend\n";
	$Genes->{$chr}{$intstart}{end} = $intend;
	$Genes->{$chr}{$intstart}{name} = $name;
       	$Genes->{$chr}{$intstart}{feature} = "intron";
	$Locations->{$chr}{$intstart}++;
      }
      else {
	$feature = "fiveprime";
	my $start = $startf - 2000;
	my $end = $startf -1;
#	print "$chr\t$name\t$strand\t$feature\t$start\t$end\n";
 	$Genes->{$chr}{$start}{end} = $end;
	$Genes->{$chr}{$start}{name} = $name;
       	$Genes->{$chr}{$start}{feature} = $feature;
	$Locations->{$chr}{$start}++;
#	print "$chr\t$name\t$strand\texon\t$startf\t$endf";
	$Genes->{$chr}{$startf}{end} = $endf;
	$Genes->{$chr}{$startf}{name} = $name;
       	$Genes->{$chr}{$startf}{feature} = "exon";
	$Locations->{$chr}{$startf}++;
#	print "\tintron\t$intstart\t$intend\n";
	$Locations->{$chr}{$startf}++;
	$Genes->{$chr}{$intstart}{end} = $intend;
	$Genes->{$chr}{$intstart}{name} = $name;
       	$Genes->{$chr}{$intstart}{feature} = "intron";
	$Locations->{$chr}{$intstart}++;
	
      }
    }
    elsif ($n == $#exonstart ) {
      if ($strand eq "-") {
	$feature = "fiveprime";
#	print "$chr\t$name\t$strand\texon\t$startf\t$endf\n";
	$Genes->{$chr}{$startf}{end} = $endf;
	$Genes->{$chr}{$startf}{name} = $name;
       	$Genes->{$chr}{$startf}{feature} = "exon";
	$Locations->{$chr}{$startf}++;
	
	my $start = $endf + 1;
	my $end = $endf + 2000 ;
#	print "$chr\t$name\t$strand\t$feature\t$start\t$end\n";
	$Genes->{$chr}{$start}{end} = $end;
	$Genes->{$chr}{$start}{name} = $name;
       	$Genes->{$chr}{$start}{feature} = $feature;
	$Locations->{$chr}{$start}++;
      }
      else {
	$feature = "threeprime";	
#	print "$chr\t$name\t$strand\texon\t$startf\t$endf\n";
	$Genes->{$chr}{$startf}{end} = $endf;
	$Genes->{$chr}{$startf}{name} = $name;
       	$Genes->{$chr}{$startf}{feature} = "exon";
	$Locations->{$chr}{$startf}++;
	
	my $end = $endf + 2000;
	my $start = $endf + 1;
#	print "$chr\t$name\t$strand\t$feature\t$start\t$end\n";
      	$Genes->{$chr}{$start}{end} = $end;
	$Genes->{$chr}{$start}{name} = $name;
       	$Genes->{$chr}{$start}{feature} = $feature;
	$Locations->{$chr}{$start}++;
      }
    }
    else {
  #    print "$chr\t$name\t$strand\texon\t$startf\t$endf";
      $Genes->{$chr}{$startf}{end} = $endf;
      $Genes->{$chr}{$startf}{name} = $name;
      $Genes->{$chr}{$startf}{feature} = "exon";
      $Locations->{$chr}{$startf}++;
      
 #     print "\tintron\t$intstart\t$intend\n";
      $Genes->{$chr}{$intstart}{end} = $intend;
      $Genes->{$chr}{$intstart}{name} = $name;
      $Genes->{$chr}{$intstart}{feature} = "intron";
      $Locations->{$chr}{$intstart}++;
    }
  }
#  print "\n";
}
}

else  {
#warn "opening gene file $genefile\n";
#156	NM_054041	chr6	-	87133852	87335775	87137094	87335442	18	87133852,87180109,87188106,87188735,87204408,87217273,87240949,87242486,87255782,87269242,87282686,87284253,87287000,87288151,87288757,87308327,87312322,87335296,	87137355,87180190,87188274,87188831,87204450,87217369,87241028,87242556,87255881,87269303,87282767,87284322,87287080,87288185,87288839,87308399,87312394,87335775,	0	Antxr1	cmpl	cmpl	0,0,0,0,0,0,2,1,1,0,0,0,1,0,2,2,2,0,


 open (GENES, $genefile) || die "Can't find gene file $genefile\n";
  while (<GENES>) {
    chomp;
    my @data = split (/\t/, $_);
    my $chr = $data[2];
    my $gstart = $data[4];
    my $strand = $data[3];
    
    my @exonstart = split (",", $data[9]);
    my @exonend = split (",", $data[10]);
    my $name = $data[12];

    
    
    #print "$chr, $gstart, $strand, $name,  @exonstart\n";
    
    for (my $n = 0; $n <=  $#exonstart; $n++) {
      my $startf = $exonstart[$n];
      my $endf = $exonend[$n];
      my $intstart = $endf + 1;
      #my $intend  = $gstart + $exonstart[$n + 1] -1;
      my $intend  = $exonstart[$n+1] -1;
      my $feature = "NA";

      
      if ($n ==0 && $n == $#exonstart ) {
	
	# print "$chr\t$name\t$strand\texon\t$startf\t$endf";
	$Genes->{$chr}{$startf}{end} = $endf;
	$Genes->{$chr}{$startf}{name} = $name;
       	$Genes->{$chr}{$startf}{feature} = "exon";
	$Locations->{$chr}{$startf}++; 
      }
      elsif ($n ==0 ) {  
	if ($strand eq "-") {
	  $feature = "threeprime";
	  my $start = $startf - 2000;
	  my $end = $startf -1;
	  #       print "$chr\t$name\t$strand\t$feature\t$start\t$end\n";
	  $Genes->{$chr}{$start}{end} = $end;
	  $Genes->{$chr}{$start}{name} = $name;
	  $Genes->{$chr}{$start}{feature} = $feature;
	  $Locations->{$chr}{$start}++;
	  #	print "$chr\t$name\t$strand\texon\t$startf\t$endf";
	  $Genes->{$chr}{$startf}{end} = $endf;
	  $Genes->{$chr}{$startf}{name} = $name;
	  $Genes->{$chr}{$startf}{feature} = "exon";
	  $Locations->{$chr}{$startf}++;
	  #	print "\tintron\t$intstart\t$intend\n";
	  $Genes->{$chr}{$intstart}{end} = $intend;
	  $Genes->{$chr}{$intstart}{name} = $name;
	  $Genes->{$chr}{$intstart}{feature} = "intron";
	  $Locations->{$chr}{$intstart}++;
	}
	else {
	  $feature = "fiveprime";
	  my $start = $startf - 2000;
	  my $end = $startf -1;
	  #	print "$chr\t$name\t$strand\t$feature\t$start\t$end\n";
	  $Genes->{$chr}{$start}{end} = $end;
	  $Genes->{$chr}{$start}{name} = $name;
	  $Genes->{$chr}{$start}{feature} = $feature;
	  $Locations->{$chr}{$start}++;
	  #	print "$chr\t$name\t$strand\texon\t$startf\t$endf";
	  $Genes->{$chr}{$startf}{end} = $endf;
	  $Genes->{$chr}{$startf}{name} = $name;
	  $Genes->{$chr}{$startf}{feature} = "exon";
	  $Locations->{$chr}{$startf}++;
	  #	print "\tintron\t$intstart\t$intend\n";
	  $Locations->{$chr}{$startf}++;
	  $Genes->{$chr}{$intstart}{end} = $intend;
	  $Genes->{$chr}{$intstart}{name} = $name;
	  $Genes->{$chr}{$intstart}{feature} = "intron";
	  $Locations->{$chr}{$intstart}++;
	  
	}
      }
      elsif ($n == $#exonstart ) {
	if ($strand eq "-") {
	  $feature = "fiveprime";
	  #	print "$chr\t$name\t$strand\texon\t$startf\t$endf\n";
	  $Genes->{$chr}{$startf}{end} = $endf;
	  $Genes->{$chr}{$startf}{name} = $name;
	  $Genes->{$chr}{$startf}{feature} = "exon";
	  $Locations->{$chr}{$startf}++;
	  
	  my $start = $endf + 1;
	  my $end = $endf + 2000 ;
	  #	print "$chr\t$name\t$strand\t$feature\t$start\t$end\n";
	  $Genes->{$chr}{$start}{end} = $end;
	  $Genes->{$chr}{$start}{name} = $name;
	  $Genes->{$chr}{$start}{feature} = $feature;
	  $Locations->{$chr}{$start}++;
	}
	else {
	  $feature = "threeprime";	
	  #	print "$chr\t$name\t$strand\texon\t$startf\t$endf\n";
	  $Genes->{$chr}{$startf}{end} = $endf;
	  $Genes->{$chr}{$startf}{name} = $name;
	  $Genes->{$chr}{$startf}{feature} = "exon";
	  $Locations->{$chr}{$startf}++;
	  
	  my $end = $endf + 2000;
	  my $start = $endf + 1;
	  #	print "$chr\t$name\t$strand\t$feature\t$start\t$end\n";
	  $Genes->{$chr}{$start}{end} = $end;
	  $Genes->{$chr}{$start}{name} = $name;
	  $Genes->{$chr}{$start}{feature} = $feature;
	  $Locations->{$chr}{$start}++;
	}
      }
      else {
	#    print "$chr\t$name\t$strand\texon\t$startf\t$endf";
	$Genes->{$chr}{$startf}{end} = $endf;
	$Genes->{$chr}{$startf}{name} = $name;
	$Genes->{$chr}{$startf}{feature} = "exon";
	$Locations->{$chr}{$startf}++;
	
	#     print "\tintron\t$intstart\t$intend\n";
	$Genes->{$chr}{$intstart}{end} = $intend;
	$Genes->{$chr}{$intstart}{name} = $name;
	$Genes->{$chr}{$intstart}{feature} = "intron";
	$Locations->{$chr}{$intstart}++;
      }
    }
  }
}


  
foreach my $chr (sort {$a cmp $b} keys %{$Locations}) {
    my $name = "NA";
    my $end = 0;
    my $start = 0;
    my $feature = "NA";
    foreach my $pos (sort {$a <=> $b } keys %{$Locations->{$chr}}) {
  #    if ($chr eq "chr16") {
#	print "$chr, $pos, $Genes->{$chr}{$pos}{name},  $Coordinates->{$chr}{$pos} \n";
 #     }
      if (exists ($Genes->{$chr}{$pos}))  {	
	$start = $pos;
	$end = $Genes->{$chr}{$pos}{end};
	$name = $Genes->{$chr}{$pos}{name};
	$feature = $Genes->{$chr}{$pos}{feature};
      }
      if (exists ( $Coordinates->{$chr}{$pos})) {     
	#next if ($end ==0);
#	next if ($start ==0);
        if ($pos >= $start && $pos <= $end) {
	  unless ($opt_R) {
	    if (length ($Coordinates->{$chr}{$pos}) > $expectedlength) {
	      print "$Coordinates->{$chr}{$pos}\t$name\t$feature\n";
	    }
	  }
	}
	else{
	  unless ($opt_K) {
	  if (length ($Coordinates->{$chr}{$pos}) > $expectedlength) {
	    print	"$Coordinates->{$chr}{$pos}\tNA\tNA\n";
	  }
	}
	}
	
      }
    }
  }

