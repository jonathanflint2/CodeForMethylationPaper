#!/usr/bin/perl
use strict;
use Getopt::Std;

our($opt_h, $opt_t,$opt_f, $opt_z, $opt_b, $opt_F, $opt_S, $opt_v, $opt_c, $opt_T, $opt_C, $opt_V, $opt_H, $opt_M, $opt_m, $opt_S);

getopts("hf:bT:c:CV:vHS:z:");
&usage if $opt_h;
sub usage {
  die "FindCpGs.pl  -c chromosome  (assumes sequence file to process is of the form chr1.fa in a given directory)
  -f <annotation file> -z compressed
 
this version includes a column with the gene name so it can be matched by other programs 
 also provides a switch to output ONLY dba and b6 (for which we have expression data)
-b only dba and b6 data
-S sequence dir
-V variant dir
-H omit het sites from the variant file
-v verbose

\n";
}

my $seqdir = defined ($opt_S) ? $opt_S :  "MOUSE.GENOME/";
my $variantdir = defined ($opt_V)?  $opt_V : "MOUSE.SEQUENCE.VARIANTS/";

my $chr = defined ($opt_c) ? $opt_c : "chr1";
my $Segments = {};

my $sequencefile = $seqdir . $chr . ".fa";

unless (-e $sequencefile) {
  die "can't find the fasta file $sequencefile\n";
}
my $variantsfile =  $variantdir . $chr . ".EightStrains.vars.txt";

unless (-e $variantsfile) {
  die "can't find the variant file $variantsfile\n";
}


my $LookUp = {};
#
#for compatibility with the names in Heffel's files - and see AddVariantsToDiff8Strains.pl -
# this alters the sort order
# 
$LookUp->{'b6'} = "ref";
$LookUp->{'aj'} = "A_J";
$LookUp->{'balb'} = "BALB_cJ";
$LookUp->{'cast'} = "CAST_EiJ";
$LookUp->{'d2'} = "DBA_2J";
$LookUp->{'fvb'} = "FVB_NJ";
$LookUp->{'pwk'} = "PWK_PhJ";
$LookUp->{'wsb'} = "WSB_EiJ";

my $Variants = {};
open (VAR, $variantsfile ) || die "can't open the variant file $variantsfile\n";
my $VHeader = {};
my $line = <VAR>;
chomp ($line);
my @data = split (/\t/, $line);
for (my $n = 0; $n <=$#data; $n++) {
  #print "$data[$n] $n\n";
  $VHeader->{$data[$n]} = $n;
}


while (<VAR>) {
  chomp;
  my @data = split(/\t/, $_);
  #chrom	pos	id	ref	alt	A_J	BALB_cJ	CAST_EiJ	DBA_2J	FVB_NJ	PWK_PhJ	WSB_EiJ
  my $seqsdp;
  my $bad = 0;
  foreach my $strain (sort keys %{$LookUp}) {
    
    # order of variants in the  seqsdp  from: foreach my $strain (sort keys %{$LookUp}) {
    # gives
    # aj balb b6 cast dba fvb pwk wsb
    
    my $varstrain = $LookUp->{$strain};
    my $col = $VHeader->{$varstrain};
    my $chr = "chr" . $data[0];
    my @seq = split (/\s+/, $data[$col]);
    
    # this generates the order 
    # print "$strain, $varstrain, $seq[0]\n";
    # print "$data[0]\t$data[1]\t$strain\t$varstrain\t$col\t$data[$col]\t$seq[0]\n";
    # Keep this for now 
    
    my $var = $seq[0];
    if ($#seq == 1) {
      if ($seq[0] ne $seq[1]) {
	$var = "H";
	$bad++;	  
      }
    }
    if ($var eq "NA") {
      $var = "N"; 
    }
   
    $Variants->{$chr}{$data[1]}{$strain} = $var;
  }
}

my $AllowedTypes = {};
$AllowedTypes->{'C'} = 'G';
$AllowedTypes->{'CG'} = 'G';
$AllowedTypes->{'A'} = 'T';
$AllowedTypes->{'AT'} = 'T';

my $typetoget = defined ($opt_T) ? $opt_T : "CG";

unless (exists($AllowedTypes->{$typetoget})) {
  die "Allowed types to get are \"C\", \"CG\" , \"AT\" or \"A\". Default is \"CG\"\n";
}


my $Queries = {};
my $Types = {};

$Types->{'CH'} = "CH";
$Types->{'CG'} = "CG";
$Types->{'CHG'} = "CH";
$Types->{'MIXED'} = "CG";
$Types->{'BOTH'} = "CG";
$Types->{'ALL'} = "ALL";
my $Pairs = {};

my $Strains = {};
if ($opt_b) {
$Strains->{'b6'}++;
$Strains->{'d2'}++;
}
else {
$Strains->{'aj'}++;
$Strains->{'b6'}++;
$Strains->{'balb'}++;
$Strains->{'cast'}++;
$Strains->{'d2'}++;
$Strains->{'fvb'}++;
$Strains->{'pwk'}++;
$Strains->{'wsb'}++;
}
my $Repeats = {};
$Repeats->{"SINE"} = "SINE";
$Repeats->{"LTR"} = "LTR";
$Repeats->{"LTR/ERVK"} = "LTR"; 
$Repeats->{"LTR/ERVK?"} = "LTR";
$Repeats->{"LTR/ERVL"} = "LTR";
$Repeats->{"LTR/ERVL-MaLR"} = "LTR"; 
$Repeats->{"LTR/ERVL?"} = "LTR";
$Repeats->{"LTR/Gypsy"} = "LTR";
$Repeats->{"LTR?"} = "LTR";
$Repeats->{"SINE/Deu"} = "SINE";
$Repeats->{"SINE/ID"} = "SINE"; 
$Repeats->{"SINE?"} = "SINE"; 
$Repeats->{"LTR/ERV1"} = "LTR";
$Repeats->{"LINE/RTE-BovB"} = "LINE";
$Repeats->{"LINE/CR1" } = "LINE";
$Repeats->{"LTR/Gypsy?"} = "LTR";
$Repeats->{"SINE/tRNA"} = "SINE"; 
$Repeats->{"LTR/ERV1?"} = "LTR";
$Repeats->{"LINE/L2"} = "LINE_L2";
$Repeats->{"LINE/RTE-X"} = "LINE";
$Repeats->{"SINE/Alu"} = "SINE"; 
$Repeats->{"LINE/L1"} = "LINE_L1"; 
$Repeats->{"SINE/MIR"} = "SINE"; 
$Repeats->{"LINE/L1?"} = "LINE"; 
$Repeats->{"SINE/B2"} = "SINE_B2"; 
$Repeats->{"LINE/Dong-R4"} = "LINE";
$Repeats->{"SINE/B4"} = "SINE_B4"; 
$Repeats->{"Simple_repeat"}= "Simple_repeat";
$Repeats->{"NA"} = "NA";

my $RegContext = {};
$RegContext->{'En-Pd'} = "Enhancer";
$RegContext->{'En-Pp'} = "Enhancer";
$RegContext->{'En-Sd'} = "Enhancer";
$RegContext->{'En-Sp'} = "Enhancer";
$RegContext->{'En-W'} = "Enhancer";
$RegContext->{'Pr-A'} = "Promoter";
$RegContext->{'Pr-B'} = "Promoter";
$RegContext->{'Pr-F'} = "Promoter";
$RegContext->{'Pr-W'} = "Promoter";
$RegContext->{'Hc-H'} = "Het";	
$RegContext->{'Hc-P'} = "Het";
$RegContext->{'Tr-I'} = "Transcription";	
$RegContext->{'Tr-P'} = "Transcription";
$RegContext->{'Tr-S'} = "Transcription";
$RegContext->{'NS'} = "NA";
$RegContext->{'enhD'} = "enhD";
$RegContext->{'enhP'} = "enhP";
$RegContext->{'prom'} = "prom";
$RegContext->{'K4m3'} = "K4m3";
$RegContext->{'CTCF'} = "CTCF";
my $Genes = {};
$Genes->{'exon'} = "exon";
$Genes->{'intron'} = "intron";
$Genes->{'fiveprime'} = "fiveprime";
$Genes->{'threeprime'} = "threeprime";
$Genes->{'NA'} = "NA";


###############
# Open annotation file
###############
my $Methylation = {};
my $Coverage = {};
my $MHeader = {};
my $Annotations = {};
#my $Heritability = {};
my $CpGIslands = {};
my $TSS = {};
my $GeneNames = {};

my $file_handle;
if ($opt_z) {
  open  $file_handle, "gzip -cd $opt_z |";
}
elsif ($opt_f)  {
  open ( $file_handle, $opt_f) || die "Can't open file $opt_f\n";
}
else {
  &usage;
}


my $head = <$file_handle>;
my $Header = {};
chomp ($head);
my @hdata = split (/\t/, $head);
for (my $n = 0; $n <= $#hdata; $n++) {
  $Header->{$hdata[$n]} = $n;  
}

unless (exists ( $Header->{'b6.cov'})) {
  die "Incorrectly formatted header\n";
}
#chr     pos     aj.mc   aj.cov  b6.mc   b6.cov  balb.mc balb.cov        cast.mc cast.cov        d2.mc   d2.cov  fvb.mc  fvb.cov pwk.mc  pwk.cov wsb.mc  wsb.cov        strand  type    total.mc        total.cov       fraction.mc     n.strains       keep    logp    z.diff  y.diff  pair    var.pos var.sdp gene    feature        ccre    chromhmm        repeat  repeat.pos      imprint y       n.logp  h2

##############
# Process the methylation file
##############
while (<$file_handle>) {
  chomp;
  my @tabdata = split(/\t/, $_);
  
  my $chr = $tabdata[$Header->{'chr'}];
  if (defined ($opt_c)) {
    unless ($chr eq $opt_c) {
      next;
    }
  }
  my $pos =    $tabdata[$Header->{'pos'}]; 
  my $pair =    $tabdata[$Header->{'pair'}];
  my $atype =  $tabdata[$Header->{'type'}];
#  my $h2 = "NA";

#if (exists ($tabdata[$Header->{'h2'}])) {
#    $h2 = $tabdata[$Header->{'h2'}];
# }

my $tss = "NA";
if (exists ($tabdata[$Header->{'tss'}])) {
    $tss = $tabdata[$Header->{'tss'}];
 }

#if ($tss ne "NA") {
#print "$chr\t$pos\t$tss\n";
#}

my $cgi = "NA";

if (exists ($tabdata[$Header->{'cpgisland'}])) {
    $cgi = $tabdata[$Header->{'cpgisland'}];
 }
elsif (exists ($tabdata[$Header->{'cgi'}])) {
    $cgi = $tabdata[$Header->{'cgi'}];
 }

my $genename = "NA";
if (exists ($tabdata[$Header->{'gene'}])) {
    $genename = $tabdata[$Header->{'gene'}];
 }
  
my $type = $Types->{$atype};
next unless ($type eq $typetoget);
  
  my $strand =  $tabdata[$Header->{'strand'}];
  ###########################
  #only take the + strand [ this may introduce a bias but it sure makes the coding easier!
  ###########################
  next if ($strand eq "-");
  
  ###########################
  #output coverages so that the pairwise comparison can define which strains are methylated by
  # downsampling coverage to make each pair equal ( see CompareMultiStrainMethylationFromCountedPairwiseComparison.pl ) 
  ###########################
  my $StrainMeth = {};
  my $StrainCov = {};
 # print "$chr\t$pos";

  foreach my $strain (keys %{$Strains} ) {
    my $cname = $strain . ".cov";
    my $mname = $strain . ".mc";
    my $cov = $tabdata[$Header->{$cname}];
    my $meth = $tabdata[$Header->{$mname}];
    # print "\t$strain\t$meth\t$cov";
    $StrainMeth->{$strain} = $meth;
    $StrainCov->{$strain} = $cov;
   
  #  print "\t here-> $StrainMeth->{$strain}\n";
  }
  
 
  #################
  # define annotations
  ################
  my $Tmp = {};
  my $foundannot = 0;
  if ( exists($Genes->{$tabdata[$Header->{'feature'}]})) {
    #&& $tabdata[$Header->{'feature'}] != "NA" ) {
    unless ($Genes->{$tabdata[$Header->{'feature'}]} eq "NA") {
      $Tmp->{$Genes->{$tabdata[$Header->{'feature'}]}}++;
      $foundannot++;
    }
  }
  if (exists($RegContext->{ $tabdata[$Header->{'chromhmm'}]})) {
    #&& $tabdata[$Header->{'chromhmm'}] != "NA" ) {
    unless ($RegContext->{ $tabdata[$Header->{'chromhmm'}]} eq "NA") {
      $Tmp->{$RegContext->{ $tabdata[$Header->{'chromhmm'}]}}++;
       $foundannot++;
    }
  }
  if (exists($RegContext->{$tabdata[$Header->{'ccre'}]}) ){
    #&& $tabdata[$Header->{'ccre'}] != "NA") {
    unless ($RegContext->{$tabdata[$Header->{'ccre'}]} eq "NA") {
      $Tmp->{$RegContext->{$tabdata[$Header->{'ccre'}]}}++;
       $foundannot++;
    }
  }
  if (exists($Repeats->{$tabdata[$Header->{'repeat'}]})) {
    #&&  $tabdata[$Header->{'repeat'}] != "NA") {
    unless ($Repeats->{$tabdata[$Header->{'repeat'}]} eq "NA" ) {
      $Tmp->{$Repeats->{$tabdata[$Header->{'repeat'}]}}++;
       $foundannot++;
    }
  }
  
  ### keep the results
  
  $Methylation->{$chr}{$pos} = $StrainMeth;
  $Coverage->{$chr}{$pos} = $StrainCov;

#  $Heritability->{$chr}{$pos} = $h2;

  $CpGIslands->{$chr}{$pos} = $cgi;
  $TSS->{$chr}{$pos} = $tss;	

#if ($TSS->{$chr}{$pos} ne "NA")  {
#print "$TSS->{$chr}{$pos} \n";
#}
  $GeneNames->{$chr}{$pos} = $genename;;

if ( $foundannot > 0) {
  #  foreach my $f (sort keys %{$Tmp}) {
  #    print "$f\n";
  #  }
    $Annotations->{$chr}{$pos} = $Tmp;
  }
  
  if ($opt_v) {
    warn "Finished reading annotations file\n";
  }
}



  
  
################
# extract the sequence
################
my @output;


open (SEQ,  $sequencefile ) || die "Cannot open sequence file $sequencefile";
my  $head = <SEQ>;
push (@output, $head);
my $line;
while (<SEQ>) {
  chomp;
  #my @data = split (/\n/, $_);
  $line.= $_
}
push (@output, $line);

#print "@output\n";


my $chr;
my $start = 0;
my $end = 2000000000;

#####################
# Print header
####################
print "chr\tpos\tcontext";
foreach my $strain (sort keys  %{$Strains}) {
  my $out = $strain . ".var";
  print "\t$out";
}

foreach my $strain (sort keys  %{$Strains}) {
  my $out = $strain . ".me";
  print "\t$out";
  my $out = $strain . ".cov";
   print "\t$out";
}

#print "\tannotations\th2\tcgi\ttss\tgene\n";

print "\tannotations\tcgi\ttss\tgene\n";

for (my $n =0; $n <= $#output; $n++) {
  chomp ($output[$n]);  
  if ($output[$n] =~ /chr/) {
    $_ = $output[$n];
    s/>//;
    s/:/ /;
    s/-/ /;
    my @info = split;
    $chr = $info[0];
    if ($#info > 0)  {
    $start = $info[1];
    $end = $info[2];
  }
    #print "$chr $start $end\n";
  }
  else {
    
    $output[$n] =~ tr/a-z/A-Z/;
    my @sequences = split (//,  $output[$n]);
    #  my $numberC =   scalar @matches;
    
        
    #start sequences 
    for (my $n =0; $n <= $#sequences; $n++) {
      my $o = $n+1;
      my $position = $start + $n ;
      my $testp = $position + 1;
      my $test2p = $position + 2;
      my $context = $sequences[$n-2] . $sequences[$n-1].  $sequences[$n] . $sequences[$n+1] .  $sequences[$n+2] ;
      my $CpG = $sequences[$n] . $sequences[$o] ;
      my $C =  $sequences[$n];
      my $found = 0;
      my $ToGetStrains = {};
      
      if ($typetoget eq $CpG) {
	
	$found++;
	# set the strains to the B6 reference
	foreach my $strain (sort keys %{$Strains}) {
	  $ToGetStrains->{$strain} =  $sequences[$n] .  $sequences[$n+1];
	}
      }
      ########################
      # check for variants and see if there is CpG in one of the other strains
      ########################
 
      if (exists ( $Variants->{$chr}{$testp} )) {
	foreach my $strain (sort keys %{$Strains}) {
	  if (exists ($Variants->{$chr}{$testp}{$strain}) )  {
	    if ($Variants->{$chr}{$testp}{$strain} eq "C") {
	      if  ( $sequences[$n+1] eq "G") {		 
	#	$ToGetStrains->{$strain} = "CG+";
		$ToGetStrains->{$strain} = "CG"; 
		$found++;
	      }
	      else {
		$ToGetStrains->{$strain} = $Variants->{$chr}{$testp}{$strain};
	      }
	    }
	    elsif  ($Variants->{$chr}{$testp}{$strain} eq "G") {
	      if ( $sequences[$n-1] eq "C") {
		#		  $ToGetStrains->{$strain} = "-CG";
		$ToGetStrains->{$strain} = "CG"; 
		$found++;
	      }
	      else {
		$ToGetStrains->{$strain} = $Variants->{$chr}{$testp}{$strain};
	      }
	    }	    
	    else {	     
		$ToGetStrains->{$strain} = "$Variants->{$chr}{$testp}{$strain}";	     
	    }	    
	  }
	}
      }
      if ($found > 0) {
	print "$chr\t$testp\t$context";	    
	foreach my $strain (sort keys %{$ToGetStrains}) {
	  print "\t$ToGetStrains->{$strain}";
	}	 
      }
      #no sequence variant so print out the sequence for each strain as B6 (as CG)
      
      
      ########################
      # check for methylation
      ########################
      if ($found > 0 ) {	
	if (exists ($Coverage->{$chr}{$testp})) {
	  foreach my $strain (sort keys %{$Strains}) {	  
	    if (exists($Methylation->{$chr}{$testp}{$strain})) {
	      print "\t$Methylation->{$chr}{$testp}{$strain}\t$Coverage->{$chr}{$testp}{$strain}";
	    }
	    else {
	      print "\tNA\tNA";
	    }
	  }	      
	}	
	else {
	  foreach my $strain (sort keys %{$Strains}) {
	    if ($ToGetStrains->{$strain} eq "CG" ) {
	      print "\t0\t0";
	    }
	    else {
	      print "\tNA\tNA";
	    }
	  }
	}
	
	if (exists($Annotations->{$chr}{$testp})) {
	  my $features;
	  foreach my $feature (sort keys %{$Annotations->{$chr}{$testp}}) {
	    #	print "\t$feature\t$Annotations->{$chr}{$testp}{$feature}";
	    $features.=  $feature . ",";
	  }
	  chop ($features);
	  print "\t$features";
	}
	else {	  
	  print "\tNA";	  
	}

#	if (exists($Heritability->{$chr}{$testp})) {
#	print "\t$Heritability->{$chr}{$testp}";
#	}
#	else {
#	print "\tNA";
#	}

if (exists($CpGIslands->{$chr}{$testp})) {
        print "\t$CpGIslands->{$chr}{$testp}";
        }
        else {
	print "\tNA";
        }

if (exists($TSS->{$chr}{$testp})) {
        print "\t$TSS->{$chr}{$testp}";
        }
        else {
	print "\tNA";
        }

if (exists($GeneNames->{$chr}{$testp})) {
        print "\t$GeneNames->{$chr}{$testp}\n";
        }
        else {
	print "\tNA\n";
        }




      }   
    }
  }
}



