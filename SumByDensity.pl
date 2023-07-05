#!/usr/bin/perl
use strict;
use Getopt::Std;

our($opt_h,$opt_p,$opt_H, $opt_m, $opt_t, $opt_d, $opt_t,$opt_f, $opt_a, $opt_z, $opt_C);

getopts("hpf:d:Hamz:Ct:");
&usage if $opt_h;
sub usage {
  die "SumByDensity.pl -f <file> 
-z compressed file
-d density (default is 10)
-m sum by density of methylated sites rather than number of CpGs
-C extract mean of all strains against outgroup Cast
sum from 0-10, 10-20 etc
-t type (CG default or CH )
-H omit header\n";
}
sub weighted_mean {
    my ($means, $subjects) = @_;

    my $total_subjects = 0;
    my $weighted_sum = 0;

    for (my $i = 0; $i < scalar @$means; $i++) {
     # print "n = $subjects->[$i], mean = $means->[$i]\n";
        $total_subjects += $subjects->[$i];
        $weighted_sum += $means->[$i] * $subjects->[$i];
    }

    my $weighted_mean = "NA";
    
      if ($total_subjects > 0) {
	$weighted_mean = $weighted_sum /$total_subjects;
      }
    
    return $weighted_mean;
}

sub calculate_mean {
    my @data = @_;
    my $sum = 0;
    my $count = scalar @data;
    foreach my $number (@data) {
      $sum += $number;
     
      }
    my $mean = "NA";
    if ($count >0) {
      $mean = $sum / $count;
    }
    # print "$mean, $sum, $count\n";
    return $mean;
}

sub calculate_similarity {
    my ($array1_ref, $array2_ref) = @_;
    my $length = scalar(@$array1_ref);
    my $same_state_count = 0;
    
    for (my $i = 0; $i < $length; $i++) {
        if ($array1_ref->[$i] == $array2_ref->[$i]) {
            $same_state_count++;
        }
    }
    
    my $probability = $same_state_count / $length;
    
    return ($same_state_count, $probability);
}

sub calculate_ratio {
    my ($array1_ref, $array2_ref) = @_;
    my $length = scalar(@$array1_ref);
    my $count_array1 = 0;
    my $count_array2 = 0;
    
    for (my $i = 0; $i < $length; $i++) {
        if ($array1_ref->[$i] == 1) {
            $count_array1++;
        }
        if ($array2_ref->[$i] == 1) {
            $count_array2++;
        }
    }
    
    my $ratio = $count_array1 / $count_array2;
     if ($ratio > 1) {
      $ratio = 1/$ratio;
    }
    return ( $ratio, $count_array1 , $count_array2);
}


sub calculate_ratio_from_counts {
    my ($both_1_count, $both_0_count, $either_1_count) = @_;
    
    # Calculate the count of 1s in the first array
    my $array1_1_count = ($both_1_count + $either_1_count) / 2;
    
    # Calculate the count of 1s in the second array
    my $array2_1_count = $both_1_count - $array1_1_count;
    
    # Calculate the ratio
    my $ratio = $array1_1_count / $array2_1_count;
    if ($ratio > 1) {
      $ratio = 1/$ratio;
    }
    return $ratio;
}

sub calculate_probability_from_counts {
    my ($both_1_count, $both_0_count, $either_1_count) = @_;
    
    my $total_count = $both_1_count + $both_0_count + $either_1_count;
    
    my $probability = ($both_1_count + $both_0_count) / $total_count;
    
    return $probability;
}




sub prop_occurrence1 {
  my ($array1_ref, $array2_ref) = @_;
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
  my $count_pq = 0;
  for my $i (0..$#array1) {
    $count_pq++ if ($array1[$i] == $p && $array2[$i] == $q);
    $count_pq++ if ($array1[$i] == $q && $array2[$i] == $p);
  }

  # Calculate the proportions
  my $prop_p_array1 = $count_p_array1 / scalar(@array1);
  my $prop_q_array1 = $count_q_array1 / scalar(@array1);
  my $prop_p_array2 = $count_p_array2 / scalar(@array2);
  my $prop_q_array2 = $count_q_array2 / scalar(@array2);
  my $prop_pq = $count_pq / scalar(@array1);
#print "\nP1\tA1\t@array1\nP1A2\t@array2\n";
#print "\nP1\t$prop_p_array1\t$prop_p_array2\t$prop_q_array1\t$prop_q_array2\t$prop_pq\n";
  return ($prop_p_array1, $prop_q_array1, $prop_p_array2, $prop_q_array2, $prop_pq);
}



sub prop_occurrence {
 my  ($array1_ref, $array2_ref) = @_;
 my $p = 1;
 my $q = 0;
 # note that the arrays will be the same size for each comparison
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
 

 my $size1 =  scalar(@array1);;
 my $size2 = scalar(@array2);
 # print "\nP2\t$prop_p_array1\t$prop_p_array2\t$prop_q_array1\t$prop_q_array2\t$prop_pp\t$prop_qq\t s1 $size1\ts2 $size2\n";
  return ($prop_p_array1,  $prop_p_array2,$prop_q_array1,$prop_q_array2, $prop_pp, $prop_qq);
}

#seg     s1      s2	n      CG.all  CG.me.0 CG.me.1 mutant.me.0     mutant.me.1
my $density = defined ($opt_d) ? $opt_d : 10;
my $type =  defined ($opt_t) ? $opt_t : "CG";
my $Features = {};
if ($type eq "CH") {
  $Features->{'CH.all'}++;
  $Features->{'CH.me.0'}++;
  $Features->{'CH.me.1'}++;
}
else {
  $Features->{'CG.all'}++;
  $Features->{'CG.me.0'}++;
  $Features->{'CG.me.1'}++;
}
$Features->{'s1.me'}++;
$Features->{'s2.me'}++;
$Features->{'mutant.me.0'}++;
$Features->{'mutant.me.1'}++;
$Features->{'ratio'}++;
#$Features->{'ratio'}++;
$Features->{'prob'}++;
$Features->{'me.diff'}++;
$Features->{'h2'}++;

my $Annotations= {};
$Annotations->{"exon"}++;
$Annotations->{"fiveprime"}++;
$Annotations->{"intron"}++;
$Annotations->{"threeprime"}++;
$Annotations->{"CTCF"}++;
$Annotations->{"Enhancer"}++;
$Annotations->{"Het"}++;
$Annotations->{"K4m3"}++;
$Annotations->{"Promoter"}++;
$Annotations->{"Transcription"}++;
$Annotations->{"enhD"}++;
$Annotations->{"enhP"}++;
$Annotations->{"prom"}++;
$Annotations->{"LINE"}++;
$Annotations->{"LINE_L1"}++;
$Annotations->{"LINE_L2"}++;
$Annotations->{"LTR"}++;
$Annotations->{"SINE"}++;
$Annotations->{"SINE_B2"}++;
$Annotations->{"SINE_B4"}++;
$Annotations->{"Simple_repeat"}++;
$Annotations->{"tss"}++;
$Annotations->{"cgi"}++;


unless ($opt_H)  {
  if ($opt_C) {
    print "seg\tscast";
    foreach my $feature (sort keys %{$Features}) {
      print "\t$feature";
    }
    foreach my $annotation (sort keys %{$Annotations}) {
      print "\t$annotation";
    }
     print "\n";
  }
  else {
    print "seg\ts1\ts2";
    foreach my $feature (sort keys %{$Features}) {
      print "\t$feature";
    }
    
    foreach my $annotation (sort keys %{$Annotations}) {
      print "\t$annotation";
    }    
    print "\n";
  } 
}


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

my $Data = {};
my $Joint = {};

my $H2 = {};
my $Header = {};
#open (FILE, $opt_f) || die "Can't open file $opt_f\n";
my $head  = <$file_handle>;
chomp ($head);
my @data = split (/\t/, $head);
for (my $n = 0; $n <= $#data; $n++) {
  $Header->{$data[$n]} = $n;
}
#foreach my $entry (sort keys %{$Header}) {
#  print "$entry\n";
#}

while (<$file_handle>) {
  chomp;
  my @data = split;
  
  
  my $chr = $data[$Header->{'chr'}];
  my $seg = $data[$Header->{'seg'}];
  my $s1 =   $data[$Header->{'s1'}];
  my $s2 =  $data[$Header->{'s2'}];
  my $n =   $data[$Header->{'n'}];
  my $m1m2 =    $data[$Header->{'m1.m2'}];
  my $me_ratio =    $data[$Header->{'me.ratio'}];
  my $me_same =    $data[$Header->{'me.same'}];
  my $seq1seq2 =    $data[$Header->{'seq1.seq2'}];
  my $seqratio =  $data[$Header->{'seq.ratio'}];
  my $seqsame =    $data[$Header->{'seq.same'}];
  my $CGall =   $data[$Header->{'CG.all'}];
  my $CGme0 =    $data[$Header->{'CG.me.0'}];
  my $CGme1 =    $data[$Header->{'CG.me.1'}];
  my $CHall =   $data[$Header->{'CH.all'}];
  my $CHme0 =    $data[$Header->{'CH.me.0'}];
  my $CHme1 =    $data[$Header->{'CH.me.1'}];
  my $bothme0 =   $data[$Header->{'both.me.0'}];
  my $bothme1 =    $data[$Header->{'both.me.1'}];
  my $mutantme0 =    $data[$Header->{'mutant.me.0'}];
  my $mutantme1 =    $data[$Header->{'mutant.me.1'}];
  my $oneme1 =   $data[$Header->{'one.me.1'}];
  my $s1me0 = $data[$Header->{'s1.me.0'}];
  my $s1me1 =  $data[$Header->{'s1.me.1'}];
  my $s2me0 =   $data[$Header->{'s2.me.0'}];
  my $s2me1 =    $data[$Header->{'s2.me.1'}];
  my $exon =    $data[$Header->{'exon'}];
  my $fiveprime =    $data[$Header->{'fiveprime'}];
  my $intron =    $data[$Header->{'intron'}];
  my $threeprime =    $data[$Header->{'threeprime'}];
  my $CTCF =    $data[$Header->{'CTCF'}];
  my $Enhancer =   $data[$Header->{'Enhancer'}];
  my $Het =  $data[$Header->{'Het'}];
  my $K4m3 = $data[$Header->{'K4m3'}];
  my $Promoter =   $data[$Header->{'Promoter'}];
  my $Transcription =   $data[$Header->{'Transcription'}];
  my $enhD =    $data[$Header->{'enhD'}];
  my $enhP =    $data[$Header->{'enhP'}];
  my $prom =    $data[$Header->{'prom'}];
  my $LINE =    $data[$Header->{'LINE'}];
  my $LINE_L1 =   $data[$Header->{'LINE_L1'}];
  my $LINE_L2 =    $data[$Header->{'LINE_L2'}];
  my $LTR =    $data[$Header->{'LTR'}];
  my $SINE =    $data[$Header->{'SINE'}];
  my $SINE_B2 =    $data[$Header->{'SINE_B2'}];
  my $SINE_B4 =    $data[$Header->{'SINE_B4'}];
  my $Simple_repeat =    $data[$Header->{'Simple_repeat'}];

  my $tss =    $data[$Header->{'tss'}];
  my $cgi = $data[$Header->{'cgi'}];
  
  my $meanh2 =    $data[$Header->{'mean.h2'}];
  my $nh2 = $data[$Header->{'n.h2'}];
  
  #print "mean h2 $meanh2, n h2 $nh2\n";
  

##############
## sum by methylation sites if $opt_m
##############
if ($opt_m) {
$n =  $CGme1;
}   
  my $d = $n;
  my $x = $d % $density;
  my $seg = $d - $x;
  if ($type eq "CG") {
    unless ($CGall eq "NA") {
      $Data->{$seg}{$s1}{$s2}{'CG.all'}+= $CGall;
      if ($s1 eq "cast") {
	$Joint->{$seg}{'CG.all'}+= $CGall;
      }
    }
    unless ($CGme0 eq "NA") {
      $Data->{$seg}{$s1}{$s2}{'CG.me.0'}+= $CGme0 ;
      if ($s1 eq "cast") {
	$Joint->{$seg}{'CG.me.0'}+= $CGme0;
      }
    }
    unless ($CGme1 eq "NA") {
      $Data->{$seg}{$s1}{$s2}{'CG.me.1'}+= $CGme1;
      if ($s1 eq "cast") {
	$Joint->{$seg}{'CG.me.1'}+= $CGme1;
      }
    }
  }
  elsif ($type eq "CH") {
    unless ($CHall eq "NA") {
      $Data->{$seg}{$s1}{$s2}{'CH.all'}+= $CHall;
      if ($s1 eq "cast") {
	$Joint->{$seg}{'CH.all'}+= $CHall;
      }
    }
    unless ($CHme0 eq "NA") {
      $Data->{$seg}{$s1}{$s2}{'CH.me.0'}+= $CHme0 ;
      if ($s1 eq "cast") {
	$Joint->{$seg}{'CH.me.0'}+= $CHme0;
      }
    }
    unless ($CHme1 eq "NA") {
      $Data->{$seg}{$s1}{$s2}{'CH.me.1'}+= $CHme1;
      if ($s1 eq "cast") {
	$Joint->{$seg}{'CH.me.1'}+= $CHme1;
      }
    }
  }
  
  unless ($mutantme0 eq "NA") {
    $Data->{$seg}{$s1}{$s2}{'mutant.me.0'}+= $mutantme0;
    if ($s1 eq "cast") {
      $Joint->{$seg}{'mutant.me.0'}+= $mutantme0;
    }
  }
  unless ($mutantme1 eq "NA") {
    $Data->{$seg}{$s1}{$s2}{'mutant.me.1'}+= $mutantme1;
    if ($s1 eq "cast") {
      $Joint->{$seg}{'mutant.me.1'}+= $mutantme1;
    }
  }
  unless ($s1me1 eq "NA") {
    $Data->{$seg}{$s1}{$s2}{'s1.me'}+= $s1me1;
     if ($s1 eq "cast") {
       $Joint->{$seg}{'s1.me'}+=$s1me1;
     }
  }
  unless ($s2me1 eq "NA") {
    $Data->{$seg}{$s1}{$s2}{'s2.me'}+= $s2me1;
    if ($s1 eq "cast") {
      $Joint->{$seg}{'s2.me'}+=$s2me1;
    }
  }
  
  unless($oneme1 eq "NA") {
    $Data->{$seg}{$s1}{$s2}{'me.diff'} += $oneme1;
    if ($s1 eq "cast") {
      $Joint->{$seg}{'me.diff'}+=$oneme1;
    }
  }
  
  unless ($bothme0 eq "NA") {
    $Data->{$seg}{$s1}{$s2}{'both.me.0'}+= $bothme0;
    if ($s1 eq "cast") {
      $Joint->{$seg}{'both.me.0'}+=$bothme0;
    }
  }
  unless ($bothme1 eq "NA") {
    $Data->{$seg}{$s1}{$s2}{'both.me.1'}+=$bothme1;
     if ($s1 eq "cast") {
       $Joint->{$seg}{'both.me.1'}+=$bothme1;
     }
  }
  unless ($oneme1 eq "NA") {
    $Data->{$seg}{$s1}{$s2}{'one.me.1'}+= $oneme1;
    if ($s1 eq "cast") {
      $Joint->{$seg}{'one.me.1'}+=$oneme1;
    }
  }
  if (exists($Header->{'n.h2'})) {
    unless ($nh2 eq "NA") {
      push (@{$H2->{$seg}{$s1}{$s2}{mean}}, $meanh2);
      push (@{$H2->{$seg}{$s1}{$s2}{n}}, $nh2);
      #  print "h2 @{$H2->{$seg}{$s1}{$s2}{mean}} n @{$H2->{$seg}{$s1}{$s2}{n}} \n";
    }
  }
   
  
  $Data->{$seg}{$s1}{$s2}{'cgi'}+=$cgi;
  $Data->{$seg}{$s1}{$s2}{'CTCF'}+=$CTCF;
  $Data->{$seg}{$s1}{$s2}{'Enhancer'}+=$Enhancer;
  $Data->{$seg}{$s1}{$s2}{'enhD'}+=$enhD;
  $Data->{$seg}{$s1}{$s2}{'enhP'}+= $enhP;
  $Data->{$seg}{$s1}{$s2}{'exon'}+=$exon;
  $Data->{$seg}{$s1}{$s2}{'fiveprime'}+=$fiveprime;
  $Data->{$seg}{$s1}{$s2}{'Het'}+=$Het;
  $Data->{$seg}{$s1}{$s2}{'intron'}+=$intron;
  $Data->{$seg}{$s1}{$s2}{'K4m3'}+=$K4m3;
  $Data->{$seg}{$s1}{$s2}{'LINE_L1'}+=$LINE_L1;
  $Data->{$seg}{$s1}{$s2}{'LINE_L2'}+=$LINE_L2;
  $Data->{$seg}{$s1}{$s2}{'LINE'}+=$LINE;
  $Data->{$seg}{$s1}{$s2}{'LTR'}+=$LTR;
  $Data->{$seg}{$s1}{$s2}{'prom'}+=$prom;
  $Data->{$seg}{$s1}{$s2}{'Promoter'}+=$Promoter;
  $Data->{$seg}{$s1}{$s2}{'Simple_repeat'}+=$Simple_repeat;
  $Data->{$seg}{$s1}{$s2}{'SINE_B2'}+=$SINE_B2;
  $Data->{$seg}{$s1}{$s2}{'SINE_B4'}+=$SINE_B4;
  $Data->{$seg}{$s1}{$s2}{'SINE'}+=$SINE;
  $Data->{$seg}{$s1}{$s2}{'threeprime'}+=$threeprime;
  $Data->{$seg}{$s1}{$s2}{'Transcription'}+=$Transcription;
  $Data->{$seg}{$s1}{$s2}{'tss'}+=$tss;

  if ($opt_C) {
    if ($s1 eq "cast") {
      $Joint->{$seg}{'cgi'}+=$cgi;
      $Joint->{$seg}{'CTCF'}+=$CTCF;
      $Joint->{$seg}{'Enhancer'}+=$Enhancer;
      $Joint->{$seg}{'enhD'}+=$enhD;
      $Joint->{$seg}{'enhP'}+= $enhP;
      $Joint->{$seg}{'exon'}+=$exon;
      $Joint->{$seg}{'fiveprime'}+=$fiveprime;
      $Joint->{$seg}{'Het'}+=$Het;
      $Joint->{$seg}{'intron'}+=$intron;
      $Joint->{$seg}{'K4m3'}+=$K4m3;
      $Joint->{$seg}{'LINE_L1'}+=$LINE_L1;
      $Joint->{$seg}{'LINE_L2'}+=$LINE_L2;
      $Joint->{$seg}{'LINE'}+=$LINE;
      $Joint->{$seg}{'LTR'}+=$LTR;
      $Joint->{$seg}{'prom'}+=$prom;
      $Joint->{$seg}{'Promoter'}+=$Promoter;
      $Joint->{$seg}{'Simple_repeat'}+=$Simple_repeat;
      $Joint->{$seg}{'SINE_B2'}+=$SINE_B2;
      $Joint->{$seg}{'SINE_B4'}+=$SINE_B4;
      $Joint->{$seg}{'SINE'}+=$SINE;
      $Joint->{$seg}{'threeprime'}+=$threeprime;
      $Joint->{$seg}{'Transcription'}+=$Transcription;
      $Joint->{$seg}{'tss'}+=$tss;
    }
  }

  
}
if ($opt_C) {
  
  foreach my $seg (sort {$a <=> $b } keys %{$Joint}) {   
    print "$seg\tcast";
    foreach my $feature (sort keys %{$Features}) {
      if (exists($Joint->{$seg}{$feature})) {
	print "\t$Joint->{$seg}{$feature}";
      }
      else {
	if ($feature eq "ratio") {
	  if (exists($Joint->{$seg}{'s1.me'}) && exists ($Joint->{$seg}{'s2.me'}))  {
	    my $ratio = "NA";
	    my $s1n = $Joint->{$seg}{'s1.me'};
	    my $s2n = $Joint->{$seg}{'s2.me'};
	    if ($s2n >= $s1n && $s2n > 0 ) {		
	      $ratio = $s1n/$s2n;
	    }
	    elsif ($s1n >= $s2n && $s1n > 0 ) {		
	      $ratio = $s2n/$s1n;
	    }
	    printf "\ts1n $s1n s2n $s2n\tRATIO %3.4f", $ratio;
	  }
	  else {
	    print "\tNA";
	  }
	}
	elsif ($feature eq "prob") {	
	  my  $bothme0 =   $Joint->{$seg}{'both.me.0'};         	      
	  my $bothme1 = $Joint->{$seg}{'both.me.1'};	      
	  my  $oneme1 =  $Joint->{$seg}{'one.me.1'};	
	  my @array1;
	  my @array2;
	 ;
	  if ($bothme1 + $bothme0 + $oneme1 > 0) {
	    my $prob_from_counts =  calculate_probability_from_counts($bothme1 , $bothme0 ,$oneme1);

	    print "\n";
	     printf "%3.4f\n", $prob_from_counts;
	   
	    print "\n";
	  }
	  else {
	    print "\tNA";
	  } 
	}
	elsif ($feature eq "h2") {
	  print "\tNA";
	}
      }
    }
          
    foreach my $annotation (sort keys %{$Annotations}) {	    
      if (exists ($Joint->{$seg}{$annotation})) {
	print "\t$Joint->{$seg}{$annotation}";
      }
      else {
	print "\t0";
      }
    }
    print "\n";
  }
}

    
else {
  
  foreach my $seg (sort {$a <=> $b } keys %{$Data}) {
    foreach my $s1 (sort keys %{$Data->{$seg}}) {
      foreach my $s2 (sort keys %{$Data->{$seg}{$s1}}) {
	    print "$seg\t$s1\t$s2";
	    foreach my $feature (sort keys %{$Features}) {
	      if (exists($Data->{$seg}{$s1}{$s2}{$feature})) {
		print "\t$Data->{$seg}{$s1}{$s2}{$feature}";
	      }
	      else {
		if ($feature eq "ratio") {
		  if (exists($Data->{$seg}{$s1}{$s2}{'s1.me'}) && exists ($Data->{$seg}{$s1}{$s2}{'s2.me'}))  {
		    my $ratio = "NA";
		    my $s1n = $Data->{$seg}{$s1}{$s2}{'s1.me'};
		    my $s2n = $Data->{$seg}{$s1}{$s2}{'s2.me'};
		    if ($s2n >= $s1n && $s2n > 0 ) {		
		      $ratio = $s1n/$s2n;
		    }
		    elsif ($s1n >= $s2n && $s1n > 0 ) {		
		      $ratio = $s2n/$s1n;
		    }
		    
		    printf "\t%3.4f", $ratio;
		  }
		  else {
		    print "\tNA";
		  }
		}
		elsif ($feature eq "prob") {
		  
		  my  $bothme0 =   $Data->{$seg}{$s1}{$s2}{'both.me.0'};         	      
		  my $bothme1 = $Data->{$seg}{$s1}{$s2}{'both.me.1'};	      
		  my  $oneme1 =  $Data->{$seg}{$s1}{$s2}{'one.me.1'};
		  
		  
		  my @array1;
		  my @array2;
		  for (my $n = 1; $n <= $bothme0 ; $n++) {
		    push (@array1, 0);
		    push (@array2, 0);
		  }
		  
		  for (my $n = 1; $n <= $bothme1 ; $n++) {
		    push (@array1, 1);
		    push (@array2, 1);
		  }
		  for (my $n = 1; $n <= $oneme1 ; $n++) {
		    push (@array1, 1);
		    push (@array2, 0);
		  }
		  if ( scalar ($#array1)  > 0 && scalar ($#array2) > 0) { 
		    my ($p_array1,  $p_array2,$q_array1,$q_array2, $prop_p, $prop_q) = prop_occurrence(\@array1, \@array2 );
		     my ($prop_p_array1, $prop_q_array1, $prop_p_array2, $prop_q_array2, $prop_pq) = prop_occurrence1(\@array1, \@array2 );;
		    my $equal = $prop_p + $prop_q;
		    printf "\t%3.4f", $equal;
		  }
		  else {
		    print "\tNA";
		  }
		  
		}
		elsif ($feature eq "h2") {
		  if (defined( $H2->{$seg}{$s1}{$s2}{mean})) {
		    my @means =  @{$H2->{$seg}{$s1}{$s2}{mean}};
		    my @subjects =  @{$H2->{$seg}{$s1}{$s2}{n}};
		    #!print "\n@array\n";
		    my $mean = weighted_mean (\@means, \@subjects);
		    if ($mean ne "NA") {
		      printf "\t%3.3f",$mean ;
		    }
		    else {
		      print "\tNA";
		    }
		  }
		  else {
		    print "\tNA";
		  }
		}
	      }
	    }
	    	    
	    foreach my $annotation (sort keys %{$Annotations}) {	    
	      if (exists ($Data->{$seg}{$s1}{$s2}{$annotation})) {
		print "\t$Data->{$seg}{$s1}{$s2}{$annotation}";
	      }
	      else {
		print "\t0";
	      }
	    }
	    print "\n";
	  }     
      
    }
    
  }
    
}
