use strict;
use warnings;

###################################
#
# Usage: perl $0 $QC_folder $sample_name $stat_exonV2_out_table
#
###################################
my $dir = $ARGV[0];
my $prefix = $ARGV[1];
my $prefix1 = "${prefix}_1";
my $prefix2 = "${prefix}_2";
my $raw1 = "$dir/${prefix1}.fq_fastqc/fastqc_data.txt";
my $raw2 = "$dir/${prefix2}.fq_fastqc/fastqc_data.txt";
my $clean1 = "$dir/${prefix1}_trimmed.fq_fastqc/fastqc_data.txt";
my $clean2 = "$dir/${prefix2}_trimmed.fq_fastqc/fastqc_data.txt";

my @qual=();
my @qual2=();
my @lines=();
my $start=0;
my $qual_start=0;
my $length_start=0;
my $mean_sum=0;
my $raw1_Q20=0;
my $raw1_Q30=0;
my $raw1_Q40=0;
my $raw2_Q20=0;
my $raw2_Q30=0;
my $raw2_Q40=0;
my $seqnum1=0;
my $length1=0;
my $gc1=0;
my $raw_Q20=0;
my $raw_Q30=0;
my $error=0;
my @clean_qual=();
my @clean_qual2=();
my $clean1_Q20=0;
my $clean1_Q30=0;
my $clean1_Q40=0;
my $clean2_Q20=0;
my $clean2_Q30=0;
my $clean2_Q40=0;
my $clean_seqnum1=0;
my $clean_gc1=0;
my $clean_Q20=0;
my $clean_Q30=0;
my $target=0;
my $clean_base=0;
my $raw_base=0;
my $raw_depth=0;
my $effective=0;
my $raw_data=0;
my $clean_data=0;
open OUT,">${prefix}\_QC.txt";
print OUT "Sample name\tRaw reads\tRaw data(G)\tRaw depth(x)\tEffective(%)\t";
print OUT "Error(%)\tQ20(%)\tQ30(%)\tGC(%)\tClean reads\tClean data(G)\t";
print OUT "Trimmed_Q20(%)\tTrimmed_Q30(%)\tTrimmed_GC(%)\n";
print OUT "$prefix\t";

if ((-e $raw2) && ($raw1 ne $raw2)){
	open IN1,$raw1;
	@lines=<IN1>;
	for (my $i=0;$i<@lines;$i++){
	
		if ($lines[$i] =~ /^Total\s+Sequences\t(\d+)/){
			$seqnum1=$1;
			#print "reads number in raw PE1 is: $seqnum1\n";
		}

		if ($lines[$i] =~ /^Sequence\s+length\t(\d+)/){
			$length1=$1;
			#print "read length of raw PE1 is: $length1\n";
		}
	
		if ($lines[$i] =~ /\%GC\t(\d+)/){
			$gc1=$1;
			#print "GC content of raw PE1 is: $gc1\n";
		}
#statistic quality of total base in raw PE1
		if ($qual_start == 1){
			if 	($lines[$i] =~ /^\d+/){
				if ($lines[$i] =~ /^\d+-\d+/){
					my @tmp = split(/\t/,$lines[$i]);
					$mean_sum+=$tmp[1];
					$mean_sum+=$tmp[1];
				}
				else{
					my @tmp = split(/\t/,$lines[$i]);
					$mean_sum+=$tmp[1];
				}
			}	
			else{
				$qual_start = 0;	
			}
		
		}
		
		#if ($lines[$i] =~ /Base\s+Mean\s+Median\s+Lower\s+Quantile/){
		if ($lines[$i] =~ /Per\s+base\s+sequence\s+quality\s+pass/){
			$qual_start = 1;
			$i++;
		}
	
#statistic Q20 Q30 percentage of total reads PE1	
		if ($start == 1){
			if ($lines[$i] =~ /^(\d+)\t(.*)/){
				$qual[$1]=$2;
			}	
			else{
				$start=0;
			}
		}
		
		if ($lines[$i] =~ /Quality\s+Count/){
			$start=1;
		}
	}
        close IN1;
        
	open IN2,$raw2;
	@lines=<IN2>;
	
	for (my $i=0;$i<@lines;$i++){
#statistic quality of total base in raw PE2
        if ($qual_start == 1){
            if 	($lines[$i] =~ /^\d+/){
                if ($lines[$i] =~ /^\d+-\d+/){
                    my @tmp = split(/\t/,$lines[$i]);
                    $mean_sum+=$tmp[1];
                    $mean_sum+=$tmp[1];
                }
                else{
                    my @tmp = split(/\t/,$lines[$i]);
                    $mean_sum+=$tmp[1];
                }
            }	
            else{
                $qual_start = 0;	
            }
            
        }
		
		if ($lines[$i] =~ /Per\s+base\s+sequence\s+quality\s+pass/){
			$qual_start=1;
			$i++;
		}
#statistic Q20 Q30 percentage of total reads PE2		
		if ($start == 1){
			if ($lines[$i] =~ /^(\d+)\t(.*)/){
				$qual2[$1]=$2;
			}	
			else{
				$start=0;
			}
		}
		
		if ($lines[$i] =~ /Quality\s+Count/){
			$start=1;
		}
	
	}

	for (my $j=0;$j<@qual;$j++){
		if($j<20){
			if($qual[$j]){
				$raw1_Q20+=$qual[$j];
				$raw1_Q30=$raw1_Q20;
			}
		}
		else {
			if($j<30){
				$raw1_Q30+=$qual[$j];
				$raw1_Q40=$raw1_Q30;
			}
			else{
				$raw1_Q40+=$qual[$j];
			}
		}
	}	
	
	for (my $j=0;$j<@qual2;$j++){
		if($j<20){
			if($qual2[$j]){
				$raw2_Q20+=$qual2[$j];
				$raw2_Q30=$raw2_Q20;
			}
		}
		else {
			if($j<30){
				$raw2_Q30+=$qual2[$j];
				$raw2_Q40=$raw2_Q30;
			}
		
			else{
				$raw2_Q40+=$qual2[$j];
			}
		}
	}

}

$raw_Q20=($raw1_Q40+$raw2_Q40-$raw1_Q20-$raw2_Q20)/($raw1_Q40+$raw2_Q40)*100;
$raw_Q30=($raw1_Q40+$raw2_Q40-$raw1_Q30-$raw2_Q30)/($raw1_Q40+$raw2_Q40)*100;
$error=10**(-0.05*$mean_sum/$length1)*100;
#print "Raw Q20: $raw_Q20\n";
#print "Raw Q30: $raw_Q30\n";
print "Error of raw data: $error\n";

if ((-e $clean2) && ($clean1 ne $clean2)){
	open IN3,$clean1;
	@lines=<IN3>;
	for (my $i=0;$i<@lines;$i++){
	
		if ($lines[$i] =~ /^Total\s+Sequences\t(\d+)/){
			$clean_seqnum1=$1;
		}
	
		if ($lines[$i] =~ /\%GC\t(\d+)/){
			$clean_gc1=$1;
		}
		
#statistic clean Q20 Q30 percentage of total reads PE1	
		if ($start == 1){
			if ($lines[$i] =~ /^(\d+)\t(.*)/){
				$clean_qual[$1]=$2;
			}	
			else{
				$start=0;
			}
		}
		
		if ($lines[$i] =~ /Quality\s+Count/){
			$start=1;
		}
#statistic clean base number of PE1		
		if ($length_start ==1){
			if ($lines[$i] =~ /^(\d+)-\d+\t(.*)/){
				$clean_base+=$1*$2;
			}	
			else{
				$length_start=0;
			#print "PE1 clean base number: $clean_base\n";
			}
		
		}
		
		if ($lines[$i] =~ /Length\s+Count/){
			$length_start=1;
		}
	}
	close IN3;
	
	
	open IN4,$clean2;
	@lines=<IN4>;
	for (my $i=0;$i<@lines;$i++){

#statistic clean Q20 Q30 percentage of total reads PE2			
		if ($start == 1){
			if ($lines[$i] =~ /^(\d+)\t(.*)/){
				$clean_qual2[$1]=$2;
			}	
			else{
				$start=0;
			}
		}
		
		if ($lines[$i] =~ /Quality\s+Count/){
			$start=1;
		}
#statistic clean base number of PE				
		if ($length_start ==1){
			if ($lines[$i] =~ /^(\d+)-\d+\t(.*)/){
				$clean_base+=$1*$2;
			}	
			else{
				$length_start=0;
				#print "PE clean base number: $clean_base\n";
			}
		
		}
		
		if ($lines[$i] =~ /Length\s+Count/){
			$length_start=1;
		}
	}
	close IN4;

	for (my $j=0;$j<@clean_qual;$j++){
		if($j<20){
			if($clean_qual[$j]){
				$clean1_Q20+=$clean_qual[$j];
				$clean1_Q30=$clean1_Q20;
			}
		}
		else {
			if($j<30){
				if($clean_qual[$j]){
					$clean1_Q30+=$clean_qual[$j];
					$clean1_Q40=$clean1_Q30;
				}
			}
			else{
				$clean1_Q40+=$clean_qual[$j];
			}
		}
	}	
	
	for (my $j=0;$j<@clean_qual2;$j++){
		if($j<20){
			if($clean_qual2[$j]){
				$clean2_Q20+=$clean_qual2[$j];
				$clean2_Q30=$clean2_Q20;
			}
		}
		else {
			if($j<30){
				$clean2_Q30+=$clean_qual2[$j];
				$clean2_Q40=$clean2_Q30;
			}
		
			else{
				$clean2_Q40+=$clean_qual2[$j];
			}
		}
	}
}

$clean_Q20=($clean1_Q40+$clean2_Q40-$clean1_Q20-$clean2_Q20)/($clean1_Q40+$clean2_Q40)*100;
$clean_Q30=($clean1_Q40+$clean2_Q40-$clean1_Q30-$clean2_Q30)/($clean1_Q40+$clean2_Q40)*100;		

#print "Clean Q20: $clean_Q20\n";
#print "Clean Q30: $clean_Q30\n";
		
open INFO,"$ARGV[2]";
while (<INFO>){
	chomp;
	if ($_ =~ /Target_region_size\(bp\):\t(\d+)/){
		$target=$1;
		print "Target region is: $target\n";
	}
}
close INFO;

$raw_base=$seqnum1*$length1*2;
$raw_data=$raw_base/1000000000;
$raw_depth=$seqnum1*$length1*2/$target;
$effective=$clean_base/$raw_base*100;
$clean_data=$clean_base/1000000000;

print OUT "$seqnum1\t";
printf OUT "%0.2f\t",$raw_data;
printf OUT "%0.2f\t",$raw_depth;
printf OUT "%0.2f\t",$effective;
printf OUT "%0.2f\t",$error;
printf OUT "%0.2f\t",$raw_Q20;
printf OUT "%0.2f\t",$raw_Q30;
printf OUT "%0.2f\t",$gc1;
print OUT "$clean_seqnum1\t";
printf OUT "%0.2f\t",$clean_data;
printf OUT "%0.2f\t",$clean_Q20;
printf OUT "%0.2f\t",$clean_Q30;
printf OUT "%0.2f\n",$clean_gc1;

