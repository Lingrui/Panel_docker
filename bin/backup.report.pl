#use Switch;
####### agrv[0] input *.candidate.txt
########agrv[1]  annotation
######## argv[2]  report.txt
####### argv[3]  pubmed.info
#############################
$n=0;
open(IN2,"<","$ARGV[1]")||die "cannot open the file: $!\n";
@anno=<IN2>;
open(INPUT,"<","$ARGV[0]")||die"cannot open the file: $!\n";
open(OUTPUT,">","$ARGV[2]")||die"cannot open the file: $!\n";
open(PUBMED,"$ARGV[3]")||die"cannot open the file: $!\n";
###read in all pubmed info
my %pmd=();
while(<PUBMED>){ 
	chomp;
	my $rl=$_;
	if($rl=~/PMID:\s(\d+)/){
		my $pmid=$1;
		$pmd{$pmid}=$rl;
	}
}

if (-s $ARGV[0]){
	my %ptitle=();
	my $pout="";
	my $pid=1;
	$out ='{"disease":"","result":[';
		while(<INPUT>){
			chomp;
			my @candidates = split(/\t/, $_);
			$disease = $candidates[1];
			$symble = $candidates[2];
			$aa_mut = $candidates[12];
			$nt_mut = $candidates[13];
			$rs_id  = $candidates[15];
			$position = $candidates[-1];
			my $pmid  = $candidates[26];
			$n++;
			foreach(@anno){
				undef @w;
				@w=split(/\t/,$_);
				$pos=join(",",$w[0],$w[1],$w[3],$w[4]);
				if ($pos eq $position){

					if ($w[-1] =~ /^0/){$het="杂合"}
					else {$het="纯合"}

					if ($w[8] ne "."){$mut=$w[8]}
					else {$mut=$w[5]}

					if ($mut =~ /exonic/) {$mut_type="外显子区"}
					elsif ($mut eq "intronic") {$mut_type="内含子区"}
					elsif ($mut eq "intergenic"){$mut_type="基因间区"}
					elsif ($mut eq "splicing"){$mut_type="可变剪切"}
					elsif ($mut eq "upstream"){$mut_type="基因上游区域"}
					elsif ($mut eq "downstream"){$mut_type="基因下游区域"}
					elsif ($mut eq "nonsynonymous SNV"){$mut_type="非同义突变"}
					elsif ($mut eq "synonymous SNV"){$mut_type="同义突变"}
					elsif ($mut eq "stopgain"){$mut_type="提前终止"}
					elsif ($mut eq "stoploss"){$mut_type="终止密码子缺失"}
					elsif ($mut eq "unknown"){$mut_type="未知"}
					elsif ($mut =~ /^ncRNA/){$mut_type="非编码RNA"}
					elsif ($mut eq "upstream;downstream"){$mut_type="转录本上下游1kb以外区域"}
					elsif ($mut eq "UTR5"){$mut_type="5′端非编码区"}
					elsif ($mut eq "UTR3"){$mut_type="3′端非编码区"}
					break;
				}
			}
			my @description_output =split (/\t/,$pmd{$pmid});
			$auther = $description_output[0];
			$magzine = $description_output[1];
			$reference_output = "$auther等在$magzine中发表论文表明，$symble基因的$nt_mut突变导致$mut_type$aa_mut，可能与$disease相关，请参阅：<br/>$description_output[-1]<br/>";
#

			$out .='{'.
				'"gene":"'. $symble . '",'.
				'"rs_id":"'. $rs_id . '",'.
				'"position":"' . $position .'",'.
				'"nucleotide":"' . $nt_mut .'",'.
				'"amino":"' . $aa_mut .'",'.
				'"heterozygous":"' . $het . '",'.
				'"type":"'. $mut_type . '",'.
				'"pathogenicity":"'. $disease . '",'.
				'"reference":" '. $reference_output . '"},';

	}#end while
	chop($out);#remove last dot 	
	$out .=']}';
	print OUTPUT "$out\n";
	}
else {
	$out ='{"disease":"","result":[]}';
	print OUTPUT "$out\n";
}





