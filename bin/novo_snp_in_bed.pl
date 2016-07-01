#!usr/bin/perl -w
#DP and dbsnp number of SNPs 
################################################################
#
# perl $0 BED file novo_SNP novo_in_BED
#
################################################################

open IN1,"$ARGV[0]" or die "could not open BED file\n";
open IN2,"$ARGV[1]" or die "could not open novo_SNP file\n";
open OUT,">$ARGV[2]" or die "could not open output file\n";

my @bed;
while (<IN1>){
	if($_ !~ /^#/){
		chomp;
#my @tmp = split(/\t/,$_);
#my $tmp_bed =join("\t",@tmp);
		push @bed,$_;
	}
}
while (<IN2>){
	if($_ !~ /^#/){
	chomp;
	my $novo_snp=$_;
	my @novo = split(/\t/,$_);
	my $chr_novo = join("","chr",$novo[0]);
	my $pos_novo = $novo[1];
	foreach (@bed){
		chomp;
		my @target = split(/\t/,$_);
#my $chr_tar = $target[0];
#my $sta_tar = $target[1];
#my $end_tar = $target[2];
		if ($chr_novo eq $target[0]){
			if (($target[1]<=$pos_novo)&&($target[2]>=$pos_novo)){
				print OUT "$novo_snp\n";
			}
		}
	}

	}
}


		
