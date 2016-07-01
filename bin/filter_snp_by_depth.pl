#!usr/bin/perl -w
#filter SNPs by depth
################################################################
#
# perl $0 novo_SNP.vcf outfile depth_cutoff
#
################################################################

open IN1,"$ARGV[0]" or die "could not open SNP file\n";
open OUT,">$ARGV[1]" or die "could not open output file\n";
my $cutoff=$ARGV[2];
my $dp=0;
while (<IN1>){
	if ($_ !~ /^#/){
		chomp;
		my @snp = split(/\t/,$_);
        if($snp[7]=~/DP=(\d+);/){
        	$dp=$1;
        	if($dp>=$cutoff){
        		print OUT join"\t",@snp;
        		print OUT "\n";
        	}
		}
    }

}




