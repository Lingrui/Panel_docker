#!usr/bin/perl -w
#DP and dbsnp number of SNPs 
################################################################
#
# perl $0 MH_SNP novo_SNP sample_name
#
################################################################

open IN1,"$ARGV[0]" or die "could not open MH_SNP file\n";
open IN2,"$ARGV[1]" or die "could not open novo_SNP file\n";
$sample=$ARGV[2];

my $mh_dbsnp=0;
my $novo_dbsnp=0;
my $mh_sum=0;
my $novo_sum=0;

open OUT, ">${sample}_tmp.txt";
print OUT "Source\tDepth\n";
while (<IN1>){
	if ($_ !~ /^#/){
		chomp;
		$mh_sum++;
		my @mh_snp = split(/\t/,$_);
		my $dp;
		if($mh_snp[2]=~/rs\d+/){
        	$mh_dbsnp++;
         }
        if($mh_snp[7]=~/DP=(\d+);/){
        	my $mh=$1;
			$dp=join("\t","mh","$mh");
			print OUT "$dp\n";
		}
    }

}

while (<IN2>){
	if ($_ !~ /^#/){
		chomp;
		$novo_sum++;
		my @novo_snp = split(/\t/,$_);
		my $dp;
		if($novo_snp[7]=~/^DP=(\d+);/){
			$dp=join("\t","novo","$1");
			print OUT "$dp\n";
		}
		if($novo_snp[7]=~/;snp138=rs(\d+);/){
			$novo_dbsnp++;
		}
	}
}
close OUT;
my $novo_ratio=$novo_dbsnp/$novo_sum*100;
my $mh_ratio=$mh_dbsnp/$mh_sum*100;
open OUT2, ">${sample}_dbsnp.txt";
print OUT2 "Sample\tnovo_sum\tnovo_dbsnp\tnovo_dbsnp_ratio(%)\tmh_sum\tmh_dbsnp\tmh_dnsnp_ratio(%)\n";
print OUT2 "$sample\t$novo_sum\t$novo_dbsnp\t";
printf OUT2 "%0.2f",$novo_ratio;
print OUT2 "\t$mh_sum\t$mh_dbsnp\t";
printf OUT2 "%0.2f\n",$mh_ratio;

my $rscipts=<<EOF;

library(ggplot2)
dp<-read.table(file="${sample}_tmp.txt",header=T,sep="\t")
png("$sample.snp_depth.png",width = 977, height = 704)
p <- ggplot(dp,aes(Depth))
p+geom_histogram(position='identity',binwidth =2,
  alpha=0.5,
  aes(y=..density..,
  fill=factor(Source))) +
  scale_x_continuous(breaks = seq(0, 200,20),limits = c(0, 200), expand = c(0, 0)) +
  stat_density(geom = 'line',
  position = 'identity',
  aes(colour=factor(Source)))+
  theme_bw()+
  ggtitle("$sample SNPs distribution \\n")
dev.off()
EOF
open RSP,">${sample}.snp_depth.R" or die "can't creat file\n";
print RSP $rscipts."\n";
close RSP;
`/usr/bin/Rscript ${sample}.snp_depth.R`;
