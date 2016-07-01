use strict;
use warnings;

############################################
#
# Usage: perl $0 sample.SNP.hg19_multianno.txt sample_name outfile
#
############################################

open ANNO,"$ARGV[0]" or die "could not open annotation file\n";
open OUT,">$ARGV[2]" or die "could not create output file\n";

my $prefix=$ARGV[1];
print OUT "基因组不同区域上SNV数目\n\n";
print OUT "Sample\texonic\tintronic\tUTR3\tUTR5\tintergenic\tncRNA_exonic\tncRNA_intronic\tupstream\tdownstream\tsplicing\tncRNA_splicing\n";
my @VAR=(0) x 11;
my @Gtype=(0) x 2;
my @EXO=(0) x 5;
print OUT "$prefix\t";
while(<ANNO>){
    chomp;
    my ($po,$type,$gt)=(split "\t")[5,8,-1];
    if($po =~ /^exonic/){
        $VAR[0]++;
        if($type=~/^synonymous/){
            $EXO[0]++;
        }elsif($type=~/^nonsynonymous/){
            $EXO[1]++;
        }elsif($type=~/stopgain/){
            $EXO[2]++;
        }elsif($type=~/stoploss/){
            $EXO[3]++;
        }elsif($type=~/unknown/){
            $EXO[4]++;
        }

    }elsif($po eq 'intronic'){
        $VAR[1]++;
    }elsif($po eq 'UTR3'){
        $VAR[2]++;
    }elsif($po eq 'UTR5'){
        $VAR[3]++;
    }elsif($po eq 'intergenic'){
        $VAR[4]++;
    }elsif($po =~ /^ncRNA_exonic/){
        $VAR[5]++;
    }elsif($po eq 'ncRNA_intronic'){
        $VAR[6]++;
    }elsif($po eq 'upstream'){
        $VAR[7]++;
    }elsif($po eq 'downstream'){
        $VAR[8]++;
    }elsif($po eq 'splicing'){
        $VAR[9]++;
    }elsif($po eq 'ncRNA_splicing'){
        $VAR[10]++;
    }
    if($gt =~ /^1/){
        $Gtype[1]++;
    }else{
        $Gtype[2]++;
    }
    $Gtype[0]++;
   

}
foreach my $one (@VAR){
    print OUT "$one\t";
}
print OUT "\n\n";
print OUT "----------------------------------------------------------------------\n\n";
print OUT "编码区上不同类型的SNV数目\n\n";

print OUT "Sample\tsynonymous_SNV\tnonsynonymous_SNV\tstopgain\tstoploss\tunknown\n$prefix\t";


foreach my $one (@EXO){
     print OUT "$one\t";
}
print OUT "\n";
my $var_list=join ",",@VAR;
my $exo_list=join ",",@EXO;
my $gt_list=$Gtype[1].",".$Gtype[2];
my $het_rate=sprintf "%0.2f",$Gtype[2]/($Gtype[1]+$Gtype[2])*100;
print OUT "----------------------------------------------------------------------\n\n";
print OUT "SNV基因型分布\n\n";
print OUT "Sample\tAll\tgenotype.Hom\tgenotype.Het\n$prefix\t";
foreach my $one (@Gtype){
    print OUT "$one\t";
}
print OUT "\n";

my $rscript=<<EOF;

library("ggplot2")
library(Cairo)

ColorMatrix<-c("#a50026","#d73027","#f46d43","#fdae61","#fee08b","#ffffbf","#d9ef8b","#a6d96a","#66bd63","#1a9850","#006837");
ColorMatrix = ColorMatrix[seq(11,1,-1)]

dt <- data.frame(values=c($var_list),labels=c("exonic","intronic","UTR3","UTR5","intergenic","ncRNA_exonic","ncRNA_intronic","upstream","downstream","splicing","ncRNA_splicing"));
pie <- ggplot(dt,aes(x="",y=values,fill=labels)) + geom_bar(stat = "identity",width = 3) + coord_polar(theta="y")
pie <- pie + scale_fill_manual(values = ColorMatrix)
pie <- pie + xlab('') + ylab('') + labs(fill="Types") + ggtitle("Genomic distribution of SNPs \\n $prefix")
pie <- pie + theme_bw() + theme(axis.text.x = element_blank(),axis.ticks = element_blank(),panel.grid  = element_blank())
CairoPNG("$prefix\_SNP_variant_statistics.png",width = 977, height = 704)
pie
dev.off()

ColorMatrix = ColorMatrix[seq(8,1,-1)]

dt2<-data.frame(values=c($exo_list),labels=c("synonymous_SNV","nonsynonymous_SNV","stopgain","stoploss","unknown"))
pie <- ggplot(dt2,aes(x="",y=values,fill=labels)) + geom_bar(stat = "identity",width = 3) + coord_polar(theta="y")
pie <- pie + scale_fill_manual(values = ColorMatrix)
pie <- pie + xlab('') + ylab('') + labs(fill="Types") + ggtitle("Exomic distribution of SNPs \\n $prefix")
pie <- pie + theme_bw() + theme(axis.text.x = element_blank(),axis.ticks = element_blank(),panel.grid  = element_blank())
CairoPNG("$prefix\_SNP_exonic_statistics.png",width = 977, height = 704)
pie
dev.off()

dt3<-data.frame(values=c($gt_list),labels=c("Hom","Het"));
his <- ggplot(dt3,aes(x=labels,y=values,fill=labels))+geom_bar(stat = "identity");
his <- his + xlab('') + ylab('') + labs(fill="Genotype") + theme_bw() +
    ggtitle("Het rate $het_rate% \\n $prefix")
pdf("$prefix\_SNP_genotype_statistics.pdf")
his
dev.off()
CairoPNG("$prefix\_SNP_genotype_statistics.png",width = 489, height = 704)
his
dev.off()

cat ("Rscript ${prefix} SNP statistic finished\n")

EOF

#open R,"/data/home/qfj/source/R-3.1.0/bin/R --vanilla --slave" or die $!;
#print R $rscript;
#close R;
open RP,">${prefix}_SNP_statistics.R" or die "could not create ${prefix}_SNP_statistics.R file\n";
print RP $rscript."\n";
close RP;
#`/usr/bin/Rscript SNP_statistics.R`;

close ANNO;
close OUT;

