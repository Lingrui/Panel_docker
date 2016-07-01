use strict;
use warnings;

###################################
#
# Usage: perl $0 table_annovar_out sample_name snp_or_indel
#
###################################


open IN,"$ARGV[0]";
my $prefix="$ARGV[1]";
my $type = "$ARGV[2]";
my @out=(0) x 2;
my @trans_novel=(0) x 2;###
my @trans=(0) x 2;###

open OUT,">$prefix\_$type\_table.txt";
print OUT "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tGeneName\tFunc\tGeneDetail\tExonicFunc\tAAChange\tesp6500si_all\t1000g2015aug_all\t1000g2015aug_eas\tsnp138\tSIFT_score\tPolyphen2_score\tMutationTaster_score\tFATHMM_score\tLR_score\tCADD\tFORMAT\tSL003\n";
while(<IN>){
	if ($_ !~ /^#/){
	chomp;
	my @tmp=split "\t";
    my $type=$tmp[7];
    my $var= join("",$tmp[3],$tmp[4]);####
    print OUT join "\t",@tmp[0..7];
    my $ot="";
    if($type=~/Gene.refGene=(.*?);/){
        $ot=$1;
        print OUT "$ot\t";
    }
    if($type=~/Func.refGene=(.*?);/){
        $ot=$1;
        print OUT "$ot\t";
    }
    if($type=~/GeneDetail.refGene=(.*?);/){
        $ot=$1;
        print OUT "$ot\t";
    }
    if($type=~/ExonicFunc.refGene=(.*?);/){
        $ot=$1;
        print OUT "$ot\t";
    }
    if($type=~/AAChange.refGene=(.*?);/){
        $ot=$1;
        print OUT "$ot\t";
    }
    if($type=~/esp6500siv2_all=(.*?);/){
        $ot=$1;
        print OUT "$ot\t";
    }
    if($type=~/1000g2015aug_all=(.*?);/){
        $ot=$1;
        print OUT "$ot\t";
    }
    if($type=~/1000g2015aug_eas=(.*?);/){
        $ot=$1;
        print OUT "$ot\t";
    }
    if($type=~/snp138=(.*?);/){
        $ot=$1;
        print OUT "$ot\t";
    }
    if($type=~/SIFT_score=(.*?);/){
        $ot=$1;
        print OUT "$ot\t";
    }
    if($type=~/Polyphen2_HDIV_score=(.*?);/){
        $ot=$1;
        print OUT "$ot\t";
    }
    if($type=~/MutationTaster_score=(.*?);/){
        $ot=$1;
        print OUT "$ot\t";
    }
    if($type=~/FATHMM_score=(.*?);/){
        $ot=$1;
        print OUT "$ot\t";
    }
    if($type=~/LR_score=(.*?);/){
        $ot=$1;
        print OUT "$ot\t";
    }
    if($type=~/CADD_raw=(.*?);/){
        $ot=$1;
        print OUT "$ot\t";
    }
    print OUT join "\t",@tmp[8..9];
    print OUT "\n";
    #check is not in dbsnp
	if($type=~/snp138=(.*?);/){
        my $check=$1;
		if($check ne '.'){
            $out[1]++;
		}else{
            $out[0]++;
            if($var eq 'AG'|$var eq 'TC'|$var eq 'GA'|$var eq 'CT'){ ###
                $trans_novel[0]++;
            }else{
                $trans_novel[1]++;
                }
        }
    }


    if($var eq 'AG'|$var eq 'TC'|$var eq 'GA'|$var eq 'CT'){
        $trans[0]++;
    }else{
        $trans[1]++;
    }
}
}

my $no_list=join ",",@out;
#my $dbSNP_rate=sprintf "%0.2f",$out[1]/($out[0]+$out[1])*100;
my $temp_dbSNP_rate=($out[0]+$out[1])?($out[1]/(($out[0]+$out[1])*100)):0;
my $dbSNP_rate=sprintf "%0.2f",$temp_dbSNP_rate;

my $tn_list=join ",",@trans_novel;

my $temp_novel_ts=$trans_novel[1]?$trans_novel[0]/$trans_novel[1]:$trans_novel[0];
my $novel_ts=sprintf "%0.2f",$temp_novel_ts;
#my $novel_ts=sprintf "%0.2f",$trans_novel[0]/$trans_novel[1];
my $temp_novel_ts_table=$trans_novel[1]?$trans_novel[0]/$trans_novel[1]:$trans_novel[0];
my $novel_ts_table=sprintf "%0.6f",$temp_novel_ts_table;
#my $novel_ts_table=sprintf "%0.6f",$trans_novel[0]/$trans_novel[1];

my $tstv_list=join ",",@trans;

#my $ts=sprintf "%0.2f",$trans[0]/$trans[1];
my $temp_ts=$trans[1]?$trans[0]/$trans[1]:$trans[0];
my $ts=sprintf "%0.2f",$temp_ts;

#my $ts_table=sprintf "%0.6f",$trans[0]/$trans[1];
my $temp_ts_table=$trans[1]?$trans[0]/$trans[1]:$trans[0];
my $ts_table=sprintf "%0.6f",$temp_ts_table;

#my $novel_rate=sprintf "%0.6f",$out[0]/($out[0]+$out[1]);
my $temp_novel_rate=($out[0]+$out[1])?($out[0]/($out[0]+$out[1])):0;
my $novel_rate=sprintf "%0.6f",$temp_novel_rate;

open OUT2,">$prefix\_$type\_tstv.txt";
print OUT2 "----------------------------------------------------------------------\n\n";
print OUT2 "新的碱基置换数目及比例\n\n";
print OUT2 "Sample\tnovel\tnovel_proportion\n";
print OUT2 "$prefix\t$out[0]\t$novel_rate\n";

if($type eq 'SNP'){
    print OUT2 "----------------------------------------------------------------------\n\n";
    print OUT2 "转换和颠换类型分布\n\n";
    print OUT2 "Sample\tnovel_ts\tnovel_ts/tv\tnovel_tv\tts\tts/tv\ttv\n";
    print OUT2 "$prefix\t$trans_novel[0]\t$novel_ts_table\t$trans_novel[1]\t$trans[0]\t$ts_table\t$trans[1]\n\n";
}


close IN;

my $rscript=<<EOF;

library("ggplot2")
library(Cairo)

ColorMatrix<-c("#a50026","#d73027","#f46d43","#fdae61","#fee08b","#ffffbf","#d9ef8b","#a6d96a","#66bd63","#1a9850","#006837");
dt<-data.frame(values=c($no_list),labels=c("Novel","In dbSNP"));
his <- ggplot(dt,aes(x=labels,y=values,fill=labels))+geom_bar(stat = "identity");
his <- his + xlab('') + ylab('') + labs(fill="Type") + theme_bw()+ggtitle("dbSNP rate $dbSNP_rate% \\n $prefix")
CairoPNG("$prefix\_$type\_in_dbSNP_statistics.png",width = 489, height = 704)
his
dev.off()

dt2<-data.frame(values=c($tn_list),labels=c("Novel Ts","Novel Tv"));
his2 <- ggplot(dt2,aes(x=labels,y=values,fill=labels))+geom_bar(stat = "identity");
his2 <- his2 + xlab('') + ylab('') + labs(fill="Type") + theme_bw()+ggtitle("Novel Ts/Tv $novel_ts \\n $prefix")
CairoPNG("$prefix\_$type\_novel_tstv.png",width = 489, height = 704)
his2
dev.off()


dt3<-data.frame(values=c($tstv_list),labels=c("Ts","Tv"));
his3 <- ggplot(dt3,aes(x=labels,y=values,fill=labels))+geom_bar(stat = "identity");
his3 <- his3 + xlab('') + ylab('') + labs(fill="Type") + theme_bw()+ggtitle("Ts/Tv $ts \\n $prefix")
CairoPNG("$prefix\_$type\_tstv.png",width = 489, height = 704)
his3
dev.off()
cat ("Rscript ${prefix} generate ${type} table finished\n")

EOF
open RP,">${prefix}_generate_${type}_table.R" or die "could not create ${prefix}_generate_SNP_InDel_table.R file\n";
print RP $rscript."\n";
close RP;
#`/usr/bin/Rscript generate_SNP_InDel_table.R`;
