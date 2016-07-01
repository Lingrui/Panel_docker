#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Cwd;

#######panel pipeline via docker#########

#basic command--------------------------------------------------
my $docker="docker run -v /usr/bin/docker:/usr/bin/docker -v /var/run/docker.sock:/var/run/docker.sock -v /data:/data -v $ENV{'RAWDATA_PATH'}:/rawdata -v $ENV{'TMP_PATH'}:/tmp -v $ENV{'LOG_PATH'}:/log -v $ENV{'REPORT_PATH'}:/report -v $ENV{'GTF_PATH'}:/gtf -v $ENV{'BED_PATH'}:/bed -v $ENV{'REF_PATH'}:/ref ";

#project information--------------------------------------------
my $reference = $ENV{"GENOME"};
my $Intervals = $ENV{"BED"};
my $samples= $ENV{"DATASET"};
my $setname= $ENV{"SETNAME"};
my $specific = $ENV{"SPECIFIC"};
my $Thread = $ENV{"THREAD"};
my $seqtype = $ENV{"ENDPOINT"};

die "Could not find reference file\n" if($reference eq "" or $Intervals eq "" or $samples eq "");
$Thread ||=2;
my $ip = "-ip 100";
my $SN = 0;

#Working directory-------------------------------------------------------------
my @sample = split ",",$samples;
my $outdir = "/report/output";
my $jsonout = "/report";
&createDir("$outdir/SNP");
&createDir("$outdir/fig_table");
&createDir("$outdir/tmp");

my $WORKSPACE_DIR = "$jsonout";
my $OUTPUT_DIR = "$outdir";
my $SNP = "$outdir/SNP";
my $fig_table = "$outdir/fig_table";
my $Tmp = "$outdir/tmp";

#Software--------------------------------------------------------------------
my $samtools="192.168.6.28:5000/weixuan/samtools:v1";
my $gatk="192.168.6.28:5000/weixuan/gatk:v3.5 java -jar /opt/GenomeAnalysisTK.jar";
my $picardtools="192.168.6.28:5000/weixuan/gatk:v3.5 java -jar /opt/picard-tools-2.0.1/picard.jar";
my $bin="/weixuan/panel/bin";
my $annovar="/weixuan/panel/annovar";

#generate runshell file----------------------------------------------------------
open RUN,">$outdir/runShell.sh" or die "could not create runShell.sh\n";

#prepare for reference genome and database-----------------------------------------
my $REF_Human = "/ref/build37/build37.fa";
my $RefDir = "/ref/build37";
my $humandb = "/ref/humandb";

my $Mills = "Mills_and_1000G_gold_standard.indels.b37plus.vcf";
my $G1000 = "1000G_phase1.indels.b37plus.vcf";
my $dbsnp = "dbsnp_138.b37plus.vcf";

foreach my $one_sample (@sample){
        my $filename_pre=basename($one_sample);
        $filename_pre=~s/.bam//ig;
   	my $PREFIX=$filename_pre;
   	open OUT,">$outdir/$PREFIX.sh" or die "could not create $PREFIX.sh\n"; ######
	my $outsh=<<EOF;

#Prepare for SNP calling-------------------------------------------------------
echo "Process 0." >>${OUTPUT_DIR}/${PREFIX}.log
$docker $samtools sort -m5G $one_sample $Tmp/$PREFIX.sorted 2>> ${OUTPUT_DIR}/${PREFIX}.log
$docker $samtools index $Tmp/$PREFIX.sorted.bam

#Test if sorted BAM exist-----------------------------------------------------
if [ ! -e $Tmp/$PREFIX.sorted.bam ];then
echo "Process ERROR error occured in sorting BAM file.">>${OUTPUT_DIR}/${PREFIX}.log
echo "Process ERROR error occured in sorting BAM file.">>${OUTPUT_DIR}/${PREFIX}.error

#SNP Calling-------------------------------------------------------------------
else echo "Start SNP calling for $PREFIX at \`date\` ......" >>${OUTPUT_DIR}/${PREFIX}.log
echo "Process 30." >>${OUTPUT_DIR}/${PREFIX}.log

# add RG and sort
$docker $picardtools AddOrReplaceReadGroups I=$Tmp/$PREFIX.sorted.bam \\
O=$Tmp/$PREFIX.sorted.head.bam SO=coordinate ID=$PREFIX PL=IONTORRENT \\
LB=LB PU=PU SM=$PREFIX VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true TMP_DIR=$Tmp/ &>> ${OUTPUT_DIR}/${PREFIX}.log

#Run Indel Realignment
echo "Start Indel Realignment for $PREFIX at \`date\` ......" >>${OUTPUT_DIR}/${PREFIX}.log
echo "Process 39." >>${OUTPUT_DIR}/${PREFIX}.log

$docker $gatk \\
-T RealignerTargetCreator -R $REF_Human -I $Tmp/$PREFIX.sorted.head.bam \\
-o $Tmp/$PREFIX.realigner.intervals -rf BadCigar -known $RefDir/$G1000 \\
-known $RefDir/$Mills -U ALLOW_SEQ_DICT_INCOMPATIBILITY -L $Intervals $ip &>>${OUTPUT_DIR}/${PREFIX}.log

$docker $gatk \\
-T IndelRealigner -R $REF_Human -I $Tmp/$PREFIX.sorted.head.bam \\
-targetIntervals $Tmp/$PREFIX.realigner.intervals \\
-o $Tmp/$PREFIX.realigned.bam -U ALLOW_SEQ_DICT_INCOMPATIBILITY \\
-model USE_READS -rf BadCigar -known $RefDir/$G1000 -known $RefDir/$Mills &>>${OUTPUT_DIR}/${PREFIX}.log
echo "Indel realignment finished at \`date\`." >>${OUTPUT_DIR}/${PREFIX}.log
echo "Process 47." >>${OUTPUT_DIR}/${PREFIX}.log

##Run Base Recalibration
$docker $gatk \\
-T BaseRecalibrator -I $Tmp/$PREFIX.realigned.bam \\
-R $REF_Human -nct $Thread \\
-o $Tmp/$PREFIX.recal_data.txt \\
-rf BadCigar -U ALLOW_SEQ_DICT_INCOMPATIBILITY -L $Intervals $ip \\
-knownSites $RefDir/$G1000 \\
-knownSites $RefDir/$Mills \\
-knownSites $RefDir/$dbsnp 2>>${OUTPUT_DIR}/${PREFIX}.log

$docker $gatk \\
-T PrintReads -I $Tmp/$PREFIX.realigned.bam \\
-R $REF_Human -nct $Thread -L $Intervals $ip \\
-BQSR $Tmp/$PREFIX.recal_data.txt \\
-o $Tmp/$PREFIX.realigned.recal.bam 2>>${OUTPUT_DIR}/${PREFIX}.log
echo "Base Recalibration finished at \`date\`." >>${OUTPUT_DIR}/${PREFIX}.log

#Test if Base Recalibration output exist-----------------------------------------------------
if [ ! -e $Tmp/$PREFIX.realigned.recal.bam ];then
echo "Process ERROR error occured in Base Recalibration.">>${OUTPUT_DIR}/${PREFIX}.log
echo "Process ERROR error occured in Base Recalibration.">>${OUTPUT_DIR}/${PREFIX}.error

# Call variants
else $docker $gatk \\
-T HaplotypeCaller -R $REF_Human -I $Tmp/$PREFIX.realigned.recal.bam -L $Intervals $ip \\
--dbsnp $RefDir/$dbsnp -o $Tmp/$PREFIX.raw.snps.indels.HaplotypeCaller.vcf \\
-stand_call_conf 30 -stand_emit_conf 30 -minPruning 3 2>>${OUTPUT_DIR}/${PREFIX}.log
echo "Call variants finished at \`date\`." >>${OUTPUT_DIR}/${PREFIX}.log
echo "Process 59." >>${OUTPUT_DIR}/${PREFIX}.log

##Filter to a call set
echo "Start filter SNP for $PREFIX at \`date\` ......" >>${OUTPUT_DIR}/${PREFIX}.log
$docker $gatk \\
-T SelectVariants -R $REF_Human -V $Tmp/$PREFIX.raw.snps.indels.HaplotypeCaller.vcf \\
-selectType SNP \\
-sn ${PREFIX} \\
-o $Tmp/$PREFIX.raw_snps.vcf 2>>${OUTPUT_DIR}/${PREFIX}.log

$docker $gatk \\
-T SelectVariants -R $REF_Human -V $Tmp/$PREFIX.raw.snps.indels.HaplotypeCaller.vcf \\
-selectType INDEL \\
-sn ${PREFIX} \\
-o $Tmp/$PREFIX.raw_indels.vcf 2>>${OUTPUT_DIR}/${PREFIX}.log

echo "Finished Variant Analysis at \`date\`...... " >>${OUTPUT_DIR}/${PREFIX}.log
echo "Process 76." >>${OUTPUT_DIR}/${PREFIX}.log

cat $Tmp/$PREFIX.raw_snps.vcf $Tmp/$PREFIX.raw_indels.vcf > $SNP/$PREFIX.SNP_InDel.vcf

#Test if variants calling output exist-----------------------------------------------------
#if [ ! -e $SNP/$PREFIX.SNP_InDel.vcf ];then
#echo "Process ERROR error occured in variants calling.">>${OUTPUT_DIR}/${PREFIX}.log
#echo "Process ERROR error occured in variants calling.">>${OUTPUT_DIR}/${PREFIX}.error

##Annotation
cd $SNP
echo "Start SNP annotation for $PREFIX at \`date\` ......" >>${OUTPUT_DIR}/${PREFIX}.log

perl $annovar/table_annovar.pl \\
$Tmp/$PREFIX.raw_snps.vcf $humandb -buildver hg19 \\
-out $PREFIX.SNP -remove \\
-protocol refGene,esp6500siv2_all,1000g2015aug_all,1000g2015aug_eas,snp138,ljb26_all \\
-operation g,f,f,f,f,f -nastring . -vcfinput &>>${OUTPUT_DIR}/${PREFIX}.log&

perl $annovar/table_annovar.pl \\
$Tmp/$PREFIX.raw_indels.vcf $humandb -buildver hg19 \\
-out $PREFIX.InDel -remove \\
-protocol refGene,esp6500siv2_all,1000g2015aug_all,1000g2015aug_eas,snp138,ljb26_all \\
-operation g,f,f,f,f,f -nastring . -vcfinput &>>${OUTPUT_DIR}/${PREFIX}.log&
wait

cat $SNP/$PREFIX.SNP.hg19_multianno.txt $SNP/$PREFIX.InDel.hg19_multianno.txt > $SNP/$PREFIX.SNP_InDel.hg19_multianno.txt
#perl $bin/vcf2mut.pl $RefDir/panel.pos.txt $SNP/$PREFIX.SNP_InDel.vcf $SNP/$PREFIX.candidates.txt 2>>${OUTPUT_DIR}/${PREFIX}.log
#perl $bin/report.pl $SNP/$PREFIX.candidates.txt $SNP/$PREFIX.SNP_InDel.hg19_multianno.txt $WORKSPACE_DIR/$PREFIX.output.json $RefDir/pubmed.info
perl $bin/vcf2mut_pathogenic.pl $RefDir/aneurysma_patho.txt $SNP/$PREFIX.SNP_InDel.vcf $SNP/$PREFIX.candidates.txt 2>>${OUTPUT_DIR}/${PREFIX}.log
perl $bin/report_pathogenic.pl $SNP/$PREFIX.candidates.txt $SNP/$PREFIX.SNP_InDel.hg19_multianno.txt $WORKSPACE_DIR/$PREFIX.output.json $RefDir/pubmed.info

#Test if annotation output exist-----------------------------------------------------
if [ ! -e $WORKSPACE_DIR/$PREFIX.output.json ];then
echo "Process ERROR error occurred in annotation.">>${OUTPUT_DIR}/${PREFIX}.log
echo "Process ERROR error occurred in annotation.">>${OUTPUT_DIR}/${PREFIX}.error

else echo "Finish SNP annotation for $PREFIX at \`date\` ......" >>${OUTPUT_DIR}/${PREFIX}.log
echo "Process 85." >>${OUTPUT_DIR}/${PREFIX}.log

###Figures and tables
cd $fig_table

#Depth of Coverage
$docker $gatk \\
-T DepthOfCoverage -L $Intervals -R $REF_Human -I $Tmp/$PREFIX.sorted.head.bam \\
-o ${fig_table}/$PREFIX -ct 1 -ct 4 -ct 20 \\
-omitBaseOutput 2>>${OUTPUT_DIR}/${PREFIX}.log
echo "Finish depth of coverage stat of $PREFIX at \`date\` ......" >>${OUTPUT_DIR}/${PREFIX}.log

perl $bin/capture_draw.pl ${PREFIX}.sample_interval_summary $PREFIX >>${OUTPUT_DIR}/${PREFIX}.log
echo "Finish capture stat of $PREFIX at \`date\` ......" >>${OUTPUT_DIR}/${PREFIX}.log
R --vanilla --slave <${PREFIX}_capture_draw.R &>>${OUTPUT_DIR}/${PREFIX}.log
echo "Finish capture draw of $PREFIX at \`date\` ......" >>${OUTPUT_DIR}/${PREFIX}.log
montage ${PREFIX}.coverage_1X.png ${PREFIX}.coverage_4X.png ${PREFIX}.coverage_20X.png \\
-tile 3x1 -geometry +0+0+0 ${PREFIX}.capture_coverage.png &>>${OUTPUT_DIR}/${PREFIX}.log

perl $bin/generate_SNP_InDel_table.pl $SNP/${PREFIX}.InDel.hg19_multianno.vcf $PREFIX InDel >>${OUTPUT_DIR}/${PREFIX}.log
perl $bin/generate_SNP_InDel_table.pl $SNP/${PREFIX}.SNP.hg19_multianno.vcf $PREFIX SNP >>${OUTPUT_DIR}/${PREFIX}.log
echo "Finish SNP and InDel table of $PREFIX at \`date\` ......" >>${OUTPUT_DIR}/${PREFIX}.log
R --vanilla --slave <${PREFIX}_generate_InDel_table.R &>>${OUTPUT_DIR}/${PREFIX}.log
R --vanilla --slave <${PREFIX}_generate_SNP_table.R &>>${OUTPUT_DIR}/${PREFIX}.log
echo "Finish SNP and InDel table draw of $PREFIX at \`date\` ......" >>${OUTPUT_DIR}/${PREFIX}.log

perl $bin/InDel_statistics.pl $SNP/${PREFIX}.InDel.hg19_multianno.txt \\
$PREFIX ${PREFIX}.InDel_statistics.txt >>${OUTPUT_DIR}/${PREFIX}.log
perl $bin/SNP_statistics.pl $SNP/${PREFIX}.SNP.hg19_multianno.txt \\
$PREFIX ${PREFIX}.SNP_statistics.txt >>${OUTPUT_DIR}/${PREFIX}.log
echo "Finish SNP and InDel stat of $PREFIX at \`date\` ......" >>${OUTPUT_DIR}/${PREFIX}.log

R --vanilla --slave <${PREFIX}_InDel_statistics.R &>>${OUTPUT_DIR}/${PREFIX}.log
R --vanilla --slave <${PREFIX}_SNP_statistics.R &>>${OUTPUT_DIR}/${PREFIX}.log
montage ${PREFIX}_SNP_genotype_statistics.png ${PREFIX}_SNP_in_dbSNP_statistics.png -tile 2x1 -geometry +0+0 ${PREFIX}_SNP_statistics.png >>${OUTPUT_DIR}/${PREFIX}.log
montage ${PREFIX}_InDel_genotype_statistics.png ${PREFIX}_InDel_in_dbSNP_statistics.png -tile 2x1 -geometry +0+0 ${PREFIX}_InDel_statistics.png >>${OUTPUT_DIR}/${PREFIX}.log
montage ${PREFIX}_SNP_tstv.png ${PREFIX}_SNP_novel_tstv.png -tile 2x1 -geometry +0+0 ${PREFIX}_tstv.png >>${OUTPUT_DIR}/${PREFIX}.log
echo "Finish SNP and InDel stat draw of $PREFIX at \`date\` ......" >>${OUTPUT_DIR}/${PREFIX}.log

echo "Process 100." >>${OUTPUT_DIR}/${PREFIX}.log
fi
fi
fi
EOF
		print OUT "$outsh\n";
		close OUT;
		print RUN "bash $outdir/$PREFIX.sh&\n";
		$SN++;
		if ($SN%15==0){
			print RUN "wait\n";
		}
}
print RUN "wait\n";
close RUN;
`bash $outdir/runShell.sh`;

chdir("$fig_table");
`bash $bin/generate_summary.sh`;

#make DIR-------------------------------------------------------------------
sub createDir{
    my ($path)=@_;
    if(! -e $path){
        `mkdir -p $path`;
    }
}


