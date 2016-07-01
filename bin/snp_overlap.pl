#!usr/bin/perl -w
#stat the overlap of two SNP dataset
################################################################
#
# perl $0 novogene_SNP MH_SNP novogene_InDel MH_InDel sample_name
#
################################################################

open IN1,"$ARGV[0]" or die "could not open novogene_SNP file\n";
open IN2,"$ARGV[1]" or die "could not open MH_SNP file\n";
open IN3,"$ARGV[2]" or die "could not open novogene_InDel file\n";
open IN4,"$ARGV[3]" or die "could not open MH_InDel file\n";
$name = $ARGV[4];
open OUT1,">${name}_overlap.snp.vcf";
open OUT2,">${name}_overlap.indel.vcf";
open OUT3,">${name}_stat.txt";
######unique part###
open OUT4,">${name}_novo.unique.snp.vcf";
open OUT5,">${name}_MH.unique.snp.vcf";
open OUT6,">${name}_novo.unique.indel.vcf";
open OUT7,">${name}_MH.unique.indel.vcf";

my @novo_snp;
my @mh_snp;
my $overlap_snp=0;
my $novo_num_snp=0;
my $mh_num_snp=0;
my $snp_ratio=0;
my @novo_indel;
my @mh_indel;
my $overlap_indel=0;
my $novo_num_indel=0;
my $mh_num_indel=0;
my $indel_ratio=0;
my @overlap_snp;
my @overlap_indel;

my %h_novo_snp=();
my $h_mh_snp=();
my %h_novo_indel=();
my $h_mh_indel=();

#######################
while (<IN1>){
	if ($_ !~ /^#/){
		chomp;
        my $temp_novo_snp=$_;
		my @tmp1_snp = split(/\t/,$_);
		my $combine_snp = join(",",$tmp1_snp[0],$tmp1_snp[1],$tmp1_snp[3],$tmp1_snp[4]);
		my $novo_pos_snp = join("","chr","$combine_snp");
		push @novo_snp,$novo_pos_snp;
        $h_novo_snp{$novo_pos_snp}=$temp_novo_snp;
	}
}
while (<IN2>){
	if ($_ !~ /^#/){
        chomp;
        my $temp_mh_snp=$_;
        $mh_num_snp++;
        my @tmp2_snp = split(/\t/,$_);
        my $mh_pos_snp = join(",",$tmp2_snp[0],$tmp2_snp[1],$tmp2_snp[3],$tmp2_snp[4]);
        push @mh_snp,$mh_pos_snp;
        $h_mh_snp{$mh_pos_snp}=$temp_mh_snp;
        foreach (@novo_snp){
            if ($_ eq $mh_pos_snp){
                push @overlap_snp,$mh_pos_snp;###
                print OUT1 join"\t",@tmp2_snp;
                print OUT1 "\n";
                $overlap_snp++;
            }
           
		}
	}
}
foreach(@novo_snp){
    chomp;
    my $in=$_;
    unless ($h_mh_snp{$in}){
        print OUT4 $h_novo_snp{$in}."\n";
    }
}

foreach(@mh_snp){
    chomp;
    my $in=$_;
    unless($h_novo_snp{$in}){
        print OUT5 $h_mh_snp{$in}."\n";
    }
}
$novo_num_snp=@novo_snp;
$snp_ratio=$overlap_snp/$novo_num_snp*100;

##############find out unique SNPs seperately##############
while (<IN3>){
    if ($_ !~ /^#/){
        chomp;
    my $temp_novo_indel=$_;
        my @tmp1_indel = split(/\t/,$_);
        my $combine_indel = join(",",$tmp1_indel[0],$tmp1_indel[1],$tmp1_indel[3],$tmp1_indel[4]);
        my $novo_pos_indel = join("","chr","$combine_indel");
        push @novo_indel,$novo_pos_indel;
        $h_novo_indel{$novo_pos_indel}=$temp_novo_indel;
    }
}

while (<IN4>){
    if ($_ !~ /^#/){
        chomp;
        my $temp_mh_indel=$_;
        $mh_num_indel++;
        my @tmp2_indel = split(/\t/,$_);
        my $mh_pos_indel = join(",",$tmp2_indel[0],$tmp2_indel[1],$tmp2_indel[3],$tmp2_indel[4]);
        push @mh_indel,$mh_pos_indel;
        $h_mh_indel{$mh_pos_indel}=$temp_mh_indel;
        foreach (@novo_indel){
            if ($_ eq $mh_pos_indel){
                print OUT2 join"\t",@tmp2_indel;
                print OUT2 "\n";
                $overlap_indel++;
            }
        }
    }
}

foreach(@novo_indel){
    chomp;
    my $in=$_;
    unless ($h_mh_indel{$in}){
        print OUT6 $h_novo_indel{$in}."\n";
    }
}

foreach(@mh_indel){
    chomp;
    my $in=$_;
    unless($h_novo_indel{$in}){
        print OUT7 $h_mh_indel{$in}."\n";
    }
}
$novo_num_indel=@novo_indel;
$indel_ratio=$overlap_indel/$novo_num_indel*100;

##########################################################



print OUT3 "Sample\tMH_SNP\tNovo_SNP\tOverlap_SNP\tRatio(%)\tMH_InDel\tNovo_InDel\tOverlap_InDel\tRatio(%)\n";
print OUT3 "$name\t";
print OUT3 "$mh_num_snp\t";
print OUT3 "$novo_num_snp\t";
print OUT3 "$overlap_snp\t";
printf OUT3 "%0.2f\t",$snp_ratio;
print OUT3 "$mh_num_indel\t";
print OUT3 "$novo_num_indel\t";
print OUT3 "$overlap_indel\t";
printf OUT3 "%0.2f\t",$indel_ratio;











	
