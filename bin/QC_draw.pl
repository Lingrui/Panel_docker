### input1 is file name without suffix fastq or fq   eg.asdf_1 for asdf_1.fq
### automatic find mate asdf_2
### input2 is the path of Fastqc output with extract file
###################################
#
# Usage: perl $0 asdf_1 fastqc_output
#
###################################
use File::Basename;
my @suffixlist = qw(.fastq .fq .tar.gz);
#my ($name,$path) = fileparse($ARGV[0]);
$name = $ARGV[0];
#$dir = "/data/anzhen/workspace/output/QC";
$dir = $ARGV[1];
#print "$name  aaa  $path\n";
$file2 = $name;
#$file2 =~ s/_R1/_R2/;
#$file2 =~ s/_1/_2/;
$file2 =~ s/(.*)_1/\1_2/;
$file1 = "$dir/$name.fq_fastqc/fastqc_data.txt";
$file2 = "$dir/$file2.fq_fastqc/fastqc_data.txt"; 
#$name =~ s/_R1|_1//;
$name =~ s/(.*)_1/\1/;
#print $name, "aaaaaa\n";
@lines =();
#print "$file1 aa $file2 \n";
if ((-e $file2) && !($file1 eq $file2)){ # if exists mate 
	open F1, $file1; #
	$col_name = "";$save =0; $last_base=0;
	@read_lines = <F1>;
	for ($i=0;$i<@read_lines;$i++){
		if ($read_lines[$i] =~ /tile\s+sequence\s+quality/){
			$last_line = $read_lines[$i-2];
			@w =split(/\s+/,$last_line);
#print $last_line."\n";
			$last_base = int($w[0]);
#print $last_base ,"aaaadfasdf\n";
			last;
		}
		if ($save ==1){
#if (/^\d+-\d+/){ #  18-19 write into same two lines
			if ($read_lines[$i]=~/^\d+-\d+/){ #  18-19 write into same two lines
#			@w = split(/\s+/,$read_lines[$i]);
#print "$read_lines[$i] a a a \n";
#				@w = split(/-/,$w[0]);
				$tmp1 = $read_lines[$i];
				$tmp1 =~ s/^(\d+)\-\d+/$1/;
				push(@lines,$tmp1);
				$tmp2 = $read_lines[$i];
				$tmp2 =~ s/^(\d+)\-(\d+)/$2/;
				push(@lines,$tmp2);

			}	
			else{push(@lines,$read_lines[$i]);}
		}
		if ($read_lines[$i] =~ /Quartile/){
			$col_name = $read_lines[$i];$save=1;
			$col_name =~ s/#//;
		}
	}
	unshift (@lines,$col_name);
	pop(@lines);
############################################## file 2 process
	open F2, $file2; #
	$col_name = "";$save =0; 
	@read_lines = <F2>;
	for ($i=0;$i<@read_lines;$i++){
		if ($read_lines[$i] =~ /tile\s+sequence\s+quality/){
#$last_line = $read_lines[$i-2];
#			@w =split(/\s+/,$last_line);
#print $last_line."\n";
#			$last_base = $w[0];
			last;
		}
		if ($save ==1){
#if (/^\d+-\d+/){ #  18-19 write into same two lines
			if ($read_lines[$i]=~/^\d+-\d+/){ #  18-19 write into same two lines
#			@w = split(/\s+/,$read_lines[$i]);
#print "$read_lines[$i] a a a \n";
#				@w = split(/-/,$w[0]);
				$tmp1 = $read_lines[$i];
				$tmp1 =~ /^(\d+)-(\d+)/ ;
				$digital = $1+$last_base;
				$tmp1 =~ s/^(\d+)-(\d+)/$digital/;
				push(@lines,$tmp1);
				$tmp2 = $read_lines[$i];
				$tmp2 =~ /^(\d+)-(\d+)/;
				$digital = $2+$last_base;
				$tmp2 =~ s/^(\d+)-(\d+)/$digital/;
				push(@lines,$tmp2);

			}	
			else{	
				$tmp1 = $read_lines[$i];
				$tmp1 =~ /^(\d+)/ ;
				$digital = $1+$last_base;
				$tmp1 =~ s/^(\d+)/$digital/;
					push(@lines,$tmp1);
			}
		}
		if ($read_lines[$i] =~ /Quartile/){
			$col_name = $read_lines[$i];$save=1;
			$col_name =~ s/#//;
		}
	}
	pop(@lines);
}
#print "$last_base aaa \n";
open TMP,">${name}_tmp.txt";
print TMP "@lines" ;#> "tmp.txt";
close TMP;
my $rscript=<<EOF;

library(ggplot2)
colour_from_novogene=rgb(144,205,180,max=255)
df = read.table("${name}_tmp.txt",header=T,fill=T,as.is=T,sep="\t")
df = df[1:2]
df\$P = 10^(-df\$Mean/10)
png("$name.error_rate.png",width = 977, height = 704)
ggplot(data=df,aes(x=Base,y=P))+
    geom_bar(colour=colour_from_novogene,stat="identity")+
	ylab("% Error rate") + xlab("Position along reads") +
    theme_bw()+
	ggtitle("Error rate distribution along reads of $name") +
    guides(fill=F)
dev.off()
png("$name.quality_score.png",width = 977, height = 704)
ggplot(data=df,aes(x=Base,y=Mean,fill=Base))+
    scale_y_continuous(limits=c(0,45))+
	ylab("Quality Score") + xlab("Position along reads") +
    theme_bw()+
	ggtitle("Quality along reads of $name") +
    geom_point(col=colour_from_novogene,size=3)+
    geom_line(col=colour_from_novogene,stat="identity")+
    guides(fill=F)

cat ("Rscript $name QC draw finished\n")

dev.off()
EOF
open RP,">${name}_QC_draw.R" or die "could not create ${name}_QC_draw.R file\n";
print RP $rscript."\n";
close RP; 
#`Rscript QC_draw.R`;
#`rm -rf tmp.txt`;
#means.barplot <- qplot(x=group, y=mean, fill=variable,
#                       data=means, geom="bar", stat="identity",
#                      position="dodge")
