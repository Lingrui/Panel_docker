use strict;
use warnings;

###################################
#
# Usage: perl $0 sample_interval_summary sample_name
#
###################################

my $path="$ARGV[0]";
my $prefix="$ARGV[1]";

my $rscript=<<EOF;

library("ggplot2")
library(Cairo)
data<-read.table("$path",header=T)

depth<-data[,3]
coverage_1X<-data[,9]
coverage_4X<-data[,10]
coverage_20X<-data[,11]

##1X coverage of capture region
first <- ggplot(data,aes(x=coverage_1X,y=..count../sum(..count..)*100))+
  geom_histogram(binwidth =5,fill = "turquoise4",bin=10)+
  scale_y_continuous(breaks = seq(0, 100,10),limits = c(0, 100), expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(0,100,10),limits =c(0,110),expand=c(0,0))+
  ylab("Percentage(%)")+
  xlab("1X Coverage(%)")+
  theme_bw()+
  ggtitle("1X in capture region \\n $prefix")
CairoPNG("$prefix.coverage_1X.png", width = 326, height = 704)
first
dev.off()

##4X coverage of capture region
second <- ggplot(data,aes(x=coverage_4X,y=..count../sum(..count..)*100))+
  geom_histogram(binwidth =5,fill = "turquoise4",bin=10)+
  scale_y_continuous(breaks = seq(0, 100,10),limits = c(0, 100), expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(0,100,10),limits =c(0,110),expand=c(0,0))+
  ylab("Percentage(%)")+
  xlab("4X Coverage(%)")+
  theme_bw()+
  ggtitle("4X in capture region \\n $prefix")
CairoPNG("$prefix.coverage_4X.png", width = 326, height = 704)
second
dev.off()

##20X coverage of capture region
third <- ggplot(data,aes(x=coverage_20X,y=..count../sum(..count..)*100))+
  geom_histogram(binwidth =5,fill = "turquoise4",bin=10)+
  scale_y_continuous(breaks = seq(0, 100,10),limits = c(0, 100), expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(0,100,10),limits =c(0,110),expand=c(0,0))+
  ylab("Percentage(%)")+
  xlab("20X Coverage(%)")+
  theme_bw()+
  ggtitle("20X in capture region \\n $prefix")
CairoPNG("$prefix.coverage_20X.png", width = 326, height = 704)
third
dev.off()

fourth <-ggplot(data, aes(x = depth,y=..density..)) +
  scale_x_continuous(breaks = seq(0, 300,50),limits = c(0, 300), expand = c(0, 0)) +
  geom_histogram(binwidth = 2,fill = "turquoise1",bin=2,colour='turquoise4')+
#  geom_line(stat="density",adjust=2, colour="violetred4")+
  theme_bw()+
  xlab("Depth in target region(X)")+
  ggtitle("Sequencing depth of capture region \\n $prefix")
CairoPNG("$prefix.capture_depth.png", width = 977, height = 704)
fourth
dev.off()

cat ("Rscript ${prefix} capture draw finished\n")

EOF
open RP,">${prefix}_capture_draw.R" or die "could not create ${prefix}_capture_draw.R file\n";
print RP $rscript."\n";
close RP;
#`/usr/bin/Rscript capture_draw.R`;

