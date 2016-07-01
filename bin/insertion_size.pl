use strict;
use warnings;

##################################
#
# Usage: perl $0 picard_insertionsize sample_name
#
##################################

open IN, "$ARGV[0]" or die "could not open insertionsize file.\n";
my $sample=$ARGV[1];
open OUT, ">${sample}_insertionsize.txt;
my $start=0;
while (<IN>){
    chomp;
    my $temp = $_;
    if ($start==1){
        if ($temp=~ /^\d+\td+/){
            print OUT $temp



