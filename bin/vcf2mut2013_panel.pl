#!/usr/bin/perl -w
#IN1=mutation list
#IN2=SNP list
open IN1, "<$ARGV[0]";
open IN2, "<$ARGV[1]";
open OUT, ">$ARGV[2]";

my %h2=();

while (<IN1>){
	chomp;
	my $rl=$_;
	my @tmp = split "\t",$rl;
	$h2{$tmp[34]}=$rl;
}

while (<IN2>){
	chomp;
	@_ = split;
	my $temp = join(",",$_[0],$_[1],$_[3],$_[4]);
#	print "$temp\n";
	if($h2{$temp}){
		print OUT $h2{$temp}."\n";
	}
}

close IN1;
close IN2;
close OUT;
