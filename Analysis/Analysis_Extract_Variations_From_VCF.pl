#!/usr/bin/perl -w
use strict;
use warnings;
#use lib 'd:/e-books';
#use beginperlbioinfo;

#The final file column names
#Variation: Variation ID
#Pos: The positions of variations in X chromosome
#Ref: Reference alleles
#Alt: Alternate alleles
#EAS: The frequency of alternate alleles in East Asians
#AF: The frequency of alternate alleles in all populations
print "variation\tPos\tRef\tAlt\tEAS\tAF\n";

#Read in ACE2 variation data from NCBI (https://www.ncbi.nlm.nih.gov/variation/tools/1000genomes/)
my $file="1000G_chrX_15575052-15624296.popvcf";
open(IN, $file) or die ("can not open $file\n");
while (my $line=<IN>) {
	chomp $line;
		my @mid=split /\t/,$line;
		#Only analysis ACE2 variantions rows 
	if ($line!~/\#/ && $mid[2]!~/;/ && $mid[2]!~/\./ 
		&& $mid[1]<=15619137  && $mid[1]>=15579156) {
		my $EAS=0;
		my $AF=0;
		my @cdx=split /:/,$mid[13];#cdx: Allele count in Chinese Dai in Xishuangbanna 
		my @chb=split /:/,$mid[15];#chb: Allele count in Han Chinese in Beijing
		my @chs=split /:/,$mid[16];#chs: Allele count in Southern Han Chinese
		my @jpt=split /:/,$mid[25];#jpt: Allele count in Japanese in Tokyo
		my @khv=split /:/,$mid[26];#khv: Allele count in Kinh in Ho Chi Minh City
		my $Asian_T=0;#Allele count in East Asians
		my $Asian_A=0;#Alternate allele count in East Asians
		$Asian_T=$cdx[0]+$chb[0]+$chs[0]+$jpt[0]+$khv[0];
		$Asian_A=$cdx[1]+$chb[1]+$chs[1]+$jpt[1]+$khv[1];

		my $total=0;#Allele count in all populations
		my $alt=0;#Alternate allele count in all populations
		my @global=split /:/,$mid[9];
		for (my $i=10 ;$i<=$#mid;$i++) {
			my @mid2=split /:/,$mid[$i];
			$total+=$mid2[0];
			$alt+=$mid2[1];
		}
		# If there is no alternate alleles in EAS, then the count of alternate alleles is assigned to 0.5
		if ($Asian_A ==0) {
			$Asian_A =0.5;
		}
		$EAS=$Asian_A/$Asian_T;
		$AF=$alt/$total;
		print "$mid[2]\tX:$mid[1]\t$mid[3]\t$mid[4]\t$EAS\t$AF\n";
	}
}
close IN;


