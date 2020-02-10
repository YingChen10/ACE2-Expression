use strict;
use warnings;
my %length;my %s;
my $in=shift;
open(IN,$in);
while(my $line=<IN>){
	chomp $line;
	my @line=split/\s+/,$line;
	my @l=split/\./,$line[0];
	$s{$l[0]}=$line[1] unless ($s{$l[0]});
	$length{$l[0]}=$line[2] unless ($length{$l[0]});
	if ($length{$l[0]} < $line[2]){$length{$l[0]}=$line[2];}
}
close IN;
foreach(sort keys %length){
	print "$_\t$s{$_}\t$length{$_}\n";
}
