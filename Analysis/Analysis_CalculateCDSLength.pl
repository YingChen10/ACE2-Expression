#calculate the length of each transcript
$in=shift;
open IN ,$in or die;
$/=">";
while (<IN>)
{
	chomp;
	next unless $_=~/\w+/;
	my @a=(split /\n/,$_,2);
	$a[1]=~s/\s+//g;
	my $n=(split /\s+/,$a[0])[3];
	my $s=(split /\s+/,$a[0])[6];
	my $l=length $a[1];
	$n=~s/gene://g;
	$s=~s/gene_symbol://g;
	print "$n\t$s\t$l\n";

}
$/="\n";
