use strict;

my $infile1 = shift;
my $infile2 = shift;
my $burnin = shift;

open(INFILE1, $infile1) or die "intput error\n";
open(INFILE2, $infile2) or die "intput error\n";

my $line = <INFILE1>;
chomp $line;

my @a = split('\t', $line);
my $nsite = (@a) - 1;

my @sitemean;
for (my $i=0; $i<$nsite; $i++)	{
	$sitemean[$i] = 0;
}

$line = <INFILE2>;

for (my $rep=1; $rep<$burnin; $rep++)	{
	my $line = <INFILE1>;
	$line = <INFILE2>;
}

my $size = 0;
my $nrep = 0;
while ((my $line1 = <INFILE1>) && (my $line2 = <INFILE2>))	{
	chomp $line1;
	my @a1 = split('\t', $line1);
	chomp $line2;
	my @a2 = split('\t', $line2);
	for (my $i=0; $i<$nsite; $i++)	{
		$sitemean[$i] = $sitemean[$i] + ($a1[$i+1] - $a2[$i+1]);
	}
	$nrep++;
}


my $grandmean = 0;
my $grandvar = 0;
for (my $i=0; $i<$nsite; $i++)	{
	$sitemean[$i] /= $nrep;
	$grandmean += $sitemean[$i];
	$grandvar += $sitemean[$i]*$sitemean[$i];
}
$grandmean /= $nsite;
$grandvar /= $nsite;
$grandvar -= $grandmean*$grandmean;

my $totmeanlogl = $nsite * $grandmean;
my $totvarlogl = $nsite * $grandvar;
my $totstdev = sqrt($totvarlogl);

print "delta log likelihood $totmeanlogl +/- $totstdev\n";

