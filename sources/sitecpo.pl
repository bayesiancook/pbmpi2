use strict;

my $infile1 = shift;
my $infile2 = shift;
my $burnin = shift;

open(INFILE1, $infile1) or die "intput error\n";
open(INFILE2, $infile2) or die "intput error\n";

for (my $rep=0; $rep<$burnin-1; $rep++)	{
	my $line = <INFILE1>;
	$line = <INFILE2>;
}

my $line1 = <INFILE1>;
chomp $line1;
my @a1 = split('\t', $line1);
my $nsite = (@a1) - 1;

my $line2 = <INFILE2>;
chomp $line2;
my @a2 = split('\t', $line2);
my $nsite2 = (@a2) - 1;

if ($nsite2 != $nsite) {
    die "error: non matching number of sites:$nsite\t$nsite2\n";
}

my @sitef1;
my @sitef2;
my @sitecpo1;
my @sitecpo2;
my @sitemean;
for (my $i=0; $i<$nsite; $i++)	{
	$sitef1[$i] = $a1[$i+1];
	$sitef2[$i] = $a2[$i+1];
	$sitecpo1[$i] = 0;
	$sitecpo2[$i] = 0;
	$sitemean[$i] = 0;
}

my $size = 0;
my $nrep = 0;
while ((my $line1 = <INFILE1>) && (my $line2 = <INFILE2>))	{
	chomp $line1;
	my @a1 = split('\t', $line1);
	chomp $line2;
	my @a2 = split('\t', $line2);
	for (my $i=0; $i<$nsite; $i++)	{
        my $tmp = $a1[$i+1] - $a2[$i+1];
		$sitemean[$i] = $sitemean[$i] + $tmp;
        $sitecpo1[$i] += exp($sitef1[$i] - $a1[$i+1]);
        $sitecpo2[$i] += exp($sitef2[$i] - $a2[$i+1]);
	}
	$nrep++;
}


my $grandmean = 0;
my $grandvar = 0;
my $logcpo1 = 0;
my $logcpo2 = 0;
for (my $i=0; $i<$nsite; $i++)	{
	$sitemean[$i] /= $nrep;
	$grandmean += $sitemean[$i];
	$grandvar += $sitemean[$i]*$sitemean[$i];

    $sitecpo1[$i] /= $nrep;
    $logcpo1 += $sitef1[$i] - log($sitecpo1[$i]);

    $sitecpo2[$i] /= $nrep;
    $logcpo2 += $sitef2[$i] - log($sitecpo2[$i]);
}

$grandmean /= $nsite;
$grandvar /= $nsite;
$grandvar -= $grandmean*$grandmean;

my $totmeanlogl = $nsite * $grandmean;
my $totvarlogl = $nsite * $grandvar;
my $totstdev = sqrt($totvarlogl);

print "delta log likelihood $totmeanlogl +/- $totstdev\n";
print "\n";
my $dcpo = $logcpo1 - $logcpo2;
print "delta log cpo: $dcpo\n";
print "\n";
print "cpo1: $logcpo1\n";
print "cpo2: $logcpo2\n";

