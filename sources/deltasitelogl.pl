use strict;

my $infile1 = shift;
my $infile2 = shift;

open (INFILE1, $infile1) or die "intput error\n";
open (INFILE2, $infile2) or die "intput error\n";

my $nsite = 0;
my $mean = 0;
my $var = 0;

while ((my $line1 = <INFILE1>) && (my $line2 = <INFILE2>))  {

    chomp $line1;
    chomp $line2;
    my @a1 = split('\t',$line1);
    my @a2 = split('\t',$line2);

    my $l1 = $a1[1];
    my $l2 = $a2[1];

    my $delta = $l1 - $l2;
    $mean += $delta;
    $var += $delta*$delta;
    $nsite += 1;
}

$mean /= $nsite;
$var /= $nsite;
$var -= $mean*$mean;

my $grandmean = $nsite*$mean;
my $grandvar = $nsite*$var;
my $grandstdev = sqrt($grandvar);

my $z = $grandmean / $grandstdev;
print "mean delta logl:\t$grandmean\t$grandstdev\t$z\n";

