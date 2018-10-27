use strict;

my $infile = shift;
my $burnin = shift;
my $index = shift;

open (INFILE, $infile) or die "input error\n";

for (my $i=0; $i<$burnin; $i++) {
    my $line = <INFILE>;
}

my $mean = 0;
my $var = 0;
my $n = 0;
foreach my $line (<INFILE>) {
    chomp $line;
    my @a = split('\t',$line);
    my $tmp = $a[$index-1];
    $mean += $tmp;
    $var += $tmp*$tmp;
    $n += 1;
}
$mean /= $n;
$var /= $n;
$var -= $mean*$mean;
my $stdev = sqrt($var);
print "mean: $mean\t$stdev\n";

