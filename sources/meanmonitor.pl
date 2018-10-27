use strict;

my $infile = shift;
my $burnin = shift;

open (INFILE, $infile) or die "input error\n";

for (my $i=0; $i<$burnin; $i++) {
    my $line = <INFILE>;
}

my $mean = 0;
my $var = 0;
my $n = 0;
foreach my $line (<INFILE>) {
    chomp $line;
    if ($line =~ /^(\S+)\s/g)   {
        $mean += $1;
        $var += $1*$1;
        $n += 1;
    }
}
$mean /= $n;
$var /= $n;
$var -= $mean*$mean;
my $stdev = sqrt($var);
print "mean: $mean\t$stdev\n";

