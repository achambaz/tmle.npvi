#!/usr/bin/perl -w
use strict;

my %hash = ();
my $file = $ARGV[0];

open (my $fh, "<", $file) or die "Can't open the file $file: ";

my $count = 0;

while (my $line =<$fh>) {
    $count++;
    chomp ($line);
    my @decomp = split("\;", $line);
    # print "@decomp\n";
    foreach (@decomp) {
	# print("\t$_\n");
	push(@{$hash{$_}}, $count); 
    }
}

while (my($key, $val) = each(%hash))
{
    # print "$key => @{$val}\n";
    print "$key\t";
    print "@{$val}\n";
}
