#!/usr/bin/perl
use strict;
use warnings;

my ($in,$out)=@ARGV;
die "perl $0 input_mcl_raw_result output_phased_mcl_result\n" if ! $out;
my $num=0;
open (F,"$in")||die"no such file: $in\n";
open (O,">$out") || die "cannot creat file: $out\n";
while (<F>) {
    chomp;
    $num++;
    print O "cluster$num:\t$_\n";
}
close F;
close O;
