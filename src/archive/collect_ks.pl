#!/usr/bin/perl
use strict;
use warnings;

my $infile=shift or die "perl $0 calculate_ks.pl_output_file\n";
my %outks;
my ($outdir,$cluster);
open (F,"$infile")||die"$!";
while (<F>) {
    chomp;
    my @a=split(/\s+/,$_);
    next unless /^NEED_COLLECTED\:/;
    ($outdir,$cluster)=($a[1],$a[2]);
    $outks{$a[3]}{$a[4]}=$a[5];
    $outks{$a[4]}{$a[3]}=$a[5]; 
}
close F;

`mkdir $outdir/ks_martix/$cluster` if (! -e "$outdir/ks_martix/$cluster");
open (OUTDIST,">$outdir/ks_martix/$cluster/ks.dist")||die "$!";
my @key1=sort keys %outks;
print OUTDIST "\t",join("\t",@key1),"\n";
for my $key1 (@key1){
    print OUTDIST "$key1\t";
    for my $key2 (@key1){
        if ($key1 eq $key2){
            print OUTDIST "0\t";
        }else{
            if (exists $outks{$key1}{$key2}){
	print OUTDIST "$outks{$key1}{$key2}\t";
            }else{
	print OUTDIST "10\t";
            }
        }
    }
    print OUTDIST "\n";
}
close OUTDIST;
