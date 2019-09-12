#!/usr/bin/perl
use strict;
use warnings;

my ($in,$out)=@ARGV;
my %h;
open (F,"$in")||die"$!";
<F>;
my $num=0;
my $gene;
while (<F>) {
    chomp;
    s/^\s+//;
    my @a=split(/\s+/,$_);
    if (/^gene\d+/){
        $gene=shift @a;
        $num++;
        $h{$gene}{num}=$num;
    }
    for my $a (@a){
        $h{$gene}{list} .= "$a\t";
    }
}
close F;

open (O,">$out")||die"$!";
my @keys=sort{$h{$a}{num} <=> $h{$b}{num}} keys %h;
print O "\t",join("\t",@keys),"\n";
for my $key (@keys){
    print O "$key\t$h{$key}{list}\n";
}
close O;

