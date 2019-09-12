#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;

my ($incds,$inpep,$tmp)=@ARGV;
die "perl $0 input_cds input_pep tmp_dir\n" if ! $tmp;

my %cds=&read_fasta("$incds");
my %pep=&read_fasta("$inpep");

open (O1,">$tmp/tmp_cds.fa");
open (O2,">$tmp/tmp_pep.fa");
open (O3,">$tmp/tmp_id.table");
my $num=0;
for my $k1 (sort keys %cds){
    next if ! exists $pep{$k1};
    $num++;
    print O1 ">gene$num\n$cds{$k1}\n";
    print O2 ">gene$num\n$pep{$k1}\n";
    print O3 "$k1\tgene$num\n";
}
close O1;
close O2;
close O3;

sub read_fasta{
    my ($tmp_in_file)=@_;
    my $tmp_fa=Bio::SeqIO->new(-format=>"fasta",-file=>"$tmp_in_file");
    my %r;
    while (my $tmp_seq=$tmp_fa->next_seq) {
        my $tmp_id=$tmp_seq -> id;
        my $tmp_seq=$tmp_seq -> seq;
        $r{$tmp_id}=$tmp_seq;
    }
    return %r;
}
