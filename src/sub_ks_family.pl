#!/usr/bin/perl
use strict;
use warnings;

my ($ksdist,$kscluster,$fltmaxks,$outdir_prefix)=@ARGV;
die "perl $0 input_ks_dist input_kscluster flt_max_ks outdir_prefix\n" if (! $outdir_prefix);

my %ks=&read_ks_dist($ksdist);
my %out;
open (F,"$kscluster")||die"$!";
while (<F>) {
    chomp;
    s/\"//g;
    my $l1=$_;
    my $l2=<F>;chomp $l2;
    my @label=split(/\s+/,$l1);
    my @group=split(/\s+/,$l2);
    my @max_group=sort{$b<=>$a} @group;
    my $max_group=$max_group[0];
    my %flt;
    for (my $i=0;$i<@label;$i++){
        $flt{$group[$i]}{$label[$i]}++;
    }
    for my $k1 (sort keys %flt){
        my @k2=sort keys %{$flt{$k1}};
        my @phasek2;
        $out{$k2[0]}="delete" if ((scalar(@k2) == 1) && (! exists $out{$k2[0]}));
        for my $k2 (@k2){
            last if exists $out{$k2};
            push @phasek2,$k2;
        }
        my $maxks=0;
        for (my $i=0;$i<@phasek2;$i++){
            for (my $j=$i+1;$j<@phasek2;$j++){
	$maxks=$ks{$phasek2[$i]}{$phasek2[$j]} if ($ks{$phasek2[$i]}{$phasek2[$j]} > $maxks);
            }
        }
        next if $maxks > $fltmaxks;
        for my $k2 (@phasek2){
            $out{$k2}="$max_group-$k1";
        }
    }
    last if scalar(keys %out) == scalar(@label);
}
close F;
my %convert;
for my $k (sort keys %out){
    my $v=$out{$k};
    next if $v eq 'delete';
    $convert{$v}{$k}++;
}

my @outk=sort keys %convert;
for (my $i=0;$i<@outk;$i++){
    my $j=$i+1;
    `mkdir -p $outdir_prefix.$j` if (! -e "$outdir_prefix.$j");
    open (O,">$outdir_prefix.$j/ks.dist");
    my @v=sort keys %{$convert{$outk[$i]}};
    print O "\t",join("\t",@v),"\n";
    for my $v (@v){
        print O "$v\t";
        for my $v2 (@v){
            print O "$ks{$v}{$v2}\t";
        }
        print O "\n";
    }
}

sub read_ks_dist{
    my ($sub_infile)=@_;
    my @sub_id;
    my %r;
    open (F,"$sub_infile")||die"$!";
    while (<F>) {
        chomp;
        my @a=split(/\s+/,$_);
        if (/^\s+/){
            @sub_id=@a;
            next;
        }
        for (my $i=1;$i<@sub_id;$i++){
            $r{$a[0]}{$sub_id[$i]}=$a[$i];
        }
    }
    close F;
    return %r;
}
