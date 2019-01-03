#!/usr/bin/perl
use strict;
use warnings;

my ($indir,$outdir)=@ARGV;
die "perl $0 inputdir outputdir\n" if (! $outdir);

my @in=<$indir/*/ks.dist>;
open (O,">$outdir/final.ks.distribution.list");
for my $in (@in){
    $in=~/(\w+\.\d+)\/ks.dist/;
    my $clustername=$1;
    my %ks=&read_ks_dist($in);
    open (F,"$in.hcluster")||die"$!";
    my $lable=<F>;chomp $lable;
    my $group=<F>;chomp $group;
    close F;
    $lable=~s/\"//g;
    my @lable=split(/\s+/,$lable);
    my @group=split(/\s+/,$group);
    my ($allks,$allnum)=(0,0);
    my %group;
    for (my $i=0;$i<@lable;$i++){
        $group{$group[$i]}{$lable[$i]}++;
    }
    for my $k1 (sort keys %{$group{1}}){
        for my $k2 (sort keys %{$group{2}}){
            my $singleks=$ks{$k1}{$k2};
            $allks += $singleks;
            $allnum ++ ;
        }
    }
    die "$in\n" if $allnum == 0;
    print O "$clustername\t",$allks/$allnum,"\n";
}
close O;


    
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
