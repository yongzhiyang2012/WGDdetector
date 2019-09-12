#!/usr/bin/perl
use strict;
use warnings;
use Cwd 'abs_path';

my ($pep,$mucle,$fastme,$Rscript)=@ARGV;
die "perl $0 imput_pep mucle protdist\n" if ! $Rscript;

my $cluastermaxnum=50;
my $basedir=abs_path($0);
$basedir=~s/split_large_gf.pl$//;
$pep=~/^(\S+)\/input.pep.file/ or die "$pep\n";
my $tmpdir=$1;
my $pwd=`pwd`;chomp $pwd;

my $align_cmd="$mucle -maxiters 1 -diags -sv -distance1 kbit20_3 -in $pep -out $pep.best.fas";
system("$align_cmd");
my $convert_cmd="$basedir/Fasta2Phylip.pl $pep.best.fas $pep.best.fas.phy";
system("$convert_cmd");
chdir("$tmpdir");
my $dist_cmd="$fastme -i input.pep.file.best.fas.phy -O outfile -p LG";
system("$dist_cmd");
my $convert_cmd2="$basedir/conver_dist.pl outfile input.pep.file.best.fas.phy.protdist";
system("$convert_cmd2");
chdir("$pwd");

## hclust
my $linenum=`wc -l $pep.best.fas.phy.protdist`;chomp $linenum;
$linenum=~s/^(\d+)\s+.*/$1/;
$linenum=$linenum-40;
open (R,">$pep.best.fas.phy.protdist.hclust.R");
print R "inputdata=read.table(\"$pep.best.fas.phy.protdist\")\n";
print R "datadist=as.dist(inputdata)\n";
print R "datahclust=hclust(datadist)\n\n";
for (my $i=2;$i<=$linenum;$i++){
    print R "out_id=cutree(datahclust,k=$i)\n";
    print R "write.table(t(labels(out_id)),file=\"$pep.best.fas.phy.protdist.hclust.out\",append = TRUE, row.names = FALSE, col.names = FALSE)\n";
    print R "write.table(t(out_id),file=\"$pep.best.fas.phy.protdist.hclust.out\",append = TRUE, row.names = FALSE, col.names = FALSE)\n";
}
close R;
my $cmd="$Rscript $pep.best.fas.phy.protdist.hclust.R";
system("$cmd");

## phase sub family
my %dist;
my $sub_cluster=0;
open (F,"$pep.best.fas.phy.protdist.hclust.out")||die"$!";
while (<F>) {
    my $l1=$_;chomp $l1;
    my $l2=<F>;chomp $l2;
    my @l1=split(/\s+/,$l1);
    my @l2=split(/\s+/,$l2);
    my %flt;
    for (my $i=0;$i<@l1;$i++){
        $flt{$l2[$i]}{$l1[$i]}++;
    }
    my $check=0;
    for my $k1 (sort keys %flt){
        my @k2=sort keys %{$flt{$k1}};
        if (scalar(@k2)<$cluastermaxnum){
            next if exists $dist{$k2[0]};
            $sub_cluster++;
            for my $k2 (@k2){
	$dist{$k2}=$sub_cluster;
            }
        }else{
            $check++;
        }
    }
    last if $check == 0;
}
close F;
open (O,">$pep.best.fas.phy.protdist.hclust.out.subGF.list");
for my $k (sort keys %dist){
    print O "$k\t$dist{$k}\n";
}
close O;

