#!/usr/bin/perl
use strict;
use warnings;
use Cwd 'abs_path';

my $basedir=abs_path($0);
$basedir=~s/\/[^\/]+$//;

my ($cds,$pep,$tmpdir,$thread_num,$muscle,$pal2nal,$outputdir)=@ARGV;
die "perl $0 cds pep tmp thread_num muscle pal2nal outputdir\n" if ! $outputdir;

my %cds=&read_fa($cds);
my %pep=&read_fa($pep);
my %list=&read_list($tmpdir);

open (SH1,">$tmpdir/run.large.split.merge_cmd.ks.sh");
open (SH2,">$tmpdir/run.large.split.merge_cmd.ks.collect.sh");

for my $k1 (sort keys %list){
    for my $k2 (sort keys %{$list{$k1}}){
        my $outcluster="$k1"."_$k2";
        my $outname="$tmpdir/align/$k1"."_$k2";
        open (O1,">$outname.input.cds.file")||die"$!";
        open (O2,">$outname.input.pep.file")||die"$!";
        for my $k3 (sort keys %{$list{$k1}{$k2}}){
            print O1 ">$k3\n$cds{$k3}\n";
            print O2 ">$k3\n$pep{$k3}\n";
        }
        close O1;
        close O2;
        print SH1 "$basedir/calculate_ks.pl $thread_num $outname.input.cds.file $outname.input.pep.file $muscle $pal2nal\n";
        print SH2 "$basedir/collect_ks.pl $outputdir $outcluster $tmpdir/align/$outcluster.input.cds.file.align.output.ks.gz\n";
    }
}
close SH1;
close SH2;


sub read_fa{
    use Bio::SeqIO;
    my ($tmpin)=@_;
    my %r;
    my $tmpfa=Bio::SeqIO->new(-format=>"fasta",-file=>"$tmpin");
    while (my $tmpseq=$tmpfa->next_seq) {
        my $id=$tmpseq->id;
        my $seq=$tmpseq->seq;
        $r{$id}=$seq;
    }
    return %r;
}
sub read_list{
    my ($tmpinputdir)=@_;
    my %r;
    my @tmpin=<$tmpinputdir/align_large/*/input.pep.file.best.fas.phy.protdist.hclust.out.subGF.list>;
    for my $tmpin (@tmpin){
        $tmpin=~/align_large\/([^\/]+)\// or die "$tmpin\n";
        my $cluster=$1;
        open (F,"$tmpin")||die"$!";
        while (<F>) {
            chomp;
            s/\"//g;
            my @a=split(/\s+/,$_);
            $r{$cluster}{$a[1]}{$a[0]}++;
        }
        close F;
    }
    return %r;
}
