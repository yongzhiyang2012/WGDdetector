#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;
use Cwd 'abs_path';

my ($cds,$pep,$cluster_file,$tmp,$outputdir,$muscle,$pal2nal,$fastme,$Rscript,$thread_num)=@ARGV;
die "perl $0 cds_file pep_file cluster_file tmp_dir outputdir muscle pal2nal fastme Rscript thread_num\n" if (! $thread_num);

`mkdir $tmp` if (! -e $tmp);
`mkdir $tmp/align` if (! -e "$tmp/align");
`mkdir $tmp/align_large` if (! -e "$tmp/align_large");

my $basedir=abs_path($0);
$basedir=~s/\/[^\/]+\/[^\/]+$//;
my $cluastermaxnum=50;
my %list=&read_cluster($cluster_file);
my %cds=&read_fasta($cds);
my %pep=&read_fasta($pep);

open (RUNKS,">$tmp/run.ks.sh")||die"$!";
open (COLLECTKS,">$tmp/collect.ks.sh")||die"$!";
open (R2,">$tmp/run.large.split.sh");
open (R3,">$tmp/run.large.split.merge_cmd.sh");
for my $k1 (sort keys %list){
    #print COLLECTKS "perl $basedir/bin/collect_ks.pl $outputdir $k1 ";
    my @gene=sort keys %{$list{$k1}};
    if (scalar(@gene) <= $cluastermaxnum){
        open (O1,">$tmp/align/$k1.input.cds.file");
        open (O2,">$tmp/align/$k1.input.pep.file");
        for my $geneid (@gene){
            print O1 ">$geneid\n$cds{$geneid}\n";
            print O2 ">$geneid\n$pep{$geneid}\n";
        }
        close O1;
        close O2;
        print RUNKS "$basedir/bin/calculate_ks.pl $thread_num $tmp/align/$k1.input.cds.file $tmp/align/$k1.input.pep.file $muscle $pal2nal\n";
        print COLLECTKS "$basedir/bin/collect_ks.pl $outputdir $k1 $tmp/align/$k1.input.cds.file.align.output.ks.gz\n";
    }else{
        `mkdir $tmp/align_large/$k1` if (! -e "$tmp/align_large/$k1");
        open (O1,">$tmp/align_large/$k1/input.cds.file");
        open (O2,">$tmp/align_large/$k1/input.pep.file");
        for my $geneid (@gene){
            print O1 ">$geneid\n$cds{$geneid}\n";
            print O2 ">$geneid\n$pep{$geneid}\n";
        }
        close O1;
        close O2;
        print R2 "$basedir/bin/split_large_gf.pl $tmp/align_large/$k1/input.pep.file $muscle $fastme $Rscript\n";
    }
}
close RUNKS;
close COLLECTKS;
close R2;
print R3 "$basedir/bin/phase.new_sub_GF_seq.pl $cds $pep $tmp $thread_num $muscle $pal2nal $outputdir\n";
close R3;

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
sub read_cluster{
    my ($tmp_in_file)=@_;
    my %r;
    open (F,"$tmp_in_file")||die"$!";
    while (<F>) {
        chomp;
        my @a=split(/\s+/,$_);
        $a[0]=~s/\:$//;
        for (my $i=1;$i<@a;$i++){
            $r{$a[0]}{$a[$i]}++;
        }
    }
    close F;
    return %r;
}
