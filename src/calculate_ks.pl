#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;
use Bio::Align::DNAStatistics;
use Cwd 'abs_path';
use Bio::AlignIO;
use threads;

my ($max_cpu,$cds,$pep,$muscle,$pal2nal)=@ARGV;
die "perl $0 max_cpu cds_file pep_file mafft_path pal2nal_path\n" if (! $pal2nal);

my $basedir=abs_path($0);
$basedir=~s/\/[^\/]+\/[^\/]+$//;

## align ##
#my $cmd_align="$mafft --quiet --thread $max_cpu --auto $pep > $pep.align";
my $cmd_align="$muscle -in $pep -out $pep.align";
system ("$cmd_align") == 0  or die "system $cmd_align failed: $?\n";
## pal2nal ##
my $cmd_pal2nal="$pal2nal $pep.align $cds -output fasta > $cds.align";
system ("$cmd_pal2nal") == 0  or die "system $cmd_pal2nal failed: $?\n";

my %seq;
my $cds_aln_obj=Bio::SeqIO->new(-format=>"fasta",-file=>"$cds.align");
while (my $seqobj=$cds_aln_obj->next_seq) {
    my $id=$seqobj->id;
    my $seq=$seqobj->seq;
    $seq{$id}=$seq;
}
my @spid=sort keys %seq;
open (O,"| gzip -c >$cds.align.output.ks.gz")||die"$!";
for (my $i=0;$i<@spid;$i++){
    for (my $j=$i+1;$j<@spid;$j++){
        my ($spid1,$spid2)=($spid[$i],$spid[$j]);
        my ($spseq1,$spseq2)=($seq{$spid1},$seq{$spid2});
        my ($tmptwoseqalign,$newseq1,$newseq2)=&two_seq_align($spseq1,$spseq2);
        next if $tmptwoseqalign < 90;
        my $Ds=`$basedir/bin/calculate_ks.single.pl \"$spid1\" \"$newseq1\" \"$spid2\" \"$newseq2\"`;
        $Ds=10 if $Ds !~ /^\d+/;
        chomp $Ds;
        print O "$spid1\t$spid2\t$Ds\n";
    }
}
close O;

sub two_seq_align{
    my ($tmpseq1,$tmpseq2)=@_;
    my ($ralign,$rseq1,$rseq2)=(0,"","");
    my @tmpseq1=split(//,$tmpseq1);
    my @tmpseq2=split(//,$tmpseq2);
    for (my $i=0;$i<@tmpseq2;$i++){
        if (($tmpseq1[$i] ne '-') && ($tmpseq2[$i] ne '-')){
            $rseq1 .= $tmpseq1[$i];
            $rseq2 .= $tmpseq2[$i];
            $ralign++;
        }
    }
    return ($ralign,$rseq1,$rseq2);
}
