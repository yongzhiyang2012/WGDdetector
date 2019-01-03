#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;
use Bio::Align::DNAStatistics;
use Cwd 'abs_path';
use Bio::AlignIO;

my ($max_cpu,$cds,$pep,$cluster_name,$outputdir,$mafft,$pal2nal)=@ARGV;
die "perl $0 max_cpu cds_file pep_file cluster_name output_dir mafft_path pal2nal_path\n" if (! $pal2nal);

`mkdir $outputdir/ks_martix` if (! -e "$outputdir/ks_martix");
`mkdir $outputdir/ks_martix/$cluster_name` if (! -e "$outputdir/ks_martix/$cluster_name");

my $basedir=abs_path($0);
$basedir=~s/\/[^\/]+\/[^\/]+$//;

## align ##
my $cmd_align="$mafft --quiet --thread $max_cpu --auto $pep > $pep.align";
system ("$cmd_align") == 0  or die "system $cmd_align failed: $?\n";
## pal2nal ##
my $cmd_pal2nal="$pal2nal $pep.align $cds -output fasta > $cds.align";
system ("$cmd_pal2nal") == 0  or die "system $cmd_pal2nal failed: $?\n";

my %ks;
my %id;
my %seq;
my $cds_aln_obj=Bio::SeqIO->new(-format=>"fasta",-file=>"$cds.align");
while (my $seqobj=$cds_aln_obj->next_seq) {
    my $id=$seqobj->id;
    my $seq=$seqobj->seq;
    $seq{$id}=$seq;
}
my @spid=sort keys %seq;
for (my $i=0;$i<@spid;$i++){
    for (my $j=$i+1;$j<@spid;$j++){
        my ($spid1,$spid2)=($spid[$i],$spid[$j]);
        my ($spseq1,$spseq2)=($seq{$spid1},$seq{$spid2});
        my ($tmptwoseqalign,$newseq1,$newseq2)=&two_seq_align($spseq1,$spseq2);
        next if $tmptwoseqalign < 90;
        my $tmpfh;
        open ($tmpfh,"echo \">$spid1\n$newseq1\n>$spid2\n$newseq2\n\"|");
        my $cds_aln_obj=Bio::AlignIO->new(-format=>"fasta",-fh=>$tmpfh);
        my $cds_aln=$cds_aln_obj->next_aln();
        my $stats = Bio::Align::DNAStatistics->new(); 
        my $results = $stats->calc_all_KaKs_pairs($cds_aln);
        my $Ds;
        for my $an (@$results){
            $Ds= $an->{D_s};
        }
        $Ds=10 if $Ds !~ /^\d+/;
        $ks{$spid1}{$spid2}=$Ds;
        $id{$spid1}++;
        $id{$spid2}++;
        close $tmpfh;
    }
}

my @id=sort keys %id;
open (O,">$outputdir/ks_martix/$cluster_name/ks.dist")||die"$!";
print O "\t";
for (my $i=0;$i<@id;$i++){
    print O "G$i\t";
}
print O "\n";
for (my $i=0;$i<@id;$i++){
    print O "G$i\t";
    for (my $j=0;$j<@id;$j++){
        my $outks;
        if ($i == $j){
            $outks=0;
        }elsif (exists $ks{$id[$i]}{$id[$j]}){
            $outks=$ks{$id[$i]}{$id[$j]};
        }elsif (exists $ks{$id[$j]}{$id[$i]}) {
            $outks=$ks{$id[$j]}{$id[$i]};
        }else{
            $outks=10;
        }
        print O "$outks\t";
    }
    print O "\n";
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
