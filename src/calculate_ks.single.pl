#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;
use Bio::Align::DNAStatistics;
use Cwd 'abs_path';
use Bio::AlignIO;

my ($spid1,$newseq1,$spid2,$newseq2)=@ARGV;
die "perl spid1 spseq1 spid2 spseq2\n" if (! $newseq2);

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

print "$Ds\n";
