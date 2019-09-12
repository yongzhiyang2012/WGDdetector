#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;
use Bio::AlignIO;

my ($infile,$outfile)=@ARGV;
die"perl $0 infastafile outphyfile\n" if (! $outfile);
my $in=Bio::AlignIO -> new (-file=>"$infile",-format=>"fasta");
my $out=Bio::AlignIO -> new(-file=>">$outfile",-format=>"phylip",-line_length=>"50",-idlength=>"15");
while (my $aln = $in->next_aln) {
    $out->write_aln($aln);
}
