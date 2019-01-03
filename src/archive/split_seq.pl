#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;

my ($cds,$pep,$clsuter_file,$line_num,$tmp,$mafft,$pal2nal,$codeml,$ctl)=@ARGV;
die "perl $0 cds pep clsuter_file line_num tmp mafft pal2nal codeml ctl\n" if (! $ctl);

`mkdir $tmp` if (! -e $tmp);
`mkdir $tmp/align` if (! -e "$tmp/align");

my %cds=&read_fasta($cds);
my %pep=&read_fasta($pep);

my $line=`head -n $line_num $clsuter_file | tail -n 1`;
chomp $line;
my @line=split(/\s+/,$line);
$line[0]=~s/\:$//;
`mkdir $tmp/align/$line[0]` if (! -e "$tmp/align/$line[0]");
open (SH,">$tmp/align/$line[0].run.sh")||die"$!";
for (my $i=1;$i<@line;$i++){
    for (my $j=$i+1;$j<@line;$j++){
        `mkdir $tmp/align/$line[0]/$i-$j` if (! -e "$tmp/align/$line[0]/$i-$j");
        open (O,">$tmp/align/$line[0]/$i-$j/cds");
        print O ">a\n$cds{$line[$i]}\n>b\n$cds{$line[$j]}\n";
        close O;
        open (O,">$tmp/align/$line[0]/$i-$j/pep");
        print O ">a\n$pep{$line[$i]}\n>b\n$pep{$line[$j]}\n";
        close O;
        print SH "cd $tmp/align/$line[0]/$i-$j ; $mafft --auto pep > pep.best.fas ; $pal2nal pep.best.fas cds -nogap -nomismatch -output paml > cds.paml; $codeml $ctl ; cd ../../../../\n";
    }
}
close SH;

sub read_fasta{
    my %r;
    my ($in)=@_;
    my $fa=Bio::SeqIO->new(-format=>"fasta",-file=>"$in");
    while (my $seq=$fa->next_seq) {
        my $id=$seq->id;
        my $seq=$seq->seq;
        $r{$id}=$seq;
    }
    return %r;
}
