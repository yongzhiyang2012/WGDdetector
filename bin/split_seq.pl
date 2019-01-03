#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;
use Cwd 'abs_path';

my ($cds,$pep,$cluster_file,$tmp)=@ARGV;
die "perl $0 cds_file pep_file cluster_file tmp_dir\n" if (! $tmp);

`mkdir $tmp` if (! -e $tmp);
`mkdir $tmp/align` if (! -e "$tmp/align");

my $basedir=abs_path($0);
$basedir=~s/\/[^\/]+\/[^\/]+$//;

my %list=&read_cluster($cluster_file);
my %cds=&read_fasta($cds);
my %pep=&read_fasta($pep);


for my $k1 (sort keys %list){
    open (O1,">$tmp/align/$k1.input.cds.file")||die"$!";
    open (O2,">$tmp/align/$k1.input.pep.file")||die"$!";
    for my $seqid (sort keys %{$list{$k1}}){
        print O1 ">$seqid\n$cds{$seqid}\n";
        print O2 ">$seqid\n$pep{$seqid}\n";
    }
    close O1;
    close O2;
}

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
