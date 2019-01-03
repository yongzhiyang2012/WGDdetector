#!/usr/bin/perl
use strict;
use warnings;

my ($pep,$tmp_dir,$outdir,$mmseqs2,$blast2graphs,$mcl,$pahsemcl,$thread_num)=@ARGV;
die "perl $0 pep tmp_dir out_dir mmseqs2 blast2graphs mcl pahse_mcl_script thread_num\n" if (! $thread_num);

`mkdir $tmp_dir` if (! -e "$tmp_dir");
`ln -s $pep $tmp_dir/input.pep.fa` if (! -e "$tmp_dir/input.pep.fa");
`mkdir $tmp_dir/tmp` if (! -e "$tmp_dir/tmp");

my $cmd;
$cmd="$mmseqs2 createdb $tmp_dir/input.pep.fa $tmp_dir/selfDB";
system ("$cmd") == 0  or die "system $cmd failed: $?\n";

$cmd="$mmseqs2 search $tmp_dir/selfDB $tmp_dir/selfDB $tmp_dir/resultDB $tmp_dir/tmp --threads $thread_num";
system ("$cmd") == 0  or die "system $cmd failed: $?\n";

$cmd="$mmseqs2 convertalis $tmp_dir/selfDB $tmp_dir/selfDB $tmp_dir/resultDB $tmp_dir/all-vs-all.blast.out.raw --format-mode 2 --threads $thread_num";
system ("$cmd") == 0  or die "system $cmd failed: $?\n";

open (F,"$tmp_dir/all-vs-all.blast.out.raw")||die"$!";
open (O,">$tmp_dir/all-vs-all.blast.out")||die"$!";
while (<F>) {
    chomp;
    my @a=split(/\s+/,$_);
    next if $a[10]>1e-10;
    print O "$_\n";
}
close F;
close O;

$cmd="$blast2graphs $tmp_dir/all-vs-all.blast.out $tmp_dir/mcl.input";
system ("$cmd") == 0  or die "system $cmd failed: $?\n";

$cmd="$mcl $tmp_dir/mcl.input_nrm_dmls_bit.abc --abc -I 1.5 -o  $tmp_dir/mcl.output";
system ("$cmd") == 0  or die "system $cmd failed: $?\n";

$cmd="$pahsemcl $tmp_dir/mcl.output $outdir/all-protein-clustering.result";
system ("$cmd") == 0  or die "system $cmd failed: $?\n";
