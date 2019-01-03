#!/usr/bin/perl
use strict;
use warnings;

my ($pep,$tmp_dir,$outdir,$makeblastdb,$blastp,$blast2graphs,$mcl,$pahsemcl,$thread_num)=@ARGV;
die "perl $0 pep tmp_dir out_dir makeblastdb blastp blast2graphs mcl pahsemcl_script thread_num\n" if (! $thread_num);

my $cmd;
`mkdir $tmp_dir` if (! -e "$tmp_dir");
`ln -s $pep $tmp_dir/input.pep.fa` if (! -e "$tmp_dir/input.pep.fa");

$cmd="$makeblastdb -in $tmp_dir/input.pep.fa -dbtype prot";
system ("$cmd") == 0  or die "system $cmd failed: $?\n";

$cmd="$blastp -db $tmp_dir/input.pep.fa -query $tmp_dir/input.pep.fa -out  $tmp_dir/all-vs-all.blast.out -evalue 1e-10 -outfmt '7 std qlen slen' -num_threads $thread_num";
system ("$cmd") == 0  or die "system $cmd failed: $?\n";

$cmd="$blast2graphs $tmp_dir/all-vs-all.blast.out $tmp_dir/mcl.input";
system ("$cmd") == 0  or die "system $cmd failed: $?\n";

$cmd="$mcl $tmp_dir/mcl.input_nrm_dmls_bit.abc --abc -I 1.5 -o  $tmp_dir/mcl.output";
system ("$cmd") == 0  or die "system $cmd failed: $?\n";

$cmd="$pahsemcl $tmp_dir/mcl.output $outdir/all-protein-clustering.result";
system ("$cmd") == 0  or die "system $cmd failed: $?\n";

