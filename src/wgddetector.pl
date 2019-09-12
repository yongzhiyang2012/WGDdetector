#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Cwd 'abs_path';
use Bio::SeqIO;

my ($input_cds,$input_pep,$output_dir,$tmp_dir,$thread_num,$clean,$run_cluster,$cluster_engine);
GetOptions(
           'input_cds=s' => \$input_cds,
           'input_pep=s' => \$input_pep,
           'output_dir=s' => \$output_dir,
           'tmp_dir=s' => \$tmp_dir,
           'thread_num=i' => \$thread_num,
           'clean=s' => \$clean,
           'run_cluster=s' => \$run_cluster,
           'cluster_engine=s' => \$cluster_engine,
          );
if ((! $input_cds) || (! $input_pep) || (! $output_dir) || (! $tmp_dir) || (! $thread_num) || (! $cluster_engine)){
    &print_help;
    exit;
}
$clean="yes" if (! $clean);
$run_cluster="yes" if (! $run_cluster);
$input_cds=abs_path($input_cds);
$input_pep=abs_path($input_pep);

`mkdir $output_dir` if (! -e "$output_dir");
`mkdir $tmp_dir` if (! -e "$tmp_dir");

my $basedir=abs_path($0);
$basedir=~s/\/[^\/]+\/[^\/]+$//;
my %par=&read_config("$basedir/config/software.config");

## phase raw input cds and pep
print `date`;
print "## phasing the raw input cds and pep\n";
my $phaseseq_cmd="$basedir/bin/phase.id.pl $input_cds $input_pep $tmp_dir";
system("$phaseseq_cmd");
$input_cds=abs_path("$tmp_dir/tmp_cds.fa");
$input_pep=abs_path("$tmp_dir/tmp_pep.fa");

## running clustering ##
print `date`;
`mkdir $output_dir/01.cluster` if (! -e "$output_dir/01.cluster");
if ($run_cluster eq "yes"){
    print "## running protein clustering ##\n\n";
    if ($cluster_engine eq 'blastp'){
        my $run_cluster_cmd = "$basedir/bin/protein_blastp_cluster.pl $input_pep $tmp_dir/01.running_cluster $output_dir/01.cluster $par{makeblastdb} $par{blastp} $par{blast2graphs} $par{mcl} $basedir/bin/phase.mcloutp2orthomcl.format.pl $thread_num";
        system("$run_cluster_cmd");
    }elsif ($cluster_engine eq 'mmseqs2'){
        my $run_cluster_cmd = "$basedir/bin/protein_mmseqs2_cluster.pl $input_pep $tmp_dir/01.running_cluster $output_dir/01.cluster $par{mmseqs2} $par{blast2graphs} $par{mcl} $basedir/bin/phase.mcloutp2orthomcl.format.pl $thread_num";
        system("$run_cluster_cmd");
    }else{
        die "give the right --cluster_engine <blastp|mmseqs2>\n";
    }
}else{
    print "## skip running protein clustering ##\n\n";
}
`rm -r $tmp_dir/01.running_cluster` if (($clean eq 'yes') && ($run_cluster eq "yes"));

## running ks estamite ##
print `date`;
`mkdir $output_dir/02.ks_estimate` if (! -e "$output_dir/02.ks_estimate");
`mkdir $tmp_dir/02.ks_estimate` if (! -e "$tmp_dir/02.ks_estimate");
print "## running ks estimating ##\n\n";
if (! -e "$output_dir/01.cluster/all-protein-clustering.result"){
    die "no such file: $output_dir/01.cluster/all-protein-clustering.result\n";
}

my $run_split_seq_cmd="$basedir/bin/split_seq.pl $input_cds $input_pep $output_dir/01.cluster/all-protein-clustering.result $tmp_dir/02.ks_estimate $output_dir/02.ks_estimate $par{muscle} $par{pal2nal} $par{fastme} $par{Rscript} $thread_num";
system("$run_split_seq_cmd");
my $parallel_run_ks_estimate_cmd="$par{parallel} -j $thread_num < $tmp_dir/02.ks_estimate/run.ks.sh";
system("$parallel_run_ks_estimate_cmd");
my $parallel_collect_ks_cmd="$par{parallel} -j $thread_num < $tmp_dir/02.ks_estimate/collect.ks.sh";
system("$parallel_collect_ks_cmd");
if (-e "$tmp_dir/02.ks_estimate/run.large.split.sh"){
    print "## running large gene families\n";
    my $third_thread_num=int($thread_num/3); $third_thread_num=1 if $third_thread_num == 0;
    my $large_gf_split="$par{parallel} -j $third_thread_num < $tmp_dir/02.ks_estimate/run.large.split.sh";
    system("$large_gf_split");
    system("sh $tmp_dir/02.ks_estimate/run.large.split.merge_cmd.sh");
    system("$par{parallel} -j $thread_num < $tmp_dir/02.ks_estimate/run.large.split.merge_cmd.ks.sh");
    system("$par{parallel} -j $thread_num < $tmp_dir/02.ks_estimate/run.large.split.merge_cmd.ks.collect.sh");
}
`rm -r $tmp_dir/02.ks_estimate` if $clean eq 'yes';

## running Hierarchical Clustering ##
print `date`;
print "## running Hierarchical Clustering ##\n\n";
`mkdir $output_dir/03.hierarchial_clustering` if (! -e "$output_dir/03.hierarchial_clustering");
`mkdir $tmp_dir/03.hierarchial_clustering` if (! -e "$tmp_dir/03.hierarchial_clustering");
open (SHHC,">$tmp_dir/03.hierarchial_clustering/01.hcluster.all.sh");
open (SUBGF,">$tmp_dir/03.hierarchial_clustering/02.sub_GF.sh");
my @outRAWksdist=`ls $output_dir/02.ks_estimate/ks_martix`;
for my $outRAWksdist (@outRAWksdist){
    chomp $outRAWksdist;
    print SHHC "$basedir/bin/hclust_ks.pl $output_dir/02.ks_estimate/ks_martix/$outRAWksdist/ks.dist all $par{Rscript}\n";
    print SUBGF "$basedir/bin/sub_ks_family.pl $output_dir/02.ks_estimate/ks_martix/$outRAWksdist/ks.dist $output_dir/02.ks_estimate/ks_martix/$outRAWksdist/ks.dist.hcluster 5 $output_dir/03.hierarchial_clustering/$outRAWksdist\n";
}
close SHHC;
close SUBGF;
my $run_Hierarchical_Clustering_cmd="$par{parallel} -j $thread_num < $tmp_dir/03.hierarchial_clustering/01.hcluster.all.sh";
system("$run_Hierarchical_Clustering_cmd");
$run_Hierarchical_Clustering_cmd="$par{parallel} -j $thread_num < $tmp_dir/03.hierarchial_clustering/02.sub_GF.sh";
system("$run_Hierarchical_Clustering_cmd");
open (TWOCLAD,">$tmp_dir/03.hierarchial_clustering/03.two_clad.sh");
my @phased_gf=<$output_dir/03.hierarchial_clustering/*/ks.dist>;
for my $phased_gf (@phased_gf){
    print TWOCLAD "$basedir/bin/hclust_ks.pl $phased_gf 2 $par{Rscript}\n";
}
close TWOCLAD;
$run_Hierarchical_Clustering_cmd="$par{parallel} -j $thread_num < $tmp_dir/03.hierarchial_clustering/03.two_clad.sh";
system("$run_Hierarchical_Clustering_cmd");
`rm -r $tmp_dir/03.hierarchial_clustering` if  $clean eq 'yes';

## collecting final ks result for plot and comapring ##
print `date`;
print "## collecting final ks result for plot and comapring ##\n\n";
`mkdir $output_dir/04.final_paralogs_ks` if (! -e "$output_dir/04.final_paralogs_ks");
system("$basedir/bin/merge_final_ks.pl $output_dir/03.hierarchial_clustering $output_dir/04.final_paralogs_ks");

## detecting ks distribution ##
#print `date`;
#print "## detecting ks distribution ##\n\n";


sub read_config{
    my ($subin)=@_;
    die "no such file: $subin\n" if (! -e "$subin");
    my %r;
    open (SUBF,"$subin");
    while (<SUBF>) {
        chomp;
        next if /^#/;
        /^(\S+)\s+=\s+(\S+)/ or die "wrong format at $subin: $_\n";
        $r{$1}=$2;
    }
    close SUBF;
    return %r;
}
sub print_help{
    print STDERR<<EOF;

wgddetector (v1.00)

Usage: wgddetector --input_cds <cds_file> --input_pep <pep_file> --output_dir <output_dir> --tmp_dir <tmp_dir> --thread_num <thread_num> --cluster_engine <blastp|mmseqs2>

Options:
        --input_cds      input CDS in fasta format
        --input_pep      input protein in fasta format
        --output_dir     the output dir, which containing the main results
        --tmp_dir        the tmp dir
        --thread_num     the thread number when runnig this script
        --cluster_engine "blastp" or "mmseqs2". when setted this as mmseqs2, the protein aligning and clustering will be faster than blast
        --clean          "yes" or "no". If selected yes, the tmp_dir will removed after finish. Default: yes
        --run_cluster    "yes" or "no". If selected yes, the protein blast and mcl clustering will generate. If selected no, the clustering file "output_dir/01.cluster/all-protein-clustering.result" must be existed

  ** The sequence IDs within the CDS and protein files must be same!
EOF
}
