#!/usr/bin/perl
use strict;
use warnings;
use threads;
use threads::shared;
use Bio::SeqIO;
use Bio::Align::DNAStatistics;
use Cwd 'abs_path';
use Bio::AlignIO;
use Bio::Align::Utilities qw(aa_to_dna_aln);

my ($max_cpu,$cds,$pep,$cluster_file,$tmp,$outputdir,$mafft)=@ARGV;
die "perl $0 max_cpu cds_file pep_file cluster_file tmp_dir output_dir mafft_path\n" if (! $mafft);

`mkdir $tmp` if (! -e $tmp);
`mkdir $tmp/align` if (! -e "$tmp/align");
`mkdir $outputdir/ks_martix` if (! -e "$outputdir/ks_martix");

my $basedir=abs_path($0);
$basedir=~s/\/[^\/]+\/[^\/]+$//;

my %list=&read_cluster($cluster_file);
my %cds=&read_fasta($cds);
my %pep=&read_fasta($pep);

my %allks :shared;

for my $k1 (sort keys %list){
    
    `mkdir $outputdir/ks_martix/$k1` if (! -e "$outputdir/ks_martix/$k1");
    my @seqid=sort keys %{$list{$k1}};
    for (my $i=0;$i<@seqid;$i++){
        my @cmd;
        for (my $j=$i+1;$j<@seqid;$j++){
            push @cmd,"G$i";
            push @cmd,"G$j";
            push @cmd,$seqid[$i];
            push @cmd,$seqid[$j];
        }
        my $jobs=0;
        my $finish_jobs=0;
        while (1) {
            while (scalar(threads->list()) < $max_cpu) {
        	last if $jobs >= scalar(@cmd)/4;
        	$jobs++;
          	my ($outGid1,$outGid2,$needid1,$needid2)=($cmd[(($jobs-1)*4)],$cmd[(($jobs-1)*4)+1],$cmd[(($jobs-1)*4)+2],$cmd[(($jobs-1)*4)+3]);
	$allks{$k1} = &share({});
	$allks{$k1}{$outGid1} = &share({});
	#$allks{$k1}{$outGid1}{$outGid2} = &share({});
	#$allks{$k1}{$outGid2} = &share({});
	#$allks{$k1}{$outGid2}{$outGid1} = &share({});
	threads->create(\&run_Nei_Gojobori_ks,$k1,$outGid1,$outGid2,$needid1,$needid2,$cds{$needid1},$pep{$needid1},$cds{$needid2},$pep{$needid2});
            }
            for my $single_job (threads->list(threads::all)){
        	if ($single_job -> is_joinable()){
        	    $single_job -> join();
        	    $finish_jobs++;
        	}
            }
            last if ($finish_jobs == (scalar(@cmd)/4));
        }
    }
}

for my $allksk1 (sort keys %allks){
    open (O,">$outputdir/ks_martix/$allksk1/ks.dist")||die"$!";
    for my $k2 (sort keys %{$allks{$allksk1}}){
        for my $k3 (sort keys %{$allks{$allksk1}{$k2}}){
            print O "$allksk1\t$k2\t$k3\t$allks{$allksk1}{$k2}{$k3}\n";
        }
    }
    close O;
    print "$outputdir/ks_martix/$allksk1/ks.dist\n";
    exit;
}


 
=cut
my $genenumlen=scalar(keys $allks{$allksk1});
    my @genenum=(0..$genenumlen);
    my @geneid;
    for my $genenum (@genenum){ push @geneid,"G$genenum";}
    print O "\t",join("\t",@geneid),"\n";
    for (my $outi=0;$outi<@geneid;$outi++){
        my $geneid1=$geneid[$outi];
        print O "$geneid1\t";
        for (my $outj=0;$outj<@geneid;$outj++){
            my $geneid2=$geneid[$outj];
            if ($geneid1 eq $geneid2){
	print O "0\t";
            }else{
	if ($outi<$outj){
	    print O "$allks{$allksk1}{$geneid1}{$geneid2}\t";
	}elsif ($outj > $outi){
	    print O "$allks{$allksk1}{$geneid2}{$geneid1}\t";
	}
            }
        }
        print O "\n";
    }
    close O;
}
=cut

sub run_Nei_Gojobori_ks{
    my ($tmp_clusterid,$tmp_outid1,$tmp_outid2,$tmp_id1,$tmp_id2,$tmp_cds_seq1,$tmp_pep_seq1,$tmp_cds_seq2,$tmp_pep_seq2)=@_;
    my %tmp_cds;
    my $tmp_read_cds_filehandle;
    open ($tmp_read_cds_filehandle,"echo \">$tmp_id1\n$tmp_cds_seq1\n>$tmp_id2\n$tmp_cds_seq2\n\"|");
    my $tmp_read_cds_obj=Bio::SeqIO->new(-fh=>$tmp_read_cds_filehandle,-format=>"fasta");
    my $tmp_read_cds_obj_single_1=$tmp_read_cds_obj->next_seq;
    $tmp_cds{$tmp_read_cds_obj_single_1 -> display_id}=$tmp_read_cds_obj_single_1;
    my $tmp_read_cds_obj_single_2=$tmp_read_cds_obj->next_seq;
    $tmp_cds{$tmp_read_cds_obj_single_2 -> display_id}=$tmp_read_cds_obj_single_2;
    close $tmp_read_cds_filehandle;
    my $tmp_mafft_result=`echo \">$tmp_id1\n$tmp_pep_seq1\n>$tmp_id2\n$tmp_pep_seq2\" | $mafft --quiet -`;
    my $tmp_pep_filehandle;
    open ($tmp_pep_filehandle,"echo \"$tmp_mafft_result\"|");
    my $tmp_pep_align=Bio::AlignIO->new(-fh=>$tmp_pep_filehandle,-format=>"fasta");
    my $tmp_pep_align_obj=$tmp_pep_align->next_aln();
    my $tmp_dna_aln = aa_to_dna_aln($tmp_pep_align_obj,\%tmp_cds);
    close ($tmp_pep_filehandle);
    my $tmp_stats = Bio::Align::DNAStatistics->new();
    my $result = $tmp_stats->calc_all_KaKs_pairs($tmp_dna_aln);
    my ($Da, $Ds, $Dn, $N, $S, $S_d, $N_d);
    for my $an (@$result){
        for my $result_key (sort keys %$an ){
            next if $result_key=~/Seq/;
            if($result_key eq "D_n"){$Dn = $an->{$result_key}};
            if($result_key eq "D_s"){$Ds = $an->{$result_key}};
            if($result_key eq "S_d"){$S_d = $an->{$result_key};}
            if($result_key eq "N_d"){$N_d = $an->{$result_key};}
            if($result_key eq "S"){$S = $an->{$result_key};}
            if($result_key eq "N"){$N = $an->{$result_key};}
        }
    }
    if($Dn !~ /\d/){$Dn = 10;}
    if($Ds !~ /\d/){$Ds = 10;}
    $allks{$tmp_clusterid}{$tmp_outid1}{$tmp_outid2}=$Ds;
    lock($allks{$tmp_clusterid}{$tmp_outid1}{$tmp_outid2});
#$allks{$tmp_clusterid}{$tmp_outid2}{$tmp_outid1}=$Ds;
}

sub randstr{
    my $maxLenth=16;
    my @a = (0..9,'a'..'z','A'..'Z');
    my $password = join '', map { $a[int rand @a] } 0..($maxLenth-1);
    return $password;
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
