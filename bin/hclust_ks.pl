#!/usr/bin/perl
use strict;
use warnings;

my ($in,$iternum,$Rscript)=@ARGV;
die "perl $0 input_ks_dist iternum Rscript\n" if (! $Rscript);

my $linenum=`wc -l $in`;chomp $linenum;
$linenum=~s/^(\d+)\s+.*/$1/;
$linenum=$linenum-1;
my $start=1;
if ($iternum ne "all"){
    $linenum=$iternum;
    $start=$iternum;
}

open (R,">$in.R")||die"$!";
print R "inputdata=read.table(\"$in\")\n";
print R "datadist=as.dist(inputdata)\n";
print R "datahclust=hclust(datadist)\n\n";
for (my $i=$start;$i<=$linenum;$i++){
    print R "out_id=cutree(datahclust,k=$i)\n";
    print R "write.table(t(labels(out_id)),file=\"$in.hcluster\",append = TRUE, row.names = FALSE, col.names = FALSE)\n";
    print R "write.table(t(out_id),file=\"$in.hcluster\",append = TRUE, row.names = FALSE, col.names = FALSE)\n";
}
close R;

my $cmd="$Rscript $in.R";
system("$cmd");
