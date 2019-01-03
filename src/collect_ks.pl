use strict;
use warnings;

my $outputdir=shift or die "perl $0 outputdir cluster_name input_ks.gz.1 input_ks.gz.2 input_ks.gz.3 ...\n";
my $cluster_name=shift or die "perl $0 outputdir cluster_name input_ks.gz.1 input_ks.gz.2 input_ks.gz.3 ...\n";
my @result=@ARGV;
die "perl $0 outputdir cluster_name input_ks.gz.1 input_ks.gz.2 input_ks.gz.3 ...\n" if scalar(@result)<1;

`mkdir -p $outputdir` if (! -e $outputdir);
`mkdir -p $outputdir/ks_martix/$cluster_name` if (! -e "$outputdir/ks_martix/$cluster_name");

my %ks;
my %id;
for my $result (@result){
    open (F,"zcat $result|")||die"$!";
    while (<F>) {
        chomp;
        my @a=split(/\s+/,$_);
        $ks{$a[0]}{$a[1]}=$a[2];
        $id{$a[0]}++;
        $id{$a[1]}++;
    }
    close F;
}

my @id=sort keys %id;
open (O,">$outputdir/ks_martix/$cluster_name/ks.dist")||die"$!";
print O "\t",join("\t",@id),"\n";
for my $k1 (@id){
    print O "$k1\t";
    for my $k2 (@id){
        my $out=10;
        if (exists $ks{$k1}{$k2}){
            $out=$ks{$k1}{$k2};
        }elsif (exists $ks{$k2}{$k1}){
            $out=$ks{$k2}{$k1};
        }elsif ($k1 eq $k2){
            $out=0;
        }
        print  O "$out\t";
    }
    print  O "\n";
}
close O;
