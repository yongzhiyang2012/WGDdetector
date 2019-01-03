use strict;
use warnings;
use Cwd 'abs_path';

my $base=abs_path($0);
$base=~s/\/setup.pl//;

`mkdir $base/bin` if (! -e "$base/bin");

my $perl_path=`which perl`; chomp $perl_path;

my @pl=`ls $base/src`;
for my $pl (@pl){
    chomp $pl;
    next unless $pl=~/pl$/;
    open (O,">$base/bin/$pl");
    open (F,"$base/src/$pl");
    my $line=0;
    while (<F>) {
        chomp;
        $line++;
        if ($line == 1){
            print O "#!$perl_path\n";
        }
        else {
            print O "$_\n";
        }
    }
    close F;
    close O;
    `chmod +x $base/bin/$pl`;
}
