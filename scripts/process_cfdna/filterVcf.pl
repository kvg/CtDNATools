#!/usr/bin/perl -w

use strict;

use Cwd;
use FindBin;
use File::Basename;
use Data::Dumper;
use lib "$FindBin::Bin/";

use ParseArgs;
use DM;

my %args = &getCommandArguments(
    'VCF' => undef,
    'OUT' => undef,
);

open(OUT, ">$args{'OUT'}");

open(VCF, $args{'VCF'});
while (my $line = <VCF>) {
    chomp($line);

    if ($line =~ /^#/) {
        print OUT "$line\n";
    } else {
        my @fields = split(/\s+/, $line);

        my ($bamDepth) = $fields[7] =~ /BAM_DEPTH=(\d+)/;

        my @gts;
        my $allow = 1;

        my @fracs;

        for (my $i = 9; $i < scalar(@fields); $i++) {
            if ($fields[$i] eq '.') {
                $allow = 0;
            } else {
                my ($GT, $AD, $DP, $GQ, $PL) = split(/:/, $fields[$i]);
                my @dps = split(/,/, $AD);

                if ($GT eq '1' && $dps[1] < 10) {
                    $allow = 0;
                } else {
                    push(@fracs, 100*$dps[1]/$bamDepth);
                }
            }
        }

        for (my $i = 1; $i < scalar(@fracs); $i++) {
            if ($fracs[$i] < $fracs[$i-1]) {
                $allow = 0;
            }
        }

        if ($allow) {
            print OUT "\t$line\n";
        }
    }
}
close(VCF);

close(OUT);
