#!/bin/perl -w
# This is perl, version 5.6.1 built for i386-linux
#
# $Id: humus.nw,v 1.53 2002/06/17 18:57:08 jabril Exp $
#
# Re-scoring SRs to produce HSP-SRs for SGP homology.
#
# USAGE: getHSPSR.pl "chrname" < SR_file.gff > HSP-SR_file.gff > stdout.report
#
use strict;
#
#
my $chr = shift;
my ($seq,$src,$ftr,$ori,$end,$sco,$str,$frm) = (0..7);

my $HSPminLEN = 1; # $ori + $HSPminLEN - 1 ==> minimum length is 1 nucleotide
my $S_CUTOFF = 26;
my $SCF = 0; # substract to tblastx scores S_CUTOFF - SCF
my $DSC = $S_CUTOFF - $SCF;
my $SHSP  = 0;    # SHSP=6 # shrink hsp by $SHSP
my $WTBX  = 0.19; # weigth of tblastx score
# my $WTBX  = 0.05; # weigth of tblastx score
# my $WTBXF = 0.30; # weigth of tblastx score first exon
# my $WTBXI = 0.20; # weigth of tblastx score internal exon
# my $WTBXT = 0.30; # weigth of tblastx score terminal exon
while (<STDIN>) {
    my @l;
    next if /^#/o;
    next if /^\s*$/o;
    chomp;
    @l = split /\s+/o, $_;
    next unless $l[$sco] > $S_CUTOFF;
    $l[$sco] = sprintf("%.6f",($l[$sco] - $DSC) * $WTBX);
    $l[$ori] += $SHSP;
    $l[$end] -= $SHSP;
    next if $l[$end] < ($l[$ori] + $HSPminLEN - 1);
    print STDOUT (join("\t", (@l)))."\n";
};
print STDERR "DONE!";

exit(0);
