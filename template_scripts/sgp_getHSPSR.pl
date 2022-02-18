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
my $PVER = sprintf("v%vd",$^V);
my $DATE = localtime;
my ($USER,$HOST);
if (defined($ENV{USER})) {
    $USER = $ENV{USER};
} else {
    chomp($USER = `whoami`);
};
if (defined($ENV{HOSTNAME})) {
    $HOST = $ENV{HOSTNAME};
} else {
    chomp($HOST = `hostname`);
};
my $host = $HOST; ###
#
my $PROG = 'getHSPSR.pl';
my $PRGVER = '0.9alpha';
#
my $chr = shift;
my ($seq,$src,$ftr,$ori,$end,$sco,$str,$frm) = (0..7);
my %g = ( "+1" => 0,  "+2" => 1,  "+3" => 2,  "+" => 3,
          "-1" => 4,  "-2" => 5,  "-3" => 6,  "-" => 7,
          ".1" => 8,  ".2" => 9,  ".3" => 10, "." => 11,
                     "all" => 12,            "sum" => [ (0) x 13 ] );
my $HSPminLEN = 1; # $ori + $HSPminLEN - 1 ==> minimum length is 1 nucleotide
my $S_CUTOFF = 26;
my $SCF = 0; # substract to tblastx scores S_CUTOFF - SCF
my $DSC = $S_CUTOFF - $SCF;
my $SHSP  = 0;    # SHSP=6 # shrink hsp by $SHSP
my $WTBX  = 0.19; # weigth of tblastx score
#my $WTBX  = 0.05; # weigth of tblastx score
my $WTBXF = 0.30; # weigth of tblastx score
my $WTBXI = 0.20; # weigth of tblastx score
my $WTBXT = 0.30; # weigth of tblastx score
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
    $g{sum}[$g{"$l[$str]$l[$frm]"}]++;
}; # while
$g{sum}[$g{"+"}] = $g{sum}[$g{"+1"}] + $g{sum}[$g{"+2"}] + $g{sum}[$g{"+3"}];
$g{sum}[$g{"-"}] = $g{sum}[$g{"-1"}] + $g{sum}[$g{"-2"}] + $g{sum}[$g{"-3"}];
$g{sum}[$g{"."}] = $g{sum}[$g{".1"}] + $g{sum}[$g{".2"}] + $g{sum}[$g{".3"}];
$g{sum}[$g{"all"}] = $g{sum}[$g{"+"}] + $g{sum}[$g{"-"}];
print STDERR "# TOTAL ".$g{sum}[$g{"all"}]." HSP-SRs on $chr: ".
             $g{sum}[$g{"+"}]." forward, ".$g{sum}[$g{"-"}]." reverse, ".
             $g{sum}[$g{"."}]." without strand.\n";
foreach my $t (qw/ 1 2 3 /) {
   printf STDERR "#\t%s : %6s\t\|\t%s : %6s\t|\t%s : %6s\n",
          "+$t",$g{sum}[$g{"+$t"}],"-$t",$g{sum}[$g{"-$t"}],
          ".$t",$g{sum}[$g{".$t"}];
}; # foreach
# foreach my $t (qw/ +1 +2 +3 -1 -2 -3 /) {
#    printf STDERR "#\t%s : %s\n",$t,$g{sum}[$g{$t}];
# }; # foreach
#
exit(0);
