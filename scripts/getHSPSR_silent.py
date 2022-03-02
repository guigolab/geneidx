# This is perl, version 5.6.1 built for i386-linux
#
# $Id: humus.nw,v 1.53 2002/06/17 18:57:08 jabril Exp $
#
# Re-scoring SRs to produce HSP-SRs for SGP homology.
#
# USAGE: getHSPSR.pl "chrname" < SR_file.gff > HSP-SR_file.gff > stdout.report
#
#
#
from sys import argv
filename = argv[1]
out_filename = argv[2]
# chr = shift;
# my ($seq,$src,$ftr,$ori,$end,$sco,$str,$frm) = (0..7);

seq,src,ftr,ori,end,sco,str,frm = range(0,8)

HSPminLEN = 1; # $ori + $HSPminLEN - 1 ==> minimum length is 1 nucleotide
S_CUTOFF = 26;
SCF = 0; # substract to tblastx scores S_CUTOFF - SCF
DSC = S_CUTOFF - SCF;
SHSP  = 0;    # SHSP=6 # shrink hsp by $SHSP
WTBX  = 0.19; # weigth of tblastx score
# WTBX  = 0.05; # weigth of tblastx score
# WTBXF = 0.30; # weigth of tblastx score first exon
# WTBXI = 0.20; # weigth of tblastx score internal exon
# WTBXT = 0.30; # weigth of tblastx score terminal exon

lines_to_write = []
with open(filename, "r") as a_file:
    for line in a_file:
        if line.startswith('#'):
            continue

        # not sure how to translate the perl line below
        # and I can't find any weirder line
        # next if /^\s*$/o;
        # if line.startswith('#'):
        #     continue

        l = line.strip().split()

        if float(l[sco]) > S_CUTOFF:
            l[sco] = round((float(l[sco]) - DSC) * WTBX, 6)
            l[ori] = float(l[ori]) + SHSP
            l[end] = float(l[end]) - SHSP
            if l[end] > (l[ori] + HSPminLEN - 1):
                print(l)
                str_l = [ str(e) for e in l ]
                # print(*l, sep = "\t")
                # for e in l:
                #     print(e)
                #     ee = str(15) + 'aaaa'
                #     print(ee)
                # # print(l)
                lines_to_write.append( '\t'.join(l) )


open(out_filename, 'wb').writelines(lines_to_write)
# f.writelines( "%s\n" % item for item in xrange(2**20) )
