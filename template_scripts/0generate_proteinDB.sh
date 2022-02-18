data="/users/rg/fcalvet/projects/geneidblastx_multispecies/data"
fa_prot="$data/fasta_proteins"
prot_dbs="$data/protein_dbs"

peptide_file="Homo_sapiens.GRCh38.pep.all.fa"
#peptide_file="Gallus_gallus.GRCg6a.pep.all.fa"
#peptide_file="Drosophila_melanogaster.r6.43.protein.all.fasta"
peptide_file_no_end=$(echo $peptide_file | rev | cut -d '.' -f2- | rev)
#species=""
species=$(echo $peptide_file | cut -d '.' -f1)


egrep 'biotype:protein_coding' $fa_prot/$peptide_file | egrep -v 'mitochondria' | sed 's/>//g' | gawk '{print $1}' | sort | uniq > $fa_prot/$species.protein_coding_Curated


FastaToTbl $fa_prot/$peptide_file | sort -k1,1  > $fa_prot/$peptide_file_no_end.tbl


join -1 1 -2 1 $fa_prot/$peptide_file_no_end.tbl $fa_prot/$species.protein_coding_Curated | TblToFasta > $fa_prot/$peptide_file_no_end.proteincoding.fasta
#(just protein-coding plus no mitochondrial proteins)


rm $fa_prot/$peptide_file_no_end.tbl;
rm $fa_prot/$species.protein_coding_Curated;
rm $fa_prot/$peptide_file;

#Format gencode proteins for Diamond:
#/software/rg/el7.2/fcamara/diamond  # now copied to ~/bin/diamond
diamond makedb --in $fa_prot/$peptide_file_no_end.proteincoding.fasta -d $prot_dbs/$peptide_file_no_end

