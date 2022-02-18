#!/bin/bash

### I have defined the path to DATA before, using this command:
### DATA="/users/rg/fcalvet/projects/HSPs/data";

#$DATA/DIAMONDBLASTX4CLUSTER_auto.sh 
#    H.sapiens
#    Homo_sapiens.GRCh38.dna_rm.primary_assembly.fa
#    $DATA/H.sapiens
#    $DATA/H.sapiens/blastx
#    /no_backup_isis/rg/fcamara/Uniprot4BLAST
#    uniref90_notransposon_201702


# source ~/.fcamara_bashrc

export SPECIES=${1:?"please type SPECIES to be annotated; BulbMite"}

export GENOME=${2:?"please select genome in fasta format wiht fa extension or without extension"}

export GENOMEDIR=${3:?"please include directory with genome file"}

export BLASTXOUT=${4:?"BLASTX outputdir"}

export DATABASE=${5:?"directory containing the biological database; /no_backup_isis/rg/fcamara/NR_database2 "}

export DATABASENAME=${6:?"name of protein database; nr_201609 "}

mkdir -p $BLASTXOUT


cat <<EOF > $BLASTXOUT/NQS.${GENOME}.$DATABASENAME
#!/bin/bash
#$ -q rg-el7,long-sl7,short-sl7
#$ -N BLAST_TAXONOMY
#$ -e $BLASTXOUT/e.${GENOME}.$DATABASENAME.log
#$ -o $BLASTXOUT/o.${GENOME}.$DATABASENAME.log
#$ -l virtual_free=30G
#$ -l h_rt=10:00:00

genome_noFA=$( basename $GENOME .fa )
db_noDMND=$( basename $DATABASENAME .dmnd )

fmt6_custom='6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qstrand qframe'

if [ ! -s  $BLASTXOUT/\$genome_noFA.\$db_noDMND.blastx.out ]; then

    echo "create DIAMOND output blastx.out";

    diamond blastx --db $DATABASE/$DATABASENAME --query $GENOMEDIR/\$genome_noFA.fa -k 0 --max-hsps 0  --outfmt \$fmt6_custom --evalue 0.0001 --threads 1 --out $BLASTXOUT/\$genome_noFA.\$db_noDMND.blastx.out


else

    echo "DIAMOND output blastx.out already created"

fi



if [ ! -s  $BLASTXOUT/\$genome_noFA.\$db_noDMND.hsp.gff ]; then

    echo "create HSPs gff file";

    gawk 'BEGIN{OFS="\t"}{if (\$13=="-"){frame=-(\$14);print \$1,"blastx","hsp",\$8,\$7,\$12,\$13,frame,\$2}else if(\$13=="+"){print \$1,"blastx","hsp",\$7,\$8,\$12,\$13,\$14,\$2}}' $BLASTXOUT/\$genome_noFA.\$db_noDMND.blastx.out | sort -k1,1 -k4,5n -k9,9 | egrep -v '^MT' | cut -d '|' -f1 > $BLASTXOUT/\$genome_noFA.\$db_noDMND.hsp.gff

else

    echo "DIAMOND output GFF already created"

fi


EOF

chmod 755 $BLASTXOUT/NQS.${GENOME}.$DATABASENAME
qsub $BLASTXOUT/NQS.${GENOME}.$DATABASENAME >> $BLASTXOUT/listJobIDs.${SPECIES}.${DATABASENAME}  2>&1


