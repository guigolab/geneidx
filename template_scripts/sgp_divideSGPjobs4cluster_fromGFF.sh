#!/bin/bash

general_scripts="/nfs/users2/rg/fcalvet/projects/scripts_to_run_Francisco";


#source /users/rg/fcamara/.bashrc


export SPECIES=${1:?"please type SPECIES to be annotated"}

export ODIR=${2:?"please type directory path of $SPECIES project"}

export GENOMEDIR=${3:?"please type directory path of the GENOME of $SPECIES"}

export GENOME=${4:?"please indicate the file name of the genome assembly for $SPECIES"}

export PARAM=${5:?"please type ${SOURCENAME} $SPECIES parameter FILE"}

export SOURCENAME=${6:?"please give name to source of evidence"}

export TYPE=${7:?"indicate if predictions are for evaluation or genomewide prediction"}

export BLASTX=${8:?"indicate name of BLASTX output file"}

export BLASTX_dir=${9:?"indicate name of the directory with the BLASTX file"}




export PERL5LIB=/users/rg/fcamara/Perlmods/:$PERL5LIB;




mkdir -p $BLASTX_dir

cd $BLASTX_dir

# echo "process DIAMOND BLASTX if not already done"

if [ ! -s  $BLASTX_dir/$SOURCENAME.$TYPE.hsp.gff ]; then

    echo "Check the previous step, HSPs gff file is empty or does not exist";

#    gawk 'BEGIN{OFS="\t"}{if ($13=="-"){frame=-($14);print $1,"blastx","hsp",$8,$7,$12,$13,frame,$2}else if($13=="+"){print $1,"blastx","hsp",$7,$8,$12,$13,$14,$2}}' $BLASTX | sort -k1,1 -k4,5n -k9,9 | egrep -v '^MT' | cut -d '|' -f1 > $SOURCENAME.$TYPE.hsp.gff

fi



# index
if [ ! -s  $GENOMEDIR/${GENOME}.i ]; then
    echo "indexing $GENOME"
    fastaindex -f $GENOMEDIR/$GENOME -i $GENOMEDIR/${GENOME}.i
fi



cd $GENOMEDIR;

if [ ! -s  $GENOMEDIR/seq_id_$GENOME.txt ]; then

    echo "creating seq_id seq_id_$GENOME.txt"

    infoseq -noheading -only -name -length $GENOMEDIR/$GENOME | gawk '{print $1"_1_"$2,$1,$2}' | shuf > $GENOMEDIR/seq_id_$GENOME.txt

    python /users/rg/fcalvet/projects/scripts_to_run_Francisco/distribute_seqs.py $GENOMEDIR/seq_id_$GENOME.txt 25;

    ls -1 seq_id_${GENOME}.txt_?? >  $GENOMEDIR/list_chunks_$GENOME.txt

    echo "created list_chunks_$GENOME.txt"

else

    echo "previously created seq_id_$GENOME.txt"
    echo "previously created list_chunks_$GENOME.txt"
    
fi






#rm -rf $ODIR/${SOURCENAME}fasta_all

mkdir -p $ODIR/${SOURCENAME}.${TYPE}.fasta_all

export QUERYDIR=$ODIR/${SOURCENAME}.${TYPE}.fasta_all



# If the directory is already present :

if [ -d  $ODIR/${SOURCENAME}.${TYPE}.2pred ]; then
    rm -rf $ODIR/${SOURCENAME}.${TYPE}.2pred
    mkdir -p  $ODIR/${SOURCENAME}.${TYPE}.2pred
else
    mkdir -p  $ODIR/${SOURCENAME}.${TYPE}.2pred
fi

export GOUTPUT=$ODIR/${SOURCENAME}.${TYPE}.2pred



# creating folder for the time files
mkdir -p  $GOUTPUT/times 


while read chunk #xaa

do 

 cat <<EOF > $GOUTPUT/NQS.$chunk.g
#!/bin/bash
#$ -q long-sl7,short-sl7
#$ -N GBTX_${SPECIES}_$chunk
#$ -e $GOUTPUT/e.${chunk}.log
#$ -o $GOUTPUT/o.${chunk}.log
#$ -l virtual_free=4G
#$ -l h_rt=4:00:00

#. /users/rg/fcamara/.bashrc

export GOUTPUT=$GOUTPUT
export QUERYDIR=$QUERYDIR

cd $GOUTPUT

while read chr chr2 len
   do
 
    fastafetch -f $GENOMEDIR/$GENOME -i $GENOMEDIR/${GENOME}.i -q "\$chr2" > $QUERYDIR/\$chr2 

    echo -e "\nrun $TYPE for \$chr2 of ${SOURCENAME}\n"

    time (egrep -w "^\$chr2" $BLASTX_dir/$SOURCENAME.$TYPE.hsp.gff > $SPECIES.\$chr2.hsp.gff) 2>> $GOUTPUT/times/time.$SOURCENAME.$TYPE.${SPECIES}.$chunk

    time (blast2gff -vg $SPECIES.\$chr2.hsp.gff > $SPECIES.\$chr2.SR.gff) 2>> $GOUTPUT/times/time.$SOURCENAME.$TYPE.${SPECIES}.$chunk

    time ($general_scripts/sgp_getHSPSR.pl "\$chr2" < $SPECIES.\$chr2.SR.gff > $SPECIES.\$chr2.HSP_SR.gff) 2>> $GOUTPUT/times/time.$SOURCENAME.$TYPE.${SPECIES}.$chunk


    rm $SPECIES.\$chr2.hsp.gff
    rm $SPECIES.\$chr2.SR.gff
    

    time (geneid -3P $PARAM -S $SPECIES.\$chr2.HSP_SR.gff $QUERYDIR/\$chr2 | sed -e 's/geneid_v1.4/geneidblastx/g' | egrep 'CDS' | sort -k4,5n  >> $SPECIES.$SOURCENAME.${chunk}.gff3) 2>> $GOUTPUT/times/time.$SOURCENAME.$TYPE.${SPECIES}.$chunk

    rm $SPECIES.\$chr2.HSP_SR.gff


    echo "### FINISHED "\$chr2; # &

    rm $QUERYDIR/\$chr2;

done < $GENOMEDIR/$chunk

echo "### FINISHED "$chunk;

EOF

   chmod 755 $GOUTPUT/NQS.$chunk.g
   qsub $GOUTPUT/NQS.$chunk.g >> $GOUTPUT/listJobIDs.${SOURCENAME}.${SPECIES}  2>&1
  # rm $GOUTPUT/NQS.$chunk.g

done < $GENOMEDIR/list_chunks_$GENOME.txt
