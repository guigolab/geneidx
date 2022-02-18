general_scripts="/nfs/users/rg/fcalvet/projects/scripts_to_run_Francisco"
data="/nfs/users/rg/fcalvet/projects/vertebrate_benchmarking/data"

# contains the protein database setted up for running BLASTx against GENCODE proteins
protein_DBS="/nfs/users/rg/fcalvet/projects/geneidblastx_multispecies/data/protein_dbs"
reference_DB_name="Vertebrates.5references"


while read species main_genome_file; do
    echo "Starting $species";

    species_dir="$data/$species.files";
    species_dir_genome="$species_dir/dna";

    species_dir_blastx_o="$species_dir/geneidblastx"
    #species_dir_blastx_o="$species_dir/geneidblastx_unmasked"
    

    if [ ! -s  $species_dir/dna/$main_genome_file.fa ]; then

        echo "unzipping genome $main_genome_file.fa.gz"
        gunzip -q $species_dir/dna/$main_genome_file.fa.gz;

    fi


    # run DIAMOND BLASTx
    mkdir -p $species_dir_blastx_o;
    mkdir -p $species_dir_blastx_o/blastx;
    
    
    # we run the genome of each specie against the human gencode proteins
    $general_scripts/sgp_DIAMONDBLASTX4CLUSTER_toGFF.sh $species $main_genome_file.fa $species_dir_genome $species_dir_blastx_o/blastx $protein_DBS $reference_DB_name

    echo "DIAMOND blastx sent";
    echo


done < $1
#done < /nfs/users/rg/fcalvet/projects/vertebrate_benchmarking/scripts/geneid_blastx_missing_files
#done < /nfs/users/rg/fcalvet/projects/vertebrate_benchmarking/scripts/downloaded_names_genomes_good_running

#done < /nfs/users/rg/fcalvet/projects/vertebrate_benchmarking/data/species_names_test

