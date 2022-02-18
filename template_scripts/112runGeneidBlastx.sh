general_scripts="/nfs/users/rg/fcalvet/projects/scripts_to_run_Francisco"
ben_scripts="/nfs/users/rg/fcalvet/projects/benchmark_abinitio/scripts"
scripts="/nfs/users/rg/fcalvet/projects/vertebrate_benchmarking/scripts"
data="/nfs/users/rg/fcalvet/projects/vertebrate_benchmarking/data"

param_file_SGP="/nfs/users/rg/fcalvet/data/params/human3iso.sgp.Hs-Mm_NEW.param"
param_file="/nfs/users/rg/fcalvet/data/params/human3isoU12.param";

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


    # run Geneid+BLASTx
    mkdir -p $species_dir_blastx_o;   
    
    $general_scripts/sgp_divideSGPjobs4cluster_fromGFF.sh $species $species_dir_blastx_o $species_dir_genome $main_genome_file.fa $param_file_SGP ${main_genome_file} $reference_DB_name ${main_genome_file}.$reference_DB_name.blastx.out $species_dir_blastx_o/blastx

    echo "GeneidBLASTX sent";
    echo


done < $1
#done < /nfs/users/rg/fcalvet/projects/vertebrate_benchmarking/scripts/geneid_blastx_missing_files
#done < /nfs/users/rg/fcalvet/projects/vertebrate_benchmarking/scripts/downloaded_names_genomes_good_running

#done < /nfs/users/rg/fcalvet/projects/vertebrate_benchmarking/data/species_names_test

