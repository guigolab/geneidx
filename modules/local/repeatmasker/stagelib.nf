process REPEATMASKER_STAGELIB {

    tag "$fasta"
    label 'process_low'

    container 'quay.io/biocontainers/repeatmasker:4.1.2.p1--pl5321hdfd78af_1'
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/repeatmasker:4.1.2.p1--pl5321hdfd78af_1':
    //     'quay.io/biocontainers/repeatmasker:4.1.2.p1--pl5321hdfd78af_1' }"

    input:
    tuple val(id), path(fasta)
    path db

    output:
    tuple val(id), path("Libraries")
    
    script:

    options = "-lib $fasta"
    
    """
    cp ${baseDir}/assets/repeatmasker/my_genome.fa .
    cp ${baseDir}/assets/repeatmasker/repeats.fa .
    cp -R /usr/local/share/RepeatMasker/Libraries .
    cp $db Libraries/Dfam.h5
    export LIBDIR=\$PWD/Libraries
    RepeatMasker $options my_genome.fa > out
    """
}
