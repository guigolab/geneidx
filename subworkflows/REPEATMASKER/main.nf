//
// Check input samplesheet and get read channels
//

include { REPEATMASKER_STAGELIB } from '../../modules/local/repeatmasker/stagelib'
include { REPEATMASKER_REPEATMASK } from '../../modules/local/repeatmasker/repeatmask'
include { FASTASPLITTER } from '../../modules/local/fastasplitter'
include { CAT_FASTA as REPEATMASKER_CAT_FASTA} from '../../modules/local/cat/fasta'
include { GUNZIP } from '../../modules/nf-core/modules/gunzip/main'

workflow REPEATMASKER {

      //tuple id, path
    take:
    assemblies
    rm_lib

    main:

    splitted_fastas = FASTASPLITTER(assemblies, params.npart_size)

    // splitted_fastas | view
    // If chunks == 1, forward - else, map each chunk to the meta hash
    splitted_fastas.branch { id,f ->
       single: f.getClass() != ArrayList
       multi: f.getClass() == ArrayList
    }.set { ch_fa_chunks }

    // ch_fa_chunks.single.view { "$it single" }
    // ch_fa_chunks.multi.view { "$it multi" }

    ch_fa_chunks.multi.flatMap { id, fastas ->
       fastas.collect { [ id, file(it)] }
    }.set { ch_chunks_split }

    // ch_chunks_split | view {it -> print "$it"}

     stage_lib = REPEATMASKER_STAGELIB(
         rm_lib,
         file(params.dummy_gff)
      )
    
    //tuple id, fasta chunk, library directory, repeatmodeler fasta 
    repeat_masker_input = 
        ch_fa_chunks.single
            .map { id,f -> [id,file(f)]}
            .mix(ch_chunks_split)
            .combine(stage_lib, by:0)
            .combine(rm_lib, by:0)
    
    masked_assemblies = REPEATMASKER_REPEATMASK(repeat_masker_input).groupTuple() | REPEATMASKER_CAT_FASTA


    masked_assemblies.view()

    // emit:
    // fasta = REPEATMASKER_CAT_FASTA.out.fasta

}