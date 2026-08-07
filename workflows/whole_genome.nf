#!/usr/bin/env nextflow

include { HYBRID_WGS } from "../subworkflows/hybrid_wgs"
include { HIFI_ONLY_WGS } from "../subworkflows/hifi_only_wgs"

workflow WHOLE_GENOME {
    
    take:
        ch_hybrid_samples
        ch_hifi_only_samples

    main:
        HYBRID_WGS ( ch_hybrid_samples )
        HIFI_ONLY_WGS ( ch_hifi_only_samples )

}
