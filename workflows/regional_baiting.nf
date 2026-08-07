#!/usr/bin/env nextflow

include { HIFI_ONLY_REGIONAL } from "../subworkflows/hifi_only_regional"
include { HYBRID_REGIONAL } from "../subworkflows/hybrid_regional"

workflow REGIONAL_BAITING {

    take:
        ch_hybrid_samples
        ch_hifi_only_samples
        ch_ref
        ch_desired_regions

    main:
        HYBRID_REGIONAL (
            ch_hybrid_samples,
            ch_ref,
            ch_desired_regions
        )

        HIFI_ONLY_REGIONAL (
            ch_hifi_only_samples,
            ch_ref,
            ch_desired_regions
        )

}
