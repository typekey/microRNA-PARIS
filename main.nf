#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

WorkflowMain.initialise(workflow, params, log)

include { PARIS2 } from './workflows/paris2'

workflow ZL {
    PARIS2 ()
}

workflow {
    ZL ()
}
