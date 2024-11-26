/*
    Useful helpers to get versioning of softwares
    Adapted from: 
    https://github.com/nf-core/rnaseq/blob/master/subworkflows/nf-core/utils_nfcore_pipeline/main.nf
*/

def print_finish_summary() {
    println ( workflow.success ? 
    """
    ===============================================================================
    Pipeline execution summary
    -------------------------------------------------------------------------------

    Run as      : ${workflow.commandLine}
    Started at  : ${workflow.start}
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    Config files: ${workflow.configFiles}
    exit status : ${workflow.exitStatus}

    --------------------------------------------------------------------------------
    ================================================================================
    """.stripIndent() : """
    Failed      : ${workflow.errorReport}
    Exit status : ${workflow.exitStatus}
    Run name    : ${workflow.runName}
    """.stripIndent()
    )
}

def print_start_params() {
    log.info """
    --------------------------------------------------------------------------------
    DGE-ANALYSIS
    
    Differential Gene Expression Analysis using paired-ended RNA sequencing data

    tonyliang19/biof501-project
    ---------------------------------------------------------------------------------
    Parameter / Metadata information of the pipeline
    =================================================================================
    Core Nextflow options
        runName             : ${workflow.runName}
        containerEngine     : ${workflow.containerEngine}
        launchDir           : ${workflow.launchDir}
        workDir             : ${workflow.workDir}
        projectDir          : ${workflow.projectDir}
        profile             : ${workflow.configFiles}

    Input/output options
        samplesheet         : ${params.samplesheet}
        outdir              : ${params.outdir}
    
    Main options
        genome              : ${params.genome}
        genome_annotation   : ${params.genome_annotation}

    Max resource options
        max_cpus            : ${params.max_cpus}
        max_memory          : ${params.max_memory}
        max_time            : ${params.max_time}
    ================================================================================ 
    For more information please see the github README at:
    https://github.com/tonyliang19/biof501-project/blob/main/README.md

    Starting the pipeline now with above options...

    """.stripIndent()
}

def getWorkflowVersion() {
    def version_string = "" as String
    if (workflow.manifest.version) {
        def prefix_v = workflow.manifest.version[0] != 'v' ? 'v' : ''
        version_string += "${prefix_v}${workflow.manifest.version}"
    }

    if (workflow.commitId) {
        def git_shortsha = workflow.commitId.substring(0, 7)
        version_string += "-g${git_shortsha}"
    }

    return version_string
} 

def processVersionsFromYAML(yaml_file) {
    def yaml = new org.yaml.snakeyaml.Yaml()
    def versions = yaml.load(yaml_file).collectEntries { k, v -> [ k.tokenize(':')[-1], v ] }
    return yaml.dumpAsMap(versions).trim()
}

def workflowVersionToYAML() {
    return """
    Workflow:
        $workflow.manifest.name: ${getWorkflowVersion()}
        Nextflow: $workflow.nextflow.version
    """.stripIndent().trim()
}

def softwareVersionsToYAML(ch_versions) {
    return ch_versions
                .unique()
                .map { version -> processVersionsFromYAML(version) }
                .unique()
                .mix(Channel.of(workflowVersionToYAML()))
}
