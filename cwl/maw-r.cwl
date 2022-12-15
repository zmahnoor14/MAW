#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

baseCommand: [ Rscript]

requirements:
  DockerRequirement:
    dockerPull: zmahnoor/maw-r:1.0.7

inputs: 
  workflow_script: File
  mzml_files:
    type: File
  gnps_rda:
    type: File
  hmdb_rda:
    type: File
  mbank_rda:
    type: File

arguments: 
    - $(inputs.workflow_script.path)
    - $(inputs.mzml_files.path)
    - $(inputs.gnps_rda.path)
    - $(inputs.hmdb_rda.path)
    - $(inputs.mbank_rda.path)
    - $(runtime.outdir)
outputs:
  results: 
    type: Directory
    outputBinding:
       glob: $(runtime.outdir)