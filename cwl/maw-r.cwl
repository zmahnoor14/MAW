#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

baseCommand: ["Rscript"]                  

requirements:
  DockerRequirement:
    dockerPull: zmahnoor/maw-r:1.0.8
  InlineJavascriptRequirement: {}

inputs: 
  workflow_script: File
  mzml_file:
    type: File
  gnps_file:
    type: File
  hmdb_file:
    type: File
  mbank_file:
    type: File
  file_id:
    type: string
  ppmx:
    type: int
  # runCamera: boolean
  collision_info: boolean
  db_name:
    type: string
  db_path:
    type: File

arguments: 
  # - Rscript
  - $(inputs.workflow_script.path)
  - $(inputs.mzml_file.path)
  - $(inputs.gnps_file.path)
  - $(inputs.hmdb_file.path)
  - $(inputs.mbank_file.path)
  - $(inputs.file_id)
  - $(inputs.ppmx)
  - |
    $(inputs.collision_info ? "TRUE" : "FALSE")
  - $(inputs.db_name)
  - $(inputs.db_path.path)
# the no output binding needed bceause the output.json will provide output
outputs:
  results:
    type: Directory

  # ms_files_isotope:
  #   type: File[]?

  # ms_files_no_isotope:
  #   type: File[]

  peaks_and_parameters:
    type: Any[]?

  msp_file:
    type: File

  ms1data:
    type: File

  gnps_dir:
    type: Directory

  hmdb_dir:
    type: Directory

  mbank_dir:
    type: Directory
  
  provenance:
    type: Directory

# $namespaces:
#   edam: http://edamontology.org/