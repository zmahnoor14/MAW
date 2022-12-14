#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

baseCommand: [ Rscript]

requirements:
  DockerRequirement:
    dockerPull: docker.io/zmahnoor/maw-r:1.0.7


  InitialWorkDirRequirement:
    listing:
      - entry: $(inputs.mzml_files)
        writable: true
inputs: 
  workflow_script: File
  mzml_files:
    type: Directory

arguments: 
    - $(inputs.workflow_script.path)
    - $(inputs.mzml_files.basename)
outputs:
  results: 
    type: Directory
    outputBinding:
       glob: $(inputs.mzml_files.basename)
  # csv_output:
  #   type:
  #     type: array
  #     items: File
  #   
  # mzml_outputs:
  #   type:
  #     type: array
  #     items: File
  #   outputBinding:
  #     glob: "*/*.mzML"
  # txt_outputs:
  #   type:
  #     type: array
  #     items: File
  #   outputBinding:
  #     glob: "*/*.txt"

#  result_files:
#    type:
#      type: array
#      items: File
