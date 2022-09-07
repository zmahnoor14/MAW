#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: Rscript
hints:
  DockerRequirement:
    dockerPull: maw_r
requirements:
  InitialWorkDirRequirement:
    listing:
      - $(inputs.mzml_files)
inputs:
  src:
    type: File
    inputBinding:
      position: 1
  mzml_files:
    type: File[]
    inputBinding:
      position: 2

outputs:
  results: Directory
  # csv_output:
  #   type:
  #     type: array
  #     items: File
  #   outputBinding:
  #     glob: "*.csv"
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
