#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

baseCommand: [ Rscript, Basic-Preprocessing-Script.r ]

hints:
  DockerRequirement:
    dockerPull: zmahnoor/maw-r-basic-preprocessing
requirements:
  InitialWorkDirRequirement:
    listing:
      - $(inputs.mzml_files)
inputs:
  mzml_files:
    type: File[]
    inputBinding: 
      position: 1

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
