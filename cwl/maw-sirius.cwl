#!/usr/bin/env cwl-runner
#docker run -ti -v /maw:/maw zmahnoor/run-sirius4 bash
#/usr/local/sirius/bin/sirius --input 1_NA_iso_NA_MS1p_182.081_SIRIUS_param.ms --output output formula --profile orbitrap --no-isotope-filter --no-isotope-score --candidates 30 --ppm-max 5 --ppm-max-ms2 15 structure --database ALL canopus 2>&1 | tee -a /tmp/sirius.log > /dev/null

cwlVersion: v1.0
class: CommandLineTool

baseCommand: [sirius]

requirements:
  DockerRequirement:
    dockerPull: docker.io/zmahnoor/run-sirius4

inputs: 
  spectrum:
    type: File
#    inputBinding: 
#      prefix: --input
  profile:
    type: string
    default: "orbitrap"
  candidates:
    type: int
    default: 30
  ppm_max:
    type: int
    default: 5
  ppm_max_ms2:
    type: int
    default: 15
  database:
    type: string
    default: "ALL"

arguments:
    - --input
    - $(inputs.spectrum.path)
    - --output
    - $(inputs.spectrum.nameroot).json
    - formula
    - --profile
    - $(inputs.profile)
    - --no-isotope-filter
    - --no-isotope-score
    - --candidates
    - $(inputs.candidates)
    - --ppm-max
    - $(inputs.ppm_max)
    - --ppm-max-ms2
    - $(inputs.ppm_max_ms2)
    - structure
    - --database
    - $(inputs.database)
    - canopus

outputs:
  results: 
    type: Directory
    outputBinding:
       glob: $(inputs.spectrum.nameroot).json
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
