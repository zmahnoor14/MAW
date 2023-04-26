#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

baseCommand: [python3]

requirements:
  DockerRequirement:
    dockerPull: zmahnoor/maw-py:1.0.7

  InitialWorkDirRequirement:
    listing:
      - entry: $(inputs.mzml_files_results)
        writable: true
inputs: 
  workflow_script: File
  mzml_files_results:
    type: Directory
  # sirius_results:
  #   type: Directory[]
  #   inputBinding:
  #     position: 3
  


arguments: 
    - $(inputs.workflow_script.path)
    - $(inputs.mzml_files_results.path)
outputs:
  results: 
    type: Directory
    outputBinding:
       glob: $(inputs.mzml_files_results.basename)
  provenance:
    type: File
    outputBinding:
       glob: provenance_python.yaml


  
