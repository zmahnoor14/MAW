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
  msp_file: File
  gnps_dir: Directory
  hmdb_dir: Directory
  mbank_dir: Directory
  metfrag_candidate_list: 
    type: File[]
  ms1data: File
  score_thresh:
    type: float
    default: 0.75
    
arguments: 
    - $(inputs.workflow_script.path)
    - cwl.output.json

outputs:
  results: 
    type: Directory
    outputBinding:
       glob: $(inputs.mzml_files_results.basename)
  provenance:
    type: File
    outputBinding:
       glob: provenance_python.yaml


  
