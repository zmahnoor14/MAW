#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

baseCommand: [python3]

requirements:
  DockerRequirement:
    dockerPull: zmahnoor/maw-py:1.0.7
  InlineJavascriptRequirement: {}
  # InitialWorkDirRequirement: {}
    # listing:
    #   - entry: $(inputs.mzml_files_results)
    #     writable: true
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
    - $(runtime.outdir)/cwl.output.json

outputs:
  results: 
    type: File
    outputBinding:
       glob: "mergedResults-with-one-Candidates.csv"
  provenance:
    type: File
    outputBinding:
       glob: provenance_python.yaml


  
