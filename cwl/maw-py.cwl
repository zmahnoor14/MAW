#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

baseCommand: [python3]

requirements:
  DockerRequirement:
    dockerPull: zmahnoor/maw-py:1.0.7
  InlineJavascriptRequirement: {}
inputs: 
  workflow_script: File
  msp_file: File
  gnps_dir: Directory
  hmdb_dir: Directory
  mbank_dir: Directory
  metfrag_candidate_list: 
    type: File[]
    inputBinding: 
      prefix: --metfrag_candidate_list
  ms1data: File
  score_thresh:
    type: float
    default: 0.75
    
arguments: 
    - $(inputs.workflow_script.path)
    - --msp_file
    - $(inputs.msp_file.path)
    - --gnps_dir
    - $(inputs.gnps_dir.path)
    - --hmdb_dir
    - $(inputs.hmdb_dir.path)
    - --mbank_dir
    - $(inputs.mbank_dir.path)
    - --ms1data
    - $(inputs.ms1data.path)
    - --score_thresh
    - $(inputs.score_thresh)

outputs:
  msp_file_df:
    type: File
    outputBinding:
       glob: "spec_results_for*.csv"
  ms1data_df:
    type: File
    outputBinding:
      glob: "MS1DATA.csv"
  # candidate_directory:
  #   type: Directory
  #   outputBinding:
  #     glob: "final_results/Candidate_Selection"
  result: 
    type: File
    outputBinding:
       glob: "mergedResults-with-one-Candidates.csv"
  # provenance:
  #   type: File
  #   outputBinding:
  #      glob: provenance_python.yaml


  
