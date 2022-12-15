#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

baseCommand: [ Rscript]

requirements:
  DockerRequirement:
<<<<<<< HEAD
    dockerPull: docker.io/zmahnoor/maw-r:1.0.7
=======
    dockerPull: zmahnoor/maw-r:1.0.7
>>>>>>> 45425f6 (only for one mzML file)

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
<<<<<<< HEAD
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
=======
       glob: $(runtime.outdir)
>>>>>>> 45425f6 (only for one mzML file)
