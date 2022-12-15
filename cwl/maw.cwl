cwlVersion: v1.0
class: Workflow


inputs: 
    python_script:
      type: File
      default:
        class: File
        path: Workflow_Python_Script_all.py
    r_script:
       type: File
       default:
         class: File
         path: Workflow_R_Script_all.r
    mzml_file: File
  
steps:
    dereplication:
        run: maw-r.cwl
        in:
            workflow_script: r_script
            mzml_file: mzml_file
        out:
            - ms_files
            - ms1_data
            - provenance
            #- parameters
    sirius:
        run: maw-sirius.cwl
        in:
            ms_file: dereplication/ms_files
            #parameter: dereplication/parameters
        scatter:
            - ms_file
            #- parameter
        out:
            [results]
    cheminformatics:
        run: maw-py.cwl
        in: 
            workflow_script: python_script
            ms1_data: dereplication/ms1_data
            sirius_results: sirius/results
        out:
            - [results]
            - provenance

outputs:
  results: 
    type: Directory
    outputSource: cheminformatics/results
