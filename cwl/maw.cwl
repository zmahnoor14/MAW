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
    mzml_files: Directory
  
steps:
    dereplication:
        run: maw-r.cwl
        in:
            workflow_script: r_script
            mzml_files: mzml_files
        out: [results]
    cheminformatics:
        run: maw-py.cwl
        in: 
            workflow_script: python_script
            mzml_files_results: dereplication/results
        out: [results]

outputs:
  results: 
    type: Directory
    outputSource: cheminformatics/results
