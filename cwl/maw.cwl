cwlVersion: v1.0
class: Workflow
doc: |
    multiple mzml files, all with the same condition

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
    mzml_files:
        type: File[]
    gnps_file:
        type: File
    hmdb_file:
        type: File
    mbank_file:
        type: File
    file_id:
        type: string
    ppmx:
        type: int
    collision_info:
        type: boolean
    db_name:
        type: string
    db_path:
        type: File
    # isotope:
    #     type: boolean
    #     default: False
  
steps:
    analysis:
        run: maw_single.cwl
        in:
            python_script: python_script
            r_script: r_script
            mzml_file: mzml_files
            gnps_file: gnps_file
            hmdb_file: hmdb_file
            mbank_file: mbank_file
            file_id: file_id
            ppmx: ppmx
            collision_info: collision_info
            db_name: db_name
            db_path: db_path
        scatter: mzml_file
        out: 
            - candidate_files
            - result

outputs:
    candidate_files:
        type: 
            type: array
            items: 
                type: array
                items: File
        outputSource: analysis/candidate_files
    result: 
        type: File[]
        outputSource: analysis/result

requirements:
    ScatterFeatureRequirement: {}
    SubworkflowFeatureRequirement: {}

