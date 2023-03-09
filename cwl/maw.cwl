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
    mzml_files:
        type: File
        #format: http://edamontology.org/format_3244
    gnps_rda:
        type: File
    hmdb_rda:
        type: File
    mbank_rda:
        type: File
    isotope:
        type: boolean
        default: False
  
steps:
    dereplication:
        run: maw-r.cwl
        in:
            workflow_script: r_script
            mzml_files: mzml_files
            gnps_rda: gnps_rda
            hmdb_rda: hmdb_rda
            mbank_rda: mbank_rda
        out:
            - ms_files
            - results
    sirius_isotope:
        run: sirius-new.cwl
        in:
            spectrum: dereplication/ms_files
            isotope: 
                default: False
            #parameter: dereplication/parameters
        scatter:
            - spectrum
            #- parameter
        out: [results]

    sirius_no_isotope:
        run: sirius-new.cwl
        in:
            spectrum: dereplication/ms_files
            isotope: 
                default: True
            #parameter: dereplication/parameters
        scatter:
            - spectrum
            #- parameter
        out: [results]

    cheminformatics:
        run: maw-py.cwl
        in: 
            workflow_script: python_script
            mzml_files_results: dereplication/results
            sirius_results: 
                source: [sirius_no_isotope/results, sirius_isotope/results]
                linkMerge: True
            
        out: [results, provenance]

outputs:
  results: 
    type: Directory
    outputSource: cheminformatics/results
  cheminformatics_prov:
    type: File
    outputSource: cheminformatics/provenance
requirements:
    ScatterFeatureRequirement: {}
