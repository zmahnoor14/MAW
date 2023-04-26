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
    mzml_file:
        type: File
        #format: http://edamontology.org/format_3244
    gnps_file:
        type: File
    hmdb_file:
        type: File
    mbank_file:
        type: File
    mzml_result:
        type: string
    file_id:
        type: string
    ppmx:
        type: int
    db_name:
        type: string
    db_path:
        type: File
    # isotope:
    #     type: boolean
    #     default: False
  
steps:
    dereplication:
        run: maw-r.cwl
        in:
            workflow_script: r_script
            mzml_file: mzml_file
            gnps_file: gnps_file
            hmdb_file: hmdb_file
            mbank_file: mbank_file
            mzml_result: mzml_result
            file_id: file_id
            ppmx: ppmx
            db_name: db_name
            db_path: db_path
        out:
            # - ms_files
            - results
            - peaks_and_parameters

    metfrag:
        run: maw-metfrag.cwl
        scatter:
            - PeakList
            - IonizedPrecursorMass
            - PrecursorIonMode
            # - IsPositiveIonMode
            - LocalDatabase
            - SampleName

        scatterMethod: dotproduct
        in:
            PeakList: 
                source: dereplication/peaks_and_parameters
                valueFrom: $(self.PeakList)
            IonizedPrecursorMass:
                source: dereplication/peaks_and_parameters
                valueFrom: $(self.IonizedPrecursorMass)
            PrecursorIonMode:
                source: dereplication/peaks_and_parameters
                valueFrom: $(self.PrecursorIonMode)
            # IsPositiveIonMode:
            #     source: dereplication/peaks_and_parameters
            #     valueFrom: $(self.IsPositiveIonMode)
            LocalDatabase:
                source: dereplication/peaks_and_parameters
                valueFrom: $(self.LocalDatabase)
            SampleName:
                source: dereplication/peaks_and_parameters
                valueFrom: $(self.SampleName)

        out: [candidate_list]

    # sirius_isotope:
    #     run: sirius-new.cwl
    #     in:
    #         spectrum: dereplication/ms_files
    #         isotope: 
    #             default: False
    #         #parameter: dereplication/parameters
    #     scatter:
    #         - spectrum
    #         #- parameter
    #     out: [results]

    # sirius_no_isotope:
    #     run: sirius-new.cwl
    #     in:
    #         spectrum: dereplication/ms_files
    #         isotope: 
    #             default: True
    #         #parameter: dereplication/parameters
    #     scatter:
    #         - spectrum
    #         #- parameter
    #     out: [results]

    cheminformatics:
        run: maw-py.cwl
        in: 
            workflow_script: python_script
            mzml_files_results: dereplication/results
            # sirius_results: 
            #     source: [sirius_no_isotope/results, sirius_isotope/results]
            #     linkMerge: True
            # candidate_list: metfrag/candidate_list
 
            
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
    StepInputExpressionRequirement: {}
