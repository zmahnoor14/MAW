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
    gnps_file:
        type: File
    hmdb_file:
        type: File
    mbank_file:
        type: File
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
    dereplication:
        run: maw-r.cwl
        in:
            workflow_script: r_script
            mzml_file: mzml_file
            gnps_file: gnps_file
            hmdb_file: hmdb_file
            mbank_file: mbank_file
            ppmx: ppmx
            collision_info: collision_info
            db_name: db_name
            db_path: db_path
        out:
            # - ms_files
            - results
            - peaks_and_parameters
            - msp_file
            - ms1data
            - gnps_dir
            - hmdb_dir
            - mbank_dir

    metfrag:
        run: maw-metfrag.cwl
        scatter:
            - PeakList
            - IonizedPrecursorMass
            - PrecursorIonMode
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
            SampleName:
                source: dereplication/peaks_and_parameters
                valueFrom: $(self.SampleName)
            LocalDatabasePath: db_path

        out: [metfrag_candidate_list]
        #out: metfrag_candidate_list 

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
            msp_file: dereplication/msp_file
            gnps_dir: dereplication/gnps_dir
            hmdb_dir: dereplication/hmdb_dir
            mbank_dir: dereplication/mbank_dir
            # sirius_results: 
            #     source: [sirius_no_isotope/results, sirius_isotope/results]
            #     linkMerge: True
            metfrag_candidate_list: metfrag/metfrag_candidate_list
            ms1data: dereplication/ms1data

        #out: [results, provenance]
        out: 
            # - msp_file_df
            # - ms1data_df
            - candidate_files
            - result

outputs:
    candidate_files:
        type: File[]
        outputSource: cheminformatics/candidate_files
    result: 
        type: File
        outputSource: cheminformatics/result

requirements:
    ScatterFeatureRequirement: {}
    StepInputExpressionRequirement: {}
    InlineJavascriptRequirement: {}

