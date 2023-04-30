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
            file_id: file_id
            ppmx: ppmx
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
            - LocalDatabasePath
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
            LocalDatabasePath:
                source: dereplication/peaks_and_parameters
                valueFrom: $(self.LocalDatabasePath)
            SampleName:
                source: dereplication/peaks_and_parameters
                valueFrom: $(self.SampleName)

        out: [metfrag_candidate_list]

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
            msp_file: 
                source: dereplication/msp_file
                valueFrom: $(self.msp_file)
            gnps_dir: 
                source: dereplication/gnps_dir
                valueFrom: $(self.gnps_dir)
            hmdb_dir: 
                source: dereplication/hmdb_dir
                valueFrom: $(self.hmdb_dir)
            mbank_dir: 
                source: dereplication/mbank_dir
                valueFrom: $(self.mbank_dir)
            # sirius_results: 
            #     source: [sirius_no_isotope/results, sirius_isotope/results]
            #     linkMerge: True
            metfrag_candidate_list: 
                source: [metfrag/metfrag_candidate_list]
            ms1data: 
                source: dereplication/ms1data
                valueFrom: $(self.ms1data)
        out: [results, provenance]

outputs:
  results: 
    type: File
    outputSource: cheminformatics/results
  cheminformatics_prov:
    type: File
    outputSource: cheminformatics/provenance
requirements:
    ScatterFeatureRequirement: {}
    StepInputExpressionRequirement: {}
    InlineJavascriptRequirement: {}

# $namespaces:
#   edam: http://edamontology.org/