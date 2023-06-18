#!/usr/bin/env cwl-runner
#docker run -ti -v /maw:/maw zmahnoor/run-sirius4 bash
#/usr/local/sirius/bin/sirius --input 1_NA_iso_NA_MS1p_182.081_SIRIUS_param.ms --output output formula --profile orbitrap --no-isotope-filter --no-isotope-score --candidates 30 --ppm-max 5 --ppm-max-ms2 15 structure --database ALL canopus 2>&1 | tee -a /tmp/sirius.log > /dev/null

cwlVersion: v1.2 #has to be 1.2 or higher version for condition
class: Workflow

inputs: 
  spectrum:
    type: File
    # default:
    #   class: File
    #   path: "*.json"
  profile:
    type: string
    default: "orbitrap"
  candidates:
    type: int
    default: 30
  ppm_max:
    type: int
    default: 5
  ppm_max_ms2:
    type: int
    default: 15
  database:
    type: string
    default: "coconut"

steps:
  
  isotopefiles:
    run: sirius_isotope.cwl
    in: 
      spectrum: spectrum
      profile: profile
      candidates: candidates
      ppm_max: ppm_max
      ppm_max_ms2: ppm_max_ms2
      database: database
    when:
      $(/isotopeNum/.test(inputs.spectrum.path))
      # glob: 
      #   "*_isotopeNum_*.ms"
    scatter:
      - spectrum
    out: [results]
      # type: Directory
      # outputBinding:
      #   glob: $(inputs.spectrum.nameroot).json

  noisotopefiles:
    run: sirius_no_isotope.cwl
    in: 
      spectrum: spectrum
      profile: profile
      candidates: candidates
      ppm_max: ppm_max
      ppm_max_ms2: ppm_max_ms2
      database: database
    when:
      $(/NA_iso/.test(inputs.spectrum.path))
      # glob: 
      #   "*_NA_iso_*.ms"
    scatter:
      - spectrum
    out: [results]
      # type: Directory
      # outputBinding:
      #   glob: $(inputs.spectrum.nameroot).json

outputs:
  results:
    type: Directory
#      outputBinding:
#       glob: $(inputs.spectrum.nameroot).json
  # WorkflowStepOutput:
  #   type: Directory
  #   outputSource: []
    outputSource:
      - isotopefiles/results
      - noisotopefiles/results
    pickValue: first_non_null
  
requirements:
  InlineJavascriptRequirement: {}
  ScatterFeatureRequirement: {}
