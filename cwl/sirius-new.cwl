cwlVersion: v1.0 
class: CommandLineTool

baseCommand: [sirius]

requirements:
  DockerRequirement:
    dockerPull: docker.io/zmahnoor/run-sirius4
  InlineJavascriptRequirement: {}

inputs: 
  spectrum:
    type: File
#    inputBinding: 
#      prefix: --input
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
  isotope:
    type: boolean
    default: False


arguments:
    - --input
    - $(inputs.spectrum.path)
    - --output
    - $(inputs.spectrum.nameroot).json
    - formula
    - --profile
    - $(inputs.profile)
    - | 
      ${
        if (inputs.isotope) { 
          return ["--no-isotope-filter",
            "--no-isotope-score"];
        } else {
          return null;
        }
        
      }
    - --candidates
    - $(inputs.candidates)
    - --ppm-max
    - $(inputs.ppm_max)
    - --ppm-max-ms2
    - $(inputs.ppm_max_ms2)
    - fingerprint
    - structure
    - --database
    - $(inputs.database)
    - compound-classes
    - write-summaries
    - --output
    - $(inputs.spectrum.nameroot).json

outputs:
  results: 
    type: Directory
    outputBinding:
       glob: $(inputs.spectrum.nameroot).json
  