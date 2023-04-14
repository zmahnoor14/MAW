cwlVersion: v1.0
class: CommandLineTool

baseCommand: ['java', '-jar', '/usr/src/myapp/MetFragCommandLine-2.5.0.jar']

requirements:
  DockerRequirement:
    dockerPull: docker.io/zmahnoor/maw-metfrag_2.5.0:1.0.4
  InlineJavascriptRequirement: {}

inputs:
    PeakListPath:
        type: File

    IonizedPrecursorMass:
        type: string

    PrecursorIonMode:
        type: int

    # IsPositiveIonMode:
    #     type: boolean
    MetFragDatabaseType:
        type: string

    LocalDatabasePath:
        type: File

    DatabaseSearchRelativeMassDeviation:
        type: int

    FragmentPeakMatchAbsoluteMassDeviation:
        type: float

    FragmentPeakMatchRelativeMassDeviation:
        type: float

    MetFragCandidateWriter:
        type: string

    SampleName:
        type: string
    ResultsPath:
        type: Directory
    MetFragPreProcessingCandidateFilter:
        type: string
    MetFragPostProcessingCandidateFilter:
        type: string
    MaximumTreeDepth:
        type: int
    NumberThreads:
        type: int

outputs:
  metfrag_results:
    type: File
    outputBinding:
      glob: "$(inputs.SampleName).csv"
