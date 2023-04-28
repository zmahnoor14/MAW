cwlVersion: v1.0
class: CommandLineTool

baseCommand: ['java', '-jar', '/usr/src/myapp/MetFragCommandLine-2.5.0.jar']

requirements:
  DockerRequirement:
    dockerPull: docker.io/zmahnoor/maw-metfrag_2.5.0:1.0.4
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
    - entryname: metfrag.inputs
      entry: |-
        PeakListPath = $(inputs.PeakList.path)
        IonizedPrecursorMass =  $(inputs.IonizedPrecursorMass)
        PrecursorIonMode = $(inputs.PrecursorIonMode)
        MetFragDatabaseType = LocalCSV
        LocalDatabasePath = $(inputs.LocalDatabasePath.path)
        DatabaseSearchRelativeMassDeviation = 5
        FragmentPeakMatchAbsoluteMassDeviation = 0.001
        FragmentPeakMatchRelativeMassDeviation = 15
        MetFragCandidateWriter = CSV
        SampleName = $(inputs.SampleName)
        ResultsPath = $(runtime.outdir)/metfrag
        MetFragPreProcessingCandidateFilter = UnconnectedCompoundFilter
        MetFragPostProcessingCandidateFilter = InChIKeyFilter
        MaximumTreeDepth = 2
        NumberThreads = $(runtime.cores)


inputs: # additional inputs for all files; make them to show certain paths
  PeakList:
      type: File
  IonizedPrecursorMass:
      type: string
  PrecursorIonMode:
      type: int
  LocalDatabasePath:
      type: File
  SampleName:
      type: string

# MetFragDatabaseType:
  #     type: string
  # DatabaseSearchRelativeMassDeviation:
  #     type: int

  # FragmentPeakMatchAbsoluteMassDeviation:
  #     type: float

  # FragmentPeakMatchRelativeMassDeviation:
  #     type: float

  # MetFragCandidateWriter:
  #     type: string

  # MetFragPreProcessingCandidateFilter:
  #     type: string
  # MetFragPostProcessingCandidateFilter:
  #     type: string
  # MaximumTreeDepth:
  #     type: int
  # NumberThreads:
  #     type: int

arguments: 
  - cwl.output.json

outputs:

  metfrag_candidate_list:
    type: File
    outputBinding:
        glob: "$(runtime.outdir)/metfrag/*.csv"
        #loadContents: true