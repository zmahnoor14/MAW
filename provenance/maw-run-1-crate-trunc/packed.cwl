{
    "$graph": [
        {
            "class": "CommandLineTool",
            "baseCommand": [
                "python3"
            ],
            "requirements": [
                {
                    "dockerPull": "zmahnoor/maw-py:1.0.7",
                    "class": "DockerRequirement"
                },
                {
                    "listing": [
                        {
                            "entry": "$(inputs.mzml_files_results)",
                            "writable": true
                        }
                    ],
                    "class": "InitialWorkDirRequirement"
                }
            ],
            "inputs": [
                {
                    "type": "Directory",
                    "id": "#maw-py.cwl/mzml_files_results"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "Directory"
                    },
                    "inputBinding": {
                        "position": 3
                    },
                    "id": "#maw-py.cwl/sirius_results"
                },
                {
                    "type": "File",
                    "id": "#maw-py.cwl/workflow_script"
                }
            ],
            "arguments": [
                "$(inputs.workflow_script.path)",
                "$(inputs.mzml_files_results.path)"
            ],
            "id": "#maw-py.cwl",
            "outputs": [
                {
                    "type": "Directory",
                    "outputBinding": {
                        "glob": "provenance-intro-artifacts"
                    },
                    "id": "#maw-py.cwl/provenance"
                },
                {
                    "type": "Directory",
                    "outputBinding": {
                        "glob": "$(inputs.mzml_files_results.basename)"
                    },
                    "id": "#maw-py.cwl/results"
                }
            ]
        },
        {
            "class": "CommandLineTool",
            "baseCommand": [
                "Rscript"
            ],
            "requirements": [
                {
                    "dockerPull": "zmahnoor/maw-r:1.0.7",
                    "class": "DockerRequirement"
                }
            ],
            "inputs": [
                {
                    "type": "File",
                    "id": "#maw-r.cwl/gnps_rda"
                },
                {
                    "type": "File",
                    "id": "#maw-r.cwl/hmdb_rda"
                },
                {
                    "type": "File",
                    "id": "#maw-r.cwl/mbank_rda"
                },
                {
                    "type": "File",
                    "id": "#maw-r.cwl/mzml_files"
                },
                {
                    "type": "File",
                    "id": "#maw-r.cwl/workflow_script"
                }
            ],
            "arguments": [
                "$(inputs.workflow_script.path)",
                "$(inputs.mzml_files.path)",
                "$(inputs.gnps_rda.path)",
                "$(inputs.hmdb_rda.path)",
                "$(inputs.mbank_rda.path)",
                "."
            ],
            "outputs": [
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputBinding": {
                        "glob": "insilico/SIRIUS/*.ms"
                    },
                    "id": "#maw-r.cwl/ms_files"
                },
                {
                    "type": "Directory",
                    "outputBinding": {
                        "glob": "prov_console"
                    },
                    "id": "#maw-r.cwl/provenance"
                },
                {
                    "type": "Directory",
                    "outputBinding": {
                        "glob": "."
                    },
                    "id": "#maw-r.cwl/results"
                }
            ],
            "id": "#maw-r.cwl"
        },
        {
            "class": "CommandLineTool",
            "baseCommand": [
                "sirius"
            ],
            "requirements": [
                {
                    "dockerPull": "docker.io/zmahnoor/run-sirius4",
                    "class": "DockerRequirement"
                }
            ],
            "inputs": [
                {
                    "type": "int",
                    "default": 30,
                    "id": "#maw-sirius.cwl/candidates"
                },
                {
                    "type": "string",
                    "default": "ALL",
                    "id": "#maw-sirius.cwl/database"
                },
                {
                    "type": "int",
                    "default": 5,
                    "id": "#maw-sirius.cwl/ppm_max"
                },
                {
                    "type": "int",
                    "default": 15,
                    "id": "#maw-sirius.cwl/ppm_max_ms2"
                },
                {
                    "type": "string",
                    "default": "orbitrap",
                    "id": "#maw-sirius.cwl/profile"
                },
                {
                    "type": "File",
                    "id": "#maw-sirius.cwl/spectrum"
                }
            ],
            "arguments": [
                "--input",
                "$(inputs.spectrum.path)",
                "--output",
                "$(inputs.spectrum.nameroot).json",
                "formula",
                "--profile",
                "$(inputs.profile)",
                "--no-isotope-filter",
                "--no-isotope-score",
                "--candidates",
                "$(inputs.candidates)",
                "--ppm-max",
                "$(inputs.ppm_max)",
                "--ppm-max-ms2",
                "$(inputs.ppm_max_ms2)",
                "structure",
                "--database",
                "$(inputs.database)",
                "canopus"
            ],
            "outputs": [
                {
                    "type": "Directory",
                    "outputBinding": {
                        "glob": "$(inputs.spectrum.nameroot).json"
                    },
                    "id": "#maw-sirius.cwl/results"
                }
            ],
            "id": "#maw-sirius.cwl"
        },
        {
            "class": "Workflow",
            "inputs": [
                {
                    "type": "File",
                    "id": "#main/gnps_rda"
                },
                {
                    "type": "File",
                    "id": "#main/hmdb_rda"
                },
                {
                    "type": "File",
                    "id": "#main/mbank_rda"
                },
                {
                    "type": "File",
                    "id": "#main/mzml_files"
                },
                {
                    "type": "File",
                    "default": {
                        "class": "File",
                        "path": "file:///mnt/tdm-dic/users/simleo/sandbox/MAW/cwl/Workflow_Python_Script_all.py"
                    },
                    "id": "#main/python_script"
                },
                {
                    "type": "File",
                    "default": {
                        "class": "File",
                        "path": "file:///mnt/tdm-dic/users/simleo/sandbox/MAW/cwl/Workflow_R_Script_all.r"
                    },
                    "id": "#main/r_script"
                }
            ],
            "steps": [
                {
                    "run": "#maw-py.cwl",
                    "in": [
                        {
                            "source": "#main/dereplication/results",
                            "id": "#main/cheminformatics/mzml_files_results"
                        },
                        {
                            "source": "#main/sirius/results",
                            "id": "#main/cheminformatics/sirius_results"
                        },
                        {
                            "source": "#main/python_script",
                            "id": "#main/cheminformatics/workflow_script"
                        }
                    ],
                    "out": [
                        "#main/cheminformatics/results"
                    ],
                    "id": "#main/cheminformatics"
                },
                {
                    "run": "#maw-r.cwl",
                    "in": [
                        {
                            "source": "#main/gnps_rda",
                            "id": "#main/dereplication/gnps_rda"
                        },
                        {
                            "source": "#main/hmdb_rda",
                            "id": "#main/dereplication/hmdb_rda"
                        },
                        {
                            "source": "#main/mbank_rda",
                            "id": "#main/dereplication/mbank_rda"
                        },
                        {
                            "source": "#main/mzml_files",
                            "id": "#main/dereplication/mzml_files"
                        },
                        {
                            "source": "#main/r_script",
                            "id": "#main/dereplication/workflow_script"
                        }
                    ],
                    "out": [
                        "#main/dereplication/ms_files",
                        "#main/dereplication/results"
                    ],
                    "id": "#main/dereplication"
                },
                {
                    "run": "#maw-sirius.cwl",
                    "in": [
                        {
                            "source": "#main/dereplication/ms_files",
                            "id": "#main/sirius/spectrum"
                        }
                    ],
                    "scatter": [
                        "#main/sirius/spectrum"
                    ],
                    "out": [
                        "#main/sirius/results"
                    ],
                    "id": "#main/sirius"
                }
            ],
            "outputs": [
                {
                    "type": "Directory",
                    "outputSource": "#main/cheminformatics/results",
                    "id": "#main/results"
                }
            ],
            "requirements": [
                {
                    "class": "ScatterFeatureRequirement"
                }
            ],
            "id": "#main"
        }
    ],
    "cwlVersion": "v1.0"
}