"""This file is a template for the configuration of the ASAB package."""

from pathlib import Path
from uuid import uuid4
from datetime import datetime

config = dict()

######################
### regular config ###
######################

projectFolder = Path("<path to the home directory of the project>").resolve()
experimentFolder = Path(__file__).resolve().parents[1]
run_ID = f"{datetime.now().strftime('%Y%m%d_%H%M%S')}_{str(uuid4())}"
runFolder = str(Path(experimentFolder).joinpath("Runs", run_ID))


simulationBalance = False
simulationArduino = False


config['projectFolder'] = projectFolder

config['experimentFolder'] = experimentFolder

config['run_ID'] = run_ID

config["runFolder"] = runFolder

config["balanceDriver"] =   {
                                "serialPort": "<COM port for the balance>",
                                "settings": {
                                    'baudrate': 115200,
                                    'bytesize': 8,
                                    'parity': 'N',
                                    'stopbits': 1,
                                    'xonxoff': False,
                                    'dsrdtr': False,
                                    'rtscts': False,
                                    'timeout': 0,
                                    'write_timeout': None,
                                    'inter_byte_timeout': None
                                    },
                                "simulated": simulationBalance
                            }

config["balance"] = {
                        "simulated": simulationBalance
                    }

config["CetoniDeviceDriver"] =  {
                                    "availableSyringes": str(Path(runFolder).joinpath(
                                                    "generated_files",
                                                    "syringes.txt"
                                                )
                                            ),
                                    "configPath": "<path to the configuration file for"
                                        "the formulation unit of ASAB>",
                                    "syringeConfig": {
                                        <dictionary stating the volume of the syringes
                                        using key e.g. "A0.0" and value e.g. "10_mL">
                                        },
                                    "valvePositionDict": str(
                                        Path(runFolder).joinpath(
                                                "generated_files",
                                                "valvePositionDict.py"
                                            )
                                        ),
                                    "flow": <float specifying the maximum
                                        total flow rate in mL/s>,
                                    "simulateBalance": simulationBalance,
                                }

config["CetoniDevice"] =    {
                                "waste": "<the node, where the waste is connected>",
                                "gas": "<the node, where the ambient can be accessed>",
                                "pathsToClean": [],
                                "measureVolumes": {
                                    "densiVisco": <volume needed for a density
                                        and viscosity measurement in mL>,
                                    "nmr": <volume needed for an NMR measurement in mL>,
                                    "electrochemistry": <volume needed for an
                                        electrochemical measurement in mL>
                                    },
                                "pickUpNode": "<node, where the formulated sample "
                                    "should be picked up>",
                                "wasteForDrying": "<the waste used, when blowing"
                                    "pressurized gas through the tubes for drying>"
                            }

config["densiViscoDriver"] = {
                            "inputFolder": <path to the directory, where inputs for 
                                density and viscosity measurements are placed>,
                            "outputFolder": <path to the directory, where outputs for 
                                density and viscosity measurements shall be placed>
                        }


config["densiVisco"] = {
                            "inputFolder": <path to the directory, where inputs for 
                                density and viscosity measurements are placed>,
                            "outputFolder": <path to the directory, where outputs for 
                                density and viscosity measurements shall be placed>
                        }

config["syringes"] =    {
                            "savePath": str(Path(runFolder).joinpath(
                                            "generated_files",
                                            "syringes.py"
                                            )
                                        ),
                        }

config["graph"] =   {
                        "pathNodes": "<path to the nodes.csv file>",
                        "pathEdges": "<path to the edges.csv file>",
                        "pathTubing": "<path to the tubing.csv file>",
                        "savePath_graph": str(Path(runFolder).joinpath(
                                            "generated_files",
                                            "graph.py"
                                            )
                                        ),
                        "saveNameValvePositionDict": str(
                                        Path(runFolder).joinpath(
                                                "generated_files",
                                                "valvePositionDict.py"
                                            )
                                        ),
                    }

config["ArduinoDriver"] =   {
                                "Arduino":  [{
                                                "ArduinoNo": <number to identify
                                                    the Arduino>,
                                                "serialPort": '<COM port for the Arduino>',
                                                "settings": {
                                                    'baudrate': 9600,
                                                    'timeout': 1
                                                    },
                                                "simulated": simulationArduino
                                            }],
                                "ArduinoRelay": [{
                                                    "relayNo": <relay number for
                                                        the Arduino>,
                                                    "ArduinoNo": <number to identify
                                                        the Arduino>,
                                                }],
                                "ArduinoValve": [{
                                                    "name": '<name for the valve>',
                                                    "positions": [<available positions for the valve
                                                         for the valve>],
                                                    "ArduinoNo": <number to identify
                                                        the Arduino>,
                                                    "relayNo": <relay number for
                                                        the Arduino>
                                                }]
                            }

config["solutionHandler"] = {
                                'stockSolutions': "<path to the stockSolutions.csv file>",
                                'chemicals': "<path to the chemicals.csv file>",
                                'units': {
                                    'molar_mass': 'g/mol',
                                    'density': 'g/cm^3',
                                    'mass': 'g'
                                    },
                                'chemicalsDict': str(Path(runFolder).joinpath(
                                            "generated_files",
                                            "chemicalsDict.py"
                                        )
                                    ),
                                'solutionsDict': str(Path(runFolder).joinpath(
                                            "generated_files",
                                            "solutionsDict.py"
                                            )
                                        )
                            }

config["PalmsensDriver"] = {
    "PalmsensSDK_python": "<path to the PalmSens Python SDK>",
    "savepath_raw": "<path, where results from EIS measurements shall be saved>"
}

config["electrochemistry"] = {
    "cell_constant": <float value of the cell constant in cm**-1>
}
