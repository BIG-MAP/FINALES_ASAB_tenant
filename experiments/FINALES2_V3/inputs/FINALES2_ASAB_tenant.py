from datetime import datetime, timedelta
from uuid import uuid4
from typing import Any, Union
import pathlib
from itertools import combinations
import json
from pathlib import Path
import pandas as pd
import numpy as np
import time
import os
import logging

from madap.echem.e_impedance import e_impedance
from madap.data_acquisition import data_acquisition as da

parent_folder = pathlib.Path(__file__).parent.resolve()
import sys
sys.path.append(parent_folder)
import configFINALES2_ASAB_tenant
conf = configFINALES2_ASAB_tenant.config
from ASAB_tenant_config import config as tenant_config

from ASAB.utility.helpers import loadVariable, saveToFile
from ASAB.utility import solutionHandler, loggingASAB
from ASAB.action import CetoniDevice_action, Palmsens_action

from FINALES2.tenants.referenceTenant import Tenant
from FINALES2.schemas import GeneralMetaData, Quantity, ServerConfig, Method
from FINALES2.user_management.classes_user_manager import User
from FINALES2.server.schemas import Request, RequestInfo

runLogger = loggingASAB.prepareMainLogger()
# Create a logger for the tenant
logger_ASAB_tenant = logging.getLogger(f"run_logger.{__name__}")
logger_ASAB_tenant.info(f"Added the path {parent_folder} to the system path")
logger_ASAB_tenant.info(f"Imported the configuration {str(conf)}.")

logger_ASAB_tenant.info("Running with the following commits and branches:\n"
    f"FINALES2_repository:\n branch: {tenant_config['FINALES_repository']['branch']}\n "
    f"commit: {tenant_config['FINALES_repository']['commit']}\n"

    f"FINALES2_schemas_repository:\n branch: {tenant_config['FINALES_schemas_repository']['branch']}\n "
    f"commit: {tenant_config['FINALES_schemas_repository']['commit']}\n"

    f"ASAB_repository:\n branch: {tenant_config['ASAB_repository']['branch']}\n "
    f"commit: {tenant_config['ASAB_repository']['commit']}\n"

    f"ASAB_experiments_repository:\n branch: {tenant_config['ASAB_experiments_repository']['branch']}\n "
    f"commit: {tenant_config['ASAB_experiments_repository']['commit']}\n"

    f"MADAP_repository:\n branch: {tenant_config['MADAP_repository']['branch']}\n "
    f"commit: {tenant_config['MADAP_repository']['commit']}\n"
)

########################################################################################
#### Helper functions ####
########################################################################################

# Copy the live_config to the runFolder
os.popen(f"copy {str(Path(conf['experimentFolder']).joinpath('live_config.py'))} {str(Path(conf['runFolder']).joinpath('inputs', 'live_config.py'))}")


def get_live_config(
        logger:logging.Logger,
        filepath:str=str(Path(conf["runFolder"]).joinpath("inputs", "live_config.py")),
        config_var:str="live_config"
        ):
    live_config = loadVariable(loadPath=filepath, variable=config_var)
    logger.info(f"Got the latest info:\n {live_config}")
    return live_config

def update_vial_status(
        vial:str,
        new_status:Union[str, None],
        logger:logging.Logger,
        config_filepath:str=str(Path(conf["runFolder"]).joinpath("inputs", "live_config.py")),
        config_var:str="live_config"
        ):
    live_config = get_live_config(logger, filepath=config_filepath, config_var=config_var)
    live_config["vial_status"][vial].append((datetime.now().strftime("%Y_%m_%d-%H_%m_%S"), new_status))
    logger.info(f"Updated status of {vial} to be {new_status}.")


########################################################################################
#### Preparing the system ####
########################################################################################

# Prepare the hardware for operation
setupData = CetoniDevice_action.prepareRun(graphSave=True, graphShow=False, refPos=True)

# Copy the initial information to the run folder
live_config_path = str(Path(conf['runFolder']).joinpath('inputs', 'live_config'))
os.popen(f"copy {Path(conf['experimentFolder']).joinpath('live_config')} {live_config_path}")
# Log the path to the live_config
logger_ASAB_tenant.info(f"Live config is saved here:\n {live_config_path}")

pumpsDict = setupData["pumps"]
valvesDict = setupData["valves"]

chemicalsDict = setupData["chemicalsDict"]
solutionsDict = setupData["solutionsDict"]

########################################################################################
#### Get limitations and vials ####
########################################################################################

def get_chemicals_limits(solutionsDict:dict, pumpsDict:dict, chemicalsDict:dict):

    solution_names = list(solutionsDict.keys())

    limitations_total = []

    all_combis_no_repeat = list(combinations(list(solution_names), len(solution_names)))

    # For each combination
    for combi in all_combis_no_repeat:
        # find the chemicals involved
        chemicals_in_combi = []
        for sol_name in combi:
            chemicals_in_combi.extend(solutionsDict[sol_name].chemicals)
        chemicals_in_combi = set(chemicals_in_combi)

        reduced_solutions_dict = solutionsDict.copy()
        for key in solutionsDict.keys():
            if key not in combi:
                del reduced_solutions_dict[key]

        combi_limits = []

        # find max. and min. amount for each chemical
        for chemical in chemicals_in_combi:
            mix_ratio_max = {f"{chemical}": 1.}

            fractions_info_max = solutionHandler.get_fractions_2(
                mixingRatio=mix_ratio_max,
                chemicalsDict=chemicalsDict,
                solutionsDict=reduced_solutions_dict,
                input_fraction="molPerMol",
                output_fraction="volPerVol"
                )

            fractions_max = fractions_info_max[2]   # TODO: include limitations due to flow of syringes
            
            # fractions_min = fractions_info_min[2]   # TODO: include limitations due to flow of syringes

            fraction_max = fractions_max.get(chemical, 0.0)
            if chemicalsDict[chemical].name == "EMC":
                fraction_min = 0.35 # based on https://iopscience.iop.org/article/10.1149/1.1393419/pdf figure 3 at 20 °C;
                # the actual acceptable fraction should be lower due to the addition of LiPF6, so this is only a rough estimate
            else:
                fraction_min = 0.0

            # assemble the limitations structure
            component_limit = {
                "chemical": {
                    "SMILES": chemicalsDict[chemical].SMILES,
                    "InChIKey": chemical
                },
                "fraction": [{"min": fraction_min, "max": fraction_max}],
                "fraction_type": ["molPerMol"]
            }
            combi_limits.append(component_limit)
        # add the structure to the list of limitations
        limitations_total.append(combi_limits)

    logger_ASAB_tenant.info(f"limitations_total: \n {limitations_total}")

    # return the list of limitations
    return limitations_total

formulation_limits = get_chemicals_limits(
    solutionsDict=solutionsDict,
    chemicalsDict=chemicalsDict,
    pumpsDict=pumpsDict
    )

def get_free_vial():
    live_config = get_live_config(logger=runLogger)
    for vial in live_config["vial_status"].keys():
        if live_config["vial_status"][vial][-1][1] is None:
            logger_ASAB_tenant.info(f"free_vial: \n {vial}")
            return vial
    logger_ASAB_tenant.info(f"free_vial: \n {vial}")
    return None

def get_temperature(T_type:str):
    live_config = get_live_config(logger=runLogger)
    if T_type == "housing":
        T = live_config["temperature_echem_cell_housing"][-1][1]
    elif T_type == "glovebox":
        T = live_config["ambient_temperature_glovebox"][-1][1]
    else:
        raise ValueError(f"Unknown type of temperature {T_type}.")
    logger_ASAB_tenant.info(f"temperature_{T_type}: \n {T}")
    return T

def get_rating_mix(mixRatio_requested:dict, mixRatio_actual:dict):
    mixRatio_requested = pd.Series(data=mixRatio_requested, index=mixRatio_requested.keys())
    mixRatio_actual = pd.Series(data=mixRatio_actual, index=mixRatio_actual.keys())
    print(
        "REQ:\n", mixRatio_requested,
        "ACTUAL\n", mixRatio_actual,
        "DIFF\n", mixRatio_requested - mixRatio_actual
    )
    deviation = np.sqrt(np.sum((mixRatio_requested - mixRatio_actual)**2))
    mix_requested = np.sqrt(np.sum((mixRatio_requested)**2))
    error_fraction = deviation/mix_requested
    print("ERROR_FRACTION_MIX\n", error_fraction)
    if error_fraction >= 1.:
        quality_rating = 0.0
    else:
        quality_rating = np.floor((1. - error_fraction) * 5.)
    print("QUALITY_RATING_MIX\n", quality_rating)
    logger_ASAB_tenant.info(f"quality_rating_mix: \n {quality_rating}")
    return quality_rating

def get_rating_EIS(info:pd.DataFrame):
    for i in info.index:
        info.at[(i, "error_conductivity / S/cm")] = abs(- conf["electrochemistry"]["cell_constant"] / (info.at[(i, "resistivity / Ohm")]**2.)) * info.at[(i, "RMSE_fit")] # This focuses only on the resistance R0 in R0-CPE1, which is relevant for the calculation of the conductivity

    average_conductivity = info["conductivity_MADAP / S/cm"].sum() / len(info["conductivity_MADAP / S/cm"])
    sum_of_squares = np.sum([(info.at[(k, "conductivity_MADAP / S/cm")] - average_conductivity)**2. for k in info.index])
    RMSE_avg = np.sqrt(sum_of_squares)
    fit_error_RMSE_avg = np.sum([((sum_of_squares)**(-0.5) * (info.at[(l, "conductivity_MADAP / S/cm")] - average_conductivity)) * info.at[(l, "error_conductivity / S/cm")] for l in info.index])
    error_avg = RMSE_avg + fit_error_RMSE_avg

    error_fraction = error_avg / average_conductivity
    print("ERROR_FRACTION_EIS\n", error_fraction)
    if error_fraction >= 1.:
        quality_rating = 0.0
    else:
        quality_rating = np.floor((1. - error_fraction) * 5.)
    print("QUALITY_RATING_EIS\n", quality_rating)
    logger_ASAB_tenant.info(f"quality_rating_EIS: \n {quality_rating}")
    return quality_rating


# ########################################################################################
# #### Functions to run the methods ####
# ########################################################################################


# utility functions to cover repeating tasks across the methods
def fraction_check(formulation:dict) -> str:
    '''This method ensures, that the fraction type for all chemicals is the same within
    the given formulation and sum to unity.'''
    master_fraction_type = formulation[0]["fraction_type"]
    fractions = []
    for component in formulation:
        if not component["fraction_type"] == master_fraction_type:
            raise ValueError(
                f"The component {component['chemical']['InChIKey']} has the "
                f"fraction type {component['fraction_type']} instead of "
                f"{master_fraction_type}.")
        fractions.append(component["fraction"])
    if not np.isclose(np.sum(fractions), 1., atol=1e-15):
        raise ValueError(
            f"The sum of the fractions is {np.sum(fractions)}. "
            "Fractions are expected to sum to unity."
        )

    logger_ASAB_tenant.info(f"Fractions are all given in {master_fraction_type}.")

    return master_fraction_type

def get_mixing_ratio(formulation:dict) -> dict:
    '''This function extracts the mixing ratio in the format compatible with ASAB from
    the definition of a formulation in the format compatible with FINALES.'''
    mixingRatio={
                component["chemical"]["InChIKey"]: component["fraction"]
                for component in formulation
                }
    logger_ASAB_tenant.info(f"Requested mixing ratio is \n {mixingRatio}")
    return mixingRatio

def mix_ratio_to_formulation(mixing_ratio:dict, fraction_type:str) -> list:
    '''This function transfers a mixing ratio from ASAB format to the FINALES
    formulation format.'''
    formulation = []
    for chemical in mixing_ratio.keys():
        chem = {
            "chemical": {
                "SMILES": chemicalsDict[chemical].SMILES,
                "InChIKey": chemical
            },
            "fraction": mixing_ratio[chemical],
            "fraction_type": fraction_type
        }
        formulation.append(chem)
    return formulation


# # define the methods to run experiments
# def density_viscosity_measurement(parameters:dict):
#     # TODO: Ensure, that the fraction type for all chemicals is the same!!!
#     mixRatio_request, mixRatio_vol, calcMix_request, flows_actual, pumpsFlows = CetoniDevice_action.mix(
#             mixingRatio={
#                 component["chemical"]["InChIKey"]: component["fraction"]
#                 for component in parameters["formulation"]
#                 },
#             fraction=parameters["formulation"][0]["fraction_type"],
#             pumps=pumpsDict,
#             valves=valvesDict
#             )
#     actualFlows = CetoniDevice_action.provideSample(
#             measurementtype="densiVisco",
#             sample_node="M1.0",
#             pumpsFlows=pumpsFlows,
#             pumps=pumpsDict,
#             valves=valvesDict
#             )
#     densiVisco_action.measure(sampleName=uuid4(), method="Lovis-DMA_MultiMeasure_20")
#     CetoniDevice_action.drainSample(
#             measurementtype="densiVisco",
#             pump="F0.0",
#             repeats=3,
#             pumps=pumpsDict,
#             valves=valvesDict
#             )
#     density_result = densiVisco_action.retrieveData()
#     return density_result


def conductivity_measurement(
    parameters:dict,
    sample_name:str,
    request_ID:str,
    method_params:dict=tenant_config["EIS"]["parameters"],
    repeats:int=3
    ) -> dict:

    eis_method = "EIS"

    fraction_type = fraction_check(formulation=parameters["formulation"])
    mixing_ratio = get_mixing_ratio(formulation=parameters["formulation"])

    # Assign a UUID to the sample and the batch
    sample_ID = f"{request_ID}_sample_{sample_name}"
    batch_ID = f"{request_ID}_batch_{uuid4()}"

    logger_ASAB_tenant.info(f"Mixing sample {sample_ID} part of batch {batch_ID} for conductivity measurement.")

    conductivity_result_info = pd.DataFrame(columns=["resistivity / Ohm", "RMSE_fit", "conductivity_MADAP / S/cm", "conductivity / S/m"])

    mixRatio_request, mixRatio_vol, calcMix_request, flows_actual, pumpsFlows = CetoniDevice_action.mix(
            mixingRatio=mixing_ratio,
            fraction=fraction_type,
            pumps=pumpsDict,
            valves=valvesDict
            )
    quality_rating_mix = get_rating_mix(mixRatio_requested=mixRatio_request, mixRatio_actual=calcMix_request)
    print("QUALITY_RATING_MIX\n", quality_rating_mix)
    
    for i in range(repeats):
        subsample_name = f"{sample_name}_{i}"
        actualFlows = CetoniDevice_action.provideSample(
                measurementtype="electrochemistry",
                sample_node="M1.0",
                pumpsFlows=pumpsFlows,
                pumps=pumpsDict,
                valves=valvesDict
                )

        print("\n\n\n FLOWS actual comparison \n", flows_actual, "\n", actualFlows)
        print("\n\n\n FLOWS actual equal \n", flows_actual == actualFlows)
        logger_ASAB_tenant.info(f"Flows actual comparison \n mix: {flows_actual} \n provideSample: {actualFlows}")

        T_cell_housing = get_temperature(T_type="housing")

        logger_ASAB_tenant.info(f"EIS measurement of sample {sample_ID} in progress...")
        Palmsens_action.measure(callback_func=None,
            sampleName=subsample_name,
            method=eis_method,
            parameters=method_params,
            save_folder=conf["PalmsensDriver"]["savepath_raw"]
        )

        logger_ASAB_tenant.info(f"Collection of data for sample {sample_ID} in progress...")
        eis_data = Palmsens_action.retrieveData(
            sampleName=subsample_name,
            method=eis_method,
            save_folder=conf["PalmsensDriver"]["savepath_raw"]
        )

        for k in ["freq_arrays", "zre_arrays", "zim_arrays"]:
            if len(eis_data[k]) > 1:
                raise ValueError(
                    f"There are {len(eis_data[k])} lists of {k} in the EIS data. "
                    "Expected only one."
                    )

        ## Analyze the data to get the ionic conductivity
        # assemble the data for MADAP
        data_for_analysis = pd.DataFrame(
            data={
                "frequency": eis_data["freq_arrays"][0],
                "real_impedance": eis_data["zre_arrays"][0],
                "neg_imaginary_impedance": eis_data["zim_arrays"][0]
            }
        )
        # MADAP negates the imaginary impedance values internally. The Palmsens SDK already,
        # does the same, so the negative of the Palmsens output needs to be passed to MADAP.
        data_for_analysis["imaginary_impedance"] = -data_for_analysis["neg_imaginary_impedance"]

        analysis_savepath = Path(
            conf["PalmsensDriver"]["savepath_raw"]
            ).joinpath("MADAP_analysis")
        analysis_savepath.mkdir(parents=True, exist_ok=True)

        analysis_plots = ["nyquist", "nyquist_fit", "bode", "residual"]

        analysis_impedance = e_impedance.EImpedance(
            da.format_data(data_for_analysis["frequency"]),
            da.format_data(data_for_analysis["real_impedance"]),
            da.format_data(data_for_analysis["imaginary_impedance"])
            )

        analysis_EIS = e_impedance.EIS(
            analysis_impedance,
            suggested_circuit = tenant_config["EIS"]["circuit"]["circuit"],
            initial_value=tenant_config["EIS"]["circuit"]["initial_guess"](data_for_analysis),
            cell_constant = conf["electrochemistry"]["cell_constant"]
            )

        logger_ASAB_tenant.info(f"Analysis of EIS data for sample {sample_ID} in progress...")
        analysis_EIS.perform_all_actions(
            str(analysis_savepath),
            plots = analysis_plots,
            optional_name = f"{subsample_name}"
        )

        # Get the conductivity value form the output of MADAP
        result_filepath = ""
        path_glob = list(analysis_savepath.joinpath('data').glob("*.json"))
        while result_filepath == "":
            time.sleep(1)
            for file in path_glob:
                if subsample_name in str(file): # Should work as long as the sample_name is a unique identifier
                    result_filepath = file
        with open(result_filepath, "r") as analysis_file:
            analyzed_data = json.load(analysis_file)

        analysis_result = pd.DataFrame(
            {
               "resistivity / Ohm": analyzed_data["Parameters"][0],
               "RMSE_fit": analyzed_data["RMSE_fit_error"],
               "conductivity_MADAP / S/cm": analyzed_data["conductivity [S/cm]"],
               "conductivity / S/m": analyzed_data["conductivity [S/cm]"]*100.
            },
            index=[i]
        )

        logger_ASAB_tenant.info(f"analysis_result: \n {analysis_result}")

        conductivity_result_info = pd.concat([conductivity_result_info, analysis_result], ignore_index=False, axis=0)

    quality_rating_EIS = get_rating_EIS(info = conductivity_result_info)
    print("QUALITY_RATING_EIS\n", quality_rating_EIS)

    quality_rating = np.floor((quality_rating_mix + quality_rating_EIS)/2.)
    print("QUALITY_RATING\n", quality_rating)

    results = {
        "sample_name": sample_name,
        "sample_ID": sample_ID,
        "mixRatio_request": mixRatio_request,
        "mixRatio_vol": mixRatio_vol,
        "actual_mix_ratio": calcMix_request,
        "flows_actual": flows_actual,
        "pumps_flows": pumpsFlows,
        "actual_flows": actualFlows,
        "conductivity": conductivity_result_info["conductivity / S/m"].to_list(),
        "batch_ID": batch_ID,
        "temperature": T_cell_housing,
        "success": True,
        "quality_rating": quality_rating
        }

    logger_ASAB_tenant.info(f"Results for sample {sample_ID}: \n {results}")

    # Save the output
    saveToFile(
        folder=str(Path(conf["runFolder"]).joinpath("Results")),
        filename=f"{sample_ID}_conductivity_two_electrodes",
        extension="py",
        data=f"result = {str(results)}")

    # Drain the sample and clean the measuring cell
    logger_ASAB_tenant.info(f"Cleaning after sample {sample_ID} in progress...")
    CetoniDevice_action.clean(
        medium="Solvent",
        pumps=["F0.0"],
        pumpsDict=pumpsDict,
        valvesDict=valvesDict,
        nodes=["electrochemistryIN"],
        repeats=1,
        goToRef=False)
    CetoniDevice_action.clean(
        medium="pressurizedGas",
        pumpsDict=pumpsDict,
        valvesDict=valvesDict,
        nodes=["electrochemistryIN"],
        repeats=1,
        goToRef=False)

    return results


def formulate(parameters:dict, request_ID: str, sample_name:str) -> dict:
    fraction_type = fraction_check(formulation=parameters["electrolyte"]["formulation"])
    mixing_ratio = get_mixing_ratio(formulation=parameters["electrolyte"]["formulation"])

    # Assign a UUID to the sample and the batch
    sample_ID = f"{request_ID}_sample_{uuid4()}"
    batch_ID = f"{request_ID}_batch_{uuid4()}"

    logger_ASAB_tenant.info(f"Mixing sample {sample_ID} part of batch {batch_ID} for electrolyte preparation.")

    T_glovebox = get_temperature(T_type="glovebox")

    mixRatio_request, mixRatio_vol, calcMix_request, flows_actual, pumpsFlows = CetoniDevice_action.mix(
        mixingRatio=mixing_ratio,
        fraction=fraction_type,
        pumps=pumpsDict,
        valves=valvesDict
        )
    
    start = datetime.now()
    actualFlows = CetoniDevice_action.provideSample(
            measurementtype=None,
            volume=parameters["volume"],
            sample_node="V7.0",
            pumpsFlows=pumpsFlows,
            pumps=pumpsDict,
            valves=valvesDict,
            endpoint=parameters["vial"]
            )
    end = datetime.now()
    duration = (end-start) - timedelta(seconds=13)
    print("DURATION\n", duration, type(duration))
    vol_estimate = duration.total_seconds() * conf["CetoniDeviceDriver"]["flow"]
    print("VOL_estimate\n", vol_estimate, type(vol_estimate))

    logger_ASAB_tenant.info(f"Estimated volume of electrolyte is {vol_estimate}.")

    print("\n\n\n FLOWS actual comparison \n", flows_actual, "\n", actualFlows)
    print("\n\n\n FLOWS actual equal \n", np.round(list(flows_actual.values()), 3) == np.round(list(actualFlows.values()),3))
    logger_ASAB_tenant.info(f"Flows actual comparison \n mix: {flows_actual} \n provideSample: {actualFlows}")
    
    mixrat = {}
    for k in flows_actual.keys():
        mixrat.update({k: actualFlows[solutionsDict[k].pump]})

    quality_rating = get_rating_mix(mixRatio_requested=mixRatio_request, mixRatio_actual=calcMix_request)
    print("QUALITY_RATING\n", quality_rating)

    results = {
        "sample_name": sample_name,
        "sample_ID": sample_ID,
        "mixRatio_request": mixRatio_request,
        "mixRatio_vol": mixRatio_vol,
        "actual_mix_ratio": calcMix_request,
        "flows_actual": flows_actual,
        "pumps_flows": pumpsFlows,
        "actual_flows": actualFlows,
        "volume": vol_estimate,
        "vial": parameters["vial"],
        "batch_ID": batch_ID,
        "temperature": T_glovebox,
        "success": True,
        "quality_rating": quality_rating
    }

    logger_ASAB_tenant.info(f"Sample {sample_ID} in vial {parameters['vial']}")

    # Save the output
    saveToFile(
        folder=str(Path(conf["runFolder"]).joinpath("Results")),
        filename=f"{sample_name}_{sample_ID}_electrolyte_flow",
        extension="py",
        data=str(results))
    
    update_vial_status(vial=parameters["vial"], new_status=f"{request_ID}", logger=runLogger)

    return results

########################################################################################
#### Methods for ASAB ####
########################################################################################

def run_ASAB(input_request:RequestInfo):
    request_ID = input_request["uuid"]
    request = Request(**input_request["request"])
    method = request.methods[0]
    parameters = request.parameters[method]

    logger_ASAB_tenant.info(f"Request information: \n request: {request_ID} \n method: {method} \n parameters: {parameters}")

    print("\n request_ID", request_ID, "\n method", method, "\n parameters", parameters)
    sample_name = str(uuid4())

    logger_ASAB_tenant.info(f"request ID: {sample_name}")

    if method=="two_electrode":
            result = conductivity_measurement(parameters=parameters, sample_name=sample_name, request_ID=request_ID)
    elif method=="flow":
        free_vial = get_free_vial()
        if free_vial is not None:
            parameters["vial"] = free_vial
            update_vial_status(vial=free_vial, new_status=f"reserved_{request_ID}", logger=runLogger)
        else:
            raise ValueError("No free vials available, mixture not possible.")  # TODO: Find a smarter way to handle this situation, so that it does not break the tenant.
        result = formulate(parameters=parameters, sample_name=sample_name, request_ID=request_ID)
    # if method=="vibrating_tube_densimetry":
    #     result = density_viscosity_measurement(parameters=parameters)
    # elif method=="rolling_ball_viscosimetry":
    #     result = density_viscosity_measurement(parameters=parameters)
    return result

def prepare_results_ASAB(request:dict, data:Any):
    logger_ASAB_tenant.info(f"Preparing results for request {request['uuid']} of sample {data['sample_name']} for posting...")
    # Assemble the chemicals_info
    chemicals_info = {}
    for chem in data["actual_mix_ratio"].keys():
        chemicals_info[chem] = {
            "name": chemicalsDict[chem].name,
            "molar_mass": (chemicalsDict[chem].molar_mass, chemicalsDict[chem].molar_mass_uncertainty),
            "density": (chemicalsDict[chem].density, chemicalsDict[chem].density_uncertainty),
            "batch": chemicalsDict[chem].batch,
            "manufacturer": chemicalsDict[chem].manufacturer,
            "manufacturing_date": chemicalsDict[chem].manufacturing_date,
        }
    # extract the technical part of the request
    request_technical = Request(**request["request"])

    formulation = mix_ratio_to_formulation(
                data["actual_mix_ratio"],
                fraction_type="molPerMol"
                )

    run_info={
                "formulation": formulation,
                "internal_reference": data["sample_name"],
                "formulation_info": {
                    "name": data["sample_name"],
                    "uuid": data["sample_ID"],
                    "preparation_date": datetime.now().strftime("%d.%m.%Y"),
                    "batch": data["batch_ID"]
                },
                "chemicals_info": chemicals_info
    }

    # prepare the data depending on the method used
    if (("two_electrode" in request_technical.methods)
        and 
        ("conductivity" in request_technical.quantity)):

        method = "two_electrode"

        output = {
            "run_info": run_info,
            "conductivity": {
                "values": data["conductivity"],
                "temperature": data["temperature"],
                "meta": {
                    "success": True,
                    "rating": data["quality_rating"]
                }
            }
        }

        actual_parameters = {
            "formulation": formulation,
            "temperature": data["temperature"]
        }


    elif (("flow" in request_technical.methods)
        and 
        ("electrolyte" in request_technical.quantity)):

        method = "flow"

        output = {
            "run_info": run_info,
            "electrolyte": {
                "formulation": formulation,
                "volume": data["volume"],
                "location": {
                    "address": "KIT, HIU, Ulm",
                    "detail_1": "ASAB",
                    "detail_2": data["vial"]}
                },
                "meta": {
                    "success": True,
                    "rating": data["quality_rating"]
                }
        }
        actual_parameters = {
            "electrolyte": {
                "formulation": formulation
            },
            "volume": data["volume"]
            }

    result = {
        "data": output,
        "quantity": request_technical.quantity,
        "method": request_technical.methods,
        "parameters": {method: actual_parameters},
        "tenant_uuid": "",
        "request_uuid": request["uuid"]
    }

    logger_ASAB_tenant.info(f"Prepared result: \n {result}")

    return result

########################################################################################
#### Method objects for the methods ####
########################################################################################

# collect all the information required for the tenant
# density_method = Method(
#     name = "vibrating_tube_densimetry",
#     quantity = "density",
#     parameters = ["formulation", "temperature"],
#     limitations = {
#     "formulation": formulation_limits,
#     "temperature": [{"min": 298.15, "max": 298.15}]
#     }
# )

# viscosity_method = Method(
#     name = "rolling_ball_viscosimetry",
#     quantity = "viscosity",
#     parameters = ["formulation", "temperature"],
#     limitations = {
#         "formulation": formulation_limits,
#         "temperature": [{"min": 298.15, "max": 298.15}]
#     }
# )

conductivity_method = Method(
    name = "two_electrode",
    quantity = "conductivity",
    parameters = ["formulation", "temperature"],
    limitations = {
        "formulation": formulation_limits,
        "temperature": [
            {
                "min": get_temperature(T_type="housing"),
                "max": get_temperature(T_type="housing")
                }
        ]
    }
)

logger_ASAB_tenant.info(f"conductivity method: \n {conductivity_method}")

electrolyte_method = Method(
    name = "flow",
    quantity = "electrolyte",
    parameters = ["formulation"],
    limitations = {
        "formulation": formulation_limits
    }
)

logger_ASAB_tenant.info(f"electrolyte method: \n {electrolyte_method}")

########################################################################################
#### Inputs for the tenant ####
########################################################################################

ASAB_quantities = {
    # "density": Quantity(
    # name = density_method.quantity,
    # methods = {density_method.name: density_method},
    # is_active = True
    # ),
    # "viscosity":  Quantity(
    # name = viscosity_method.quantity,
    # methods = {viscosity_method.name: viscosity_method},
    # is_active = True
    # ),
    "conductivity": Quantity(
    name = conductivity_method.quantity,
    methods = {conductivity_method.name: conductivity_method},
    is_active = True
    ),
    "electrolyte": Quantity(
    name = electrolyte_method.quantity,
    methods = {electrolyte_method.name: electrolyte_method},
    is_active = True
    )
}

########################################################################################
#### Assembling the tenant ####
########################################################################################

ASAB_tenant = Tenant(
    general_meta = GeneralMetaData(**tenant_config["general_meta"]),
    sleep_time_s = 5,
    quantities = ASAB_quantities,
    tenant_config = (
        f"ASAB tenant configuration:\n {str(conf)}"
        f"\nThe specific config is:\n {str(tenant_config)}"
    ),
    run_method = run_ASAB,
    prepare_results = prepare_results_ASAB,
    FINALES_server_config = ServerConfig(**tenant_config["ServerConfig"]),
    end_run_time = tenant_config["end_run_time"],
    operators = [User(**tenant_config["operator"])],
    tenant_user = User(**tenant_config["tenant_user"]),
    tenant_uuid = "71e8e983f58c4ab1937f929fa641e2d8"
)

logger_ASAB_tenant.info(f"general metadata: \n {ASAB_tenant.general_meta}")
logger_ASAB_tenant.info(f"ASAB quantities: \n {ASAB_tenant.quantities}")
logger_ASAB_tenant.info(f"ASAB tenant configuration: {ASAB_tenant.tenant_config}")
logger_ASAB_tenant.info(f"ASAB FINALES Server configuration: \n {ASAB_tenant.FINALES_server_config}")
logger_ASAB_tenant.info(f"ASAB end run time: \n {ASAB_tenant.end_run_time}")

for user in ASAB_tenant.operators + [ASAB_tenant.tenant_user]:
    # Avoid logging the password
    user_copy = user.model_copy().__dict__
    del user_copy["password"]
    logger_ASAB_tenant.info(f"ASAB {user_copy['username']} (omitting password): \n {user_copy}")



if __name__ == "__main__":

    ASAB_tenant.run()
