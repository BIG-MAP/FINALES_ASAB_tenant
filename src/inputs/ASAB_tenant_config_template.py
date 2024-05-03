"""This file is a template for the configuration of the ASAB tenant for FINALES."""

import pandas as pd
from datetime import datetime, timedelta

def initial_guess(data:pd.DataFrame):
    return ["<initial guess for the selected equivalent circuit>"]

config = {}

config["FINALES_repository"] = {
    "branch": "<used branch>",
    "commit": "<used commit>"
}

config["FINALES_schemas_repository"] = {
    "branch": "<used branch>",
    "commit": "<used commit>"
}

config["ASAB_repository"] = {
    "branch": "<used branch>",
    "commit": "<used commit>"
}

config["ASAB_experiments_repository"] = {
    "branch": "<used branch>",
    "commit": "<used commit>"
}

config["MADAP_repository"] = {
    "branch": "<used branch>",
    "commit": "<used commit>"
}

config["EIS"] = {
    "parameters": {
        "<parameters for an EIS measurement as in the PalmSens SDK version 5.9>"
        },
    "circuit": {
        "circuit": "<equivalent circuit as accepted by the MADAP package>",
        "initial_guess": initial_guess
    }
}

config["general_meta"] = {
    "name": "ASAB_tenant",
    "description": "This tenant can formulate and analyze liquid electrolytes."
}

config["ServerConfig"] = {
    "host": "<host IP address>",
    "port": "<host port>"
}

config["end_run_time"] = datetime.now() + timedelta(minutes=1440) # 24 h

config["operator"] = {
    "username": "<operator username>",
    "password": "<operator password>",
    "usergroups": ["<operator usergroups>"],
}

config["tenant_user"] = {
    "username": "<tenant username>",
    "password": "<tenant password>",
    "usergroups": ["<tenant usergroups>"],
}
