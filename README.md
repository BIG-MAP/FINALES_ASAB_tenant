# FINALES_ASAB_tenant

**A**utonomous **S**ynthesis and **A**nalysis of **B**attery electrolytes Tenant

The tenant interfacing the ASAB system with FINALES.


# Related Documents and Links to FINALES

Documents related to this project and its broader context can be found on the respective Wiki page of this project: [https://github.com/BIG-MAP/FINALES2/wiki/Links](https://github.com/BIG-MAP/FINALES2/wiki/Links)

Links to FINALES:

1. FINALES latest version on Github
[https://github.com/BIG-MAP/FINALES2](https://github.com/BIG-MAP/FINALES2)

1. FINALES used version 1.1.0 on Zenodo
[10.5281/zenodo.10987727](10.5281/zenodo.10987727)

1. Schemas used with FINALES version 1.1.0
[https://github.com/BIG-MAP/FINALES2_schemas](https://github.com/BIG-MAP/FINALES2_schemas)


# Description

The ASAB tenant provides the capability to formulate electrolytes from up to 6 stock solutions, which may either be used in the assembly of coin cells or for analysis of their ionic conductivity. The measurement of ionic conductivity is achieved by means of electrochemical impedance spectroscopy (EIS). This tenant uses the version 1.1.0 of the Modular and Autonomous Data Analysis Platform (MADAP) ([https://doi.org/10.5281/zenodo.8220661](https://doi.org/10.5281/zenodo.8220661), [https://github.com/fuzhanrahmanian/MADAP](https://github.com/fuzhanrahmanian/MADAP)) for the automated analysis of the obtained EIS data.

# Installation

To install the tenant, please follow the steps listed below:

1. Clone this repository
1. Install the packages reported in the requirements.txt
1. Clone the repository of FINALES version 1.1.0 [https://github.com/BIG-MAP/FINALES2](https://github.com/BIG-MAP/FINALES2)
1. Clone the repository of the ASAB laboratory automation package version 1.0.1 [https://github.com/Helge-Stein-Group/ASAB](https://github.com/Helge-Stein-Group/ASAB)
(There is a package called ASAB available from pipy, but this is not the laboratory automation package)
1. Install the FINALES and ASAB packages by switching to the respective directory and
running `pip install . `
1. Adjust the configuration files to match your setup and replace all the placeholders
enclosed in `<>`.
The files to manually alter are:
    - `ASAB_tenant_config_template.py`
    - `configFINALES2_ASAB_tenant_template.py`
    - `edges.csv`
    - `nodes.csv`
    - `stockSolutions.csv`
    - `tubing.csv`
1. Alter the bottom lines of the script `FINALES2_ASAB_tenant.py` to read
    ```
    if __name__ == "__main__":
        ASAB_tenant.tenant_object_to_json()
        # ASAB_tenant.run()
    ```
    and run the script. This will update the `ASAB_tenant_tenant.json` file, which you
    can then send to your FINALES server administrator.
1. Once the server administrator provides you with a tenant UUID, provide it in the
respective key in the instantiation of the ASAB tenant in the `FINALES2_ASAB_tenant.py`
file and alter the bottom lines to read
    ```
    if __name__ == "__main__":
        # ASAB_tenant.tenant_object_to_json()
        ASAB_tenant.run()
    ```

You are now all set to use the tenant once the FINALES instance and all the hardware,
for which it is configured, is running.

# Usage

To start the tenant run the script `FINALES2_ASAB_tenant.py`.

# Acknowledgements

This project received funding from the European Union’s [Horizon 2020 research and innovation programme](https://ec.europa.eu/programmes/horizon2020/en) under grant agreement [No 957189](https://cordis.europa.eu/project/id/957189) (BIG-MAP).
The authors acknowledge BATTERY2030PLUS, funded by the European Union’s Horizon 2020 research and innovation program under grant agreement no. 957213.
This work contributes to the research performed at CELEST (Center for Electrochemical Energy Storage Ulm-Karlsruhe) and was co-funded by the German Research Foundation (DFG) under Project ID 390874152 (POLiS Cluster of Excellence).
