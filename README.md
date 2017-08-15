# phmnl_statistical_worfklow_api
Run the statistical workflow published and maintained by Workflow4Metabolomics within the PhenoMeNal EU project. The workflow processes a Metabolights study (http://www.ebi.ac.uk/metabolights/) and generates univariate and multivariate statistics. For more information please visit to https://github.com/phnmnl/phenomenal-h2020/wiki/Sacurine-statistical-workflow.
The script uses wft4galaxy (http://wft4galaxy.readthedocs.io) as API accessing a Galaxy instance with an installed PhenoMeNal e-infrastructure.

## Description
This script processes all assays and runs the statistics for all factor values found in the investigation file of the given study.

## Requirements
- R >= 3.4.0
- docker (wft4galaxy)
- python (isatab2w4m)

## Preparations
- checkout the reporsitory to your local machine<br>
```bash
git clone https://github.com/c-ruttkies/phmnl_statistical_worfklow_api.git
cd phmnl_statistical_worfklow_api
mkdir output studies
```

- checkout isatab2w4m scripts<br>
```bash
git clone -b develop https://github.com/workflow4metabolomics/mtbls-dwnld
```

- install wft4galaxy-docker (visit http://wft4galaxy.readthedocs.io/installation.html#id2)

- download MTBLS studies of interest from the Metabolights ftp mirror (ftp://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/)
- only the isatab files are needed
```
ls studies/MTBLS404
a_sacurine.txt  audit  i_sacurine.txt  m_sacurine.txt  s_sacurine.txt
```

## Usage
- display help message by simply running the script without arguments
- at the end you can find the example command for processing a single MTBLS study
- make sure to adapt galaxy_url and galaxy_key to point to your PhenoMeNal Galaxy instance

```bash
./scripts/run_statistical_workflow.r 
```

```bash
Usage: ./run_statistical_workflow.r isatab2w4m=... study_path=... study_name=... ga_file_template=... output_path=... log_file=... galaxy_url=... galaxy_key=...


Arguments:

        isatab2w4m - path to isatab2w4m script (https://github.com/workflow4metabolomics/mtbls-dwnld/tree/develop)

        study_path - local path containing MTBLS studies

        study_name - name of the study to be processed (e.g. MTBLS404)

        ga_file_template - Galaxy workflow file for the statistical workflow

        output_path - local path used to store the result files

        log_file - path to csv file used for log messages

        galaxy_url - url to Galaxy server with installed PhenoMeNal on which to run the workflow

        galaxy_key - API key to access the Galaxy server


Description:

        Runs the statistical worfklow (https://github.com/phnmnl/phenomenal-h2020/wiki/Sacurine-statistical-workflow) using wft4galaxy.
        Make sure that wft4galaxy-docker is available on your machine. Visit http://wft4galaxy.readthedocs.io/installation.html#id2 for installation instructions.
        The R package Risa needs to be installed. Visit https://bioconductor.org/packages/release/bioc/html/Risa.html for installation instructions.


Example: 

        ./scripts/run_statistical_workflow.r isatab2w4m="mtbls-dwnld/isatab2w4m" \
                study_path="studies" \
                study_name="MTBLS404" \
                ga_file_template="template/w4m-sacurine-statistics.ga" \
                output_path="output" \
                log_file="MTBLS404_logs.csv" \
                galaxy_url="http://127.0.0.1:30700" \
                galaxy_key="ce30966564d6a42b22f951ca023081ed"
```
