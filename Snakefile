# Define common parameters
# NOTE: ALIGNN is supported but results seem total wrong!
# from datetime import datetime
# date_str = datetime.now().strftime("%d%b%Y")
# DATABASE = f"Results_{date_str}.json"
DATABASE = "SMA_Results.json"


# Function to get environment file
def get_env_file(model):
    env = f"env/{model.lower()}.yml"
    if model == "M3GNet":
        return env
    elif model == "CHGNet":
        return env
    else:
        return "env/base.yml"


rule create_db:
    input:
        create="src/scripts/NewASEDatabase.py",
    output:
        db="src/data/" + DATABASE,
        create="src/data/COMPLETED_TASKS/created.database.done",
    conda:
        "env/base.yml"
    shell:
        "python src/scripts/NewASEDatabase.py {DATABASE}; touch {output.create}"


rule minimize_and_calculate_eos:
    input:
        minimize="src/scripts/MinimizeStructure.py",
        script="src/scripts/CalculateEOS.py",
        db=ancient("src/data/" + DATABASE),
        create="src/data/COMPLETED_TASKS/created.database.done",
    output:
        done=touch(
            "src/data/COMPLETED_TASKS/{chemsys}_{structure}_{model}.min_eos.done"
        ),
    threads: 3
    resources:
        mem_gb=7,
    conda:
        lambda wildcards: get_env_file(wildcards.model)
    params:
        runner="src/scripts/CalculationRunner.py",
        chemsys=lambda wildcards: wildcards.chemsys,
        structure=lambda wildcards: wildcards.structure,
        model=lambda wildcards: wildcards.model,
    shell:
        "python {params.runner} --calc_type min_eos --dbname {DATABASE} --chemsys {params.chemsys} --structure {params.structure} --model {params.model}"


rule calculate_phonons:
    input:
        config="src/scripts/Config.py",
        script="src/scripts/CalculatePhonons.py",
        db=ancient("src/data/" + DATABASE),
        create="src/data/COMPLETED_TASKS/created.database.done",
        mineos="src/data/COMPLETED_TASKS/{chemsys}_{structure}_{model}.min_eos.done",
    output:
        done=touch(
            "src/data/COMPLETED_TASKS/{chemsys}_{structure}_{model}.phonons.done"
        ),
    threads: 12
    resources:
        mem_gb=30,
    conda:
        lambda wildcards: get_env_file(wildcards.model)
    params:
        runner="src/scripts/CalculationRunner.py",
        chemsys=lambda wildcards: wildcards.chemsys,
        structure=lambda wildcards: wildcards.structure,
        model=lambda wildcards: wildcards.model,
    shell:
        "python {params.runner} --calc_type phonons --dbname {DATABASE} --chemsys {params.chemsys} --structure {params.structure} --model {params.model}"


rule calculate_elastic:
    input:
        script="src/scripts/CalculateElastic.py",
        db=ancient("src/data/" + DATABASE),
        create="src/data/COMPLETED_TASKS/created.database.done",
        mineos="src/data/COMPLETED_TASKS/{chemsys}_{structure}_{model}.min_eos.done",
    output:
        done=touch(
            "src/data/COMPLETED_TASKS/{chemsys}_{structure}_{model}.elastic.done"
        ),
    threads: 1
    resources:
        mem_gb=15,
    conda:
        lambda wildcards: get_env_file(wildcards.model)
    params:
        runner="src/scripts/CalculationRunner.py",
        chemsys=lambda wildcards: wildcards.chemsys,
        structure=lambda wildcards: wildcards.structure,
        model=lambda wildcards: wildcards.model,
    shell:
        "python {params.runner} --calc_type elastic --dbname {DATABASE} --chemsys {params.chemsys} --structure {params.structure} --model {params.model}"


import os


# NOTE: We use this to bypass running if we change a script in a rule above
# but don't want to rerun if all the cached database and checkpoint files
# are present for the cache_db rul. This means if you want to tigger the
# calculation rules you need to delete the checkpoint files.
def is_cache_db_ready():
    dependencies = [
        "src/data/CACHED/" + DATABASE,
        "src/data/COMPLETED_TASKS/niti.database.aggregated.done",
        "src/data/COMPLETED_TASKS/ptti.database.aggregated.done",
        "src/data/COMPLETED_TASKS/nialco.database.aggregated.done",
    ]
    return all(os.path.exists(dep) for dep in dependencies)


if not is_cache_db_ready():

    include: "NiTi.Snakefile"
    include: "PtTi.Snakefile"
    include: "NiAlCo.Snakefile"


# NOTE: This is the final database after all calculations/simulations.
# All the plots that leverage the ASE JSON database use this file path.
# Also Zenodo caching to sandbox seems broken. I can create a record but
# and showyourwork validates but when it tries to create a draft it fails.
rule cache_db:
    input:
        db="src/data/" + DATABASE,
        niti_done="src/data/COMPLETED_TASKS/niti.database.aggregated.done",
        ptti_done="src/data/COMPLETED_TASKS/ptti.database.aggregated.done",
        alnico_done="src/data/COMPLETED_TASKS/nialco.database.aggregated.done",
    output:
        dbc="src/data/CACHED/" + DATABASE,
    #    cache: True
    shell:
        "mkdir -p src/data/CACHED && cp src/data/{DATABASE} src/data/CACHED/"


rule generate_niti_bz_appendix:
    input:
        script="src/scripts/AppendixBZ.py",
    output:
        figure_ibz_B2="src/tex/figures/B2_BrillouinZonePointsSampled.png",
        figure_ibz_B19="src/tex/figures/B19_BrillouinZonePointsSampled.png",
        figure_ibz_B19P="src/tex/figures/B19P_BrillouinZonePointsSampled.png",
        figure_ibz_BCO="src/tex/figures/BCO_BrillouinZonePointsSampled.png",
        table_qpoints_B2="src/tex/output/B2_SpecialSymmetryPointsBZ.tex",
        table_qpoints_B19="src/tex/output/B19_SpecialSymmetryPointsBZ.tex",
        table_qpoints_B19P="src/tex/output/B19P_SpecialSymmetryPointsBZ.tex",
        table_qpoints_BCO="src/tex/output/BCO_SpecialSymmetryPointsBZ.tex",
    conda:
        "env/base.yml"
    script:
        "src/scripts/AppendixBZ.py"


rule show_phonon_config_appendix:
    input:
        config="src/scripts/Config.py",
        script="src/scripts/WriteCodeToTex.py",
    output:
        code_listing="src/tex/output/Phonons_ConfigSettingsCode.tex",
    params:
        code="Config.py",
        funcname="get_phonon_config",
        output="Phonons_ConfigSettingsCode.tex",
    shell:
        "python src/scripts/WriteCodeToTex.py {params.code} {params.funcname} {params.output}"
