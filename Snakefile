# Define common parameters
# NOTE: ALIGNN is supported but results seem total wrong!
NiTi_MODELS = ["Mutter","Zhong","Ko","M3GNet","CHGNet","MACE","ALIGNN"]
NiTi_STRUCTURES = ["B2","B19","B19P","BCO"]
PtTi_MODELS = ["Kim","M3GNet","CHGNet","MACE"]
PtTi_STRUCTURES = ["B2","B19"]

#from datetime import datetime
#date_str = datetime.now().strftime("%d%b%Y")
#DATABASE = f"Results_{date_str}.json"
DATABASE = f"Results_06Jan2024.json"

# Function to get environment file
def get_env_file(model):
    return f"env/{model.lower()}.yml"

rule create_db:
    input:
        create="src/scripts/NewASEDatabase.py",
    output:
        db="src/data/" + DATABASE,
        create="src/data/COMPLETED_TASKS/created.database.done"
    conda:
        "env/ase.yml"
    shell:
        "python src/scripts/NewASEDatabase.py {DATABASE}; touch {output.create}"

rule minimize_and_calculate_eos:
    input:
        minimize="src/scripts/MinimizeStructure.py",
        script="src/scripts/CalculateEOS.py",
        db=ancient("src/data/" + DATABASE),
	create="src/data/COMPLETED_TASKS/created.database.done"
    output:
        done=touch("src/data/COMPLETED_TASKS/{chemsys}_{structure}_{model}.min_eos.done")
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
        config="src/scripts/Config.py", # phonon settings, this gets changed to converge so commented
        script="src/scripts/CalculatePhonons.py",
        db=ancient("src/data/" + DATABASE),
	create="src/data/COMPLETED_TASKS/created.database.done",
        mineos="src/data/COMPLETED_TASKS/{chemsys}_{structure}_{model}.min_eos.done"
    output:
        done=touch("src/data/COMPLETED_TASKS/{chemsys}_{structure}_{model}.phonons.done")
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
        mineos="src/data/COMPLETED_TASKS/{chemsys}_{structure}_{model}.min_eos.done"
    output:
        done=touch("src/data/COMPLETED_TASKS/{chemsys}_{structure}_{model}.elastic.done")
    conda:
        lambda wildcards: get_env_file(wildcards.model)
    params:
        runner="src/scripts/CalculationRunner.py",
        chemsys=lambda wildcards: wildcards.chemsys,
        structure=lambda wildcards: wildcards.structure,
        model=lambda wildcards: wildcards.model,
    shell:
        "python {params.runner} --calc_type elastic --dbname {DATABASE} --chemsys {params.chemsys} --structure {params.structure} --model {params.model}"

include: "NiTi.Snakefile"

include: "PtTi.Snakefile"
