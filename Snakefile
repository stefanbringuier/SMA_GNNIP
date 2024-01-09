# Define common parameters
# NOTE: ALIGNN is supported but results seem total wrong!
NiTi_MODELS = ["Mutter","Zhong","Ko","Kouvasi","DeepMD","M3GNet","CHGNet","MACE","ALIGNN"]
NiTi_STRUCTURES = ["B2","B19","B19P","BCO"]
PtTi_MODELS = ["Kim","M3GNet","CHGNet","MACE","ALIGNN"]
PtTi_STRUCTURES = ["B2","B19"]

#from datetime import datetime
#date_str = datetime.now().strftime("%d%b%Y")
#DATABASE = f"Results_{date_str}.json"
DATABASE = "SMA_Results_07Jan2024.json"

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
    threads: 2
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
        #config="src/scripts/Config.py", # phonon settings, commented because will run everything, better to remove *phonons.done chk file
        #script="src/scripts/CalculatePhonons.py",
        db=ancient("src/data/" + DATABASE),
	create="src/data/COMPLETED_TASKS/created.database.done",
        mineos="src/data/COMPLETED_TASKS/{chemsys}_{structure}_{model}.min_eos.done"
    output:
        done=touch("src/data/COMPLETED_TASKS/{chemsys}_{structure}_{model}.phonons.done")
    threads: 4
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
        #script="src/scripts/CalculateElastic.py",
        db=ancient("src/data/" + DATABASE),
	create="src/data/COMPLETED_TASKS/created.database.done",
        mineos="src/data/COMPLETED_TASKS/{chemsys}_{structure}_{model}.min_eos.done"
    output:
        done=touch("src/data/COMPLETED_TASKS/{chemsys}_{structure}_{model}.elastic.done")
    conda:
        lambda wildcards: get_env_file(wildcards.model)
    threads: 4
    params:
        runner="src/scripts/CalculationRunner.py",
        chemsys=lambda wildcards: wildcards.chemsys,
        structure=lambda wildcards: wildcards.structure,
        model=lambda wildcards: wildcards.model,
    shell:
        "python {params.runner} --calc_type elastic --dbname {DATABASE} --chemsys {params.chemsys} --structure {params.structure} --model {params.model}"

include: "NiTi.Snakefile"

include: "PtTi.Snakefile"



