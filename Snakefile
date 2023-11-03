MODELS = ["Zhong", "Mutter", "M3GNet", "CHGNet", "MACE", "ALIGNN"]
NiTi_STRUCTURES = ["B2", "B19P", "B19", "B33", "Pbcm", "B32", "R_Phase"]
NiTi_EOS_DATABASE = "NiTi_EOS.json"
NiTi_PHONON_DATABASE = "NiTi_Phonons.json"

def get_env_file(model):
    return f"env/{model.lower()}.yml"

rule run_all_niti_eos:
    input:
        "src/data/" + NiTi_EOS_DATABASE,
        "src/tex/figures/NiTi_EOS_Comparison.png"
    cache:
        True

rule add_niti_to_eos_db:
    input:
        "src/scripts/GenerateEOS.py"  # corrected the path
    output:
        touch("src/data/COMPLETED_TASKS/{structure}_{model}.eos.done")
    conda:
        lambda wildcards: get_env_file(wildcards.model)  # use a lambda to pass the wildcards to get_env_file
    params:
        dbname=NiTi_EOS_DATABASE,
        structure=lambda wildcards: wildcards.structure,
        model=lambda wildcards: wildcards.model
    shell:
        "python src/scripts/GenerateEOS.py --dbname {params.dbname} --structure {params.structure} --model {params.model}"

rule aggregate_niti_eos_db:
    input:
        expand("src/data/COMPLETED_TASKS/{structure}_{model}.eos.done", structure=NiTi_STRUCTURES, model=MODELS)
    output:
        touch("src/data/COMPLETED_TASKS/eos.database.done")

rule add_metadata_to_niti_eos_db:
    input:
        "src/scripts/SetMetadata.py",
	"src/data/" + NiTi_EOS_DATABASE,
	"src/data/COMPLETED_TASKS/eos.database.done"
    output:
        touch("src/data/COMPLETED_TASKS/eos.metadata.done")
    conda:
        "env/ase.yml"
    script:
        "src/scripts/SetMetadata.py"
    
rule plot_niti_eos:
    input:
        "src/scripts/GenerateEOS.py",
        "src/scripts/PlotNiTiEOS.py",
        "src/data/COMPLETED_TASKS/eos.database.done",
    output:
        "src/tex/figures/NiTi_EOS_Comparison.png"
    conda:
        "env/ase.yml"
    script:
        "src/scripts/PlotNiTiEOS.py"