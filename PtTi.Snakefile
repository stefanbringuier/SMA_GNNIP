### PART OF MAIN Snakefile ###
### THIS IS ALL PtTi RULES ###
PtTi_CHEMSYS="PtTi"
PtTi_PROCESS_MODELS = ["Kim","M3GNet","CHGNet","MACE"]

# Rule to aggregate the PtTi to database
rule aggregate_ptti_db:
    input:
        db="src/data/" + DATABASE,
	create="src/data/COMPLETED_TASKS/created.database.done",
        minimize="src/scripts/MinimizeStructure.py",
        eos="src/scripts/CalculateEOS.py",
        phonons="src/scripts/CalculatePhonons.py",
        elastic="src/scripts/CalculateElastic.py",
        mineos_calc=expand("src/data/COMPLETED_TASKS/{chemsys}_{structure}_{model}.min_eos.done", chemsys=PtTi_CHEMSYS, structure=PtTi_STRUCTURES, model=PtTi_PROCESS_MODELS),
        phonons_calc=expand("src/data/COMPLETED_TASKS/{chemsys}_{structure}_{model}.phonons.done", chemsys=PtTi_CHEMSYS, structure=PtTi_STRUCTURES, model=PtTi_PROCESS_MODELS),
        elastic_calc=expand("src/data/COMPLETED_TASKS/{chemsys}_{structure}_{model}.elastic.done", chemsys=PtTi_CHEMSYS, structure=PtTi_STRUCTURES, model=PtTi_PROCESS_MODELS),
    output:
        aggregated="src/data/COMPLETED_TASKS/ptti.database.aggregated.done"
    shell:
        "touch {output.aggregated}"

# NOTE: My intent is to cache the database, but I think this is not needed
# I believe the purpose behind caching is to store intermediate results between
# rules/workflows. For example if you need to fetch a database and had rules
# that kept doing this, then you could cache that fetch. I think!.
# What I was tyring to do was create the action of storing the databse once its
# created because I'm frequentyly changeing whats being added to the database.
rule cache_ptti_db:
    input:
        db="src/data/" + DATABASE,
        done="src/data/COMPLETED_TASKS/ptti.database.aggregated.done"
    output:
        dbc = "src/data/CACHED/" + DATABASE
    cache:
        True
    shell:
        "mkdir -p src/data/CACHED && cp src/data/{DATABASE} src/data/CACHED/"

# Rule for plotting PtTi Cohesive Energy
rule plot_ptti_ecoh:
    input:
        data = "src/data/" + DATABASE,
        script="src/scripts/PlotCohesiveEnergy.py",
        aggregated="src/data/COMPLETED_TASKS/ptti.database.aggregated.done"
    output:
        figure="src/tex/figures/PtTi_CohesiveEnergyPlot.png"
    threads:
        1
    conda:
        "env/ase.yml"
    params:
        chemsys="PtTi"
    shell:
        "python src/scripts/PlotCohesiveEnergy.py {DATABASE} {params.chemsys}"
