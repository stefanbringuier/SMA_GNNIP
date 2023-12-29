### PART OF MAIN Snakefile ###
### THIS IS ALL PtTi RULES ###

# Rule to aggregate the NiTi to database
rule aggregate_ptti_db:
    input:
        db="src/data/" + DATABASE,
	create="src/data/COMPLETED_TASKS/created.database.done",
        minimize="src/scripts/MinimizeStructure.py",
        eos="src/scripts/CalculateEOS.py",
        phonons="src/scripts/CalculatePhonons.py",
        elastic="src/scripts/CalculateElastic.py",
        mineos_calc=expand("src/data/COMPLETED_TASKS/{structure}_{model}.min_eos.done", chemsys="PtTi", structure=PtTi_STRUCTURES, model=PtTi_MODELS),
        phonons_calc=expand("src/data/COMPLETED_TASKS/{structure}_{model}.phonons.done", chemsys="PtTi", structure=PtTi_STRUCTURES, model=PtTi_MODELS),
        elastic_calc=expand("src/data/COMPLETED_TASKS/{structure}_{model}.elastic.done", chemsys="PtTi", structure=PtTi_STRUCTURES, model=PtTi_MODELS),
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
