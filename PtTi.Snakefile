### PART OF MAIN Snakefile ###
### THIS IS ALL PtTi RULES ###
PtTi_CHEMSYS="PtTi"
PtTi_PROCESS_MODELS = ["Kim","M3GNet","CHGNet","MACE"]

# Rule to aggregate the PtTi to database
rule aggregate_ptti_db:
    input:
        db="src/data/" + DATABASE,
	create="src/data/COMPLETED_TASKS/created.database.done",
        #minimize="src/scripts/MinimizeStructure.py",
        #eos="src/scripts/CalculateEOS.py",
        #phonons="src/scripts/CalculatePhonons.py",
        #elastic="src/scripts/CalculateElastic.py",
        mineos_calc=expand("src/data/COMPLETED_TASKS/{chemsys}_{structure}_{model}.min_eos.done", chemsys=PtTi_CHEMSYS, structure=PtTi_STRUCTURES, model=PtTi_PROCESS_MODELS),
        phonons_calc=expand("src/data/COMPLETED_TASKS/{chemsys}_{structure}_{model}.phonons.done", chemsys=PtTi_CHEMSYS, structure=PtTi_STRUCTURES, model=PtTi_PROCESS_MODELS),
        elastic_calc=expand("src/data/COMPLETED_TASKS/{chemsys}_{structure}_{model}.elastic.done", chemsys=PtTi_CHEMSYS, structure=PtTi_STRUCTURES, model=PtTi_PROCESS_MODELS),
    output:
        aggregated="src/data/COMPLETED_TASKS/ptti.database.aggregated.done"
    shell:
        "touch {output.aggregated}"

# Rule for plotting PtTi Cohesive Energy
rule plot_ptti_ecoh:
    input:
        data = "src/data/" + DATABASE,
        script="src/scripts/PlotCohesiveEnergy.py",
        aggregated="src/data/COMPLETED_TASKS/ptti.database.aggregated.done"
    output:
        figure="src/tex/figures/PtTi_CohesiveEnergyPlot.png"
    conda:
        "env/ase.yml"
    params:
        chemsys="PtTi"
    shell:
        "python src/scripts/PlotCohesiveEnergy.py {DATABASE} {params.chemsys}"
