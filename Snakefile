# Define common parameters
# NOTE: ALIGNN is supported but seems total wrong results!
NiTi_MODELS = ["Mutter","Zhong","Ko","M3GNet","CHGNet","MACE"]
NiTi_STRUCTURES = ["B2","B19","B19P","BCO"]
PtTi_MODELS = [""]
PtTi_STRUCTURES = [""]

DATABASE = "Results.json"

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

# Original rule that does everything
# rule calculate:
#     input:
#         db=ancient("src/data/" + DATABASE),
#         script="src/scripts/Generate.py",
# 	create="src/data/COMPLETED_TASKS/created.database.done",
#         minimize="src/scripts/MinimizeStructure.py",
#         eos="src/scripts/CalculateEOS.py",
#         #config="src/scripts/Config.py", # phonon settings, this gets changed to converge so commented
#         phonons="src/scripts/CalculatePhonons.py",
#         elastic="src/scripts/CalculateElastic.py"
#     output:
#         done=touch("src/data/COMPLETED_TASKS/{structure}_{model}.done")
#     conda:
#         lambda wildcards: get_env_file(wildcards.model)
#     params:
#         structure=lambda wildcards: wildcards.structure,
#         model=lambda wildcards: wildcards.model,
#         n_strain_phonons=13
#     shell:
#         "python {input.script} --dbname {DATABASE} --structure {params.structure} --model {params.model} --nph {params.n_strain_phonons}"

rule minimize_and_calculate_eos:
    input:
        minimize="src/scripts/MinimizeStructure.py",
        script="src/scripts/CalculateEOS.py",
        db=ancient("src/data/" + DATABASE),
	create="src/data/COMPLETED_TASKS/created.database.done"
    output:
        done=touch("src/data/COMPLETED_TASKS/{structure}_{model}.min_eos.done")
    conda:
        lambda wildcards: get_env_file(wildcards.model)
    params:
        runner="src/scripts/CalculationRunner.py",
        structure=lambda wildcards: wildcards.structure,
        model=lambda wildcards: wildcards.model,
    shell:
        "python {params.runner} --calc_type min_eos --dbname {DATABASE} --structure {params.structure} --model {params.model}"

rule calculate_phonons:
    input:
        config="src/scripts/Config.py", # phonon settings, this gets changed to converge so commented
        script="src/scripts/CalculatePhonons.py",
        db=ancient("src/data/" + DATABASE),
	create="src/data/COMPLETED_TASKS/created.database.done",
        mineos="src/data/COMPLETED_TASKS/{structure}_{model}.min_eos.done"
    output:
        done=touch("src/data/COMPLETED_TASKS/{structure}_{model}.phonons.done")
    conda:
        lambda wildcards: get_env_file(wildcards.model)
    params:
        runner="src/scripts/CalculationRunner.py",
        structure=lambda wildcards: wildcards.structure,
        model=lambda wildcards: wildcards.model,
    shell:
        "python {params.runner} --calc_type phonons --dbname {DATABASE} --structure {params.structure} --model {params.model}"

rule calculate_elastic:
    input:
        script="src/scripts/CalculateElastic.py",
        db=ancient("src/data/" + DATABASE),
	create="src/data/COMPLETED_TASKS/created.database.done",
        mineos="src/data/COMPLETED_TASKS/{structure}_{model}.min_eos.done"
    output:
        done=touch("src/data/COMPLETED_TASKS/{structure}_{model}.elastic.done")
    conda:
        lambda wildcards: get_env_file(wildcards.model)
    params:
        runner="src/scripts/CalculationRunner.py",
        structure=lambda wildcards: wildcards.structure,
        model=lambda wildcards: wildcards.model,
    shell:
        "python {params.runner} --calc_type elastic --dbname {DATABASE} --structure {params.structure} --model {params.model}"

### THIS IS ALL NiTi Rules ###
           
# Rule to aggregate the NiTi to database
rule aggregate_niti_db:
    input:
        db="src/data/" + DATABASE,
	create="src/data/COMPLETED_TASKS/created.database.done",
        minimize="src/scripts/MinimizeStructure.py",
        eos="src/scripts/CalculateEOS.py",
        phonons="src/scripts/CalculatePhonons.py",
        elastic="src/scripts/CalculateElastic.py",
        mineos_calc=expand("src/data/COMPLETED_TASKS/{structure}_{model}.min_eos.done", structure=NiTi_STRUCTURES, model=NiTi_MODELS),
        phonons_calc=expand("src/data/COMPLETED_TASKS/{structure}_{model}.phonons.done", structure=NiTi_STRUCTURES, model=NiTi_MODELS),
        elastic_calc=expand("src/data/COMPLETED_TASKS/{structure}_{model}.elastic.done", structure=NiTi_STRUCTURES, model=NiTi_MODELS),
    output:
        aggregated="src/data/COMPLETED_TASKS/niti.database.aggregated.done"
    shell:
        "touch {output.aggregated}"

# NOTE: My intent is to cache the database, but I think this is not needed
# I believe the purpose behind caching is to store intermediate results between
# rules/workflows. For example if you need to fetch a database and had rules
# that kept doing this, then you could cache that fetch. I think!.
# What I was tyring to do was create the action of storing the databse once its
# created because I'm frequentyly changeing whats being added to the database.
rule cache_niti_db:
    input:
        db="src/data/" + DATABASE,
        done="src/data/COMPLETED_TASKS/niti.database.aggregated.done"
    output:
        dbc = "src/data/CACHED/" + DATABASE
    cache:
        True
    shell:
        "mkdir -p src/data/CACHED && cp src/data/{DATABASE} src/data/CACHED/"


# Rule for plotting NiTi Cohesive Energy
rule plot_niti_ecoh:
    input:
        data = "src/data/" + DATABASE,
        script="src/scripts/PlotCohesiveEnergy.py",
        aggregated="src/data/COMPLETED_TASKS/niti.database.aggregated.done"
    output:
        figure="src/tex/figures/NiTi_CohesiveEnergyPlot.png"
    conda:
        "env/ase.yml"
    params:
        chemsys="NiTi"
    shell:
        "python src/scripts/PlotCohesiveEnergy.py {DATABASE} {params.chemsys}"

# Rule for plotting NiTi EOS
rule plot_niti_eos:
    input:
        data = "src/data/" + DATABASE,
        script="src/scripts/PlotEOS.py",
        aggregated="src/data/COMPLETED_TASKS/niti.database.aggregated.done"
    output:
        figure="src/tex/figures/NiTi_EquationOfStates.png"
    conda:
        "env/ase.yml"
    params:
        chemsys="NiTi"
    shell:
        "python src/scripts/PlotEOS.py {DATABASE} {params.chemsys}"

# Rule for generating NiTi equilibrium table
rule generate_niti_equil_table:
    input:
        script="src/scripts/GenerateEquilibriumTable.py",
        aggregated="src/data/COMPLETED_TASKS/niti.database.aggregated.done"
    output:
        table="src/tex/output/Table_NiTi_Equilibrium_Structures.tex"
    conda:
        "env/ase.yml"
    params:
        chemsys="NiTi"
    shell:
        "python src/scripts/GenerateEquilibriumTable.py {DATABASE} {params.chemsys}"


# Rule for plotting NiTi phonons
rule plot_niti_phonons:
    input:
        plotphonons="src/scripts/PlotPhonons.py",
        aggregated="src/data/COMPLETED_TASKS/niti.database.aggregated.done"
    output:
        figure_all_models_B2="src/tex/figures/NiTi_B2_ModelsPhononBandstructures.png",
        figure_all_models_B19P="src/tex/figures/NiTi_B19P_ModelsPhononBandstructures.png",
        figure_all_models_B19="src/tex/figures/NiTi_B19_ModelsPhononBandstructures.png",
        figure_all_models_BCO="src/tex/figures/NiTi_BCO_ModelsPhononBandstructures.png",
        figure_strain_mutter_B2="src/tex/figures/NiTi_Mutter_B2_StrainsPhononBandstructures.png",
        figure_strain_zhong_B2="src/tex/figures/NiTi_Zhong_B2_StrainsPhononBandstructures.png",
        figure_strain_M3GNet_B2="src/tex/figures/NiTi_M3GNet_B2_StrainsPhononBandstructures.png",
        figure_strain_CHGNet_B2="src/tex/figures/NiTi_CHGNet_B2_StrainsPhononBandstructures.png",
        figure_strain_MACE_B2="src/tex/figures/NiTi_MACE_B2_StrainsPhononBandstructures.png",
    conda:
        "env/ase.yml"
    params:
        num_strains=6,
        chemsys="NiTi"
    shell:
        "python src/scripts/PlotPhonons.py \
        {DATABASE} \
        {params.num_strains} \
        {params.chemsys} \
        --models {NiTi_MODELS} \
        --structures {NiTi_STRUCTURES}"

# # Rule for generating NiTi M-Mode Gruneisen parameters.
# rule generate_niti_m_mode_gruneisen:
#     input:
#         script="src/scripts/GenerateMGruneisenTable.py",
#         aggregated="src/data/COMPLETED_TASKS/niti.database.aggregated.done"
#     output:
#         table="src/tex/output/Table_NiTi_M_ModeGruneisen.tex"
#     conda:
#         "env/ase.yml"
#     shell:
#         "python src/scripts/GenerateMGruneisenTable.py {DATABASE}"

rule generate_niti_bz_appendix:
    input:
        script="src/scripts/AppendixBZ.py",
    output:
        figure_ibz_B2 = "src/tex/figures/B2_BrillouinZonePointsSampled.png",
        figure_ibz_B19 = "src/tex/figures/B19_BrillouinZonePointsSampled.png",
        figure_ibz_B19P = "src/tex/figures/B19P_BrillouinZonePointsSampled.png",
        figure_ibz_BCO = "src/tex/figures/BCO_BrillouinZonePointsSampled.png",
        table_qpoints_B2 = "src/tex/output/B2_SpecialSymmetryPointsBZ.tex",
        table_qpoints_B19 ="src/tex/output/B19_SpecialSymmetryPointsBZ.tex",
        table_qpoints_B19P ="src/tex/output/B19P_SpecialSymmetryPointsBZ.tex",
        table_qpoints_BCO ="src/tex/output/BCO_SpecialSymmetryPointsBZ.tex",
    conda:
        "env/ase.yml"
    script:
        "src/scripts/AppendixBZ.py"
