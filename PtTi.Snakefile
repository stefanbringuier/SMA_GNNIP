### PART OF MAIN Snakefile ###
### THIS IS ALL PtTi RULES ###
PtTi_CHEMSYS = "PtTi"
PtTi_PROCESS_MODELS = ["Kim", "M3GNet", "CHGNet", "MACE"]
PtTi_STRUCTURES = ["B2", "B19"]


# Rule to aggregate the PtTi to database
rule aggregate_ptti_db:
    input:
        db="src/data/" + DATABASE,
        create="src/data/COMPLETED_TASKS/created.database.done",
        mineos_calc=expand(
            "src/data/COMPLETED_TASKS/{chemsys}_{structure}_{model}.min_eos.done",
            chemsys=PtTi_CHEMSYS,
            structure=PtTi_STRUCTURES,
            model=PtTi_PROCESS_MODELS,
        ),
        phonons_calc=expand(
            "src/data/COMPLETED_TASKS/{chemsys}_{structure}_{model}.phonons.done",
            chemsys=PtTi_CHEMSYS,
            structure=PtTi_STRUCTURES,
            model=PtTi_PROCESS_MODELS,
        ),
        elastic_calc=expand(
            "src/data/COMPLETED_TASKS/{chemsys}_{structure}_{model}.elastic.done",
            chemsys=PtTi_CHEMSYS,
            structure=PtTi_STRUCTURES,
            model=PtTi_PROCESS_MODELS,
        ),
    output:
        aggregated="src/data/COMPLETED_TASKS/ptti.database.aggregated.done",
    shell:
        "touch {output.aggregated}"


# Rule for plotting PtTi Cohesive Energy
rule plot_ptti_ecoh:
    input:
        data="src/data/" + DATABASE,
        script="src/scripts/PlotCohesiveEnergy.py",
        aggregated="src/data/COMPLETED_TASKS/ptti.database.aggregated.done",
    output:
        figure="src/tex/figures/PtTi_CohesiveEnergyPlot.png",
    conda:
        "env/base.yml"
    params:
        chemsys="PtTi",
    shell:
        "python src/scripts/PlotCohesiveEnergy.py {DATABASE} {params.chemsys}"


# Rule for plotting PtTi EOS
rule plot_ptti_eos:
    input:
        data="src/data/CACHED/" + DATABASE,
        script="src/scripts/PlotEOS.py",
        aggregated="src/data/COMPLETED_TASKS/ptti.database.aggregated.done",
    output:
        figure="src/tex/figures/PtTi_EquationOfStates.png",
    conda:
        "env/base.yml"
    params:
        chemsys="PtTi",
    shell:
        "python src/scripts/PlotEOS.py {DATABASE} {params.chemsys}"


# Rule for plotting PtTi phonons
rule plot_ptti_phonons:
    input:
        data="src/data/CACHED/" + DATABASE,
        plotphonons="src/scripts/PlotPhonons.py",
        aggregated="src/data/COMPLETED_TASKS/ptti.database.aggregated.done",
    output:
        figure_all_models_B2="src/tex/figures/PtTi_B2_ModelsPhononBandstructures.png",
        figure_all_models_B19P="src/tex/figures/PtTi_B19P_ModelsPhononBandstructures.png",
        figure_all_models_B19="src/tex/figures/PtTi_B19_ModelsPhononBandstructures.png",
        figure_all_models_BCO="src/tex/figures/PtTi_BCO_ModelsPhononBandstructures.png",
        figure_strain_mutter_B2="src/tex/figures/PtTi_Mutter_B2_StrainsPhononBandstructures.png",
        figure_strain_zhong_B2="src/tex/figures/PtTi_Zhong_B2_StrainsPhononBandstructures.png",
        figure_strain_ko_B2="src/tex/figures/PtTi_Ko_B2_StrainsPhononBandstructures.png",
        figure_strain_M3GNet_B2="src/tex/figures/PtTi_M3GNet_B2_StrainsPhononBandstructures.png",
        figure_strain_CHGNet_B2="src/tex/figures/PtTi_CHGNet_B2_StrainsPhononBandstructures.png",
        figure_strain_MACE_B2="src/tex/figures/PtTi_MACE_B2_StrainsPhononBandstructures.png",
    conda:
        "env/base.yml"
    params:
        num_strains=5,
        chemsys="PtTi",
    shell:
        "python src/scripts/PlotPhonons.py \
        {DATABASE} \
        {params.num_strains} \
        {params.chemsys} \
        --models {PtTi_PROCESS_MODELS} \
        --structures {PtTi_STRUCTURES}"


# # Rule for generating PtTi M-Mode Gruneisen parameters.
rule generate_ptti_m_mode_gruneisen:
    input:
        data="src/data/CACHED/" + DATABASE,
        script="src/scripts/GruneisenParameters.py",
        aggregated="src/data/COMPLETED_TASKS/ptti.database.aggregated.done",
    output:
        table="src/tex/output/Table_PtTi_M_ModeGruneisen.tex",
        figure="src/tex/figures/Plot_PtTi_M_ModeGruneisen.png",
    conda:
        "env/base.yml"
    params:
        chemsys="PtTi",
        structure="B2",
        qpoint="M",
    shell:
        "python src/scripts/GruneisenParameters.py {DATABASE} {params.chemsys} {params.structure} {params.qpoint}"


# Rule for generating PtTi equilibrium table
rule generate_ptti_elastic_table:
    input:
        data="src/data/CACHED/" + DATABASE,
        script="src/scripts/GenerateElasticTable.py",
        aggregated="src/data/COMPLETED_TASKS/ptti.database.aggregated.done",
    output:
        table="src/tex/output/Table_PtTi_Elastic_Constants.tex",
    conda:
        "env/base.yml"
    params:
        chemsys="PtTi",
    shell:
        "python src/scripts/GenerateElasticTable.py {DATABASE} {params.chemsys}"


# Appendix: Rule for generating PtTi equilibrium table
rule generate_ptti_equil_table:
    input:
        data="src/data/CACHED/" + DATABASE,
        script="src/scripts/GenerateEquilibriumTable.py",
        aggregated="src/data/COMPLETED_TASKS/ptti.database.aggregated.done",
    output:
        table="src/tex/output/Table_PtTi_Equilibrium_Structures.tex",
    conda:
        "env/base.yml"
    params:
        chemsys="PtTi",
    shell:
        "python src/scripts/GenerateEquilibriumTable.py {DATABASE} {params.chemsys}"
