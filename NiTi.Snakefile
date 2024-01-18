### PART OF MAIN Snakefile     ###
### THIS IS ALL NiTi RULES     ###
NiTi_CHEMSYS = "NiTi"
NiTi_PROCESS_MODELS = ["Mutter", "Zhong", "Ko", "M3GNet", "CHGNet", "MACE"]


# Rule to aggregate the NiTi to database
rule aggregate_niti_db:
    input:
        db="src/data/" + DATABASE,
        create="src/data/COMPLETED_TASKS/created.database.done",
        #minimize="src/scripts/MinimizeStructure.py",
        #eos="src/scripts/CalculateEOS.py",
        #phonons="src/scripts/CalculatePhonons.py",
        #elastic="src/scripts/CalculateElastic.py",
        mineos_calc=expand(
            "src/data/COMPLETED_TASKS/{chemsys}_{structure}_{model}.min_eos.done",
            chemsys=NiTi_CHEMSYS,
            structure=NiTi_STRUCTURES,
            model=NiTi_PROCESS_MODELS,
        ),
        phonons_calc=expand(
            "src/data/COMPLETED_TASKS/{chemsys}_{structure}_{model}.phonons.done",
            chemsys=NiTi_CHEMSYS,
            structure=NiTi_STRUCTURES,
            model=NiTi_PROCESS_MODELS,
        ),
        elastic_calc=expand(
            "src/data/COMPLETED_TASKS/{chemsys}_{structure}_{model}.elastic.done",
            chemsys=NiTi_CHEMSYS,
            structure=NiTi_STRUCTURES,
            model=NiTi_PROCESS_MODELS,
        ),
    output:
        aggregated="src/data/COMPLETED_TASKS/niti.database.aggregated.done",
    shell:
        "touch {output.aggregated}"


# NOTE: System environment needs to have xvfb-run or remove if display server exist
rule visualize_niti_structures:
    input:
        structures="src/scripts/Structures.py",
        script="src/scripts/VisualizeStructures.py",
    output:
        figure="src/tex/figures/NiTi_VisualizedStructures.png",
    threads: 4
    conda:
        "env/ovito.yml"
    params:
        chemsys="NiTi",
    shell:
        "xvfb-run python src/scripts/VisualizeStructures.py {params.chemsys} {NiTi_STRUCTURES}"


# Rule for plotting NiTi Cohesive Energy
rule plot_niti_ecoh:
    input:
        data="src/data/CACHED/" + DATABASE,
        script="src/scripts/PlotCohesiveEnergy.py",
        aggregated="src/data/COMPLETED_TASKS/niti.database.aggregated.done",
    output:
        figure="src/tex/figures/NiTi_CohesiveEnergyPlot.png",
    threads: 1
    conda:
        "env/ase.yml"
    params:
        chemsys="NiTi",
    shell:
        "python src/scripts/PlotCohesiveEnergy.py {DATABASE} {params.chemsys}"


# Rule for plotting NiTi EOS
rule plot_niti_eos:
    input:
        data="src/data/CACHED/" + DATABASE,
        script="src/scripts/PlotEOS.py",
        aggregated="src/data/COMPLETED_TASKS/niti.database.aggregated.done",
    output:
        figure="src/tex/figures/NiTi_EquationOfStates.png",
    conda:
        "env/ase.yml"
    params:
        chemsys="NiTi",
    shell:
        "python src/scripts/PlotEOS.py {DATABASE} {params.chemsys}"


# Rule for plotting NiTi phonons
rule plot_niti_phonons:
    input:
        data="src/data/CACHED/" + DATABASE,
        plotphonons="src/scripts/PlotPhonons.py",
        aggregated="src/data/COMPLETED_TASKS/niti.database.aggregated.done",
    output:
        figure_all_models_B2="src/tex/figures/NiTi_B2_ModelsPhononBandstructures.png",
        figure_all_models_B19P="src/tex/figures/NiTi_B19P_ModelsPhononBandstructures.png",
        figure_all_models_B19="src/tex/figures/NiTi_B19_ModelsPhononBandstructures.png",
        figure_all_models_BCO="src/tex/figures/NiTi_BCO_ModelsPhononBandstructures.png",
        figure_strain_mutter_B2="src/tex/figures/NiTi_Mutter_B2_StrainsPhononBandstructures.png",
        figure_strain_zhong_B2="src/tex/figures/NiTi_Zhong_B2_StrainsPhononBandstructures.png",
        figure_strain_ko_B2="src/tex/figures/NiTi_Ko_B2_StrainsPhononBandstructures.png",
        figure_strain_M3GNet_B2="src/tex/figures/NiTi_M3GNet_B2_StrainsPhononBandstructures.png",
        figure_strain_CHGNet_B2="src/tex/figures/NiTi_CHGNet_B2_StrainsPhononBandstructures.png",
        figure_strain_MACE_B2="src/tex/figures/NiTi_MACE_B2_StrainsPhononBandstructures.png",
    conda:
        "env/ase.yml"
    params:
        num_strains=5,
        chemsys="NiTi",
    shell:
        "python src/scripts/PlotPhonons.py \
        {DATABASE} \
        {params.num_strains} \
        {params.chemsys} \
        --models {NiTi_PROCESS_MODELS} \
        --structures {NiTi_STRUCTURES}"


# # Rule for generating NiTi M-Mode Gruneisen parameters.
rule generate_niti_m_mode_gruneisen:
    input:
        data="src/data/CACHED/" + DATABASE,
        script="src/scripts/GruneisenParameters.py",
        aggregated="src/data/COMPLETED_TASKS/niti.database.aggregated.done",
    output:
        table="src/tex/output/Table_NiTi_M_ModeGruneisen.tex",
        figure="src/tex/figures/Plot_NiTi_M_ModeGruneisen.png",
    conda:
        "env/ase.yml"
    params:
        chemsys="NiTi",
        structure="B2",
        qpoint="M",
    shell:
        "python src/scripts/GruneisenParameters.py {DATABASE} {params.chemsys} {params.structure} {params.qpoint}"


# Rule for generating NiTi equilibrium table
rule generate_niti_elastic_table:
    input:
        data="src/data/CACHED/" + DATABASE,
        script="src/scripts/GenerateElasticTable.py",
        aggregated="src/data/COMPLETED_TASKS/niti.database.aggregated.done",
    output:
        table="src/tex/output/Table_NiTi_Elastic_Constants.tex",
    conda:
        "env/ase.yml"
    params:
        chemsys="NiTi",
    shell:
        "python src/scripts/GenerateElasticTable.py {DATABASE} {params.chemsys}"


# Appendix: Rule for generating NiTi equilibrium table
rule generate_niti_equil_table:
    input:
        data="src/data/CACHED/" + DATABASE,
        script="src/scripts/GenerateEquilibriumTable.py",
        aggregated="src/data/COMPLETED_TASKS/niti.database.aggregated.done",
    output:
        table="src/tex/output/Table_NiTi_Equilibrium_Structures.tex",
    conda:
        "env/ase.yml"
    params:
        chemsys="NiTi",
    shell:
        "python src/scripts/GenerateEquilibriumTable.py {DATABASE} {params.chemsys}"
