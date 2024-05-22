### PART OF MAIN Snakefile ###
### THIS IS ALL Nialco RULES ###
NiAlCo_CHEMSYS = "NiAlCo"
NiAlCo_PROCESS_MODELS = ["Pun", "M3GNet", "CHGNet", "MACE"]
NiAlCo_STRUCTURES = ["L21P"]


# Rule to aggregate the NiAlCo to database
rule aggregate_nialco_db:
    input:
        db="src/data/" + DATABASE,
        create="src/data/COMPLETED_TASKS/created.database.done",
        mineos_calc=expand(
            "src/data/COMPLETED_TASKS/{chemsys}_{structure}_{model}.min_eos.done",
            chemsys=NiAlCo_CHEMSYS,
            structure=NiAlCo_STRUCTURES,
            model=NiAlCo_PROCESS_MODELS,
        ),
        # phonons_calc=expand(
        #     "src/data/COMPLETED_TASKS/{chemsys}_{structure}_{model}.phonons.done",
        #     chemsys=NiAlCo_CHEMSYS,
        #     structure=NiAlCo_STRUCTURES,
        #     model=NiAlCo_PROCESS_MODELS,
        # ),
        elastic_calc=expand(
            "src/data/COMPLETED_TASKS/{chemsys}_{structure}_{model}.elastic.done",
            chemsys=NiAlCo_CHEMSYS,
            structure=NiAlCo_STRUCTURES,
            model=NiAlCo_PROCESS_MODELS,
        ),
    output:
        aggregated="src/data/COMPLETED_TASKS/nialco.database.aggregated.done",
    shell:
        "touch {output.aggregated}"


# Rule for plotting NiAlCo Cohesive Energy
# rule plot_nialco_ecoh:
#     input:
#         data="src/data/" + DATABASE,
#         script="src/scripts/PlotCohesiveEnergy.py",
#     output:
#         figure="src/tex/figures/NiAlCo_CohesiveEnergyPlot.png",
#     conda:
#         "env/base.yml"
#     params:
#         chemsys="NiAlCo",
#     shell:
#         "python src/scripts/PlotCohesiveEnergy.py {DATABASE} {params.chemsys}"
# # Rule for plotting Nialco EOS
# rule plot_nialco_eos:
#     input:
#         data="src/data/CACHED/" + DATABASE,
#         script="src/scripts/PlotEOS.py",
#     output:
#         figure="src/tex/figures/NiAlCo_EquationOfStates.png",
#     conda:
#         "env/base.yml"
#     params:
#         chemsys="NiAlCo",
#     shell:
#         "python src/scripts/PlotEOS.py {DATABASE} {params.chemsys}"
# # Rule for plotting NiAlCo phonons don't plot strained phonons.
# rule plot_nialco_phonons:
#     input:
#         data="src/data/CACHED/" + DATABASE,
#         plotphonons="src/scripts/PlotPhonons.py",
#     #        aggregated="src/data/COMPLETED_TASKS/nialco.database.aggregated.done",
#     output:
#         figure_all_models_B2="src/tex/figures/NiAlCo_L2_1P_ModelsPhononBandstructures.png",
#     conda:
#         "env/base.yml"
#     params:
#         num_strains=5,
#         chemsys="NiAlCo",
#     shell:
#         "python src/scripts/PlotPhonons.py \
#         {DATABASE} \
#         --skip-strains \
#         {params.chemsys} \
#         --models {Nialco_PROCESS_MODELS} \
#         --structures {Nialco_STRUCTURES}"
# Rule for generating NiAlCo equilibrium table
# rule generate_nialco_elastic_table:
#     input:
#         data="src/data/CACHED/" + DATABASE,
#         script="src/scripts/GenerateElasticTable.py",
#     output:
#         table="src/tex/output/Table_NiAlCo_Elastic_Constants.tex",
#     conda:
#         "env/base.yml"
#     params:
#         chemsys="NiAlCo",
#     shell:
#         "python src/scripts/GenerateElasticTable.py {DATABASE} {params.chemsys}"
# # Appendix: Rule for generating NiAlCo equilibrium table
# rule generate_nialco_equil_table:
#     input:
#         data="src/data/CACHED/" + DATABASE,
#         script="src/scripts/GenerateEquilibriumTable.py",
#     output:
#         table="src/tex/output/Table_NiAlCo_Equilibrium_Structures.tex",
#     conda:
#         "env/base.yml"
#     params:
#         chemsys="NiAlCo",
#     shell:
#         "python src/scripts/GenerateEquilibriumTable.py {DATABASE} {params.chemsys}"
