# Define common parameters
MODELS = ["Mutter","Zhong","M3GNet","CHGNet","MACE","ALIGNN"]
NiTi_STRUCTURES = ["B2","B19","B19P","BCO"] #["BCO","B33", "Pbcm", "B32","R-Phase]
NiTi_DATABASE = "NiTi_Structures.json"

# Function to get environment file
def get_env_file(model):
    return f"env/{model.lower()}.yml"

# Rule to create the NiTi database
rule create_niti_db:
    input:
        create="src/scripts/NewASEDatabase.py",
        script="src/scripts/Generate.py",
        eos="src/scripts/CalculateEOS.py",
        phonons="src/scripts/CalculatePhonons.py"        
    output:
        db="src/data/" + NiTi_DATABASE,
        create="src/data/COMPLETED_TASKS/created.database.done"
    conda:
        "env/ase.yml"
    shell:
        "python src/scripts/NewASEDatabase.py {NiTi_DATABASE}; touch {output.create}"

rule add_niti_to_db:
    input:
        db=ancient("src/data/" + NiTi_DATABASE),
        script="src/scripts/Generate.py",
	create="src/data/COMPLETED_TASKS/created.database.done",
        eos="src/scripts/CalculateEOS.py",
        phonons="src/scripts/CalculatePhonons.py"
    output:
        done=touch("src/data/COMPLETED_TASKS/{structure}_{model}.done")
    conda:
        lambda wildcards: get_env_file(wildcards.model)
    params:
        structure=lambda wildcards: wildcards.structure,
        model=lambda wildcards: wildcards.model,
        nph=13
    shell:
        "python {input.script} --dbname {NiTi_DATABASE} --structure {params.structure} --model {params.model} --nph {params.nph}"

# Rule to aggregate the NiTi database
rule aggregate_niti_db:
    input:
        done=expand("src/data/COMPLETED_TASKS/{structure}_{model}.done", structure=NiTi_STRUCTURES, model=MODELS),
        db="src/data/" + NiTi_DATABASE,
	create="src/data/COMPLETED_TASKS/created.database.done",
        eos="src/scripts/CalculateEOS.py",
        phonons="src/scripts/CalculatePhonons.py"        
    output:
        aggregated="src/data/COMPLETED_TASKS/niti.database.aggregated.done"
    shell:
        "touch {output.aggregated}"

rule cache_niti_db:
    input:
        db="src/data/" + NiTi_DATABASE,
        done="src/data/COMPLETED_TASKS/niti.database.aggregated.done"
    output:
        dbc = "src/data/CACHED/" + NiTi_DATABASE
    cache:
        True
    shell:
        "mkdir -p src/data/CACHED && cp src/data/{NiTi_DATABASE} src/data/CACHED/"
        
# Rule for plotting NiTi EOS
rule plot_niti_eos:
    input:
        script="src/scripts/PlotNiTiEOS.py",
        aggregated="src/data/COMPLETED_TASKS/niti.database.aggregated.done"
    output:
        figure="src/tex/figures/NiTi_EOS_Comparison.png"
    conda:
        "env/ase.yml"
    shell:
        "python src/scripts/PlotNiTiEOS.py {NiTi_DATABASE}"

# Rule for generating NiTi equilibrium table
rule generate_niti_equil_table:
    input:
        script="src/scripts/GenerateNiTiEquilTable.py",
        aggregated="src/data/COMPLETED_TASKS/niti.database.aggregated.done"
    output:
        table="src/tex/output/Table_NiTi_Equilibrium_Structures.tex"
    conda:
        "env/ase.yml"
    shell:
        "python src/scripts/GenerateNiTiEquilTable.py {NiTi_DATABASE}"


# Rule for plotting All models NiTi phonons
rule plot_niti_all_model_phonons:
    input:
        phonons="src/scripts/PlotPhonons.py",
        script="src/scripts/PlotNoStrainPhonons.py",
        aggregated="src/data/COMPLETED_TASKS/niti.database.aggregated.done"
    output:
        figure_all_models_B2="src/tex/figures/B2_combined_models_phonon_bandstructures.png",
        figure_all_models_B19P="src/tex/figures/B19P_combined_models_phonon_bandstructures.png",
        figure_all_models_B19="src/tex/figures/B19_combined_models_phonon_bandstructures.png",
        figure_all_models_BCO="src/tex/figures/BCO_combined_models_phonon_bandstructures.png"
    conda:
        "env/ase.yml"
    shell:
        "python src/scripts/PlotNoStrainPhonons.py {NiTi_DATABASE} {NiTi_STRUCTURES}"

# Rule for plotting NiTi EAM phonons
rule plot_niti_eam_phonons:
    input:
        phonons="src/scripts/PlotPhonons.py",
        script="src/scripts/PlotEAMPhonons.py",
        aggregated="src/data/COMPLETED_TASKS/niti.database.aggregated.done"
    output:
        figure_mutter_B2="src/tex/figures/Mutter_B2_combined_strain_phonon_bandstructures.png",
        figure_zhong_B2="src/tex/figures/Zhong_B2_combined_strain_phonon_bandstructures.png"
    conda:
        "env/ase.yml"
    params:
        npoints=5
    shell:
        "python src/scripts/PlotEAMPhonons.py {NiTi_DATABASE} {params.npoints} {NiTi_STRUCTURES}"

# Rule for plotting NiTi GNN phonons
rule plot_niti_gnn_phonons:
    input:
        phonons="src/scripts/PlotPhonons.py",
        script="src/scripts/PlotGNNPhonons.py",
        aggregated="src/data/COMPLETED_TASKS/niti.database.aggregated.done"
    output:
        figure_m3gnet_B2="src/tex/figures/M3GNet_B2_combined_strain_phonon_bandstructures.png",
#        figure_m3gnet_B19P="src/tex/figures/M3GNet_B19P_combined_strain_phonon_bandstructures.png",
#        figure_m3gnet_B19="src/tex/figures/M3GNet_B19_combined_strain_phonon_bandstructures.png",
#        figure_m3gnet_BCO="src/tex/figures/M3GNet_BCO_combined_strain_phonon_bandstructures.png",
        figure_chgnet_B2="src/tex/figures/CHGNet_B2_combined_strain_phonon_bandstructures.png",
#        figure_chgnet_B19P="src/tex/figures/CHGNet_B19P_combined_strain_phonon_bandstructures.png",
#        figure_chgnet_B19="src/tex/figures/CHGNet_B19_combined_strain_phonon_bandstructures.png",
#        figure_chgnet_BCO="src/tex/figures/CHGNet_BCO_combined_strain_phonon_bandstructures.png",
        figure_mace_B2="src/tex/figures/MACE_B2_combined_strain_phonon_bandstructures.png",
#        figure_mace_B19P="src/tex/figures/MACE_B19P_combined_strain_phonon_bandstructures.png",
#        figure_mace_B19="src/tex/figures/MACE_B19_combined_strain_phonon_bandstructures.png",
#        figure_mace_BCO="src/tex/figures/MACE_BCO_combined_strain_phonon_bandstructures.png",
    conda:
        "env/ase.yml"
    params:
        npoints=5
    shell:
        "python src/scripts/PlotGNNPhonons.py {NiTi_DATABASE} {params.npoints} {NiTi_STRUCTURES}"

# Rule for generating NiTi M-Mode Gruneisen parameters.
rule generate_niti_m_mode_gruneisen:
    input:
        script="src/scripts/GenerateMGruneisenTable.py",
        aggregated="src/data/COMPLETED_TASKS/niti.database.aggregated.done"
    output:
        table="src/tex/output/Table_NiTi_M_ModeGruneisen.tex"
    conda:
        "env/ase.yml"
    shell:
        "python src/scripts/GenerateMGruneisenTable.py {NiTi_DATABASE}"
