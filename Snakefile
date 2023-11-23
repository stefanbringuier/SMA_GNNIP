# Define common parameters
MODELS = ["Mutter", "Zhong","M3GNet","CHGNet","MACE","ALIGNN"]
NiTi_STRUCTURES = ["B2"] #["B2","B19","B19P","B33", "Pbcm", "B32"]
NiTi_DATABASE = "NiTi_Structures.json"

# Function to get environment file
def get_env_file(model):
    return f"env/{model.lower()}.yml"

# Rule to create the NiTi database
rule create_niti_db:
    input:
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
    cache:
        True
    shell:
        "touch {output.aggregated}"

# Rule for plotting NiTi EOS
rule plot_niti_eos:
    input:
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
        figure_all_models_b2="src/tex/figures/B2_combined_models_phonon_bandstructures.png"
    conda:
        "env/ase.yml"
    shell:
        "python src/scripts/PlotNoStrainPhonons.py {NiTi_DATABASE} {NiTi_STRUCTURES}"

# Rule for plotting NiTi EAM phonons
rule plot_niti_eam_phonons:
    input:
        aggregated="src/data/COMPLETED_TASKS/niti.database.aggregated.done"
    output:
        figure_mutter_b2="src/tex/figures/Mutter_B2_combined_strain_phonon_bandstructures.png",
        figure_zhong_b2="src/tex/figures/Zhong_B2_combined_strain_phonon_bandstructures.png"
    conda:
        "env/ase.yml"
    params:
        npoints=5
    shell:
        "python src/scripts/PlotEAMPhonons.py {NiTi_DATABASE} {params.npoints} {NiTi_STRUCTURES}"

# Rule for plotting NiTi GNN phonons
rule plot_niti_gnn_phonons:
    input:
        aggregated="src/data/COMPLETED_TASKS/niti.database.aggregated.done"
    output:
        figure_m3gnet_b2="src/tex/figures/M3GNet_B2_combined_strain_phonon_bandstructures.png",
        figure_chgnet_b2="src/tex/figures/CHGNet_B2_combined_strain_phonon_bandstructures.png",
        figure_mace_b2="src/tex/figures/MACE_B2_combined_strain_phonon_bandstructures.png",
    conda:
        "env/ase.yml"
    params:
        npoints=5
    shell:
        "python src/scripts/PlotGNNPhonons.py {NiTi_DATABASE} {params.npoints} {NiTi_STRUCTURES}"
