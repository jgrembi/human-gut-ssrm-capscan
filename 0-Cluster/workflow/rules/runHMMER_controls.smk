
rule runProdigal:
	input: join(config["coAssemblyControlDir"],"{control}/megahit/filtered_contigs_1kb.fasta")
	output:
		geneCoords=join(config["prodigalGeneCoordDir"],"{control}_geneCoord.out"),
		proteinSeqs=join(config["prodigalProteinSeqDir"], "{control}.faa")
	conda:
		"../envs/prodigal.yml"
	shell:
		"""
		prodigal -i {input} -o {output.geneCoords} -a {output.proteinSeqs}
		"""

rule runHMMER:
	input:
		join(config["prodigalProteinSeqDir"],"{control}.faa")
	output:
		hmmOut_bact = join(config["hmmOutDir"],"{control}_dsrAB_bact.hmm.out"),
		domOut_bact = join(config["domtblDir"],"{control}_dsrAB_bact.domtblout"),
		hmmOut_arch = join(config["hmmOutDir"],"{control}_dsrAB_arch.hmm.out"),
		domOut_arch = join(config["domtblDir"],"{control}_dsrAB_arch.domtblout"),
		msa_bact= join(config["msaDir"],"sto/{control}_dsrAB_bact.sto"),
		msa_arch= join(config["msaDir"],"sto/{control}_dsrAB_arch.sto")
	params:
		hmm_profile_bact=config["dsrAB_bact_HMMProfile"],
		hmm_profile_arch=config["dsrAB_arch_HMMProfile"]
	conda:
		"../envs/hmmer.yml"
	shell:
		"""
		hmmsearch -o {output.hmmOut_bact} --domtblout {output.domOut_bact} -A {output.msa_bact} {params.hmm_profile_bact} {input}
		hmmsearch -o {output.hmmOut_arch} --domtblout {output.domOut_arch} -A {output.msa_arch} {params.hmm_profile_arch} {input}
		"""