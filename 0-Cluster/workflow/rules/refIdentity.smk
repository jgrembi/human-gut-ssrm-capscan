# pull fasta entries for the reference leaves
rule pullRefLeaves:
  input:
		OTUs_to_use = join(config["refIdentityDir"],"OTUs_to_ref.txt"),
		refFASTA = config["refMSA"]
	output:
	  join(config["refIdentityDir"],"OTUs_to_ref.faa")
	resources:
		mem_mb=4000,
		time=60
	conda:
		"../envs/biopython.yml"
	params:
		scripts_dir=config["scriptsDir"]
	shell:
		"""
		python3 {params.scripts_dir}/selectOTUSeqsInCluster.py {input.refFASTA} {input.OTUs_to_use} {output}
		"""

# two exact hit duplicates that mess up the pullQuerySeqs rule. 
rule removeTwoDups:
	input:
		join(config["cleanHitsDir"],"compiled_dsrAB_hits_scoreThreshold_assemblies.faa")
	output:
		join(config["refIdentityDir"],"compiled_dsrAB_hits_scoreThreshold_assemblies_filtered.faa")
	resources:
		mem_mb=4000,
		time=60
	conda:
		"../envs/biopython.yml"
	shell:
		"""
		awk '
		BEGIN {{ drop[">ORJZ01.17253.1_1/2-182"]; drop[">OWAH01.53685.1_1/1-179"]; drop[">OQPP01.20025.1_1/64-390"]  }}
		/^>/ {{ skip = ($0 in drop); if (skip) delete drop[$0] }}
		!skip
		' {input} > {output}
		"""

rule pullQuerySeqs:
	input:
		queries = join(config["refIdentityDir"],"hits_all.txt"),
		fasta = join(config["refIdentityDir"],"compiled_dsrAB_hits_scoreThreshold_assemblies_filtered.faa")
	output:
		join(config["refIdentityDir"],"hits_all.faa")
	resources:
		mem_mb=4000,
		time=60
	conda:
		"../envs/biopython.yml"
	params:
		scripts_dir=config["scriptsDir"]
	shell:
		"""
		python3 {params.scripts_dir}/selectOTUSeqsInCluster.py {input.fasta} {input.queries} {output}
		"""

rule pullQuerySeqs_capsule:
	input:
		queries = join(config["refIdentityDir"],"hits_capsule.txt"),
		fasta = join(config["cleanHitsDir"],"StoolCapsule_arch_bact_compiled_dsrAB_hits_scoreThreshold.faa")
	output:
		join(config["refIdentityDir"],"hits_capsule.faa")
	resources:
		mem_mb=4000,
		time=60
	conda:
		"../envs/biopython.yml"
	params:
		scripts_dir=config["scriptsDir"]
	shell:
		"""
		python3 {params.scripts_dir}/selectOTUSeqsInCluster.py {input.fasta} {input.queries} {output}
		"""

rule makeBlastDatabase:
	input: join(config["refIdentityDir"],"OTUs_to_ref.faa")
	output:
		join(config["refIdentityDir"],"db/OTUs_to_ref.pdb"),
		join(config["refIdentityDir"],"db/OTUs_to_ref.phr"),
		join(config["refIdentityDir"],"db/OTUs_to_ref.pin"),
		join(config["refIdentityDir"],"db/OTUs_to_ref.pjs"),
		join(config["refIdentityDir"],"db/OTUs_to_ref.pot"),
		join(config["refIdentityDir"],"db/OTUs_to_ref.psq"),
		join(config["refIdentityDir"],"db/OTUs_to_ref.ptf"),
		join(config["refIdentityDir"],"db/OTUs_to_ref.pto")
	resources:
		mem_mb=4000,
		time=60
	params:
		join(config["refIdentityDir"],"db/OTUs_to_ref")
	conda:
		"../envs/blast.yml"
	shell:
		"""
		makeblastdb -in {input} -out {params} -dbtype prot
		"""

rule blastHitsToOTUs:
	input:
		query = join(config["refIdentityDir"],"hits_all.faa"),
		pdb = join(config["refIdentityDir"],"db/OTUs_to_ref.pdb"),
		phr = join(config["refIdentityDir"],"db/OTUs_to_ref.phr"),
		pin = join(config["refIdentityDir"],"db/OTUs_to_ref.pin"),
		pjs = join(config["refIdentityDir"],"db/OTUs_to_ref.pjs"),
		pot = join(config["refIdentityDir"],"db/OTUs_to_ref.pot"),
		psq = join(config["refIdentityDir"],"db/OTUs_to_ref.psq"),
		ptf = join(config["refIdentityDir"],"db/OTUs_to_ref.ptf"),
		pto = join(config["refIdentityDir"],"db/OTUs_to_ref.pto")
	output: join(config["refIdentityDir"],"out/NCBI_blast.txt")
	resources:
		mem_mb=4000,
		time=240
	threads: 8
	conda:
		"../envs/blast.yml"
	params:
		format = "10 qseqid sseqid pident qcovhsp evalue",
		db = join(config["refIdentityDir"],"db/OTUs_to_ref")
	shell:
		"""
		blastp -db {params.db} -query {input.query} -out {output} -max_hsps 1 -num_threads {threads} -outfmt '{params.format}'
		"""

rule blastHitsToOTUs_capsule:
	input:
		query = join(config["refIdentityDir"],"hits_capsule.faa"),
		pdb = join(config["refIdentityDir"],"db/OTUs_to_ref.pdb"),
		phr = join(config["refIdentityDir"],"db/OTUs_to_ref.phr"),
		pin = join(config["refIdentityDir"],"db/OTUs_to_ref.pin"),
		pjs = join(config["refIdentityDir"],"db/OTUs_to_ref.pjs"),
		pot = join(config["refIdentityDir"],"db/OTUs_to_ref.pot"),
		psq = join(config["refIdentityDir"],"db/OTUs_to_ref.psq"),
		ptf = join(config["refIdentityDir"],"db/OTUs_to_ref.ptf"),
		pto = join(config["refIdentityDir"],"db/OTUs_to_ref.pto")
	output: join(config["refIdentityDir"],"out/capsule_blast.txt")
	resources:
		mem_mb=4000,
		time=240
	threads: 8
	conda:
		"../envs/blast.yml"
	params:
		format = "10 qseqid sseqid pident qcovhsp evalue",
		db = join(config["refIdentityDir"],"db/OTUs_to_ref")
	shell:
		"""
		blastp -db {params.db} -query {input.query} -out {output} -max_hsps 1 -num_threads {threads} -outfmt '{params.format}'
		"""
