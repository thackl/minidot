PATH := bin:util/minimap:$(PATH)

.PHONY: all clean samples

all:
	@echo "minidot doesn't require building. Use"
	@echo "make minimap      # to install minimap in /util"
	@echo "make sample-XXX   # to run samples: virus, prochlorococcus, arabidopsis"

clean:
	-rm -fr util
	-rm -r samples

minimap:
	hash minimap || $(MAKE) util/minimap

util/minimap:
	mkdir -p util
	cd util && git clone https://github.com/lh3/minimap
	cd util/minimap && make

## prochlorococcus sample
PROC=samples/prochlorococcus
PMED4=$(PROC)/MED4.fa
PSB=$(PROC)/SB.fa
PNATLA2=$(PROC)/NATLA2.fa
PLG=$(PROC)/LG.fa
PMIT9313=$(PROC)/MIT9313.fa

sample-prochlorococcus: minimap $(PROC) $(PMED4) $(PSB) $(PNATLA2) $(PLG) $(PMIT9313)
	@echo "\nRunning minidot"
	minidot -s -o $(PROC)/prochlorococcus.pdf $(PMED4) $(PSB) $(PNATLA2) $(PLG) $(PMIT9313)

$(PROC):
	mkdir -p $(PROC)

$(PMED4):
	curl -# ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000011465.1_ASM1146v1/GCF_000011465.1_ASM1146v1_genomic.fna.gz | gunzip > $(PMED4)

$(PSB):
	curl -# ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000760115.1_ASM76011v1/GCF_000760115.1_ASM76011v1_genomic.fna.gz | gunzip > $(PSB)

$(PNATLA2):
	curl -# ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000760155.1_ASM76015v1/GCF_000760155.1_ASM76015v1_genomic.fna.gz | gunzip > $(PNATLA2)

$(PLG):
	curl -# ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000012465.1_ASM1246v1/GCF_000012465.1_ASM1246v1_genomic.fna.gz | gunzip > $(PLG)

$(PMIT9313):
	curl -# ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000011485.1_ASM1148v1/GCF_000011485.1_ASM1148v1_genomic.fna.gz | gunzip > $(PMIT9313)

## arabidopsis sample
ARAB=samples/arabidopsis
ATAL=$(ARAB)/A.thaliana.fa
ALYR=$(ARAB)/A.lyrata.fa

sample-arabidopsis: minimap $(ARAB) $(ATAL) $(ALYR)
	@echo "\nRunning minidot"
	minidot -o $(ARAB)/arabidopsis.pdf $(ALYR) $(ATAL)

$(ARAB):
	mkdir -p $(ARAB)
$(ATAL):
	curl -# ftp://ftp.ensemblgenomes.org/pub/plants/current/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.31.dna.genome.fa.gz | gunzip > $(ATAL)
$(ALYR):
	curl -# ftp://ftp.ensemblgenomes.org/pub/plants/current/fasta/arabidopsis_lyrata/dna/Arabidopsis_lyrata.v.1.0.31.dna.genome.fa.gz | gunzip > $(ALYR)



