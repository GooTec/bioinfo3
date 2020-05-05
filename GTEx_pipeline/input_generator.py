import json 

filelist_f = '/data/dataset/PDAC/PDAC_RNAseq_fastq/filelist.txt'

filelist_f = open(filelist_f,'r')
filelist = filelist_f.readlines()

for i in range(len(filelist)):
    filelist[i] = filelist[i][:8]

filelist = set(flielist)

for PDAC_ID in filelist:
    PDAC_ID = 'PDAC-004'
    file_dir = "/data/dataset/PDAC/PDAC_RNAseq_fastq/"
    file_name1 = '-R_1.fastq.gz'
    file_name2 = '-R_2.fastq.gz'

    fastq1 = file_dir + PDAC_ID + file_name1
    fastq2 = file_dir + PDAC_ID + file_name2

    input_dict = {
    "rnaseq_pipeline.memory": 40.0,
    "rnaseq_pipeline.num_threads": 10,
    "rnaseq_pipeline.fastq1": fastq1,
    "rnaseq_pipeline.fastq2": fastq2,
    "rnaseq_pipeline.disk_space": 20,
    "rnaseq_pipeline.prefix": PDAC_ID,
    "rnaseq_pipeline.star_index": "/home/twgoo/genome_ref/star_index_oh75.tar.gz",
    "rnaseq_pipeline.rsem_reference": "/home/twgoo/genome_ref/rsem_reference.tar.gz",
    "rnaseq_pipeline.genes_gtf": "/home/twgoo/genome_ref/genecode.v26.genes.gtf"
    }

    with open(PDAC_ID+'.json', 'w') as fp:
        json.dump(input_dict, fp)