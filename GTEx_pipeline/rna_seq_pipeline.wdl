task star {
    # Input Files 
    
    # Paired end read로 sequencing을 진행한 경우, 1,2 개의 파일이 생성됨. 
    File fastq1 
    File? fastq2
    
    String prefix
    # STAR software를 이용해 STAR index 파일을 먼저 생성해주어야함. 
    # Ensemble DB에서 본인이 원하는 종의 Reference genome assembly와 Annotation File을 이용해 생성가능
    # 저는 이번 분석에서 Humon genome reference 38 version을 이용하였습니다.
    # Hg38 reference genome assembly file link: ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/
    # Hg38 genome annotation file link: ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens

    # STAR index를 통해 생성된 파일들을 압축해 File 변수로 전달해주어야함. 
    File star_index

    # STAR options
    Int? outFilterMultimapNmax
    Int? alignSJoverhangMin
    Int? alignSJDBoverhangMin
    Int? outFilterMismatchNmax
    Float? outFilterMismatchNoverLmax
    Int? alignIntronMin
    Int? alignIntronMax
    Int? alignMatesGapMax
    String? outFilterType
    Float? outFilterScoreMinOverLread
    Float? outFilterMatchNminOverLread
    Int? limitSjdbInsertNsj
    String? outSAMstrandField
    String? outFilterIntronMotifs
    String? alignSoftClipAtReferenceEnds
    String? quantMode
    String? outSAMattrRGline
    String? outSAMattributes
    File? varVCFfile
    String? waspOutputMode
    Int? chimSegmentMin
    Int? chimJunctionOverhangMin
    String? chimOutType
    Int? chimMainSegmentMultNmax
    Int? chimOutJunctionFormat
    File? sjdbFileChrStartEnd

    # Runtime Variables 
    Float memory
    Int disk_space
    Int num_threads

    command {
        set -euo pipefail

        if [[ ${fastq1} == *".tar" || ${fastq1} == *".tar.gz" ]]; then
            tar -xvvf ${fastq1}
            fastq1_abs=$(for f in *_1.fastq*; do echo "$(pwd)/$f"; done | paste -s -d ',')
            fastq2_abs=$(for f in *_2.fastq*; do echo "$(pwd)/$f"; done | paste -s -d ',')
            if [[ $fastq1_abs == *"*_1.fastq*" ]]; then  # no paired-end FASTQs found; check for single-end FASTQ
                fastq1_abs=$(for f in *.fastq*; do echo "$(pwd)/$f"; done | paste -s -d ',')
                fastq2_abs=''
            fi
        else
            # make sure paths are absolute
            fastq1_abs=${fastq1}
            fastq2_abs=${fastq2}
            if [[ $fastq1_abs != /* ]]; then
                fastq1_abs=$PWD/$fastq1_abs
                fastq2_abs=$PWD/$fastq2_abs
            fi
        fi

        echo "FASTQs:"
        echo $fastq1_abs
        echo $fastq2_abs

        # extract index
        echo $(date +"[%b %d %H:%M:%S] Extracting STAR index")
        mkdir star_index
        tar -xvvf ${star_index} -C star_index --strip-components=1

        mkdir star_out
        # placeholders for optional outputs
        touch star_out/${prefix}.Aligned.toTranscriptome.out.bam
        touch star_out/${prefix}.Chimeric.out.sorted.bam
        touch star_out/${prefix}.Chimeric.out.sorted.bam.bai
        touch star_out/${prefix}.ReadsPerGene.out.tab  # run_STAR.py will gzip

        /src/run_STAR.py \
            star_index $fastq1_abs $fastq2_abs ${prefix} \
            --output_dir star_out \
            ${"--outFilterMultimapNmax " + outFilterMultimapNmax} \
            ${"--alignSJoverhangMin " + alignSJoverhangMin} \
            ${"--alignSJDBoverhangMin " + alignSJDBoverhangMin} \
            ${"--outFilterMismatchNmax " + outFilterMismatchNmax} \
            ${"--outFilterMismatchNoverLmax " + outFilterMismatchNoverLmax} \
            ${"--alignIntronMin " + alignIntronMin} \
            ${"--alignIntronMax " + alignIntronMax} \
            ${"--alignMatesGapMax " + alignMatesGapMax} \
            ${"--outFilterType " + outFilterType} \
            ${"--outFilterScoreMinOverLread " + outFilterScoreMinOverLread} \
            ${"--outFilterMatchNminOverLread " + outFilterMatchNminOverLread} \
            ${"--limitSjdbInsertNsj " + limitSjdbInsertNsj} \
            ${"--outSAMstrandField " + outSAMstrandField} \
            ${"--outFilterIntronMotifs " + outFilterIntronMotifs} \
            ${"--alignSoftClipAtReferenceEnds " + alignSoftClipAtReferenceEnds} \
            ${"--quantMode " + quantMode} \
            ${"--outSAMattrRGline " + outSAMattrRGline} \
            ${"--outSAMattributes " + outSAMattributes} \
            ${"--varVCFfile " + varVCFfile} \
            ${"--waspOutputMode " + waspOutputMode} \
            ${"--chimSegmentMin " + chimSegmentMin} \
            ${"--chimJunctionOverhangMin " + chimJunctionOverhangMin} \
            ${"--chimOutType " + chimOutType} \
            ${"--chimMainSegmentMultNmax " + chimMainSegmentMultNmax} \
            ${"--chimOutJunctionFormat " + chimOutJunctionFormat} \
            ${"--sjdbFileChrStartEnd " + sjdbFileChrStartEnd} \
            --threads ${num_threads}

        
        
        rm -rf star_index
    }


    output {
        File bam_file = "star_out/${prefix}.Aligned.sortedByCoord.out.bam"
        File bam_index = "star_out/${prefix}.Aligned.sortedByCoord.out.bam.bai"
        File transcriptome_bam = "star_out/${prefix}.Aligned.toTranscriptome.out.bam"
        File chimeric_junctions = "star_out/${prefix}.Chimeric.out.junction.gz"
        File chimeric_bam_file = "star_out/${prefix}.Chimeric.out.sorted.bam"
        File chimeric_bam_index = "star_out/${prefix}.Chimeric.out.sorted.bam.bai"
        File read_counts = "star_out/${prefix}.ReadsPerGene.out.tab.gz"
        File junctions = "star_out/${prefix}.SJ.out.tab.gz"
        File junctions_pass1 = "star_out/${prefix}._STARpass1/${prefix}.SJ.pass1.out.tab.gz"
        Array[File] logs = ["star_out/${prefix}.Log.final.out", "star_out/${prefix}.Log.out", "star_out/${prefix}.Log.progress.out"]
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_rnaseq:V9"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
    }

    meta {
        author: "Francois Aguet"
    }
}


task markduplicates {

    File input_bam
    String prefix
    Int? max_records_in_ram
    Float? sorting_collection_size_ratio

    Float memory
    Int java_memory = floor(memory - 0.5)
    Int disk_space
    Int num_threads

    String output_bam = sub(basename(input_bam), "\\.bam$", ".md.bam")

    command {
        set -euo pipefail
        python3 -u /src/run_MarkDuplicates.py ${input_bam} ${prefix} \
            --memory ${java_memory} \
            ${"--max_records_in_ram " + max_records_in_ram} \
            ${"--sorting_collection_size_ratio " + sorting_collection_size_ratio}
        samtools index ${output_bam}
    }

    output {
        File bam_file = "${output_bam}"
        File bam_index = "${output_bam}.bai"
        File metrics = "${prefix}.marked_dup_metrics.txt"
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_rnaseq:V9"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
    }

    meta {
        author: "Francois Aguet"
    }
}

task rsem {

    File transcriptome_bam
    File rsem_reference
    String prefix

    Float memory
    Int disk_space
    Int num_threads

    Int? max_frag_len
    String? estimate_rspd
    String? is_stranded
    String? paired_end

    command {
        set -euo pipefail
        mkdir rsem_reference
        tar -xvvf ${rsem_reference} -C rsem_reference --strip-components=1

        /src/run_RSEM.py \
            ${"--max_frag_len " + max_frag_len} \
            ${"--estimate_rspd " + estimate_rspd} \
            ${"--is_stranded " + is_stranded} \
            ${"--paired_end " + paired_end} \
            --threads ${num_threads} \
            rsem_reference ${transcriptome_bam} ${prefix}
        gzip *.results

        rm -rf rsem_reference
    }

    output {
        File genes="${prefix}.rsem.genes.results.gz"
        File isoforms="${prefix}.rsem.isoforms.results.gz"
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_rnaseq:V9"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
    }

    meta {
        author: "Francois Aguet"
    }
}

task rnaseqc2 {

    File bam_file
    File genes_gtf
    String sample_id
    String? strandedness 
    File? intervals_bed
    String? flags

    Float memory
    Int disk_space
    Int num_threads

    command {
        set -euo pipefail
        echo $(date +"[%b %d %H:%M:%S] Running RNA-SeQC 2")
        touch ${sample_id}.fragmentSizes.txt
        rnaseqc ${genes_gtf} ${bam_file} . -s ${sample_id} ${"--bed " + intervals_bed} ${"--stranded " + strandedness} -vv ${flags}
        echo "  * compressing outputs"
        gzip *.gct
        echo $(date +"[%b %d %H:%M:%S] done")
    }

    output {
        File gene_tpm = "${sample_id}.gene_tpm.gct.gz"
        File gene_counts = "${sample_id}.gene_reads.gct.gz"
        File exon_counts = "${sample_id}.exon_reads.gct.gz"
        File metrics = "${sample_id}.metrics.tsv"
        File insertsize_distr = "${sample_id}.fragmentSizes.txt"
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_rnaseq:V9"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
    }

    meta {
        author: "Francois Aguet"
    }
}

workflow star_workflow {
    String prefix
    Int disk_space
    Float memory
    Int num_threads
    File fastq1 
    File? fastq2
    File star_index
    File rsem_reference
    File genes_gtf
    
    call star{
        input : fastq1= fastq1,
            fastq2 = fastq2,
            star_index=star_index,
            prefix=prefix, 
            memory=memory,
            num_threads=num_threads, 
            disk_space=disk_space
    }

    call markduplicates {
        input: input_bam=star.bam_file, 
            prefix=prefix, 
            memory=memory,
            num_threads=num_threads, 
            disk_space=disk_space
    }
    call rsem{
        input: transcriptome_bam=star.transcriptome_bam,
            rsem_reference=rsem_reference,
            prefix=prefix,
            memory=memory,
            disk_space=disk_space,
            num_threads=num_threads
    }

    call rnaseqc2{
        input: bam_file=markduplicates.bam_file,
            sample_id =prefix,
            genes_gtf=genes_gtf,
            memory=memory,
            disk_space=disk_space,
            num_threads=num_threads
    }
}
