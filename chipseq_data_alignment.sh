

### Create working directories
    base=/Users/aayushojha/Documents/chipseq
    dir_data=/Users/aayushojha/Documents/chipseq/data
    dir_reference=/Users/aayushojha/Documents/chipseq/data/reference
    dir_software=/Users/aayushojha/Documents/sequencing_projects/software

    cd ${dir_data}
    

### Download samples
        cd ${dir_data}
            mkdir fastq
            dir_fastq=/Users/aayushojha/Documents/chipseq/data/fastq
            cd ${dir_fastq}
                wget ....
    

### Download reference genome
        cd ${dir_data}
            mkdir reference
                dir_reference=/Users/aayushojha/Documents/chipseq/data/reference
                cd ${dir_reference}
                    # Go to UCSC and downlaod the grch38 reference genome an
                    # save hg38.fa.gz here.
                    # Unzip it to get hg38.fa

                    ref_genome=${dir_reference}/hg38.fa

                    cat ${ref_genome} | less -S
                    # press 'q' to Quit seeing contents of ref_genome. 




### READ MAPPING USING bowtie2
    #1. INSTALL bowtie2
        dir_software=/Users/aayushojha/Documents/sequencing_projects/software
        cd ${dir_software}
            mkdir bowtie2
            cd ${dir_software}/bowtie2
                conda install bowtie2 -c bioconda


    #2. Download index
        # The aligners usually needs an index file that is created using the genome reference files
        # (fasta) we downloaded from UCSC, NCBI, ENSEMBL or GENCODE.
        # For bowtie2, there is a pre-built index file for us to download.
        # Scroll down the page https://bowtie-bio.sourceforge.net/bowtie2/index.shtml and on the left
        # you will see the index to download.
        # We will download the GRCH38 (aka, hg38 without decoy).
        # click it, it will download a GRCh38_noalt_as.zip file of 3.5G.
        # unzip it and place the folder to the data/reference folder.
        # If there is no pre-built index, take a look at this nf-core/references.


    #3. Running bowtie2
        # make sure you are in the data folder
        cd data/

        bowtie2 -x reference/GRCh38_noalt_as/GRCh38_noalt_as -U fastq/YAP.fastq.gz -S fastq/YAP1.sam --threads 8  --no-mixed --no-discordant -k 1
        #This is done to check if everthying is working.
        
        # Once everything is working, we can run the code in loop for
        # every fastq.gz file in fastq directory. 
        # The outputs will be saved in the fastq dir as .sam file. 
        for fq in fastq/*fastq.gz
        do
            bowtie2 -x reference/GRCh38_noalt_as/GRCh38_noalt_as -U $fq -S ${fq/fastq.gz/sam} --threads 8  --no-mixed --no-discordant -k 1
        done


### SAM to BAM
    #1. INSTALL samtools FROM: https://github.com/samtools/samtools
        mkdir /Users/aayushojha/Documents/sequencing_projects/software/samtools
        cd /Users/aayushojha/Documents/sequencing_projects/software/samtools
        wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 -O samtools.tar.bz2
        tar -xjvf samtools.tar.bz2
        cd samtools-{version}
        make
        sudo make prefix=/Users/aayushojha/Documents/sequencing_projects/software/samtools install

    # 2. Loading samtools 
        file_samtools=/Users/aayushojha/Documents/sequencing_projects/software/samtools/bin/samtools

    #3. Running samtools for multiple sam files at once to convert them 
        #into their respective bam output
        for sam in fastq/*sam
        do
            ${file_samtools} view -bS $sam > ${sam/sam/bam}
        done
    
### BAM to SORTED BAM + Index
    #1. Converting BAM to Sorted BAM for multiple files at once
        for bam in fastq/*bam
        do
            ${file_samtools} sort -@ 4 -o ${bam/bam/sorted.bam} $bam
        done

      
    #2. remove all the unsorted bam 
        ls fastq/*bam 
        ls fastq/*bam | grep -v sorted.bam | xargs rm

    #3. Indexing the sorted bam file 
        for bam in fastq/*sorted.bam
        do
            ${file_samtools} index $bam
        done



### Deeptools 

    #1. Install deeptools 
        # For my mac, the deeptools did not install in previous environment (reproduce_figure),
        # so I created a completely new environment to install and run deeptools
        # New env named "deeptools_env"
        conda create -n deeptools_env -c bioconda -c conda-forge deeptools
        conda activate deeptools_env    

    #2. Check if deeptools is running by using bamCoverage for a single file

        # bamCoverage is a tool from deepTools that converts BAM files (binary alignment files
        # from sequencing data) into BigWig files for visualization and analysis.
        # It calculates read coverage across the genome and allows for  normalization to account
        # for sequencing depth differences.

        cd ${dir_fastq}
        bamCoverage --bam YAP.sorted.bam --normalizeUsing RPKM --extendReads 200 -o YAP1.bw
    
    #3. Run bamCoverage for multiple files at once
        for bam in *sorted.bam
        do
            bamCoverage --bam $bam --normalizeUsing RPKM --extendReads 200 -o ${bam/sorted.bam/bw}
        done
    

### Peak calling
       
    #1. Install macs
        # MACS is the most popular peak caller for ChIPseq.
        # MACS is now MACS3! https://macs3-project.github.io/MACS/docs/INSTALL.html
        pip install --upgrade macs3

    #2. Check if the callpeak subcommand works or not. 
        #The MACS3 callpeak subcommands have many paramters and
        # you want to read https://macs3-project.github.io/MACS/docs/callpeak.html

        macs3 callpeak -t YAP.sorted.bam -c IgG.sorted.bam -f BAM -n YAP -g hs --outdir YAP_peak

    #3. Use callpeak to find the peaks for all samples at once. 
        for bam in *sorted.bam
        do
            if [[ "$bam" != "IgG.sorted.bam" ]]; then
                macs3 callpeak -t $bam -c IgG.sorted.bam -f BAM -n ${bam%.sorted.bam} -g hs --outdir ${bam/.sorted.bam/_peak}
            fi
        done
    
    #4. Check how many peaks we get for each transcription factor
        find . -name "*Peak"  | xargs wc -l
            # 11512 ./TEAD4_peak/TEAD4_peaks.narrowPeak
            # 10719 ./TAZ_peak/TAZ_peaks.narrowPeak
            #     9807 ./YAP_peak/YAP_peaks.narrowPeak
            # 32038 total


### Converting .bed file to .bw file for histone modification files
    dir_data=/Users/aayushojha/Documents/chipseq/data
    cd ${dir_data}/public_data

    # Download bedGraphToBigWig from UCSC
        wget http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.arm64/bedGraphToBigWig
        chmod +x bedGraphToBigWig

    # Download .chrom.sizes file from UCSC
        wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
        # if hg38 needed:
        # wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes

    # Load bedGraphToBigWig:
        bgtobw=/Users/aayushojha/Documents/chipseq/data/public_data/bedGraphToBigWig
        hg19_chrom_sizes=/Users/aayushojha/Documents/chipseq/data/public_data/hg19.chrom.sizes
        # hg38_chrom_sizes=/Users/aayushojha/Documents/chipseq/data/public_data/hg38.chrom.sizes
                    
    
        # Use bedGraphToBigWig to convert the bedGraph files to bigWig files
            ${bgtobw} GSM1204470_MDAMB231.H3K4me1_1.hg19.tags.bedGraph hg19.chrom.sizes H3K4me1.bw
            ${bgtobw} GSM1204472_MDAMB231.H3K4me3_1.hg19.tags.bedGraph hg19.chrom.sizes H3K4me3.bw
            ${bgtobw} GSM1204474_MDAMB231.H3K27Ac_1.hg19.tags.bedGraph hg19.chrom.sizes H3K27ac.bw

    # If you have multiple bed files, you can sort them first and then convert them to bigwig files 
        # sort -k1,1 -k2,2n H3K4me1.bed > H3K4me1.sorted.bedGraph
        # awk '{print $1, $2, $3, $5}' H3K4me1.sorted.bedGraph > H3K4me1.final.bedGraph
        # ${bgtobw} H3K4me1.final.bedGraph hg38.chrom.sizes H3K4me1.bw

        # sort -k1,1 -k2,2n H3K4me3.bed > H3K27ac.sorted.bedGraph
        # awk '{print $1, $2, $3, $5}' H3K27ac.sorted.bedGraph > H3K27ac.final.bedGraph
        # ${bgtobw} H3K27ac.final.bedGraph hg38.chrom.sizes H3K27ac.bw

        # sort -k1,1 -k2,2n H3K27ac.bed > H3K27ac.sorted.bedGraph
        # awk '{print $1, $2, $3, $5}' H3K27ac.sorted.bedGraph > H3K27ac.final.bedGraph
        # ${bgtobw} H3K27ac.final.bedGraph hg38.chrom.sizes H3K27ac.bw

