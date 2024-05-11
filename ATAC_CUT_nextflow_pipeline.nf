#!/usr/bin/env nextflow

nextflow.enable.dsl=2


// Define script inputs
params.dataPath = ""
params.fastqFile = ""
params.bamFile = ""
params.sequencingSummaryFile = ""
params.referenceGenome = ""
params.condaEnv3in1 = ""
params.condaEnvPycoQC = ""
params.condaEnvNgsPlot = ""
params.condaEnvR = ""
params.condaEnvSicer = ""
params.experimentType = ""
params.barplot = ""
params.piechart = ""

println("params.dataPath value: " + params.dataPath)
println("Expanded dataPath value: " + params.dataPath)

// Check that required parameters are provided
if (!params.dataPath || !params.sequencingSummaryFile ||
    !params.referenceGenome || !params.condaEnv3in1 || !params.condaEnvPycoQC ||
    !params.condaEnvNgsPlot || !params.condaEnvSicer || !params.experimentType) {
    log.error "Missing required parameters! Please provide all necessary parameters."
    exit 1
}

if (params.fastqFile) {
    ch_fastqFile = Channel.fromPath(params.fastqFile)
} else if (params.bamFile) {
    ch_bamFile = Channel.fromPath(params.bamFile)
} else {
    error "No valid input file provided. Please specify either a FASTQ or BAM file."
}

ch_dataPath = Channel.fromPath(params.dataPath)
//ch_fastqFile = Channel.fromPath(params.fastqFile)
//ch_bamFile = Channel.fromPath(params.bamFile)
ch_sequencingSummaryFile = Channel.fromPath(params.sequencingSummaryFile)
ch_referenceGenome = Channel.fromPath(params.referenceGenome)
ch_condaEnvPycoQC = Channel.fromPath(params.condaEnvPycoQC)
ch_condaEnv3in1 = Channel.fromPath(params.condaEnv3in1)
ch_condaEnvNgsPlot = Channel.fromPath(params.condaEnvNgsPlot)
ch_condaEnvR = Channel.fromPath(params.condaEnvR)
ch_condaEnvSicer = Channel.fromPath(params.condaEnvSicer)
ch_experimentType = Channel.fromPath(params.experimentType)
ch_barplot = Channel.fromPath(params.barplot)
ch_piechart = Channel.fromPath(params.barplot)

//Define processes

//Process to create to directories
process createDirectories {
    
    input:
    
    path dataPath
    
    output:
    
    file "create_directories_done.txt" 
    
    script:
    """
    mkdir -p ${dataPath}/01_minionqc
    mkdir -p ${dataPath}/02_pycoQC
    mkdir -p ${dataPath}/03_fastqc
    mkdir -p ${dataPath}/04_porechop
    mkdir -p ${dataPath}/05_minimap
    mkdir -p ${dataPath}/06_samtools
    mkdir -p ${dataPath}/07_sicer2
    mkdir -p ${dataPath}/08_ngsplot
    mkdir -p ${dataPath}/09_homer
    mkdir -p ${dataPath}/10_bamtobed

    echo 'done' > create_directories_done.txt
    """
   }
   
//Process to create directories when starting from bam  File directly
process bamcreateDirectories {

	input:
	
	path dataPath
	
	output:
	
	file "create_directories_done.txt"
	
	script:
	"""
	mkdir -p ${dataPath}/02_pycoQC
	mkdir -p ${dataPath}/06_samtools
    	mkdir -p ${dataPath}/07_sicer2
    	mkdir -p ${dataPath}/08_ngsplot
    	mkdir -p ${dataPath}/09_homer
	mkdir -p ${dataPath}/10_bamtobed

	echo 'done' > create_directories_done.txt
	"""
	}
	
// Run MinionQC
process runMinionQC {

        input:
        
        file sequencingSummaryFile 
        path dataPath
        env condaEnv3in1
        file directories
	
        script:
        
        """
        export PATH=${params.condaEnv3in1}/bin:$PATH
        Rscript ${params.condaEnv3in1}/bin/MinIONQC.R -i ${sequencingSummaryFile} -o ${dataPath}/01_minionqc
        exit
        """
    }

// Run PycoQC
process runPycoQC {

        input:
        
        file sequencingSummaryFile 
        path dataPath
        env condaEnvPycoQC
        file directories
	
        script:
        
        """
	export PATH=${params.condaEnvPycoQC}/bin:$PATH
        pycoQC -f ${sequencingSummaryFile} -o ${dataPath}/02_pycoQC/pycoQC.html
        exit
        """
        
    }

// Run FastQC
process runFastQC {

        input:
        
        file fastqFile 
        path dataPath
        env condaEnv3in1
        file directories

        script:
        
        """
        export PATH=${params.condaEnv3in1}/bin:$PATH 
        fastqc -t 10 -o ${dataPath}/03_fastqc ${fastqFile}
        exit
        """
        
    }

// Run porechop_abi on the fastq.gz file

process runPorechop {
	
        input:
        file fastqFile 
        env condaEnv3in1
        path dataPath
        file directories
	
	output:
	path "${dataPath}/04_porechop/barcode.trim.fastq" 
	
        script:
        """
	export PATH=${params.condaEnv3in1}/bin:$PATH 
        porechop_abi -abi -i ${fastqFile} -o ${dataPath}/04_porechop/barcode.trim.fastq
        exit
        """
    }
    
// Use minimap2 to align the sequencing data to the reference genome
process runMinimap2 {

        input:
        
        path trimmedFastq 
        file referenceGenome
        path dataPath
        env condaEnv3in1
        
        output:
	path "${dataPath}/05_minimap/alignment.sam"

        script:
        
        """
        export PATH=${params.condaEnv3in1}/bin:$PATH 
        minimap2 -ax map-ont ${referenceGenome} ${trimmedFastq} > ${dataPath}/05_minimap/alignment.sam
        exit
        """
    }

// Convert SAM to BAM and sort
process convertSamToBam {

        input:
        
        path samFile
        path dataPath
        env condaEnv3in1 
	
	output:
	path "${dataPath}/06_samtools/sorted_alignment.bam" 
	
        script:
        
        """
        export PATH=${params.condaEnv3in1}/bin:$PATH
        samtools view -bS ${samFile} | samtools sort -o ${dataPath}/06_samtools/sorted_alignment.bam
        exit 
        """
    }
    
// sort Bam file (if the input is a bam file directly
process sortBamFile {

	input:
	
	file bamFile
	path dataPath
	env condaEnv3in1
	
	output:
	path "${dataPath}/06_samtools/sorted_alignment.bam"
	
	script:
	
	"""
	export PATH=${params.condaEnv3in1}/bin:$PATH
	samtools sort ${bamFile} -o ${dataPath}/06_samtools/sorted_alignment.bam
	exit
	"""
	
	}
	
	 
// index bam file
process indexBam {

	input:
	
	path sortedBamFile 

	path dataPath
	env condaEnv3in1
	
	output:
	path "${dataPath}/06_samtools/sorted_alignment.bam.bai" 
	
	script:
	
	"""
	export PATH=${params.condaEnv3in1}/bin:$PATH
	samtools index ${sortedBamFile} -o ${dataPath}/06_samtools/sorted_alignment.bam.bai
	exit
	"""
	
}

// Run PycoQC again with sorted BAM file
process runPycoQCPost {

        input:
        
        file sequencingSummaryFile 
        path sortedBamFile 
        path indexedBamFile
        path dataPath
        env condaEnv3in1
        
        output:
        
        path "${dataPath}/02_pycoQC/post_pycoQC.html"

        script:
        
        """
        export PATH=${params.condaEnvPycoQC}/bin:$PATH
        pycoQC -f ${sequencingSummaryFile} -o ${dataPath}/02_pycoQC/post_pycoQC.html --bam_file ${sortedBamFile}
        exit
        """
    }

// Run NGSplot based on experiment type
process runNgsPlot {

        input:
        
        file sortedBamFile 
	path dataPath
        env condaEnvNgsPlot
	
	
        script:
        
        """
        export PATH=${params.condaEnvNgsPlot}/bin:$PATH
        if [ "${params.experimentType}" == "Cut&Tag" ]; then
            ${params.condaEnvNgsPlot}/bin/ngs.plot.r -G hg38 -R genebody -L 20000 -FL 700 -C ${sortedBamFile} -O ${dataPath}/08_ngsplot/genebodyoutput -T Cut-and-Tag
            ${params.condaEnvNgsPlot}/bin/ngs.plot.r -G hg38 -R tss -L 20000 -FL 700 -C ${sortedBamFile} -O ${dataPath}/08_ngsplot/tssoutput -T Cut-and-Tag
        elif [ "${params.experimentType}" == "ATAC-Seq" ]; then
            ${params.condaEnvNgsPlot}/bin/ngs.plot.r -G hg38 -R genebody -L 5000 -FL 400 -C ${sortedBamFile} -O ${dataPath}/08_ngsplot/genebodyoutput -T ATAC-Seq
            ${params.condaEnvNgsPlot}/bin/ngs.plot.r -G hg38 -R tss -L 5000 -FL 400 -C ${sortedBamFile} -O ${dataPath}/08_ngsplot/tssoutput -T ATAC-Seq
        else
            echo "Invalid experiment type selected. Please select either ATAC-Seq or Cut&Tag"
            exit 1
        fi
        exit
        """
    }

// Filter mitochondrial contamination from the sorted alignment BAM file
process filterMitochondrialContamination {

        input:
        
        path sortedBamFile 
        path dataPath
        env condaEnv3in1
        path indexedBamFile
	
	output:
	
	path "${dataPath}/06_samtools/filtered_alignment.bam"
	
        script:
        
        """
        export PATH=${params.condaEnv3in1}/bin:$PATH
        samtools idxstats ${sortedBamFile} | cut -f 1 | grep -v chrM | xargs samtools view -b ${sortedBamFile} > ${dataPath}/06_samtools/filtered_alignment.bam
        exit
        """
    }

// Convert BAM to BED
process filteredBamToBed {

        input:

        path sortedBamFile 
	
        path dataPath
        env condaEnv3in1

        output:

        path "${dataPath}/07_sicer2/filtered_alignment.bed"

        script:

        """
        export PATH=${params.condaEnv3in1}/bin:$PATH
        bamToBed -i ${sortedBamFile} > ${dataPath}/07_sicer2/filtered_alignment.bed
        exit
        """
    }



// Perform sicer2 peak calling
process Sicerpeakcalling {

        input:

        path BedFile
        path dataPath
        env condaEnvSicer

        output:

        path "${dataPath}/07_sicer2/peaks.bed"

        script:

        """
        cd ${dataPath}/07_sicer2
        export PATH=${params.condaEnvSicer}/bin:$PATH
        sicer -t ${BedFile} -s hg38
        awk 'BEGIN{OFS="\\t"}{print \$1,\$2,\$3}' filtered_alignment-W200-G600.scoreisland > peaks.bed
        exit
        """
    }


// Perform peak annotation with HOMER
process runHomerPeakAnnotation {

        input:
        
        path peaksbed
        path dataPath
        env condaEnv3in1
        
        output:
        path "${dataPath}/09_homer/statisticspeaks.txt"
        
        script:
        
        """
        export PATH=${params.condaEnv3in1}/bin:$PATH
        annotatePeaks.pl ${peaksbed} hg38 > ${dataPath}/09_homer/annotatedpeaks.txt -annStats ${dataPath}/09_homer/statisticspeaks.txt
        exit
        """
    }

// Convert BAM to BED
process BamToBed {

        input:

        path sortedBamFile
        path dataPath
        env condaEnv3in1

        output:

        path "${dataPath}/10_bamtobed/alignment.bed" 

        script:

        """
        export PATH=${params.condaEnv3in1}/bin:$PATH
        bamToBed -i ${sortedBamFile} > ${dataPath}/10_bamtobed/alignment.bed
        exit
        """
    }

    
//make pie chart with annotated peaks from statisticspeaks.txt file

process generatePieChart {

	input:
	
	path statisticspeaksFile
	path dataPath
	path piechart
	env condaEnvR
	
	output:
	path "${dataPath}/09_homer/pie_chart.png"
	
	script:
	"""
	export PATH=${params.condaEnvR}/bin:$PATH
	Rscript ${dataPath}/PieChart.R ${statisticspeaksFile} ${dataPath}/09_homer/pie_chart.png
	exit
	"""
}
    
//Create chromosome barplot
process generateBarplot {

	input:
	
	path unfilteredbedFile
	path dataPath
	path barplot
	env condaEnvR
	
	output: 
	
	path "${dataPath}/10_bamtobed/chromosome_plot.png"
	
	script:
	"""
	export PATH=${params.condaEnvR}/bin:$PATH
	Rscript ${dataPath}/Barplot_Chr.R ${unfilteredbedFile} ${dataPath}/10_bamtobed/chromosome_plot.png
	exit
	"""
}

    
// Define the workflow

workflow {
	
	// Declare all variables that will hold process outputs at the top level of the workflow
	def directories, trimmedFastq, samFile, sortedBamFile, indexedBamFile, unfilteredbedfile

	if (params.fastqFile) {
	
		directories = createDirectories (ch_dataPath)
		runMinionQC(ch_sequencingSummaryFile, ch_dataPath, ch_condaEnvPycoQC, directories)
		runPycoQC(ch_sequencingSummaryFile, ch_dataPath, ch_condaEnvPycoQC, directories)
		runFastQC(ch_fastqFile, ch_dataPath, ch_condaEnv3in1, directories)
		trimmedFastq = runPorechop(ch_fastqFile, ch_condaEnv3in1, ch_dataPath, directories)
		samFile = runMinimap2(trimmedFastq, ch_referenceGenome, ch_dataPath, ch_condaEnv3in1)
		sortedBamFile = convertSamToBam(samFile, ch_dataPath, ch_condaEnv3in1)

	} else if (params.bamFile) {
		
		directories = bamcreateDirectories (ch_dataPath)
		sortedBamFile = sortBamFile(ch_bamFile, ch_dataPath, ch_condaEnv3in1)
	
	}
			
		
	indexedBamFile = indexBam(sortedBamFile, ch_dataPath, ch_condaEnv3in1)
	
	runPycoQCPost(ch_sequencingSummaryFile, sortedBamFile, indexedBamFile, ch_dataPath, ch_condaEnv3in1)
	
	runNgsPlot(sortedBamFile, ch_dataPath, ch_condaEnvNgsPlot)
	
	filteredBamFile = filterMitochondrialContamination(sortedBamFile, ch_dataPath, ch_condaEnv3in1, indexedBamFile)
	
	BedFile = filteredBamToBed(filteredBamFile, ch_dataPath, ch_condaEnv3in1)

	peaksbed = Sicerpeakcalling(BedFile, ch_dataPath, ch_condaEnvSicer)
	
	statisticspeaksFile = runHomerPeakAnnotation(peaksbed, ch_dataPath, ch_condaEnv3in1)
	
	unfilteredbedfile = BamToBed(sortedBamFile, ch_dataPath, ch_condaEnv3in1)

	generatePieChart(statisticspeaksFile, ch_dataPath, ch_piechart, ch_condaEnvR)
	
	generateBarplot(unfilteredbedfile, ch_dataPath, ch_barplot, ch_condaEnvR)
    
}

