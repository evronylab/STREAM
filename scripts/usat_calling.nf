// enable DSL2 syntax
nextflow.enable.dsl=2


/*

  CHANNELS

*/

// load sample metadata and FASTQ paths
if( params.input_type == 'fastq' ) {
  samples_ch =
  Channel
    .fromPath( params.samples )
    .splitCsv( sep: '\t' )
// load sample metadata and CRAM paths
} else if( params.input_type == 'cram' ) {
  Channel
    .fromPath( params.samples )
    .splitCsv( sep: '\t' )
    .multiMap { it ->
      all: it
      cram: it[2]
      index: it[3]
    }
    .set {samples_ch}
} else {
  println 'No input file type provided'
}


/*

  WORKFLOW

*/

workflow {

  if( params.input_type == 'fastq' ) {
    CONCAT_FASTQ( samples_ch )
    FASTQC( CONCAT_FASTQ.out.fastq )
    ALIGN_BAM( CONCAT_FASTQ.out.fastq )
    REMOVE_DUPS( ALIGN_BAM.out.bam )
    SORT_BAM( REMOVE_DUPS.out.unmarked )
    if( params.usats_probes_picard == 'false') {
      println 'Not performing hybrid capture quality analysis, since no probes specified'
    }
    else {
      PICARD( SORT_BAM.out.cram )
    }
    SPLIT_REGIONS()
    HIPSTR( SORT_BAM.out.hipstr.toList(), SPLIT_REGIONS.out.usats.flatten(), SORT_BAM.out.index.collect() )
    CONCAT_VCF( HIPSTR.out.small_vcf.collect(), HIPSTR.out.index.collect() )
    SPLIT_VCF( SORT_BAM.out.cram, CONCAT_VCF.out.joint_vcf )
    EXPANSION_HUNTER( SORT_BAM.out.cram )
    GANGSTR( SORT_BAM.out.cram )
    BEDTOOLS_COVERAGE( SORT_BAM.out.cram )
    if( params.supp_panel_picard == 'false') {
      println 'Not performing supplementary panel genotyping, since no panel specified'
    }
    else if( params.supp_probes_picard == 'false' ) {
      SUPP_VARIANTS( SORT_BAM.out.cram )
      println 'Not performing supplementary panel capture quality analysis, since no probes specified'
    }
    else {
      SUPP_VARIANTS( SORT_BAM.out.cram )
      SUPP_PICARD( SORT_BAM.out.cram )
    }
    QUERY_HIPSTR( SPLIT_VCF.out.vcf )
    QUERY_EH( EXPANSION_HUNTER.out.vcf )
    QUERY_GANGSTR( GANGSTR.out.vcf )
    JOIN_CALLS( SORT_BAM.out.cram.join(QUERY_HIPSTR.out.table).join(QUERY_EH.out.table).join(QUERY_GANGSTR.out.table).join(BEDTOOLS_COVERAGE.out.table) )
    JOIN_SAMPLES( JOIN_CALLS.out.rds.collect() )
  }
  
  else if( params.input_type == 'cram' ) {
    if( params.usats_probes_picard == 'false') {
      println 'Not performing hybrid capture quality analysis, since no probes specified'
    }
    else {
      PICARD( samples_ch.all )
    }
    SPLIT_REGIONS()
    HIPSTR( samples_ch.cram.toList(), SPLIT_REGIONS.out.usats.flatten(), samples_ch.index.collect() )
    CONCAT_VCF( HIPSTR.out.small_vcf.collect(), HIPSTR.out.index.collect() )
    SPLIT_VCF( samples_ch.all, CONCAT_VCF.out.joint_vcf )
    EXPANSION_HUNTER( samples_ch.all )
    GANGSTR( samples_ch.all )
    BEDTOOLS_COVERAGE( samples_ch.all )
    if( params.supp_panel_picard == 'false') {
      println 'Not performing supplementary panel genotyping, since no panel specified'
    }
    else if( params.supp_probes_picard == 'false' ) {
      SUPP_VARIANTS( samples_ch.all )
      println 'Not performing supplementary panel capture quality analysis, since no probes specified'
    }
    else {
      SUPP_VARIANTS( samples_ch.all )
      SUPP_PICARD( samples_ch.all )
    }
    QUERY_HIPSTR( SPLIT_VCF.out.vcf )
    QUERY_EH( EXPANSION_HUNTER.out.vcf )
    QUERY_GANGSTR( GANGSTR.out.vcf )
    JOIN_CALLS( samples_ch.all.join(QUERY_HIPSTR.out.table).join(QUERY_EH.out.table).join(QUERY_GANGSTR.out.table).join(BEDTOOLS_COVERAGE.out.table) )
    JOIN_SAMPLES( JOIN_CALLS.out.rds.collect() )  
  }
  
  else {
    println 'Invalid input file type'
  }
}


/*

  PROCESSES

*/

// concatenate separate FASTQ files

process CONCAT_FASTQ {

  time '4h'
  memory '32 GB'

  input:
    // load FASTQ files
    tuple val( sampleID ), val( sex ), val( read1files ), val( read2files ), val( trio ), val( sampleType ), val(opticalDistance)

  output:
    // output channel with concatenated FASTQ files
    tuple val( sampleID ), val( sex ), path( "${sampleID}.R1.fastq.gz" ), path( "${sampleID}.R2.fastq.gz" ), val( trio ), val( sampleType ), val(opticalDistance), emit: fastq

  script:
  """
  IFS="," read -a read1array <<< ${read1files}

  IFS="," read -a read2array <<< ${read2files}

  cat \${read1array[@]} > ${sampleID}.R1.fastq.gz

  cat \${read2array[@]} > ${sampleID}.R2.fastq.gz

  """
}


// QC raw sequencing data

process FASTQC {

  time '6h'
  memory '32 GB'

  publishDir("${params.results}/${sampleID}_results", mode: 'copy')

  input:
    // load concatenated FASTQ files
    tuple val( sampleID ), val( sex ), path( read1 ), path( read2 ), val( trio ), val( sampleType ), val(opticalDistance)

  output:
    // save HTML reports
    tuple val( sampleID ), path("*_fastqc.html")

  script:
  """
  module load ${params.fastqc}

  # make results directory if it doesn't already exist
  mkdir -p \$(dirname ${params.results}/${sampleID}_results)
  
  # run fastqc on read 1 and read 2 FASTQ files
  fastqc ${read1} ${read2}

  """
}


// Align reads to reference genome

process ALIGN_BAM {

  time '30h'
  cpus 16
  memory '32 GB'

  input:
    // load concatenated FASTQ files
    tuple val( sampleID ), val( sex ), path( read1 ), path( read2 ), val( trio ), val( sampleType ), val(opticalDistance)

  output:
    // output channel for REMOVE_DUPS
    tuple val( sampleID ), val( sex ), path( "${sampleID}.bam" ), val( trio ), val( sampleType ), val(opticalDistance), emit: bam

  script:
  """
  module load ${params.bwa}

  bwa mem -R "@RG\\tID:${sampleID}\\tSM:${sampleID}\\tLB:${sampleID}\\tPL:ILLUMINA" -t 16 ${params.reference} ${read1} ${read2} > ${sampleID}.bam
  
  """
}


// Remove optical duplicates from BAM

process REMOVE_DUPS {

  time '10h'
  cpus 2
  memory '32 GB'

  publishDir("${params.results}/${sampleID}_results", mode: 'copy', pattern: "*.markdup.metrics.txt")

  input:
    // load BAM file
    tuple val( sampleID ), val( sex ), path( bam ), val( trio ), val( sampleType ), val(opticalDistance)

  output:
    // output channel for SORT_BAM
    tuple val( sampleID ), val( sex ), path( "${sampleID}.unmarked.bam" ), val( trio ), val( sampleType ), emit: unmarked
    // save duplication metrics
    path( "${sampleID}.markdup.metrics.txt" )

  script:
  """
  module load ${params.gatk}
  
  # Remove optical duplicates from BAM
  java -jar ${params.picard} MarkDuplicates \\
    I=${bam} \\
    O=${sampleID}.markdup.bam \\
    REMOVE_DUPLICATES=false \\
    METRICS_FILE=${sampleID}.markdup.metrics.txt \\
    OPTICAL_DUPLICATE_PIXEL_DISTANCE=${opticalDistance} \\
    ASSUME_SORT_ORDER=queryname \\
    CLEAR_DT=false \\
    REMOVE_SEQUENCING_DUPLICATES=true \\
    READ_NAME_REGEX=[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9\\-]+:[0-9]:\\([0-9]+\\):\\([0-9]+\\):\\([0-9]+\\).* \\
    TMP_DIR=${params.gatk_tmp}
    
  # Unmark duplicate reads
  gatk UnmarkDuplicates \
    -I ${sampleID}.markdup.bam \
    -O ${sampleID}.unmarked.bam
  """
}


// Sort BAM file and output as CRAM

process SORT_BAM {

  time '10h'
  cpus 2
  memory '32 GB'

  publishDir("${params.results}/${sampleID}_results", mode: 'copy')

  input:
    // load FASTQ files
    tuple val( sampleID ), val( sex ), path( unmarked ), val( trio ), val( sampleType )

  output:
    // output channel for most processes
    tuple val( sampleID ), val( sex ), path( "${sampleID}.cram" ), path( "${sampleID}.cram.crai" ), val( trio ), val( sampleType ), emit: cram

    // output channels for HipSTR
    // enables HipSTR joint calling
    path( "${sampleID}.cram" ), emit: hipstr
    path( "${sampleID}.cram.crai" ), emit: index

  script:
  """
  module load ${params.samtools}
  
  # Sort BAM  
  java -jar ${params.picard} SortSam \\
    I=${unmarked} \\
    O=${sampleID}.unmarked.sorted.bam \\
    SORT_ORDER=coordinate \\
    TMP_DIR=${params.gatk_tmp}; samtools index ${sampleID}.unmarked.sorted.bam
  
  # Convert to CRAM
  samtools view -F 2304 -C -T ${params.reference} --write-index -O cram -@4 -o ${sampleID}.cram##idx##${sampleID}.cram.crai ${sampleID}.unmarked.sorted.bam
  """
}


// Check hybrid capture quality

process PICARD {

  time '8h'
  memory '32 GB'

  publishDir("${params.results}/${sampleID}_results", mode: 'copy')

  input:
    // load CRAM
    // CRAM should have duplicates removed
    tuple val( sampleID ), val( sex ), path( cram ), path( index ), val( trio ), val( sampleType )

  output:
    // save Picard CollectHsMetrics report
    tuple val( sampleID ), path( "${sampleID}_picard_usats.txt" ), path( "${sampleID}_picard_probes.txt" )

  script:
  """
  module load ${params.jdk}
  # run Picard with microsatellites as "target" coordinates
  # check zero-coverage targets from this file
  java -jar ${params.picard} CollectHsMetrics \\
        I=${cram} \\
        O=${sampleID}_picard_usats.txt \\
        R=${params.reference} \\
        BAIT_INTERVALS=${params.usats_probes_picard} \\
        TARGET_INTERVALS=${params.usats_panel_picard}

  # run Picard with probes as "target" coordinates
  # check capture quality metrics from this file
  java -jar ${params.picard} CollectHsMetrics \\
        I=${cram} \\
        O=${sampleID}_picard_probes.txt \\
        R=${params.reference} \\
        BAIT_INTERVALS=${params.usats_probes_picard} \\
        TARGET_INTERVALS=${params.usats_probes_picard}

  """
}


// Divide HipSTR input file into smaller files
// Will use to submit smaller, faster jobs

process SPLIT_REGIONS {

  time '20m'

  output:
    // channel divided input files
    path( "split_regions_*" ), emit: usats

  script:
  """
  # use numerical suffixes to avoid file name overlap
  split -n l/${params.split} --numeric-suffixes --additional-suffix=.txt ${params.hipstr_usats} split_regions_

  """
}


// Run HipSTR jointly on all CRAMs using divided list of microsatellites

process HIPSTR {

  time '10h'
  memory '32 GB'
  
  publishDir("${params.results}/hipstr_aln_viz", pattern: "${params.outputbasename}_*.aln.viz.gz", mode: 'copy')

  input:
    // load CRAMs
    // must be just CRAMs to make comma separated list for HipSTR
    path( cram )

    // load microsatellite lists
    path( usats )

    // load all CRAM indexes
    path( index )

  output:
    // save HipSTR visualization files
    path( "${params.outputbasename}_*.aln.viz.gz" )
  
    // channel for VCFs to be concatenated
    path( "array_${usats}_hipstr.vcf.gz" ), emit: small_vcf

    // channel for VCF indexes
    path( "array_${usats}_hipstr.vcf.gz.csi" ), emit: index

  script:
  """
  module load ${params.bcftools}
  
  # Do joint calling on all samples

  # Make comma separated list of CRAMs
  BAMS=\$(echo ${cram} | sed 's/ /,/g')
  
  # Check if all samples are male
  CHECK_SEX=\$(cut -f 2 ${params.samples} | sort | uniq)
  
  if [ "\$CHECK_SEX" = "M" ]; then
  
   # run HipSTR with halpoid X/Y calling if all samples are male
    ${params.hipstr} \
      --bams \${BAMS} \
      --fasta ${params.reference} \
      --regions ${usats} \
      --str-vcf array_${usats}_hipstr.vcf.gz \
      --viz-out ${params.outputbasename}_${usats}.aln.viz.gz \
      --min-reads 20 \
      --max-str-len 150 \
      --no-rmdup \
      --haploid-chrs chrX,chrY \
      --output-filters; bcftools index array_${usats}_hipstr.vcf.gz
  
  else

    # run HipSTR with diploid calling on all chromosomes if not all samples are male
    ${params.hipstr} \
      --bams \${BAMS} \
      --fasta ${params.reference} \
      --regions ${usats} \
      --str-vcf array_${usats}_hipstr.vcf.gz \
      --viz-out ${params.outputbasename}_${usats:14}.aln.viz.gz \
      --min-reads 20 \
      --max-str-len 150 \
      --no-rmdup \
      --output-filters; bcftools index array_${usats}_hipstr.vcf.gz
      
  fi

  """
}


// Combine HipSTR VCFs into one VCF containing all samples and loci

process CONCAT_VCF {

  time '3h'

  publishDir("${params.results}", mode: 'copy')

  input:
    // load all HipSTR VCFs
    // must be just VCFs to output as list in script
    path( small_vcf )

    // load all HipSTR VCF indexes
    path( index )

  output:
    // channel for combined VCF
    // save combined VCF
    tuple path( "${params.outputbasename}_hipstr.vcf.gz" ), path( "${params.outputbasename}_hipstr.vcf.gz.csi" ), emit: joint_vcf

  script:
  """
  module load ${params.bcftools}

  # concatenate all HipSTR VCFs and index combined VCF
  bcftools concat --allow-overlaps ${small_vcf} | bcftools sort --output ${params.outputbasename}_hipstr.vcf.gz -Oz -;
  bcftools index ${params.outputbasename}_hipstr.vcf.gz

  """
}


// Split combined HipSTR VCF by sample

process SPLIT_VCF {

  time '2h'

  publishDir("${params.results}/${sampleID}_results", mode: 'copy')

  input:
    // load sample metadata
    tuple val( sampleID ), val( sex ), path( cram ), path( index ), val( trio ), val( sampleType )

    // load combined HipSTR VCF
    tuple path( joint_vcf ), path( index )

  output:
    // channel for individual VCFs
    // save individual VCFs
    tuple val( sampleID ), path( "${sampleID}_hipstr.vcf.gz" ), emit: vcf

  script:
  """
  module load ${params.bcftools}

  # split combined VCF by sample ID
  bcftools view -Oz -s ${sampleID} ${joint_vcf} > ${sampleID}_hipstr.vcf.gz

  """
}


// Run Expansion Hunter on each sample individually

process EXPANSION_HUNTER {

  time '40h'

  cpus 14

  memory '128 GB'

  publishDir("${params.results}/${sampleID}_results", mode: 'copy')

  input:
    // load sample metadata and CRAM
    tuple val( sampleID ), val( sex ), path( cram ), path( index ), val( trio ), val( sampleType )

  output:
    // channel for Expansion Hunter VCFs
    // save output VCFs
    tuple val( sampleID ), path( "${sampleID}_eh.vcf.gz" ), path( "${sampleID}_eh.vcf.gz.csi" ), emit: vcf

  script:
  """

  module load ${params.eh}
  module load ${params.bcftools}

  # The JSON file from convert_to_EH.sh specifies the regions to be analyzed (called a "variant catalog" by EH)
  # Must convert "M/F" sex designation to "male/female" for Expansion Hunter

  # define sex variable to fit Expansion Hunter syntax
  if [ ${sex} == "M" ]
  then
    eh_sex="male"
  else
    eh_sex="female"
  fi

  # run Expansion Hunter
  # zip and index the output VCF
  ExpansionHunter --reads ${cram} \
      --reference ${params.reference} \
      --variant-catalog ${params.eh_usats} \
      --output-prefix ${sampleID}_eh \
      --sex \${eh_sex} \
      --threads 12; bgzip -f ${sampleID}_eh.vcf; bcftools index ${sampleID}_eh.vcf.gz

  """
}


// Run GangSTR on each sample individually

process GANGSTR {

  time '10h'

  cpus 4

  memory '48 GB'

  publishDir("${params.results}/${sampleID}_results", mode: 'copy')

  input:
    // load sample metadata and CRAM
    tuple val( sampleID ), val( sex ), path( cram ), path( index ), val( trio ), val( sampleType )

  output:
    // channel for GangSTR VCFs
    // save output VCFs
    tuple val( sampleID ), path( "${sampleID}_gangstr.vcf.gz" ), path( "${sampleID}_gangstr.vcf.gz.csi" ), emit: vcf

  script:
  """
  module load ${params.gangstr}
  module load ${params.bcftools}

  # run GangSTR
  # zip and index the output VCF
  GangSTR \
      --bam ${cram} \
      --ref ${params.reference} \
      --regions ${params.gangstr_usats} \
      --out ${sampleID}_gangstr \
      --bam-samps ${sampleID} \
      --samp-sex ${sex} \
      --min-sample-reads 20 \
      --nonuniform \
      --frrweight 0 \
      --spanweight 0 \
      --flankweight 0; bgzip -f ${sampleID}_gangstr.vcf; bcftools index ${sampleID}_gangstr.vcf.gz
  """
}


// Run bedtools coverage on each sample to get spanning read depth

process BEDTOOLS_COVERAGE {

  time '3h'

  publishDir("${params.results}/${sampleID}_results", mode: 'copy')

  input:
    // load CRAM
    tuple val( sampleID ), val( sex ), path( cram ), path( index ), val( trio ), val( sampleType )

  output:
    // channel for bedtools coverage file
    // save coverage file
    tuple val( sampleID ), path( "${sampleID}_bedtools.txt" ), emit: table

  script:
  """
  module load ${params.bedtools}

  # Counts number of reads fully spanning microsatellite locus

  # Make column titles
  echo -e "name\tB_depth" > ${sampleID}_bedtools.txt

  # get spanning read depth at each locus
  bedtools coverage -sorted -g ${params.reference_index} -f 1.0 -a ${params.bedtools_usats} -b ${cram} | awk -e '{print \$4"\\t"\$5}' >> ${sampleID}_bedtools.txt

  # final tab-delimited output: [name] [depth]

  """
}


// Genotype variants from supplementary panel if present

process SUPP_VARIANTS {

  time '8h'
  cpus 2
  memory '32 GB'

  publishDir("${params.results}/${sampleID}_results", mode: 'copy')

  input:
    // load CRAM
    tuple val( sampleID ), val( sex ), path( cram ), path( index ), val( trio ), val( sampleType )

  output:
    // save variant calls
    tuple val( sampleID ), path( "${sampleID}.GATK.g.vcf.gz" ), path( "${sampleID}.GATK.g.vcf.gz.tbi" ), path( "${sampleID}.GATK.vcf.gz" ), path( "${sampleID}.GATK.vcf.gz.tbi" )

  script:
  """
  module load ${params.gatk}
  
  gatk --java-options "-Xmx20g" HaplotypeCaller \
     -I ${cram} -O ${sampleID}.GATK.g.vcf.gz \
     -R ${params.reference} \
     --intervals ${params.supp_panel_picard} \
     -ERC GVCF \
     -G StandardAnnotation \
     -G StandardHCAnnotation \
     -G AS_StandardAnnotation \
     --max-reads-per-alignment-start 0 \
     --interval-padding 50
   
   gatk GenotypeGVCFs \
     -R ${params.reference} \
     --intervals ${params.supp_panel_picard} \
     -V ${sampleID}.GATK.g.vcf.gz \
     -O ${sampleID}.GATK.vcf.gz
  
  """
}


// Check hybrid capture quality of supplementary panel

process SUPP_PICARD {

  time '8h'

  memory '32 GB'

  publishDir("${params.results}/${sampleID}_results", mode: 'copy')

  input:
    // load CRAM
    tuple val( sampleID ), val( sex ), path( cram ), path( index ), val( trio ), val( sampleType )

  output:
    // save Picard CollectHsMetrics report
    tuple val( sampleID ), path( "${sampleID}_picard_supp_panel.txt" ), path( "${sampleID}_picard_supp_probes.txt" )

  script:
  """
  module load ${params.jdk}
  # run Picard CollectHsMetrics on supplementary panel coordinates
  java -jar ${params.picard} CollectHsMetrics \\
        I=${cram} \\
        O=${sampleID}_picard_supp_panel.txt \\
        R=${params.reference} \\
        BAIT_INTERVALS=${params.supp_probes_picard} \\
        TARGET_INTERVALS=${params.supp_panel_picard}
        
  # run Picard CollectHsMetrics on supplementary panel coordinates
  # probes as targets
  java -jar ${params.picard} CollectHsMetrics \\
        I=${cram} \\
        O=${sampleID}_picard_supp_probes.txt \\
        R=${params.reference} \\
        BAIT_INTERVALS=${params.supp_probes_picard} \\
        TARGET_INTERVALS=${params.supp_probes_picard}
        
  """
}


/*

VCF DATA EXTRACTION
  
*/


/*
HipSTR fields: https://hipstr-tool.github.io/HipSTR/
  ID - name of locus (for joining)
  FORMAT/GT - indicates whether or not the locus has been genotyped (0/1/2 = genotyped, "." = missed)
    Since HipSTR does joint calling, the genotype number may be higher than 2 if more than two genotypes are seen amongst all     samples
  FORMAT/Q - quality of the call (from manual: "Posterior probability of unphased genotype")
  FORMAT/GB - From manual: "Base pair differences of genotype from reference"
    allele 1|allele 2 (ex: 0|0) (for genotyped alleles)
  FORMAT/DP - From manual: "Number of valid reads used for sample’s genotype"
  FORMAT/DSTUTTER - From manual: "Total number of reads with a stutter indel in the STR region"
  FORMAT/DFLANKINDEL - From manual: "Total number of reads with an indel in the regions flanking the STR"
  FORMAT/ALLREADS - gives number of spanning reads supporting a certain bp difference based on alignment (for all observed bp differences)
  FORMAT/MALLREADS - gives number of spanning reads supporting a certain bp difference based on maximum likelihood haplotype (for all observed bp differences)
  FORMAT/GLDIFF - From manual: "Difference in likelihood between the reported and next best genotypes"
    Note: HipSTR only uses spanning reads to genotype
    bp diff 1|# of supporting reads;bp diff 2|# of supporting reads (for all observed bp differences) (ex: -2|2;0|111;2|1)
*/


// Extract relevant fields from HipSTR VCF and save in table

process QUERY_HIPSTR {

  time '1h'

  input:
    // load HipSTR VCF
    tuple val( sampleID ), path( vcf )

  output:
    // channel for HipSTR data table
    tuple val( sampleID ), path( "${sampleID}_hipstr_table.txt" ), emit: table

  script:
  """
  module load ${params.bcftools}

  # Extract fields needed for quality filtering from HipSTR's VCF
  echo -e "name\tH_GT\tH_GB\tH_Q\tH_DP\tH_DSTUTTER\tH_DFLANKINDEL\tH_ALLREADS\tH_MALLREADS\tH_GLDIFF" > ${sampleID}_hipstr_table.txt
  bcftools query -f '%ID\t[%GT\t%GB\t%Q\t%DP\t%DSTUTTER\t%DFLANKINDEL\t%ALLREADS\t%MALLREADS\t%GLDIFF]\n' ${sampleID}_hipstr.vcf.gz >> ${sampleID}_hipstr_table.txt

  """
}


/*
ExpansionHunter fields: https://github.com/Illumina/ExpansionHunter and https://support-docs.illumina.com/SW/DRAGEN_v310/Content/SW/RepeatGenotyping.htm
  INFO/VARID - name of locus (for joining)
  INFO/REF - number of repeats in reference genome
  FORMAT/GT - indicates whether or not the locus has been genotyped (0/1/2 = genotyped, "." = missed)
  FORMAT/REPCN - the number of repeats in the allele (for genotyped alleles)
      allele 1/allele 2 (ex: 14/15)
  FORMAT/SO - From manual: "Type of reads that support the allele; can be SPANNING, FLANKING, or INREPEAT meaning that the reads span, flank, or are fully contained     in the repeat" (for genotyped alleles)
      allele 1/allele 2 (ex: SPANNING/SPANNING)
  FORMAT/ADSP - From manual: "Number of spanning reads consistent with the allele" (for genotyped alleles)
      allele 1/allele 2 (ex: 93/67)
  FORMAT/ADFL - From manual: "Number of flanking reads consistent with the allele" (for genotyped alleles)
    Note: "flanking" reads in EH are reads that end in the repeat on one side
    allele 1/allele 2
  FORMAT/ADIR - From manual: "Number of in-repeat reads consistent with the allele" (for genotyped alleles)
    Note: "in-repeat" reads are reads that fully consist of the repeat (the repeat is longer than the read length)
    allele 1/allele 2
  FORMAT/LC - the average coverage of the locus
*/

// Extract relevant fields from ExpansionHunter VCF

process QUERY_EH {

  time '1h'

  input:
    // load ExpansionHunter VCF
    tuple val( sampleID ), path( vcf ), path( index )

  output:
    // channel for ExpansionHunter data table
    tuple val( sampleID ), path( "${sampleID}_eh_table.txt" ), emit: table

  script:
  """
  module load ${params.bcftools}

  # Extracts fields needed for quality filtering from ExpansionHunter's VCF
  echo -e "name\tE_REF\tE_GT\tE_REPCN\tE_SO\tE_ADSP\tE_ADFL\tE_ADIR\tE_LC" > ${sampleID}_eh_table.txt
  bcftools query -f '%INFO/VARID\t%INFO/REF\t[%GT\t%REPCN\t%SO\t%ADSP\t%ADFL\t%ADIR\t%LC]\n' ${sampleID}_eh.vcf.gz >> ${sampleID}_eh_table.txt

  """
}


/*
GangSTR fields: https://github.com/gymreklab/GangSTR
  CHR and POS - chromosome and start position (for joining)
  INFO/REF - number of repeats in reference genome
  FORMAT/GT - indicates whether or not the locus has been genotyped (0/1/2 = genotyped, "." = missed)
  FORMAT/Q - quality score
  FORMAT/REPCN - From manual: "Genotype given in number of copies of the repeat motif" (for genotyped alleles)
      allele 1,allele 2 (ex: 13,14)
  FORMAT/REPCI - From manual: "95% Confidence interval for each allele based on bootstrapping"
  FORMAT/DP - From manual: "Read Depth (number of informative reads)"
  FORMAT/RC - number of reads in each class (enclosing, spanning, FRR, flanking)
    # of enclosing reads,# of spanning reads,#number of FRR reads,# of flanking reads
    Note: see GangSTR paper for read type definitions
  FORMAT/ENCLREADS - gives number of spanning reads supporting a certain repeat count (from manual: "Summary of reads in enclosing class in | separated key-value pairs. Keys are number of copies and values show number of reads with that many copies.") (for all observed bp differences)
      Note: "enclosing" reads in GangSTR are defined the same way as "spanning" reads in HipSTR and ExpansionHunter 
      repeat count 1,# of supporting reads|repeat count 2,# of supporting reads (ex: 6,2|13,39|14,38)
  FORMAT/FLNKREADS - gives number of flanking reads supporting a certain repeat count (from manual: "Summary of reads in flanking class in | separated key-value pairs. Keys are number of copies and values show number of reads with that many copies.")
      Note: "flanking" reads in GangSTR are reads in which the repeat is at the end of the read.
      repeat count 1,# of supporting reads|repeat count 2,# of supporting reads
*/

// Extract relevant fields from GangSTR VCF

process QUERY_GANGSTR {

  time '2h'

  input:
    tuple val( sampleID ), path( vcf ), path( index )

  output:
    tuple val( sampleID ), path( "${sampleID}_gangstr_table.txt" ), emit: table

  script:
  """
  module load ${params.bcftools}

  # Extracts fields needed for quality filtering from GangSTR's VCF
  echo -e "chr\tstart\tG_REF\tG_GT\tG_Q\tG_REPCN\tG_REPCI\tG_DP\tG_RC\tG_ENCLREADS\tG_FLNKREADS" > ${sampleID}_gangstr_table.txt
  bcftools query -f '%CHROM\t%POS\t%INFO/REF\t[%GT\t%Q\t%REPCN\t%REPCI\t%DP\t%RC\t%ENCLREADS\t%FLNKREADS]\n' ${sampleID}_gangstr.vcf.gz >> ${sampleID}_gangstr_table.txt

  """
}


// Join panel metadata with data tables from HipSTR, Expansion Hunter, GangSTR, and bedtools

process JOIN_CALLS {

  time '1h'

  input:
    // load data tables from HipSTR, Expansion Hunter, GangSTR, and bedtools
    tuple val( sampleID ), val( sex ), path( cram ), path( index ), val( trio ), val( sampleType ), path( hipstr ), path( expansionhunter ), path( gangstr ), path( bedtools )

  output:
    // channel for sample data frames
    path( "${sampleID}.rds" ), emit: rds

  script:
  """
  module load ${params.r}
  Rscript ${params.script_path}/join_str_calls_nextflow.R ${params.panel} ${sampleID} ${sex} ${trio} ${sampleType} ${hipstr} ${expansionhunter} ${gangstr} ${bedtools}

  """
}


// Combine all sample data frames

process JOIN_SAMPLES {

  time '2h'
  memory '32 GB'

  publishDir("${params.results}", mode: 'copy')

  input:
    // load all sample data frames
    path( RDS )

  output:
    // save final R dataframe
    path( "${params.outputbasename}.rds" )

  script:
  """
  module load ${params.r}
  # save sample TSV to results directory
  cp ${params.samples} ${params.results}
  # save configuration file to results directory
  cp `echo ${workflow.configFiles} | sed -e 's/\\[//' -e 's/\\]//'` ${params.results}

  # arguments: [RDS files] [basename for joint dataframe]
  # Nextflow automatically outputs all RDS files in space-separated list
  Rscript ${params.script_path}/join_str_samples_nextflow.R ${RDS} ${params.outputbasename}

  """
}
