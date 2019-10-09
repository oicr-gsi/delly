version 1.0

workflow delly {
input {
    # If we are in somatic mode, normal file follows tumor file in the input array
    Array[File]+ inputBams
    String? outputFileNamePrefix = ""
}

String? sampleID = if outputFileNamePrefix=="" then basename(inputBams[0], ".bam") else outputFileNamePrefix
# If we see more than one (two) bams switch to somatic mode
scatter (f in inputBams) {
    call dupmarkBam { input: inputBam = f}
}

String callType = if length(inputBams) == 1 then "unpaired" else "somatic"

scatter (m in ["DEL", "DUP", "INV", "INS", "BND"]) {
 call runDelly { input: inBams = dupmarkBam.outputBam, inBai = dupmarkBam.outputBai, dellyMode = m, callType = callType, sampleName = sampleID }
}

# Go on with merging and zipping/indexing
call mergeAndZip as mergeAndZipALL { input: inputVcfs = select_all(runDelly.outVcf), inputTbis = select_all(runDelly.outTbi), sampleName = sampleID, callType = callType}

# Go on with processing somatic - filtered files
if (callType == "somatic") {
 call mergeAndZip as mergeAndZipFiltered { input: inputVcfs = select_all(runDelly.outVcf_filtered), inputTbis = select_all(runDelly.outTbi_filtered), sampleName = sampleID, callType = callType, prefix = "_filtered"}
}

meta {
  author: "Peter Ruzanov"
  email: "peter.ruzanov@oicr.on.ca"
  description: "StructuralVariation 2.0"
}

output {
  File? mergedIndex = mergeAndZipALL.dellyMergedTabixIndex
  File? mergedVcf   = mergeAndZipALL.dellyMergedVcf
  File? mergedFilteredIndex = mergeAndZipFiltered.dellyMergedTabixIndex
  File? mergedFilteredVcf   = mergeAndZipFiltered.dellyMergedVcf
}

}

# ==========================================
#  TASK 1 of 3: mark duplicates with picard
# ==========================================
task dupmarkBam {
input {
	File   inputBam
        Int?   jobMemory  = 20
        String? modules = "java/8 picard/2.19.2" 
}

parameter_meta {
 inputBam: "Input .bam file"
 jobMemory: "memory allocated for Job"
 modules: "Names and versions of modules for picard-tools and java"
}

command <<<
 java -Xmx~{jobMemory-4}G -jar $PICARD_ROOT/picard.jar MarkDuplicates \
                              TMP_DIR=picardTmp \
                              ASSUME_SORTED=true \
                              VALIDATION_STRINGENCY=LENIENT \
                              OUTPUT="~{basename(inputBam, '.bam')}_dupmarked.bam" \
                              INPUT=~{inputBam} \
                              CREATE_INDEX=true \
                              METRICS_FILE="~{basename(inputBam)}.mmm"
>>>

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
}

output {
  File outputBam = "~{basename(inputBam, '.bam')}_dupmarked.bam"
  File outputBai = "~{basename(inputBam, '.bam')}_dupmarked.bai"
}
}

# ================================
#  TASK 2 of 3: run delly
# ================================
task runDelly {
input { 
        # We may have 1 or 2 files here, tumor always first
        Array[File]+ inBams
        Array[File]+ inBai
        String dellyMode
        String? sampleName = "SAMPLE"
        String excludeList
        String? refFasta = "$HG19_ROOT/hg19_random.fa"
        String? callType = "unpaired"
        String? modules = "delly/0.8.1 bcftools/1.9 tabix/0.2.6 hg19/p13"
        Int? mappingQuality = 30
        Int? jobMemory = 10
}

parameter_meta {
 inBams: "Input .bam files"
 inBai: "Input .bai files"
 dellyMode: "Mode specifying type of call"
 sampleName: "Normally passed from workflow block, prefix for making output files"
 excludeList: "List of regions to exclude (telomeres and centromeres)"
 refFasta: "reference assembly file"
 callType: "unpaired or somatic"
 mappingQuality: "defines quality threshold for reads to use in calling SVs"
 jobMemory: "memory allocated for Job"
 modules: "Names and versions of modules for picard-tools and java"
}

command <<<
delly call -t ~{dellyMode} \
      -x ~{excludeList} \
      -o "~{sampleName}.~{dellyMode}.~{callType}.bcf" \
      -q ~{mappingQuality} \
      -g ~{refFasta} \
         ~{sep=' ' inBams}
if [ "~{callType}" == "somatic" ]; then
   echo "Somatic mode requested, will run delly filtering for somatic SVs"
   bcftools view "~{sampleName}.~{dellyMode}.~{callType}.bcf" | grep ^# | tail -n 1 | \
   sed 's/.*FORMAT\t//' | awk -F "\t" '{print $1"\ttumor";print $2"\tcontrol"}' > samples.tsv
   delly filter -f somatic -o "~{sampleName}.~{dellyMode}.~{callType}_filtered.bcf" -s samples.tsv \
                              "~{sampleName}.~{dellyMode}.~{callType}.bcf"
   bcftools view "~{sampleName}.~{dellyMode}.~{callType}_filtered.bcf" | \
   bgzip -c > "~{sampleName}.~{dellyMode}.~{callType}_filtered.vcf.gz"
fi

bcftools view "~{sampleName}.~{dellyMode}.~{callType}.bcf" | bgzip -c > "~{sampleName}.~{dellyMode}.~{callType}.vcf.gz"

if [ -e "~{sampleName}.~{dellyMode}.~{callType}.vcf.gz" ]; then
   tabix -p vcf "~{sampleName}.~{dellyMode}.~{callType}.vcf.gz"
fi
if [ -e "~{sampleName}.~{dellyMode}.~{callType}_filtered.vcf.gz" ]; then
   tabix -p vcf "~{sampleName}.~{dellyMode}.~{callType}_filtered.vcf.gz"
fi
>>>

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
}

output {
  File? outVcf = "~{sampleName}.~{dellyMode}.~{callType}.vcf.gz"
  File? outTbi = "~{sampleName}.~{dellyMode}.~{callType}.vcf.gz.tbi"
  File? outVcf_filtered = "~{sampleName}.~{dellyMode}.~{callType}_filtered.vcf.gz"
  File? outTbi_filtered = "~{sampleName}.~{dellyMode}.~{callType}_filtered.vcf.gz.tbi"
}
}


# =====================================================
#  TASK 3 of 3: merge vcf files and bgzip/index results
# =====================================================
task mergeAndZip {
input {
        Array[File] inputVcfs
        Array[File] inputTbis
        String? sampleName = "SAMPLE"
        String? callType = "unpaired"
        String? modules = "vcftools/0.1.16 tabix/0.2.6"
        String? prefix = ""
	Int? jobMemory = 10
}

parameter_meta {
 inputVcfs: "Input .bam files"
 inputTbis: "Input .bai files"
 sampleName: "Normally passed from workflow block, prefix for making output files"
 callType: "unpaired or somatic"
 modules: "Names and versions of modules for picard-tools and java"
 prefix: "parameter to use when we need to append _filtered to the file's name"
 jobMemory: "memory allocated for Job"
}


command <<<
       vcf-concat ~{sep=' ' inputVcfs} | vcf-sort | bgzip -c > "~{sampleName}.~{callType}~{prefix}.delly.merged.vcf.gz"
       tabix -p vcf "~{sampleName}.~{callType}~{prefix}.delly.merged.vcf.gz"
>>>

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
}

output {
  File? dellyMergedVcf        = "~{sampleName}.~{callType}~{prefix}.delly.merged.vcf.gz"
  File? dellyMergedTabixIndex = "~{sampleName}.~{callType}~{prefix}.delly.merged.vcf.gz.tbi"
}
}

