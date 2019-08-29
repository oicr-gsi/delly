version 1.0

workflow dellyWorkflow {
input {
    # If we are in somatic mode, normal file follows tumor file in the input array
    Array[File] inputBams
    String      sampleID
}

# If we see more than one (two) bams switch to somatic mode
scatter (f in inputBams) {
    call dupmarkBam { input: inputBam = f}
}

String callType = if length(inputBams) == 1 then "unpaired" else "somatic"

scatter (m in ["DEL", "DUP", "INV", "TRA"]) {
 call runDelly { input: inBams = dupmarkBam.outputBam, inBai = dupmarkBam.outputBai, dellyMode = m, callType = callType, sampleName = sampleID }
 if (runDelly.validOutput) {
     File delly_out = runDelly.outVcf
     File delly_idx = runDelly.outTbi
 }
}

Array[File?] delly_out_maybes = delly_out
Array[File?] delly_tbi_maybes = delly_idx
Array[File] delly_out_valids  = select_all(delly_out_maybes)
Array[File] delly_tbi_valids  = select_all(delly_tbi_maybes)

# Go on with merging and zipping/indexing
call mergeAndZip { input: inputVcfs = delly_out_valids, inputTbis = delly_tbi_valids, sampleName = sampleID, callType = callType}

}

# ==========================================
#  TASK 1 of 3: imark duplicates with picard
# ==========================================
task dupmarkBam {
input {
	File   inputBam
        Int?   jobMemory  = 20
        Int?   javaMemory = 12
        String? modules = "picard/1.72" 
}

command <<<
 java -Xmx~{javaMemory}G -jar $PICARDROOT/MarkDuplicates.jar \
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
#  TASK 2 of 3: index bam file
# ================================
task runDelly {
input { 
        # We may have 1 or 2 files here, tumor always first
        Array[File] inBams
        Array[File] inBai
        String dellyMode
        String sampleName
        String? excludeList = "/.mounts/labs/PDE/data/reference/hg19/delly/human.hg19.excl.tsv"
        String? refFasta = "/.mounts/labs/PDE/data/reference/hg19_random/fasta/UCSC/hg19_random.fa"
        String? callType = "unpaired"
        String? modules = "delly/0.8.1 bcftools-1.7/1.7 tabix/0.2.6"
        Int? mappingQuality = 30
        Int? jobMemory = 10
}

command <<<
VALID_TAG="false"
delly call -t ~{dellyMode} \
      -x ~{excludeList} \
      -o "~{sampleName}.~{dellyMode}.~{callType}.bcf" \
      -q ~{mappingQuality} \
      -g ~{refFasta} \
         ~{sep=' ' inBams}
bcftools view "~{sampleName}.~{dellyMode}.~{callType}.bcf" | bgzip -c > "~{sampleName}.~{dellyMode}.~{callType}.vcf.gz"

if [ -e "~{sampleName}.~{dellyMode}.~{callType}.vcf.gz" ]; then
   bcftools view "~{sampleName}.~{dellyMode}.~{callType}.bcf" | bgzip -c > "~{sampleName}.~{dellyMode}.~{callType}.vcf.gz"
   tabix -p vcf "~{sampleName}.~{dellyMode}.~{callType}.vcf.gz"
   VALID_TAG="true"
else
   touch "~{sampleName}.~{dellyMode}.~{callType}.vcf.gz"
   touch "~{sampleName}.~{dellyMode}.~{callType}.vcf.gz.tbi"
fi
echo $VALID_TAG 1>&2
>>>

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
}

output {
  Boolean validOutput = read_boolean(stderr())
  File outVcf = "~{sampleName}.~{dellyMode}.~{callType}.vcf.gz"
  File outTbi = "~{sampleName}.~{dellyMode}.~{callType}.vcf.gz.tbi"
}
}


# =====================================================
#  TASK 3 of 3: merge vcf files and bgzip/index results
# =====================================================
task mergeAndZip {
input {
        Array[File] inputVcfs
        Array[File] inputTbis
        String sampleName
        String? callType = "unpaired"
        String? modules = "vcftools/0.1.10 tabix/0.2.6"
	Int? jobMemory = 10
}

command <<<
       vcf-merge ~{sep=' ' inputVcfs} | bgzip -c > "~{sampleName}.~{callType}.delly.merged.vcf.gz"
       tabix -p vcf "~{sampleName}.~{callType}.delly.merged.vcf.gz"
>>>

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
}

output {
  File dellyMergedVcf       = "~{sampleName}.~{callType}.delly.merged.vcf.gz"
  File delyMergedTabixIndex = "~{sampleName}.~{callType}.delly.merged.vcf.gz.tbi"
}
}

