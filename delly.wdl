version 1.0

struct GenomeResources {
    String rundelly_module
    String rundelly_fasta
    String rundelly_exclude_list
}

workflow delly {

input {
  # If we are in somatic mode, normal file follows tumor file in the input array
  File inputTumor
  File? inputNormal
  Boolean markdup = true
  String outputFileNamePrefix
  String reference
}

Map[String,GenomeResources] resources = {
  "hg19": {
    "rundelly_module": "delly/0.9.1 bcftools/1.9 tabix/0.2.6 hg19/p13 hg19-delly/1.0",
    "rundelly_fasta": "$HG19_ROOT/hg19_random.fa",
    "rundelly_exclude_list": "$HG19_DELLY_ROOT/human.hg19.excl.tsv"
  },
   "hg38": {
    "rundelly_module": "delly/0.9.1 bcftools/1.9 tabix/0.2.6 hg38/p12 hg38-delly/1.0",
    "rundelly_fasta": "$HG38_ROOT/hg38_random.fa",
    "rundelly_exclude_list": "$HG38_DELLY_ROOT/human.hg38.excl.tsv"
   }
}

Array[File] inputBams= select_all([inputTumor,inputNormal])
String sampleID = outputFileNamePrefix
String callType = if length(inputBams) == 1 then "unmatched" else "somatic"

# If we see more than one (two) bams switch to somatic mode
scatter (f in inputBams) { 
  call dupmarkBam { input: inputBam = f, dedup = if markdup then "dedup" else "nomark"}
} 

scatter (m in ["DEL", "DUP", "INV", "INS", "BND"]) {
  call runDelly { input: inBams = dupmarkBam.outputBam, inBai = dupmarkBam.outputBai, dellyMode = m, callType = callType, sampleName = sampleID, modules = resources [ reference ].rundelly_module, refFasta = resources [ reference ].rundelly_fasta, excludeList = resources [ reference ].rundelly_exclude_list}
}

# Go on with merging and zipping/indexing
call mergeAndZip as mergeAndZipALL { input: inputVcfs = select_all(runDelly.outVcf), inputTbis = select_all(runDelly.outTbi), sampleName = sampleID, callType = callType, prefix = "_all"}

# Go on with processing somatic - filtered files
if (callType == "somatic") {
 call mergeAndZip as mergeAndZipFiltered { input: inputVcfs = select_all(runDelly.outVcf_filtered), inputTbis = select_all(runDelly.outTbi_filtered), sampleName = sampleID, callType = callType, prefix = "_filtered"}
}

parameter_meta {
  inputTumor: "Tumor input .bam file."
  inputNormal: "Normal input .bam file."
  markdup: "A switch between marking duplicate reads and indexing with picard."
  outputFileNamePrefix: "Output prefix to be used with result files."
  reference: "the reference genome for input sample"
}

meta {
  author: "Peter Ruzanov"
  email: "peter.ruzanov@oicr.on.ca"
  description: "The Delly workflow produces a set of vcf files with different types of structural variant calls: Translocations, Deletions, Inversions and Duplications It uses .bam files as input. The below graph describes the process:\n![delly flowchart](docs/delly-wf.png)\n### Preprocessing\nThe expected inputs for the DELLY tool are aligned sequence (bam format), properly sorted and indexed, with marked duplicates. \n### Mark duplicates\nPicard Tools MarkDuplicates is used to flag reads as PCR or optical duplicates and is activated be default.  If providing bam files with duplicates marked, this can be disabled.\n```\n java -jar MarkDuplicates.jar\n INPUT=sample.bam\n OUTPUT=sample.dedup.bam    \n METRICS_FILE=sample.metrics\n```\n### Detect deletions\n```\ndelly\n-t DEL\n-x excludeList.tsv\n-o sample.jumpy.bam\n-q 0\n-g hn19.fa\nsample.bam\n```\n### Detect tandem duplications\n```\ndelly\n-t DUP\n-x excludeList.tsv\n-o sample.jumpy.bam\n-q 0\n-g hn19.fa\nsample.bam\n```\n### Detect inversions\n```\ndelly\n-t INV\n-x excludeList.tsv\n-o sample.jumpy.bam\n-q 0\n-g hn19.fa\nsample.bam\nDetect translocations\n```\n### Detecting translocations\n```\ndelly\n-t TRA\n-x excludeList.tsv\n-o sample.jumpy.bam\n-q 0\n-g hn19.fa\nsample.bam\n```\n### Post-processing\nEach DELLY tool produces several files, which will all need to be merged together after the chromosomes are finished processing. The output format is described on the DELLY webpage. The merging script may require a small parser to combine the output from multiple runs in together.\nMerge DELLY results with vcftools"
  dependencies: [
      {
        name: "picard/2.19.2",
        url: "https://master.dl.sourceforge.net/project/picard/picard-tools/1.89/picard-tools-1.89.zip"
      },
      {
        name: "java/8",
        url: "https://github.com/AdoptOpenJDK/openjdk8-upstream-binaries/releases/download/jdk8u222-b10/OpenJDK8U-jdk_x64_linux_8u222b10.tar.gz"
      },
      {
        name: "delly/0.9.1",
        url: "https://github.com/dellytools/delly/releases/download/v0.9.1/delly_v0.9.1_linux_x86_64bit"
      },
      {
        name: "bcftools/1.9",
        url: "https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2"
      },
      {
        name: "tabix/0.2.6",
        url: "https://sourceforge.net/projects/samtools/files/tabix/tabix-0.2.6.tar.bz2"
      },
      {
        name: "vcftools/0.1.16",
        url: "https://github.com/vcftools/vcftools/archive/v0.1.16.tar.gz"
      }
    ]
    output_meta: {
      mergedVcf: "vcf file containing all structural variant calls",
      mergedIndex: "tabix index of the vcf file containing all structural variant calls",
      mergedFilteredVcf: "filtered vcf file containing structural variant calls",
      mergedFilteredIndex: "tabix index of the filtered vcf file containing structural variant calls",
      mergedFilteredPassVcf: "filtered vcf file containing PASS structural variant calls",
      mergedFilteredPassIndex: "tabix index of the filtered vcf file containing PASS structural variant calls"
    }
}

output {
  File mergedIndex = mergeAndZipALL.dellyMergedTabixIndex
  File mergedVcf   = mergeAndZipALL.dellyMergedVcf
  File? mergedFilteredIndex = mergeAndZipFiltered.dellyMergedTabixIndex
  File? mergedFilteredVcf   = mergeAndZipFiltered.dellyMergedVcf
  File? mergedFilteredPassIndex = mergeAndZipFiltered.dellyMergedPassTabixIndex
  File? mergedFilteredPassVcf = mergeAndZipFiltered.dellyMergedPassVcf
}

}

# ==========================================
#  TASK 1 of 3: mark duplicates with picard
# ==========================================
task dupmarkBam {
input {
  File inputBam
  Int jobMemory = 20
  Int timeout   = 20
  String dedup = "dedup"
  String modules = "java/8 picard/2.19.2"
}

parameter_meta {
 inputBam: "Input .bam file"
 jobMemory: "memory allocated for Job"
 dedup: "A switch between marking duplicate reads and indexing with picard"
 modules: "Names and versions of modules for picard-tools and java"
 timeout: "Timeout in hours"
}

command <<<
 set -eu -o pipefail
 if [ "~{dedup}" == "dedup" ]; then
  java -Xmx~{jobMemory-8}G -jar $PICARD_ROOT/picard.jar MarkDuplicates \
                                TMP_DIR=picardTmp \
                                ASSUME_SORTED=true \
                                VALIDATION_STRINGENCY=LENIENT \
                                OUTPUT="~{basename(inputBam, '.bam')}_dupmarked.bam" \
                                INPUT=~{inputBam} \
                                CREATE_INDEX=true \
                                METRICS_FILE="~{basename(inputBam)}.mmm"
 else
  ln -s ~{inputBam} ~{basename(inputBam)}
  java -Xmx~{jobMemory-8}G -jar $PICARD_ROOT/picard.jar BuildBamIndex \
                              VALIDATION_STRINGENCY=LENIENT \
                              INPUT=~{basename(inputBam)} \
                              OUTPUT="~{basename(inputBam, '.bam')}.bai"
 fi
>>>

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  timeout: "~{timeout}"
} 

output {
  File outputBam = if "~{dedup}" == "dedup" then "~{basename(inputBam, '.bam')}_dupmarked.bam" else "~{basename(inputBam)}"
  File outputBai = if "~{dedup}" == "dedup" then "~{basename(inputBam, '.bam')}_dupmarked.bai" else "~{basename(inputBam, '.bam')}.bai"
}
}

# ================================
#  TASK 2 of 3: run delly
# ================================
task runDelly {
input { 
  Array[File]+ inBams
  Array[File]+ inBai
  String dellyMode
  String sampleName
  String excludeList
  String refFasta
  String callType = "unmatched"
  String modules
  Int mappingQuality = 30
  Int jobMemory = 16
  Int timeout = 20
}

parameter_meta {
 inBams: "Input .bam files"
 inBai: "Input .bai files"
 dellyMode: "Mode specifying type of call"
 sampleName: "Normally passed from workflow block, prefix for making output files"
 excludeList: "List of regions to exclude (telomeres and centromeres)"
 refFasta: "reference assembly file"
 callType: "unmatched or somatic"
 mappingQuality: "defines quality threshold for reads to use in calling SVs"
 jobMemory: "memory allocated for Job"
 timeout: "Timeout in hours"
 modules: "Names and versions of modules for picard-tools and java"
}

command <<<
set -eu -o pipefail
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
  timeout: "~{timeout}"
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
  String sampleName = "SAMPLE"
  String callType = "unmatched"
  String modules = "bcftools/1.9 vcftools/0.1.16 tabix/0.2.6"
  String prefix = ""
  Int variantSupport = 0
  Int jobMemory = 10
}

parameter_meta {
 inputVcfs: "Input .bam files"
 inputTbis: "Input .bai files"
 sampleName: "Normally passed from workflow block, prefix for making output files"
 callType: "unmatched or somatic"
 modules: "Names and versions of modules for picard-tools and java"
 prefix: "parameter to use when we need to append _filtered to the file's name"
 variantSupport: "Paired-end support for structural variants, in pairs. Default is 10"
 jobMemory: "memory allocated for Job"
}


command <<<
  set -eu -o pipefail
  vcf-concat ~{sep=' ' inputVcfs} | vcf-sort | bgzip -c > "~{sampleName}.~{callType}~{prefix}.vcf.gz"
  tabix -p vcf "~{sampleName}.~{callType}~{prefix}.vcf.gz"
  if [ -e ~{sampleName}.~{callType}_filtered.vcf.gz ]; then
    bcftools view -i "%FILTER='PASS' & INFO/PE>~{variantSupport}" ~{sampleName}.~{callType}~{prefix}.vcf.gz -Oz -o ~{sampleName}.~{callType}~{prefix}.pass.vcf.gz
    tabix -p vcf ~{sampleName}.~{callType}~{prefix}.pass.vcf.gz
  fi
>>>

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
}

output {
  File dellyMergedVcf        = "~{sampleName}.~{callType}~{prefix}.vcf.gz"
  File dellyMergedTabixIndex = "~{sampleName}.~{callType}~{prefix}.vcf.gz.tbi"
  File? dellyMergedPassVcf    = "~{sampleName}.~{callType}~{prefix}.pass.vcf.gz"
  File? dellyMergedPassTabixIndex =  "~{sampleName}.~{callType}~{prefix}.pass.vcf.gz.tbi"
}
}
