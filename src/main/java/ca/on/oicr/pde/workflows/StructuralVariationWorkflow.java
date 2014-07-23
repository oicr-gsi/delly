package ca.on.oicr.pde.workflows;

import ca.on.oicr.pde.utilities.workflows.OicrWorkflow;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import net.sourceforge.seqware.common.util.Log;
import net.sourceforge.seqware.pipeline.workflowV2.model.Job;
import net.sourceforge.seqware.pipeline.workflowV2.model.SqwFile;
import net.sourceforge.seqware.pipeline.workflowV2.model.Workflow;

public class StructuralVariationWorkflow extends OicrWorkflow {

    private String dataDir;
    private String samtools;
    private String picard;
    private String delly;
    private String tabixDir;
    private String vcftoolsDir;
    private String queue;
    private String sampleName;
    private String excludeList;
    private String refFasta;
    private String mappingQuality;
    private Boolean manualOutput;
    private String[] inputBamFiles;
    private Map<String, SqwFile> inputFiles;
    private final String[] DELLY_TYPES = {"DEL", "DUP", "INV", "TRA"};  // Subject to change in future DELLY versions?
    private final String DEDUP_BAM_SUFFIX  = ".dupmarked.bam";
    private final String FILTR_BAM_SUFFIX  = ".dupmarked.filtered.bam";
    private final String VCF_MERGED_SUFFIX = ".delly.merged.vcf";
    
    

    public StructuralVariationWorkflow() {
        super();
        this.sampleName = "NA";
        this.inputFiles = new HashMap<String, SqwFile>();
    }

    /**
     * getInputFile will return a path to an input SqwFile stored in the private
     * inputFiles HashMap
     *
     * @param name
     * @return path to a SqwlFile stored in inputFiles HashMap
     */
    private SqwFile getInputFile(String name) {
        return inputFiles.get(name);
    }

    /**
     * createInputFile is a service function that would create a new SqwFile,
     * set it's type to 'Input' and add the path to the private inputFile
     * HashMap
     *
     * @param name
     * @param workingPath
     * @return Input SqwFile that would be identified by name variable passed to
     * this function
     */
    protected SqwFile createInputFile(String name, String workingPath) {
        SqwFile file = new SqwFile();
        file.setForceCopy(false);
        file.setIsInput(true);
        file.setSourcePath(workingPath);
        inputFiles.put(name, file);
        return file;
    }

    @Override
    public Map<String, SqwFile> setupFiles() {

        try {
            // TODO setup picard, delly and samtools
            if (getProperty("ref_fasta") == null) {
                Logger.getLogger(StructuralVariationWorkflow.class.getName()).log(Level.SEVERE, "The ref_fasta param was null! You need reference fasta file to run delly!");
                return (null);
            } else {
                this.refFasta = getProperty("ref_fasta");
            }

            if (getProperty("exclude_list") == null) {
                Logger.getLogger(StructuralVariationWorkflow.class.getName()).log(Level.SEVERE, "The exclude-list param was null! You need exclude list to run delly!");
                return (null);
            } else {
                this.excludeList = getProperty("exclude_list");
            }

            if (getProperty("sample_name") != null) {
                String sampleTmp = getProperty("sample_name").toString();
                if (!sampleTmp.isEmpty()) {
                    this.sampleName = sampleTmp;
                }
            }

            this.mappingQuality = getProperty("mapping_quality");
            if (null == this.mappingQuality || this.mappingQuality.isEmpty()) {
                this.mappingQuality = "0";
            }

            /**
             * Assume that if we have a tool path starting with "/" it would
             * contain /bin directory and we just append this value to the
             * WorkflowBaseDir, otherwise construct the path as if we only have
             * the tools directory name under bin
             */
            if (getProperty("samtools") == null) {
                Logger.getLogger(StructuralVariationWorkflow.class.getName()).log(Level.SEVERE, "The samtools param was null! We need samtools!");
                return (null);
            } else {
                this.samtools = getProperty("samtools").startsWith("/") ? getWorkflowBaseDir() + getProperty("samtools") : getWorkflowBaseDir() + "/bin/" + getProperty("samtools");
            }

            if (getProperty("picard") == null) {
                Logger.getLogger(StructuralVariationWorkflow.class.getName()).log(Level.SEVERE, "The picard param was null! We need picard!");
                return (null);
            } else {
                this.picard = getProperty("picard").startsWith("/") ? getWorkflowBaseDir() + getProperty("picard") : getWorkflowBaseDir() + "/bin/" + getProperty("picard");
            }

            if (getProperty("delly") == null) {
                Logger.getLogger(StructuralVariationWorkflow.class.getName()).log(Level.SEVERE, "The delly param was null! We need delly!");
                return (null);
            } else {
                this.delly = getProperty("delly").startsWith("/") ? getWorkflowBaseDir() + getProperty("delly")
                        : getWorkflowBaseDir() + "/bin/" + getProperty("delly") + "/" + getProperty("delly") + "_linux_x86_64bit";
            }
            
            if (getProperty("vcftools") == null) {
                Logger.getLogger(StructuralVariationWorkflow.class.getName()).log(Level.SEVERE, "vcftools is not set, we need it to call vcftools correctly");
                return (null);
            } else {
                this.vcftoolsDir = getProperty("vcftools").startsWith("/") ? getWorkflowBaseDir() + getProperty("vcftools")
                        : getWorkflowBaseDir() + "/bin/" + getProperty("vcftools");
            }

            if (getProperty("tabix") == null) {
                Logger.getLogger(StructuralVariationWorkflow.class.getName()).log(Level.SEVERE, "tabix is not set, we need it to call tabix correctly");
                return (null);
            } else {
                this.tabixDir = getProperty("tabix").startsWith("/") ? getWorkflowBaseDir() + getProperty("tabix")
                        : getWorkflowBaseDir() + "/bin/" + getProperty("tabix");
            }

            // INPUTS (bam files and their ids, multiple files WILL BE MERGED)
            if (getProperty("input_bams") == null) {
                Logger.getLogger(StructuralVariationWorkflow.class.getName()).log(Level.SEVERE, "The input_bams param was null! We need some bam files to work on!");
                return (null);
            } else {
                inputBamFiles = getProperty("input_bams").split(",");
            }

            manualOutput = Boolean.valueOf(getOptionalProperty("manual_output", "false"));
            queue = getOptionalProperty("queue", "");

            // iterate over inputs
            int fileIndex = 0;
            for (String filePath : inputBamFiles) {
                fileIndex++;
                Log.stdout("CREATING FILE: bam_inputs_" + fileIndex);
                SqwFile file = this.createInputFile("bam_inputs_" + fileIndex, filePath);
                file.setType("application/bam");
                file.setIsInput(true);
            }

        } catch (Exception e) {
            Logger.getLogger(StructuralVariationWorkflow.class.getName()).log(Level.SEVERE, null, e);
        }

        return this.getFiles();
    }

    @Override
    public void setupDirectory() {
        try {
            String datadir = getProperty("data_dir");
            if (null == datadir) {
                datadir = "data/";
            }
            this.addDirectory(datadir);
            this.dataDir = datadir.endsWith("/") ? datadir : datadir + "/";
        } catch (Exception e) {
            Logger.getLogger(StructuralVariationWorkflow.class.getName()).log(Level.WARNING, null, e);
        }

    }

    @Override
    public void buildWorkflow() {
        try {
            Workflow workflow = this.getWorkflow();
            Job bamJob = workflow.createBashJob("aggregate_bam_input");
            bamJob.setCommand("sleep 1");
            for (int i = 1; i <= inputBamFiles.length; i++) {
                bamJob.addFile(getInputFile("bam_inputs_" + i));
            }
            // Picard job
            Job picardJob = workflow.createBashJob("picard_deduplicate");
            picardJob.setCommand(
                    getWorkflowBaseDir() + "/bin/" + getProperty("bundled_jre") + "/bin/java"
                    + " -Xmx" + getProperty("picard_memory") + "M"
                    + " -jar " + this.picard + "/MarkDuplicates.jar"
                    + " TMP_DIR=" + this.dataDir + "picardTmp"
                    + " ASSUME_SORTED=true"
                    + " VALIDATION_STRINGENCY=LENIENT"
                    + " OUTPUT=" + this.dataDir + this.sampleName + DEDUP_BAM_SUFFIX
                    + " METRICS_FILE=" + this.dataDir + this.sampleName + ".bam.mmm");
            /* 
             + "REMOVE_DUPLICATES=true "
             + "CREATE_INDEX=true "); */

            for (int f = 0; f < this.inputBamFiles.length; f++) {
                picardJob.getCommand().addArgument(" INPUT=" + this.inputBamFiles[f]);
            }
            picardJob.setMaxMemory("16000");
            picardJob.addParent(bamJob);
            if (!this.queue.isEmpty()) {
                picardJob.setQueue(this.queue);
            }

            // Filtering job - 1280 means Remove Reads which are:
            // * not primary alignment
            // * are PCR or optical duplicates
            Job filterJob = workflow.createBashJob("samtools_filtering");
            filterJob.setCommand(this.samtools
                    + "/samtools view -b -F 1280 "
                    + this.dataDir + this.sampleName + DEDUP_BAM_SUFFIX
                    + " > " + this.dataDir + this.sampleName + FILTR_BAM_SUFFIX);

            filterJob.setMaxMemory("3000");
            filterJob.addParent(picardJob);
            if (!this.queue.isEmpty()) {
                filterJob.setQueue(this.queue);
            }

            // Indexing job - we needs this since DELLY uses index file
            Job indexJob = workflow.createBashJob("samtools_indexed");
            indexJob.setCommand(this.samtools
                    + "/samtools index "
                    + this.dataDir + this.sampleName + FILTR_BAM_SUFFIX);

            indexJob.setMaxMemory("3000");
            indexJob.addParent(filterJob);
            if (!this.queue.isEmpty()) {
                indexJob.setQueue(this.queue);
            }

            // DELLY jobs (4 different types)
            List <String> indexedVcfs = new ArrayList<String>();
            List <Job>    tabixJobs   = new ArrayList<Job>();
            
            for (String type : this.DELLY_TYPES) {
                Job dellyJob = workflow.createBashJob("delly_job_" + type);
                String outputFile = this.dataDir + this.sampleName + "." + type + ".vcf";
                dellyJob.setCommand(this.delly
                        + " -t " + type
                        + " -x " + this.excludeList
                        + " -o " + outputFile
                        + " -q " + this.mappingQuality
                        + " -g " + this.refFasta + " "
                        + this.dataDir + this.sampleName + FILTR_BAM_SUFFIX);
                dellyJob.setMaxMemory(getProperty("delly_memory"));
                dellyJob.addParent(indexJob);
                if (!this.queue.isEmpty()) {
                    dellyJob.setQueue(this.queue);
                }
                SqwFile dellyCalls = this.createOutputFile(outputFile, "text/vcf", this.manualOutput);
                dellyJob.addFile(dellyCalls);

                // tabix job for indexing (needed for merging by vcftools)
                Job tabixJob = this.getWorkflow().createBashJob("bgzip_idx_" + type);
                tabixJob.setCommand(this.tabixDir + "/bgzip -c "
                        + outputFile + " > "
                        + outputFile + ".gz;"
                        + this.tabixDir + "/tabix -p vcf "
                        + outputFile + ".gz");
                indexedVcfs.add(outputFile + ".gz");
                tabixJob.setMaxMemory("2000");
                tabixJob.addParent(dellyJob);
                if (!this.queue.isEmpty()) {
                    tabixJob.setQueue(this.queue);
                }
                tabixJobs.add(tabixJob);
                
            }

            // Concatenate all indexed vcf.gz files
            String mergeThese = "";
            for(String vcf : indexedVcfs) {
                mergeThese = mergeThese.concat(" " + vcf);
            }
            
            // Merging (vcftoolsDir) Job - make one file
            Job mergeJob = this.getWorkflow().createBashJob("vcf_merge");
            mergeJob.setCommand(this.vcftoolsDir + "/bin/vcf-merge"
                + mergeThese + " > "
                + this.dataDir + this.sampleName + VCF_MERGED_SUFFIX);
            
            SqwFile dellyVcf = this.createOutputFile(this.dataDir + this.sampleName + VCF_MERGED_SUFFIX, "text/vcf", this.manualOutput);
            mergeJob.setMaxMemory("3000");
            mergeJob.addFile(dellyVcf);
            for (Job tj : tabixJobs) {
                mergeJob.addParent(tj);
            }
            if (!this.queue.isEmpty()) {
                    mergeJob.setQueue(this.queue);
            }
            
        
        } catch (Exception e) {
            Logger.getLogger(getClass().getName()).log(Level.SEVERE, null, e);
        }

    }
}
