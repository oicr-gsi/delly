package ca.on.oicr.pde.workflows;

import ca.on.oicr.pde.utilities.workflows.OicrWorkflow;
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
    private String queue;
    private String sampleName;
    private String excludeList;
    private String refFasta;
    private String mappingQuality;
    private Boolean manualOutput;
    private String[] inputBamFiles;
    private Map<String, SqwFile> inputFiles;
    private final String[] DELLY_TYPES = { "DEL", "DUP", "INV", "TRA" };

    public StructuralVariationWorkflow() {
        super();
        this.sampleName = "NA";
        this.manualOutput = null;
        this.inputBamFiles = null;
    }

    private SqwFile getInputFile(String name) {
        return inputFiles.get(name);
    }

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

            if (getProperty("exclude-list") == null) {
                Logger.getLogger(StructuralVariationWorkflow.class.getName()).log(Level.SEVERE, "The exclude-list param was null! You need exclude list to run delly!");
                return (null);
            } else {
                this.excludeList = getProperty("exclude-list");
            }

            if (getProperty("sample_name") != null) {
                String sampleTmp = getProperty("sample_name").toString();
                if (!sampleTmp.isEmpty()) {
                    this.sampleName = sampleTmp;
                }
            }
            
            this.mappingQuality = getProperty("mapping_quality");
            if ( null == this.mappingQuality )
                 this.mappingQuality = "0";

            if (getProperty("samtools") == null) {
                Logger.getLogger(StructuralVariationWorkflow.class.getName()).log(Level.SEVERE, "The samtools param was null! We need samtools!");
                return (null);
            } else {
                this.samtools = getProperty("samtools").startsWith("/") ? getWorkflowBaseDir() + "bin/" + getProperty("samtools") : getWorkflowBaseDir() + "/" + getProperty("samtools");
            }

            if (getProperty("picard") == null) {
                Logger.getLogger(StructuralVariationWorkflow.class.getName()).log(Level.SEVERE, "The picard param was null! We need picard!");
                return (null);
            } else {
                this.picard = getProperty("picard").startsWith("/") ? getWorkflowBaseDir() + "bin/" + getProperty("picard") : getWorkflowBaseDir() + "/" + getProperty("picard");
            }

            if (getProperty("delly") == null) {
                Logger.getLogger(StructuralVariationWorkflow.class.getName()).log(Level.SEVERE, "The delly param was null! We need delly!");
                return (null);
            } else {
                this.delly = getProperty("delly").startsWith("/") ? getWorkflowBaseDir() + "bin/" + getProperty("delly") : getWorkflowBaseDir() + "/" + getProperty("delly");
            }

            // INPUTS (bam files and their ids, multiple files will be merged)
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
            //Picard job
            Job picardJob = workflow.createBashJob("pcard_deduplicate");
            picardJob.setCommand(
                    getWorkflowBaseDir() + "/bin/" + getProperty("bundled_jre") + "/bin/java "
                    + "-Xmx" + getProperty("picard_memory") + "M "
                    + "-jar " + this.picard + "/MarkDuplicates.jar "
                    + "TMP_DIR=" + this.dataDir + "picardTmp "
                    + "OUTPUT=" + this.dataDir + this.sampleName + "deduplicated.bam "
                    + "METRICS_FILE=" + this.dataDir + this.sampleName + ".bam.mmm ");
            // + "ASSUME_SORTED=true "
            // + "VALIDATION_STRINGENCY=LENIENT "
            // + "REMOVE_DUPLICATES=true "
            // + "CREATE_INDEX=true ");

            for (int f = 0; f < this.inputBamFiles.length; f++) {
                picardJob.getCommand().addArgument("INPUT=" + this.inputBamFiles[f]); // this.dataDir + filtered_bam);
            }
            picardJob.setMaxMemory("16000");
            picardJob.addParent(bamJob);
            if (!this.queue.isEmpty()) {
                picardJob.setQueue(this.queue);
            }

            //Samtools job
            Job filter_job = workflow.createBashJob("samtools_filtering");
            filter_job.setCommand(this.samtools
                    + "/samtools view -b -F 1280 "
                    + this.dataDir + this.sampleName + "deduplicated.bam "
                    + "> " + this.dataDir + this.sampleName + "filtered.bam");

            filter_job.setMaxMemory("3000");
            filter_job.addParent(picardJob);
            if (!this.queue.isEmpty()) {
                filter_job.setQueue(this.queue);
            }

            //DELLY jobs
            for(String type : this.DELLY_TYPES) {
            Job delly_job = workflow.createBashJob("delly_job");
            String outputFile = this.dataDir + this.sampleName + "." + type + ".vcf";
            delly_job.setCommand(this.delly
                    + "/delly -t " + type 
                    + " -x " + this.excludeList
                    + " -o " + outputFile
                    + " -q " + this.mappingQuality
                    + " -g " + this.refFasta + " "
                    + this.dataDir + this.sampleName + "filtered.bam");
            delly_job.setMaxMemory("6000");
            delly_job.addParent(picardJob);
            if (!this.queue.isEmpty()) {
                delly_job.setQueue(this.queue);
            }
            SqwFile dellyCalls = this.createOutputFile(outputFile, "text/vcf", this.manualOutput);
            delly_job.addFile(dellyCalls);
            }

        } catch (Exception e) {
            Logger.getLogger(getClass().getName()).log(Level.SEVERE, null, e);
        }

    }
}
