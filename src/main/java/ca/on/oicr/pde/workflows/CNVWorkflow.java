package ca.on.oicr.pde.workflows;

import ca.on.oicr.pde.utilities.workflows.OicrWorkflow;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import net.sourceforge.seqware.common.util.Log;
import net.sourceforge.seqware.pipeline.workflowV2.model.Job;
import net.sourceforge.seqware.pipeline.workflowV2.model.SqwFile;
import net.sourceforge.seqware.pipeline.workflowV2.model.Workflow;

public class CNVWorkflow extends OicrWorkflow {

    //Versions of tools
    private String bicseqVersion;
    private String freecVersion;
    private String samtoolsVersion;
    
    //FREEC
    private String chrLengthFile = "";
    private String templateType;
    private String freecVarCoeff = "";
            
    //References
    private boolean manualOutput;
    private boolean doSort = true;
    private String  queue;

    //Data
    private String[] normal;
    private String[] tumor;
    private String[] localInputNormalFiles;
    private String[] localInputTumorFiles;
    private String[] normalBases;
    private String[] tumorBases;
    private String[] chromData; // if we have this, do mutect split/merge that should be faster

    //Misc
    private String dataDir;
    private String tempDir;
    private int bicseqInterval;
    private int bicseqSpread;

    private boolean doCrosscheck = false;
    private static final String BICSEQ_I_DEFAULT = "150";
    private static final String BICSEQ_S_DEFAULT = "20";
    private static final String FREEC_CV_DEFAULT = "0.5";
    /**
     * 
     * 
     *
     */
    
    @Override
    public Map<String, SqwFile> setupFiles() {

        try {

            String nextProp = getProperty("input_files_normal");
            if (nextProp == null) {
                Logger.getLogger(CNVWorkflow.class.getName()).log(Level.SEVERE, "input_files_normal is not set, we need at least one bam file");
                return (null);
            } else {
                this.normal = nextProp.split(",");
            }

            nextProp = getProperty("force_crosscheck");
            if (nextProp == null) {
                Logger.getLogger(CNVWorkflow.class.getName()).log(Level.WARNING, "force_crosscheck is not set, will do only one normal vs N tumors or pairwise N vs T");
            } else {
                this.doCrosscheck = nextProp.isEmpty() || nextProp.equalsIgnoreCase("false") ? false : true;
            }

            nextProp = getProperty("input_files_tumor");
            if (nextProp == null) {
                Logger.getLogger(CNVWorkflow.class.getName()).log(Level.SEVERE, "input_files_tumor is not set, we need at least one bam file");
                return (null);
            } else {
                this.tumor = nextProp.split(",");
                if (this.tumor.length != this.normal.length) {
                    if (this.normal.length == 1) {
                        this.doCrosscheck = true;
                    } else {
                        Logger.getLogger(CNVWorkflow.class.getName()).log(Level.SEVERE, "numbers for normal and tumor bam files are not the same, crosscheck isn't forced - "
                                + "check your .ini file");
                        return (null);
                    }
                }
            }

            nextProp = getProperty("chromosome_length");
            if (nextProp == null) {
                Logger.getLogger(CNVWorkflow.class.getName()).log(Level.WARNING, "chromosome_length is not set, mutect won't split by chromosome, slower option");
                this.chromData = null;
            } else {
                this.chromData = nextProp.split(",");
            }

            if (getProperty("queue") == null) {
                Logger.getLogger(CNVWorkflow.class.getName()).log(Level.WARNING, "Queue not set, will run on a queue assigned by sge");
                this.queue = "";
                return (null);
            } else {
                this.queue = getProperty("queue");
            }

            nextProp = getProperty("manual_output");
            if (nextProp == null) {
                Logger.getLogger(CNVWorkflow.class.getName()).log(Level.WARNING, "manual_output is not set, will put the file into automatically generated dir");
            } else {
                this.manualOutput = nextProp.isEmpty() || nextProp.equalsIgnoreCase("false") ? false : true;
            }

            
            // =====================Application Versions
            if (getProperty("bicseq_version") == null) {
                Logger.getLogger(CNVWorkflow.class.getName()).log(Level.SEVERE, "bicseq_version is not set, we need it to call BICseq correctly");
                return (null);
            } else {
                this.bicseqVersion = getProperty("bicseq_version");
            }

            if (getProperty("freec_version") == null) {
                Logger.getLogger(CNVWorkflow.class.getName()).log(Level.SEVERE, "freec_version is not set, we need it to call FREEC correctly");
                return (null);
            } else {
                this.freecVersion = getProperty("freec_version");
            }
            
            if (getProperty("samtools_version") == null) {
                Logger.getLogger(CNVWorkflow.class.getName()).log(Level.SEVERE, "samtools_version is not set, we need it to call samtools correctly");
                return (null);
            } else {
                this.samtoolsVersion = getProperty("samtools_version");
            }
            
            nextProp = getProperty("freec_var_coefficient");
            if (nextProp == null) {
                this.freecVarCoeff = FREEC_CV_DEFAULT;
                Logger.getLogger(CNVWorkflow.class.getName()).log(Level.WARNING, "freec_var_coefficient is not set, will use default value " + FREEC_CV_DEFAULT);
            } else {
                this.freecVarCoeff = nextProp.isEmpty() ? FREEC_CV_DEFAULT : getProperty("freec_var_coefficient");;
            }
           
            //=============A special flag that determines if we need to sort/index
            nextProp = getProperty("do_sort");
            if (nextProp == null || nextProp.isEmpty()) {
                Logger.getLogger(CNVWorkflow.class.getName()).log(Level.WARNING, "do_sort is not set, will deduce it from the names of input files");
                this.doSort = !this.normal[0].contains("sorted");
            } else {
                this.doSort = nextProp.isEmpty() || nextProp.equalsIgnoreCase("false") ? false : true;
            }

            // Register input files and set up local files
            this.localInputNormalFiles = new String[this.normal.length];
            this.localInputTumorFiles = new String[this.tumor.length];
            this.normalBases = new String[this.normal.length];
            this.tumorBases = new String[this.tumor.length];
            this.bicseqInterval = Integer.valueOf(getOptionalProperty("biqseq_interval", BICSEQ_I_DEFAULT));
            this.bicseqSpread   = Integer.valueOf(getOptionalProperty("biqseq_spread",   BICSEQ_S_DEFAULT));

            String[] types = {"normal", "tumor"};
            for (String type : types) {
                int listSize = type.equals("normal") ? this.normal.length : this.tumor.length;
                for (int fileIndex = 0; fileIndex < listSize; fileIndex++) {
                    String bamBasename;
                    /*
                     * vcftools analyzes filename, and if .gz is present, it will treat the file az bgzipped
                     * So, in our case we need to make sure that we don't have .gz in our files' names
                     */
                    if (type.equals("normal")) {
                        bamBasename = this.normal[fileIndex].substring(this.normal[fileIndex].lastIndexOf("/") + 1, this.normal[fileIndex].lastIndexOf(".bam"));
                        this.normalBases[fileIndex] = bamBasename.replaceAll(".gz.", ".");
                        this.normalBases[fileIndex] = this.normalBases[fileIndex].replaceAll(".fastq.annotated", "");
                    } else {
                        bamBasename = this.tumor[fileIndex].substring(this.tumor[fileIndex].lastIndexOf("/") + 1, this.tumor[fileIndex].lastIndexOf(".bam"));
                        this.tumorBases[fileIndex] = bamBasename.replaceAll(".gz.", ".");
                        this.tumorBases[fileIndex] = this.tumorBases[fileIndex].replaceAll(".fastq.annotated", "");
                    }

                    Log.stdout("CREATING FILE: input_bam_" + fileIndex + "_" + type);

                    SqwFile file = this.createFile("input_bam_" + fileIndex + "_" + type);
                    file.setType("application/bam");
                    file.setIsInput(true);
                    file.setForceCopy(false);

                    if (type.equals("normal")) {
                        file.setSourcePath(this.normal[fileIndex]);
                        this.localInputNormalFiles[fileIndex] = !this.doSort ? file.getProvisionedPath() : this.dataDir + bamBasename + ".sorted.bam";
                    } else {
                        file.setSourcePath(this.tumor[fileIndex]);
                        this.localInputTumorFiles[fileIndex] = !this.doSort ? file.getProvisionedPath() : this.dataDir + bamBasename + ".sorted.bam";
                    }
                }
            }// finished registering input bam files

            return this.getFiles();

        } catch (Exception ex) {
            Logger.getLogger(CNVWorkflow.class.getName()).log(Level.SEVERE, null, ex);
            return (null);
        }
    }

    @Override
    public void setupDirectory() {
        try {
            this.dataDir = getProperty("data_dir");
            this.tempDir = "tempfiles/";
            if (this.dataDir == null || this.dataDir.isEmpty()) {
                this.dataDir = "data/";
            }
            if (!this.dataDir.endsWith("/")) {
                this.dataDir = this.dataDir.concat("/");
            }
            this.addDirectory(this.tempDir);
            this.addDirectory(this.dataDir);

        } catch (Exception e) {
            Logger.getLogger(CNVWorkflow.class.getName()).log(Level.WARNING, null, e);
        }
    }

    @Override
    public void buildWorkflow() {

        try {
            
            
          // TODO - crosscheck configured and routed here, each job should accept a pair of normal/tumor bam or just one bam
          // for reference-free analysis
          for (int n = 0; n < this.normal.length; n++) {
               for (int t = 0; t < this.tumor.length; t++) {
               //TODO need to think how to configure this for crosscheck properly if we don't have reference
               launchBicSeq(this.localInputNormalFiles[n],
                            this.localInputTumorFiles[t], n);
                   
               }
          }
          
          // if this.templateType == "WG"
          //==============Launch Supported Analyses for Whole Genome Data
               //launchBicSeq(this.localInputNormalFiles[n],
               //             this.localInputTumorFiles[t], n);
          // launchHMMcopy();
          // launchSupportedAnalyses("WG");
          
          
          // else if this.templateType == "EX"
          //==============Launch Supported Analyses for Whole Exome Data
          // launch VarSeq2();
          // launch FREEC();
          // launchSupportedAnalyses("EX);


        } catch (Exception e) {
            Logger.getLogger(getClass().getName()).log(Level.SEVERE, null, e);
        }
    }

    /**
     * BICseq configuring/launching
     */
    private void launchBicSeq(String inputNormal, String inputTumor, int id) {

        // Job convertJob and create configFile
        Job convertJob = this.getWorkflow().createBashJob("bicseq_prepare");
            
        String configFile = "bicseq_config." + id + ".conf";
        convertJob.setCommand(getWorkflowBaseDir() + "/dependencies/configureBICseq.pl"
                            + " --input-normal " + inputNormal
                            + " --input-tumor "  + inputTumor
                            + " --outdir " + this.dataDir
                            + " --config-file " + configFile
                            + " --samtools " + getWorkflowBaseDir() + "/bin/BICseq-" + this.bicseqVersion
                            + "/PERL_pipeline/BICseq_" + this.bicseqVersion + "/SAMgetUnique/samtools-0.1.7a_getUnique-0.1.1/samtools");
        convertJob.setMaxMemory("4000");
        Log.stdout("Created BICseq convert Job");
        
        
        // Launch BICSeq, provision results
        // PERL_pipeline/BICseq_1.1.2/BIC-seq/BIC-seq.pl --I 150,20 /u/pruzanov/Data/CNVtools/BICseq/test1.config /scratch2/users/pruzanov/Data/CNVTOOLS/BIC-seq.hn.test1 \"InitialTest\"
        String resultDir = "BICseq_out_" + id;
        Job launchJob = this.getWorkflow().createBashJob("bicseq_launch");
        String resultID = this.makeBasename(inputNormal, ".bam") + ".vs." 
                        + this.makeBasename(inputTumor,  ".bam");
        
        launchJob.setCommand(getWorkflowBaseDir() + "/dependencies/launchBICseq.pl"
                           + " --config-file " + configFile
                           + " --outdir " + this.dataDir + resultDir
                           + " --bigseq-interval " + this.bicseqInterval
                           + " --bicseq-spread "   + this.bicseqSpread
                           + " --result-id " + resultID
                           + " --biqseq " + getWorkflowBaseDir() + "/bin/BICseq-" + this.bicseqVersion
                           + "/PERL_pipeline/BICseq_" + this.bicseqVersion + "/BIC-seq/BIC-seq.pl");
        launchJob.setMaxMemory("6000");
        launchJob.addParent(convertJob);
        Log.stdout("Created BICseq launch Job");
    }
    
    /**
     * HMMcopy configuring/launching
     */
    private void launchHMMcopy() {
        // take inputs, convert into .seq
        
        // compose config file
        
        // Launch HMMcopy scripts, provision results
    }
    
    
    /**
     * Varscan configuring/launching
     */
    private void launchVarscan() {
        
    }
    
    /**
     * FREEC configuring/launching
     */
    private void launchFREEC(String inputNormal, String inputTumor, int id) {
        
        // Job convertJob and create configFile
        Job freecJob = this.getWorkflow().createBashJob("freec_launch");
            
        String configFile = "freec_config." + id + ".conf";
        freecJob.setCommand(getWorkflowBaseDir() + "/dependencies/launchFREEC.pl"
                            + " --input-normal " + inputNormal
                            + " --input-tumor "  + inputTumor
                            + " --lenfile " + this.chrLengthFile
                            + " --id " + id
                            + " --freec " + getWorkflowBaseDir() + "/bin/FREEC-" + this.freecVersion + "/freec"
                            + " --data-type " + this.templateType
                            + " --outdir " + this.dataDir
                            + " --samtools " + getWorkflowBaseDir() + "/bin/samtools-" + this.samtoolsVersion + "/samtools");
                            
        if (!this.freecVarCoeff.isEmpty()) {
         freecJob.getCommand().addArgument(" var-coefficient " + this.freecVarCoeff);
        }              
                
        freecJob.setMaxMemory("6000");
        
        Log.stdout("Created FREEC launch Job");
    }
    
    /**
     * Utility function
     * 
     * @param path
     * @param extension
     * @return 
     */
    private String makeBasename(String path, String extension) {
        return path.substring(path.lastIndexOf("/"), path.lastIndexOf(extension));
    }
}
