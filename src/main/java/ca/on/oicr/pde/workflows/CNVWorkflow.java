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
    private String mutectVersion;
    private String strelkaVersion;
    private String samtoolsVersion;
    private String vcftoolsVersion;
    private String bicseqVersion;
    private String picardVersion;
    
    //References
    private String referenceFasta;
    private String bundledJRE;
    private String cosmicVcf;
    private String dbSNPvcf;
    private boolean manualOutput;
    private boolean doSort = true;
    private String queue;

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
    private String alignerSoftware;
    private int bicseqInterval;
    private int bicseqSpread;

    private boolean doCrosscheck = false;
    private static final String BICSEQ_I_DEFAULT = "150";
    private static final String BICSEQ_S_DEFAULT = "20";

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

            nextProp = getProperty("jre-version");
            if (nextProp == null) {
                Logger.getLogger(CNVWorkflow.class.getName()).log(Level.SEVERE, "jre-version is not set, we need it to run picard");
                return (null);
            } else {
                this.bundledJRE = nextProp;
            }

            if (getProperty("reference_fasta") == null) {
                Logger.getLogger(CNVWorkflow.class.getName()).log(Level.WARNING, "reference_fasta is not set!");
            } else {
                this.referenceFasta = getProperty("reference_fasta");
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
            nextProp = getProperty("mutect_version");
            if (nextProp == null) {
                Logger.getLogger(CNVWorkflow.class.getName()).log(Level.WARNING, "MuTect version not set");
                this.mutectVersion = "";
                return (null);
            } else {
                this.mutectVersion = nextProp;
            }

            nextProp = getProperty("strelka_version");
            if (nextProp == null) {
                Logger.getLogger(CNVWorkflow.class.getName()).log(Level.WARNING, "Strelka version not set");
                this.strelkaVersion = "";
                return (null);
            } else {
                this.strelkaVersion = nextProp;
            }

            nextProp = getProperty("samtools_version");
            if (nextProp == null) {
                Logger.getLogger(CNVWorkflow.class.getName()).log(Level.WARNING, "samtools version not set");
                this.samtoolsVersion = "";
                return (null);
            } else {
                this.samtoolsVersion = nextProp;
            }

            nextProp = getProperty("vcftools_version");
            if (nextProp == null) {
                Logger.getLogger(CNVWorkflow.class.getName()).log(Level.SEVERE, "vcftools_version is not set, we need it to call vcftools correctly");
                return (null);
            } else {
                this.vcftoolsVersion = getProperty("vcftools_version");
            }

            if (getProperty("bicseq_version") == null) {
                Logger.getLogger(CNVWorkflow.class.getName()).log(Level.SEVERE, "bicseq_version is not set, we need it to call BICseq correctly");
                return (null);
            } else {
                this.bicseqVersion = getProperty("bicseq_version");
            }

            if (getProperty("picard_version") == null) {
                Logger.getLogger(CNVWorkflow.class.getName()).log(Level.SEVERE, "picard_version is not set, we need it to call tabix correctly");
                return (null);
            } else {
                this.picardVersion = getProperty("picard_version");
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
               }
          }
          
          // if this.templateType == "WG"
          //==============Launch Supported Analyses for Whole Genome Data
          // launchBicSeq();
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
    private void launchBicSeq(SqwFile inputNormal, SqwFile inputTumor, int id) {

        // Job convertJob and create configFile
        Job convertJob = this.getWorkflow().createBashJob("bicseq_prepare");
            
        String configFile = "bicseq_config_" + id;// + TODO
        convertJob.setCommand(getWorkflowBaseDir() + "/dependencies/configureBICseq.pl"
                            + " --input-normal " + inputNormal.getProvisionedPath()
                            + " --input-tumor "  + inputTumor.getProvisionedPath()
                            + " --outdir " + this.dataDir
                            + " --config-file " + configFile
                            + " --samtools " + getWorkflowBaseDir() + "/bin/BICseq-" + this.bicseqVersion
                            + "/PERL_pipeline/BICseq_" + this.bicseqVersion + "/SAMgetUnique/samtools-0.1.7a_getUnique-0.1.1/samtools");
        convertJob.setMaxMemory("4000");
        Log.stdout("Created BICseq convert Job");
        
        
        // Launch BICSeq, provision results
        //PERL_pipeline/BICseq_1.1.2/BIC-seq/BIC-seq.pl --I 150,20 /u/pruzanov/Data/CNVtools/BICseq/test1.config /scratch2/users/pruzanov/Data/CNVTOOLS/BIC-seq.hn.test1 \"InitialTest\"
        String resultDir = "BICseq_out";// + TODO
        Job launchJob = this.getWorkflow().createBashJob("bicseq_launch");
        launchJob.setCommand(getWorkflowBaseDir() + "/dependencies/launchBICseq.pl"
                           + " --input-config " + configFile
                           + " --outdir " + this.dataDir + resultDir
                           + " --bigseq-interval" + this.bicseqInterval
                           + " --bicseq-spread"   + this.bicseqSpread
                           + " --biqseq " + getWorkflowBaseDir() + "/bin/BICseq-" + this.bicseqVersion
                           + "/PERL_pipeline/BICseq_" + this.bicseqVersion);
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
    private void launchFREEC() {
        
    }
}
