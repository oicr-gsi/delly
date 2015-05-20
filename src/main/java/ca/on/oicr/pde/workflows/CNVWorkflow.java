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

public class CNVWorkflow extends OicrWorkflow {

    //Versions of tools
    private String bicseqVersion;
    private String freecVersion;
    private String varscanVersion;
    private String samtoolsVersion;
    private String hmmcopyVersion;
    private String rVersion;
    
    //FREEC
    private String chrLengthFile = "";
    private String templateType;
    private String freecVarCoeff = "";
            
    //References
    private boolean manualOutput;
    private boolean doSort = true;
    private String  queue;
    private String  refFasta;
    private String  chromSizeFile;

    //Data
    private String[] normal;
    private String[] tumor;
    private String[] localInputNormalFiles;
    private String[] localInputTumorFiles;
    private String[] normalBases;
    private String[] tumorBases;

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

            this.normal   = getProperty("input_files_normal").split(",");
            this.tumor    = getProperty("input_files_tumor").split(",");
            this.refFasta = getProperty("reference_fasta");
            this.chromSizeFile = getProperty("chromsize_file");
            
            this.queue = getProperty("queue");
            
            if (this.tumor.length != this.normal.length) {
                    if (this.normal.length == 1) {
                        this.doCrosscheck = true;
                    } else {
                        Logger.getLogger(CNVWorkflow.class.getName()).log(Level.SEVERE, "numbers for normal and tumor bam files are not the same, crosscheck isn't forced - "
                                + "check your .ini file");
                        return (null);
                    }
            }

            // =====================Application Versions           
            this.bicseqVersion = getProperty("bicseq_version");
            this.freecVersion = getProperty("freec_version");
            this.samtoolsVersion = getProperty("samtools_version");
            this.varscanVersion  = getProperty("varscan_version");
            this.samtoolsVersion = getProperty("samtools_version");
            this.hmmcopyVersion = getProperty("hmmcopy_version");
            this.rVersion = getProperty("R_version");
            this.freecVarCoeff = getOptionalProperty("freec_var_coefficient", FREEC_CV_DEFAULT);
           
            //=============A special flag that determines if we need to sort/index
            String sortFlag = getProperty("do_sort");
            if (sortFlag == null || sortFlag.isEmpty()) {
                Logger.getLogger(CNVWorkflow.class.getName()).log(Level.WARNING, "do_sort is not set, will deduce it from the names of input files");
                this.doSort = !this.normal[0].contains("sorted");
            } else {
                this.doSort = sortFlag.isEmpty() || sortFlag.equalsIgnoreCase("false") ? false : true;
            }
            
            String manualFlag = getProperty("manual_output");
            if (manualFlag == null) {
                Logger.getLogger(CNVWorkflow.class.getName()).log(Level.WARNING, "manual_output is not set, will put the file into automatically generated dir");
            } else {
                this.manualOutput = manualFlag.isEmpty() || manualFlag.equalsIgnoreCase("false") ? false : true;
            }

            // Register input files and set up local files
            this.localInputNormalFiles = new String[this.normal.length];
            this.localInputTumorFiles  = new String[this.tumor.length];
            this.normalBases           = new String[this.normal.length];
            this.tumorBases            = new String[this.tumor.length];
            this.bicseqInterval        = Integer.valueOf(getOptionalProperty("biqseq_interval", BICSEQ_I_DEFAULT));
            this.bicseqSpread          = Integer.valueOf(getOptionalProperty("biqseq_spread",   BICSEQ_S_DEFAULT));

            String[] types = {"normal", "tumor"};
            for (String type : types) {
                int listSize = type.equals("normal") ? this.normal.length : this.tumor.length;
                for (int fileIndex = 0; fileIndex < listSize; fileIndex++) {
                    String bamBasename;
                    /* vcftools analyzes filename, and if .gz is present, it will treat the file az bgzipped
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
            } // finished registering input bam files

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
          
          List<Job> sortJobs = new ArrayList<Job>();

          if (this.doSort) {
          String[] types = {"normal", "tumor"};
            for (String type : types) {
                int listSize = type.equals("normal") ? this.normal.length : this.tumor.length;
                for (int fileIndex = 0; fileIndex < listSize; fileIndex++) {
                    String bamBasename = "";
                    String filePath = "";
                if (type.equals("normal")) {
                  bamBasename = this.makeBasename(this.localInputNormalFiles[fileIndex], ".bam");
                  filePath = this.normal[fileIndex];
                } else {
                  bamBasename = this.makeBasename(this.localInputTumorFiles[fileIndex], ".bam");
                  filePath = this.tumor[fileIndex];
                }
                  
                Job jobSamSort = this.getWorkflow().createBashJob("index_sort");
                jobSamSort.setCommand(getWorkflowBaseDir() + "/bin/samtools-" + this.samtoolsVersion + "/samtools sort "
                                    + filePath + " "
                                    + this.dataDir + bamBasename); // localIndexed file
                jobSamSort.setMaxMemory("9000");
                 if (!this.queue.isEmpty()) {
                     jobSamSort.setQueue(this.queue);
                 }
                sortJobs.add(jobSamSort);
                }
           }
          }

          for (int n = 0; n < this.normal.length; n++) {
               for (int t = 0; t < this.tumor.length; t++) {
               //TODO need to think how to configure this for crosscheck properly if we don't have reference
               if (this.templateType.equals("WG")) {
                 // LAUNCH BICseq
                 launchBicSeq(this.localInputNormalFiles[n],
                              this.localInputTumorFiles[t], n, sortJobs);
                 // LAUNCH HMMcopy
                 // launchHMMcopy();
                 
                 // LAUNCH FREEC
                 launchFREEC(this.localInputNormalFiles[n],
                             this.localInputTumorFiles[t], n, null);
                 // launchSupportedAnalyses("WG");
               } else if (this.templateType.equals("EX")) {
                 // LAUNCH FREEC
                 launchFREEC(this.localInputNormalFiles[n],
                             this.localInputTumorFiles[t], n, sortJobs);
                 // LAUNCH Varscan
                 launchVarscan(this.localInputNormalFiles[n],
                               this.localInputTumorFiles[t], n, null);
                 
                 // launchSupportedAnalyses("EX);
               } else {
                   throw new RuntimeException("Unsupported template type, workflow will terminate!");
               }
                   
              }
          }
          
          // if this.templateType == "WG"
          //==============Launch Supported Analyses for Whole Genome Data
               //launchBicSeq(this.localInputNormalFiles[n],
               //             this.localInputTumorFiles[t], n);
          
          
          
          // else if this.templateType == "EX"
          //==============Launch Supported Analyses for Whole Exome Data
          // launch VarSeq2();
          // launch FREEC();
          


        } catch (Exception e) {
            Logger.getLogger(getClass().getName()).log(Level.SEVERE, null, e);
        }
    }

    /**
     * BICseq configuring/launching
     */
    private void launchBicSeq(String inputNormal, String inputTumor, int id, List<Job> parents) {

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
        if (parents != null) {
            for (Job p : parents) {
                convertJob.addParent(p);
            }
        }
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
    private void launchHMMcopy(String inputNormal, String inputTumor, int id, List<Job> parents) {
        
        String[] allInputs = {inputNormal, inputTumor};
        String outputDir = this.dataDir + "HMMcopy." + id + "/";
        List<Job> setupJobs = new ArrayList<Job>();
        String normalMpb;
        String tumorMpb;
        
        // Inputs converted into .wig and then - .bw format
        for (String inFile : allInputs) {
            // =============== indexing =========================================
            Job indexJob = this.getWorkflow().createBashJob("hmmcopy_index");
            indexJob.setCommand(getWorkflowBaseDir() + "/bin/HMMcopy-" + this.hmmcopyVersion + "/bin/readCounter -b " + inFile);
            indexJob.setMaxMemory("4000");
            
            if (parents != null) {
             for (Job p : parents) {
                indexJob.addParent(p);
             }
            }
            
            //============ converting to wig format =============================
            Job convertJob = this.getWorkflow().createBashJob("hmmcopy_convert");
            convertJob.setCommand(getWorkflowBaseDir() + "/bin/HMMcopy-" + this.hmmcopyVersion + "/bin/readCounter " + inFile
                                + " > " + outputDir + this.makeBasename(inFile, ".bam") + "_reads.wig");
            convertJob.setMaxMemory("4000");
            convertJob.addParent(indexJob);
            
            //================== convert to bigwig format =======================
            Job convertBwJob = this.getWorkflow().createBashJob("hmmcopy_bw_convert");
            convertBwJob.setCommand(getWorkflowBaseDir() + "/bin/HMMcopy-" + this.hmmcopyVersion + "/util/bigwig/wigToBigWig"
                                  + " -clip " + outputDir + this.makeBasename(inFile, ".bam") + "_reads.wig "
                                  + this.chrLengthFile +  " " + outputDir + this.makeBasename(inFile, ".bam") + ".bw");

            convertBwJob.setMaxMemory("4000");
            convertBwJob.addParent(convertJob);
            
            //===================calculate mappability ==========================
            Job mappableJob = this.getWorkflow().createBashJob("hmmcopy_mappability");
            mappableJob.setCommand(getWorkflowBaseDir() + "/bin/HMMcopy-" + this.hmmcopyVersion + "/bin/mapCounter " 
                                 + outputDir + this.makeBasename(inFile, ".bam") + ".bw > "
                                 + outputDir + this.makeBasename(inFile, ".bam") + ".mappable.wig");
            
            if (inFile.equals(inputNormal)) {
              normalMpb = outputDir + this.makeBasename(inFile, ".bam") + ".mappable.wig";
            } else {
              tumorMpb  = outputDir + this.makeBasename(inFile, ".bam") + ".mappable.wig";
            }
            
            mappableJob.setMaxMemory("4000");
            mappableJob.addParent(convertBwJob);
            
            setupJobs.add(mappableJob);
        }
  
        // cg content calculation for bins:

        // bin/gcCounter <FASTA reference> > gc.wig

        //  average mappability for bins

        // bin/mapCounter <BigWig file> > map.wig
        // Launch HMMcopy scripts, provision results
    }
         
       
    
    
    
    /**
     * Varscan configuring/launching
     */
    private void launchVarscan(String inputNormal, String inputTumor, int id, List<Job> parents) {
        
        Job varscanJob = this.getWorkflow().createBashJob("bicseq_prepare");   
        String outputDir = this.dataDir + "Varscan2." + id + "/"; 
        
        varscanJob.setCommand(getWorkflowBaseDir() + "/dependencies/launchVarscan2.pl"
                            + " --input-normal " + inputNormal
                            + " --input-tumor "  + inputTumor
                            + " --output-dir "   + outputDir
                            + " --rlibs-dir "    + getWorkflowBaseDir() + "/bin"
                            + " --ref-fasta "    + refFasta
                            + " --java "         + getWorkflowBaseDir() + "/bin/jre" + getProperty("jre-version") + "/bin/java"
                            + " --varscan "      + getWorkflowBaseDir() + "/bin/VarScan.v" + varscanVersion
                            + " --id "           + id
                            + " --samtools "     + getWorkflowBaseDir() + "/bin/samtools-" 
                                                 + this.samtoolsVersion + "/samtools");
        varscanJob.setMaxMemory("6000");
        if (parents != null) {
            for (Job p : parents) {
                varscanJob.addParent(p);
            }
        }
        Log.stdout("Created Varscan launch Job");
    }
    
    /**
     * FREEC configuring/launching
     */
    private void launchFREEC(String inputNormal, String inputTumor, int id, List<Job> parents) {
        
        // Job convertJob and create configFile
        Job freecJob = this.getWorkflow().createBashJob("freec_launch");
        String outputDir = this.dataDir + "FREEC." + id + "/";      
        freecJob.setCommand(getWorkflowBaseDir() + "/dependencies/launchFREEC.pl"
                            + " --input-normal " + inputNormal
                            + " --input-tumor "  + inputTumor
                            + " --lenfile " + this.chrLengthFile
                            + " --id " + id
                            + " --freec " + getWorkflowBaseDir() + "/bin/FREEC-" + this.freecVersion + "/freec"
                            + " --data-type " + this.templateType
                            + " --outdir " + outputDir
                            + " --samtools " + getWorkflowBaseDir() + "/bin/samtools-" + this.samtoolsVersion + "/samtools");
                            
        if (!this.freecVarCoeff.isEmpty()) {
         freecJob.getCommand().addArgument(" --var-coefficient " + this.freecVarCoeff);
        }              
                
        freecJob.setMaxMemory("6000");
        if (parents != null) {
            for (Job p : parents) {
                freecJob.addParent(p);
            }
        }
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
