package ca.on.oicr.pde.workflows;

import ca.on.oicr.pde.utilities.workflows.OicrWorkflow;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import net.sourceforge.seqware.common.util.Log;
import net.sourceforge.seqware.pipeline.workflowV2.model.Command;
import net.sourceforge.seqware.pipeline.workflowV2.model.Job;
import net.sourceforge.seqware.pipeline.workflowV2.model.SqwFile;

public class CNVWorkflow extends OicrWorkflow {

    //Versions of tools
    private String bicseqVersion;
    private String freecVersion;
    private String varscanVersion;
    private String samtoolsVersion;
    private String hmmcopyVersion;
    private String rModule;
    private String rLibDir;
    
    //FREEC
    private String chrLengthFile = "";
    private String templateType;
    private String freecVarCoeff = "";
    private String freecWindow;
            
    //References
    private boolean manualOutput;
    private boolean doSort = true;
    private String  queue;
    private String  refFasta;
    private String  targetFile = "";
    private String[] supportedChromosomes;
    
    //HMMcopy readCounter
    private String readCounterChromosomes;
    private String readCounterWindow;

    //HMMcopy
    private String refGCfile;
    private String refMAPfile;
    
    /**   Varscan Filtering parameters:
     *    Defaults are set by the software, not this workflow
	--min-coverage	Minimum read depth at a position to make a call [8]
	--amp-threshold	Lower bound for log ratio to call amplification [0.25]
	--del-threshold	Upper bound for log ratio to call deletion (provide as positive number) [0.25]
	--min-region-size	Minimum size (in bases) for a region to be counted [10]
	--recenter-up	Recenter data around an adjusted baseline > 0 [0]
	--recenter-down	Recenter data around an adjusted baseline < 0 [0]
        
     */
    
    //Varscan
    private String varscanMinCoverage;
    private String varscanDelCoverage;
    private String varscanMinRegion;
    private String varscanRecenterUp;
    private String varscanRecenterDown;
    private String varscanPvalueThreshold;
    private String varscanJavaXmx;
    
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

    private boolean skipFlag;
    private static final String BICSEQ_I_DEFAULT         = "150";
    private static final String BICSEQ_S_DEFAULT         = "20";
    private static final String FREEC_CV_DEFAULT         = "0.05";
    private static final String FREEC_WINDOW_DEFAULT_EX  = "500";
    private static final String FREEC_WINDOW_DEFAULT_WG  = "50000";
    private static final boolean DEFAULT_SKIP_IF_MISSING = true;  // Conditional provisioning
    private static final String  FREEC_PREFIX  = "freec_";
    private static final String HMMCOPY_PREFIX = "hmmcopy_";
    private static final String BICSEQ_PREFIX  = "bicseq_";
    private static final String RLIBDIR_BASE   = "CNV.R_modules-";
    private static final String VARSCAN_PREFIX = "varscan_";
    private final static String WG           = "WG";
    private final static String PVALUE         = "0.05";
    private static final String VARSCAN_JAVA_MEM = "4";
    
    
    /**
     * 
     * 
     *
     */
    
    @Override
    public Map<String, SqwFile> setupFiles() {

        try {

            this.normal        = getProperty("input_files_normal").split(",");
            this.tumor         = getProperty("input_files_tumor").split(",");
            this.refFasta      = getProperty("reference_fasta");
            this.refGCfile     = getProperty("reference_gc");
            this.refMAPfile    = getProperty("reference_map");
            this.chrLengthFile = getProperty("reference_len_file");
            this.templateType  = getProperty("template_type");
            this.supportedChromosomes = getProperty("supported_chromosomes").split(",");

            this.readCounterWindow = getProperty("readcounter_window");
            this.readCounterChromosomes = getProperty("supported_chromosomes");

            String fileSkipFlag = this.getOptionalProperty("skip_missing_files", Boolean.toString(DEFAULT_SKIP_IF_MISSING));
            this.varscanPvalueThreshold = this.getOptionalProperty("varscan_pvalue", PVALUE);
            this.varscanJavaXmx = this.getOptionalProperty("varscan_java_xmx", VARSCAN_JAVA_MEM);
            this.skipFlag = Boolean.valueOf(fileSkipFlag);
            
            //allow all template types other than WG (given that their have valid interval file associated)
            if (!this.templateType.equals(WG)) { 
                this.targetFile  = getProperty("target_file");
                this.freecWindow = getOptionalProperty("freec_window", FREEC_WINDOW_DEFAULT_EX);
            } else {
                this.freecWindow     = getOptionalProperty("freec_window", FREEC_WINDOW_DEFAULT_WG);
            }
            
            this.queue = getProperty("queue");
            
            if (this.tumor.length != this.normal.length) {
                    if (this.normal.length != 1) {
                        Logger.getLogger(CNVWorkflow.class.getName()).log(Level.SEVERE, "numbers for normal and tumor bam files are not the same, crosscheck isn't forced - "
                                + "check your .ini file");
                        return (null);
                    }
            }

            //=====================Application Versions           
            this.bicseqVersion   = getProperty("bicseq_version");
            this.freecVersion    = getProperty("freec_version");
            this.samtoolsVersion = getProperty("samtools_version");
            this.varscanVersion  = getProperty("varscan_version");
            this.samtoolsVersion = getProperty("samtools_version");
            this.hmmcopyVersion  = getProperty("hmmcopy_version");
            this.rModule         = getProperty("R_module");
            this.rLibDir         = RLIBDIR_BASE + getProperty("Rlibs_version");
            this.freecVarCoeff   = getOptionalProperty("freec_var_coefficient", FREEC_CV_DEFAULT);
            this.varscanMinCoverage = getOptionalProperty("varscan_min_coverage", null);
            this.varscanDelCoverage = getOptionalProperty("varscan_del_coverage", null);
            this.varscanMinRegion   = getOptionalProperty("varscan_min_region", null);
            this.varscanRecenterUp  = getOptionalProperty("varscan_recenter_up", null);
            this.varscanRecenterUp  = getOptionalProperty("varscan_recenter_down", null);
            
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
               /**
                * TODO need to think how to configure this for crosscheck properly 
                * if we don't have reference (TBD next iteration)
                */
               if (this.templateType.equals(WG)) {
                 // LAUNCH BICseq
                 launchBicSeq(this.localInputNormalFiles[n],
                              this.localInputTumorFiles[t], n + 1, sortJobs);
                 // LAUNCH HMMcopy
                 launchHMMcopy(this.localInputNormalFiles[n],
                               this.localInputTumorFiles[t], n + 1, sortJobs);
                 
                 // LAUNCH FREEC
                 launchFREEC(this.localInputNormalFiles[n],
                             this.localInputTumorFiles[t], n + 1, null);                
               } else if (!this.targetFile.isEmpty()) {
                 // LAUNCH FREEC
                 launchFREEC(this.localInputNormalFiles[n],
                             this.localInputTumorFiles[t], n + 1, sortJobs);
                 // LAUNCH Varscan
                 launchVarscan(this.localInputNormalFiles[n],
                               this.localInputTumorFiles[t], n + 1, null);
                 
               } else {
                   throw new RuntimeException("Unsupported template type, workflow will terminate!");
               }
                   
              }
          }
          
          // Summary job may be added in a next release

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
        String resultDir = "BICseq_out_" + id + "/";
        Job launchJob = this.getWorkflow().createBashJob("bicseq_launch");
        String resultID = BICSEQ_PREFIX + this.makeBasename(inputNormal, ".bam") + ".vs." 
                                        + this.makeBasename(inputTumor,  ".bam");
        
        launchJob.setCommand("module load " + this.rModule + ";"
                           + getWorkflowBaseDir() + "/dependencies/launchBICseq.pl"
                           + " --config-file " + this.dataDir + configFile
                           + " --outdir " + this.dataDir + resultDir
                           + " --bicseq-interval " + this.bicseqInterval
                           + " --bicseq-spread "   + this.bicseqSpread
                           + " --result-id " + resultID
                           + " --bicseq " + getWorkflowBaseDir() + "/bin/BICseq-" + this.bicseqVersion
                           + "/PERL_pipeline/BICseq_" + this.bicseqVersion + "/BIC-seq/BIC-seq.pl");
        launchJob.setMaxMemory("6000");
        launchJob.addParent(convertJob);
        Log.stdout("Created BICseq launch Job");
        
        // Provision files normal.vs.tumor.bicseg, normal.vs.tumor.png, normal.vs.tumor.wig
        SqwFile bicseqSegFile = createOutputFile(this.dataDir + resultDir + resultID + ".bicseg", "text/plain", this.manualOutput);
        bicseqSegFile.setSkipIfMissing(skipFlag);
        bicseqSegFile.getAnnotations().put("variation_calling_algorithm", "BICseq " + this.bicseqVersion);
        launchJob.addFile(bicseqSegFile);
        
        SqwFile bicseqPngFile = createOutputFile(this.dataDir + resultDir + resultID + ".png",    "image/png",  this.manualOutput);
        bicseqPngFile.setSkipIfMissing(skipFlag);
        bicseqPngFile.getAnnotations().put("variation_calling_algorithm", "BICseq " + this.bicseqVersion);
        launchJob.addFile(bicseqPngFile);
        
        SqwFile bicseqWigFile = createOutputFile(this.dataDir + resultDir + resultID + ".wig",    "text/plain", this.manualOutput);
        bicseqWigFile.setSkipIfMissing(skipFlag);
        bicseqWigFile.getAnnotations().put("variation_calling_algorithm", "BICseq " + this.bicseqVersion);
        launchJob.addFile(bicseqWigFile);
    }
    
    /**
     * HMMcopy configuring/launching
     */
    private void launchHMMcopy(String inputNormal, String inputTumor, int id, List<Job> parents) {
        
        String[] allInputs = {inputNormal, inputTumor};
        String outputDir = this.dataDir + "HMMcopy." + id + "/";
        List<Job> setupJobs = new ArrayList<Job>();
       
        // Inputs converted into .wig and then - .bw format
        for (String inFile : allInputs) {
        // =============== indexing =========================================
        Job indexJob = this.getWorkflow().createBashJob("hmmcopy_index");
        Command indexCommand = indexJob.getCommand();
        indexCommand.addArgument("mkdir -p " + outputDir + ";");
        indexCommand.addArgument(getWorkflowBaseDir() + "/bin/HMMcopy-" + this.hmmcopyVersion + "/bin/readCounter");
        indexCommand.addArgument("--window");
        indexCommand.addArgument(readCounterWindow);
        indexCommand.addArgument("--chromosome");
        indexCommand.addArgument("\"" + readCounterChromosomes + "\"");
        indexCommand.addArgument("--build");
        indexCommand.addArgument(inFile);
        indexJob.setMaxMemory("4000");
            
        if (parents != null) {
          for (Job p : parents) {
             indexJob.addParent(p);
          }
        }
            
        //============ converting to wig format =============================
        Job convertJob = this.getWorkflow().createBashJob("hmmcopy_convert");
        Command convertCommand = convertJob.getCommand();
        convertCommand.addArgument(getWorkflowBaseDir() + "/bin/HMMcopy-" + this.hmmcopyVersion + "/bin/readCounter");
        convertCommand.addArgument("--window");
        convertCommand.addArgument(readCounterWindow);
        convertCommand.addArgument("--chromosome");
        convertCommand.addArgument("\"" + readCounterChromosomes + "\"");
        convertCommand.addArgument(inFile);
        convertCommand.addArgument(">");
        convertCommand.addArgument(this.makeBasename(inFile, ".bam") + "_reads.wig");
        convertJob.setMaxMemory("4000");
        convertJob.addParent(indexJob);
            
        setupJobs.add(convertJob);
        }
       
       // Launch HMMcopy scripts, provision results
       //================== run HMMcopy ====================================
        Job hmmJob = this.getWorkflow().createBashJob("hmmcopy_launch");
        String resultID = HMMCOPY_PREFIX + this.makeBasename(inputNormal, ".bam") + ".vs." 
                                         + this.makeBasename(inputTumor,  ".bam");
        hmmJob.setCommand("module load " + this.rModule + ";"
                        + getWorkflowBaseDir() + "/dependencies/launchHMMcopy.pl "
                        + " --r-libdir "     + getWorkflowBaseDir() + "/bin/" + this.rLibDir
                        + " --normal-wig "   + this.makeBasename(inputNormal, ".bam") + "_reads.wig "
                        + " --tumor-wig "    + this.makeBasename(inputTumor, ".bam") + "_reads.wig "
                        + " --cg-file "      + this.refGCfile
                        + " --map-file "     + this.refMAPfile
                        + " --hmm-script "   + getWorkflowBaseDir() + "/dependencies/run_HMMcopy.r"
                        + " --output-base "  + outputDir + resultID);

        hmmJob.setMaxMemory("6000");
        for (Job p : setupJobs) {
            hmmJob.addParent(p);
        }
        
        Log.stdout("Created HMMcopy launch Job");
        
        // Provision .seg, .tsv, .bias_plot.png, .c_plot.chr*.png, .s_plot.chr*.png
        SqwFile hmmcopySegFile = createOutputFile(outputDir + resultID + ".seg", "text/plain", this.manualOutput);
        hmmcopySegFile.setSkipIfMissing(skipFlag);
        hmmcopySegFile.getAnnotations().put("variation_calling_algorithm", "HMMcopy " + this.hmmcopyVersion);
        hmmJob.addFile(hmmcopySegFile);
        
        SqwFile hmmcopyTsvFile = createOutputFile(outputDir + resultID + ".tsv", "text/plain", this.manualOutput);
        hmmcopyTsvFile.setSkipIfMissing(skipFlag);
        hmmcopyTsvFile.getAnnotations().put("variation_calling_algorithm", "HMMcopy " + this.hmmcopyVersion);
        hmmJob.addFile(hmmcopyTsvFile);
        
        SqwFile hmmcopyBiasPlotFile = createOutputFile(outputDir + resultID + ".png", "image/png", this.manualOutput);
        hmmcopyBiasPlotFile.setSkipIfMissing(skipFlag);
        hmmcopyBiasPlotFile.getAnnotations().put("variation_calling_algorithm", "HMMcopy " + this.hmmcopyVersion);
        hmmJob.addFile(hmmcopyBiasPlotFile);
        
        for(String chrom : this.supportedChromosomes) {

            SqwFile hmmcopyCPlotFile = createOutputFile(outputDir + resultID + ".c_plot." + chrom + ".png", "image/png", this.manualOutput);
            hmmcopyCPlotFile.setSkipIfMissing(skipFlag);
            hmmcopyCPlotFile.getAnnotations().put("variation_calling_algorithm", "HMMcopy " + this.hmmcopyVersion);
            hmmJob.addFile(hmmcopyCPlotFile);
            
            SqwFile hmmcopySPlotFile = createOutputFile(outputDir + resultID + ".s_plot." + chrom + ".png", "image/png", this.manualOutput);
            hmmcopySPlotFile.setSkipIfMissing(skipFlag);
            hmmcopySPlotFile.getAnnotations().put("variation_calling_algorithm", "HMMcopy " + this.hmmcopyVersion);
            hmmJob.addFile(hmmcopySPlotFile);
   
        }
    }
         

    /**
     * Varscan configuring/launching
     */
    private void launchVarscan(String inputNormal, String inputTumor, int id, List<Job> parents) {
        
        Job varscanJob = this.getWorkflow().createBashJob("launch_varscan");   
        String outputDir = this.dataDir + "Varscan2." + id + "/"; 
        String resultID = VARSCAN_PREFIX + this.makeBasename(inputNormal, ".bam") + ".vs." 
                                         + this.makeBasename(inputTumor,  ".bam");
        varscanJob.setCommand("module load " + this.rModule + ";"
                            + getWorkflowBaseDir() + "/dependencies/launchVarscan2.pl"
                            + " --input-normal " + inputNormal
                            + " --input-tumor "  + inputTumor
                            + " --output-dir "   + outputDir
                            + " --r-libdir "     + getWorkflowBaseDir() + "/bin/" + this.rLibDir
                            + " --ref-fasta "    + refFasta
                            + " --java "         + getWorkflowBaseDir() + "/bin/jre" + getProperty("jre-version") + "/bin/java"
                            + " --varscan "      + getWorkflowBaseDir() + "/bin/VarScan.v" + varscanVersion + ".jar"
                            + " --id "           + resultID
                            + " --xmxmem "       + this.varscanJavaXmx
                            + " --p-value "      + this.varscanPvalueThreshold
                            + " --samtools "     + getWorkflowBaseDir() + "/bin/samtools-" 
                                                 + this.samtoolsVersion + "/samtools");
        if (null != this.varscanMinCoverage) {
            varscanJob.getCommand().addArgument(" --min-coverage " + this.varscanMinCoverage);
        }
        
        if (null != this.varscanDelCoverage) {
            varscanJob.getCommand().addArgument(" --del-coverage " + this.varscanDelCoverage);
        }
        
        if (null != this.varscanMinRegion) {
            varscanJob.getCommand().addArgument(" --min-region-size " + this.varscanMinRegion);
        }
                            
        if (null != this.varscanRecenterUp) {
            varscanJob.getCommand().addArgument(" --recenter-up " + this.varscanRecenterUp);
        }
        
        if (null != this.varscanRecenterDown) {
            varscanJob.getCommand().addArgument(" --recenter-down "+ this.varscanRecenterDown);
        }

        varscanJob.setMaxMemory("12000");
        if (parents != null) {
            for (Job p : parents) {
                varscanJob.addParent(p);
            }
        }
        Log.stdout("Created Varscan launch Job");
        
        // Provision files .copynumber, .copynumber.segmented, .copynumber.filtered
        //                 .copynumber.filtered.s_plot.png, .copynumber.filtered.s_plot.png       
        
        SqwFile varscanCopyFile = createOutputFile(outputDir + resultID + ".copynumber", 
                                                   "text/plain", this.manualOutput);
        varscanCopyFile.setSkipIfMissing(skipFlag);
        varscanCopyFile.getAnnotations().put("variation_calling_algorithm", "Varscan " + this.varscanVersion);
        varscanJob.addFile(varscanCopyFile);
        
        SqwFile varscanCopySegFile = createOutputFile(outputDir + resultID + ".copynumber.segmented",
                                                      "text/plain", this.manualOutput);
        varscanCopySegFile.setSkipIfMissing(skipFlag);
        varscanCopySegFile.getAnnotations().put("variation_calling_algorithm", "Varscan " + this.varscanVersion);
        varscanJob.addFile(varscanCopySegFile);
        
        SqwFile varscanCopyFilteredFile = createOutputFile(outputDir + resultID + ".copynumber.filtered",
                                                           "text/plain", this.manualOutput);
        varscanCopyFilteredFile.setSkipIfMissing(skipFlag);
        varscanCopyFilteredFile.getAnnotations().put("variation_calling_algorithm", "Varscan " + this.varscanVersion);
        varscanJob.addFile(varscanCopyFilteredFile);
     
        SqwFile varscanWPlotFile = createOutputFile(outputDir + resultID + ".copynumber.filtered.w_plot.png",
                                                    "image/png", this.manualOutput);
        varscanWPlotFile.setSkipIfMissing(skipFlag);
        varscanWPlotFile.getAnnotations().put("variation_calling_algorithm", "Varscan " + this.varscanVersion);
        varscanJob.addFile(varscanWPlotFile);
        
        SqwFile varscanSPlotFile = createOutputFile(outputDir + resultID + ".copynumber.filtered.s_plot.png",
                                                    "image/png", this.manualOutput);
        varscanSPlotFile.setSkipIfMissing(skipFlag);
        varscanSPlotFile.getAnnotations().put("variation_calling_algorithm", "Varscan " + this.varscanVersion);
        varscanJob.addFile(varscanSPlotFile);

    }
    
    /**
     * FREEC configuring/launching
     */
    private void launchFREEC(String inputNormal, String inputTumor, int id, List<Job> parents) {
        
        // Job convertJob and create configFile
        Job freecJob = this.getWorkflow().createBashJob("freec_launch");
        String outputDir = this.dataDir + "FREEC." + id + "/";
        String resultID  = FREEC_PREFIX + this.makeBasename(inputTumor,  ".bam") + ".bam";
        freecJob.setCommand("module load " + this.rModule + ";"
                            + getWorkflowBaseDir() + "/dependencies/launchFREEC.pl"
                            + " --r-libdir "     + getWorkflowBaseDir() + "/bin/" + this.rLibDir
                            + " --input-normal " + inputNormal
                            + " --input-tumor "  + inputTumor
                            + " --lenfile "      + this.chrLengthFile
                            + " --id "           + id
                            + " --freec " + getWorkflowBaseDir() + "/bin/FREEC-" + this.freecVersion + "/freec"
                            + " --data-type " + this.templateType
                            + " --outdir "    + outputDir
                            + " --prefix "    + FREEC_PREFIX
                            + " --samtools "  + getWorkflowBaseDir() + "/bin/samtools-" + this.samtoolsVersion + "/samtools");
                            
        if (!this.freecVarCoeff.isEmpty()) {
         freecJob.getCommand().addArgument(" --var-coefficient " + this.freecVarCoeff);
        }
        if (!this.templateType.equals(WG)) {
         freecJob.getCommand().addArgument(" --target-file " + this.targetFile);
         freecJob.getCommand().addArgument(" --window "      + this.freecWindow);
        }
                
        freecJob.setMaxMemory("16000");
        if (parents != null) {
            for (Job p : parents) {
                freecJob.addParent(p);
            }
        }
        Log.stdout("Created FREEC launch Job");
        
        // Provision [tumor bam]_CNVs.p.value.txt, *_ratio.BedGraph, *_ratio_noNA.txt.png, *_sample.cpn, *_control.cpn
        SqwFile freecCNVFile = createOutputFile(outputDir + resultID + "_CNVs.p.value.txt", 
                                                "text/plain", this.manualOutput);
        freecCNVFile.setSkipIfMissing(skipFlag);
        freecCNVFile.getAnnotations().put("variation_calling_algorithm", "FREEC " + this.freecVersion);
        freecJob.addFile(freecCNVFile);
        
        SqwFile freecBedGraphFile = createOutputFile(outputDir + resultID + "_ratio.BedGraph",
                                                     "text/bed", this.manualOutput);
        freecBedGraphFile.setSkipIfMissing(skipFlag);
        freecBedGraphFile.getAnnotations().put("variation_calling_algorithm", "FREEC " + this.freecVersion);
        freecJob.addFile(freecBedGraphFile);
        
        SqwFile freecRatioPlotFile = createOutputFile(outputDir + resultID + "_ratio_noNA.txt.png",
                                                      "image/png", this.manualOutput);
        freecRatioPlotFile.setSkipIfMissing(skipFlag);
        freecRatioPlotFile.getAnnotations().put("variation_calling_algorithm", "FREEC " + this.freecVersion);
        freecJob.addFile(freecRatioPlotFile);
        
        // Raw (copy number profile) files
        SqwFile freecSampleCpnFile = createOutputFile(outputDir + resultID + "_sample.cpn",
                                                      "text/plain", this.manualOutput);
        freecSampleCpnFile.setSkipIfMissing(skipFlag);
        freecSampleCpnFile.getAnnotations().put("variation_calling_algorithm", "FREEC " + this.freecVersion);
        freecJob.addFile(freecSampleCpnFile);
        
        SqwFile freecControlCpnFile = createOutputFile(outputDir + this.makeBasename(inputNormal,  ".bam") + ".bam" 
                                                                           + "_control.cpn", "text/plain", this.manualOutput);
        freecControlCpnFile.setSkipIfMissing(skipFlag);
        freecControlCpnFile.getAnnotations().put("variation_calling_algorithm", "FREEC " + this.freecVersion);
        freecJob.addFile(freecControlCpnFile);
    }
    
    /**
     * Utility function
     * 
     * @param path
     * @param extension
     * @return 
     */
    private String makeBasename(String path, String extension) {
        return path.substring(path.lastIndexOf("/") + 1, path.lastIndexOf(extension));
    }

}
