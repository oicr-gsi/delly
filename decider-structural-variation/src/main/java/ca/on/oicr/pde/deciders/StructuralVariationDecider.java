/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ca.on.oicr.pde.deciders;

import java.io.File;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.*;
import net.sourceforge.seqware.common.hibernate.FindAllTheFiles.Header;
import net.sourceforge.seqware.common.module.FileMetadata;
import net.sourceforge.seqware.common.module.ReturnValue;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import net.sourceforge.seqware.common.util.maptools.MapTools;

/**
 *
 * @author pruzanov@oicr.on.ca
 */
public class StructuralVariationDecider extends OicrDecider {
    private final Logger logger = LogManager.getLogger(StructuralVariationDecider.class);
    private SimpleDateFormat format = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss.S");
    private Map<String, BeSmall> fileSwaToSmall;
    private String templateTypeFilter = "";
    private String picard_memory = "6000";
    private String delly_memory  = "8000";
    private String sampleName     = "";
    private String refFasta = "";
    private String queue = " ";
    private String excludeList = "";
    private String mappingQuality = " ";
    private String callMode = "";
    private List<String> duplicates;
    private final static String BAM_METATYPE = "application/bam";
    private static final String ALIGNER_TOKEN = "file.aligner";
    private final static String UNMATCHED = "unmatched";
    private final static String SOMATIC  = "somatic";
 
    public StructuralVariationDecider() {
        super();
        fileSwaToSmall  = new HashMap<String, BeSmall>();
        
        parser.acceptsAll(Arrays.asList("ini-file"), "Optional: the location of the INI file.").withRequiredArg();
        
        parser.accepts("picard-memory","Optional. Set the memory allocated to java heap "
                + "when running Picard merge collapse step, the default is 8000").withRequiredArg();
        parser.accepts("delly-memory","Optional. Set the memory allocated to java heap "
                + "when running Delly step, the default is 4000").withRequiredArg();
        parser.accepts("queue", "Optional: Override the default queue setting (production) setting it to something else").withRequiredArg();
        parser.accepts("ref-fasta", "Optional: the path to reference fasta used by Delly "
	        + " Default: /.mounts/labs/PDE/data/reference/hg19_random/fasta/UCSC/hg19_random.fa").withRequiredArg();
        parser.accepts("template-type", "Optional: may be used as a filter to restrict the analysis to a certain template type").withRequiredArg();
        parser.accepts("mode", "Required: Either somatic or unmatched").withRequiredArg();
        parser.accepts("exclude-list", "Optional: the path to file with targets not used by Delly "
	        + " Default: /.mounts/labs/PDE/data/reference/hg19/delly/human.hg19.excl.tsv").withRequiredArg();
        parser.accepts("mapping-quality", "Optional: parameter used by Delly, not set by default ").withRequiredArg();
    }

    @Override
    public ReturnValue init() {
        logger.info("INIT");
	this.setMetaType(Arrays.asList(BAM_METATYPE)); // It seems that this does not have effect 
                        
        //Handle .ini file - we accept only memeory size allocated to different steps
        if (options.has("ini-file")) {
            File file = new File(options.valueOf("ini-file").toString());
            if (file.exists()) {
                String iniFile = file.getAbsolutePath();
                Map<String, String> iniFileMap = new HashMap<String, String>();
                MapTools.ini2Map(iniFile, iniFileMap);
   
		this.delly_memory = iniFileMap.get("delly_memory");
		this.picard_memory = iniFileMap.get("picard_memory");
            } else {
                 logger.error("The given INI file does not exist: " + file.getAbsolutePath());
                System.exit(1);
            }
	}

	//Group by sample if no other grouping selected
        if (this.options.has("group-by")) {
            logger.info("group-by parameter passed, but this decider does not allow overriding the default grouping, which is group-by FILE_SWA");
        }
        
       	if (this.options.has("root-sample-name")) {
            this.sampleName   = options.valueOf("root-sample-name").toString();
	}
        
        if (this.options.has("queue")) {
            this.queue   = options.valueOf("queue").toString();
	}
        
        if (this.options.has("template-type")) {
            this.templateTypeFilter = options.valueOf("template-type").toString();
            logger.info("Setting template type is not necessary, however if set the decider will run the workflow only on this type of data");
	}
        
        if (this.options.has("picard-memory")) {
            this.picard_memory = options.valueOf("picard-memory").toString();
            logger.info("Setting heap size for Picard step to " + this.picard_memory + " Megabytes as requested");
	}
        
        if (this.options.has("delly-memory")) {
            this.picard_memory = options.valueOf("delly-memory").toString();
            logger.info("Setting memory amount for Delly step to " + this.delly_memory + " Megabytes as requested");
	}

        if (this.options.has("ref-fasta")) {
            this.refFasta = options.valueOf("ref-fasta").toString();
            logger.info("Changing ref-fasta file to " + this.refFasta);
        }
        
        if (this.options.has("exclude-list")) {
            this.excludeList = options.valueOf("exclude-list").toString();
            logger.info("Changing exclude-list file to " + this.excludeList);
        }
        
        if (this.options.has("mapping-quality")) {
            this.mappingQuality = options.valueOf("mapping-quality").toString();
            logger.info("Changing mapping-quality " + this.mappingQuality);
        }

        if (this.options.has("mode")) {
             this.callMode = options.valueOf("mode").toString();
              if (!this.callMode.equals(SOMATIC) && !this.callMode.equals(UNMATCHED)) {
                 logger.error("Mode needs to be specified as either unmatched or somatic");
                 return new ReturnValue(ReturnValue.INVALIDPARAMETERS);                 
              }
        } else {
            logger.error("Mode is a required parameter and needs to be specified");
            return new ReturnValue(ReturnValue.INVALIDPARAMETERS);
        }

        //allows anything defined on the command line to override the defaults here.
        ReturnValue val = super.init();
        return val;
    }

    /**
     * Final check
     * @param commaSeparatedFilePaths
     * @param commaSeparatedParentAccessions
     * @return    */
    @Override
    protected ReturnValue doFinalCheck(String commaSeparatedFilePaths, String commaSeparatedParentAccessions) {
        String[] filePaths = commaSeparatedFilePaths.split(",");
        boolean haveNorm = false;
        boolean haveTumr = false;
        int tumorFiles  = 0;
        int normalFiles = 0;
        int groupIds = 0;

        // Check for duplicate file names and exclude them from analysis
        this.duplicates = detectDuplicates(commaSeparatedFilePaths);
        
        if (this.callMode.equals(UNMATCHED) && filePaths.length > 0) {
            return super.doFinalCheck(commaSeparatedFilePaths, commaSeparatedParentAccessions);
        }
        
        for (String p : filePaths) {
             if (null != this.duplicates && this.duplicates.contains(p)) {
                logger.error("File [" + p + "] has a name that cannot be disambiguated in current set, will skip it");
                continue;
            }
             
            for (BeSmall bs : fileSwaToSmall.values()) {
                if (!bs.getPath().equals(p)) {
                    continue;
                }
                String tt = bs.getTissueType();
                //GP-598: make sure we have group_ids for tumor files if there is more than one
                if (!tt.isEmpty() && tt.equals("R")) {
                    haveNorm = true;
                    normalFiles++;
                } else if (!tt.isEmpty()) {
                    haveTumr = true;
                    tumorFiles++;
                    String gi = bs.getGroupID();
                    if (null != gi && !gi.equals("NA")) {
                        groupIds++;
                    }
                }
            }
        }

        if (tumorFiles > 1 && (tumorFiles != groupIds)) {
            logger.error("We have multiple tumor files but not all of them have group id value, WON'T RUN");
            return new ReturnValue(ReturnValue.INVALIDPARAMETERS);
        }
        
        if (normalFiles > 1) {
            logger.error("We have multiple normal files but we can only run the analysis with one normal, WON'T RUN");
            return new ReturnValue(ReturnValue.INVALIDPARAMETERS);
        }

        if (haveNorm && haveTumr) {
            return super.doFinalCheck(commaSeparatedFilePaths, commaSeparatedParentAccessions);
        }

        String absent = haveNorm ? "Tumor" : "Normal";
        logger.error("Data for " + absent + " tissue are not available, WON'T RUN");
        return new ReturnValue(ReturnValue.INVALIDPARAMETERS);
    }

    @Override
    protected boolean checkFileDetails(ReturnValue returnValue, FileMetadata fm) {
        logger.info("CHECK FILE DETAILS:" + fm.getFilePath());
	this.sampleName      = returnValue.getAttribute(Header.SAMPLE_NAME.getTitle());
        String tissueType    = returnValue.getAttribute(Header.SAMPLE_TAG_PREFIX.getTitle() + "geo_tissue_type");
        String tissueOrigin  = returnValue.getAttribute(Header.SAMPLE_TAG_PREFIX.getTitle() + "geo_tissue_origin");
        if (null != tissueType && null != tissueOrigin) {
            String tissueOriTyp = tissueOrigin + "_" + tissueType;
            this.sampleName = this.sampleName.substring(0, this.sampleName.lastIndexOf(tissueOriTyp) + tissueOriTyp.length());
        }
        
        String currentTtype    = returnValue.getAttribute(Header.SAMPLE_TAG_PREFIX.getTitle() + "geo_library_source_template_type");
        // Filter the data of a different template type if filter is specified
        if (!this.templateTypeFilter.isEmpty() && !this.templateTypeFilter.equalsIgnoreCase(currentTtype)) {
            logger.error("Template type " + currentTtype + " does not pass filter which is set to " + this.templateTypeFilter);
            return false;
        }
        return super.checkFileDetails(returnValue, fm);
    }

    @Override
    public Map<String, List<ReturnValue>> separateFiles(List<ReturnValue> vals, String groupBy) {
        // get files from study
        Map<String, ReturnValue> iusDeetsToRV = new HashMap<String, ReturnValue>();
        // Override the supplied group-by value
        for (ReturnValue currentRV : vals) {
            // Ensure that we have files of application/bam metatype, it looks like if we
            // don't BeSmall will be comparing dates for files of different metatypes, 
            // which is not what we desire (maybe done on BeSmall level, but we are doing it here)
            // TODO : if modified BeSmall takes care of MIME-TYPE, remove this check
            boolean metatypeOK = false;
            
            for (int f = 0; f < currentRV.getFiles().size(); f++) {
               try {
                 if (currentRV.getFiles().get(f).getMetaType().equals(BAM_METATYPE))
                     metatypeOK = true;
               } catch (Exception e) {
                 logger.error("Error checking a file");
               }
            }
            if (!metatypeOK)
                continue; // Go to the next value
            
            BeSmall currentSmall = new BeSmall(currentRV);
            fileSwaToSmall.put(currentRV.getAttribute(groupBy), currentSmall);
            
            //make sure you only have the most recent single file for each
            //sequencer run + lane + barcode + meta-type
            String fileDeets = currentSmall.getIusDetails();
            Date currentDate = currentSmall.getDate();
            
            //if there is no entry yet, add it
            if (iusDeetsToRV.get(fileDeets) == null) {
                logger.debug("Adding file " + fileDeets + " -> \n\t" + currentSmall.getPath());
                iusDeetsToRV.put(fileDeets, currentRV);
            } 
            //if there is an entry, compare the current value to the 'old' one in
            //the map. if the current date is newer than the 'old' date, replace
            //it in the map
            else {
                ReturnValue oldRV = iusDeetsToRV.get(fileDeets);
                
                BeSmall oldSmall = fileSwaToSmall.get(oldRV.getAttribute(Header.FILE_SWA.getTitle()));
                Date oldDate = oldSmall.getDate();
                if (currentDate.after(oldDate)) {
                    logger.debug("Adding file " + fileDeets + " -> \n\t" + currentSmall.getDate()
                            + "\n\t instead of file "
                            + "\n\t" + oldSmall.getDate());
                    iusDeetsToRV.put(fileDeets, currentRV);
                } else {
                    logger.debug("Disregarding file " + fileDeets + " -> \n\t" + currentSmall.getDate()
                            + "\n\tas older than duplicate sequencer run/lane/barcode in favour of "
                            + "\n\t" + oldSmall.getDate());
                    logger.debug(currentDate + " is before " + oldDate);
                }
            }
        }

        // only use those files that entered into the iusDeetsToRV
        // since it's a map, only the most recent values
        List<ReturnValue> newValues = new ArrayList<ReturnValue>(iusDeetsToRV.values());       
        Map<String, List<ReturnValue>> map = new HashMap<String, List<ReturnValue>>();

        // group files according to the designated header (e.g. sample SWID)
        for (ReturnValue r : newValues) {
            String currVal = fileSwaToSmall.get(r.getAttribute(Header.FILE_SWA.getTitle())).getGroupByAttribute();
            
            List<ReturnValue> vs = map.get(currVal);
            if (vs == null) {
                vs = new ArrayList<>();
            }
            vs.add(r);
            map.put(currVal, vs);
        }

        return map;
    }

    
   @Override
    public ReturnValue customizeRun(WorkflowRun run) {
        String inputFile          = "";
        StringBuilder inputTumors = new StringBuilder();
                 
        for (FileAttributes atts : run.getFiles()) {
            String p = atts.getPath();
            for (BeSmall bs : fileSwaToSmall.values()) {
                if (!bs.getPath().equals(p)) {
                    continue;
                }

                if (this.callMode.equals(UNMATCHED)) {
                    if (inputTumors.length() != 0) {
                        inputTumors.append(",");
                    }
                    inputTumors.append(p);
                    continue;
                }
                    
                String tt = bs.getTissueType();
                if (!tt.isEmpty() && tt.equals("R")) {
                    inputFile = p;
                } else if (!tt.isEmpty()) {
                    if (inputTumors.length() != 0) {
                        inputTumors.append(",");
                    }
                    inputTumors.append(p);
                }
            }
        }
        
        if (this.callMode.equals(UNMATCHED)) {
            run.addProperty("input_bams", inputTumors.toString());
        } else {
            run.addProperty("input_bams", inputFile);
            run.addProperty("input_tumors", inputTumors.toString());
        }
               
        run.addProperty("picard_memory", this.picard_memory);
        run.addProperty("delly_memory", this.delly_memory);
        run.addProperty("data_dir", "data");
        run.addProperty("queue", this.queue);
        run.addProperty("sample_name", this.sampleName);
        run.addProperty("mapping_quality", this.mappingQuality);
        run.addProperty("mode", this.callMode);
        
        if (this.refFasta != null && !this.refFasta.isEmpty())
            run.addProperty("ref_fasta", this.refFasta);
        
        if(this.excludeList != null && !this.excludeList.isEmpty())
            run.addProperty("exclude_list", this.excludeList);       
                              
        return new ReturnValue();
    }

   
   public static void main(String args[]){
 
        List<String> params = new ArrayList<>();
        params.add("--plugin");
        params.add(StructuralVariationDecider.class.getCanonicalName());
        params.add("--");
        params.addAll(Arrays.asList(args));
        System.out.println("Parameters: " + Arrays.deepToString(params.toArray()));
        net.sourceforge.seqware.pipeline.runner.PluginRunner.main(params.toArray(new String[params.size()]));
         
    }
   
   /*
    * Utility function for detecting duplicate file names
   */
   public static List<String> detectDuplicates(String commaSeparatedFilePaths) {
       
       String [] filePaths = commaSeparatedFilePaths.split(",");
       List<String> list    = new ArrayList<>();
       List<String> checker = new ArrayList<>();
       
       for (String path : filePaths) {
           String baseName = makeBasename(path, ".bam");
           
           if (checker.contains(baseName) && !list.contains(path)) {
               list.add(path);
           } else {
               checker.add(baseName);
           }
       }

       return list.isEmpty() ? null : list;
       
   }
   
   /**
     * Utility function
     * 
     * @param path
     * @param extension
     * @return 
     */
    public static String makeBasename(String path, String extension) {
        return path.substring(path.lastIndexOf("/") + 1, path.lastIndexOf(extension));
    }
   
   private class BeSmall {

        private Date   date = null;
        private String iusDetails = null;
        private String groupID = null;
        private String groupByAttribute = null;
        private String tissueType = null;
        private String path = null;
        
        public BeSmall(ReturnValue rv) {
            try {
                date = format.parse(rv.getAttribute(Header.PROCESSING_DATE.getTitle()));
            } catch (ParseException ex) {
                logger.error("Bad date!", ex);
                ex.printStackTrace();
            }
            FileAttributes fa = new FileAttributes(rv, rv.getFiles().get(0));
            //TODO check if metatype is doing what it is supposed to do (adding MIME type here)
            iusDetails = fa.getLibrarySample() + fa.getSequencerRun() + fa.getLane() + fa.getBarcode() + fa.getMetatype();
            tissueType = fa.getLimsValue(Lims.TISSUE_TYPE);
            //groupByAttribute = fa.getDonor() + ":" + fa.getLimsValue(Lims.TISSUE_ORIGIN) + ":" + fa.getLimsValue(Lims.LIBRARY_TEMPLATE_TYPE);
            groupID    = fa.getLimsValue(Lims.GROUP_ID);
            if (null == groupID || groupID.isEmpty()) {
                groupID = "NA";
            }
            
            
            StringBuilder gba = new StringBuilder();

            if (callMode.equals(UNMATCHED)) {
                gba.append(fa.getLibrarySample());
                
                if (null != fa.getLimsValue(Lims.TISSUE_TYPE))
                 gba.append(":").append(fa.getLimsValue(Lims.TISSUE_TYPE));
            
                if (null != fa.getLimsValue(Lims.TISSUE_PREP))
                 gba.append(":").append(fa.getLimsValue(Lims.TISSUE_PREP));
            
                if (null != fa.getLimsValue(Lims.TISSUE_REGION))
                 gba.append(":").append(fa.getLimsValue(Lims.TISSUE_REGION));
            } else {
                gba.append(fa.getDonor());
            }
            
            if (null != fa.getLimsValue(Lims.LIBRARY_TEMPLATE_TYPE))
             gba.append(":").append(fa.getLimsValue(Lims.LIBRARY_TEMPLATE_TYPE));
            
            String aligner = rv.getAttribute(ALIGNER_TOKEN);
            if (null != aligner && !aligner.isEmpty()) {
                gba.append(":").append(aligner);
            }
            
                                           
            groupByAttribute = gba.toString();
            path = rv.getFiles().get(0).getFilePath() + "";
        }

        public Date getDate() {
            return date;
        }

        public void setDate(Date date) {
            this.date = date;
        }

        public String getGroupByAttribute() {
            return groupByAttribute;
        }

        public void setGroupByAttribute(String groupByAttribute) {
            this.groupByAttribute = groupByAttribute;
        }
              
        public String getIusDetails() {
            return iusDetails;
        }

        public void setIusDetails(String iusDetails) {
            this.iusDetails = iusDetails;
        }

        public String getPath() {
            return path;
        }

        public void setPath(String path) {
            this.path = path;
        }
        
        public String getTissueType() {
            return tissueType;
        }

        /**
         * @return the groupID
         */
        public String getGroupID() {
            return groupID;
        }

        /**
         * @param groupID the groupID to set
         */
        public void setGroupID(String groupID) {
            this.groupID = groupID;
        }
    }
    
}
