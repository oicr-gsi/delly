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
import net.sourceforge.seqware.common.util.Log;
import net.sourceforge.seqware.common.util.maptools.MapTools;


/**
 *
 * @author pruzanov@oicr.on.ca
 */
public class StructuralVariationDecider extends OicrDecider {

    private SimpleDateFormat format = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss.S");
    private Map<String, BeSmall> fileSwaToSmall;
    
    private String templateTypeFilter = "";
    private String picard_memory = "6000";
    private String delly_memory  = "8000";
    private String sampleName     = "";
    private String output_prefix  = "./";
    private String output_dir = "seqware-results";
    private String refFasta = "";
    private String queue = " ";
    private String excludeList = "";
    private String mappingQuality = " ";
    private String manualOutput = "false";
    private final static String BAM_METATYPE = "application/bam";
 
    public StructuralVariationDecider() {
        super();
        fileSwaToSmall  = new HashMap<String, BeSmall>();
        parser.acceptsAll(Arrays.asList("ini-file"), "Optional: the location of the INI file.").withRequiredArg();
        
        parser.accepts("picard-memory","Optional. Set the memory allocated to java heap "
                + "when running Picard merge collapse step, the default is 8000").withRequiredArg();
        parser.accepts("delly-memory","Optional. Set the memory allocated to java heap "
                + "when running Delly step, the default is 4000").withRequiredArg();
        parser.accepts("output-path", "Optional: the path where the files should be copied to "
                + "after analysis. Corresponds to output-prefix in INI file. Default: ./").withRequiredArg();
        parser.accepts("output-folder", "Optional: the name of the folder to put the output into relative to "
	        + "the output-path. Corresponds to output-dir in INI file. Default: seqware-results").withRequiredArg();
        parser.accepts("queue", "Optional: Override the default queue setting (production) setting it to something else").withRequiredArg();
        parser.accepts("ref-fasta", "Optional: the path to reference fasta used by Delly "
	        + " Default: /.mounts/labs/PDE/data/reference/hg19_random/fasta/UCSC/hg19_random.fa").withRequiredArg();
        parser.accepts("exclude-list", "Optional: the path to file with targets not used by Delly "
	        + " Default: /.mounts/labs/PDE/data/reference/hg19/delly/human.hg19.excl.tsv").withRequiredArg();
        parser.accepts("mapping-quality", "Optional: parameter used by Delly, not set by default ").withRequiredArg();
    }

    @Override
    public ReturnValue init() {
        Log.stdout("INIT");
	this.setMetaType(Arrays.asList(BAM_METATYPE)); // It seems that this does not have effect 
        this.setGroupingStrategy(Header.FILE_SWA); // Need to check if we need it to be so
                
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
                Log.stderr("The given INI file does not exist: " + file.getAbsolutePath());
                System.exit(1);
            }
	}

	//Group by sample if no other grouping selected
        if (this.options.has("group-by")) {
            Log.stderr("group-by parameter passed, but this decider does not allow overriding the default grouping, which is group-by FILE_SWA");
        }
        
       	if (this.options.has("root-sample-name")) {
            this.sampleName   = options.valueOf("root-sample-name").toString();
	}
        
        if (this.options.has("queue")) {
            this.queue   = options.valueOf("queue").toString();
	}
        
        if (this.options.has("template-type")) {
            this.templateTypeFilter = options.valueOf("template-type").toString();
            Log.stdout("Setting template type is not necessary, however if set the decider will run the workflow only on this type of data");
	}
        
        if (this.options.has("picard-memory")) {
            this.picard_memory = options.valueOf("picard-memory").toString();
            Log.stdout("Setting heap size for Picard step to " + this.picard_memory + " Megabytes as requested");
	}
        
        if (this.options.has("delly-memory")) {
            this.picard_memory = options.valueOf("delly-memory").toString();
            Log.stdout("Setting memory amount for Delly step to " + this.delly_memory + " Megabytes as requested");
	}

        if (this.options.has("ref-fasta")) {
            this.refFasta = options.valueOf("ref-fasta").toString();
            Log.stdout("Changing ref-fasta file to " + this.refFasta);
        }
        
        if (this.options.has("exclude-list")) {
            this.excludeList = options.valueOf("exclude-list").toString();
            Log.stdout("Changing exclude-list file to " + this.excludeList);
        }
        
        if (this.options.has("mapping-quality")) {
            this.mappingQuality = options.valueOf("mapping-quality").toString();
            Log.stdout("Changing mapping-quality " + this.mappingQuality);
        }
        
        if (this.options.has("output-path")) {
             this.output_prefix = options.valueOf("output-path").toString();
              if (!this.output_prefix.endsWith("/")) {
                 this.output_prefix += "/";
              }
        }
        
        if (this.options.has("output-folder")) {
            this.output_dir = options.valueOf("output-folder").toString();
	}
        // according to Mei, setting manual_output to true if ouput_prefix is not "./" should not be used, SEQPROD want randomly 
        // named directory in seqware-results dir tree
        if (this.options.has("manual-output") && this.options.valueOf("output-folder").toString().equalsIgnoreCase("true")) {
            this.manualOutput = "true";
	}
        //allows anything defined on the command line to override the defaults here.
        ReturnValue val = super.init();
        return val;
    }


    @Override
    protected boolean checkFileDetails(ReturnValue returnValue, FileMetadata fm) {
        Log.stdout("CHECK FILE DETAILS:" + fm.getFilePath());
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
            Log.stderr("Template type " + currentTtype + " does not pass filter which is set to " + this.templateTypeFilter);
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
                 Log.stderr("Error checking a file");
                 continue;
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
                Log.debug("Adding file " + fileDeets + " -> \n\t" + currentSmall.getPath());
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
                    Log.debug("Adding file " + fileDeets + " -> \n\t" + currentSmall.getDate()
                            + "\n\t instead of file "
                            + "\n\t" + oldSmall.getDate());
                    iusDeetsToRV.put(fileDeets, currentRV);
                } else {
                    Log.debug("Disregarding file " + fileDeets + " -> \n\t" + currentSmall.getDate()
                            + "\n\tas older than duplicate sequencer run/lane/barcode in favour of "
                            + "\n\t" + oldSmall.getDate());
                    Log.debug(currentDate + " is before " + oldDate);
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
                vs = new ArrayList<ReturnValue>();
            }
            vs.add(r);
            map.put(currVal, vs);
        }

        return map;
    }

    
   @Override
    public ReturnValue customizeRun(WorkflowRun run) {
        String inputFiles    = "";
                
        for (FileAttributes atts : run.getFiles()) {
            if (!inputFiles.isEmpty()) {
                inputFiles += ",";
            }
            inputFiles += atts.getPath();
        }

        FileAttributes fa = run.getFiles()[0];
        String tubeId = fa.getLimsValue(Lims.TUBE_ID);
               
        run.addProperty("input_bams", inputFiles);
        run.addProperty("picard_memory", this.picard_memory);
        run.addProperty("delly_memory", this.delly_memory);
        run.addProperty("data_dir", "data");
	run.addProperty("output_prefix",this.output_prefix);
	run.addProperty("output_dir", this.output_dir);
        run.addProperty("queue", this.queue);
        run.addProperty("sample_name", this.sampleName);
        run.addProperty("manual_output", this.manualOutput);
        run.addProperty("mapping_quality", this.mappingQuality);
        
        if (!this.refFasta.isEmpty())
            run.addProperty("ref_fasta", this.refFasta);
        
        if(!this.excludeList.isEmpty())
            run.addProperty("exclude_list", this.excludeList);       
                
        String g_id = fa.getLimsValue(Lims.GROUP_ID);
        if (null != g_id && !g_id.equalsIgnoreCase("NA")) {
          run.addProperty("group_id", g_id);
        } else {
          run.addProperty("group_id", "NA");    
        }
        
        String g_id_desc = fa.getLimsValue(Lims.GROUP_DESC);
        if (null != g_id_desc && !g_id_desc.equalsIgnoreCase("NA")) {
          run.addProperty("group_id_description", g_id_desc);
        } else {
          run.addProperty("group_id_description", "NA");    
        }
        
        if (null != tubeId && !tubeId.equalsIgnoreCase("NA")) {
          run.addProperty("external_name", tubeId);
        } else {
          run.addProperty("external_name", "NA");    
        }
        
        return new ReturnValue();
    }

   
   public static void main(String args[]){
 
        List<String> params = new ArrayList<String>();
        params.add("--plugin");
        params.add(StructuralVariationDecider.class.getCanonicalName());
        params.add("--");
        params.addAll(Arrays.asList(args));
        System.out.println("Parameters: " + Arrays.deepToString(params.toArray()));
        net.sourceforge.seqware.pipeline.runner.PluginRunner.main(params.toArray(new String[params.size()]));
         
    }
   
   private class BeSmall {

        private Date   date = null;
        private String iusDetails = null;
        private String groupByAttribute = null;
        private String path = null;
        
        public BeSmall(ReturnValue rv) {
            try {
                date = format.parse(rv.getAttribute(Header.PROCESSING_DATE.getTitle()));
            } catch (ParseException ex) {
                Log.error("Bad date!", ex);
                ex.printStackTrace();
            }
            FileAttributes fa = new FileAttributes(rv, rv.getFiles().get(0));
            //TODO check if metatype is doing what it is supposed to do (adding MIME type here)
            iusDetails = fa.getLibrarySample() + fa.getSequencerRun() + fa.getLane() + fa.getBarcode() + fa.getMetatype();
            
            //groupByAttribute = fa.getDonor() + ":" + fa.getLimsValue(Lims.TISSUE_ORIGIN) + ":" + fa.getLimsValue(Lims.LIBRARY_TEMPLATE_TYPE);
            groupByAttribute = fa.getLibrarySample();
            
            if (null != fa.getLimsValue(Lims.TISSUE_TYPE))
                 groupByAttribute = groupByAttribute.concat(":" + fa.getLimsValue(Lims.TISSUE_TYPE));
            
            if (null != fa.getLimsValue(Lims.TISSUE_PREP))
                 groupByAttribute = groupByAttribute.concat(":" + fa.getLimsValue(Lims.TISSUE_PREP));
            
            if (null != fa.getLimsValue(Lims.TISSUE_REGION))
                 groupByAttribute = groupByAttribute.concat(":" + fa.getLimsValue(Lims.TISSUE_REGION)); 
            
            if (null != fa.getLimsValue(Lims.GROUP_ID)) 
                groupByAttribute = groupByAttribute.concat(":" + fa.getLimsValue(Lims.GROUP_ID));
            
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
    }
    
}
