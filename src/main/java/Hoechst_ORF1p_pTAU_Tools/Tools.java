package Hoechst_ORF1p_pTAU_Tools;


import Hoechst_ORF1p_pTAU_Tools.Cellpose.CellposeSegmentImgPlusAdvanced;
import Hoechst_ORF1p_pTAU_Tools.Cellpose.CellposeTaskSettings;
import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.measure.Measurements;
import ij.measure.ResultsTable;
import ij.plugin.Duplicator;
import ij.plugin.RGBStackMerge;
import ij.plugin.ZProjector;
import ij.plugin.filter.ParticleAnalyzer;
import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import javax.swing.ImageIcon;
import loci.common.services.ServiceException;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom2.Object3DInt;
import mcib3d.geom2.Object3DPlane;
import mcib3d.geom2.Objects3DIntPopulation;
import mcib3d.geom2.Objects3DIntPopulationComputation;
import mcib3d.geom2.measurements.MeasureIntensity;
import mcib3d.geom2.measurements.MeasureVolume;
import mcib3d.geom2.measurementsPopulation.MeasurePopulationColocalisation;
import mcib3d.image3d.ImageHandler;
import org.apache.commons.io.FilenameUtils;


/**
 * @author ORION-CIRB
 */
public class Tools {
    private final ImageIcon icon = new ImageIcon(this.getClass().getResource("/Orion_icon.png"));
    private final String helpUrl = "https://github.com/orion-cirb/Hoechst_ORF1p_pTAU.git";
    
    public String[] channelNames = {"Hoechst", "ORF1p", "pTAU"};
    public Calibration cal;
    public double pixelVol;
    
    // Cellpose
    private String cellposeEnvDirPath = IJ.isWindows()? System.getProperty("user.home")+"\\miniconda3\\envs\\CellPose" : "/opt/miniconda3/envs/cellpose";
    public final String cellposeModelPath = IJ.isWindows()? System.getProperty("user.home")+"\\.cellpose\\models\\" : "";
    // Nuclei detection
    public String cellposeNucModel = "cyto";
    public int cellposeNucDiam = 100;
    public double cellposeNucStitchThresh = 0.5;
    public double minNucVol = 200;
    public double maxNucVol = 6000;
    // pTAU detection
    public String cellposePtauModel = "cyto2_pTAU";
    public int cellposePtauDiam = 140;
    public double cellposePtauStitchThresh = 0.5;
    public double minPtauVol = 1000;
    public double maxPtauVol = 20000;
    
    
    
    /**
     * Display a message in the ImageJ console and status bar
     */
    public void print(String log) {
        System.out.println(log);
        IJ.showStatus(log);
    }
    
    
    /**
     * Check that needed modules are installed
     */
    public boolean checkInstalledModules() {
        ClassLoader loader = IJ.getClassLoader();
        try {
            loader.loadClass("mcib3d.geom.Object3D");
        } catch (ClassNotFoundException e) {
            IJ.log("3D ImageJ Suite not installed, please install from update site");
            return false;
        }
        return true;
    }
    
    
    /**
     * Find images extension
     */
    public String findImageType(File imagesFolder) {
        String ext = "";
        String[] files = imagesFolder.list();
        for (String name : files) {
            String fileExt = FilenameUtils.getExtension(name);
            switch (fileExt) {
               case "nd" :
                   ext = fileExt;
                   break;
                case "czi" :
                   ext = fileExt;
                   break;
                case "lif"  :
                    ext = fileExt;
                    break;
                case "ics" :
                    ext = fileExt;
                    break;
                case "ics2" :
                    ext = fileExt;
                    break;
                case "lsm" :
                    ext = fileExt;
                    break;
                case "tif" :
                    ext = fileExt;
                    break;
                case "tiff" :
                    ext = fileExt;
                    break;
            }
        }
        return(ext);
    }
    
    
    /**
     * Find images in folder
     */
    public ArrayList<String> findImages(String imagesFolder, String imageExt) {
        File inDir = new File(imagesFolder);
        String[] files = inDir.list();
        if (files == null) {
            System.out.println("No image found in " + imagesFolder);
            return null;
        }
        ArrayList<String> images = new ArrayList();
        for (String f : files) {
            String fileExt = FilenameUtils.getExtension(f);
            if (fileExt.equals(imageExt) && !f.startsWith("."))
                images.add(imagesFolder + File.separator + f);
        }
        Collections.sort(images);
        return(images);
    }
    
    
    /**
     * Find image calibration
     */
    public void findImageCalib(IMetadata meta) {
        cal = new Calibration();
        cal.pixelWidth = meta.getPixelsPhysicalSizeX(0).value().doubleValue();
        cal.pixelHeight = cal.pixelWidth;
        if (meta.getPixelsPhysicalSizeZ(0) != null)
            cal.pixelDepth = meta.getPixelsPhysicalSizeZ(0).value().doubleValue();
        else
            cal.pixelDepth = 1;
        cal.setUnit("microns");
        System.out.println("XY calibration = " + cal.pixelWidth + ", Z calibration = " + cal.pixelDepth);
    }
    
    
    /**
     * Find channels name
     * @throws loci.common.services.DependencyException
     * @throws loci.common.services.ServiceException
     * @throws loci.formats.FormatException
     * @throws java.io.IOException
     */
    public String[] findChannels(String imageName, IMetadata meta, ImageProcessorReader reader) throws loci.common.services.DependencyException, ServiceException, FormatException, IOException {
        int chs = reader.getSizeC();
        String[] channels = new String[chs];
        String imageExt =  FilenameUtils.getExtension(imageName);
        switch (imageExt) {
            case "nd" :
                for (int n = 0; n < chs; n++) 
                {
                    if (meta.getChannelID(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelName(0, n);
                }
                break;
            case "nd2" :
                for (int n = 0; n < chs; n++) 
                {
                    if (meta.getChannelID(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelName(0, n);
                }
                break;
            case "lif" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null || meta.getChannelName(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelName(0, n);
                break;
            case "czi" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelFluor(0, n);
                break;
            case "ics" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelExcitationWavelength(0, n).value().toString();
                break;
            case "ics2" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelExcitationWavelength(0, n).value().toString();
                break;   
            default :
                for (int n = 0; n < chs; n++)
                    channels[n] = Integer.toString(n);
        }
        return(channels);         
    }
    
    
    /**
     * Generate dialog box
     */
    public String[] dialog(String[] channels) {
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        gd.setInsets​(0, 120, 0);
        gd.addImage(icon);
        
        gd.addMessage("Channels", Font.getFont("Monospace"), Color.blue);
        int index = 0;
        for (String chName: channelNames) {
            gd.addChoice(chName+" : ", channels, channels[index]);
            index++;
        }
        
        gd.addMessage("Nuclei detection", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("Min nucleus volume (µm3): ", minNucVol);
        gd.addNumericField("Max nucleus volume (µm3): ", maxNucVol);
        
        gd.addMessage("pTAU cells detection", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("Min cell volume (µm3): ", minPtauVol);
        gd.addNumericField("Max cell volume (µm3): ", maxPtauVol);
        
        gd.addMessage("Image calibration", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("XY pixel size (µm): ", cal.pixelWidth);
        gd.addNumericField("Z pixel depth (µm):", cal.pixelDepth);
        gd.addHelp(helpUrl);
        gd.showDialog();
        
        String[] chChoices = new String[channelNames.length];
        for (int n = 0; n < chChoices.length; n++) {
            chChoices[n] = gd.getNextChoice();
        }
        
        minNucVol = gd.getNextNumber();
        maxNucVol = gd.getNextNumber();
        minPtauVol = gd.getNextNumber();
        maxPtauVol = gd.getNextNumber();
        cal.pixelWidth = cal.pixelHeight = gd.getNextNumber();
        cal.pixelDepth = gd.getNextNumber();
        pixelVol = cal.pixelWidth*cal.pixelWidth*cal.pixelDepth;
        
        if (gd.wasCanceled())
            chChoices = null;
        
        return(chChoices);
    }
    
    
    /**
     * Flush and close an image
     */
    public void closeImg(ImagePlus img) {
        img.flush();
        img.close();
    }
    
    
    /**
     * Look for all 3D cells in a z-stack: 
     * - apply CellPose in 2D slice by slice 
     * - let CellPose reconstruct cells in 3D using the stitch threshold parameter
     */
    public Objects3DIntPopulation cellposeDetection(ImagePlus img, String cellposeModel, int diameter, double stitchThreshold, double volMin, double volMax) throws IOException{
        ImagePlus imgIn = new Duplicator().run(img);

        // Define CellPose settings
        CellposeTaskSettings settings = new CellposeTaskSettings(cellposeModel, 1, diameter, cellposeEnvDirPath);
        settings.setStitchThreshold(stitchThreshold);
        settings.useGpu(true);
       
        // Run CellPose
        CellposeSegmentImgPlusAdvanced cellpose = new CellposeSegmentImgPlusAdvanced(settings, imgIn);
        ImagePlus imgOut = cellpose.run();
        imgOut.setCalibration(cal);
       
        // Get cells as a population of objects and filter them
        Objects3DIntPopulation pop = new Objects3DIntPopulation(ImageHandler.wrap(imgOut));
        Objects3DIntPopulation popFilter = new Objects3DIntPopulationComputation(pop).getExcludeBorders(ImageHandler.wrap(img), false);
        filterDetectionsByZ(popFilter);
        filterDetectionsBySize(popFilter, volMin, volMax);
        popFilter.resetLabels();
        System.out.println(popFilter.getNbObjects() + " detections remaining after filtering (" + (pop.getNbObjects()-popFilter.getNbObjects()) + " filtered out)");
        
        closeImg(imgIn);
        closeImg(imgOut);
        return(popFilter);
    }
    
    
    /**
     * Remove objects present in only one z slice
     */
    public void filterDetectionsByZ(Objects3DIntPopulation pop) {
        pop.getObjects3DInt().removeIf(p -> (p.getBoundingBox().zmax == p.getBoundingBox().zmin));           
    }
    
    
    /**
     * Remove objects with size < min and size > max
     */
    public void filterDetectionsBySize(Objects3DIntPopulation pop, double min, double max) {
        pop.getObjects3DInt().removeIf(p -> (new MeasureVolume(p).getVolumeUnit() < min) || (new MeasureVolume(p).getVolumeUnit() > max));
    }
        
    
    /**
     * Measure image intensity background noise:
     * Z projection over min intensity + read median intensity
     * @param img
     */
    public double measureBackgroundNoise(ImagePlus img) {
      ImagePlus imgProj = doZProjection(img, ZProjector.MIN_METHOD);
      double bg = imgProj.getProcessor().getStatistics().median;
      System.out.println("Background (median of the min projection) = " + bg);
      closeImg(imgProj);
      return(bg);
    }
    
    
    /**
     * Z-project a stack
     */
    public ImagePlus doZProjection(ImagePlus img, int param) {
        ZProjector zproject = new ZProjector();
        zproject.setMethod(param);
        zproject.setStartSlice(1);
        zproject.setStopSlice(img.getNSlices());
        zproject.setImage(img);
        zproject.doProjection();
       return(zproject.getProjection());
    }
    
    
    /**
     * Colocalize nuclei with a population of cells
     */
    public ArrayList<Cell> colocalization(Objects3DIntPopulation nucleiPop, Objects3DIntPopulation cellPop) {
        ArrayList<Cell> cells = new ArrayList<Cell>();
        if (nucleiPop.getNbObjects() > 0) {
            MeasurePopulationColocalisation coloc = new MeasurePopulationColocalisation(nucleiPop, cellPop);
            float label = 1;
            
            for (Object3DInt nucleus: nucleiPop.getObjects3DInt()) {
                Cell newCell = new Cell(nucleus);
             
                for (Object3DInt cell: cellPop.getObjects3DInt()) {
                    double colocVal = coloc.getValueObjectsPair(nucleus, cell);
                    if (colocVal > 0.25*nucleus.size()) {
                        newCell.setPtau(cell);
                        cellPop.removeObject(cell);
                        break;
                    }
                }
               
                newCell.setLabel(label);
                cells.add(newCell);
                label++;
            }
        }
        return(cells);
    }
    

    /**
     * Compute and save PV and PNN cells parameters
     */
    public void writeCellsParameters(ArrayList<Cell> cells, ImagePlus imgOrf1p, double bgOrf1p) { 
        ImageHandler imhOrf1p = ImageHandler.wrap(imgOrf1p);
        
        for (Cell cell : cells) {
            // Nucleus
            Object3DInt nucleus = cell.nucleus;
            double nucVol = new MeasureVolume(nucleus).getVolumeUnit();
            double nucCirc = computeNucleusCircularity(nucleus, imgOrf1p);
            double nucOrf1pIntMean = new MeasureIntensity(nucleus, imhOrf1p).getValueMeasurement(MeasureIntensity.INTENSITY_AVG)-bgOrf1p;
            double nucOrf1pIntSd = new MeasureIntensity(nucleus, imhOrf1p).getValueMeasurement(MeasureIntensity.INTENSITY_SD);
            cell.setNucParams(nucVol, nucCirc, nucOrf1pIntMean, nucOrf1pIntSd);
        }
    }
    
    
    public double computeNucleusCircularity(Object3DInt nucleus, ImagePlus img) {
        List<Object3DPlane> planes = nucleus.getObject3DPlanes();
        Object3DPlane middlePlane = planes.get((int)(0.5*((double)planes.size()-1)));
        ImageHandler planeImh = ImageHandler.wrap(IJ.createImage("", "8-bit black", img.getWidth(), img.getHeight(), img.getNSlices()));
        middlePlane.drawObject(planeImh, 255);
        ImagePlus planeImg = planeImh.getImagePlus();
        planeImg.setZ(middlePlane.getZPlane()+1);
        IJ.setAutoThreshold(planeImg, "Default dark no-reset");
        
        ResultsTable rt = new ResultsTable();
        ParticleAnalyzer pa = new ParticleAnalyzer(ParticleAnalyzer.CLEAR_WORKSHEET+ParticleAnalyzer.LIMIT, Measurements.SHAPE_DESCRIPTORS, rt, 0, Double.MAX_VALUE);
        pa.analyze(planeImg​);
        double circ = rt.getValue("Circ.", 0);
        
        closeImg(planeImg);
        return(circ);
    }
    
    
    /**
     * Save detected cells in image
     */
    public void drawResults(ArrayList<Cell> cells, ImagePlus imgHoechst, ImagePlus imgPtau, String imgName) {
        ImageHandler imgObj1 = ImageHandler.wrap(imgHoechst).createSameDimensions();
        ImageHandler imgObj2 = imgObj1.createSameDimensions();
        
        for(Cell cell: cells) {
            int label = cell.params.get("label").intValue();
            cell.nucleus.drawObject(imgObj1, label);
            if (cell.isPtau) {
                cell.ptau.drawObject(imgObj2, label);
            }
            
        }
        
        ImagePlus[] imgColors = {imgObj1.getImagePlus(), null, imgObj2.getImagePlus(), imgHoechst, imgPtau};
        ImagePlus imgObjects = new RGBStackMerge().mergeHyperstacks(imgColors, true);
        imgObjects.setCalibration(cal);
        FileSaver ImgObjectsFile = new FileSaver(imgObjects);
        ImgObjectsFile.saveAsTiff(imgName + ".tif");
        
        closeImg(imgObj1.getImagePlus());
        closeImg(imgObj2.getImagePlus());
        closeImg(imgObjects);
    }
} 
