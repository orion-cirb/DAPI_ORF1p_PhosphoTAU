import Hoechst_ORF1p_pTAU_Tools.Cell;
import Hoechst_ORF1p_pTAU_Tools.Tools;
import ij.*;
import ij.plugin.PlugIn;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.logging.Level;
import java.util.logging.Logger;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.common.services.ServiceFactory;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.formats.services.OMEXMLService;
import loci.plugins.BF;
import loci.plugins.util.ImageProcessorReader;
import loci.plugins.in.ImporterOptions;
import mcib3d.geom2.Objects3DIntPopulation;
import org.apache.commons.io.FilenameUtils;
import org.scijava.util.ArrayUtils;


/**
 * Detect Hoechst nuclei and measure their intensity in ORF1p channel
 * Detect pTAU cells and colocalize them with Hoechst nuclei
 * @author ORION-CIRB
 */
public class Hoechst_ORF1p_pTAU implements PlugIn {
    
    Tools tools = new Tools();   
    
    public void run(String arg) {
        try {
            if ((!tools.checkInstalledModules())) {
                return;
            } 
            
            String imageDir = IJ.getDirectory("Choose directory containing image files...");
            if (imageDir == null) {
                return;
            }
            
            // Find images with fileExt extension
            String fileExt = tools.findImageType(new File(imageDir));
            ArrayList<String> imageFiles = tools.findImages(imageDir, fileExt);
            if (imageFiles == null) {
                IJ.showMessage("Error", "No images found with " + fileExt + " extension");
                return;
            }
            
            // Create OME-XML metadata store of the latest schema version
            ServiceFactory factory;
            factory = new ServiceFactory();
            OMEXMLService service = factory.getInstance(OMEXMLService.class);
            IMetadata meta = service.createOMEXMLMetadata();
            ImageProcessorReader reader = new ImageProcessorReader();
            reader.setMetadataStore(meta);
            reader.setId(imageFiles.get(0));
            
            // Find image calibration
            tools.findImageCalib(meta);
            
            // Find channel names
            String[] channels = tools.findChannels(imageFiles.get(0), meta, reader);
            
            // Generate dialog box
            String[] channelChoices = tools.dialog(channels);
            if (channels == null) {
                IJ.showStatus("Plugin canceled");
                return;
            }
            
            // Create output folder
            String outDirResults = imageDir + File.separator+ "Results_" + new SimpleDateFormat("yyyy-MM-dd_HH-mm-ss").format(new Date()) + File.separator;
            File outDir = new File(outDirResults);
            if (!Files.exists(Paths.get(outDirResults))) {
                outDir.mkdir();
            }
            
            // Write header in results file
            String header = "Image name\tNucleus label\tNucleus volume (µm3)\tNucleus circularity\t" +
                     "ORF1p bg\tNucleus ORF1p bg-corr int mean\tNucleus ORF1p int sd\tis pTAU?\n";
            FileWriter fwCellsResults = new FileWriter(outDirResults + "results.xls", false);
            BufferedWriter cellsResults = new BufferedWriter(fwCellsResults);
            cellsResults.write(header);
            cellsResults.flush();

            for (String f : imageFiles) {
                String rootName = FilenameUtils.getBaseName(f);
                tools.print("--- ANALYZING IMAGE " + rootName + " ------");
                reader.setId(f);
                
                ImporterOptions options = new ImporterOptions();
                options.setId(f);
                options.setSplitChannels(true);
                options.setQuiet(true);
                options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                
                // Open Hoechst channel
                tools.print("- Analyzing Hoechst channel -");
                int indexCh = ArrayUtils.indexOf(channels, channelChoices[0]);
                ImagePlus imgHoechst = BF.openImagePlus(options)[indexCh];
                // Detect Hoechst nuclei with CellPose
                Objects3DIntPopulation nucleiPop = tools.cellposeDetection(imgHoechst, tools.cellposeNucModel, tools.cellposeNucDiam, tools.cellposeNucStitchThresh, tools.minNucVol, tools.maxNucVol);

                // Open ORF1p channel
                tools.print("- Analyzing ORF1p channel -");
                indexCh = ArrayUtils.indexOf(channels, channelChoices[1]);
                ImagePlus imgOrf1p = BF.openImagePlus(options)[indexCh];
                // Measure ORF1p channel background
                double bgOrf1p = tools.measureBackgroundNoise(imgOrf1p);

                // Open pTAU channel
                tools.print("- Analyzing pTAU channel -");
                indexCh = ArrayUtils.indexOf(channels, channelChoices[2]);
                ImagePlus imgPtau = BF.openImagePlus(options)[indexCh];
                // Detect pTAU cells with CellPose
                Objects3DIntPopulation ptauPop = tools.cellposeDetection(imgPtau, tools.cellposeModelPath+tools.cellposePtauModel, tools.cellposePtauDiam, tools.cellposePtauStitchThresh, tools.minPtauVol, tools.maxPtauVol);
               
                tools.print("- Colocalizing nuclei with pTAU cells -");
                ArrayList<Cell> cells = tools.colocalization(nucleiPop, ptauPop);
                
                tools.print("- Measuring cells parameters -");
                tools.writeCellsParameters(cells, imgOrf1p, bgOrf1p);
                
                // Draw results
                tools.print("- Saving results -");
                tools.drawResults(cells, imgHoechst, imgPtau, outDirResults+rootName);
                
                // Write results
                for (Cell cell : cells) {
                    cellsResults.write(rootName+"\t"+cell.params.get("label")+"\t"+cell.params.get("nucVol")+"\t"+cell.params.get("nucCirc")+
                        "\t"+bgOrf1p+"\t"+cell.params.get("nucOrf1pIntMean")+"\t"+cell.params.get("nucOrf1pIntSd")+"\t"+cell.isPtau+"\n");
                    cellsResults.flush();
                }
                
                tools.closeImg(imgHoechst);
                tools.closeImg(imgOrf1p);
                tools.closeImg(imgPtau);
            }
        } catch (IOException | DependencyException | ServiceException | FormatException ex) {
            Logger.getLogger(Hoechst_ORF1p_pTAU.class.getName()).log(Level.SEVERE, null, ex);
        }
        tools.print("--- Analysis done ---");
    }    
}    
