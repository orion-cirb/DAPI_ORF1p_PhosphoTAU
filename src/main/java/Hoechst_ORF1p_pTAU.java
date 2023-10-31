
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
import java.util.HashMap;
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
 * Detect Hoechst nuclei, ORF1p cells and pTAU cells and compute their colocalization
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
            String header = "Image name\tNuc label\tNuc volume\tNuc circularity\tNuc ORF1p int mean\tNuc pTAU int mean\t" +
                     "is ORF1p?\tNeuN volume\tNeuN int mean\tNeuN ORF1p int mean\t" +
                     "is ORF1p?\tORF1p volume\tORF1p int mean\tORF1p NeuN int mean\n";
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
                tools.print("- Analyzing DAPI channel -");
                int indexCh = ArrayUtils.indexOf(channels, channelChoices[0]);
                ImagePlus imgHoechst = BF.openImagePlus(options)[indexCh];
                // Detect Hoechst nuclei with CellPose
                Objects3DIntPopulation dapiPop = tools.cellposeDetection(imgHoechst, true, tools.cellposeNucModel, 1, tools.cellposeNucDiam, tools.cellposeNucStitchThresh, tools.minNucVol, tools.maxNucVol);

                // Open NeuN channel
                tools.print("- Analyzing NeuN channel -");
                indexCh = ArrayUtils.indexOf(channels, channelChoices[1]);
                ImagePlus imgNeun = BF.openImagePlus(options)[indexCh];
                // Detect NeuN cells with CellPose
                Objects3DIntPopulation neunPop = tools.cellposeDetection(imgNeun, true, tools.cellposeNeunModel, 1, tools.cellposeNeunDiam, tools.cellposeNeunStitchThresh, tools.minNeunVol, tools.maxNeunVol);
                
                // Open ORF1p channel
                tools.print("- Analyzing ORF1p channel -");
                indexCh = ArrayUtils.indexOf(channels, channelChoices[2]);
                ImagePlus imgOrf1p = BF.openImagePlus(options)[indexCh];
                // Detect NeuN cells with CellPose
                Objects3DIntPopulation orf1pPop = tools.cellposeDetection(imgOrf1p, true, tools.cellposeOrf1pModel, 1, tools.cellposeOrf1pDiam, tools.cellposeOrf1pStitchThresh, tools.minOrf1pVol, tools.maxOrf1pVol);
               
                tools.print("- Colocalizing nuclei with NeuN and ORF1p cells -");
                ArrayList<Cell> cells = tools.colocalization(dapiPop, neunPop, orf1pPop);
                
                tools.print("- Computing backgrounds statistics -");
                HashMap<String, Double> bgStats = tools.getBackgroundStats(imgHoechst, imgNeun, imgOrf1p);
                
                tools.print("- Measuring cells parameters -");
                tools.writeCellsParameters(cells, imgHoechst, imgNeun, imgOrf1p, bgStats);
                
                // Draw results
                tools.print("- Saving results -");
                tools.drawResults(imgHoechst, cells, parentFolder, rootName, outDirResults);
                
                // Write results
                for (Cell cell : cells) {
                    cellsResults.write(parentFolder+"\t"+rootName+"\t"+cell.params.get("label")+"\t"+cell.params.get("nucVol")+"\t"+cell.params.get("nucCirc")+
                        "\t"+cell.params.get("nucIntMean")+"\t"+cell.params.get("nucIntSd")+"\t"+cell.params.get("nucIntMeanNeun")+
                        "\t"+cell.params.get("nucIntSdNeun")+"\t"+cell.params.get("nucIntMeanOrf1p")+"\t"+cell.params.get("nucIntSdOrf1p")+
                        "\t"+(cell.neun!=null)+"\t"+cell.params.get("neunVol")+"\t"+cell.params.get("neunIntMean")+"\t"+cell.params.get("neunIntSd")+
                        "\t"+cell.params.get("neunIntMeanOrf1p")+"\t"+cell.params.get("neunIntSdOrf1p")+"\t"+cell.params.get("neunCytoVol")+
                        "\t"+cell.params.get("neunCytoIntMean")+"\t"+cell.params.get("neunCytoIntSd")+"\t"+cell.params.get("neunCytoIntMeanOrf1p")+
                        "\t"+cell.params.get("neunCytoIntSdOrf1p")+"\t"+(cell.orf1p!=null)+"\t"+cell.params.get("orf1pVol")+"\t"+cell.params.get("orf1pIntMean")+
                        "\t"+cell.params.get("orf1pIntSd")+"\t"+cell.params.get("orf1pIntMeanNeun")+"\t"+cell.params.get("orf1pIntSdNeun")+
                        "\t"+cell.params.get("orf1pCytoVol")+"\t"+cell.params.get("orf1pCytoIntMean")+"\t"+cell.params.get("orf1pCytoIntSd")+
                        "\t"+cell.params.get("orf1pCytoIntMeanNeun")+"\t"+cell.params.get("orf1pCytoIntSdNeun")+"\n");
                    cellsResults.flush();
                }
                
                tools.closeImg(imgHoechst);
                tools.closeImg(imgNeun);
                tools.closeImg(imgOrf1p);
            }
        } catch (IOException | DependencyException | ServiceException | FormatException ex) {
            Logger.getLogger(Hoechst_ORF1p_pTAU.class.getName()).log(Level.SEVERE, null, ex);
        }
        tools.print("--- Analysis done ---");
    }    
}    
