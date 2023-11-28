package Hoechst_ORF1p_pTAU_Tools;

import java.util.HashMap;
import mcib3d.geom2.Object3DInt;

/**
 * @author ORION-CIRB
 */
public class Cell {
    
    // Nucleus
    public Object3DInt nucleus;
    // pTAU cell
    public boolean isPtau;
    public Object3DInt ptau;
    // Parameters
    public HashMap<String, Double> params;
    

    public Cell(Object3DInt nucleus) {
        this.nucleus = nucleus;
        this.ptau = null;
        this.isPtau = false;
        this.params = new HashMap<>();
    }
    
    public void setLabel(double label) {
        params.put("label", label);
    }
    
    public void setNucParams(double nucVol, double nucCirc, double nucOrf1pIntMean, double nucOrf1pIntSd) {
        params.put("nucVol", nucVol);
        params.put("nucCirc", nucCirc);
        params.put("nucOrf1pIntMean", nucOrf1pIntMean);
        params.put("nucOrf1pIntSd", nucOrf1pIntSd);
    }
    
    public void setPtau(Object3DInt ptau) {
        this.ptau = ptau;
        this.isPtau = true;
    }
}
