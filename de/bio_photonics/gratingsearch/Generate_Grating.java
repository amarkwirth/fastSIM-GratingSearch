/*
This file is part of fastSIM Grating Search, ported to Java.

This code is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or
(at your option) any later version.

The code is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with the code.  If not, see <http://www.gnu.org/licenses/>
*/

package de.bio_photonics.gratingsearch;

import java.util.List;
import java.util.ArrayList;

import ij.IJ;
import ij.gui.GenericDialog;
import ij.ImagePlus;
import java.io.File;
import javax.swing.JFileChooser;

import org.fairsim.linalg.Vec2d;
import org.fairsim.fiji.DisplayWrapper;
import org.fairsim.fiji.ImageVector;

public class Generate_Grating implements ij.plugin.PlugIn {

    String mktitle(String prefix, int saveNr, double NA, double [] resImp, double [] wavelength, int r, int wl, int ang, int pha) {
        String title = String.format("%s", prefix);
        title += "_"+String.format("%02d",   saveNr)+"_";
        title += String.format("_NA%1.2f", NA);
        title += "_res";
        for (int re=0; re<resImp.length; re++) {
            if(re>0) title+="x";
            title += String.format("%1.2f", resImp[re]);
        }
        title += String.format("_wl%.0f", wavelength[wl]);
        title += String.format("_ang%d", ang);
        title += String.format("_pha%d", pha);
        return title;
    }
    

    String mkfolder(String prefix, int saveNr, double NA, double [] resImp, double [] wavelength) {
        String folder = String.format("%s", prefix);
        folder += "_"+String.format("%02d", saveNr)+"_";
        folder += "_na";
        folder += String.format("%1.2f", NA);
        folder += "_res";
        for (int re=0; re<resImp.length; re++) {
            if(re>0) folder+="x";
            folder += String.format("%1.2f", resImp[re]);
        }
        folder += "_wl";
        for (int wl=0; wl<wavelength.length; wl++) {
            if(wl>0) folder+="x";
            folder+= String.format("%.0f", wavelength[wl]);
        }
        return folder;
    }

    public void run(String arg) {

    if (arg.equals("clear")) {
        IJ.log("clearing stored gratings");
        IJ.setProperty("de.bio_photonics.gratingsearch.phaseNumber",null);
        IJ.setProperty("de.bio_photonics.gratingsearch.lastGratings",null);
        return;
    }

    // check / retrieve parameters
    if ( IJ.getProperty("de.bio_photonics.gratingsearch.phaseNumber")==null ||
        IJ.getProperty("de.bio_photonics.gratingsearch.angNumber")==null ||
        IJ.getProperty("de.bio_photonics.gratingsearch.resImp")==null ||
        IJ.getProperty("de.bio_photonics.gratingsearch.wl")==null ||
        IJ.getProperty("de.bio_photonics.gratingsearch.na")==null ||
        IJ.getProperty("de.bio_photonics.gratingsearch.lastGratings")==null ||
        IJ.getProperty("de.bio_photonics.gratingsearch.width")==null ||
        IJ.getProperty("de.bio_photonics.gratingsearch.height")==null ||
        IJ.getProperty("de.bio_photonics.gratingsearch.prefix")==null ||
        IJ.getProperty("de.bio_photonics.gratingsearch.resImp")==null 
        ) {
        IJ.log("No pattern information stored, search for pattern first!");
        return;
    }

    // TODO: Check for cast errors?
    int nrPhases = (Integer)IJ.getProperty("de.bio_photonics.gratingsearch.phaseNumber");

    Object a= IJ.getProperty("de.bio_photonics.gratingsearch.lastGratings");
    List tmp = ( a instanceof List )?((List)a):(null);
    List<Grating [][]> gratList = new ArrayList<Grating [][]>();

    for ( Object b : tmp ) 
        if ( b instanceof Grating[][] ) gratList.add((Grating [][])b);    
    

    if (gratList !=null && gratList.size()==0) {
        IJ.log("No pattern in list, run search again");
        return;
    }

    // show dialog
    GenericDialog gd = new GenericDialog("Pattern generation");
    gd.addNumericField(String.format("Set nr [0 - %d]",gratList.size()-1),0,0);
    gd.addCheckbox("Save all", false);

    gd.showDialog();
    if (gd.wasCanceled())
        return;
    
    // get and check parameters
    int nr = (int)gd.getNextNumber();
    Boolean saveAll = gd.getNextBoolean();
    String prefix = (String)IJ.getProperty("de.bio_photonics.gratingsearch.prefix");

    int width = (Integer)IJ.getProperty("de.bio_photonics.gratingsearch.width");
    int height = (Integer)IJ.getProperty("de.bio_photonics.gratingsearch.height");
    final String path;
    final JFileChooser fc = new JFileChooser();
        fc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
        fc.setAcceptAllFileFilterUsed(false);
        fc.setDialogTitle("choose directory");
        int returnVal = fc.showSaveDialog(null);
        if(returnVal == JFileChooser.APPROVE_OPTION) {
        path = fc.getSelectedFile().getAbsolutePath();
    }
        else path="";
    

        // generate pattern
        DisplayWrapper img = new DisplayWrapper(width, height, "Pattern");

    int saveStart = 0; int saveStop = gratList.size()-1;
    if(saveAll.equals(false)){
        saveStart = nr;
        saveStop = nr;
    }
        for(int saveNr = saveStart; saveNr <= saveStop; saveNr++) {
            Grating [][] gr = gratList.get(saveNr);
            
            int nrAngles            = (Integer)IJ.getProperty("de.bio_photonics.gratingsearch.angNumber") ;
            double [] resImp        =  (double[])IJ.getProperty("de.bio_photonics.gratingsearch.resImp");
            double [] wavelength    = (double[])IJ.getProperty("de.bio_photonics.gratingsearch.wl");
            double NA            = (Double)IJ.getProperty("de.bio_photonics.gratingsearch.na");

            for(int r = 0; r<resImp.length; r++) {
                for(int wl=0; wl<wavelength.length; wl++) {
                    for (int ang =0; ang<nrAngles; ang++) {
                        for (int pha =0; pha<nrPhases; pha++) {
                            double phase = pha*Math.PI*2/nrPhases;
                            ImageVector pttr = ImageVector.create(width,height);
                            Grating gri=gr[ang+r*nrAngles][wl];
                            gri.writeToVector( pttr , phase );
                            
                            String title = mktitle(prefix, saveNr, NA, resImp, wavelength, r, wl, ang, pha);
                            String folder = mkfolder(prefix, saveNr, NA, resImp, wavelength);
                            
                            System.out.print("gr["+(wl+r*wavelength.length)+"]["+ang+"] -> ");
                            System.out.println(title);
                            img.addImage( pttr, title);
                            if(returnVal == JFileChooser.APPROVE_OPTION) {
                                ImagePlus tmpIP = new ImagePlus("test", pttr.img());
                                String directory= path+File.separator+folder+File.separator;
                                new File(directory).mkdirs();
                                String absoluteFilename = directory+title+".bmp";
                                IJ.save(tmpIP, absoluteFilename);
                            }
                        }
                    }
                }
            }
            if(saveAll.equals(false)){
                img.display();
            }
        }
    }
}
