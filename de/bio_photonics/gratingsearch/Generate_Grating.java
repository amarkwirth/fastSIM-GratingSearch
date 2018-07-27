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
import javax.swing.JFileChooser;

import org.fairsim.linalg.Vec2d;
import org.fairsim.fiji.DisplayWrapper;
import org.fairsim.fiji.ImageVector;

public class Generate_Grating implements ij.plugin.PlugIn {

    public void run(String arg) {

	if (arg=="clear") {
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
// 	gd.addNumericField("SLM width",1280,0);
// 	gd.addNumericField("SIM height",1024,0);
// 	gd.addStringField("file prefix","automatic", 30);
	gd.showDialog();
	if (gd.wasCanceled())
	    return;
    
	// get and check parameters
	int nr = (int)gd.getNextNumber();
// 	int width   = (int)gd.getNextNumber();
// 	int height  = (int)gd.getNextNumber();
	String prefix = (String)IJ.getProperty("de.bio_photonics.gratingsearch.prefix");

	int width = (Integer)IJ.getProperty("de.bio_photonics.gratingsearch.width");
	int height = (Integer)IJ.getProperty("de.bio_photonics.gratingsearch.height");
//     OpenDialog directory = new ij.io.OpenDialog("Choose a directory"); 
// 	String path = directory.getDirectory();
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

//         int test1 = 0;
//         int test2 = 0;
//         for ( Grating [] gr : gratList.get(0)) {
//         test1++;
//             for ( Grating gri : gr ) {
//                 test2++;
//                 System.out.println("test1:" + test1 + " test2:" + test2);
//             }
//         }
        
        Grating [][] gr = gratList.get(nr);
        
        int nrAngles            = (Integer)IJ.getProperty("de.bio_photonics.gratingsearch.angNumber") ;
        double [] resImp        =  (double[])IJ.getProperty("de.bio_photonics.gratingsearch.resImp");
        double [] wavelength    = (double[])IJ.getProperty("de.bio_photonics.gratingsearch.wl");
        double [] NA            = (double[])IJ.getProperty("de.bio_photonics.gratingsearch.na");

        for(int r = 0; r<resImp.length; r++) {
            for(int wl=0; wl<wavelength.length; wl++) {
                for (int ang =0; ang<nrAngles; ang++) {
                    for (int pha =0; pha<nrPhases; pha++) {
                        double phase = pha*Math.PI*2/nrPhases;
                        ImageVector pttr = ImageVector.create(width,height);
                        System.out.println("gr["+(wl+r*wavelength.length)+"]["+ang+"]");
                        Grating gri=gr[wl+r*wavelength.length][ang];
                        gri.writeToVector( pttr , phase );
                        String title = String.format("%s", prefix);
                        title += String.format("_NA:%1.2f", NA[r]);
                        title += String.format("_res:%1.2f", resImp[r]);
                        title += String.format("_wl:%.0f", wavelength[wl]);
                        title += String.format("_ang:%d", ang);
                        title += String.format("_pha:%d", pha);
                        img.addImage( pttr, title);
                        if(returnVal == JFileChooser.APPROVE_OPTION) {
                            ImagePlus tmpIP = new ImagePlus("test", pttr.img());
                            IJ.save(tmpIP, String.format("%s/%s_NA%1.2f_res%1.2f_wl%.0f_ang%d_pha%d.bmp", path, prefix, NA[r], resImp[r], wavelength[wl], ang, pha));
                        }
                    }
                }
            }
        }
//         img.display();
    }
}
