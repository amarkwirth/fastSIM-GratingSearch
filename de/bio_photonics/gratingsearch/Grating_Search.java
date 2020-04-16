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

import java.util.ArrayList;
import java.util.List;

import org.fairsim.fiji.ImageVector;
import org.fairsim.linalg.Cplx;
import org.fairsim.linalg.Vec2d;
import org.fairsim.linalg.Transforms;
import org.fairsim.sim_algorithm.SimUtils;

import org.fairsim.fiji.DisplayWrapper;

import org.fairsim.utils.ImageDisplay;
import org.fairsim.utils.Tool;

import org.fairsim.utils.SimpleMT;

import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.IJ;
import java.util.Collections;

/**
 * Java port of fastSIM SLM grating search algorithm. The original MATLAB code
 * can be found here, please also cite their publication if you use this
 * software to create SIM gratings:
 *
 * https://github.com/nanoimaging/fastSIM_GratingSearchforSLM
 */
public class Grating_Search implements ij.plugin.PlugIn {

    // search value space
    final int axmax = 30, aymax = 30;
    final int bxmin = 2, bymin = 2;
    final int bxmax = 30, bymax = 30;

    // fourier space pxl
    final int fsPxl = 512;

    /**
     * Calculated possible gratings by iterating parameter space. Each parameter
     * set and resulting grating is tested to meet conditions: (a) grating
     * period within given limits, (b) grating is shiftable by a given number of
     * equi-distant phase steps.
     *
     * @param gratPerMin Minimal grating period, in #pixels
     * @param gratPerMax Maximal grating period, in #pixels
     * @param phaseSteps Number of phase steps
     * @param wavelength wavelength of the coherent illumination light, in nm 
     * @return List of gratings satisfying the conditions
     *
     */
    public List<Grating> calcGrat(
        final double gratPerMin,
        final double gratPerMax,
        final int phaseSteps,
        final double wavelength) {

        // List of grating candidates
        List<Grating> candidates = new ArrayList<Grating>(1000);

        // outer loop set: loop grating sizes
        for (int ax = 0; ax < axmax; ax++) {
            Tool.tell(String.format(
                "DONE: %d / %d, #candidates %d", ax, axmax,
                candidates.size()));
            for (int ay = -aymax; ay < aymax; ay++) {

                // inner loop set: loop grating shifts
                for (int bx = bxmin; bx < bxmax; bx++) {
                    for (int by = -bxmax; by < bymax; by++) {

                        if (Math.abs(by) < bymin) {
                            continue;
                        }

                        // create the grating, auto-computes parameters
                        Grating current = new Grating(ax, ay, bx, by, wavelength);

                        // check if the grating period is within bounds
                        if (current.gratPer < gratPerMin
                            || current.gratPer > gratPerMax) {
                            continue;
                        }

                        // check if n equid. pha. can be shifted horizontal
                        if (current.testStepHorizontal(phaseSteps)) {
                            current.shiftDir = 0;
                            fuzzyAdd(current, candidates);
                        }
                        // check if n equid. pha. can be shifted vertical
                        if (current.testStepVertical(phaseSteps)) {
                            current.shiftDir = 1;
                            fuzzyAdd(current, candidates);
                        }
                    }
                }
            }
        }

        return candidates;
    }

    /**
     * Add a grating only if no grating with similar direction is already
     * present
     */
    private static boolean fuzzyAdd(Grating a, List<Grating> candidates) {

        // Tolerance when looking for equal direction, in rad
        final double tolDir = 1e-8;

        for (Grating b : candidates) {
            if (Math.abs(b.gratDir - a.gratDir) < tolDir) {
                if (Math.abs(b.gratPer - a.gratPer) < tolDir) {
                    return false;
                }
            }
        }

        candidates.add(a);
        return true;
    }

    /**
     * Helper: returns the set with the lowest average euclidean distance
     */
    private static Grating[] lowestAvrEuclDist(Grating[] a, List< Grating[]> b) {
        double minAvr = Double.MAX_VALUE;
        Grating[] ret = null;

        for (Grating[] c : b) {
            double avr = 0;
            for (int l = 0; l < a.length; l++) {
                avr += Grating.scaledDistance(a[l], c[l]);
            }
            if (avr < minAvr) {
                ret = c;
                minAvr = avr;
            }
        }
        return ret;
    }

    /**
     * Select combinations of grating orientations. From a list of possible
     * gratings, generate subsets with equidistant angles, i.e. for n
     * directions, angle between gratings is approx. pi/n.
     *
     * @param candidates List of grating candidates
     * @param candidates_2 a second List of grating candidates
     * @param Tirf boolean: only find candidates with matching patterns ind "candidates_2"
     * @param nrDirs	Number of directions (typ. 3, 5, ...)
     * @param tolDirection Allowed deviation from ideal angle, in rad
     * @return List of matching grating combinations
     *
     */
    public List<Grating[]> selectDirs(
        List< Grating> candidates,
        List< Grating> candidates_2,
        boolean Tirf,
        final int nrDirs,
        final double tolDirection) {

        double sector = Math.PI / nrDirs;

        int count = 0;
        int maxCount = candidates.size();

        List<Grating[]> ret = new ArrayList<Grating[]>();

        // outer loop, select only dirs in the first sector
        for (Grating c : candidates) {

            count++;
            if (count % 100 == 0) {
                ij.IJ.showStatus(String.format(" %d / %d tested, found: %d ",
                    count, maxCount, ret.size()));
            }

            if (!c.testDirection(0, sector)) {
                continue;
            }

            // create the array
            Grating[] col = new Grating[nrDirs];
            Grating[] colTirf = new Grating[nrDirs];
            col[0] = c;
            colTirf[0] = c;

            // see if the other directions can be found within limits
            for (int i = 1; i < nrDirs; i++) {
                double minDist = tolDirection;
                for (Grating d : candidates) {
                    double dist = Math.abs(
                        Grating.wrapPi(d.gratDir) - sector * i
                        - Grating.wrapPi(c.gratDir));
                    if (dist < minDist) {
                        col[i] = d;
                        minDist = dist;
                    }
                }
                if (Tirf) {
                    double minDistTirf = tolDirection;
                    for (Grating d : candidates_2) {
                        double dist = Math.abs(
                            Grating.wrapPi(d.gratDir) - sector * i
                            - Grating.wrapPi(c.gratDir));
                        if (dist < minDistTirf) {
                            colTirf[i] = d;
                            minDistTirf = dist;
                        }
                    }
                }
            }

            // check if a complete set was found
            boolean complete = true;
            for (Grating e : col) {
                if (e == null) {
                    continue;
                }
            }
            for (Grating e : colTirf) {
                if (e == null) {
                    continue;
                }
            }

            ret.add(col);
        }
        return ret;
    }

    /**
     * For a combination of gratings, check if a mask could block unwanted
     * orders.
     *
     * @param candidates List of gratings to check
     * @param illum	Illumination vector to apply
     * @param maxUnwanted Maximum ratio of amplitude passing through mask
     * @param maskSize	Size of mask, in pixel
     * @param	outputAlsoFailed Add also the gratings that failed unwanted order
     * testing
     * @param	spatial	Receives spatial images of illuminated SLM, may be null
     * @param fourier	Receives Fourier images of illuminated SLM, may be null
     * @param name	String to add to the displayed image caption
     *
     */
    public boolean fourierCheck(Grating[] candidates, Vec2d.Real illum,
        double maxUnwated, int maskSize, boolean outputAlsoFailed,
        ImageDisplay spatial, ImageDisplay fourier, String name) {

        Vec2d.Real sumFreq = Vec2d.createReal(fsPxl, fsPxl);
        Vec2d.Real[] gratFreq = Vec2d.createArrayReal(candidates.length, fsPxl, fsPxl);

        // Compute the gratings Fourier space and sum them up
        for (int i = 0; i < candidates.length; i++) {
            Vec2d.Cplx tmp = Vec2d.createCplx(fsPxl, fsPxl);
//            System.out.println("i/candidates.length = "+i+"/"+candidates.length);
            candidates[i].writeToVector(tmp);	// get grating as vector 
            
            tmp.times(illum);			// apply illumination vector
            
            Transforms.fft2d(tmp, false);	// transform
            
            Transforms.swapQuadrant(tmp);	// quadrant-swap FFT result 
            
            gratFreq[i].copyMagnitude(tmp);
            
            sumFreq.add(gratFreq[i]);			// add to full vector
        }

        // Compute wanted and unwanted contributions
        double[] wanted = new double[candidates.length];
        double[] unwanted = new double[candidates.length];
        boolean ok = true;
        String magResult = "";

        for (int i = 0; i < candidates.length; i++) {

            double[] mPos = candidates[i].peakPos(fsPxl);

            wanted[i]
                = sumRegion(gratFreq[i], (int) (fsPxl / 2 + mPos[0]), (int) (fsPxl / 2 + mPos[1]), 8)
                + sumRegion(gratFreq[i], (int) (fsPxl / 2 - mPos[0]), (int) (fsPxl / 2 - mPos[1]), 8);

            for (int j = 0; j < candidates.length; j++) {
                if (i == j) {
                    continue;
                }
                unwanted[i]
                    += sumRegion(gratFreq[j], (int) (fsPxl / 2 + mPos[0]), (int) (fsPxl / 2 + mPos[1]), maskSize)
                    + sumRegion(gratFreq[j], (int) (fsPxl / 2 - mPos[0]), (int) (fsPxl / 2 - mPos[1]), maskSize);
            }

            //Tool.trace(String.format("%8.5f / %8.5f -> %8.5f", wanted[i], unwanted[i], 
            //unwanted[i]/wanted[i]));
            if ((unwanted[i] / wanted[i]) > maxUnwated) {
                ok = false;
            }
            magResult += String.format("%1d: %5.4f ", i, unwanted[i] / wanted[i]);

            // store unwanted / wanted
            candidates[i].unwantedMod = unwanted[i] / wanted[i];

        }

        // store spectrum (if display is set != null)
        if ((fourier != null) && (outputAlsoFailed || ok)) {

            // markers for positions
            ImageDisplay.Marker[] maskRings = new ImageDisplay.Marker[candidates.length * 4];
            for (int i = 0; i < candidates.length; i++) {
                double[] mPos = candidates[i].peakPos(fsPxl);
                maskRings[i * 4 + 0]
                    = new ImageDisplay.Marker(fsPxl / 2 + mPos[0], fsPxl / 2 + mPos[1], 8, 8, false);
                maskRings[i * 4 + 1]
                    = new ImageDisplay.Marker(fsPxl / 2 - mPos[0], fsPxl / 2 - mPos[1], 8, 8, false);
                maskRings[i * 4 + 2]
                    = new ImageDisplay.Marker(fsPxl / 2 + mPos[0], fsPxl / 2 + mPos[1], maskSize, maskSize, false);
                maskRings[i * 4 + 3]
                    = new ImageDisplay.Marker(fsPxl / 2 - mPos[0], fsPxl / 2 - mPos[1], maskSize, maskSize, false);
            }
            //fourier.addImage(magnitude( sumFreq ), 
            //name+" "+magResult+((ok)?("ok"):("!!not OK")), maskRings);
            fourier.addImage(sumFreq,
                name + " " + magResult + ((ok) ? ("ok") : ("!!not OK")), maskRings);
        }

        // store pattern (if display is set != null)
        if ((spatial != null) && (outputAlsoFailed || ok)) {
            for (int i = 0; i < candidates.length; i++) {
                Vec2d.Real img = Vec2d.createReal(fsPxl, fsPxl);
                candidates[i].writeToVector(img, 0);
                img.times(illum);
                spatial.addImage(img, name + " ang: " + i + " " + candidates[i]);
            }
        }

        return ok;
    }

    /**
     * sum up the magnitude within a sub-region (mask) of a vector
     */
    double sumRegion(Vec2d.Cplx vec, int xp, int yp, int size) {

        double sum = 0;

        for (int y = Math.max(0, yp - size); y < Math.min(vec.vectorHeight() - 1, yp + size); y++) {
            for (int x = Math.max(0, xp - size); x < Math.min(vec.vectorWidth() - 1, xp + size); x++) {
                sum += vec.get(x, y).abs();
            }
        }

        return sum;
    }

    /**
     * sum up the magnitude within a sub-region (mask) of a vector
     */
    double sumRegion(Vec2d.Real vec, int xp, int yp, int size) {

        double sum = 0;

        for (int y = Math.max(0, yp - size); y < Math.min(vec.vectorHeight() - 1, yp + size); y++) {
            for (int x = Math.max(0, xp - size); x < Math.min(vec.vectorWidth() - 1, xp + size); x++) {
                sum += vec.get(x, y);
            }
        }

        return sum;
    }

    /**
     * Create a Gaussian illumination profile.
     */
    public Vec2d.Real createGaussIllum(double fwhm, int size) {

        final double sigma = fwhm / 2.355;
        Vec2d.Real vec = Vec2d.createReal(size, size);

        for (int y = 0; y < size; y++) {
            for (int x = 0; x < size; x++) {

                double dist = Math.hypot(y - size / 2., x - size / 2.);
                double fac = Math.exp(-(dist * dist) / (2 * sigma * sigma));
                vec.set(x, y, (float) fac);
            }
        }
        return vec;
    }

    /**
     * Convolve a Fourier spectrum with the residual SLM structure
     */
    public void convSLMstructure(Vec2d.Cplx spec) {

        Transforms.fft2d(spec, true);
        Vec2d.Cplx residual = Vec2d.createCplx(spec);

        final int w = spec.vectorWidth();
        final int h = spec.vectorWidth();

        // structure directly from the original matlab script
        residual.set(0, h / 2 - 1, Cplx.Float.one().mult(0.08f));   // top
        residual.set(w - 1, h / 2 - 1, Cplx.Float.one().mult(0.08f));   // bottom

        residual.set(w / 4 - 1, 0, Cplx.Float.one().mult(0.03f));   // top left
        residual.set(w * 3 / 4 - 1, 0, Cplx.Float.one().mult(0.03f));   // bottom left

        residual.set(w / 4 - 1, h / 2 - 1, Cplx.Float.one().mult(0.03f));   // top left
        residual.set(w * 3 / 4 - 1, h / 2 - 1, Cplx.Float.one().mult(0.03f));   // bottom left

        residual.set(w / 2 - 1, h / 2 - 1, Cplx.Float.one());    // zero order

        // TODO: check if this is what matlabs conv2 would do
        Transforms.swapQuadrant(residual);
        Transforms.fft2d(residual, true);

        spec.times(residual);
        Transforms.fft2d(spec, false);

    }

    /**
     * Helper function, returns the magnitude of an input vector
     */
    Vec2d.Real magnitude(Vec2d.Cplx in) {
        Vec2d.Real ret = Vec2d.createReal(in);
        ret.copyMagnitude(in);
        return ret;
    }

    // ugly trick to allow for array of generic typhes
    private interface ListOfGratings extends List<Grating[]> {
    };

    
    public Grating[][] outputGratingSets(List<List<Grating[]>> otherWavelength, Grating[] currentGrating, double[] wavelength, int nr_dir, Vec2d.Real gaussProfile, double max_unwanted, int mask_size, DisplayWrapper imgSpatial, DisplayWrapper imgFourier, List<Grating[][]> resultList, String nameString) {
//        Tool.trace("otherWavelength.size() = " + otherWavelength.size());
        //                 String res = "1";
        //                 for (int ch = 1; ch < channels.length; ch++) {
        //                     res += " " + otherWavelength.get(ch - 1).size();
        //                 }
        
        //                 Tool.trace("found a full combination, adding it to final results:" + res);
        
        Grating[][] fullSet = new Grating[wavelength.length][];
        
        fullSet[0] = currentGrating;
        for (int ch = 1; ch < wavelength.length; ch++) {
            fullSet[ch] = lowestAvrEuclDist(currentGrating, otherWavelength.get(ch - 1));
        }
        
        // run the fourier check again, just so we can show the output
        for (int ch = 0; ch < wavelength.length; ch++) {
            fourierCheck(fullSet[ch], gaussProfile, max_unwanted, mask_size,
                         true, imgSpatial, imgFourier, String.format("set: %2d wl:%3.0f %s", resultList.size(), wavelength[ch], nameString ));
        }
        
        // resort (wavelength last, for output)
        Grating[][] fullSetResorted = new Grating[nr_dir][wavelength.length];
        for (int ch = 0; ch < wavelength.length; ch++) {
            for (int d = 0; d < nr_dir; d++) {
                fullSetResorted[d][ch] = fullSet[ch][d];
            }
        }
        return fullSetResorted;
    }
    
    
    boolean matchOtherWavelengths(List< List< Grating [] >> otherWavelength, double [] wavelength, int nr_dir, Grating [] currentGrating, double max_euclDist, int max_candidates, int countOk, Vec2d.Real gaussProfile, double max_unwanted, int mask_size, List<Grating[][]> resultList, List<List<Grating>> all, String iterationString) {
        boolean hasCandidatesForAll = true;
        // for each additional channel ...
        for (int ch = 1 ; ch < wavelength.length; ch++ ) {
            Tool.tell("matchOtherWavelengths " + iterationString + " "+ ch + "/" + (wavelength.length-1) + " " + " ok: "+ resultList.size());
            // ... create a list of candidates for that wavelength by...
            List< List<Grating> > directionCandidates = new ArrayList<List <Grating>>(nr_dir);
            
            // ... for every orientation: find all gratings that are in o.k.
            // eucl. distance to the main wavelength grating ...
            for (int d=0; d<nr_dir; d++) {
                Grating main = currentGrating[d];
                List<Grating> thisDir = new ArrayList<Grating>();
                for ( Grating cand : all.get(ch) ) {
                    //                        if(Math.abs(Grating.wrapPi(cand.gratDir) - Grating.wrapPi(main.gratDir)) > max_angle)
                    //                            continue;
                    if ( Grating.scaledDistance( cand, main ) > max_euclDist)
                        continue;
                    thisDir.add(cand);
                }
                directionCandidates.add(thisDir);
            }
            
            // ... check if there are combinations available for all directions ..
            boolean candidatesAvailableForAllDir = true;
            for (int d=0; d<nr_dir; d++) {
                if ( directionCandidates.get(d).size() == 0 )
                    candidatesAvailableForAllDir=false;
            }
            if (!candidatesAvailableForAllDir) {
                return false;
            }
            
            // ... create some random combinations of these ... 
            Grating [][] combinations = new Grating[max_candidates][nr_dir];
            for (int i=0; i<max_candidates; i++) {
                for (int d=0; d<nr_dir; d++) {
                    List<Grating> forThisDir = directionCandidates.get(d);
                    combinations[i][d] = forThisDir.get( (int)
                    (Math.random() * forThisDir.size()) );
                }
            }
            
            // ... and Fourier-check these combinations
            List< Grating []> secondaryList = new ArrayList< Grating [] >();
            countOk =0;
            for (int i=0; i<max_candidates; i++) {
                boolean ok = fourierCheck( combinations[i], gaussProfile, max_unwanted, mask_size, false, null, null, "i:"+i);
                if (ok) {
                    countOk++;
                    secondaryList.add( combinations[i] );
//                     Tool.trace("countOk = " + countOk);
                }
            }
            
            // if we did not find any candidates, set the bool to false
            if (countOk==0)
                return false;
            
            otherWavelength.add(secondaryList);
//            Tool.trace("otherWavelength.add(secondaryList) -> otherWavelength.size() = " + otherWavelength.size());
        }
        return true;
    }
            
    /**
     * Run the full calculation.
     *
     * @param wavelength list of wavelength to calculate for
     * @param gratMin	for each wavelength, minimal grating period, in #pxl
     * @param gratMax	for each wavelength, maximal grating period, in #pxl
     * @param phases	Number of phases
     * @param nr_dir	Number of pattern orientations
     * @param max_angle Maximal deviation from ideal (pi/n) angle distribution
     * of pattern direction, in rad
     * @param mask_size Size of holes in mask, in #pxl
     * @param max_unwanted Maximal contribution of unwanted orders
     * @param output_failed If to output also failed pattern
     * @param max_candidates How many candidates to calculate in first step
     * @return List of matching gratings
     *
     */
    public List<Grating[][]> calculate(
        double[] wavelength,
        boolean Tirf,
        double resImpAvr[], double NA,
        double[][] gratMin, double[][] gratMax,
        int phases, int nr_dir, double max_angle, int mask_size,
        double max_unwanted, double max_euclDist,
        boolean output_failed, final int max_candidates) {
        
        int countModOkAdded = 0;
        int tirf = Tirf ? 1 : 0;
        

//      resultList[#result][#angles][#wl] #angles doubled for tirf
        List<Grating[][]> resultList = new ArrayList<Grating[][]>();

        SimpleMT.useParallel(true);
        
        Tool.trace("-- Compute all candidates --");

        // for each wavelength, create a list of candidates
        List< List<Grating>> all = new ArrayList< List<Grating>>();
        List< List<Grating>> allTirf = new ArrayList< List<Grating>>();
        for (int ch = 0; ch < wavelength.length; ch++) {
            List< Grating> candidate = calcGrat(gratMin[0][ch], gratMax[0][ch], phases, wavelength[ch]);
            all.add(candidate);
            Tool.trace(String.format("NA %1.2f resImp %1.2f %5.0f nm : %d candidates", NA, resImpAvr[0], wavelength[ch], candidate.size()));
            if (Tirf) {
                List< Grating> candidateTirf = calcGrat(gratMin[1][ch], gratMax[1][ch], phases, wavelength[ch]);
                allTirf.add(candidateTirf);
                Tool.trace(String.format("NA %1.2f resImp %1.2f %5.0f nm : %d candidates", NA, resImpAvr[1], wavelength[ch], candidateTirf.size()));
            } else {
                allTirf = all;
            }
        }
        
        double [] channels = new double[wavelength.length*(1+1*tirf)];
        for (int ch=0; ch<wavelength.length; ch++) {
            channels[ch] = wavelength[ch];
            if(Tirf) {
                channels[ch+wavelength.length] = wavelength[ch];
                all.add(allTirf.get(ch));
            }
        }

        // for the first, main wavelength, get a list of matching gratings...
        Tool.trace("-- Compute direction combinations (main wavelength) --");
        List<Grating[]> dirs = selectDirs(all.get(0), allTirf.get(0), Tirf, nr_dir, max_angle);
        List<Grating[]>  dirsTirf = new ArrayList<Grating[]>();
        Tool.trace("   grating pairs for main wavelength: " + dirs.size());
        if(Tirf) {
            dirsTirf = selectDirs(allTirf.get(0),all.get(0), Tirf, nr_dir, max_angle);
            Tool.trace("   grating pairs for main wavelength (TIR): " + dirsTirf.size());
        }

        // run these pairs through fourier checking ...
        Tool.trace("-- Main wavelength: compute wanted vs. unwanted orders --");
        List<Grating[]> modOk = new ArrayList<Grating[]>();

        Vec2d.Real gaussProfile = createGaussIllum(fsPxl / 2.2, fsPxl);
        int countOk = 0;
        for (int i = 0; i < dirs.size(); i++) {

            Tool.tell("Main WL: FFT set " + i + "/" + dirs.size() + " ok: " + countOk + "/" + max_candidates);

            boolean ok = fourierCheck(dirs.get(i), gaussProfile, max_unwanted, mask_size, false, null, null, "i:" + i);

            if (ok) {
                countOk++;
                modOk.add(dirs.get(i));
            }
            if (countOk >= max_candidates) {
                break;
            }
        }

        Tool.trace(" main wavelength grating sets passing Fourier check: " + countOk);
        Tool.tell(" main wavelength grating sets passing Fourier check: " + countOk);

        DisplayWrapper imgSpatial = new DisplayWrapper(fsPxl, fsPxl, "Spatial");
        DisplayWrapper imgFourier = new DisplayWrapper(fsPxl, fsPxl, "Fourier");

        // now, try to find matching secondary and tertiary wavelength ...
        Tool.trace("-- find matching candidates at additional wavelengths or resolution settings--");

        int j=0;
        for ( Grating [] currentGrating : modOk ) {	    // loop each candidate
            j++;
            Tool.trace(j + "/" + modOk.size());
            String iterationString = ("addWL for mainWL " + j + "/" + modOk.size());
            //find patterns with matching angle and second pattern size
            List< List< Grating [] >> otherWavelength = new ArrayList< List< Grating []> >();
            List< List< Grating [] >> otherWavelengthTirf = new ArrayList< List< Grating []> >();
            boolean hasCandidatesForAll=true;
            boolean hasCandidatesForAllTirf=true;
            hasCandidatesForAll = matchOtherWavelengths(otherWavelength, wavelength, nr_dir, currentGrating, max_euclDist, max_candidates, countOk, gaussProfile, max_unwanted, mask_size, resultList, all, iterationString);
            if(!hasCandidatesForAll)
                continue;
            Tool.trace(j + "/" + modOk.size() + " has candidates for all angles"); 

                //create List of Gratings[nr_dirs] that match modOk[nr_dirs]
            List<Grating[]> modOkTirf = new ArrayList<Grating[]>();
            if(Tirf) {
                int countOkTirf = 0;
                for (int i = 0; i < dirsTirf.size(); i++) {
                    boolean matchDir=true;
                    for (int dir=0; dir < dirsTirf.get(i).length; dir++) {
                        if(Math.abs(Grating.wrapPi(currentGrating[dir].gratDir) - Grating.wrapPi(dirsTirf.get(i)[dir].gratDir) ) > max_angle)
                            matchDir = false;
                    }
                    if(matchDir == false)
                        continue;
                    boolean ok = fourierCheck(dirsTirf.get(i), gaussProfile, max_unwanted, mask_size, false, null, null, "i:" + i);
                    if (ok) {
                        countOkTirf++;
                    }
                        countOkTirf++;
                        modOkTirf.add(dirsTirf.get(i));
//                         if (countOkTirf >= max_candidates) {
//                             break;
//                         }
                }
                int k=0;
                int countResAdded=0;
                for (Grating [] currentGratingTirf : modOkTirf) {
                    k++;
//                     Tool.trace(j + "/" + modOk.size() + " " + k + "/" + modOkTirf.size());
                    iterationString = ("main WL "+j + "/" + modOk.size() + " addWL for mainTIRWL " + k + "/" + modOkTirf.size());
                    hasCandidatesForAllTirf = matchOtherWavelengths(otherWavelengthTirf, wavelength, nr_dir, currentGratingTirf, max_euclDist, max_candidates, countOkTirf, gaussProfile, max_unwanted, mask_size, resultList, allTirf, iterationString);
                   
                    // output some success
                    if (hasCandidatesForAll && hasCandidatesForAllTirf) {
//                         Tool.trace(j + "/" + modOk.size() + " has candidates for all angles, TIRF");
                        Grating[][] tmpResultPlain = outputGratingSets( otherWavelength, currentGrating, wavelength, nr_dir, gaussProfile, max_unwanted, mask_size, imgSpatial, imgFourier, resultList, "main");
                        Grating[][] tmpResultTirf = outputGratingSets( otherWavelengthTirf, currentGratingTirf, wavelength, nr_dir, gaussProfile, max_unwanted, mask_size, imgSpatial, imgFourier, resultList, "extra");
                        Grating[][] tmpResult = new Grating[nr_dir*2][Math.min(tmpResultPlain[0].length, tmpResultTirf[0].length)];
//                        Tool.trace("tmpResult "+tmpResult.length+ " "+tmpResult[0].length);
//                        Tool.trace("tmpResult "+tmpResultTirf.length+" "+tmpResultTirf[0].length);
                        for(int x=0; x<nr_dir; x++) {
                            for(int y=0; y<tmpResult[0].length; y++) {
                                tmpResult[x][y] = tmpResultPlain[x][y];
                            }
                            for(int y=0; y<tmpResult[0].length; y++) {
                                tmpResult[x+nr_dir][y] = tmpResultTirf[x][y];
                            }
                        }
//                        Tool.trace("tmpResult "+tmpResult.length+ " "+tmpResult[0].length);
                        resultList.add(tmpResult);
                        countResAdded++;
                        Tool.trace(j + "/" + modOk.size() + " has candidates for all angles, TIR. Adding result " +countResAdded+ ". Total: "+resultList.size()+"/"+max_candidates);
//                         if(countResAdded >= max_candidates) break;
//                         if(countModOkAdded >= 1) break;
                    } else {
//                         Tool.trace("no luck");
                    }
                }
                
            } else {
                if (hasCandidatesForAll) {
                        resultList.add(outputGratingSets(otherWavelength, currentGrating, wavelength,  nr_dir, gaussProfile, max_unwanted, mask_size, imgSpatial, imgFourier, resultList, ""));
                    }
            }
            if(resultList.size()>=max_candidates)
                break;
        }

        imgSpatial.display();
        imgFourier.display();

        Tool.trace("    Number of full sets of pattern " + resultList.size());

        /*
	// check if unwanted orders can be blocked by masking
	Tool.trace("-- Compute wanted vs. unwanted orders --");
	
	Vec2d.Real gaussProfile = createGaussIllum( fsPxl/2.2, fsPxl);
	int countOk=0;
	for (int i=0; i< dirs.size() ; i++) {
	    
	    if (i%10==0) Tool.tell(" FFT "+i+"/"+dirs.size());
	    
	    boolean ok = true;
	    
	    for (int ch=0; ch<wavelength.length; ch++) {
		ok &= fourierCheck( dirs.get(i)[ch], gaussProfile, max_unwanted, mask_size, 
		    output_failed, imgSpatial, imgFourier, "i:"+i); 
	    }
	    
	    if (ok) {
		countOk++;
		resultList.add( dirs.get(i) );
	    }
	}
	
	Tool.trace("    Number of pattern after modulation check: "+ countOk);
         */
        return resultList;
    }

    // ---- Start methods ----
    /**
     * ImageJ plugin run method
     */
    @Override
    public void run(String arg) {

        // redirect log output to FIJIs log
        Tool.setLogger(new Tool.Logger() {
            @Override
            public void writeTrace(String w) {
                ij.IJ.log(w);
            }

            @Override
            public void writeShortMessage(String w) {
                ij.IJ.showStatus(w);
            }

            @Override
            public void writeError(String w, boolean err) {
                ij.IJ.log("ERR: " + err);
            }

        });

        // generate simple GUI
        GenericDialog gd = new GenericDialog("SLM pattern search");
        gd.addMessage("System parameters");
        String[] items = {"SXGA-3DM", "DLP6500FYE", "other"};
        gd.addRadioButtonGroup("SLM", items, 3, 1, "SXGA-3DM");
        gd.addStringField("    other", "SLM", 30);
        gd.addNumericField("    other: pixel size", 13.62, 2, 6, "µm");
        gd.addNumericField("    other: pixels X  ", 1280, 0);
        gd.addNumericField("    other: pixels Y  ", 1024, 0);
        gd.addMessage("              ");
        gd.addNumericField("                SLM scale factor", (300), 2);
        gd.addNumericField("               objective lens NA", 1.47, 2);
        gd.addNumericField("     refractive index n (sample)", 1.36, 2);
        gd.addNumericField("resolution enhancement average  ", 1.80, 3);
        gd.addNumericField("resolution enhancement range +- ", 0.005, 3);

        gd.addCheckbox(    "create additional TIR pattern", true);
//         gd.addNumericField("        objective lens NA or sample n", 1.36, 2);
//         gd.addNumericField("       resolution enhancement average", 2.00, 4);
//         gd.addNumericField("      resolution enhancement range +-", 0.002, 4);

        gd.addMessage("Wavelength to analyse");
        gd.addNumericField("main wavelength", 488, 0, 6, "nm");
        gd.addCheckbox("use additional wavelength 1", true);
        gd.addNumericField(" add wavelength 1", 568, 0, 6, "nm");
        gd.addCheckbox("use additional wavelength 2", true);
        gd.addNumericField("add wavelength 2", 647, 0, 6, "nm");

        gd.addMessage("Pattern parameters");
        gd.addNumericField("#pattern_directions", 3, 0);
        gd.addNumericField("max_deviation_ideal_angle", 2, 1, 6, "deg");
        gd.addNumericField("#phases", 3, 0);
        gd.addNumericField("max_eucl_dist(approx. pxl)", .05, 2);
        gd.addMessage("Modulation");
        gd.addNumericField("max_unwanted_modulation", 0.015, 3);
        gd.addNumericField("mask_size", 15, 0, 6, "pxl");
        gd.addCheckbox("Output_also_failed", false);
        gd.addMessage("Cancel");
        gd.addNumericField("max_nr_candidates", 200, 0);

        gd.showDialog();
        if (gd.wasCanceled()) {
            return;
        }

        // ---- get parameters ----
        final double pxlSize, slmPxlSize, slmScale;
        final int slmPxlX, slmPxlY;
        final String slmType = gd.getNextRadioButton(), prefixSlm;

        if (slmType.equals("SXGA-3DM")) {
            prefixSlm = "SXGA3DM";
            gd.getNextString();
            slmPxlSize = 13.62;
            gd.getNextNumber();
            slmPxlX = 1280;
            gd.getNextNumber();
            slmPxlY = 1024;
            gd.getNextNumber();
            slmScale = gd.getNextNumber();
            pxlSize = 1000. * slmPxlSize / slmScale;
        } else if (slmType.equals("DLP6500FYE")) {
            prefixSlm = "DLP6500";
            gd.getNextString();
            slmPxlSize = 7.56;
            gd.getNextNumber();
            slmPxlX = 1920;
            gd.getNextNumber();
            slmPxlY = 1080;
            gd.getNextNumber();
            slmScale = gd.getNextNumber();
            pxlSize = 1000. * slmPxlSize / slmScale;
        } else {
            prefixSlm = gd.getNextString();
            slmPxlSize = gd.getNextNumber();
            slmPxlX = (int) gd.getNextNumber();
            slmPxlY = (int) gd.getNextNumber();
            slmScale = gd.getNextNumber();
            pxlSize = 1000. * slmPxlSize / slmScale;
        }
        double refInd;
        double objNA;
        double[] resImpAvr = new double[2];
        double[] resImpRange = new double[2];
        boolean Tirf;
        double[] wavelength_gui = new double[3];
        boolean[] wavelength_gui_switch = new boolean[3];
        double[] wavelength;

        objNA = gd.getNextNumber();
        refInd = gd.getNextNumber();
        if(objNA < refInd) {
            refInd=objNA;
        }
        resImpAvr[0] = gd.getNextNumber();
        resImpRange[0] = gd.getNextNumber();
        Tirf = gd.getNextBoolean();
        int tirf = Tirf ? 1 : 0;
        resImpAvr[1] = 1.5+(objNA)/refInd/2.;
        resImpRange[1] =resImpRange[0]/2;

        // TODO: there must be a nicer way to code this
        {
            // copy all wavelength and switches
            int count_active = 1;
            for (int i = 0; i < 3; i++) {
                wavelength_gui[i] = gd.getNextNumber();
                if (i != 0) {
                    wavelength_gui_switch[i] = gd.getNextBoolean();
                    if (wavelength_gui_switch[i] == true) {
                        count_active++;
                    }
                }
            }

            // create the final array 
            wavelength = new double[count_active];

            wavelength[0] = wavelength_gui[0];
            int count_pos = 1;
            for (int i = 1; i < 3; i++) {
                if (wavelength_gui_switch[i] == true) {
                    wavelength[count_pos++] = wavelength_gui[i];
                }
            }
        }

        final int nrDirs = (int) gd.getNextNumber();
        final double maxAngleDev =  gd.getNextNumber() * Math.PI/180.;
        final int nrPhases = (int) gd.getNextNumber();
        final double maxEuclDist = gd.getNextNumber();

        final double maxUnwMod = gd.getNextNumber();
        final int maskSize = (int) gd.getNextNumber();
        final boolean outputFailed = gd.getNextBoolean();
        final int maxCandidates = (int) gd.getNextNumber();
        final String prefix = String.format("%s_%da%dp", prefixSlm, nrDirs, nrPhases);

        // 1 - calculate the ranges for all from resolution enhancement
        //
        // This is pttr/2 = lambda / ( 2 * NA * (resImp-1) * pxlSize)
        // where
        //  pttr    : pattern spacing, in pxls
        //  lambda  : exitation wavelength
        //  NA	    : objective's NA
        //  resImp  : factor of resolution improvement, e.g. 1.75x
        //  pxlSize : projected pixel size 
        //
        double resImpMin[] = new double[2];
        double resImpMax[] = new double[2];
        resImpMin[0] = resImpAvr[0] - resImpRange[0];
        resImpMin[1] = resImpAvr[1] - resImpRange[1];
        resImpMax[0] = resImpAvr[0] + resImpRange[0];
        resImpMax[1] = resImpAvr[1] + resImpRange[1];


        IJ.log(String.format("Searching pattern for res. improvement %7.4f -- %7.4f",
            resImpMin[0], resImpMax[0]));
        IJ.log(String.format("Searching pattern for res. improvement %7.4f -- %7.4f",
            resImpMin[1], resImpMax[1]));

        IJ.log(String.format("Objective %7.4f NA, proj. SLM pixels %5.1f nm", refInd, pxlSize));

        double[][] gratMin = new double[2][3];
        double[][] gratMax = new double[2][3];

        for (int ch = 0; ch < wavelength.length; ch++) {
            gratMin[0][ch] = 2 * wavelength[ch] / (2 * refInd * (resImpMax[0] - 1) * pxlSize);
            gratMin[1][ch] = 2 * wavelength[ch] / (2 * refInd * (resImpMax[1] - 1) * pxlSize);
            gratMax[0][ch] = 2 * wavelength[ch] / (2 * refInd * (resImpMin[0] - 1) * pxlSize);
            gratMax[1][ch] = 2 * wavelength[ch] / (2 * refInd * (resImpMin[1] - 1) * pxlSize);

            IJ.log(String.format("Search range %5.0f nm: %7.4f --- %7.4f pxl", wavelength[ch],
                gratMin[0][ch], gratMax[0][ch]));
                IJ.log(String.format("Search range %5.0f nm: %7.4f --- %7.4f pxl", wavelength[ch],
                gratMin[1][ch], gratMax[1][ch]));

        }

        // run the actual calculation
        List<Grating[][]> res = calculate(wavelength, Tirf,
            resImpAvr, refInd,
            gratMin, gratMax, nrPhases,
            nrDirs, maxAngleDev, maskSize,
            maxUnwMod, maxEuclDist, outputFailed, maxCandidates);

        // store the result, if any
        if (res.size() > 0) {
            ij.IJ.setProperty("de.bio_photonics.gratingsearch.phaseNumber", nrPhases);
            ij.IJ.setProperty("de.bio_photonics.gratingsearch.angNumber", nrDirs );
            ij.IJ.setProperty("de.bio_photonics.gratingsearch.resImp", resImpAvr);
            ij.IJ.setProperty("de.bio_photonics.gratingsearch.wl", wavelength);
            ij.IJ.setProperty("de.bio_photonics.gratingsearch.na", refInd);
            ij.IJ.setProperty("de.bio_photonics.gratingsearch.lastGratings", res);
            ij.IJ.setProperty("de.bio_photonics.gratingsearch.width", slmPxlX);
            ij.IJ.setProperty("de.bio_photonics.gratingsearch.height", slmPxlY);
            ij.IJ.setProperty("de.bio_photonics.gratingsearch.prefix", prefix);
            ij.IJ.setProperty("de.bio_photonics.gratingsearch.resImp", resImpAvr);
        }

        // find the lowest modulation
        double minTransmission = Double.MAX_VALUE;
        int minPos = 0;
        
//      [#][maxUnwantedMod|biggestAngleDevBetween(Non)Tirf)|AvgRes|spread|AngTirfRes|spead]
        double resTab[][] = new double[res.size()][7];
        
        for (int i = 0; i < res.size(); i++) { //for all results
            // get the highest unwanted modulation from the current grating set
            double maxHere = Double.MIN_VALUE;
            double avgRes[] = {0,0};
            double minRes[] = {Double.MAX_VALUE, Double.MAX_VALUE};
            double maxRes[] = {Double.MIN_VALUE, Double.MIN_VALUE};
            // res[#result][#angles][#wl] #angles doubled for tirf
            Grating[][] grs = res.get(i); //pick result
            for (int dir=0; dir<nrDirs*(1+tirf); dir++) {    //for all angles
                Grating[] gr = grs[dir];
                for (int ch=0; ch<gr.length; ch++) {    //for all wavelengths
                    Grating g = gr[ch];
                    if (g.unwantedMod > maxHere) {
                        maxHere = g.unwantedMod;
                    }
                    double curRes = 1 + wavelength[ch] / g.gratPer / refInd / pxlSize;
                    avgRes[dir/nrDirs] += curRes;
                    if(curRes<minRes[dir/nrDirs]) {
                        minRes[dir/nrDirs] = curRes;
                    }
                    if(curRes>maxRes[dir/nrDirs]) {
                        maxRes[dir/nrDirs] = curRes;
                    }
                }
            }
            avgRes[0] = avgRes[0]/(nrDirs*wavelength.length);
            avgRes[1] = avgRes[1]/(nrDirs*wavelength.length);
            resTab[i][0] = maxHere;
            resTab[i][1] = 0;
            resTab[i][2] = avgRes[0];
            resTab[i][3] = maxRes[0]-minRes[0];
            resTab[i][4] = avgRes[1];
            resTab[i][5] = maxRes[1]-minRes[1];
            resTab[i][6] = i;
            
            // see if this is the minimum in the list
            if (maxHere < minTransmission) {
                minTransmission = maxHere;
                minPos = i;
            }

            IJ.log("Pattern set " + i + " highest unwanted modulation: " + maxHere);
        }

        java.util.Arrays.sort(resTab, new java.util.Comparator<double[]>() {
            public int compare(double[] a, double[] b) {
                return Double.compare(a[0], b[0]);
            }
        });
        
        // output position of lowest unwanted modulation 
        IJ.log("Lowest unwanted modulation at pos " + minPos);
        for (int r=0; r<resTab.length; r++) {
            IJ.log((int)resTab[r][6]+" "+resTab[r][0]+" "+resTab[r][1]+" "+resTab[r][2]+" "+resTab[r][3]+" "+resTab[r][4]+" "+resTab[r][5] );
        }
    }

    /**
     * main method
     */
    public static void main(String[] args) {

        Grating_Search gs = new Grating_Search();
        new ij.ImageJ(ij.ImageJ.EMBEDDED);

        gs.run("");

        /*
	double gratMin = Double.parseDouble( args[0] );
	double gratMax = Double.parseDouble( args[1] );
	int nrPhases   = Integer.parseInt( args[2] );
	int nrDirs     = Integer.parseInt( args[3] );

	gs.calculate(gratMin, gratMax, nrPhases, 
	    nrDirs, 3./180*Math.PI, 20, 0.02, false, 400); */
    }

}
