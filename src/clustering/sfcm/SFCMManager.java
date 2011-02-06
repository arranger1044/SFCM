/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package clustering.sfcm;

import ij.IJ;
import ij.process.ByteProcessor;
import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.Duplicator;
import ij.process.ColorProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageConverter;
import ij.process.StackConverter;
import java.io.IOException;
import java.io.RandomAccessFile;
import util.ColorManager;
import util.ColorSpaceConversion;
import vectorLib.VectorProcessor;

/**
 * This class isolates the Spatial Fuzzy C-Means Algorithm (SFC) making it
 * independent from using the ImageJ API. It takes the responsibility to run it
 * and to communicate the results to the SFCMPlugin class.
 * It also provides a way to encode the input image in a different color space
 * and to decode the clustered image in several ways.
 * The class conforms to the @link{ClusteringDelegate} interface, implementing the
 * delegation pattern and allowing the delegating algorithm to communicate with
 * the plugin.
 * Optionally it can be instantiated by the @link{TestManager} to execute some tests
 * @see SFCM
 * @see ClusteringDelegate
 * @see SFCMPlugin
 * @see TestManager
 */
public class SFCMManager implements ClusteringDelegate{

   /**
    * Boolean used to check if a random sequence has to be seeded
    * @see #randomizationSeed
    */
   private boolean seedEnabled = true;
   /**
    * The integer used to generate a random sequence
    * Default value is 48
    * @see Random
    */
   private int randomizationSeed = 48;
   /**
    * The iteration error upper limit
    * Default value is 0.0001
    */
   private double tolerance = 0.0001;
   /**
    * The number of cluster in which we want to partition the image
    * Default value is 4
    */
   private int numberOfClusters = 4;
   /**
    * A boolean indicating whether the clustered image shall be visualized in
    * gray scale.
    * Default value is true
    * @see #encodeClusteredImageInGray(vectorLib.VectorProcessor, int[][], float[][])
    */
   private boolean grayScaleVisualization = true;
   /**
    * A boolean indicating whether the clustered image shall be visualized
    * using random RGB colors.
    * Default value is false
    * @see #encodeClusteredImageInRGB(vectorLib.VectorProcessor, int[][], float[][], int[], int)
    */
   private boolean randomRGBVisualization = false;
   /**
    * A boolean indicating whether the clustered image shall be visualized
    * using an @link{ImageStack} of binary slices.
    * Default value is false
    * @see #encodeClusteredImageBinaryStack(vectorLib.VectorProcessor, int[][], float[][])
    */
   private boolean binaryStackVisualization = false;
   /**
    * A boolean indicating whether the clustered image shall be visualized
    * using an RGB image that allows only the cluster center colors.
    * Default value is false
    * @see #encodeClusteredImageWithClusterCenterColors(vectorLib.VectorProcessor, int[][], float[][]) 
    */
   private boolean clusterCenterColorsVisualization = false;
   /**
    * An array of strings listing the possible initialization modes for the cluster
    * center matrix and the cluster membership matrix;
    * Items are : "<b>random V (Forgy)</b>", "<b>K-Means++</b>", "<b>random U</b>"
    * @see #initializationMode
    */
   private static final String[] initModes = {"random V (Forgy)", "K-Means++", "random U"};
   /**
    * A string containing the current initialization mode string.
    * Default value is <b>K-Means++</b>
    * @see #initModes
    */
   private String initializationMode = "K-Means++";
   /**
    * A string array listing the valid color conversion space strings.
    * Values are: "<b>None</b>", "<b>XYZ</b>", "<b>L*a*b*</b>", "<b>HSB</b>"
    * @see #colorSpace
    */
   private static final String[] colorSpaces = {"None", "XYZ", "L*a*b*", "HSB"};
   /**
    * A string containing the current color conversion mode
    * Default value is "<b>None</b>"
    * @see #colorSpaces
    */
   private String colorSpace = "None";
   /**
    * An integer holding the image type of the currently opened image, according
    * to ImageJ types.
    * @see #run;
    */
   private int imageType;
   /**
    * A boolean telling whether the algorithm status updates are to be printed on
    * the system console.
    * Default value is false
    */
   private boolean printOnConsole = true;
   /**
    * The fuzzyness parameter (<code>m</code>) as a double
    * Default value is 2.0
    */
   private double fuzzyness = 2.0f;
   /**
    * A boolean indicating whether the clustered image shall be visualized
    * using an @link{ImageStack} whose slices represents the fuzzy memberships of
    * a pixel to a specified cluster
    * Default value is false
    * @see #encodeClusteredImageFuzzyStack(vectorLib.VectorProcessor, float[][], float[][])
    */
   private boolean fuzzyStackVisualization = false;
   /**
    * The maximum number of allowed iteration for the SFCM algorithm
    * @see SFCM
    */
   private long maxIterations = 400;
   /**
    * An array of strings listing the allowed stopping criterions for the algorithm
    * Listed values are: "<b>Frobenius Norm on U</b>", "<b>Frobenius Norm on V</b>",
    * "<b>Max Norm on U</b>", "<b>Max Norm on V</b>"
    * @see SFCM
    */
   private static final String[] stopCriterions = { "Frobenius Norm on U",
                                                   "Frobenius Norm on V", "Max Norm on U", "Max Norm on V"};
   /**
    * THe string containing the currently chosen value for the stop criterion.
    * Default value is "<b>Max Norm on U</b>"
    */
   private String stopCriterion = "Max Norm on U";
   /**
    * The radius for the window representing the neighborhood for a image pixel
    * Default value is 2;
    */
   private int windowRadius = 2;
   /**
    * A double used to weight the cluster membership function in the SFCM algorithm.
    * Default value is 1.0
    * @see SFCM
    */
   private double membershipWeight = 1.0;
   /**
    * A double used to weight the spatial function in the SFCM algorithm.
    * Default value is 0.0
    * @see SFCM
    */
   private double spatialFunctionWeight = 0.0;
   /**
    * A string array containing the allowed spatial function strings that identify
    * the kind of the spatial function used in the SFCM algorithm
    * Listed values are: "<b>Likeliest Cluster</b>", "<b>Weightiest Cluster</b>"
    * @see SFCM
    * @see #spatialFunction
    */
   private static final String[] spatialFunctions = {"Likeliest Cluster", "Weightiest Cluster"};
   /**
    * A string containing the currently selected spatial function kind
    * Default value is "<b>Likeliest Cluster</b>"
    * @see #spatialFunctions
    */
   private String spatialFunction = "Likeliest Cluster";
   /**
    * A reference to the file to which, while in testing mode, the algorithm status
    * updates are redirected. If <code>null</code> there is no need to write on file
    * Default value is <code>null</code>
    */
   private RandomAccessFile RAF = null;
   /**
    * A boolean indicating if the class has been instantiated in order to run a
    * batch of tests.
    * Default value is false
    * @see TestManager
    */
   private boolean testing = false;
   /**
    * A boolean indicating if the class has been instantiated in order to run a
    * Validation test.
    * Default value is false
    * @see TestManager
    * @see ClusteringValidity
    */
   private boolean validation = false;

   /**
    * Public getter for <code>seedEnabled</code>
    * @return the boolean value set
    * @see #seedEnabled
    */
   public boolean getSeedEnabled(){
       return seedEnabled;
   }

   /**
    * Public setter for <code>seedEnabled</code>
    * @param enabled the new boolean value to set
    * @see #seedEnabled
    */
   public void setSeedEnabled(boolean enabled){
       seedEnabled = enabled;
   }

   /**
    * Public getter for <code>randomizationSeed</code>
    * @return the currently used value
    * @see #randomizationSeed
    */
   public int getRandomizationSeed() {
        return randomizationSeed;
   }

   /**
    * Public setter for <code>randomizationSeed</code>
    * @param randomizationSeed a new randomization seed
    * @see #randomizationSeed
    */
   public void setRandomizationSeed(final int randomizationSeed) {
        this.randomizationSeed = randomizationSeed;
   }

   /**
    * Public getter for <code>numberOfClusters</code>
    * @return the currently used value
    * @see #numberOfClusters
    */
   public int getNumberOfClusters() {
       return numberOfClusters;
   }

   /**
    * Public setter for <code>numberOfClusters</code>
    * @param numberOfClusters a new number of clusters
    * @see #numberOfClusters
    */
   public void setNumberOfClusters(final int numberOfClusters) {
       this.numberOfClusters = numberOfClusters;
   }

   /**
    * Public getter for <code>tolerance</code>
    * @return the currently used value
    * @see #tolerance
    */
   public double getTolerance() {
       return tolerance;
   }

   /**
    * Public setter for <code>tolerance</code>
    * @param tolerance a new tolerance value
    * @see #tolerance
    */
   public void setTolerance(final float tolerance) {
       this.tolerance = tolerance;
   }

   /**
    * Public getter for <code>grayScaleVisualization</code>
    * @return the currently used boolean value
    * @see #grayScaleVisualization
    */
   public boolean getGrayScaleVisualization(){
       return grayScaleVisualization;
   }

   /**
    * Public setter for <code>grayScaleVisualization</code>
    * @param visualization the new boolean value
    * @see #grayScaleVisualization
    */
   public void setGrayScaleVisualization(boolean visualization){
       grayScaleVisualization = visualization;
   }

   /**
    * Public getter for <code>randomRGBVisualization</code>
    * @return the currently used boolean value
    * @see #randomRGBVisualization
    */
   public boolean getRandomRGBVisualization(){
       return randomRGBVisualization;
   }

   /**
    * Public setter for <code>randomRGBVisualization</code>
    * @param visualization the new boolean value
    * @see #randomRGBVisualization
    */
   public void setRandomRGBVisualization(boolean visualization){
       randomRGBVisualization = visualization;
   }

   /**
    * Public getter for <code>binaryStackVisualization</code>
    * @return the currently used boolean value
    * @see #binaryStackVisualization
    */
   public boolean getBinaryStackVisualization(){
       return binaryStackVisualization;
   }

   /**
    * Public setter for <code>binaryStackVisualization</code>
    * @param visualization the new boolean value
    * @see #binaryStackVisualization
    */
   public void setBinaryStackVisualization(boolean visualization){
       binaryStackVisualization = visualization;
   }

   /**
    * Public getter for <code>clusterCenterColorsVisualization</code>
    * @return the currently used boolean value
    * @see #clusterCenterColorsVisualization
    */
   public boolean getClusterCenterColorsVisualization(){
       return clusterCenterColorsVisualization;
   }

   /**
    * Public setter for <code>clusterCenterColorsVisualization</code>
    * @param visualization the new boolean value
    * @see #clusterCenterColorsVisualization
    */
   public void setClusterCenterColorsVisualization(boolean visualization){
       clusterCenterColorsVisualization = visualization;
   }

   /**
    * Public getter for <code>initializationMode</code>
    * @return the currently used initialization mode string
    * @see #initializationMode
    */
   public String getInizializationMode(){
       return initializationMode;
   }

   /**
    * Public setter for <code>initializationMode</code>
    * @param initMode a new initialization mode string
    * @see #initializationMode
    */
   public void setInitializationMode(String initMode){
       initializationMode = initMode;
   }

   /**
    * Public getter for <code>initModes</code>
    * @return the string array for the initialization modes
    * @see #initModes
    */
   public String[] getInitModes(){
       return initModes;
   }

   /**
    * Public getter for <code>colorSpaces</code>
    * @return the string array listing the color spaces
    * @see #colorSpaces
    */
   public String[] getColorSpaces(){
       return colorSpaces;
   }

   /**
    * Public getter for <code>colorSpace</code>
    * @return the currently used color space string
    * @see #colorSpace
    */
   public String getColorSpace(){
       return colorSpace;
   }

   /**
    * Public setter for <code>colorSpace</code
    * @param space the new color space string
    * @see #colorSpace
    */
   public void setColorSpace(String space){
       colorSpace = space;
   }

   /**
    * Public getter for <code>fuzzyness</code>
    * @return the currently used double value
    * @see #fuzzyness
    */
   public double getFuzzyness(){
       return fuzzyness;
   }

   /**
    * Public setter for <code>fuzzyness</code>
    * @param m the new double value
    * @see #fuzzyness
    */
   public void setFuzzyness(double m){
       fuzzyness = m;
   }

   /**
    * Public setter for the <code>RAF</code>
    * @param testFile the new file reference
    * @see #RAF
    */
   public void setTestFilePointer(RandomAccessFile testFile){
       RAF = testFile;
   }

   /**
    * Public getter for <code>fuzzyStackVisualization</code>
    * @return the currently used boolean value
    * @see #fuzzyStackVisualization
    */
   public boolean getFuzzyStackVisualization(){
       return fuzzyStackVisualization;
   }

   /**
    * Public setter for <code>fuzzyStackVisualization</code>
    * @param visualization the new boolean value
    * @see #fuzzyStackVisualization
    */
   public void setFuzzyStackVisualization(boolean visualization){
       fuzzyStackVisualization = visualization;
   }

   /**
    * Validates <code>fuzzyness</code>. Allowed values in (1.0, +inf)
    * @return a boolean value indicating whether the current value is valid or not
    * @see #fuzzyness
    */
   public boolean validateFuzzyness(){
       boolean validFuzzyness = (fuzzyness > 1f);

       return validFuzzyness;
   }

   /**
    * Validates <code>numberOfClusters</code>. Allowed values in [2, +inf)
    * @return a boolean value indicating whether the current value is valid or not
    * @see #numberOfClusters
    */
   public boolean validateClusterNumber(){
       return  numberOfClusters > 1;
   }

   /**
    * Validates <code>tolerance</code>. Allowed values in (0, +inf)
    * @return a boolean value indicating whether the current value is valid or not
    * @see #tolerance
    */
   public boolean validateTolerance(){
       return tolerance > 0;
   }

   /**
    * Public getter for <code>maxIterations</code>
    * @return the currently used long value
    * @see #maxIterations
    */
   public long getMaxIterations(){
       return maxIterations;
   }

   /**
    * Public setter for <code>maxIterations</code>
    * @param iterations the new max iterations value
    * @see #maxIterations
    */
   public void setMaxIterations(long iterations){
       maxIterations = iterations;
   }

   /**
    * Public getter for <code>stopCriterion</code>
    * @return the currently used string
    * @see #stopCriterion
    */
   public String getStopCriterion(){
       return stopCriterion;
   }

   /**
    * Public getter for <code>stopCriterions</code>
    * @return the allowed strings as stop criterions
    * @see #stopCriterions
    */
   public String[] getStopCriterions(){
       return stopCriterions;
   }

   /**
    * Public setter for <code>stopCriterion</code>
    * @param criterion the new stop criterion string
    * @see #stopCriterion
    */
   public void setStopCriterion(String criterion){
       stopCriterion = criterion;
   }

   /**
    * Public getter for <code>testing</code>
    * @return the currently used boolean value
    * @see #testing
    */
   public boolean isTesting(){
       return testing;
   }

   /**
    * Public setter for <code>testing</code>
    * @param test the new boolean value
    * @see #testing
    */
   public void setTesting(boolean test){
       testing = test;
   }

   /**
    * Public getter for <code>windowRadius</code>
    * @return the currently used integer value
    * @see #windowRadius
    */
   public int getWindowRadius(){
       return windowRadius;
   }

   /**
    * Public setter for <code>windowRadius</code>
    * @param radius the new integer value
    * @see #windowRadius
    */
   public void setWindowRadius(int radius){
       windowRadius = radius;
   }

   /**
    * Validates <code>windowRadius</code>. Allowed values in [0, +inf)
    * @return a boolean value indicating whether the current value is valid or not
    * @see #windowRadius
    */
   public boolean validateWindowsRadius(){
       return windowRadius >= 0;
   }

   /**
    * Public getter for <code>membershipWeight</code>
    * @return the currently used double value
    * @see #membershipWeight
    */
   public double getMembershipWeight(){
       return membershipWeight;
   }

   /**
    * Public setter for <code>membershipWeight</code>
    * @param weight the new membership weight value
    * @see #membershipWeight
    */
   public void setMembershipWeight(double weight){
       membershipWeight = weight;
   }

   /**
    * Validates <code>membershipWeight</code>. Allowed values in [0.0, +inf)
    * @return a boolean value indicating whether the current value is valid or not
    * @see #membershipWeight
    */
   public boolean validateMembershipWeight(){
       return membershipWeight >= 0.0;
   }

   /**
    * Public getter for <code>spatialFunctionWeight</code>
    * @return the currently used value for the spatial function weight
    * @see #spatialFunctionWeight
    */
   public double getSpatialFunctionWeight(){
       return spatialFunctionWeight;
   }

   /**
    * Public setter for <code>spatialFunctionWeight</code>
    * @param weight a new double for the spatial function weight
    * @see #spatialFunctionWeight
    */
   public void setSpatialFunctionWeight(double weight){
       spatialFunctionWeight = weight;
   }

   /**
    * Validates <code>spatialFunctionWeight</code>. Allowed values in [0.0, +inf)
    * @return a boolean value indicating whether the current value is valid or not
    * @see #spatialFunctionWeight
    */
   public boolean validateSpatialFunctionWeight(){
       return spatialFunctionWeight >= 0.0;
   }

   /**
    * Public getter for <code>spatialFunctions</code>
    * @return the string array listing the possible spatial functions kinds
    * @see #spatialFunctions
    */
   public String[] getSpatialFunctions(){
       return spatialFunctions;
   }

   /**
    * Public getter for <code>spatialFunction</code>
    * @return the currently used spatial function string
    * @see #spatialFunction
    */
   public String getSpatialFunction(){
       return spatialFunction;
   }

   /**
    * Public setter for <code>spatialFunction</code>
    * @param function the new spatial function value
    * @see #spatialFunction
    */
   public void setSpatialFunction(String function){
       spatialFunction = function;
   }

   /**
    * Public setter for <code>validation</code>
    * @param a boolean value, true in the case we should run a validation test
    * @see #validation
    */
   public void setValidation(boolean validity){
        validation = validity;
   }

   /**
    * The main method of the class. It encapsulates the @link{SFCM} <b>run</b>.
    * It takes an @link{ImagePlus} and converts it into a @link{VectorProcessor},
    * provides the conversion into another color space and extracts the pixel in
    * a data matrix form in order to run the SFCM algorithm. Once the algorithm
    * terminated, it encodes new images based on the clustering partition obtained
    * and the visualization mode specified.
    * @param img the currently opened image in ImageJ
    * @return and array of @link{ImageStack} containing the encoded clustered images
    * @see SFCM
    * @see SFCMPlugin
    * @see VectorProcessor
    * @see #computeInitializationMode()
    * @see #computeSpatialFunctionMode()
    * @see #computeStopCriterion()
    * @see #convertMatrixColorSpace(float[][], java.lang.String)
    * @see #convertToColorSpace(ij.ImagePlus, java.lang.String)
    */
   public ImageStack[] run(ImagePlus img)
    {
        ImageStack[] imgArray = null;
        imageType = img.getType();

        VectorProcessor vp = convertToColorSpace(img, colorSpace);

        float [][] imageData = vp.getPixels();

        int initMode = computeInitializationMode();
        int stopCrit = computeStopCriterion();
        int sFunction = computeSpatialFunctionMode();

        /* Calling the clustering algorithm */
        final long startTime = System.currentTimeMillis();
        Object[] resultMatrixes = SFCM.run(imageData, numberOfClusters, tolerance,
                                             seedEnabled,
                                             randomizationSeed, initMode, this,
                                             fuzzyness, maxIterations, stopCrit,
                                             windowRadius, membershipWeight,
                                             spatialFunctionWeight, sFunction,
                                             vp.getWidth(), testing);
        final long endTime = System.currentTimeMillis();
        final String message = "Clustering completed in " + (endTime - startTime) + " ms.";
        updateStatus(message);

        int [][] clusterMemberships = (int[][]) resultMatrixes[0];
        float [][] clusterCenters = (float[][]) resultMatrixes[1];

        /* If a validation test is being run we must compute some validity indexes*/
        if (validation) {
            float[][] membershipMatrix = (float[][]) resultMatrixes[2];
            double vpc = ClusteringValidity.bezdekPartitionCoefficient(membershipMatrix);
            double vpe = ClusteringValidity.partitionEntropyIndex(membershipMatrix);
            double vxb = ClusteringValidity.compactnessAndSeparationMetric(imageData,
                    membershipMatrix,
                    clusterCenters,
                    fuzzyness);
            try {
                testingValidation(fuzzyness, membershipWeight, spatialFunctionWeight, 
                        windowRadius, initializationMode, spatialFunction,
                        numberOfClusters, colorSpace, vpc, vpe, vxb);
            } catch (IOException ex) {
                ex.printStackTrace();
            }
        }

        Object[] stacks = new Object[5];
        
        int visualizationModes = 0;

        /* Computing the encodings according the visualizaton modes selected by
           the user*/
        if (grayScaleVisualization)
        {
            ImageStack imgsbp = new ImageStack(vp.getWidth(), vp.getHeight());
            imgsbp.addSlice("GrayScale", encodeClusteredImageInGray(vp, clusterMemberships,
                                                                    clusterCenters));
            stacks[visualizationModes] = imgsbp;
            visualizationModes++;
        }

        if (randomRGBVisualization)
        {
            ImageStack imgsbp2 = new ImageStack(vp.getWidth(), vp.getHeight());
            imgsbp2.addSlice("RGB", encodeClusteredImageInRGB(vp, clusterMemberships,
                                                              clusterCenters, null,
                                                              randomizationSeed));
            stacks[visualizationModes] = imgsbp2;
            visualizationModes++;
        }

        if (binaryStackVisualization)
        {
            stacks[visualizationModes] = encodeClusteredImageBinaryStack(vp, clusterMemberships,
                                                                         clusterCenters);
            visualizationModes++;
        }

        if (clusterCenterColorsVisualization)
        {
            float [][] computedClusterCenters = convertMatrixColorSpace(clusterCenters,
                                                                        colorSpace);
            stacks[visualizationModes] = encodeClusteredImageWithClusterCenterColors(vp, clusterMemberships,
                                                                                     computedClusterCenters);
            visualizationModes++;
        }

        if (fuzzyStackVisualization)
        {
            stacks[visualizationModes] = encodeClusteredImageFuzzyStack(vp, (float[][]) resultMatrixes[2],
                                                                         clusterCenters);
            visualizationModes++;
        }

        imgArray = new ImageStack[visualizationModes];
        for (int i = 0; i < visualizationModes; i++)
        {
            imgArray[i] = (ImageStack)stacks[i];
        }

        return imgArray;
    }


   /**
    * Converts the @link{ImagePlus} received in input into the color space specified
    * by the input string. If the string equals "<b>None</b>" then no color space
    * conversion is applied and the image is converted into a float stack.
    * To convert the color space it calls some utility methods from
    * @link{ColorSpaceConversion}.
    * @param imp the image to convert
    * @param colorSpace the string specifying the color space for the conversion.
    * Allowed values are: "<b>None</b>", "<b>XYZ</b>", "<b>L*a*b*</b>", "<b>HSB</b>"
    * @return a VectorProcessor representing the converted image
    * @see ColorSpaceConversion
    * @see #convertToFloatStack(ij.ImagePlus)
    * @see #convertToHSBFloatStack(ij.ImagePlus)
    */
   private VectorProcessor convertToColorSpace(ImagePlus imp, String colorSpace){

       VectorProcessor vp = null;

       if (colorSpace.equals("XYZ"))
       {
            vp = ColorSpaceConversion.rgbToXYZVectorProcessor((ColorProcessor)imp.getProcessor());
       }
       else if (colorSpace.equals("L*a*b*"))
       {
            vp = ColorSpaceConversion.rgbToLabVectorProcessor((ColorProcessor)imp.getProcessor());
       }
       else if (colorSpace.equals("HSB"))
       {
            final ImagePlus stack = convertToHSBFloatStack(imp);
            ImageStack ims = stack.getStack();
            System.out.println("Stack slices " + ims.getSize());
            vp = new VectorProcessor(ims);
//            float [] hsb = vp.get(3, 44);
//
//            System.out.println("HS "+ hsb[0] + " " + hsb[1] + " " + hsb[2]);
       }
       else // if (colorSpace.equals("None"))
       {
            final ImagePlus stack = convertToFloatStack(imp);
            ImageStack ims = stack.getStack();

            System.out.println("Stack slices " + ims.getSize());
            vp = new VectorProcessor(ims);
       }

       return vp;
   }

   /**
    * Computes an integer encoding the stop criterion type currently chosen.
    * It throws an exception if no valid criterion is chosen.
    * @return an integer representing the stop criterion
    * @see #stopCriterion
    * @see #stopCriterions
    */
   private int computeStopCriterion(){

       int stopCrit = -1;

       for (int i = 0; i < stopCriterions.length; i++)
       {
           if(stopCriterion.equals(stopCriterions[i]))
           {
               stopCrit = i;
           }
       }

       if(stopCrit == -1)
       {
          throw new IllegalArgumentException("Invalid Stop Criterion");
       }

       return stopCrit;
   }

   /**
    * Computes an integer encoding the initialization mode currently chosen.
    * It throws an exception if no valid mode is chosen.
    * @return an integer representing the init mode
    * @see #initModes
    * @see #initializationMode
    */
   private int computeInitializationMode(){

        int initMode = -1;

        for (int i = 0; i < initModes.length; i++)
        {
            if (initializationMode.equals(initModes[i]))
            {
                initMode = i;
            }
        }

        if (initMode == -1)
        {
            throw new IllegalArgumentException("Invalid Initialization Mode");
        }
        return initMode;
   }

   /**
    * Computes an integer encoding the spatial function kind currently chosen.
    * It throws an exception if no valid function kind is chosen.
    * @return an integer representing the spatial function
    * @see #spatialFunction
    * @see #spatialFunctions
    */
   private int computeSpatialFunctionMode(){

       int sFunction = -1;

       for (int i = 0; i < spatialFunctions.length; i++)
       {
           if(spatialFunction.equals(spatialFunctions[i]))
           {
               sFunction = i;
           }
       }

       if(sFunction == -1)
       {
           throw new IllegalArgumentException("Invalid Spatial Function Choice");
       }
       return sFunction;
   }

   /**
    * Encodes the input images into a new one representing the clustering partition
    * resulted from the algorithm by using a gray scaled image
    * @param vp the @link{VectorProcessor} representing the original image
    * @param clusterMemberships the clustering partition matrix
    * @param clusterCenters the cluster center matrix
    * @return a @link{ByteProcessor} containing the encoded image
    */
   private ByteProcessor encodeClusteredImageInGray(VectorProcessor vp,
                                                     int [][] clusterMemberships,
                                                     float [][] clusterCenters){

        final ByteProcessor dest = new ByteProcessor(vp.getWidth(), vp.getHeight());
        final VectorProcessor.PixelIterator iterator = vp.pixelIterator();
        final int nClusters = clusterCenters.length;
        while (iterator.hasNext())
        {
            final float[] v = iterator.next();
            for (int j = 0; j < nClusters; j++)
            {
                if (clusterMemberships[iterator.getOffset()][j] == 1)
                {
                    dest.putPixel(iterator.getX(), iterator.getY(), j);
                }
            }
        }
        return dest;
    }

   /**
    * Encodes the input images into a new one representing the clustering partition
    * resulted from the algorithm by using a stack whose slices are binary images
    * representing the membership of a pixel to a specified cluster (each slice
    * represents a different cluster)
    * @param vp the @link{VectorProcessor} representing the original image
    * @param clusterMemberships the clustering partition matrix
    * @param clusterCenters the cluster center matrix
    * @return an @link{ImageStack} containing the encoded binary stack
    */
   private ImageStack encodeClusteredImageBinaryStack(VectorProcessor vp,
                                                      int [][] clusterMemberships,
                                                      float [][] clusterCenters){

        final ImageStack dest = new ImageStack(vp.getWidth(), vp.getHeight());
        final VectorProcessor.PixelIterator iterator = vp.pixelIterator();
        final int nClusters = clusterCenters.length;

        ByteProcessor[] binaryClusters = new ByteProcessor[nClusters];
        for (int i = 0; i < nClusters; i++)
        {
            binaryClusters[i] = new ByteProcessor(vp.getWidth(), vp.getHeight());
        }

        while (iterator.hasNext())
        {
            final float[] v = iterator.next();
            for (int j = 0; j < nClusters; j++)
            {
                if (clusterMemberships[iterator.getOffset()][j] == 1)
                {
                    binaryClusters[j].putPixel(iterator.getX(), iterator.getY(), 0);
                }
                else
                {
                    binaryClusters[j].putPixel(iterator.getX(), iterator.getY(), 255);
                }
            }
        }

        for (int i = 0; i < nClusters; i++)
        {
            dest.addSlice("", binaryClusters[i]);
        }

        return dest;
    }

   /**
    * Encodes the input images into a new one representing the <i>fuzzy</i> clustering partition
    * resulted from the algorithm by using a stack whose slices are gray scaled images
    * representing the fuzzy membership of a pixel to a specified cluster (each slice
    * represents a different cluster)
    * @param vp the @link{VectorProcessor} representing the original image
    * @param clusterMemberships the fuzzy clustering membership matrix
    * @param clusterCenters the cluster center matrix
    * @return an @link{ImageStack} containing the encoded fuzzy stack
    */
   private ImageStack encodeClusteredImageFuzzyStack(VectorProcessor vp,
                                                      float [][] clusterMemberships,
                                                      float [][] clusterCenters){

        final ImageStack dest = new ImageStack(vp.getWidth(), vp.getHeight());
        final VectorProcessor.PixelIterator iterator = vp.pixelIterator();
        final int nClusters = clusterCenters.length;

        ByteProcessor[] fuzzyClusters = new ByteProcessor[nClusters];
        for (int i = 0; i < nClusters; i++)
        {
            fuzzyClusters[i] = new ByteProcessor(vp.getWidth(), vp.getHeight());
        }

        while (iterator.hasNext())
        {
            final float[] v = iterator.next();
            for (int j = 0; j < nClusters; j++)
            {
                float inverseValue = 1.0f - clusterMemberships[iterator.getOffset()][j];
                float grayValue = (int) (inverseValue * 255);
                fuzzyClusters[j].setf(iterator.getX(), iterator.getY(), grayValue);
            }
        }

        for (int i = 0; i < nClusters; i++)
        {
            dest.addSlice("Fuzzy", fuzzyClusters[i]);
        }

        return dest;
    }

    /**
     * Encodes an image into a RGB colored images according to the clustering results,
     * it assigns each cluster a random generated color using @link{ColorManger}
     * and each pixel will assume the color of the cluster it belongs to.
     * @param vp the original image as a @link{VectorProcessor}
     * @param clusterMemberships the cluster membership matrix
     * @param clusterCenters the cluster center matrix
     * @param mixingColor an optional color for mixing the random color sequence
     * @param randomizationSeed an integer used to generate the random colors
     * @return a @link{ColorProcessor} encoding the clustered image
     * @see ColorManager
     */
    private ColorProcessor encodeClusteredImageInRGB(VectorProcessor vp, int [][] clusterMemberships,
                                                     float [][] clusterCenters, int [] mixingColor,
                                                     int randomizationSeed){

        final ColorProcessor dest = new ColorProcessor(vp.getWidth(), vp.getHeight());
        final VectorProcessor.PixelIterator iterator = vp.pixelIterator();
        final int nClusters = clusterCenters.length;

        /* Creating an array of random RGB colors, one for each cluster */
        int [][] clusterColors = new int [nClusters][4];
        for (int i = 0; i < nClusters; i++)
        {
            clusterColors[i] = ColorManager.randomMixedRGBColor(mixingColor);
            System.out.println("C" + i + " " + clusterColors[i][0] +
                    " " + clusterColors[i][1] + " " + clusterColors[i][2]);
        }

        while (iterator.hasNext())
        {
            final float[] v = iterator.next();
            for (int j = 0; j < nClusters; j++)
            {
                if (clusterMemberships[iterator.getOffset()][j] == 1)
                {

                    dest.putPixel(iterator.getX(), iterator.getY(), clusterColors[j]);

                }
            }
        }

        System.out.println(dest.getColor(4, 33));
        return dest;
    }

    /**
     * Encodes an image into a RGB colored images according to the clustering results,
     * each pixel will assume the RGB color of the centroid whose cluster it belongs to.
     * @param vp the original image as a @link{VectorProcessor}
     * @param clusterMemberships the cluster membership matrix
     * @param clusterCenters the cluster center matrix
     * @return an @link{ImageStack} encoding the RGB clustered image
     */
    private ImageStack encodeClusteredImageWithClusterCenterColors(VectorProcessor vp,
                                                                   int [][] clusterMemberships,
                                                                   float [][] clusterCenters){

        final ImageStack dest = new ImageStack(vp.getWidth(), vp.getHeight());
        final VectorProcessor.PixelIterator iterator = vp.pixelIterator();
        final int nClusters = clusterCenters.length;

        for (int i = 0; i < vp.getNumberOfValues(); i++)
        {
            dest.addSlice("Band i", new FloatProcessor(vp.getWidth(), vp.getHeight()));
        }

        final Object[] pixels = dest.getImageArray();

        while (iterator.hasNext())
        {
            final float[] v = iterator.next();
            for (int j = 0; j < nClusters; j++)
            {
                if (clusterMemberships[iterator.getOffset()][j] == 1)
                {
                    for (int k = 0; k < vp.getNumberOfValues(); k++)
                    {
                        ((float[]) pixels[k])[iterator.getOffset()] = clusterCenters[j][k];
                    }
                }
            }
        }

        return dest;
    }

    /**
     * Converts a float sample x features matrix into the RGB color space according
     * to the current color space expressed in the input string. Allowed values
     * are: "<b>None</b>", "<b>XYZ</b>", "<b>L*a*b*</b>", "<b>HSB</b>". If
     * "<b>None</b>", than no conversion is applied. It is called when it is
     * necessary to encode the clustered matrix into a RGB image using the cluster
     * centroids colors.
     * @param dataFeatureMatrix the matrix representing the image converted
     * @param colorSpace the color space from which convert the matrix back
     * @return a float matrix resulting from the conversion
     * @see #encodeClusteredImageWithClusterCenterColors(vectorLib.VectorProcessor, int[][], float[][]) 
     */
    private float [][] convertMatrixColorSpace(float [][] dataFeatureMatrix, String colorSpace){

        float [][] convertedMatrix = null;

        if (colorSpace.equals("L*a*b*"))
        {
            convertedMatrix = ColorSpaceConversion.labToRGBMatrix(dataFeatureMatrix);
        }
        else if (colorSpace.equals("XYZ"))
        {
            convertedMatrix = ColorSpaceConversion.xyzToRGBMatrix(dataFeatureMatrix);
        }
        else if (colorSpace.equals("HSB"))
        {
            convertedMatrix = ColorSpaceConversion.hsbToRGBMatrix(dataFeatureMatrix);
        }
//        else if (colorSpace.equals("YCrCb"))
//        {
//        }
        else // if (colorSpace.equals("None"))
        {
            convertedMatrix  = dataFeatureMatrix;
        }

        return convertedMatrix;
    }

    /**
     * Convert image to a stack of FloatProcessors.
     *
     * @param src image to convert.
     * @return float stack.
     */
    private static ImagePlus convertToFloatStack(final ImagePlus src) {

        final ImagePlus dest = new Duplicator().run(src);

        // Remember scaling setup
        final boolean doScaling = ImageConverter.getDoScaling();

        try {
            // Disable scaling
            ImageConverter.setDoScaling(false);

            if (src.getType() == ImagePlus.COLOR_RGB)
            {
                if (src.getStackSize() > 1)
                {
                    throw new IllegalArgumentException("Unsupported image type: "
                            + "RGB with more than one slice.");
                }

                final ImageConverter converter = new ImageConverter(dest);
                converter.convertToRGBStack();
            }

            if (dest.getStackSize() > 1)
            {
                final StackConverter converter = new StackConverter(dest);
                converter.convertToGray32();

            }
            else
            {
                final ImageConverter converter = new ImageConverter(dest);
                converter.convertToGray32();
            }

            return dest;
        }
        finally
        {
            // Restore original scaling option
            ImageConverter.setDoScaling(doScaling);
        }
    }

    /**
     * Converts an @link{ImagePlus} into a float stack representing the HSB
     * conversion of the original image. Used when the HSB color space is chosen
     * for representing the input image.
     * @param src the original Image to convert
     * @return the converted ImagePlus containing the float stack
     * @see #convertToColorSpace(ij.ImagePlus, java.lang.String)
     */
    private static ImagePlus convertToHSBFloatStack(final ImagePlus src) {

        final ImagePlus dest = new Duplicator().run(src);

        // Remember scaling setup
        final boolean doScaling = ImageConverter.getDoScaling();

        try {
            // Disable scaling
            ImageConverter.setDoScaling(false);

            if (src.getType() == ImagePlus.COLOR_RGB)
            {
                if (src.getStackSize() > 1)
                {
                    throw new IllegalArgumentException("Unsupported image type: "
                            + "RGB with more than one slice.");
                }

                final ImageConverter converter = new ImageConverter(dest);
                converter.convertToHSB();
            }

            if (dest.getStackSize() > 1)
            {
                final StackConverter converter = new StackConverter(dest);
                converter.convertToGray32();

            }
            else
            {
                final ImageConverter converter = new ImageConverter(dest);
                converter.convertToGray32();
            }

            return dest;
        }
        finally
        {
            // Restore original scaling option
            ImageConverter.setDoScaling(doScaling);
        }
    }

    /**
     * Implements the ClusteringDelegate method of the same signature.
     * Updates the ImageJ status showing a message that contains the algorithm
     * iteration number and the iteration error
     * @param V the optional cluster center matrix (not used)
     * @param U the optional cluster membership matrix (not used)
     * @param nIteration the current iteration
     * @param error the iteration error
     * @see ClusteringDelegate
     * @see IJ
     */
    @Override
    public void updateStatus(float[][] V, int[][] U, long nIteration, float error) {

        final String message = "Fuzzy C-Means iteration " + nIteration + ", error: " + error;
        IJ.showStatus(message);
    }

    /**
     * Implements the ClusteringDelegate method of the same signature.
     * Updates the ImageJ status showing a textual message and optionally
     * printing it on the console
     * @param message the textual message to display
     * @see #printOnConsole
     * @see ClusteringDelegate
     */
    @Override
    public void updateStatus(String message) {

        IJ.showStatus(message);
        if (printOnConsole)
        {
            System.out.println(message);
        }
    }

    /**
     * Implements the ClusteringDelegate method of the same signature.
     * Updates the ImageJ status showing a message containing the algorithm
     * updated status and information on the matrixes and the iteration errors
     * obtained by checking the convergence. Optionally displays the message on
     * the system console and writes the results on a file if testing
     * @param V the optional cluster center matrix (not used)
     * @param U the optional cluster membership matrix (not used)
     * @param nIteration the current iteration
     * @param errorJ the current error in the objective function
     * @param errorU the current error computed on the matrix U
     * @param errorV the current error computed on the matrix V
     * @throws IOException exception thrown while trying to write on file
     * @see #RAF
     * @see ClusteringDelegate
     * @see #printOnConsole
     */
    @Override
    public void updateStatus(float[][] V, int[][] U, long nIteration, float errorJ,
                             float errorU, float errorV) throws IOException{
        final String message = "Fuzzy C-Means iteration " + nIteration
                                + ", errorJ: " + errorJ + ", errorU: " + errorU + ", errorV: " + errorV;
        IJ.showStatus(message);
        if (printOnConsole)
        {
            System.out.println(message);
        }

        String csvMessage = nIteration
                                + "," + errorJ + "," + errorU + "," + errorV + "\n";

        if (RAF != null)
        {
            RAF.writeBytes(csvMessage);
        }
    }

    /**
     * Writes on the <code>RAF</code> file the sequence of parameters when a validation
     * test is run after the @link{SFCM} algorithm finished.
     * @param m the fuzzyness value used
     * @param p the membership weight value used
     * @param q the spatial function weight value used
     * @param rad the window radius used
     * @param init the initialization criterion used
     * @param spatialFunction the spatial function chosen
     * @param k the used number of cluster
     * @param colorSpace the color space used
     * @param vpc the partition coefficient computed for validating the algorithm
     * @param vpe the partition entropy value computed for validation
     * @param vxb the xie and benn's index computed for validation
     * @throws IOException an exception if it it not posisble to write on the file
     * @see #RAF
     * @see #validation
     * @see #run
     * @see #TestManager
     * @see ClusterValidity
     */
    public void testingValidation(double m, double p, double q, int rad, String init,
                                  String spatialFunction, int k, String colorSpace,
                                  double vpc, double vpe, double vxb) throws IOException {
        final String message = k + "," + m + "," + p + "," + q + "," + rad + ","
                + init + "," + spatialFunction + "," + colorSpace + ","
                + vpc + "," + vpe + "," + vxb + "\n";

        RAF.writeBytes(message);
    }

    public String getConfigurationString(){
        String config = " k:" + numberOfClusters + " m:" + fuzzyness + " rnd:" +
                randomizationSeed + " init:" + initializationMode + " p:" + 
                membershipWeight + " q:" + spatialFunctionWeight + " r:" + windowRadius
                + " cs:" + colorSpace + " stc:" + stopCriterion;
        return config;
    }
}
