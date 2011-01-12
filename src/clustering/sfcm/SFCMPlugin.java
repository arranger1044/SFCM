/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package clustering.sfcm;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import ij.process.ImageConverter;
import ij.process.ImageProcessor;
import ij.process.StackConverter;

/**
 *
 * @author valerio
 */
public class SFCMPlugin implements PlugIn{

    public static final String RESULTS_WINDOW_TITLE = " Spatial Fuzzy C-Means cluster centers";

    private static boolean showCentroidImage;
    private static boolean sendToResultTable;

    private static final boolean APPLY_LUT = false;
    private static final boolean AUTO_BRIGHTNESS = true;


    private static final String TITLE = "Spatial Fuzzy C-Means Clustering";
    private static final String ABOUT = "" +
            "Spatial Fuzzy C-Means in ImageJ, built upon the k-means plugin found" +
            " in the ij-plugin-toolkit";

    /* Creating an array of illegal Image types */
    private static final int[] illegalImages = {ImagePlus.COLOR_256};
    private SFCMManager sFCMM;

    @Override
    public void run(String arg) {

        /* Showing the About Message */
        if ("about".equalsIgnoreCase(arg))
        {
            IJ.showMessage("About " + TITLE, ABOUT);
            return;
        }

        /* Getting the current Image */
        final ImagePlus imp = IJ.getImage();

        /* Validate it against the illegal Image Type List */
        if (!validateImage(imp))
        {
            IJ.error(TITLE, "Image type not supported");
            return;
        }

        /* Create an instance of KMeansManager */
        sFCMM = new SFCMManager();

        /* Create an instance of Generic Dialog */
        GenericDialog configDialog = new GenericDialog("Fuzzy C-Means Configuration");
        /* and configure it according to the defaults in KMeansManager */
        configDialog = configureDialog(sFCMM, configDialog);
        /* Show Dialog */
        configDialog.showDialog();
        if (configDialog.wasCanceled())
        {
            return;
        }

        /* Configuring the KMeansManager */
        getConfigurationFromDialog(configDialog, sFCMM);

        /* validate input parameters */
        String validParameters = validateInputParameters(sFCMM);
        if(validParameters != null)
        {
            IJ.error(TITLE, validParameters);
            return;
        }

        /* Validate the correct color space conversion */
        if (!validateColorSpace(imp, sFCMM.getColorSpace()))
        {
            IJ.error(TITLE, "Invalid Color Space Conversion");
            return;
        }

        /* Calling the clustering algorithm */
        final long startTime = System.currentTimeMillis();
        ImageStack[] imgArray = sFCMM.run(imp);
        final long endTime = System.currentTimeMillis();
        final String message = "Clustering and Visualization completed in " +
                                                      (endTime - startTime) +
                                                                     " ms.";
        IJ.showStatus(message);

        /* Show result images */
        String type = null;
        for (int i = 0; i < imgArray.length; i++)
        {
            ImagePlus r = null;

            //System.out.println(i + " len " + imgArray.length + " " + imgArray[i].getSize());
            ImageProcessor IP = imgArray[i].getProcessor(1);
            int slices = imgArray[i].getSize();
            //System.out.println("class " + IP.getClass());
            if (IP.getClass() == ij.process.ByteProcessor.class && slices == 1)
            {
                //System.out.println("Entrato");
                adjustBrightness(IP, sFCMM.getNumberOfClusters());
                r = new ImagePlus("Clusters in Gray Scale", IP);
            }
            else if (IP.getClass() == ij.process.FloatProcessor.class)
            {

                r = encodeRGBImageFromStack(imp.getType(), imgArray[i]);

            }
            else
            {
                r = new ImagePlus("Clusters", imgArray[i]);
            }

            //System.out.println(type);

            r.show();
        }
    }

    private boolean validateImage(ImagePlus img){

        boolean validImage = true;
        for (int i = 0; i < illegalImages.length && validImage; i++)
        {
            if (img.getType() == illegalImages[i])
            {
                validImage = false;
            }
        }

        return validImage;
    }

    private boolean validateColorSpace(ImagePlus img, String colorSpace){

        boolean validColorSpace = true;
        if (img.getStackSize() > 1 && !colorSpace.equals("None"))
        {
            validColorSpace = false;
        }
        else if (img.getStackSize() == 1 && !(img.getType() == ImagePlus.COLOR_RGB) && !colorSpace.equals("None"))
        {
            validColorSpace = false;
        }
        else
        {
            validColorSpace = true;
        }

        return validColorSpace;
    }

    private String validateInputParameters(SFCMManager sFCMM){
        
        String errorMessage = null;

        if (!sFCMM.validateFuzzyness())
        {
            errorMessage = "Invalid Fuzzyness value!";
        }

        if (!sFCMM.validateClusterNumber())
        {
            if (errorMessage != null)
            {
                errorMessage += "\nInvalid Number of Clusters!";
            }
            else
            {
                errorMessage = "\nInvalid Number of Clusters!";
            }
        }

        if (!sFCMM.validateTolerance())
        {
            if (errorMessage != null)
            {
                errorMessage += "\nInvalid Tolerance value!";
            }
            else
            {
                errorMessage = "\nInvalid Tolerance value!";
            }
        }

        if (!sFCMM.validateMembershipWeight())
        {
            if (errorMessage != null)
            {
                errorMessage += "\nInvalid Membership Weight value!";
            }
            else
            {
                errorMessage = "\nInvalid Membership Weight value!";
            }
        }

        if (!sFCMM.validateSpatialFunctionWeight())
        {
            if (errorMessage != null)
            {
                errorMessage += "\nInvalid Spatial Function Weight value!";
            }
            else
            {
                errorMessage = "\nInvalid Spatial Function Weight value!";
            }
        }

        if (!sFCMM.validateWindowsRadius())
        {
            if (errorMessage != null)
            {
                errorMessage += "\nInvalid Window Radius value!";
            }
            else
            {
                errorMessage = "\nInvalid Window Radius value!";
            }
        }
        return errorMessage;
    }

    private GenericDialog configureDialog(SFCMManager sFCMM, GenericDialog dialog){

        dialog.addNumericField("Number_of_clusters", sFCMM.getNumberOfClusters(), 0);
        dialog.addNumericField("Max_number_of_iterations", sFCMM.getMaxIterations(), 0);
        dialog.addChoice("Stop_criterion", sFCMM.getStopCriterions(), sFCMM.getStopCriterion());
        dialog.addNumericField("Tolerance_threshold", sFCMM.getTolerance(), 8);
        dialog.addChoice("Initialization Mode", sFCMM.getInitModes(), sFCMM.getInizializationMode());
//        dialog.addCheckbox("Random_clusters_initialization", sFCMM.getRandomInitialization());
//        dialog.addCheckbox("KMeans++_clusters_initialization", sFCMM.getKMeansPlusPlusInitialization());
        dialog.addNumericField("Randomization_seed", sFCMM.getRandomizationSeed(), 0);
        dialog.addNumericField("Fuzzyness_parameter", sFCMM.getFuzzyness(), 1);
        dialog.addNumericField("Window_radius", sFCMM.getWindowRadius(), 0);
        dialog.addNumericField("Membership_weight", sFCMM.getMembershipWeight(), 1);
        dialog.addChoice("Spatial_function", sFCMM.getSpatialFunctions(), sFCMM.getSpatialFunction());
        dialog.addNumericField("Spatial_function_weight", sFCMM.getSpatialFunctionWeight(), 1);
        dialog.addChoice("Color_Space_Conversion", sFCMM.getColorSpaces(), sFCMM.getColorSpace());
        dialog.addCheckbox("Show_clusters_as_centroid_value", sFCMM.getClusterCenterColorsVisualization());
        dialog.addCheckbox("Show_clusters_as_random_RGB", sFCMM.getRandomRGBVisualization());
        dialog.addCheckbox("Show_clusters_as_gray_levels", sFCMM.getGrayScaleVisualization());
        dialog.addCheckbox("Show_clusters_as_binary_stack", sFCMM.getBinaryStackVisualization());
        dialog.addCheckbox("Show_clusters_as_fuzzy_stack", sFCMM.getFuzzyStackVisualization());
        return dialog;
    }

    private void getConfigurationFromDialog(GenericDialog dialog, SFCMManager sFCMM){
        sFCMM.setNumberOfClusters((int) Math.round(dialog.getNextNumber()));
        sFCMM.setMaxIterations(Math.round(dialog.getNextNumber()));
        sFCMM.setStopCriterion(dialog.getNextChoice());
        sFCMM.setTolerance((float) dialog.getNextNumber());
        sFCMM.setInitializationMode(dialog.getNextChoice());
//        sFCMM.setRandomInitialization(dialog.getNextBoolean());
//        sFCMM.setKMeansPlusPlusInitialization(dialog.getNextBoolean());
        sFCMM.setRandomizationSeed((int) Math.round(dialog.getNextNumber()));
        sFCMM.setFuzzyness(dialog.getNextNumber());
        sFCMM.setWindowRadius((int) Math.round(dialog.getNextNumber()));
        sFCMM.setMembershipWeight(dialog.getNextNumber());
        sFCMM.setSpatialFunction(dialog.getNextChoice());
        sFCMM.setSpatialFunctionWeight(dialog.getNextNumber());
        sFCMM.setColorSpace(dialog.getNextChoice());
        sFCMM.setClusterCenterColorsVisualization(dialog.getNextBoolean());
        sFCMM.setRandomRGBVisualization(dialog.getNextBoolean());
        sFCMM.setGrayScaleVisualization(dialog.getNextBoolean());
        sFCMM.setBinaryStackVisualization(dialog.getNextBoolean());
        sFCMM.setFuzzyStackVisualization(dialog.getNextBoolean());
    }

    private void adjustBrightness(ImageProcessor IP, int nClusters){
        if (AUTO_BRIGHTNESS)
        {
            IP.setMinAndMax(0, nClusters);
        }
    }

    private ImagePlus encodeRGBImageFromStack(final int originalImageType, final ImageStack centroidValueStack){

        final boolean doScaling = ImageConverter.getDoScaling();
        try
        {
            ImageConverter.setDoScaling(false);
            final ImagePlus cvImp = new ImagePlus("Cluster centroid values", centroidValueStack);

            if (centroidValueStack.getSize() > 1)
            {
                final StackConverter stackConverter = new StackConverter(cvImp);

                switch (originalImageType)
                {
                    case ImagePlus.COLOR_RGB:
                        stackConverter.convertToGray8();
                        final ImageConverter imageConverter = new ImageConverter(cvImp);
                        imageConverter.convertRGBStackToRGB();
                        break;
                    case ImagePlus.GRAY8:
                        stackConverter.convertToGray8();
                        break;
                    case ImagePlus.GRAY16:
                        stackConverter.convertToGray16();
                        break;
                    case ImagePlus.GRAY32:
                        // No action needed
                        break;
                    default:
                        throw new IllegalArgumentException("Unsupported input image type: " + originalImageType);
                }
            }
            else
            {
                final ImageConverter converter = new ImageConverter(cvImp);
                // Convert image back to original type
                switch (originalImageType)
                {
                    case ImagePlus.COLOR_RGB:
                        throw new IllegalArgumentException("Internal error: RGB image cannot have a single band.");
                    case ImagePlus.GRAY8:
                        converter.convertToGray8();
                        break;
                    case ImagePlus.GRAY16:
                        converter.convertToGray16();
                        break;
                    case ImagePlus.GRAY32:
                        // No action needed
                        break;
                    default:
                        throw new IllegalArgumentException("Unsupported input image type: " + originalImageType);
                }
            }

            return cvImp;
        }
        finally
        {
            ImageConverter.setDoScaling(doScaling);
        }
    }

}
