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
 * This class communicates with the user by implementing the ImageJ PlugIn interface
 * It provides an graphical user interface to let the user configure the parameters
 * of the @link{SFCM} algorithm. Then, it configure the @{SFCMManager} instance
 * that will run the algorithm and in the end will display the resulted images
 * encoded by the manager instance.
 * @see SFCM
 * @see SFCMManager
 * @see #run(java.lang.String) 
 */
public class SFCMPlugin implements PlugIn{

    /**
     * The title of the resulting window
     */
    public static final String RESULTS_WINDOW_TITLE = " Spatial Fuzzy C-Means cluster centers";
    /**
     * The title displayed by the plugin
     */
    private static final String TITLE = "Spatial Fuzzy C-Means Clustering";
    /**
     * The plugin <b>about</b> message
     */
    private static final String ABOUT = "" +
            "Spatial Fuzzy C-Means in ImageJ, built upon the k-means plugin found" +
            " in the ij-plugin-toolkit";
    /**
     * An array of illegal Image types
     */
    private static final int[] illegalImages = {ImagePlus.COLOR_256};
    /**
     * The instance of the @link{SFCMManager}
     */
    private static SFCMManager staticSFCMM = new SFCMManager();

    /**
     *  Overridden method of the one provided by the @link{PlugIn} interface in ImageJ
     *  It takes the currently selected image, validates the image,
     *  validates the user inputs and instantiates the @link{SFCMManager} class
     *  which manages to run the @link{SFCM} algorithm, in the end it presents
     *  the clustered images to the user according to the selected visualization method.
     *  @param arg an optional string used to determine if the plugin is run for
     *  the <b>about</b> help
     *  @see #validateColorSpace(ij.ImagePlus, java.lang.String)
     *  @see #validateImage(ij.ImagePlus)
     *  @see #validateInputParameters(clustering.sfcm.SFCMManager)
     *  @see #configureDialog(clustering.sfcm.SFCMManager, ij.gui.GenericDialog)
     */
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

        SFCMManager sFCMM = null;
        /* Create an instance of KMeansManager */
        sFCMM = new SFCMManager();


        /* Create an instance of Generic Dialog */
        GenericDialog configDialog = new GenericDialog("Spatial Fuzzy C-Means Configuration");
        /* and configure it according to the defaults in KMeansManager */
        configDialog = configureDialog(staticSFCMM, configDialog);
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
                adjustBrightness(IP, sFCMM.getNumberOfClusters());
                r = new ImagePlus("Gray Scale" + sFCMM.getConfigurationString(), IP);
            }
            else if (IP.getClass() == ij.process.FloatProcessor.class)
            {

                r = encodeRGBImageFromStack(imp.getType(), imgArray[i]);

            }
            else
            {
                r = new ImagePlus("Clusters" + sFCMM.getConfigurationString(), imgArray[i]);
            }

            //System.out.println(type);

            r.show();
        }

        staticSFCMM = sFCMM;
    }

    /**
     * Validates the image type and in case of an invalid image it returns false.
     * Invalid image types are <b>ImagePlus.COLOR_256</b> as we can't work with those.
     * @param img the @link{ImagePlus}, input of the plugin
     * @return true if the image is valid , false otherwise
     */
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

    /**
     * Validates the color space conversion chosen in input.
     * If the user chose to convert a single sliced, non @link{ImagePlus.COLOR_RGB}
     * the input is invalid, if the chosen image is a multi slice image and the
     * user wants to convert it, the input is invalid, otherwise it is valid
     * in all the other case we will have a good input parameter
     * @param img the @link{ImagePlus}, input of the plugin
     * @param colorSpace the color space chosen for conversion allowed values are:
     * "<b>None</b>", "<b>XYZ</b>", "<b>L*a*b*</b>", "<b>HSB</b>"
     * @return true if parameters are valid otherwise false.
     */
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

    /**
     * Validates the user input parameters following the following constraints and buils
     * an error message. Calls @link{SFCMManager} to validate each value
     * 1 The fuzzyness value should be positive real number and strictly greater than 1
     * 2 The number of cluster should be greater than 1, it is a integer
     * 3 The tolerance should be greater than 0 , it can be a real number
     * 4 The membership weight should be a positive real number greater than 0
     * 5 The spatial function weight should be a positive real number greater than 0
     * 6 The radius of the window should be and integer and should be strictly greater then 0 
     * @param sFCMM the @link{SFCMManager} instace that runs the algorithm
     * @return An error message telling the user which parameter has been wrongly entered.
     * @see SFCMManager
     */
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

    /**
     * Configures an instance of @link{GenericDialog} adding the gui components
     * to let the user configure the algorithm parameters. To see a complete list
     * of parameters and their explanation , see @link{SFCMManager}
     * @param sFCMM the @link{SFCMManager} instace that runs the algorithm
     * @param the @link{GenericDialog} instance used to show the parameters
     * @return a configured dialog
     */
    private GenericDialog configureDialog(SFCMManager sFCMM, GenericDialog dialog){

        dialog.addCheckbox("Enable_random_seeding", sFCMM.getSeedEnabled());
        dialog.addNumericField("Randomization_seed", sFCMM.getRandomizationSeed(), 0);
        dialog.addNumericField("Number_of_clusters", sFCMM.getNumberOfClusters(), 0);
        dialog.addNumericField("Max_number_of_iterations", sFCMM.getMaxIterations(), 0);
        dialog.addChoice("Stop_criterion", sFCMM.getStopCriterions(), sFCMM.getStopCriterion());
        dialog.addNumericField("Tolerance_threshold", sFCMM.getTolerance(), 8);
        dialog.addChoice("Initialization Mode", sFCMM.getInitModes(), sFCMM.getInizializationMode());
//        dialog.addCheckbox("Random_clusters_initialization", sFCMM.getRandomInitialization());
//        dialog.addCheckbox("KMeans++_clusters_initialization", sFCMM.getKMeansPlusPlusInitialization());

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

    /**
     * This method takes all the input parameters from the interface dialog
     * and sest the proper variables in the @link{SFCMManager} instance which
     * will run the algorithm.
     * @param dialog the @link{GenericDialog} instance used to show the parameters
     * @param sFCMM the @link{SFCMManager} instace that runs the algorithm
     */
    private void getConfigurationFromDialog(GenericDialog dialog, SFCMManager sFCMM){
        sFCMM.setSeedEnabled(dialog.getNextBoolean());
        sFCMM.setRandomizationSeed((int) Math.round(dialog.getNextNumber()));
        sFCMM.setNumberOfClusters((int) Math.round(dialog.getNextNumber()));
        sFCMM.setMaxIterations(Math.round(dialog.getNextNumber()));
        sFCMM.setStopCriterion(dialog.getNextChoice());
        sFCMM.setTolerance((float) dialog.getNextNumber());
        sFCMM.setInitializationMode(dialog.getNextChoice());
//        sFCMM.setRandomInitialization(dialog.getNextBoolean());
//        sFCMM.setKMeansPlusPlusInitialization(dialog.getNextBoolean());

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

    /**
     * Adjusts the brightness of the pixels of a gray scaled image representing
     * a clustered image in order to make the cluster colors more distinctive.
     * @param IP the @link{ImageProcessor} representing the gray scale image
     * @param nClusters the number of clusters
     * @see #run(java.lang.String) 
     */
    private void adjustBrightness(ImageProcessor IP, int nClusters){

        IP.setMinAndMax(0, nClusters);
    }

    /**
     * Converts a stack into an RGB image or in a stack of gray of 8bit or 16bit or 32bit.
     * It is used to compute back a layered image obtained from @link{SFCMManager}
     * in order to visualize it
     * @param originalImageType integer rapresenting the type of the original image
     * @param centroidValueStack the stack of the clustered image
     * @return the new Image econded in the stack.
     */
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
