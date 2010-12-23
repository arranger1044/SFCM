/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package clustering.kmeans;

/**
 *
 * @author valerio
 */

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import ij.process.ImageConverter;
import ij.process.ImageProcessor;
import ij.process.StackConverter;




public class KMeansPlugin implements PlugIn{

    public static final String RESULTS_WINDOW_TITLE = "k-means cluster centers";

    private static boolean showCentroidImage;
    private static boolean sendToResultTable;

    private static final boolean APPLY_LUT = false;
    private static final boolean AUTO_BRIGHTNESS = true;


    private static final String TITLE = "k-means Clustering";
    private static final String ABOUT = "" +
            "Refactoring and Generalization of the k-means Plugin for ImageJ found" +
            " in the ij-plugin-toolkit";

    /* Creating an array of illegal Image types */
    private static final int[] illegalImages = {ImagePlus.COLOR_256};
    private KMeansManager KMM;

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
        KMM = new KMeansManager();

        /* Create an instance of Generic Dialog */
        GenericDialog configDialog = new GenericDialog("K-means Configuration");
        /* and configure it according to the defaults in KMeansManager */
        configDialog = configureDialog(KMM, configDialog);
        /* Show Dialog */
        configDialog.showDialog();
        if (configDialog.wasCanceled())
        {
            return;
        }

        /* Configuring the KMeansManager */
        getConfigurationFromDialog(configDialog, KMM);

        /* Calling the clustering algorithm and keeping the time */
        final long startTime = System.currentTimeMillis();
        //ImageProcessor[] imgArray = KMM.run(imp);
        ImageStack[] imgArray = KMM.run(imp);
        final long endTime = System.currentTimeMillis();

        // Apply default color map
//        if (APPLY_LUT)
//        {
//        //    bp.setColorModel(defaultColorModel());
//        }
//        if (AUTO_BRIGHTNESS)
//        {
//            for (int i = 0; i < imgArray.length; i++)
//            {
//                imgArray[i].setMinAndMax(0, KMM.getNumberOfClusters());
//            }
//        }

        /* Show result images */
        String type = null;
        for (int i = 0; i < imgArray.length; i++)
        {
            ImagePlus r = null;

//            if (imgArray[i] instanceof ByteProcessor)
//            {
//                type = "in Gray Scale";
//                System.out.println(((ByteProcessor)imgArray[i]).getPixel(4, 33));
//                r = new ImagePlus("Clusters " + type, (ByteProcessor)imgArray[i]);
//            }
//            else if (imgArray[i] instanceof ColorProcessor)
//            {
//                type = "in random RGB Colors";
//                System.out.println(((ColorProcessor)imgArray[i]).getColor(4, 33));
//                imgArray[i].convertToRGB();
//                System.out.println(((ColorProcessor)imgArray[i]).getColor(4, 33));
//                r = new ImagePlus("Clusters " + type, (ColorProcessor)imgArray[i]);
//            }

            System.out.println(i + "len " + imgArray.length + " " + imgArray[i].getSize());
            ImageProcessor IP = imgArray[i].getProcessor(1);
            int slices = imgArray[i].getSize();
            System.out.println("class " + IP.getClass());
            if (IP.getClass() == ij.process.ByteProcessor.class && slices == 1)
            {
                System.out.println("Entrato");
                adjustBrightness(IP, KMM.getNumberOfClusters());
                r = new ImagePlus("Clusters in Gray Scale", IP);
            }
            else if (IP.getClass() == ij.process.FloatProcessor.class)
            {
                final boolean doScaling = ImageConverter.getDoScaling();
                try
                {
                    ImageConverter.setDoScaling(false);
                    r = new ImagePlus("Clusters", imgArray[i]);
                    final StackConverter stackConverter = new StackConverter(r);
                    stackConverter.convertToGray8();
                    final ImageConverter imageConverter = new ImageConverter(r);
                    imageConverter.convertRGBStackToRGB();
                }
                finally
                {
                    ImageConverter.setDoScaling(doScaling);
                }

            }
            else
            {
                r = new ImagePlus("Clusters", imgArray[i]);
            }
            
            //System.out.println(type);
            
            r.show();
        }

        IJ.showStatus("Clustering completed in " + (endTime - startTime) + " ms.");

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

    private GenericDialog configureDialog(KMeansManager KMM, GenericDialog dialog){
        
        dialog.addNumericField("Number_of_clusters", KMM.getNumberOfClusters(), 0);
        dialog.addNumericField("Cluster_center_tolerance", KMM.getTolerance(), 8);
        dialog.addCheckbox("Enable_randomization_seed", KMM.isRandomizationSeedEnabled());
        dialog.addNumericField("Randomization_seed", KMM.getRandomizationSeed(), 0);
        dialog.addCheckbox("Show_clusters_as_centroid_value", showCentroidImage);
        return dialog;
    }

    private void getConfigurationFromDialog(GenericDialog dialog, KMeansManager KMM){
        KMM.setNumberOfClusters((int) Math.round(dialog.getNextNumber()));
        KMM.setTolerance((float) dialog.getNextNumber());
        KMM.setRandomizationSeedEnabled(dialog.getNextBoolean());
        KMM.setRandomizationSeed((int) Math.round(dialog.getNextNumber()));
        showCentroidImage = dialog.getNextBoolean();

    }

    private void adjustBrightness(ImageProcessor IP, int nClusters){
        if (AUTO_BRIGHTNESS)
        {
            IP.setMinAndMax(0, nClusters);
        }
    }
}

