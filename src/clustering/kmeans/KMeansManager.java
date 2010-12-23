/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package clustering.kmeans;

/**
 *
 * @author valerio
 */
import ij.process.ByteProcessor;
import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.Duplicator;
import ij.process.ColorProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageConverter;
import ij.process.ImageProcessor;
import ij.process.StackConverter;
import util.ColorManager;
import vectorLib.VectorProcessor;


public class KMeansManager {

   private int randomizationSeed = 48;
   private boolean randomizationSeedEnabled = true;
   private double tolerance = 0.0001;
   private int numberOfClusters = 4;

   public int getRandomizationSeed() {
        return randomizationSeed;
    }

    public void setRandomizationSeed(final int randomizationSeed) {
        this.randomizationSeed = randomizationSeed;
    }

    /**
     * If <code>true</code>, random number generator will be initialized with a
     * <code>randomizationSeed</code>. If <code>false</code> random number generator will be
     * initialized using 'current' time.
     *
     * @return {@code true} when randomization seed is enabled.
     * @see #getRandomizationSeed()
     */
    public boolean isRandomizationSeedEnabled() {
        return randomizationSeedEnabled;
    }

    public void setRandomizationSeedEnabled(final boolean randomizationSeedEnabled) {
        this.randomizationSeedEnabled = randomizationSeedEnabled;
    }

    public int getNumberOfClusters() {
        return numberOfClusters;
    }

    public void setNumberOfClusters(final int numberOfClusters) {
        this.numberOfClusters = numberOfClusters;
    }

    /**
     * Return tolerance used to determine cluster centroid distance. This tolerance is used to
     * determine if a centroid changed location between iterations.
     *
     * @return cluster centroid location tolerance.
     */
    public double getTolerance() {
        return tolerance;
    }

    public void setTolerance(final float tolerance) {
        this.tolerance = tolerance;
    }

    public ImageStack[] run(ImagePlus img)
    {
        ImageStack[] imgArray = null;

        final ImagePlus stack = convertToFloatStack(img);
        ImageStack ims = stack.getStack();

        VectorProcessor vp = new VectorProcessor(ims);

        float [][] imageData = vp.getPixels();

        Object[] resultMatrixes = KMeans.run(imageData, numberOfClusters, tolerance, randomizationSeed);

        int [][] clusterMemberships = (int[][]) resultMatrixes[0];
        float [][] clusterCenters = (float[][]) resultMatrixes[1];

        imgArray = new ImageStack[4];

        ImageStack imgsbp = new ImageStack(vp.getWidth(), vp.getHeight());
        imgsbp.addSlice("GrayScale", encodeClusteredImageInGray(vp, clusterMemberships, clusterCenters));
        imgArray[0] = imgsbp;

        ImageStack imgsbp2 = new ImageStack(vp.getWidth(), vp.getHeight());
        imgsbp2.addSlice("RGB", encodeClusteredImageInRGB(vp, clusterMemberships, clusterCenters, null, randomizationSeed));
        imgArray[1] = imgsbp2;

        imgArray[2] = encodeClusteredImageBinaryStack(vp, clusterMemberships, clusterCenters);

        imgArray[3] = encodeClusteredImageWithClusterCenterColors(vp, clusterMemberships, clusterCenters);
//        System.out.println("GR " + ((ByteProcessor)imgArray[0]).getPixel(4, 33));
//        imgArray[1] = encodeClusteredImageInRGB(vp, clusterMemberships, clusterCenters, null, randomizationSeed);
//        System.out.println("RGC " + ((ColorProcessor)imgArray[1]).getColor(4, 33));

        return imgArray;
    }



    private ByteProcessor encodeClusteredImageInGray(VectorProcessor vp, int [][] clusterMemberships, float [][] clusterCenters){

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

   private ImageStack encodeClusteredImageBinaryStack(VectorProcessor vp, int [][] clusterMemberships, float [][] clusterCenters){

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
            System.out.println("C" + i + " " + clusterColors[i][0] + " " + clusterColors[i][1] + " " + clusterColors[i][2]);
        }

        while (iterator.hasNext())
        {
            final float[] v = iterator.next();
            for (int j = 0; j < nClusters; j++)
            {
                if (clusterMemberships[iterator.getOffset()][j] == 1)
                {
                    //System.out.println("C" + j + " " + clusterColors[j][0] + " " + clusterColors[j][1] + " " + clusterColors[j][2]);
                    dest.putPixel(iterator.getX(), iterator.getY(), clusterColors[j]);
                    //dest.putPixel(iterator.getX(), iterator.getY(), clusterColors[j][0]*clusterColors[j][1]*clusterColors[j][2]);
                }
            }
        }

        System.out.println(dest.getColor(4, 33));
        return dest;
    }

    private ImageStack encodeClusteredImageWithClusterCenterColors(VectorProcessor vp, int [][] clusterMemberships, float [][] clusterCenters){

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
                    throw new IllegalArgumentException("Unsupported image type: RGB with more than one slice.");
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

}
