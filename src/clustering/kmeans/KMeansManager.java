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
import ij.process.ImageConverter;
import ij.process.StackConverter;
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

    public ByteProcessor[] run(ImagePlus img)
    {
        ByteProcessor[] bpArray = null;

        final ImagePlus stack = convertToFloatStack(img);
        ImageStack ims = stack.getStack();

        VectorProcessor vp = new VectorProcessor(ims);

        float [][] imageData = vp.getPixels();

        Object[] resultMatrixes = KMeans.run(imageData, numberOfClusters, tolerance, randomizationSeed);

        int [][] clusterMemberships = (int[][]) resultMatrixes[0];
        float [][] clusterCenters = (float[][]) resultMatrixes[1];

        bpArray = new ByteProcessor[1];

        bpArray[0] = encodeClusteredImage(vp, clusterMemberships, clusterCenters);

        return bpArray;
    }

    private ByteProcessor encodeClusteredImage(VectorProcessor vp, int [][] clusterMemberships, float [][] clusterCenters){

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
