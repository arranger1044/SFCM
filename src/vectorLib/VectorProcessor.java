/*
 * Image/J Plugins
 * Copyright (C) 2002-2010 Jarek Sacha
 * Author's email: jsacha at users dot sourceforge dot net
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 * Latest release available at http://sourceforge.net/projects/ij-plugins/
 */
package vectorLib;

import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.ProgressBar;
import ij.plugin.Duplicator;
import ij.process.ColorProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageConverter;
import ij.process.StackConverter;
import util.Validate;

import java.awt.*;


/**
 * Represents vector valued image. Value at each pixel in the image is a vector of floating point
 * numbers.
 *
 * @author Jarek Sacha
 */
public class VectorProcessor {

    private final int width;
    private final int height;
    private final int numberOfValues; /* number of the features, component to describe a pixel */
    private final float[][] pixels; /* the rows are the pixels, the column are the features */
    private Rectangle roi;
    // TODO: use net.sf.ij_plugins.util.progress instead of ij.gui.ProgressBar for more flexibility.
    private ProgressBar progressBar;


    /**
     * Constructor that builds a VectorProcessor up given the width and height
     * of the ROI and the number of stacked images from which it has been created.
     * It gets called by the constructor that takes an ImageStack as parameter
     * @link{public VectorProcessor(final ImageStack stack)}
     *
     * @param width the width of the image ROI
     * @param height the height of the image ROI
     * @param numberOfValues the depth of the stack from which it has been created
     * i.e. the number of values (components, features) that describe each pixel
     */
    public VectorProcessor(final int width, final int height, final int numberOfValues) {
        this.width = width;
        this.height = height;
        this.numberOfValues = numberOfValues;
        pixels = new float[width * height][numberOfValues];
        roi = new Rectangle(0, 0, width, height);
    }


    /**
     * Constructor that builds a VectorProcessor up given a @link{ColorProcessor}
     * It calls the constructor @link{public VectorProcessor(final ImagePlus imp)}
     * passing as parameters an empty string and the @link{ColorProcessor}
     *
     * @param cp the @link{ColorProcessor} used to build up the VectorProcessor
     */
    public VectorProcessor(final ColorProcessor cp) {
        this(new ImagePlus("", cp));
    }


    /**
     * Constructor that builds a VectorProcessor up given an @link{ImagePlus} object
     * It calls the constructor @link{public VectorProcessor(final ImageStack stack)}
     * after converting the @link{ImagePlus} object into a @link{FloatProcessor}s stack
     * 
     * @param imp the @link{ImagePlus} object used to get the stack
     */
    public VectorProcessor(final ImagePlus imp) {
        this(convertToFloatStack(imp));
    }


    /**
     * Constructor that builds up a VectorProcessor froma a stack of
     * {@link FloatProcessor}s
     * This is the constructor that is being used in @link{KMeans}
     *
     * @param stack a stack of {@link FloatProcessor}s.
     */
    public VectorProcessor(final ImageStack stack) {
        /* Calls the constructor
         @link{VectorProcessor(final int width, final int height, final int numberOfValues)}*/
        this(stack.getWidth(), stack.getHeight(), stack.getSize());

        /* Gets the images in the stack and puts them into an array of Object*/
        final Object[] slices = stack.getImageArray();
        /* For each image (= component) in the stack */
        for (int i = 0; i < numberOfValues; ++i)
        {
            /* get the pixels as an array of float */
            final float[] values = (float[]) slices[i];
            /* for each pixel */
            for (int j = 0; j < values.length; j++)
            {
                /* store the value */
                pixels[j][i] = values[j];
            }
        }
    }


    /**
     * Public getter for the processor width
     * @return width of the image.
     */
    public int getWidth() {
        return width;
    }


    /**
     * Public getter for the processor height
     * @return height of the image.
     */
    public int getHeight() {
        return height;
    }


    /**
     * Public getter for the number of features of the processor
     * @return number of values at each pixel in the image.
     */
    public int getNumberOfValues() {
        return numberOfValues;
    }


    /**
     * Gives direct access to pixel values in the image.
     * The first index is the pixel number (between 0
     * and width*height-1), the second index references within each pixel value.
     *
     * @return reference to the array containing pixel values in the image.
     */
    public float[][] getPixels() {
        return pixels;
    }


    /**
     * Public getter to the roi internal reference
     * @return region of interest within the image.
     */
    public Rectangle getRoi() {
        return roi;
    }


    /**
     * Public setter of the internal roi reference
     * @param roi new ROI to be set
     * @see #getRoi()
     */
    public void setRoi(final Rectangle roi) {
        this.roi = roi;
    }


    /**
     * Returns the progress bar reference
     * @return
     */
    public ProgressBar getProgressBar() {
        return progressBar;
    }


    /**
     * Sets the progress bar reference
     * @param progressBar a reference to a @link{ProgressBar}
     */
    public void setProgressBar(final ProgressBar progressBar) {
        this.progressBar = progressBar;
    }


    /**
     * Returns an iterator through the pixels in the processor
     * @return pixel value iterator.
     * @see PixelIterator
     */
    public PixelIterator pixelIterator() {
        return new VectorProcessor.PixelIterator();
    }


    /**
     * returns an interator for a 3x4 Neighborhood for the pixels in the processor
     * @return pixel value iterator.
     * @see Neighborhood3x3
     * @see Iterator
     */
    public Iterator iterator() {
        return new VectorProcessor.Iterator();
    }


    /**
     * Convert VectorProcessor to an array of {@link FloatProcessor}'s.
     *
     * @return this VectorProcessor represented as an array of {@link FloatProcessor}'s
     * @see #toFloatStack()
     */
    public FloatProcessor[] toFloatProcessors() {
        final FloatProcessor[] r = new FloatProcessor[numberOfValues];
        for (int i = 0; i < numberOfValues; ++i)
        {
            final FloatProcessor fp = new FloatProcessor(width, height);
            final float[] values = (float[]) fp.getPixels();
            for (int j = 0; j < values.length; j++)
            {
                values[j] = pixels[j][i];
            }
            r[i] = fp;
        }

        return r;
    }


    /**
     * Convert VectorProcessor to ImagePlus with FloatProcessor stack.
     *
     * @return ImagePlus representation of this object.
     * @see #toFloatProcessors()
     */
    public ImagePlus toFloatStack() {
        final ImageStack stack = new ImageStack(width, height);
        final FloatProcessor[] fps = toFloatProcessors();
        for (int i = 0; i < fps.length; i++)
        {
            stack.addSlice("band " + i, fps[i]);
        }
        return new ImagePlus("From VectorProcessor", stack);
    }


    /**
     * Creates an @link{ImageStack} by converting an @link{ImagePlus} object
     * Duplicate of a method in @link{KmeansClusteringPlugin}
     *
     * @param src the @link{ImagePlus} object to convert
     * @return an @link{ImageStack} object containing the
     */
    private static ImageStack convertToFloatStack(final ImagePlus src) {

        final ImagePlus imp = duplicate(src);

        // Remember scaling setup
        final boolean doScaling = ImageConverter.getDoScaling();

        try
        {
            // Disable scaling
            ImageConverter.setDoScaling(false);

            if (imp.getType() == ImagePlus.COLOR_RGB)
            {
                if (imp.getStackSize() > 1)
                {
                    throw new RuntimeException("Unsupported image type: stack of COLOR_RGB");
                }

                final ImageConverter converter = new ImageConverter(imp);
                converter.convertToRGBStack();
            }

            if (imp.getStackSize() > 1)
            {
                final StackConverter converter = new StackConverter(imp);
                converter.convertToGray32();
            } 
            else
            {
                final ImageConverter converter = new ImageConverter(imp);
                converter.convertToGray32();
            }

            // FIXME: make sure that there are no memory leaks
//            imp.flush();
            return imp.getStack();

        } 
        finally
        {
            // Restore original scaling option
            ImageConverter.setDoScaling(doScaling);
        }
    }


    /**
     * Duplicates an @link{ImagePlus} object using the @link{Duplicator} class of ImageJ
     * 
     * @param imp the @link{ImagePlus} object to duplicate
     * @return the duplicated the @link{ImagePlus} object
     */
    private static ImagePlus duplicate(final ImagePlus imp) {
        final Duplicator duplicator = new Duplicator();
        return duplicator.run(imp);
    }


    /**
     * Return pixel value at coordinates (<code>x</code>, <code>y</code>).
     * It calls @link{float[] get(final int x, final int y, float[] dest)} giving
     * <code>x</code>, <code>y</code> as the first two parameters and <code>null</code> 
     * as the destination
     *
     * @param x x is the width
     * @param y y is the height
     * @return pixel value.
     */
    public float[] get(final int x, final int y) {
        return get(x, y, null);
    }


    /**
     * Return pixel value at coordinates (<code>x</code>, <code>y</code>).
     * Use {@code dest} to store the value.
     *
     * @param x    x is the width
     * @param y    y is the height
     * @param dest array to store pixel value, can be {@code null}. If null it wll be instantiated
     * @return pixel value. 
     */
    public float[] get(final int x, final int y, float[] dest) {

        /* If the values of x and y are not legal */
        if (x < 0 || x >= width || y < 0 || y >= height)
        {
            throw new IllegalArgumentException("Value of coordinates (x,y)=(" + x + "," + y
                    + ") is out of range [" + width + "," + height + "].");
        }

        /* If dest is null, instantiate it */
        if (dest == null)
        {
            dest = new float[numberOfValues];
        } 
        else
        {
            if (dest.length != numberOfValues)
            {
                throw new IllegalArgumentException("Invalid length of array dest.");
            }
        }

        /* Computing the row position in pixels and copying the array into dest */
        final int offset = x + y * width;
        final float[] v = pixels[offset];
        System.arraycopy(v, 0, dest, 0, v.length);
        return dest;
    }


    /**
     * Set value of pixel value at coordinates (<code>x</code>, <code>y</code>).
     *
     * @param x x is the width
     * @param y y is the height
     * @param v is the pixel value (and array fo floating features)
     */
    public void set(final int x, final int y, float[] v) {

        /* If the values of x and y are not legal */
        if (x < 0 || x >= width || y < 0 || y >= height)
        {
            throw new IllegalArgumentException("Value of coordinates (x,y) is out of range.");
        }

        /* Check if the value vector v is not null and has a proper lenght */
        Validate.argumentNotNull(v, "v");

        if (v.length != numberOfValues)
        {
            throw new IllegalArgumentException("Invalid size of argument 'v' expecting " + numberOfValues
                    + ", got " + v.length + ".");
        }

        /* Computing the row position in pixels and copying the array into v */
        final int offset = x + y * width;
        final float[] s = pixels[offset];
        System.arraycopy(v, 0, s, 0, v.length);
    }


    /**
     * Duplicate method for the VectorProcessor class
     * Calls the constructor @link{VectorProcessor(final int width, final int height, final int numberOfValues)},
     * then copies the pixel values
     *
     * @return and object that is the duplicate of the instance on which it is called
     */
    public VectorProcessor duplicate() {
        
        /* Create a new VectorProcessor */
        final VectorProcessor r = new VectorProcessor(this.width, this.height, this.numberOfValues);
        /* Clones the ROI */
        r.roi = (Rectangle) (roi != null ? roi.clone() : null);
        // TODO: ignore progress bar?
        r.progressBar = null;

        /* For each pixel copy its features */
        for (int i = 0; i < pixels.length; ++i)
        {
            System.arraycopy(pixels[i], 0, r.pixels[i], 0, numberOfValues);
        }

        return r;
    }


    /**
     * Class that implements @link{java.util.Iterator} in order to create an Iterator
     * that can iterate over the pixels.
     */
    public class PixelIterator implements java.util.Iterator<float[]> {

        final int xMin = roi.x;
        final int xMax1 = roi.x + roi.width - 1;
        final int rowOffset = width;
        final int yMin = roi.y;
        final int yMax1 = roi.y + roi.height - 1;
        int x = roi.x - 1;
        int y = roi.y;


        private PixelIterator() {
        }


        /**
         * Returns the x coordinate of the current pixel
         * @return x coordinate of the current pixel
         */
        public int getX() {
            if (x < xMin || x > xMax1)
            {
                throw new IllegalStateException("Illegal value of x, " + x + ".");
            }
            return x;
        }


        /**
         * Returns the y coordinate of the current pixel
         * @return y coordinate of the current pixel
         */
        public int getY() {
            if (y < yMin || y > yMax1)
            {
                throw new IllegalStateException("Illegal value of y, " + y + ".");
            }
            return y;
        }


        /**
         * Returns whether there is another pixel to iterate through
         * @return a boolean that tells if there are more pixel to iterate through
         */
        @Override
        public boolean hasNext() {
            return x < xMax1 || y < yMax1;
        }


        /**
         * Returns the next pixel to iterate through
         * @return the next pixel value represented ad a float array (its features)
         */
        @Override
        public float[] next() {
            // Update center location
            if (x < xMax1)
            {
                ++x;
            } 
            else
            {
                if (y < yMax1)
                {
                    x = xMin;
                    ++y;
                }

                if (progressBar != null)
                {
                    progressBar.show(y - yMin, yMax1 - yMin);
                }
            }

            final int offset = x + y * width;

            return pixels[offset];
        }


        /**
         * Not supported.
         */
        @Override
        public void remove() {
            throw new UnsupportedOperationException("Method remove() not supported.");
        }


        /**
         * Computes the row value for the current pixel in the matrix <code>pixels</code>
         * @return the row index in <code>pixels</code>
         */
        public int getOffset() {
            return x + y * width;
        }
    }

    
    /**
     * Represents 3x3 neighborhood. the center pixel is <code>p5</code>. Pixels <code>p1</code> to
     * <code>p3</code> are in the top row, <code>p4</code> to <code>p6</code> in the middle, and
     * <code>p7</code> to <code>p9</code> in the bottom of the neighborhood.
     */
    public static class Neighborhood3x3 {

        float[] p1,
                p2,
                p3,
                p4,
                p5,
                p6,
                p7,
                p8,
                p9;
        int x,
            y,
            offset;
    }

    /**
     * Iterator over 3x3 neighborhood of vector valued pixels.
     */
    public class Iterator implements java.util.Iterator<Neighborhood3x3> {

        final int xMin = Math.max(roi.x, 1);
        final int xMax = Math.min(roi.x + roi.width, width - 1) - 1;
        final int rowOffset = width;
        final int yMin = Math.max(roi.y, 1);
        final int yMax = Math.min(roi.y + roi.height, height - 1) - 1;
        int x = xMin - 1;
        int y = yMin;
        final Neighborhood3x3 neighborhood3x3 = new Neighborhood3x3();


        private Iterator() {
        }


        @Override
        public boolean hasNext() {
            return x < xMax || y < yMax;
        }


        @Override
        public Neighborhood3x3 next() {
            // Update center location
            if (x < xMax)
            {
                ++x;
            } 
            else
            {
                if (y < yMax)
                {
                    x = xMin;
                    ++y;
                }

                if (progressBar != null)
                {
                    progressBar.show(y - yMin, yMax - yMin);
                }
            }
            
            int offset = x + y * width;

            // Update neighbourhood information
            neighborhood3x3.p1 = pixels[offset - rowOffset - 1];
            neighborhood3x3.p2 = pixels[offset - rowOffset];
            neighborhood3x3.p3 = pixels[offset - rowOffset + 1];

            neighborhood3x3.p4 = pixels[offset - 1];
            neighborhood3x3.p5 = pixels[offset];
            neighborhood3x3.p6 = pixels[offset + 1];

            neighborhood3x3.p7 = pixels[offset + rowOffset - 1];
            neighborhood3x3.p8 = pixels[offset + rowOffset];
            neighborhood3x3.p9 = pixels[offset + rowOffset + 1];

            neighborhood3x3.x = x;
            neighborhood3x3.y = y;
            neighborhood3x3.offset = offset;

            return neighborhood3x3;
        }


        /**
         * Not supported.
         */
        @Override
        public void remove() {
            throw new UnsupportedOperationException("Method remove() not supported.");
        }

    }
}
