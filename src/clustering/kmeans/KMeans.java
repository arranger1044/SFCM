/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package clustering.kmeans;

import java.awt.Point;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 *
 * @author valerio
 */
public class KMeans {

    private KMeans(){
    }

    public static Object[] run (float [][] X, int k, double tolerance,
                                int randomSeed){

        float [][] V = new float[k][X[0].length];
        int [][] U = new int [X.length][k];

        V = initializeClusterCenterMatrix(X, V, k, randomSeed);
        U = computeClusterMembership(X, U, V, k);

        return clusterize(X, U, V, k, tolerance);
    }

    private static float [][] initializeClusterCenterMatrix(float [][] X,
                                                            float [][] V,
                                                            int k,
                                                            int randomSeed){

        V = kMeansPlusPlusInitialization(X, V, k, randomSeed);
        return V;
    }

    private static float[][] kMeansPlusPlusInitialization(float [][] X, float [][] V, int k, int randomSeed) {
        final Random random = createRandom(randomSeed);

        //final int nbClusters = config.getNumberOfClusters();
//        final int width = vp.getWidth();
//        final int height = vp.getHeight();
//        final int nbPixels = width * height;
        final int nbPixels = X.length;
        final int nFeatures = X[0].length;
//        // Cluster centers
//        final List<float[]> centers = new ArrayList<float[]>();
        // Location of pixels used as cluster centers
        final List<Integer> centerLocation = new ArrayList<Integer>(); // contiene i centri già calcolati
        int clusterCreated = 0;
        // Choose one center uniformly at random from among pixels
        {
//            final Point p = toPoint(random.nextInt(nbPixels), width);
//            centerLocation.add(p);
//            centers.add(vp.get(p.x, p.y));
            final int p = random.nextInt(nbPixels);
            centerLocation.add(p);
            System.arraycopy(X[p], 0, V[0], 0, nFeatures);
            //V[0] = X[p];
            clusterCreated++;
        }

        final double[] dp2 = new double[nbPixels];
        while (clusterCreated < k) {
            //assert centers.size() == centerLocation.size();

            // For each data point p compute D(p), the distance between p and the nearest center that
            // has already been chosen.
            double sum = 0;
            //final float[][] centersArray = centers.toArray(new float[centers.size()][]);
            for (int offset = 0; offset < nbPixels; offset++) {  
                                                                 
                //final Point p = toPoint(offset, width);

                // Test that this is not a repeat of already selected center 
//                if (centerLocation.contains(p)) {
//                    continue;
//                }
                if (centerLocation.contains(offset)) {
                    continue;
                }

                // Distance to closest cluster
                //final float[] v = vp.get(p.x, p.y);
                final float[] v = X[offset];
//                final int cci = closestCluster(v, centersArray);
//                final double d = distance(v, centersArray[cci]);
                final int cci = closestCluster(v, V, clusterCreated);
                final double d = distance(v, V[cci]);
                sum += d * d;
                dp2[offset] = sum;  // distribuzione cumulata   vettore di pixel nel quale metti la somma = d*d del pixel corrente e di tutti i precedenti.
            }


            // Add one new data point at random as a new center, using a weighted probability distribution where
            // a point p is chosen with probability proportional to D(p)^2
            final double r = random.nextDouble() * sum;
            for (int offset = 0; offset < nbPixels; offset++) {
                //final Point p = toPoint(offset, width);

                // Test that this is not a repeat of already selected center
//                if (centerLocation.contains(p)) {
//                    continue;
//                }
                if (centerLocation.contains(offset)) {
                    continue;
                }

                if (dp2[offset] >= r)
                {
                    centerLocation.add(offset);
//                    final float[] v = vp.get(p.x, p.y);
//                    centers.add(v);
                    System.arraycopy(X[offset], 0, V[clusterCreated], 0, nFeatures);
                    clusterCreated++;
                    break;
                }
            }
        }

        //return centers.toArray(new float[centers.size()][]);
        return V;
    }

    private static int [][] computeClusterMembership(float [][] X, int [][] U, float [][] V, int k){

        for (int i = 0; i < X.length; i++)
        {
            final int c = closestCluster(X[i], V, k);
            for (int j = 0; j < k; j++)
            {
                if(j == c)
                {
                    U[i][j] = 1;
                }
                else
                {
                    U[i][j] = 0;
                }
            }
        }
        return U;
    }

    private static Random createRandom(int randomSeed) {
//        return config.isRandomizationSeedEnabled()
//                ? new Random(config.getRandomizationSeed())
//                : new Random();
        return new Random(randomSeed);
    }

    /**
     * Return index of the closest cluster to point <code>x</code>.
     *
     * @param x              point features.
     * @param clusterCenters cluster centers features.
     * @return index of the closest cluster
     */
    private static int closestCluster(final float[] x, final float[][] clusterCenters, int k) {

        double minDistance = Double.MAX_VALUE;
        int closestCluster = -1;

        for (int i = 0; i < k; i++)
        {
            final float[] clusterCenter = clusterCenters[i];
            final double d = distance(clusterCenter, x);
            if (d < minDistance)
            {
                minDistance = d;
                closestCluster = i;
            }
        }

        return closestCluster;
    }


    /**
     * Distance between points <code>a</code> and <code>b</code>.
     * This is the squared rooted version
     * @param a first point.
     * @param b second point.
     * @return distance.
     */
    private static double distance(final float[] a, final float[] b) {
        double sum = 0;
        for (int i = 0; i < a.length; i++) {
            final double d = a[i] - b[i];
            sum += d * d;
        }
        return Math.sqrt(sum);
    }


    private static Object[] clusterize(float [][] X, int [][] U, float [][] V, int k, double tolerance){
        boolean converged = false;
        long count = 0;

        Object[] resultMatrixes = new Object[2];
        final int nFeatures = X[0].length;

        while (!converged)
        {
            /* Creates an array of @link{MeanElement}s */
            final MeanElement[] newClusterMeans = new MeanElement[k];
            for (int i = 0; i < newClusterMeans.length; i++)
            {
                newClusterMeans[i] = new MeanElement(nFeatures);
            }

            /* computing all U can be redundant for K-Means */
            for (int i = 0; i < X.length; i++)
            {
                final int c = closestCluster(X[i], V, k);
                newClusterMeans[c].add(X[i]);
            }

            /* Check for convergence */
            float distanceSum = 0;
            for (int i = 0; i < k; i++)
            {
                final float[] clusterCenter = V[i];
                final float[] newClusterCenter = newClusterMeans[i].mean();
                //clusterCenters[i] = newClusterCenter;
                distanceSum += distance(clusterCenter, newClusterCenter);
            }

            converged = distanceSum < tolerance;

            for (int i = 0; i < k; i++)
            {
                V[i] = newClusterMeans[i].mean();
            }

            ++count;

            final String message = "k-means iteration " + count + ", cluster error: " + distanceSum;
            System.out.println(message);

        }

        U = computeClusterMembership(X, U, V, k);
        resultMatrixes[0] = U;
        resultMatrixes[1] = V;
        return resultMatrixes;
    }

//    private static Point toPoint(final int offset, final int width) {
//        final int y = offset / width;
//        final int x = offset - y * width;
//        return new Point(x, y);
//    }

    private static final class MeanElement {

        private final double[] sum;
        private int count;


        public MeanElement(final int elementSize) {
            sum = new double[elementSize];
        }


        public void add(final float[] x) {

            if (x.length != sum.length)
            {
                throw new java.lang.IllegalArgumentException("Invalid element size, got " + x.length + ", expecting" + sum.length);
            }

            for (int i = 0; i < x.length; i++)
            {
                sum[i] += x[i];
            }
            ++count;
        }


        public float[] mean() {

            final float[] r = new float[sum.length];
            for (int i = 0; i < r.length; i++)
            {
                r[i] = (float) (sum[i] / count);
            }

            return r;
        }
    }

}
