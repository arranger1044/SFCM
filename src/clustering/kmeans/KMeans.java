/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package clustering.kmeans;

import java.util.HashSet;
import java.util.Random;

/**
 *
 * @author valerio
 */
public class KMeans {

    private KMeans(){
    }

    public static Object[] run (float [][] X, int k, double tolerance,
                                int randomSeed, int initMode, boolean seedEnabled,
                                ClusteringDelegate delegate){

        float [][] V = new float[k][X[0].length];
        int [][] U = new int [X.length][k];

        V = initializeClusterCenterMatrix(X, V, k, seedEnabled, randomSeed,
                                          initMode, delegate);
        U = computeClusterMembership(X, U, V, k);

        return clusterize(X, U, V, k, tolerance, delegate);
    }

    private static float [][] initializeClusterCenterMatrix(float [][] X,
                                                            float [][] V,
                                                            int k,
                                                            boolean seedEnabled,
                                                            int randomSeed,
                                                            int initMode,
                                                            ClusteringDelegate delegate){

        final long startTime = System.currentTimeMillis();

        switch(initMode)
        {
            case 0:
                V = randomInitialization(X, V, k, seedEnabled, randomSeed);
                System.out.println("Random");
                break;
            case 1:
                V = kMeansPlusPlusInitialization(X, V, k, seedEnabled, randomSeed);
                System.out.println("K-Means++");
                break;
            default:
                break;
        }

        final long endTime = System.currentTimeMillis();
        //System.out.println("Initialization done in " + (endTime - startTime) + "ms");
        String message = "Initialization done in " + (endTime - startTime) + "ms";
        delegate.updateStatus(message);
        
        return V;
    }

    private static float[][] randomInitialization(float [][] X, float [][] V, int k,
                                                  boolean seedEnabled, int randomSeed) {

        final Random random = createRandom(seedEnabled, randomSeed);
        final int nbPixels = X.length;
        final int nFeatures = X[0].length;

        final HashSet<Integer> centerLocations = new HashSet<Integer>();
        int clusterCreated = 0;

        while (clusterCreated < k)
        {
            final int clusterCandidate = random.nextInt(nbPixels);

            /* Check if it has already been extracted */
            if (!centerLocations.contains(clusterCandidate))
            {
                /* Let's copy the pixel values in the centroid matrix */
                System.arraycopy(X[clusterCandidate], 0, V[clusterCreated], 0, nFeatures);
                centerLocations.add(clusterCandidate);
                clusterCreated++;
            }
        }
        return V;
    }

    private static float[][] kMeansPlusPlusInitialization(float [][] X, float [][] V, int k,
                                                          boolean seedEnabled, int randomSeed) {

        final Random random = createRandom(seedEnabled, randomSeed);
        final int nbPixels = X.length;
        final int nFeatures = X[0].length;

        // Location of pixels used as cluster centers
        final HashSet<Integer> centerLocation = new HashSet<Integer>();
        int clusterCreated = 0;
        // Choose one center uniformly at random from among pixels
        {
            final int p = random.nextInt(nbPixels);
            centerLocation.add(p);
            System.arraycopy(X[p], 0, V[0], 0, nFeatures);
            clusterCreated++;
        }

        final double[] dp2 = new double[nbPixels];
        while (clusterCreated < k) {
            // For each data point p compute D(p), the distance between p and the nearest center that
            // has already been chosen.
            double sum = 0;
            
            for (int offset = 0; offset < nbPixels; offset++) {  

                if (centerLocation.contains(offset)) {
                    continue;
                }
                // Distance to closest cluster

                final float[] v = X[offset];
                final int cci = closestCluster(v, V, clusterCreated);
                final double d = euclideanDistance(v, V[cci]);
                sum += d * d;
                dp2[offset] = sum;
            }


            /* Add one new data point at random as a new center, using a weighted probability distribution where
             * a point p is chosen with probability proportional to D(p)^2 */
            final double r = random.nextDouble() * sum;
            for (int offset = 0; offset < nbPixels; offset++) {

                /* Test that this is not a repeat of already selected center */
                if (centerLocation.contains(offset)) {
                    continue;
                }

                if (dp2[offset] >= r)
                {
                    centerLocation.add(offset);
                    System.arraycopy(X[offset], 0, V[clusterCreated], 0, nFeatures);
                    clusterCreated++;
                    break;
                }
            }
        }

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

    private static Random createRandom(boolean seedEnabled, int randomSeed) {
//        return config.isRandomizationSeedEnabled()
//                ? new Random(config.getRandomizationSeed())
//                : new Random();
        return (seedEnabled ?
            new Random(randomSeed)
            : new Random());
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
            final double d = euclideanDistanceSqr(clusterCenter, x);
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
    private static double euclideanDistance(final float[] a, final float[] b) {
        double sum = 0;
        for (int i = 0; i < a.length; i++) {
            final double d = a[i] - b[i];
            sum += d * d;
        }
        //return sum;
        return Math.sqrt(sum);
    }

    /**
     * Distance between points <code>a</code> and <code>b</code>.
     * @param a first point.
     * @param b second point.
     * @return distance.
     */
    private static double euclideanDistanceSqr(final float[] a, final float[] b) {
        double sum = 0;
        for (int i = 0; i < a.length; i++) {
            final double d = a[i] - b[i];
            sum += d * d;
        }
        return sum;
        //return Math.sqrt(sum);
    }


    private static Object[] clusterize(float [][] X, int [][] U, float [][] V, 
                                       int k, double tolerance, ClusteringDelegate delegate){
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
                distanceSum += euclideanDistanceSqr(clusterCenter, newClusterCenter);
            }

            converged = distanceSum < tolerance;

            for (int i = 0; i < k; i++)
            {
                V[i] = newClusterMeans[i].mean();
            }

            ++count;

            delegate.updateStatus(null, null, count, distanceSum);

        }

        U = computeClusterMembership(X, U, V, k);
        resultMatrixes[0] = U;
        resultMatrixes[1] = V;
        return resultMatrixes;
    }

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
