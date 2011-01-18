/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package clustering.sfcm;

import java.io.IOException;
import java.util.HashSet;
import java.util.Random;

/**
 * This class encapsulate the Spatial Fuzzy C-Means as proposed in: <i>Chuang,
 * Tzeng, Chen Fuzzy c-means clustering with spatial information for image
 * segmentation. Computerized Medical Imaging and Graphics 30 (2006) 9–15</i>.
 * It implements the Singleton pattern and provides as the only public class
 * method @link{#run}. It also implements the delegation pattern, being free from
 * ImageJ API it can be  more easily reused
 */
public class SFCM {

    /**
     * The class constructor is private in order to prevent the class from being
     * instantiated
     */
    private SFCM(){
    }

    /**
     * The only one method that is visible to the other classes,
     * it is used for start of the algorithm: it inizializaes the matrixes U and V,
     * run the SFCM on them and then pack them in an Object array that returns
     * to the caller.
     * <p>
     * @param X data matrix, a samples x features matrix
     * @param k number of clusters
     * @param tolerance threashold used for checking the convergence of the algorithm
     * @param randomSeed the randomization seed used for generating random numbers
     * @param initMode an integer that specifies the initialization mode for the
     * matrixes <code>U</code> and <code>V</code>.
     * <b>0</b> - firstly initializes the centroid matrix <code>V</code>
     * randomly <i>E. Forgy,
     * “Cluster Analysis of Multivariate Data: Effi- ciency vs. Interpretability
     * of Classification,” Biometrics, vol. 21, pp. 768, 1965<i/>.
     * <b>1</b> - Firstly initializes the centroid matrix <code>V</code> accoding to the
     * <i>K-Means++</i> criterion @link{http://en.wikipedia.org/wiki/K-means++}
     * <b>2</b> - Firstly initializes the cluster membership matrix <code>U</code> randomly
     * @param delegate is the object which is delegate to output partial results
     * of the algorithm.
     * @param m is the fuzzyness parameter
     * @param iterations is the max number of iterations possible for the algorithm
     * @param stopCriterion is an integer representing one of the possible stop criterions:
     * <b>0</b> - it computes the Frobenius norm on the difference between the actual
     * membership matrix <code>U</code> and the previous one from the last iteration
     * <b>1</b> - it computes the Frobenius norm on the difference between the actual
     * centroid matrix <code>V</code> and the previous one from the last iteration
     * <b>2</b> - it computes the max norm on the difference between the actual
     * membership matrix <code>U</code> and the previous one from the last iteration
     * <b>3</b> - it computes the max norm on the difference between the actual
     * centroid matrix <code>V</code> and the previous one from the last iteration
     * @param r is the radius of the window used to compute the spatial information
     * @param p is the weight for the membership values
     * @param q is the weight for the chosen spatial function <code>h</code>
     * @param spatialFunction is an integer representing the choosen spatial function:
     * <b>0</b> - Computes the spatial function as the sum of the membership values
     * for the samples in the neighborhood according to the same cluster
     * <b>1</b> - Computes the spatial function as the number of clusters membership
     * in the neighborhood that are the greatest value for the neighbor according to
     * that cluster
     * @param width the width of the matrix form, as the samples can be viewed when
     * computing the neightborhood for each one
     * @param testing  a boolean value that tells whether the algorithm is run for testing
     * @return a vector of object that has three elements: a defuzzyfied memberhip matrix,
     * the cluster centroid matrix and the cluster membership matrix
     *
     * @see SFCMManager.run()
     * @see ClusteringDelegate
     * @see #defuzzyfyClusterMemberships(float[][])
     */
    public static Object[] run (float [][] X, int k, double tolerance,
                                int randomSeed, int initMode,
                                ClusteringDelegate delegate,
                                double m, long iterations, int stopCriterion,
                                int r, double p, double q, int spatialFunction,
                                int width,
                                boolean testing){

        float [][] V = new float[k][X[0].length];
        float [][] U = new float [X.length][k];

        initializeMatrixes(X, V, U, k, randomSeed, initMode, m);

        Object[] clusteredMatrixes = clusterize(X, U, V, k, tolerance, delegate, 
                                                m, iterations, stopCriterion, r,
                                                p, q, spatialFunction, width,
                                                testing);
        Object[] matrixes = new Object[3];
        matrixes[0] = defuzzyfyClusterMemberships((float [][])clusteredMatrixes[0]);
        matrixes[1] = clusteredMatrixes[1];
        matrixes[2] = clusteredMatrixes[0];
        return matrixes;
    }

    /**
     * This method initialize the matrixes <code>U</code> and <code>V</code> according
     * with the initialization mode specified.
     * @param X data matrix
     * @param V cluster center matrix to be initialized
     * @param U cluster membership matrix to be initialized
     * @param k number of clusters
     * @param randomSeed randomization seed used to generate a random sequence
     * @param initMode an integer that specifies the initialization mode for the
     * matrixes <code>U</code> and <code>V</code>.
     * <b>0</b> - firstly initializes the centroid matrix <code>V</code>
     * randomly <i>E. Forgy,
     * “Cluster Analysis of Multivariate Data: Effi- ciency vs. Interpretability
     * of Classification,” Biometrics, vol. 21, pp. 768, 1965<i/>.
     * <b>1</b> - Firstly initializes the centroid matrix <code>V</code> accoding to the
     * <i>K-Means++</i> criterion @link{http://en.wikipedia.org/wiki/K-means++}
     * <b>2</b> - Firstly initializes the cluster membership matrix <code>U</code> randomly
     * @param delegate is the object which is delegate to output partial results
     * of the algorithm.
     * @param m the fuzziness parameter
     * @return an Object array which containes the initialized <code>U</code>
     * and <code>V</code> matrixes
     * @see #randomInitialization(float[][], float[][], int, int)
     * @see #updateClusterMembershipMatrix(float[][], float[][], float[][], double, double[][])
     * @see #kMeansPlusPlusInitialization(float[][], float[][], int, int)
     * @see #initializeClusterMembershipRandom(float[][], int, int)
     */
    private static Object [] initializeMatrixes(float [][] X, float [][] V, 
                                                float [][] U, int k, 
                                                int randomSeed, int initMode,
                                                double m){
        
        Object[] initedMatrixes = new Object[2];
        double [][] D = null;

        switch(initMode)
        {
            case 0:
                V = randomInitialization(X, V, k, randomSeed);
                D = euclideanDistanceMatrix(X, V, D, m);
                U = updateClusterMembershipMatrix(X, U, V, m, D);
                System.out.println("Random V");
                break;
            case 1:
                V = kMeansPlusPlusInitialization(X, V, k, randomSeed);
                //printMatrix(V);
                D = euclideanDistanceMatrix(X, V, D, m);
                //printMatrix(D);
                U = updateClusterMembershipMatrix(X, U, V, m, D);
                System.out.println("\n\nU\n\n");
                //printMatrix(U);
                System.out.println("K-Means++");
                
                break;
            case 2:
                U = initializeClusterMembershipRandom(U, k, randomSeed);
                System.out.println("Random U");
                break;
            default:
                break;
        }

        initedMatrixes[0] = U;
        initedMatrixes[1] = V;
        return initedMatrixes;
    }


    /**
     * Utility method for printing a float matrix on console
     * @param A the matrix to be printed
     */
    private static void printMatrix(float [][] A){
        for (int i = 0; i < A.length; i++)
        {
            for (int j = 0; j < A[0].length; j ++)
            {
                System.out.print(A[i][j] + " ");
            }
            System.out.println("");
        }
    }

    /**
     * Utility method for printing a double matrix on console
     * @param A the matrix to be printed
     */
    private static void printMatrix(double [][] A){
        for (int i = 0; i < A.length; i++)
        {
            for (int j = 0; j < A[0].length; j ++)
            {
                System.out.print(A[i][j] + " ");
            }
            System.out.println("");
        }
    }

    /**
     * Initializes the cluster membership matrix <code>U</code> randomly by generating
     * random values in [0, 1] and normalizing them in order to have the sum of
     * membership values for a pixel equal one.
     * @param U cluster membership matrix to be initialized
     * @param k number of clusters
     * @param randomSeed randomization seed used to generate a random sequence
     * @return the initialized cluster membership matrix
     * @see #initializeMatrixes(float[][], float[][], float[][], int, int, int, double)
     */
    private static float [][] initializeClusterMembershipRandom(float [][] U, int k, int randomSeed){

        int nClusters = U[0].length;
        final Random random = createRandom(randomSeed);
        for (int i = 0; i < U.length; i++)
        {
            float sum = 0;
            for (int j = 0; j < nClusters; j++)
            {
                U[i][j] =  random.nextFloat();
                sum += U[i][j];
            }
            for (int j = 0; j < nClusters; j++)
            {
                U[i][j] /= sum;
            }
        }
        //printMatrix(U);
        return U;
    }

    /**
     * Initializes the cluster center matrix<code>V</code> randomly by selecting
     * sample from X and checking not to select an already selected sample
     * @param X data matrix
     * @param V cluster centroid matrix to be initialized
     * @param k number of clusters
     * @param randomSeed randomization seed used to generate a random sequence
     * @return the initialized cluster matrix
     * @see #initializeMatrixes(float[][], float[][], float[][], int, int, int, double)
     */
    private static float[][] randomInitialization(float [][] X, float [][] V, int k, int randomSeed) {

        final Random random = createRandom(randomSeed);
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

    /**
     * Initializes the cluster center matrix<code>V</code> according to the
     * <i>K-means++</i> initialization criterion. The first centroid is selected
     * randomly from the samples and the remaining ones are selected according
     * to a weighted probability distribution on the distance from the samples and
     * the already selected centroids.
     * @link{http://en.wikipedia.org/wiki/K-means++}
     * @param X the data matrix
     * @param V the centroid matrix to be initialized
     * @param k the number of clusters
     * @param randomSeed randomization seed used to generate a random sequence
     * @return the initialized cluster center matrix
     * @see #initializeMatrixes(float[][], float[][], float[][], int, int, int, double)
     */
    private static float[][] kMeansPlusPlusInitialization(float [][] X, float [][] V, int k, int randomSeed) {

        final Random random = createRandom(randomSeed);
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


            /* Add one new data point at random as a new center, using a weighted
             * probability distribution where
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

    /**
     * Creates a <i>Random</i> instance using a seed
     * @param randomSeed an integer used to generate the Random instance
     * @return the created instance
     * @see #initializeMatrixes(float[][], float[][], float[][], int, int, int, double)
     */
    private static Random createRandom(int randomSeed) {
//        return config.isRandomizationSeedEnabled()
//                ? new Random(config.getRandomizationSeed())
//                : new Random();
        return new Random(randomSeed);
    }

    /**
     * Returns the index of the closest cluster to point <code>x</code>.
     * Called in the K-Means++ initialization of the cluster centers matrix
     * @param x  a sample as an array of float.
     * @param clusterCenters cluster centers matrix.
     * @param k number of clusters
     * @return index of the closest cluster
     * @see #kMeansPlusPlusInitialization(float[][], float[][], int, int)
     */
    private static int closestCluster(final float[] x, final float[][] clusterCenters, int k) {

        double minDistance = Double.MAX_VALUE;
        int closestCluster = -1;

        for (int i = 0; i < k; i++)
        {
            final float[] clusterCenter = clusterCenters[i];
            final double d = euclideanDistance(clusterCenter, x);
            if (d < minDistance)
            {
                minDistance = d;
                closestCluster = i;
            }
        }

        return closestCluster;
    }


    /**
     * Euclidean distance  computed between points <code>a</code> and <code>b</code>.
     * This is the squared rooted version
     * @param a first point.
     * @param b second point.
     * @return the euclidean distance computed.
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
     * Computes a distance matrix starting from two sample x features matrixes
     * It is used to precompute all the distances between the samples and the centroids
     * and it is used to update the cluster membership matrix <code>U</code>:
     * The distance matrix holds one more column used to precompute the sum of the
     * inverse exponentiated distances between all the centroids and a sample. In this way
     * the update rule for the <code>U</code> matrix is more efficiently computed
     * @param A sample x feature matrix
     * @param B sample x feature matrix
     * @param D distance matrix, if null it will be instantiated
     * @param m fuzzyness parameter
     * @return the instantiated distance matrix
     * @see #updateClusterMembershipMatrix(float[][], float[][], float[][], double, double[][])
     */
    private static double [][] euclideanDistanceMatrix(final float [][] A, final float [][] B,
                                                            double [][] D, double m){

        if (D == null)
        {
            System.out.println("Distance Matrix allocated");
            D = new double[A.length][B.length + 1];
        }

        double exp = (1 / (m - 1));
        int nDims = A[0].length;

        for (int i = 0; i < A.length; i++)
        {
            D[i][B.length] = 0.0;
            for (int j = 0; j < B.length; j++)
            {
                double sum = 0;
                for (int k = 0; k < nDims; k++)
                {
                    final double d = A[i][k] - B[j][k];
                    sum += d * d;
                }

                //D[i][j] = (float)Math.sqrt(sum);
                D[i][j] = sum;
                //System.out.println("D"+ i + j + " " + D[i][j]);
                if (D[i][j] == 0.0)
                {
                    D[i][B.length] = -1;
                }
                else
                {
                    if (D[i][B.length] != -1)
                    {
                        D[i][B.length] += Math.pow(1 / D[i][j], exp);
                    }
                }
            }
        }
        return D;
    }

    /**
     * Computes the exponential memberhip matrix <code>U</code>.
     * It stores the result in a matrix <code>Um</code> conformed to <code>U</code>.
     * If the objective function shall be not computed anymore this method can be
     * deprecated and the the exponentiation computed when necessary in the update
     * rule for the centroid matrix <code>V</code>
     * @param U cluster membership matrix
     * @param Um exponentiated cluster membership matrix, if null it will be instantiated
     * @param m the fuzzyness parameter
     * @return the exponentiated Um matrix
     * @see #updateClusterCenterMatrix(float[][], float[][], float[][])
     */
    private static float [][] computeExponentialMembership(float [][] U, float [][] Um, double m){

        int nClusters = U[0].length;
        if (Um == null)
        {
            System.out.println("Allocating U^m matrix");
            Um = new float[U.length][nClusters];
        }

        for (int i = 0; i < U.length; i++)
        {
            for (int j = 0; j < nClusters; j++)
            {
                Um[i][j] = (float) Math.pow(U[i][j], m);
                if ( Float.isNaN(Um[i][j]))
                {
                    throw new IllegalArgumentException("Nan found "+ U[i][j] + " m " + m);
                }
            }
        }
        return Um;
    }

    /**
     * Computes the objective function <b>J</b> as proposed in the <i>Fuzzy C-Means</i>
     * algorithm. It is used as a possible criterion to check the convergence in the
     * testing mode (when the difference between the current value and the previous
     * iteration one is sent to the delegate)
     * @param X the data matrix.
     * @param Um the exponentiated cluster membership matrix.
     * @param V the cluster center matrix.
     * @param D the distance matrix.
     * @return the value of objective function
     * @see #clusterize(float[][], float[][], float[][], int,
     * double, clustering.sfcm.ClusteringDelegate, double, long,
     * int, int, double, double, int, int, boolean)
     */
    private static float computeObjectiveFunction(float [][] X, float [][] Um, float [][]V, double [][] D){

        float objF = 0;
        for(int i = 0; i < X.length; i++)
        {
             for(int j = 0; j < V.length; j++)
               {
                //float distancePixelToCluster = euclideanDistance(i, j, D);
                double distancePixelToCluster = D[i][j];
                //System.out.println("Distance " + distancePixelToCluster);
                objF += distancePixelToCluster * Um[i][j];
                //System.out.println("obj " + objF);
               }
        }
        return objF;
    }

    /**
     * Updates the cluster membership matrix <code>U</code> according to the
     * <i>Fuzzy C-Means</i> algorithm.
     * Each value of <code>U</code> is updated as the division between the sum of
     * the exponentiated inverse distances between the sample and all the centroids
     * and the exponentiated distance between the sample and the current centroid
     * It uses the precomputed matrix D in order to speed the computation up.
     * @param X the data matrix
     * @param U the cluster membership matrix to be updated
     * @param V the cluster centroid matrix
     * @param m the fuzzyness parameter used in the exponentiation
     * @param D the distance matrix with an additional column with precomputed values
     * @return the updated  cluster membership matrix
     * @see #clusterize(float[][], float[][], float[][], int, double,
     * clustering.sfcm.ClusteringDelegate, double, long, int, int, double,
     * double, int, int, boolean)
     * @see #euclideanDistanceMatrix(float[][], float[][], double[][], double)
     */
    private static float[][] updateClusterMembershipMatrix(float [][] X, float [][] U, float [][] V,
                                                           double m, double [][] D){

        int nClusters = V.length;
        double exp = (1.0 / (m - 1.0));
        for(int i = 0; i < X.length; i++)
        {
            int count = 0;

            for(int j = 0; j < nClusters; j++)
            {
                    //float num = euclideanDistance(i, j, D);
                    double distance = D[i][j];
                    if(distance != 0.0)
                    {
//                                            //  sum of distances from this data point to all clusters.
//                        float sumTerms = 0;
//                        for(int k = 0; k < V.length; k++)
//                        {
//                            //float thisDistance = euclideanDistance(i,k, D);
//                            float thisDistance = D[i][k];
//
//                            //sumTerms += Math.pow(num / thisDistance, (2f / (m - 1f)));
//                            sumTerms += Math.pow(num / thisDistance, (1f / (m - 1f)));
////                            if (thisDistance == 0.0f)
////                            {
////                                System.out.println("Ouch " + i + " " + k + " " + j + " " +
////                                        Math.pow(num / thisDistance, (2f / (m - 1f))));
////                            }
//
////                                             if ( Float.isNaN(thisDistance))
////                    {
////                        throw new IllegalArgumentException("thisDistance"+ num + " m " + m);
////                    }
////                                             if ( Float.isNaN(sumTerms))
////                    {
////                        throw new IllegalArgumentException("sumterms "+ thisDistance + "
////                            d " + num+ " j " + j + " k " + k);
////                    }
//                        }
//
//                        //sumTerms = (float)(Math.pow(num, (1f / (m - 1f)))) / D[i][V.length];
//                        //System.out.println(D[i][V.length] + " " + sumTerms + " " + (float)(Math.pow(num, (1f / (m - 1f)))));
//                        U[i][j] = (1f / sumTerms);
                        double denominatorSum = D[i][nClusters];
                        if (denominatorSum != -1)
                        {
                            
                            double numerator = Math.pow(distance, exp);
                            //U[i][j] = 1 / (numerator / denominator);
                            U[i][j] = (float) (1f / (numerator * denominatorSum));
                            if ( Float.isNaN(U[i][j]))
                            {
                                throw new IllegalArgumentException("Nan found " + " m " + m + " num " + numerator
                                        + " den " + denominatorSum + " dis " + distance + " exp " + exp);
                            }
                        }
                        else
                        {
                            U[i][j] = 0;
                        }
                    }
                    else
                    {
                        count++;
//                        float sum = 0;
//                        for (int h = 0; h < V.length; h++)
//                        {
//                            if (h != j)
//                                sum += U[i][h];
//                        }
//                        System.out.println("Sum " + sum);
                        U[i][j] = 1.0f;
                    }

                    if ( Float.isNaN(U[i][j]))
                    {
                        throw new IllegalArgumentException("Nan found " + " m " + m);
                    }
            }

            if (count > 1)
            {
                for(int j = 0; j < V.length; j++)
                {
                    if (U[i][j] == 1.0f)
                    {
                        U[i][j] = 1f / count;
                    }
                }
            }
        }
        
        return U;

    }

    /**
     * Updates the cluster center matrix <code>V</code> according to the <i>Fuzzy
     * C-Means</i> algorithm.
     * Each element of  <code>V</code> is updated considering the sum of the samples
     * multiplied by their exponentiated membership and normalized against the
     * sum of the exponentiated memberships.
     * It uses the already exponentiated matrix <code>Um</code>
     * @param X data matrix
     * @param Um exponentiated cluster membership matrix
     * @param V cluster center matrix to update
     * @return the updated cluster center matrix
     * @see #clusterize(float[][], float[][], float[][], int, double,
     * clustering.sfcm.ClusteringDelegate, double, long, int, int, double,
     * double, int, int, boolean)
     * @see #computeExponentialMembership(float[][], float[][], double)
     */
    private static float[][] updateClusterCenterMatrix(float [][] X, float [][] Um,
                                                       float [][] V){


        double numerator = 0;
        double denominator = 0;

        int nClusters = V.length;
        int nFeatures = X[0].length;

        for (int i = 0; i < nClusters; i++)
        {

            for (int j = 0; j < nFeatures; j++)
            {
                numerator = 0;
                denominator = 0;
                for(int k = 0; k < X.length; k++)
                {
                    /// trovare tutti i pixel di 1 determinato cluster k
                    float x[] = X[k];
                    numerator += Um[k][i] * x[j];
                    denominator += Um[k][i];
                    //System.out.println("Um" + k + i + " " + Um[k][i] + " x " + x[j]);
                    //System.out.println("N " + numerator + " D " + denominator);
                }
                /// cluster x caratteristiche. ovvero V
                V[i][j] = (float) (numerator / denominator);
                if ( Float.isNaN(V[i][j]))
                {
                    for(int k = 0; k < X.length; k++)
                    {
                        //System.out.println(Um[k][i]);
                    }
                    throw new IllegalArgumentException("nume "+ numerator + " denom "
                            + denominator + "i " + i);
                }
                //System.out.println("V" + i + j + "n " + numerator + " d " + denominator  );
            }
        }
        return V;
    }

    /**
     * Computes the defuzzyfication of the cluster membership matrix <code>U</code>
     * by using a max operator on the membership values of a sample against all the
     * clusters.
     * It instantiates a new matrix of integers containing only 1 or 0 (hard partition)
     * @param U the cluster membership matrix
     * @return the hard partition matrix
     * @see #run(float[][], int, double, int, int,
     * clustering.sfcm.ClusteringDelegate, double, long, int, int, double,
     * double, int, int, boolean)
     */
    private static int [][] defuzzyfyClusterMemberships(float [][] U){

        int nClusters = U[0].length;
        int [][] defU = new int [U.length][nClusters];

        for (int i = 0; i < U.length; i++)
        {

            float max = -1;
            int maxPos = -1;
            for(int j = 0; j < nClusters; j++)
            {
                if ( U[i][j] > max)
                {
                    max = U[i][j];
                    maxPos = j;
                }
            }
            defU[i][maxPos] = 1;
        }
        
        return defU;
    }

    /**
     * Implements a not shallow copy between two matrixes. It is used to save the
     * previous iteration matrixes <code>U</code> and <code>V</code>
     * @param A the source matrix
     * @param B the destination matrix
     * @see System.arraycopy
     * @see #clusterize(float[][], float[][], float[][], int, double,
     * clustering.sfcm.ClusteringDelegate, double, long, int, int, double,
     * double, int, int, boolean)
     */
    private static void matrixCopy(float [][] A, float [][] B){

        int nDim = A[0].length;
        for (int i = 0; i < A.length; ++i)
        {
            System.arraycopy(A[i], 0, B[i], 0, nDim);
        }
    }

    /**
     * Computes the Frobenius norm of the difference matrix between two input matrixes.
     * They must have the same dimensions.
     * It is used as a possible criterion to check convergence.
     * @param A a float matrix
     * @param B a float matrix
     * @return the computed  Frobenius norm
     * @see #checkConvergence(float[][], float[][], float[][], float[][], int)
     */
    private static float frobeniusNorm(float [][] A, float [][] B){
        
        float d = 0;
        float sum = 0;
        int nDim = A[0].length;
        for(int i = 0; i < A.length; i++)
        {
            for(int j = 0; j < nDim; j++)
            {
                 d = A[i][j] - B[i][j];
                 sum += d*d;
            }
        }
        return (float) Math.sqrt(sum);
    }


    /**
     * Computes the max norm of the difference matrix between two input matrixes.
     * They must have the same dimensions.
     * It is used as a possible criterion to check convergence.
     * @param A a float matrix
     * @param B a float matrix
     * @return the computed max norm
     * @see #checkConvergence(float[][], float[][], float[][], float[][], int)
     */
    private static float maxNorm(float [][] A, float [][] B){

        float max = 0;
        float abs = 0;
        int nDim = A[0].length;
        for(int i = 0; i < A.length; i++)
        {
            for(int j = 0; j < nDim; j++)
            {
                abs = Math.abs(A[i][j] - B[i][j]);
                if (abs > max)
                {
                    max = abs;
                }
            }
        }
        return max;
    }

    /**
     * Checks for the convergence of the algorithm. The stopping criterion is
     * specified in input.
     * @param U the cluster membership matrix
     * @param oldU the previous iteration cluster membership matrix
     * @param V the cluster center matrix
     * @param oldV the previous iteration cluster center matrix
     * @param criterion is an integer representing one of the possible stop criterions:
     * <b>0</b> - it computes the Frobenius norm on the difference between the actual
     * membership matrix <code>U</code> and the previous one
     * <b>1</b> - it computes the Frobenius norm on the difference between the actual
     * centroid matrix <code>V</code> and the previous one
     * <b>2</b> - it computes the max norm on the difference between the actual
     * membership matrix <code>U</code> and the previous one
     * <b>3</b> - it computes the max norm on the difference between the actual
     * centroid matrix <code>V</code> and the previous one
     * @return the computed difference according to the chosen criterion
     * @see #clusterize(float[][], float[][], float[][], int, double,
     * clustering.sfcm.ClusteringDelegate, double, long, int, int, double,
     * double, int, int, boolean)
     * @see #frobeniusNorm(float[][], float[][])
     * @see #maxNorm(float[][], float[][]) 
     */
    private static float checkConvergence(float [][] U, float [][] oldU,
                                          float [][] V, float [][] oldV, int criterion){
        float diff = 0;

        switch(criterion)
        {
            case 0:
                diff = frobeniusNorm(U, oldU);
                break;
            case 1:
                diff = frobeniusNorm(V, oldV);
                break;
            case 2:
                diff = maxNorm(U, oldU);
                break;
            case 3:
                diff = maxNorm(V, oldV);
                break;
            default:
                break;
        }
        return diff;
    }

    /**
     * Updates the cluster membership matrix <code>U</code> according to the
     * <i>Spatial Fuzzy C-Means</i> algorithm implemented. It updates each element of
     * <code>U</code> by weightening it exponentiating it and multiplying it by the
     * exponentiation of a spatial function which in determines the weight of the
     * neighborhood for the current sample.
     * In order to consider a neighborhood for a sample we reinterpret the sample
     * set as a matrix and consider a neighbourhood a window of a certain radius
     * centered on the sample (which is a matrix element).
     * Different kinds of spatial function are proposed and have to be specified as
     * input parameters.
     * @param U the cluster membership matrix to update
     * @param r the radius of the window representing the neighbourhood for a sample
     * @param p the exponent used to weight the membership values
     * @param q the exponent used to weight the spatial function
     * @param spatialFunctionType an integer representing the choosen spatial function:
     * <b>0</b> - Computes the spatial function as the sum of the membership values
     * for the samples in the neighborhood according to the same cluster
     * <b>1</b> - Computes the spatial function as the number of clusters membership
     * in the neighborhood that are the greatest value for the neighbor according to
     * that cluster
     * @param cols the number of column in the matrix representation of the samples
     * @param xy an integer matrix used to precompute the matrix indices of the matrix form
     * @param offset an integer matrix used to precompute the array indices of
     * the array form of the sample
     * @param uPhQ a support matrix used to store the computed value of the update
     * @return the updated cluster membership matrix
     * @see #clusterize(float[][], float[][], float[][], int, double,
     * clustering.sfcm.ClusteringDelegate, double, long, int, int, double,
     * double, int, int, boolean)
     */
    private static float [][] updateMembershipsWithSpatialInformation(float [][] U, int r,
                                                                      double p, double q,
                                                                      int spatialFunctionType,
                                                                      int cols, int [][] xy,
                                                                      int [][] offset,
                                                                      double [][] uPhQ){

        //float[][] Uupdate = new float[U.length][U[0].length];
        int nClusters = U[0].length;
        if(uPhQ == null)
        {
            uPhQ = new double[U.length][nClusters + 1];
        }

        int rows = U.length / cols;

        for(int i = 0; i < U.length; i++)
        {
            //double numerator = 0;
            uPhQ[i][nClusters] = 0;

            for(int j = 0; j < nClusters; j++)
            {

                int row = xy[i][0];
                int col = xy[i][1];
                double h = 0;

                for(int ry = -r; ry <= r; ry++)
                {
                    int y = row + ry;
                    if (y >= 0 && y < rows)
                    {
                        for (int rx = -r; rx <= r; rx++)
                        {
                            int x = col + rx;
                            if (x >= 0 && x < cols)
                            {
                                int elem = offset[y][x];

                                if(spatialFunctionType == 0)
                                {
                                    h += U[elem][j];
                                }
                                else
                                {
                                    boolean exit = false;
                                    for(int k = 0; k < nClusters && !exit ; k++)
                                    {
                                        if(U[elem][j] < U[elem][k])
                                        {
                                            exit = true;
                                        }
                                    }

                                    if(exit == false)
                                    {
                                        h += 1;
                                    }
                                }
                            }
                        }
                    }
                }

//                if (h == 0.0)
//                    System.out.println("ciautistico");
                double uPTimeshQ = Math.pow(U[i][j], p) *
                                   Math.pow(h, q);
                uPhQ[i][j] = uPTimeshQ;
                uPhQ[i][nClusters] += uPTimeshQ;

//                int row = i / cols;
//                int col = i - row * cols;
//
//                numerator = (float) ((Math.pow(U[i][j], p)) *
//                        Math.pow(functionH(U, r, row, col, j, rows, cols, spatialFunctionType), q));
//
//                float denominator = 0;
//
//                for (int k = 0; k < U[0].length; k++)
//                {
//                    denominator += ((Math.pow(U[i][k], p)) *
//                            (Math.pow(functionH(U, r, row, col, k, rows, cols, spatialFunctionType), q)));
//                }
//
//                Uupdate[i][j] = numerator / denominator;
            }
        }

        for(int i = 0; i < U.length; i++)
        {
            double denominatorSum = uPhQ[i][nClusters];
            for(int j = 0; j < nClusters; j++)
            {
                U[i][j] = (float) (uPhQ[i][j] / denominatorSum);
            }
        }
        //U = Uupdate;
        return U;
    }

    private static float [][] opt1UpdateMembershipsWithSpatialFunctionH(float [][] U, int rad,
                                                                          double p, double q,
                                                                          int cols, double [] supC,
                                                                          double [][] uPhQ){
        int nClusters = U[0].length;
        if(uPhQ == null)
        {
            uPhQ = new double[U.length][nClusters + 1];
        }
        if (supC == null)
        {
            supC = new double[cols];
        }

        int rows = U.length / cols;

        for(int i = 0; i < U.length; i++)
        {
            uPhQ[i][nClusters] = 0;
        }

        double uPTimeshQ  = 0;
        for(int k = 0; k < nClusters; k++)
        {
            for (int i = 0; i < rows; ++i)
            {
                double h = 0;
                int qFirst = 0;
                int qLast = 0;
                for (int rx = 0; rx <= rad; rx++)
                {
                    if (rx < cols)
                    {
                        double mem = 0;
                        for (int ry =- rad; ry <= rad; ry++)
                        {
                            int y = i + ry;

                            if (y >= 0 && y < rows)
                            {
                                //mem += mat[y][rx];
                                mem += U[y * cols + rx][k];
                            }
                        }
                        supC[qLast] = mem;
                        qLast++;
                        h += mem;
                    }
                }
                int j = 0, col, y;
                //System.out.println("i " + i + " j " + j + "h " + h);
                uPTimeshQ = Math.pow(U[i * cols][k], p) *
                            Math.pow(h, q);
                uPhQ[i * cols][k] = uPTimeshQ;
                uPhQ[i * cols][nClusters] += uPTimeshQ;

                qFirst = 0;
                qLast = rad + 1;
                for (j = 1; j < cols; j++)
                {
                    col = j + rad;
                    if (j - rad > 0)
                    {
                        h -= supC[qFirst];
                        qFirst++;
                    }
                    if (j + rad < cols)
                    {
                        double mem = 0;
                        for (int ry =- rad; ry <= rad; ry++)
                        {
                            y = i + ry;

                            if (y >= 0 && y < rows)
                            {
                                //mem += mat[y][col];
                                mem += U[y * cols + col][k];

                            }
                        }
                        supC[qLast] = mem;
                        qLast++;
                        h += mem;

                    }
                    //System.out.println("i " + i + " j " + j + "h " + h);
                    int offset = i * cols + j;
                    uPTimeshQ = Math.pow(U[offset][k], p) *
                                Math.pow(h, q);
                    uPhQ[offset][k] = uPTimeshQ;
                    uPhQ[offset][nClusters] += uPTimeshQ;
                }
            }
        }
        
        for(int i = 0; i < U.length; i++)
        {
            double denominatorSum = uPhQ[i][nClusters];
            for(int j = 0; j < nClusters; j++)
            {
                U[i][j] = (float) (uPhQ[i][j] / denominatorSum);
            }
        }
        return U;
    }

    private static float [][] opt1UpdateMembershipsWithSpatialFunctionG(float [][] U, int rad,
                                                                        double p, double q,
                                                                        int cols, double [] supC,
                                                                        double [][] uPhQ,
                                                                        short [][] G){
        int nClusters = U[0].length;
        if(uPhQ == null)
        {
            uPhQ = new double[U.length][nClusters + 1];
        }
        if (supC == null)
        {
            supC = new double[cols];
        }
        if (G == null)
        {
            G = new short[U.length][nClusters];
        }

        /* Computing the G matrix just once */
        for (int i = 0; i < G.length; ++i)
        {
            double max = U[i][0];
            int maxPos = 0;
            for (int j = 1; j < nClusters; ++j)
            {
                G[i][j] = 0;
                if(U[i][j] > max)
                {
                    max = U[i][j];
                    maxPos = j;
                }
            }
            G[i][maxPos] = 1;
        }
        int rows = U.length / cols;

        for(int i = 0; i < U.length; i++)
        {
            uPhQ[i][nClusters] = 0;
        }

        double uPTimeshQ  = 0;
        for(int k = 0; k < nClusters; k++)
        {
            for (int i = 0; i < rows; ++i)
            {
                double h = 0;
                int qFirst = 0;
                int qLast = 0;
                for (int rx = 0; rx <= rad; rx++)
                {
                    if (rx < cols)
                    {
                        double mem = 0;
                        for (int ry =- rad; ry <= rad; ry++)
                        {
                            int y = i + ry;

                            if (y >= 0 && y < rows)
                            {
                                //mem += mat[y][rx];
                                mem += G[y * cols + rx][k];
                            }
                        }
                        supC[qLast] = mem;
                        qLast++;
                        h += mem;
                    }
                }
                int j = 0, col, y;
                //System.out.println("i " + i + " j " + j + "h " + h);
                uPTimeshQ = Math.pow(U[i * cols][k], p) *
                            Math.pow(h, q);
                uPhQ[i * cols][k] = uPTimeshQ;
                uPhQ[i * cols][nClusters] += uPTimeshQ;

                qFirst = 0;
                qLast = rad + 1;
                for (j = 1; j < cols; j++)
                {
                    col = j + rad;
                    if (j - rad > 0)
                    {
                        h -= supC[qFirst];
                        qFirst++;
                    }
                    if (j + rad < cols)
                    {
                        double mem = 0;
                        for (int ry =- rad; ry <= rad; ry++)
                        {
                            y = i + ry;

                            if (y >= 0 && y < rows)
                            {
                                //mem += mat[y][col];
                                mem += G[y * cols + col][k];

                            }
                        }
                        supC[qLast] = mem;
                        qLast++;
                        h += mem;

                    }
                    //System.out.println("i " + i + " j " + j + "h " + h);
                    int offset = i * cols + j;
                    uPTimeshQ = Math.pow(U[offset][k], p) *
                                Math.pow(h, q);
                    uPhQ[offset][k] = uPTimeshQ;
                    uPhQ[offset][nClusters] += uPTimeshQ;
                }
            }
        }

        for(int i = 0; i < U.length; i++)
        {
            double denominatorSum = uPhQ[i][nClusters];
            for(int j = 0; j < nClusters; j++)
            {
                U[i][j] = (float) (uPhQ[i][j] / denominatorSum);
            }
        }
        return U;
    }

    private static float [][] opt2UpdateMembershipsWithSpatialFunctionH(float [][] U, int rad,
                                                                          double p, double q,
                                                                          int cols, double [][] sumU,
                                                                          double [][] uPhQ){
        int nClusters = U[0].length;

        int rows = U.length / cols;
        
        if(uPhQ == null)
        {
            uPhQ = new double[U.length][nClusters + 1];
        }
        if (sumU == null)
        {
            sumU = new double[U.length][nClusters];
        }

        for(int i = 0; i < U.length; i++)
        {
            uPhQ[i][nClusters] = 0;
        }

        int height = rows - 1;
        int width = cols - 1;
        int yMin, yMax, xMin, xMax, y1w, y2w, yw, ywx = 0;
        double s1, s2, s3, s4, tot;
        /* Building up the integral image in sumU*/
        
        for (int k = 0; k < nClusters; ++k)
        {
            double s = 0;
            for (int x = 0; x < cols; ++x)
            {
                s += U[x][k];
                sumU[x][k] = s;
            }

            for (int y = 1; y < rows; ++y )
            {
                yw = y * cols;
                sumU[yw][k] = sumU[yw - cols][k] + U[yw][k];
                for (int x = 1; x < cols; ++x )
                {
                    ywx = yw + x;
                    sumU[ywx][k] = sumU[ywx - cols][k] + sumU[ywx - 1][k]
                                 + U[ywx][k] - sumU[ywx - cols - 1][k];
                }
            }

            for ( int y = 0; y <= height; ++y )
            {
                yMin = Math.max(-1, y - rad - 1);
                yMax = Math.min(height, y + rad);
                y1w = yMin * cols;
                y2w = yMax * cols;
                for (int x = 0; x <= width; ++x)
                {
                    xMin = Math.max(-1, x - rad - 1);
                    xMax = Math.min(width, x + rad);
                    if (y1w < 0 && xMin < 0)
                    {
                        s1 = 0;
                        s2 = 0;
                        s3 = 0;
                    }
                    else if (y1w < 0)
                    {
                        s1 = 0;
                        s2 = 0;
                        s3 = sumU[y2w + xMin][k];

                    }
                    else if (xMin < 0)
                    {
                        s1 = 0;
                        s3 = 0;
                        s2 = sumU[y1w + xMax][k];
                    }
                    else
                    {
                        s1 = sumU[y1w + xMin][k];
                        s2 = sumU[y1w + xMax][k];
                        s3 = sumU[y2w + xMin][k];
                    }

                    s4 = sumU[y2w + xMax][k];
                    tot = s1 + s4 - s2 - s3;
                    int offset = y * cols + x;
                    double uPTimeshQ = Math.pow(U[offset][k], p) *
                                       Math.pow(tot, q);
                    uPhQ[offset][k] = uPTimeshQ;
                    uPhQ[offset][nClusters] += uPTimeshQ;
//                        System.out.print(s1 + "+" + s4 + "-" +
//                                        s2 + "-" + s3 + "=");
//                        System.out.println( tot);
                }
            }
        }
//        for (int k = 0; k < nClusters; ++k)
//        {
//
//        }


        for(int i = 0; i < U.length; i++)
        {
            double denominatorSum = uPhQ[i][nClusters];
            for(int j = 0; j < nClusters; j++)
            {
                U[i][j] = (float) (uPhQ[i][j] / denominatorSum);
            }
        }
        return U;
    }


    private static float [][] opt2UpdateMembershipsWithSpatialFunctionG(float [][] U, int rad,
                                                                          double p, double q,
                                                                          int cols, int [][] sumU,
                                                                          double [][] uPhQ,
                                                                          short [][] G){
        int nClusters = U[0].length;

        int rows = U.length / cols;

        if(uPhQ == null)
        {
            uPhQ = new double[U.length][nClusters + 1];
        }
        if (sumU == null)
        {
            sumU = new int[U.length][nClusters];
        }

        if (G == null)
        {
            G = new short[U.length][nClusters];
        }

        /* Computing the G matrix just once */
        for (int i = 0; i < G.length; ++i)
        {
            double max = U[i][0];
            int maxPos = 0;
            for (int j = 1; j < nClusters; ++j)
            {
                G[i][j] = 0;
                if(U[i][j] > max)
                {
                    max = U[i][j];
                    maxPos = j;
                }
            }
            G[i][maxPos] = 1;
        }

        for(int i = 0; i < U.length; i++)
        {
            uPhQ[i][nClusters] = 0;
        }

        int height = rows - 1;
        int width = cols - 1;
        int yMin, yMax, xMin, xMax, y1w, y2w, yw, ywx = 0;
        int s1, s2, s3, s4, tot;
        /* Building up the integral image in sumU*/

        for (int k = 0; k < nClusters; ++k)
        {
            int s = 0;
            for (int x = 0; x < cols; ++x)
            {
                s += G[x][k];
                sumU[x][k] = s;
            }

            for (int y = 1; y < rows; ++y )
            {
                yw = y * cols;
                sumU[yw][k] = sumU[yw - cols][k] + G[yw][k];
                for (int x = 1; x < cols; ++x )
                {
                    ywx = yw + x;
                    sumU[ywx][k] = sumU[ywx - cols][k] + sumU[ywx - 1][k]
                                 + G[ywx][k] - sumU[ywx - cols - 1][k];
                }
            }

            for ( int y = 0; y <= height; ++y )
            {
                yMin = Math.max(-1, y - rad - 1);
                yMax = Math.min(height, y + rad);
                y1w = yMin * cols;
                y2w = yMax * cols;
                for (int x = 0; x <= width; ++x)
                {
                    xMin = Math.max(-1, x - rad - 1);
                    xMax = Math.min(width, x + rad);
                    if (y1w < 0 && xMin < 0)
                    {
                        s1 = 0;
                        s2 = 0;
                        s3 = 0;
                    }
                    else if (y1w < 0)
                    {
                        s1 = 0;
                        s2 = 0;
                        s3 = sumU[y2w + xMin][k];

                    }
                    else if (xMin < 0)
                    {
                        s1 = 0;
                        s3 = 0;
                        s2 = sumU[y1w + xMax][k];
                    }
                    else
                    {
                        s1 = sumU[y1w + xMin][k];
                        s2 = sumU[y1w + xMax][k];
                        s3 = sumU[y2w + xMin][k];
                    }

                    s4 = sumU[y2w + xMax][k];
                    tot = s1 + s4 - s2 - s3;
                    int offset = y * cols + x;
                    double uPTimeshQ = Math.pow(U[offset][k], p) *
                                       Math.pow(tot, q);
                    uPhQ[offset][k] = uPTimeshQ;
                    uPhQ[offset][nClusters] += uPTimeshQ;
                }
            }
        }

        for(int i = 0; i < U.length; i++)
        {
            double denominatorSum = uPhQ[i][nClusters];
            for(int j = 0; j < nClusters; j++)
            {
                U[i][j] = (float) (uPhQ[i][j] / denominatorSum);
            }
        }
        return U;
    }

    /**
     * Precomputes in a samples x 2 matrix the matrix form indices associated
     * to each array form indices.
     * @param offset the number of samples
     * @param width the number of columns of the matrix form
     * @return the instantiated precomputed matrix
     * @see #updateMembershipsWithSpatialInformation(float[][], int, double,
     * double, int, int, int[][], int[][], double[][])
     */
    private static int [][] precomputeMatrixForm(int offset, int width){
        int [][] matrixForm = new int [offset][2];
        for (int i = 0; i < offset; i++)
        {
            matrixForm[i][0] = i / width;
            matrixForm[i][1] = i - matrixForm[i][0] * width;
        }
        return matrixForm;
    }

    /**
     * Precomputes in a rows x columns matrix the array form indices associated
     * to each pair of matrix form indices.
     * @param rows the number of rows in matrix form
     * @param columns the number of columns in matrix form
     * @return the instantiated and precomputed matrix
     * @see #updateMembershipsWithSpatialInformation(float[][], int, double,
     * double, int, int, int[][], int[][], double[][])
     */
    private static int [][] precomputeArrayForm(int rows, int columns){
        int [][] arrayForm = new int [rows][columns];
        for(int i = 0; i < rows; i++)
        {
            for (int j = 0; j < columns; j++)
            {
                arrayForm[i][j] = j + i * columns;
            }
        }
        return arrayForm;
    }

    /**
     * Clusterize is the main method that orchestrates all the steps of the
     * <i>Spatial Fuzzy C-Means</i>algorithm.
     * It updates the cluster center matrix <code>V</code>,
     * computes the distances between the samples and the cluster centers,
     * updates the cluster membership matrix <code>U</code>,
     * and updates it again according to an inputted spatial function
     * @param X data matrix (samples x features)
     * @param U clsuter membership matrix (samples x clusters)
     * @param V cluster centers matrix (clusters x features)
     * @param k number of clusters
     * @param tolerance threshold used to check the convergence between iterations
     * @param delegate the object to which we delegate to update each iteration status
     * @param m  the fuzzyness parameter
     * @param iterations the max number of allowed iterations
     * @param stopCriterion is an integer representing one of the possible stop criterions:
     * <b>0</b> - it computes the Frobenius norm on the difference between the actual
     * membership matrix <code>U</code> and the previous one from the last iteration
     * <b>1</b> - it computes the Frobenius norm on the difference between the actual
     * centroid matrix <code>V</code> and the previous one from the last iteration
     * <b>2</b> - it computes the max norm on the difference between the actual
     * membership matrix <code>U</code> and the previous one from the last iteration
     * <b>3</b> - it computes the max norm on the difference between the actual
     * centroid matrix <code>V</code> and the previous one from the last iteration
     * @param r is the radius of the window representing the neighbourhood for a sample
     * @param p is the weight for the membership value update
     * @param q is the weight for the spatial function
     * @param spatialFunction is an integer representing the choosen spatial function:
     * <b>0</b> - Computes the spatial function as the sum of the membership values
     * for the samples in the neighborhood according to the same cluster
     * <b>1</b> - Computes the spatial function as the number of clusters membership
     * in the neighborhood that are the greatest value for the neighbor according to
     * that cluster
     * @param width the width of the matrix form, as the samples can be viewed when
     * computing the neightborhood for each one
     * @param testing a boolean value that tells whether the algorithm is run for testing
     * @return an array of Objects containing the cluster membership matrix and
     * the cluster centroid matrix
     * @see #updateClusterCenterMatrix(float[][], float[][], float[][])
     * @see #checkConvergence(float[][], float[][], float[][], float[][], int)
     * @see #computeExponentialMembership(float[][], float[][], double)
     * @see #euclideanDistanceMatrix(float[][], float[][], double[][], double)
     * @see #updateClusterMembershipMatrix(float[][], float[][], float[][],
     * double, double[][])
     * @see #updateMembershipsWithSpatialInformation(float[][], int, double,
     * double, int, int, int[][], int[][], double[][])
     */
    private static Object[] clusterize(float [][] X, float [][] U, float [][] V,
                                       int k, double tolerance, ClusteringDelegate delegate,
                                       double m, long iterations, int stopCriterion,
                                       int r, double p, double q, int spatialFunction,
                                       int width, boolean testing){
        boolean converged = false;
        long count = 0;

        Object[] resultMatrixes = new Object[2];
        final int nFeatures = X[0].length;
        
        double [][] D = null;
        float [][] Um = null;
        Um = computeExponentialMembership(U, Um, m);
        System.out.println("\n\nUM\n\n");
        //printMatrix(Um);

        int [][] xy = precomputeMatrixForm(X.length, width);
        int [][] offset = precomputeArrayForm(X.length / width, width);
        double [][] uPhQ = null;
        double [] colCached = null;
        double [][] integralU = null;
        short [][] gCached = null;
        int [][] integralG = null;

        float [][] oldU = null;
        if (stopCriterion == 0 || stopCriterion == 2 || testing)
        {
            oldU = new float [U.length][k];
        }

        float [][] oldV = null;
        if (stopCriterion == 1 || stopCriterion == 3 || testing)
        {
            oldV = new float [V.length][nFeatures];
        }

        float oldJ = 0; //computeObjectiveFunction(X, U, V, D);
        float distance = Float.MAX_VALUE;

        while (count < iterations && distance > tolerance)
        {
//            System.out.println("it "+ count);

            V = updateClusterCenterMatrix(X, Um, V);
            //printMatrix(V);
//            if (count < 2){
//                System.out.println("\n\nV\n\n " + count);
//                printMatrix(V);
//            }

            D = euclideanDistanceMatrix(X, V, D, m);
//            if (count < 2){
//                System.out.println("\n\nD\n\n " + count);
//                printMatrix(D);
//            }
            
            U = updateClusterMembershipMatrix(X, U, V, m, D);
//            if (count < 2){
//                System.out.println("\n\nU\n\n " + count);
//                printMatrix(U);
//            }


            if ((r != 0) && (p != 1.0 || q != 0.0))
            {
//                U = updateMembershipsWithSpatialInformation(U, r, p, q, spatialFunction,
//                                                        width, xy, offset, uPhQ);
                
                if (spatialFunction == 0)
                {
                    U = opt1UpdateMembershipsWithSpatialFunctionH(U, r, p, q, width,
                                                                    colCached, uPhQ);
                }
                else if (spatialFunction == 1)
                {
                    U = opt1UpdateMembershipsWithSpatialFunctionG(U, r, p, q, width,
                                                                    colCached, uPhQ,
                                                                    gCached);
                }

//                if (spatialFunction == 0)
//                {
//                    U = opt2UpdateMembershipsWithSpatialFunctionH(U, r, p, q, width,
//                                                                  integralU, uPhQ);
//                }
//                else if (spatialFunction == 1)
//                {
//                    U = opt2UpdateMembershipsWithSpatialFunctionG(U, r, p, q, width,
//                                                                  integralG, uPhQ,
//                                                                  gCached);
//                }

            }

//                        if(count == 0){
//                for (k = 0; k < X.length; k++){
//                    System.out.println(U[k][2]);
//                }
//            }
            
            Um = computeExponentialMembership(U, Um, m);

            if (!testing)
            {
                distance = checkConvergence(U, oldU, V, oldV, stopCriterion);
            }
            

//            if (diffJ < tolerance)
//            {
//                converged = true;
//            }
//            else
//            {

//            }

            ++count;

            if (testing)
            {
                float diffJ, diffU, diffV;
                float newJ = computeObjectiveFunction(X, Um, V, D);
                diffJ = Math.abs(newJ - oldJ);
                oldJ = newJ;
                diffU = maxNorm(U, oldU);
                diffV = maxNorm(V, oldV);
                try
                {
                //delegate.updateStatus(null, null, count, diffJ);
                    delegate.updateStatus(null, null, count, diffJ, diffU, diffV);
                }
                catch (IOException ex)
                {
                    ex.printStackTrace();
                }
            }
            else
            {
                delegate.updateStatus(null, null, count, distance);
            }

            if (stopCriterion == 0 || stopCriterion == 2 || testing)
            {
                matrixCopy(U, oldU);
            }

            if (stopCriterion == 1 || stopCriterion == 3 || testing)
            {
                matrixCopy(V, oldV);
            }
            

        }

        //U = computeClusterMembership(X, U, V, k);
        resultMatrixes[0] = U;
        resultMatrixes[1] = V;
        return resultMatrixes;
    }
}
