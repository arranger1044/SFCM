/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package clustering.sfcm;

import java.io.IOException;
import java.util.HashSet;
import java.util.Random;

/**
 * This class encapsulates the Spatial Fuzzy C-Means as proposed in: <i>Chuang,
 * Tzeng, Chen Fuzzy c-means clustering with spatial information for image
 * segmentation. Computerized Medical Imaging and Graphics 30 (2006) 9–15</i>.
 * It implements the Singleton pattern and provides as the only public class
 * method @link{run}. It also implements the delegation pattern, being free from
 * ImageJ API it can be  more easily reused
 */
public class SFCM {

    /**
     * A double value used for correcting numerical error when the spatial function
     * computed while updating the cluster membership matrix <code>U</code> equals 0
     * for each samples for a given cluster.
     * @see #opt1UpdateMembershipsWithSpatialFunctionG(float[][], int, double, double, int, double[], double[][], short[][]) 
     */
    private static double zeroCorrection = 1E-10;

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
     * @param seedEnabled boolean true if a random sequence is generated from a seed
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
     * @see SFCMManager
     * @see ClusteringDelegate
     * @see #defuzzyfyClusterMemberships(float[][])
     */
    public static Object[] run (float [][] X, int k, double tolerance, boolean seedEnabled,
                                int randomSeed, int initMode,
                                ClusteringDelegate delegate,
                                double m, long iterations, int stopCriterion,
                                int r, double p, double q, int spatialFunction,
                                int width,
                                boolean testing){

        /* Allocating U and V matrixes */
        float [][] V = new float[k][X[0].length];
        float [][] U = new float [X.length][k];

        /* Initializing U and V according to 'initMode' */
        initializeMatrixes(X, V, U, k, seedEnabled, randomSeed, initMode, m);

        /* Starting the clustering algorithm */
        Object[] clusteredMatrixes = clusterize(X, U, V, k, tolerance, delegate, 
                                                m, iterations, stopCriterion, r,
                                                p, q, spatialFunction, width,
                                                testing);

        /* Returning the matrixes in an array of objects*/
        Object[] matrixes = new Object[3];
        /* defuzzyfying the U matrix */
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
     * @param seedEnabled boolean true if a random sequence is generated from a seed
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
     * @param m the fuzzyness parameter
     * @return an Object array which containes the initialized <code>U</code>
     * and <code>V</code> matrixes
     * @see #randomInitialization(float[][], float[][], int, int)
     * @see #updateClusterMembershipMatrix(float[][], float[][], float[][], double, double[][])
     * @see #kMeansPlusPlusInitialization(float[][], float[][], int, int)
     * @see #initializeClusterMembershipRandom(float[][], int, int)
     */
    private static Object [] initializeMatrixes(float [][] X, float [][] V, 
                                                float [][] U, int k, boolean seedEnabled,
                                                int randomSeed, int initMode,
                                                double m){
        
        Object[] initedMatrixes = new Object[2];
        double [][] D = null;

        switch(initMode)
        {
            case 0:
                /* Random initialization on V first, then U is updated */
                V = randomInitialization(X, V, k, seedEnabled, randomSeed);
                D = euclideanDistanceMatrix(X, V, D, m);
                U = updateClusterMembershipMatrix(X, U, V, m, D);
                System.out.println("Random V");
                break;
            case 1:
                /* K-Means++ initialization on V, then U is updated */
                V = kMeansPlusPlusInitialization(X, V, k, seedEnabled, randomSeed);
                D = euclideanDistanceMatrix(X, V, D, m);
                U = updateClusterMembershipMatrix(X, U, V, m, D);
                System.out.println("K-Means++");
                break;
            case 2:
                /* Random initialization on U, no need to update V */
                U = initializeClusterMembershipRandom(U, k, seedEnabled, randomSeed);
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
     * @param seedEnabled boolean true if a random sequence is generated from a seed
     * @param randomSeed randomization seed used to generate a random sequence
     * @return the initialized cluster membership matrix
     * @see #initializeMatrixes(float[][], float[][], float[][], int, int, int, double)
     */
    private static float [][] initializeClusterMembershipRandom(float [][] U, int k,
                                                                boolean seedEnabled,
                                                                int randomSeed){

        int nClusters = U[0].length;
        final Random random = createRandom(seedEnabled, randomSeed);
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
        return U;
    }

    /**
     * Initializes the cluster center matrix <code>V</code> randomly by selecting
     * samples from X and checking not to pick an already selected sample
     * @param X data matrix
     * @param V cluster centroid matrix to be initialized
     * @param k number of clusters
     * @param seedEnabled boolean true if a random sequence is generated from a seed
     * @param randomSeed randomization seed used to generate a random sequence
     * @return the initialized cluster matrix
     * @see #initializeMatrixes(float[][], float[][], float[][], int, int, int, double)
     */
    private static float[][] randomInitialization(float [][] X, float [][] V, int k,
                                                  boolean seedEnabled, int randomSeed) {

        final Random random = createRandom(seedEnabled, randomSeed);
        final int nbPixels = X.length;
        final int nFeatures = X[0].length;

        final HashSet<Integer> centerLocations = new HashSet<Integer>();
        int clusterCreated = 0;

        /* We need to pick k clusters */
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
     * Initializes the cluster center matrix <code>V</code> according to the
     * <i>K-means++</i> initialization criterion. The first centroid is selected
     * randomly from the samples and the remaining ones are selected according
     * to a weighted probability distribution on the squared distance from the samples
     * and the nearest centroid, chosen from thealready selected centroids.
     * @link{http://en.wikipedia.org/wiki/K-means++}
     * @param X the data matrix
     * @param V the centroid matrix to be initialized
     * @param k the number of clusters
     * @param seedEnabled boolean true if a random sequence is generated from a seed
     * @param randomSeed randomization seed used to generate a random sequence
     * @return the initialized cluster center matrix
     * @see #initializeMatrixes(float[][], float[][], float[][], int, int, int, double)
     */
    private static float[][] kMeansPlusPlusInitialization(float [][] X, float [][] V, int k,
                                                          boolean seedEnabled, int randomSeed) {

        final Random random = createRandom(seedEnabled, randomSeed);
        final int nbPixels = X.length;
        final int nFeatures = X[0].length;

        /* Location of data used as cluster centers */
        final HashSet<Integer> centerLocation = new HashSet<Integer>();
        int clusterCreated = 0;
        /* Choose one center uniformly at random from among samples */
        {
            final int p = random.nextInt(nbPixels);
            centerLocation.add(p);
            System.arraycopy(X[p], 0, V[0], 0, nFeatures);
            clusterCreated++;
        }

        final double[] dp2 = new double[nbPixels];
        while (clusterCreated < k) {
            /* For each data point p compute D(p), the distance between p and
               the nearest center that has already been chosen. */
            double sum = 0;

            for (int offset = 0; offset < nbPixels; offset++) {

                /* Let's not pick already chosen centers */
                if (centerLocation.contains(offset)) {
                    continue;
                }
                /* Distance to closest cluster */
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
     * @param seedEnabled boolean true if a random sequence is generated from a seed
     * @param randomSeed an integer used to generate the Random instance
     * @return the created instance
     * @see #initializeMatrixes(float[][], float[][], float[][], int, int, int, double)
     */
    private static Random createRandom(boolean seedEnabled, int randomSeed) {
//        return config.isRandomizationSeedEnabled()
//                ? new Random(config.getRandomizationSeed())
//                : new Random();
        return (seedEnabled ?
            new Random(randomSeed)
            : new Random());
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

                D[i][j] = sum;

                /* If the distance is zero then we must 'mark' the last column element
                   with -1 */
                if (D[i][j] == 0.0)
                {
                    D[i][B.length] = -1;
                }
                else
                {
                    /* If the distance is not zero and the last column element has not
                     already been marked, let's store there the precomputed sum */
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
//                if ( Float.isNaN(Um[i][j]))
//                {
//                    throw new IllegalArgumentException("Nan found "+ U[i][j] + " m " + m);
//                }
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
                double distancePixelToCluster = D[i][j];
                objF += distancePixelToCluster * Um[i][j];
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
                /* The distance is stored in the Distance Matrix D */
                    double distance = D[i][j];
                    /* If the distance is not zero we can continue the computation */
                    if(distance != 0.0)
                    {
                        /* The precomputed denominator sum is stored in the last column
                           element of the matrix D */
                        double denominatorSum = D[i][nClusters];
                        /* If the sum has not been marked as -1 then we can compute
                           the U[i][j] element safely */
                        if (denominatorSum != -1)
                        {
                            double numerator = Math.pow(distance, exp);
                            U[i][j] = (float) (1f / (numerator * denominatorSum));
//                            if ( Float.isNaN(U[i][j]))
//                            {
//                                throw new IllegalArgumentException("Nan found " + " m " + m + " num " + numerator
//                                        + " den " + denominatorSum + " dis " + distance + " exp " + exp);
//                            }
                        }
                        else
                        {
                            U[i][j] = 0; /* We assign zero otherwise */
                        }
                    }
                    else
                    {
                        count++;
                        U[i][j] = 1.0f; /* If the distance is zero we assign 1 for now and
                                            count how many times for a row it becomes 1*/
                    }

//                    if ( Float.isNaN(U[i][j]))
//                    {
//                        throw new IllegalArgumentException("Nan found " + " m " + m);
//                    }
            }

            /* If for this row of U we assigned 1.0 to more than one sample,
               we must normalize*/
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

        double denominator = 0;

        int nClusters = V.length;
        int nFeatures = X[0].length;
        double [] numerator = null;

        for (int i = 0; i < nClusters; i++)
        {
            /* Creating an array which will store the intermediate values for each
               feature */
            numerator = new double[nFeatures];
            denominator = 0;
            for(int k = 0; k < X.length; k++)
            {
                float x[] = X[k];
                double um = Um[k][i];
                for (int j = 0; j < nFeatures; j++)
                {
                    numerator[j] += um * x[j];
                }

                denominator += um;
            }

            /* Updating each element of V cycling through the features*/
            for (int j = 0; j < nFeatures; j++)
            {
                 V[i][j] = (float) (numerator[j] / denominator);
//                 if ( Float.isNaN(V[i][j]))
//                 {
//                    for(int k = 0; k < X.length; k++)
//                    {
//                        //System.out.println(Um[k][i]);
//                    }
//                    throw new IllegalArgumentException("nume "+ numerator + " denom "
//                            + denominator + "i " + i);
//                 }
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

                double uPTimeshQ = Math.pow(U[i][j], p) *
                                   Math.pow(h, q);
                uPhQ[i][j] = uPTimeshQ;
                uPhQ[i][nClusters] += uPTimeshQ;

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
     * Computes the update rule of the membership matrix <code>U</code> according
     * to the <i>Spatial Fuzzy C-Mean</i> algorithm.
     * It implements the spatial function as the sum of the membership values for
     * the data belonging to the same neighborhood.
     * An optimization is implemented when computing the neighbourhood window.
     * An array, with length equal to the
     * number of the matrix columns, is used as a queue to compute che sum of the
     * window column values and to efficiently and incrementally move the window
     * along the columns. For each row the array is initially recomputed from scratch
     * @param U the cluster membership matrix to update
     * @param rad the window radius
     * @param p the membership weight
     * @param q the spatial function weight
     * @param cols the number of columns in matrix form
     * @param supC the array acting as a queue, if null it will be allocated
     * @param uPhQ a support matrix used to store intermediate computed values,
     * if null it will be allocated
     * @return the updated cluster membership matrix
     */
    private static float [][] opt1UpdateMembershipsWithSpatialFunctionH(float [][] U, int rad,
                                                                          double p, double q,
                                                                          int cols, double [] supC,
                                                                          double [][] uPhQ){
        int nClusters = U[0].length;

        /* Allocating the null data structures */
        if(uPhQ == null)
        {
            uPhQ = new double[U.length][nClusters + 1];
        }
        if (supC == null)
        {
            supC = new double[cols];
        }

        int rows = U.length / cols;

        /* Zeroing the uPhQ support matrix */
        for(int i = 0; i < U.length; i++)
        {
            uPhQ[i][nClusters] = 0;
        }

        double uPTimeshQ  = 0;

        /* For each column of U*/
        for(int k = 0; k < nClusters; k++)
        {
            /* Moving along rows but considering the data sample in matrix form
               rows * cols */
            for (int i = 0; i < rows; ++i)
            {
                double h = 0;
                /* For each new row of the matrix form we recompute supC starting
                 from zero */
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
                                mem += U[y * cols + rx][k];
                            }
                        }
                        /* Adding a new column */
                        supC[qLast] = mem;
                        qLast++;
                        h += mem;
                    }
                }
                int j = 0, col, y;

                /* We can already compute the summed window value for the first
                 element of the row in matrix form */
                uPTimeshQ = Math.pow(U[i * cols][k], p) *
                            Math.pow(h, q);
                uPhQ[i * cols][k] = uPTimeshQ;
                uPhQ[i * cols][nClusters] += uPTimeshQ;

                /* Now we cycle through the columns*/
                qFirst = 0;
                qLast = rad + 1;
                for (j = 1; j < cols; j++)
                {
                    col = j + rad;
                    /* If we can, we remove the first column in supC */
                    if (j - rad > 0)
                    {
                        h -= supC[qFirst];
                        qFirst++;
                    }
                    /* and if we have not reached the last usable column we add
                     another value to supC at the end of the queue*/
                    if (j + rad < cols)
                    {
                        double mem = 0;
                        for (int ry =- rad; ry <= rad; ry++)
                        {
                            y = i + ry;

                            if (y >= 0 && y < rows)
                            {
                                mem += U[y * cols + col][k];
                            }
                        }
                        supC[qLast] = mem;
                        qLast++;
                        h += mem;
                    }

                    /* In the end we compute the exponentiated value and store it */
                    int offset = i * cols + j;
                    uPTimeshQ = Math.pow(U[offset][k], p) *
                                Math.pow(h, q);
                    uPhQ[offset][k] = uPTimeshQ;
                    uPhQ[offset][nClusters] += uPTimeshQ;
                }
            }
        }

        /* Updating U after all the values have been computed */
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
     * Computes the update rule of the membership matrix <code>U</code> according
     * to the <i>Spatial Fuzzy C-Mean</i> algorithm.
     * It implements the spatial function as the number of the cluster membership values that
     * are the greates for every pixel in the neighborhood. So it will compute
     * a defuzzyfied version of U and apply the computation on it.
     * An optimization is implemented when computing the neighbourhood window.
     * An array, with length equal to the
     * number of the matrix columns, is used as a queue to compute che sum of the
     * window column values and to efficiently and incrementally move the window
     * along the columns. For each row the array is initially recomputed from scratch.
     * @param U the cluster membership matrix to update
     * @param rad the window radius
     * @param p the membership weight
     * @param q the spatial function weight
     * @param cols the number of columns in matrix form
     * @param supC the array acting as a queue, if null it will be allocated
     * @param uPhQ a support matrix used to store intermediate computed values,
     * if null it will be allocated
     * @param G a support matrix used as
     * @return the updated cluster membership value
     */
    private static float [][] opt1UpdateMembershipsWithSpatialFunctionG(float [][] U, int rad,
                                                                        double p, double q,
                                                                        int cols, double [] supC,
                                                                        double [][] uPhQ,
                                                                        short [][] G){
        int nClusters = U[0].length;

        /* Allocating the needed support data structures if null */
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

        /* Allocating a support array for checking if a column of G goes to 0*/
        int [] g0count = new int [nClusters];

        /* Computing the G matrix just once as a defuzzyfied version of U */
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

        /* Zeroing the uPhQ support matrix */
        for(int i = 0; i < U.length; i++)
        {
            uPhQ[i][nClusters] = 0;
        }

        double uPTimeshQ  = 0;
        /* For each column of U*/
        for(int k = 0; k < nClusters; k++)
        {
            /* Moving along rows but considering the data sample in matrix form
               rows * cols */
            for (int i = 0; i < rows; ++i)
            {
                int h = 0;
                /* For each new row of the matrix form we recompute supC starting
                   from zero */
                int qFirst = 0;
                int qLast = 0;
                for (int rx = 0; rx <= rad; rx++)
                {
                    if (rx < cols)
                    {
                        int mem = 0;
                        for (int ry =- rad; ry <= rad; ry++)
                        {
                            int y = i + ry;

                            if (y >= 0 && y < rows)
                            {
                                mem += G[y * cols + rx][k];
                            }
                        }
                        /* Adding a new column */
                        supC[qLast] = mem;
                        qLast++;
                        h += mem;
                    }
                }
                int j = 0, col, y;
                /* If h goes to zero we must keep track of it */
                if (h == 0)
                {
                    g0count[k]++;
                }
                /* We can already compute the summed window value for the first
                 element of the row in matrix form */
                uPTimeshQ = Math.pow(U[i * cols][k], p) *
                            Math.pow(h, q);
                uPhQ[i * cols][k] = uPTimeshQ;
                uPhQ[i * cols][nClusters] += uPTimeshQ;

                qFirst = 0;
                qLast = rad + 1;
                
                /* Now we cycle through the columns*/
                for (j = 1; j < cols; j++)
                {
                    col = j + rad;
                    /* If we can, we remove the first column in supC */
                    if (j - rad > 0)
                    {
                        h -= supC[qFirst];
                        qFirst++;
                    }
                    /* and if we have not reached the last usable column we add
                       another value to supC at the end of the queue*/
                    if (j + rad < cols)
                    {
                        int mem = 0;
                        for (int ry =- rad; ry <= rad; ry++)
                        {
                            y = i + ry;

                            if (y >= 0 && y < rows)
                            {
                                mem += G[y * cols + col][k];
                            }
                        }
                        supC[qLast] = mem;
                        qLast++;
                        h += mem;

                    }
                    /* Again we must check if h equals zero */
                    if (h == 0)
                    {
                        g0count[k]++;
                    }
                    int offset = i * cols + j;
                    uPTimeshQ = Math.pow(U[offset][k], p) *
                                Math.pow(h, q);
                    uPhQ[offset][k] = uPTimeshQ;
                    uPhQ[offset][nClusters] += uPTimeshQ;
                }
            }
        }

        /* Checking the support array to see if a column is completly zero valued */
        for(int j = 0; j < nClusters; j++)
        {
            if(g0count[j] == U.length)
            {
                /* The column is zero valued, we must correct uPhQ */
                for (int i = 0; i < U.length; i++)
                {
                    uPhQ[i][j] = zeroCorrection;
                }
            }
        }

        /* Updating U after all the values have been computed */
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

        /* Allocating a support array for checking if a column of G goes to 0*/
        int [] g0count = new int [nClusters];

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
                    if (tot == 0)
                    {
                        g0count[k]++;
                    }
                    int offset = y * cols + x;
                    double uPTimeshQ = Math.pow(U[offset][k], p) *
                                       Math.pow(tot, q);
                    uPhQ[offset][k] = uPTimeshQ;
                    uPhQ[offset][nClusters] += uPTimeshQ;
                }
            }
        }

        /* Checking the support array to see if a column is completly zero valued */
        for(int j = 0; j < nClusters; j++)
        {
            if(g0count[j] == U.length)
            {
                /* The column is zero valued, we must correct uPhQ */
                for (int i = 0; i < U.length; i++)
                {
                    uPhQ[i][j] = zeroCorrection;
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
     * Deprecated.
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
     * Deprecated.
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

//        int [][] xy = precomputeMatrixForm(X.length, width);
//        int [][] offset = precomputeArrayForm(X.length / width, width);

        /* Support data structures */
        double [][] uPhQ = null;
        double [] colCached = null;
        double [][] integralU = null;
        short [][] gCached = null;
        int [][] integralG = null;

        /* Allocating oldU and oldV only if testing or an appropriate stop
         criterion has been chosen */
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

        /* Cycling till the number of iterations gets past the max number allowed
         or the convergence distance value is greater than a threshold tolerance */
        while (count < iterations && distance > tolerance)
        {

            /* Updating the cluster centroid matrix */
            V = updateClusterCenterMatrix(X, Um, V);

            /* Computing the distance between the new centroids and the data */
            D = euclideanDistanceMatrix(X, V, D, m);

            /* Updating the cluster membership matrix */
            U = updateClusterMembershipMatrix(X, U, V, m, D);

            /* If p == 1 and q == 0 OR r == 0 there is no need to call the
             cluster membership matrix update spatial update since it will be
             a classic Fuzzy C-Means version */
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

            /* Exponentiating the cluster membership matrix */
            Um = computeExponentialMembership(U, Um, m);

            /* Checking the convergence , only if we are not testing */
            if (!testing)
            {
                distance = checkConvergence(U, oldU, V, oldV, stopCriterion);
            }

            ++count;

            /* If testing we must compute the convergence criterion on U, V and also J */
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

            /* Updating oldU and oldV selectively according to the stop criterion */
            if (stopCriterion == 0 || stopCriterion == 2 || testing)
            {
                matrixCopy(U, oldU);
            }

            if (stopCriterion == 1 || stopCriterion == 3 || testing)
            {
                matrixCopy(V, oldV);
            }
            

        }

        resultMatrixes[0] = U;
        resultMatrixes[1] = V;
        return resultMatrixes;
    }
}
