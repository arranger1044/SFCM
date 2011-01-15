/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package clustering.sfcm;

import java.io.IOException;

/**
 * Interface which must be implemented by a delegate class that will be asked to
 * use the ImageJ api in order to update the status of the clustering algorithm
 * (SFCM). ie the delegating actor.
 * @see SFCM
 * @see SFCMManager
 */
public interface ClusteringDelegate {

    /**
     * Provides a way to communicate the current status of the clustering algorithm
     * via the current values in the cluster center matrix, cluster membership
     * matrix, current iteration and error from the previous iteration
     * @param V the cluster center matrix
     * @param U the cluster membership matrix
     * @param nIteration the current iteration number
     * @param error the current iteration error
     */
    public void updateStatus(float [][] V, int [][] U, long nIteration, float error);

    /**
     * Provides a way to communicate the current status of the clustering algorithm
     * via a message
     * @param message a string containing a message describing the current status
     */
    public void updateStatus(String message);

    /**
     * Provides a way to communicate the current status of the clustering algorithm
     * via a more detailed view on its current iteration internal variable.
     * It is used when testing. It tries to write to file
     * @param V the cluster center matrix
     * @param U the cluster memebership matrix
     * @param nIteration the current iteration
     * @param errorJ the current iteration error of the objective function
     * @param errorU the current iteration error of the memebership matrix
     * @param errorV the current iterration error of the cluster center matrix
     * @throws IOException exception cause by trying to write to file
     * @see TestManager
     */
    public void updateStatus(float [][] V, int [][] U, long nIteration, 
                             float errorJ, float errorU, float errorV) throws IOException;

}
