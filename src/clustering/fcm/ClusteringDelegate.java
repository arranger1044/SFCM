/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package clustering.fcm;

import java.io.IOException;

/**
 *
 * @author valerio
 */
public interface ClusteringDelegate {

    public void updateStatus(float [][] V, int [][] U, long nIteration, float error);

    public void updateStatus(String message);

    public void updateStatus(float [][] V, int [][] U, long nIteration, 
                             float errorJ, float errorU, float errorV) throws IOException;

}
