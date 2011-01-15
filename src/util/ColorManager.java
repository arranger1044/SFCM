/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package util;

import java.util.Random;

/**
 * Provides some random color generators, used in SFCMManager to display
 * the clustered image. Implements the Singleton pattern
 * @see SFCMManager
 */
public class ColorManager {

    /**
     * ColorManager shall not be instantiated
     */
    private ColorManager(){}

    /**
     * Generates a random color in the RGB color space
     * @return an integer array representing the rgb components of the generated color
     */
    static public int [] randomRGBColor(){
        
        Random rnd = new Random();
        int red = rnd.nextInt(256);
        int green = rnd.nextInt(256);
        int blue = rnd.nextInt(256);

        int [] color = new int [3];
        color[0] = red;
        color[1] = green;
        color[2] = blue;

        return color;
    }

    /**
     * Generates a random color in the RGB color space using a mixing color to
     * blend the hues
     * @param mixingColor an integer array representing the rgb components of
     * the mixing color
     * @return an integer array representing the rgb components of the generated color
     */
    static public int [] randomMixedRGBColor(int [] mixingColor){

        Random rnd = new Random();
        int red = rnd.nextInt(256);
        int green = rnd.nextInt(256);
        int blue = rnd.nextInt(256);

        if (mixingColor != null)
        {
            red = (red + mixingColor[0]) / 2;
            green = (green + mixingColor[1]) / 2;
            blue = (blue + mixingColor[2]) / 2;
        }

        int [] color = new int [4];
        color[0] = red;
        color[1] = green;
        color[2] = blue;
        color[3] = 255;

        return color;
    }
}
