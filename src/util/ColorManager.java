/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package util;

import java.util.Random;

/**
 *
 * @author valerio
 */
public class ColorManager {

    private ColorManager(){}

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
