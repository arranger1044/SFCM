/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package clustering.fcm;
import ij.ImagePlus;
import java.io.File;
import java.io.FilenameFilter;
import java.io.RandomAccessFile;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;

/**
 *
 * @author valerio
 */
public class TestManager extends junit.framework.TestCase{

    private String [] image;
    private String path;
    private int [] randomizationSeed = {48, 1000};
    private int [] numberOfClusters = {4, 12, 20};
    private final String[] initModes = {"random V (Forgy)", "K-Means++", "random U"};
    private float [] fuzzyness = {1.1f, 2, 5};

    public TestManager(final java.lang.String test) throws IOException{
         super(test);

         path = new File(".").getCanonicalPath();
         System.out.println(path);
    }

    public void startTest(RandomAccessFile filePointer, float m, int c,
                                           int randomSeed, String image, String initmode) throws IOException{


            filePointer.writeBytes("file name : " + image + "\n");
            filePointer.writeBytes("fuzziness : " + m + "\n");
            filePointer.writeBytes("number of clusters " + c + "\n");
            filePointer.writeBytes("initmode " + initmode + "\n");
            filePointer.writeBytes("random seed " + randomSeed + "\n\n");

            filePointer.writeBytes("ITERATION,DIFFJ,DIFFU,DIFFV\n");
    }

    public void endTest(RandomAccessFile filePointer) throws IOException{

            filePointer.close();

    }

    public void runTest(ImagePlus imp, RandomAccessFile filePointer, float m, int c,
                        int randomSeed, String initmode){

            FCMManager fcmmInstance = new FCMManager();
            fcmmInstance.setTesting(true);
            fcmmInstance.setNumberOfClusters(c);
            fcmmInstance.setFuzzyness(m);
            fcmmInstance.setRandomizationSeed(randomSeed);
            fcmmInstance.setInitializationMode(initmode);

            fcmmInstance.setTestFilePointer(filePointer);
            fcmmInstance.setBinaryStackVisualization(false);
            fcmmInstance.setClusterCenterColorsVisualization(false);
            fcmmInstance.setGrayScaleVisualization(false);
            fcmmInstance.setRandomRGBVisualization(false);
            
            fcmmInstance.run(imp);
    }

    private String getCurrentDateTime() {

        DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd HH.mm.ss");
        Date date = new Date();
        return dateFormat.format(date);
    }

    private String[] getImagesInTestDir(File dir){

        FilenameFilter filter = new FilenameFilter() {
            @Override
            public boolean accept(File dir, String name) {

                return !name.startsWith(".");
            }
        };

        String[] files = dir.list(filter);

        System.out.println(files.length);
        return files;
    }

    public void test(){


        ij.io.Opener openfile = new ij.io.Opener();

        File file = null;
        RandomAccessFile filePointer = null;

        image =  getImagesInTestDir(new File("./test/data"));

        for(int i = 0; i < image.length; i++)
        {
            for(int j = 0; j < fuzzyness.length; j++)
            {
                for(int k = 0; k < numberOfClusters.length; k++)
                {
                    for(int p = 0; p < randomizationSeed.length; p++)
                    {
                        for(int u = 0; u < initModes.length; u++)
                        {
                            try
                            {
                                System.out.println(image[i] + " fuzzyness:" +  fuzzyness[j] +
                                        " k:" + numberOfClusters[k] + " rnd:" + randomizationSeed[p] +
                                        " init:" + initModes[u]);
                                ImagePlus imp = null;
                                imp = openfile.openImage(this.path + "/test/data/" + image[i]);
                                file = new File(path + "/test/log/" + getCurrentDateTime() + " " + image[i] +
                                        "." + fuzzyness[j] + "." + numberOfClusters[k] + "." +
                                        randomizationSeed[p] + "." + initModes[u] +".log");
                                filePointer = new RandomAccessFile(file, "rw");
                                this.startTest(filePointer, fuzzyness[j], numberOfClusters[k],
                                        randomizationSeed[p], image[i], initModes[u]);
                                this.runTest(imp, filePointer, fuzzyness[j], numberOfClusters[k],
                                        randomizationSeed[p], initModes[u]);
                                this.endTest(filePointer);
                            }
                            catch (IOException ex)
                            {
                                ex.printStackTrace();
                            }
                        }

                    }
                }
            }
        }
    }
}

