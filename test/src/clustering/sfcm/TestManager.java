/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package clustering.sfcm;
import ij.ImagePlus;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FilenameFilter;
import java.io.RandomAccessFile;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author valerio
 */
public class TestManager extends junit.framework.TestCase{

    private String [] image;
    private String path;
    private int [] randomizationSeed = {48, 1000};
    private int [] numberOfClusters = {4, 5, 6};
    private final String[] initModes = {"random V (Forgy)", "K-Means++", "random U"};
    private double [] fuzzyness = {1.1f, 2, 5};
    private double [] p = {1.0, 2.0};
    private double [] q = {0.0, 1.0, 2.0};
    private int [] r = {1, 2, 5};
    private String [] spatialFunction = {"Likeliest Cluster", "Weightiest Cluster"};
    private final String[] colorSpaces = {"None", "XYZ", "L*a*b*", "HSB"};

    public void startValidationTest(RandomAccessFile filePointer) throws IOException
    {
            filePointer.writeBytes("K,M,P,Q,R,INIT,SF,CS,VPC,VPE,VXB\n");
    }

    public TestManager(final java.lang.String test) throws IOException{
         super(test);

         path = new File(".").getCanonicalPath();
         System.out.println(path);
    }


    public void endTest(RandomAccessFile filePointer) throws IOException{

            filePointer.close();

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



    public void test01Validation()
    {
        ij.io.Opener openfile = new ij.io.Opener();

        File file = null;
        RandomAccessFile filePointer = null;

        image =  getImagesInTestDir(new File("./test/data"));

        for(int i = 0; i < image.length; i++)
        {
            file = new File(path + "/test/log/" + getCurrentDateTime() + " " + image[i] + ".log");
            ImagePlus imp = null;
            imp = openfile.openImage(this.path + "/test/data/" + image[i]);
            try {
                filePointer = new RandomAccessFile(file, "rw");
                this.startValidationTest(filePointer);

                for(int j = 0; j < fuzzyness.length; j++)
                {
                for(int k = 0; k < numberOfClusters.length; k++)
                {
                    for(int t = 0; t < p.length; t++)
                    {
                        for(int u = 0; u < initModes.length; u++)
                        {
                            for(int g = 0; g < q.length; g++)
                            {
                                for(int rad = 0; rad < r.length; rad++)
                                {
                                    for(int y = 0; y < spatialFunction.length; y++)
                                    {
                                        for(int z = 0; z < colorSpaces.length; z++)
                                        {
                                            this.runValidationTest(imp, filePointer, fuzzyness[j], numberOfClusters[k],
                                                initModes[u],p[t],q[g],r[rad],colorSpaces[z],spatialFunction[y]);

                                        }
                                   }
                                }
                            }
                        }

                    }
                }
            }

            this.endTest(filePointer);

            } catch (FileNotFoundException ex) {
                Logger.getLogger(TestManager.class.getName()).log(Level.SEVERE, null, ex);
            }
            catch (java.io.IOException ex) {
                Logger.getLogger(TestManager.class.getName()).log(Level.SEVERE, null, ex);
            }

        }
    }
    private void runValidationTest(ImagePlus imp, RandomAccessFile filePointer, double m, int c,
                       String initmode,double p, double q, int r, String colorSpace, String sf){

            SFCMManager sfcmmInstance = new SFCMManager();
            sfcmmInstance.setTesting(false);
            sfcmmInstance.setValidation(true);
            sfcmmInstance.setNumberOfClusters(c);
            sfcmmInstance.setFuzzyness(m);
            sfcmmInstance.setRandomizationSeed(79);
            sfcmmInstance.setInitializationMode(initmode);
            sfcmmInstance.setMembershipWeight(p);
            sfcmmInstance.setSpatialFunction(sf);
            sfcmmInstance.setSpatialFunctionWeight(q);
            sfcmmInstance.setWindowRadius(r);
            sfcmmInstance.setColorSpace(colorSpace);
            sfcmmInstance.setTestFilePointer(filePointer);
            sfcmmInstance.setBinaryStackVisualization(false);
            sfcmmInstance.setClusterCenterColorsVisualization(false);
            sfcmmInstance.setGrayScaleVisualization(false);
            sfcmmInstance.setRandomRGBVisualization(false);
            sfcmmInstance.run(imp);

    }
}
