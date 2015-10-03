package org.jgrasstools.hortonmachine.models.hm;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;

import org.geotools.coverage.grid.GridCoverage2D;
import org.jgrasstools.gears.io.rasterreader.RasterReader;
import org.jgrasstools.gears.libs.monitor.PrintStreamProgressMonitor;
import org.jgrasstools.hortonmachine.modules.statistics.cb.OmsCb;
import org.jgrasstools.hortonmachine.utils.HMTestCase;

/* Test for the Cb module.
 * 
 * @author Andrea Antonello (www.hydrologis.com)
 */
public class Cb_for_meMATTEO extends HMTestCase {

	public void testCb() throws Exception {

		String pathh = "/Users/giuseppeformetta/Desktop/VAlmazia/PIO/cell"; //to the location

		String files;
		File folder = new File(pathh);
		File[] listOfFiles = folder.listFiles();
		String[] vet=new String[listOfFiles.length-1];
		for (int i = 1; i < listOfFiles.length; i++) {
			System.out.println(listOfFiles[i].getName());
			//vet[i-1]=java.util.Arrays.toString(listOfFiles[i].getName().split(".asc")).toString();
			//System.out.println(vet[i-1]);

		}

		 String path = new File(
		 "/Users/giuseppeformetta/Desktop/VAlmazia/VENETO/cell/basinnumbered") //path to basin number 
		 .getAbsolutePath();
		 RasterReader reader = new RasterReader();
			reader.file = path;
			reader.fileNovalue = -9999.0;
			reader.geodataNovalue = Double.NaN;
			reader.process();
			GridCoverage2D rasterbacini = reader.outRaster;

		for (int ii = 1; ii < listOfFiles.length; ii++) {
			
			

			String path2 = listOfFiles[ii].toString();
			RasterReader reader2 = new RasterReader();
			reader2.file = path2;
			reader2.fileNovalue = -9999.0;
			reader2.geodataNovalue = Double.NaN;
			reader2.process();
			
			GridCoverage2D h2cd = reader2.outRaster;

			PrintStreamProgressMonitor pm = new PrintStreamProgressMonitor(
					System.out, System.out);

			OmsCb cb = new OmsCb();
			cb.pBins = 0;
			cb.pFirst = 1;
			cb.pLast = 4;
			cb.inRaster1 = rasterbacini;
			cb.inRaster2 = h2cd;
			cb.pm = pm;

			cb.process();
			String a="/Users/giuseppeformetta/Desktop/venetoInputdata/"; //to write
			double[][] moments = cb.outCb;
			System.out.println(a.concat(listOfFiles[ii].getName()));
			String ppp=a.concat(listOfFiles[ii].getName());

			FileWriter Rstatfile = new FileWriter(ppp);
			PrintWriter errestat = new PrintWriter(Rstatfile);
			for (int i = 0; i < (moments.length); i++) {
				for (int j = 0; j < (moments[0].length); j++) {

					errestat.print(moments[i][j] + " ");
					System.out.print(moments[i][j] + " ");

				}
				errestat.println();
				System.out.println();
			}

			Rstatfile.close();
		}

		System.out.print("ciao");

	}

}
