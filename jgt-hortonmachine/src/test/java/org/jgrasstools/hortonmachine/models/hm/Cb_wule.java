package org.jgrasstools.hortonmachine.models.hm;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;

import org.geotools.coverage.grid.GridCoverage2D;
import org.jgrasstools.gears.io.rasterreader.OmsRasterReader;
import org.jgrasstools.gears.libs.monitor.PrintStreamProgressMonitor;
import org.jgrasstools.hortonmachine.modules.statistics.cb.OmsCb;
import org.jgrasstools.hortonmachine.utils.HMTestCase;

/* Test for the Cb module.
 * 
 * @author Andrea Antonello (www.hydrologis.com)
 */
public class Cb_wule extends HMTestCase {

	public void testCb() throws Exception {

		String pathh = "/Users/administrator/Documents/Posina_Loc3003/rainfallmap/rainfalmap/cell/";

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
		 "/Users/administrator/Documents/Posina_Loc3003/rainfallmap/rainfalmap/cell/numberedsubbasin")
		 .getAbsolutePath();
		 OmsRasterReader reader = new OmsRasterReader();
			reader.file = path;
			reader.fileNovalue = -9999.0;
			reader.geodataNovalue = Double.NaN;
			reader.process();
			GridCoverage2D rasterbacini = reader.outRaster;

		for (int ii = 1; ii < listOfFiles.length; ii++) {
			
			

			String path2 = listOfFiles[ii].toString();
			OmsRasterReader reader2 = new OmsRasterReader();
			reader2.file = path2;
			reader2.fileNovalue = -9999.0;
			reader2.geodataNovalue = Double.NaN;
			reader2.process();
			
			GridCoverage2D rainfall = reader2.outRaster;

			PrintStreamProgressMonitor pm = new PrintStreamProgressMonitor(
					System.out, System.out);

			OmsCb cb = new OmsCb();
			cb.pBins = 100;
			cb.pFirst = 1;
			cb.pLast = 4;
			cb.inRaster1 = rasterbacini; //Users/administrator/Documents/Posina_Loc3003/rainfallmap/rainfallmoment
			cb.inRaster2 = rainfall;
			cb.pm = pm;
			cb.process();
			
			int c = 0;
			FileWriter Rstatfile = new FileWriter(
					"/Users/administrator/Documents/Posina_Loc3003/rainfallmap/rainfallmoment"
							+c+".csv");
			
			//int iii = 0;
			//String a="/Users/administrator/Documents/Posina_Loc3003/rainfallmap/rainfallmoment"; //to write
			double[][] moments = cb.outCb;
			//System.out.println(a.concat(listOfFiles[ii].getName()));
			//String ppp=a.concat(listOfFiles[ii].getName());

		//	FileWriter Rstatfile = new FileWriter(ppp);
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
