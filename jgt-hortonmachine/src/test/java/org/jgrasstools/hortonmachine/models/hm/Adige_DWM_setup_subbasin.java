package org.jgrasstools.hortonmachine.models.hm;

import java.io.File;

import org.geotools.coverage.grid.GridCoverage2D;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.jgrasstools.gears.io.rasterreader.RasterReader;
import org.jgrasstools.gears.io.rasterwriter.OmsRasterWriter;
import org.jgrasstools.gears.io.shapefile.OmsShapefileFeatureWriter;
import org.jgrasstools.gears.libs.modules.Variables;
import org.jgrasstools.gears.modules.r.cutout.OmsCutOut;
import org.jgrasstools.gears.modules.v.vectorize.OmsVectorizer;
import org.jgrasstools.hortonmachine.modules.basin.basinshape.OmsBasinShape;
import org.jgrasstools.hortonmachine.modules.demmanipulation.markoutlets.OmsMarkoutlets;
import org.jgrasstools.hortonmachine.modules.demmanipulation.wateroutlet.OmsExtractBasin;
import org.jgrasstools.hortonmachine.modules.demmanipulation.wateroutlet.OmsWateroutlet;
import org.jgrasstools.hortonmachine.modules.geomorphology.aspect.OmsAspect;
import org.jgrasstools.hortonmachine.modules.geomorphology.draindir.OmsDrainDir;
import org.jgrasstools.hortonmachine.modules.geomorphology.flow.OmsFlowDirections;
import org.jgrasstools.hortonmachine.modules.geomorphology.slope.OmsSlope;
import org.jgrasstools.hortonmachine.modules.network.extractnetwork.OmsExtractNetwork;
import org.jgrasstools.hortonmachine.modules.network.netnumbering.OmsNetNumbering;
import org.jgrasstools.hortonmachine.modules.network.networkattributes.OmsNetworkAttributesBuilder;
import org.jgrasstools.hortonmachine.utils.HMTestCase;

/* Test for the Cb module.
 * 
 * @author Andrea Antonello (www.hydrologis.com)
 */
public class Adige_DWM_setup_subbasin extends HMTestCase {

	public void testCb() throws Exception {

		

	
			
	/*		OmsPitfiller pit = new OmsPitfiller();
			pit.inElev = EUdem_projected2_32632;
			pit.process();
			GridCoverage2D pitcoverage = pit.outPit;
			OmsRasterWriter.writeRaster("/Users/administrator/Documents/Adige_project/FILL_DEM.asc", pitcoverage);
			*/
			String path3 = new File(
					 "/Users/administrator/Documents/Adige_project/FILL_DEM.asc") 
					 .getAbsolutePath();
					 RasterReader reader3 = new RasterReader();
						reader3.file = path3;
						reader3.fileNovalue = -9999.0;
						reader3.geodataNovalue = Double.NaN;
						reader3.process();
						GridCoverage2D pitcoverage = reader3.outRaster;
						
			
			
			OmsFlowDirections  flowdir = new OmsFlowDirections();
			flowdir.inPit = pitcoverage;
			flowdir.process();
			GridCoverage2D flowcoverage = flowdir.outFlow;
			
			OmsSlope slope = new OmsSlope();
			slope.inPit = pitcoverage;
			slope.inFlow = flowcoverage;
			slope.process();
			GridCoverage2D slopecoverage = slope.outSlope;
			
					
			OmsDrainDir draindir = new OmsDrainDir();
			draindir.inPit = pitcoverage;
			draindir.inFlow = flowcoverage; 
			draindir.pLambda = 1;
			draindir.process();
			GridCoverage2D draindircoverage = draindir.outFlow;
			GridCoverage2D tcacoverage = draindir.outTca;
		//	OmsRasterWriter.writeRaster("/Users/administrator/Documents/Adige_project/tca_test_1.asc", tcacoverage);
			
			
			
			OmsMarkoutlets markout = new OmsMarkoutlets();
			markout.inFlow = draindircoverage;
			markout.process();
			GridCoverage2D markoutcoverage = markout.outFlow;
			
			
			OmsExtractNetwork extractnet = new OmsExtractNetwork();
			extractnet.inFlow = markoutcoverage;
			extractnet.inTca = tcacoverage;
			//extractnet.inSlope = slopecoverage;
			extractnet.pThres = 25;
			extractnet.pMode = Variables.TCA;
			//extractnet.pExp = 0.5;
			extractnet.process();
			GridCoverage2D networkcoverage = extractnet.outNet;
		//	OmsRasterWriter.writeRaster("/Users/administrator/Documents/Adige_project/net_A50.asc", networkcoverage);
			
			
			OmsWateroutlet outlet = new OmsWateroutlet();
			outlet.inFlow = markoutcoverage;
			outlet.pEast =  661741.6993;
			outlet.pNorth = 5109593.766;
			outlet.process();
			
			GridCoverage2D mybasin = outlet.outBasin;
			//GridCoverage2D basinareacoverage = outlet.outArea;
		
		
			
			OmsCutOut pitcutout = new OmsCutOut();
			pitcutout.inMask = mybasin;
			pitcutout.inRaster = pitcoverage;
			pitcutout.process();
			GridCoverage2D pitcoverage2 = pitcutout.outRaster;
			
			OmsCutOut draindircutout = new OmsCutOut();
			draindircutout.inMask = mybasin;
			draindircutout.inRaster = draindircoverage;
			draindircutout.process();
			GridCoverage2D draindircoverage2 = draindircutout.outRaster;
			
			
			
			OmsCutOut Markoutcutout = new OmsCutOut();
			Markoutcutout.inMask = mybasin;
			Markoutcutout.inRaster = markoutcoverage;
			Markoutcutout.process();
			GridCoverage2D markoutcoverage2 = Markoutcutout.outRaster;
			
			OmsCutOut networkcutout = new OmsCutOut();
			networkcutout.inMask = mybasin;
			networkcutout.inRaster = networkcoverage;
			networkcutout.process();
			GridCoverage2D networkcoverage2 = networkcutout.outRaster;
			OmsRasterWriter.writeRaster("/Users/administrator/Documents/Adige_project/net_Avisio_2.asc", networkcoverage2);
			
			
			OmsCutOut tcacutout = new OmsCutOut();
			tcacutout.inMask = mybasin;
			tcacutout.inRaster = tcacoverage;
			tcacutout.process();
			GridCoverage2D tcacoverage2 = tcacutout.outRaster;
		//	OmsRasterWriter.writeRaster("/Users/administrator/Documents/Adige_project/tca_test.asc", tcacutout.outRaster);
			
		
			
					
			
			OmsNetworkAttributesBuilder netbuilder = new OmsNetworkAttributesBuilder();
			netbuilder.inDem = pitcoverage2;
			netbuilder.inFlow = markoutcoverage2;
			netbuilder.inNet = networkcoverage2;
			netbuilder.inTca = tcacoverage2;
			netbuilder.process();
			
			GridCoverage2D hacknetcoverage = netbuilder.outHack;;
			
			SimpleFeatureCollection netbiuldercoverage = netbuilder.outNet;
			OmsShapefileFeatureWriter.writeShapefile("/Users/administrator/Documents/Adige_project/net_atribute_Avisio.shp", netbiuldercoverage, pm);
			
			
			OmsExtractBasin extractbasin = new OmsExtractBasin();
			extractbasin.inFlow = markoutcoverage;
			extractbasin.inNetwork = netbiuldercoverage;
			extractbasin.pEast = 661741.6993;
			extractbasin.pNorth = 5109593.766;
			extractbasin.pSnapbuffer = 100;
			extractbasin.process();
			
			OmsShapefileFeatureWriter.writeShapefile("/Users/administrator/Documents/Adige_project/maskshape_Avisio.shp", extractbasin.outVectorBasin, pm);
			
					
			
	       
			OmsNetNumbering netnumber = new OmsNetNumbering();
			netnumber.inFlow = markoutcoverage2;
			netnumber.inNet = networkcoverage2;
			//netnumber.inTca = tcacoverage;
			//netnumber.pThres = 10000;
			netnumber.process();
			
			GridCoverage2D netnumercoverage = netnumber.outNetnum;
			OmsRasterWriter.writeRaster("/Users/administrator/Documents/Adige_project/netnumber_Avisio.asc", netnumercoverage);
			GridCoverage2D subbasincoverage = netnumber.outBasins;
			OmsRasterWriter.writeRaster("/Users/administrator/Documents/Adige_project/subbasin_Avisio.asc", subbasincoverage);
			
			
		
			
			OmsVectorizer vectorizer = new OmsVectorizer();
			vectorizer.inRaster = netnumercoverage;
			vectorizer.process();
		    OmsShapefileFeatureWriter.writeShapefile("/Users/administrator/Documents/Adige_project/netnumshapefeature_Avisio.shp", vectorizer.outVector, pm);
			
			
			OmsBasinShape basinshape = new OmsBasinShape();
			basinshape.inElev = pitcoverage2;
			basinshape.inBasins = subbasincoverage;
			basinshape.process();
			
			SimpleFeatureCollection basinshapefeature = basinshape.outBasins;
			OmsShapefileFeatureWriter.writeShapefile("/Users/administrator/Documents/Adige_project/basinshape_Avisio.shp", basinshapefeature, pm);
			
	
			
			OmsCutOut aspectcutout = new OmsCutOut();
			aspectcutout.inMask = mybasin;
			aspectcutout.inRaster = pitcoverage2;
			aspectcutout.process();
			GridCoverage2D dem4aspectcoverage = aspectcutout.outRaster;
			
			OmsAspect aspect = new OmsAspect();
			aspect.inElev = dem4aspectcoverage;
			aspect.process();
			OmsRasterWriter.writeRaster("/Users/administrator/Documents/Adige_project/aspect_Avisio.asc", dem4aspectcoverage);
			
			
		

	}

}