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
import org.jgrasstools.hortonmachine.modules.demmanipulation.pitfiller.OmsPitfiller;
import org.jgrasstools.hortonmachine.modules.demmanipulation.wateroutlet.OmsWateroutlet;
import org.jgrasstools.hortonmachine.modules.geomorphology.aspect.OmsAspect;
import org.jgrasstools.hortonmachine.modules.geomorphology.draindir.OmsDrainDir;
import org.jgrasstools.hortonmachine.modules.geomorphology.flow.OmsFlowDirections;
import org.jgrasstools.hortonmachine.modules.geomorphology.slope.OmsSlope;
import org.jgrasstools.hortonmachine.modules.network.extractnetwork.OmsExtractNetwork;
import org.jgrasstools.hortonmachine.modules.network.netnumbering.OmsNetNumbering;
import org.jgrasstools.hortonmachine.modules.network.networkattributes.OmsNetworkAttributesBuilder;
import org.jgrasstools.hortonmachine.utils.HMTestCase;


public class UBN_DWM_set_up extends HMTestCase {

	public void testCb() throws Exception {
		
		String path = new File(
				// "/Users/administrator/Documents/PHDResearch/UBN_Geomorphometry/dem90.asc"
				"/Users/administrator/Dropbox/BLUENILE_project/20150909142639_1688175544.tif") 
				 .getAbsolutePath();
				 RasterReader reader = new RasterReader();
					reader.file = path;
					reader.fileNovalue = -32768.0;
					reader.geodataNovalue = Double.NaN;
					reader.process();
					GridCoverage2D dem90 = reader.outRaster;

   String path2 = new File(
		 "/Users/administrator/Documents/PHDResearch/UBN_Geomorphometry/UBN_dem_90.asc") //path to basin number 
		 .getAbsolutePath();
		 RasterReader reader2 = new RasterReader();
			reader2.file = path2;
			reader2.fileNovalue = -9999.0;
			reader2.geodataNovalue = Double.NaN;
			reader2.process();
			GridCoverage2D dem500 = reader2.outRaster;
			
			
			OmsPitfiller pit = new OmsPitfiller();
			pit.inElev = dem500;
			pit.process();
			GridCoverage2D pitcoverage = pit.outPit;
			
			
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
			
			OmsMarkoutlets markout = new OmsMarkoutlets();
			markout.inFlow = draindircoverage;
			markout.process();
			GridCoverage2D markoutcoverage = markout.outFlow;
			
			
			OmsExtractNetwork extractnet = new OmsExtractNetwork();
			extractnet.inFlow = markoutcoverage;
			extractnet.inTca = tcacoverage;
			//extractnet.inSlope = slopecoverage;
			extractnet.pThres = 1500;
			extractnet.pMode = Variables.TCA;
			//extractnet.pExp = 0.5;
			extractnet.process();
			GridCoverage2D networkcoverage = extractnet.outNet;
		
			
			OmsWateroutlet outlet = new OmsWateroutlet();
			outlet.inFlow = markoutcoverage;
			outlet.pEast = 711427.4284;
			outlet.pNorth = 1242179.6221;
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
			OmsRasterWriter.writeRaster("/Users/administrator/Documents/PHDResearch/UBN_Geomorphometry/net_test.asc", networkcoverage2);
			
			
			OmsCutOut tcacutout = new OmsCutOut();
			tcacutout.inMask = mybasin;
			tcacutout.inRaster = tcacoverage;
			tcacutout.process();
			GridCoverage2D tcacoverage2 = tcacutout.outRaster;
			OmsRasterWriter.writeRaster("/Users/administrator/Documents/PHDResearch/UBN_Geomorphometry/tca_test.asc", tcacutout.outRaster);
			
		
			
			
			OmsNetworkAttributesBuilder netbuilder = new OmsNetworkAttributesBuilder();
			netbuilder.inDem = pitcoverage2;
			netbuilder.inFlow = markoutcoverage2;
			netbuilder.inNet = networkcoverage2;
			netbuilder.inTca = tcacoverage2;
			netbuilder.process();
			
			GridCoverage2D hacknetcoverage = netbuilder.outHack;;
			
			SimpleFeatureCollection netbiuldercoverage = netbuilder.outNet;
			OmsShapefileFeatureWriter.writeShapefile("/Users/administrator/Documents/PHDResearch/UBN_Geomorphometry/net_atribute_test.shp", netbiuldercoverage, pm);
			
			
			
	       
			OmsNetNumbering netnumber = new OmsNetNumbering();
			netnumber.inFlow = markoutcoverage2;
			netnumber.inNet = networkcoverage2;
			//netnumber.inTca = tcacoverage;
			//netnumber.pThres = 10000;
			netnumber.process();
			
			GridCoverage2D netnumercoverage = netnumber.outNetnum;
			OmsRasterWriter.writeRaster("/Users/administrator/Documents/PHDResearch/UBN_Geomorphometry/netnumber_test.asc", netnumercoverage);
			GridCoverage2D subbasincoverage = netnumber.outBasins;
			OmsRasterWriter.writeRaster("/Users/administrator/Documents/PHDResearch/UBN_Geomorphometry/subbasin_test.asc", subbasincoverage);
			
			
		
			
			OmsVectorizer vectorizer = new OmsVectorizer();
			vectorizer.inRaster = netnumercoverage;
			vectorizer.process();
		    OmsShapefileFeatureWriter.writeShapefile("/Users/administrator/Documents/PHDResearch/UBN_Geomorphometry/netnumshapefeature_test.shp", vectorizer.outVector, pm);
			
			
			OmsBasinShape basinshape = new OmsBasinShape();
			basinshape.inElev = pitcoverage2;
			basinshape.inBasins = subbasincoverage;
			basinshape.process();
			
			SimpleFeatureCollection basinshapefeature = basinshape.outBasins;
			OmsShapefileFeatureWriter.writeShapefile("/Users/administrator/Documents/PHDResearch/UBN_Geomorphometry/basinshape_test.shp", basinshapefeature, pm);
			
			OmsCutOut aspectcutout = new OmsCutOut();
			aspectcutout.inMask = mybasin;
			aspectcutout.inRaster = dem90;
			aspectcutout.process();
			GridCoverage2D dem4aspectcoverage = aspectcutout.outRaster;
			
			OmsAspect aspect = new OmsAspect();
			aspect.inElev = dem4aspectcoverage;
			aspect.process();
			OmsRasterWriter.writeRaster("/Users/administrator/Documents/PHDResearch/UBN_Geomorphometry/aspect_test.asc", dem4aspectcoverage);
			
			
		

	}

}