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
public class Adige_DWM_setup extends HMTestCase {

	public void testCb() throws Exception {

		
		/*String path = new File(
			//	 "/Users/administrator/Documents/Adige_project/Adige_loc32632/SRTM_DEM90.tif"
				"/Users/administrator/Downloads/20150915104540_1337763757.tif"
				) 
				 .getAbsolutePath();
				 RasterReader reader = new RasterReader();
					reader.file = path;
					reader.fileNovalue = -32768.0;
					reader.geodataNovalue = Double.NaN;
					reader.process();
					GridCoverage2D dem90 = reader.outRaster;
					
					OmsRasterResolutionResampler resampler = new OmsRasterResolutionResampler();
					resampler.inGeodata = dem90;
					resampler.pInterpolation = "NEAREST_NEIGHTBOUR";
					resampler.pXres = 200.0;
					resampler.pYres = 200.0;
					resampler.process();
					GridCoverage2D dem150 = resampler.outGeodata;
					
					OmsRasterWriter.writeRaster("/Users/administrator/Documents/Adige_project/srtm_150.asc", dem150);
					*/
					
					
					
					
					
					

		/*	OmsRasterReprojector rasterreproject = new OmsRasterReprojector();
			rasterreproject.inRaster = dem30;
			rasterreproject.pCode = "EPSG:3003";
			rasterreproject.pEast = 1764836.690743859;
			rasterreproject.pWest = 1575029.7252891343;
			rasterreproject.pNorth = 5238637.317618476;
			rasterreproject.pSouth = 5020234.09290608;
			rasterreproject.pInterpolation = "INTERP_NEAREST";
			rasterreproject.process();
			GridCoverage2D dem30_reprojected = rasterreproject.outRaster;*/
			
			String path2 = new File(
					// "/Users/administrator/Documents/Adige_project/adige_location4258/adige_EUdem/cell/EU_DEM_inProj32632.asc"
					 "/Users/administrator/Documents/Adige_project/adige_location4258/adige_EUdem/cell/EU_DEM_inProj32632_2.asc") 
					 .getAbsolutePath();
					 RasterReader reader2 = new RasterReader();
						reader2.file = path2;
						reader2.fileNovalue = -9999.0;
						reader2.geodataNovalue = Double.NaN;
						reader2.process();
						GridCoverage2D EUdem_projected2_32632 = reader2.outRaster;
						OmsRasterWriter.writeRaster("/Users/administrator/Documents/Adige_project/EUdem_projected2_32632.asc", EUdem_projected2_32632);
						
		     		  /* OmsShapefileFeatureReader shapefileReader = new OmsShapefileFeatureReader();
						shapefileReader.file = "/Users/administrator/Documents/Adige_project/Adige_loc32632/maskshape4cutout.shp";
						shapefileReader.readFeatureCollection();
						SimpleFeatureCollection shapefileReaderFC = shapefileReader.geodata;
						
						OmsRasterVectorIntersector rasvecInter = new OmsRasterVectorIntersector();
						rasvecInter.inVector = shapefileReaderFC;
						rasvecInter.inRaster = EUdem_projected2_32632;
						rasvecInter.process();
						GridCoverage2D rastercutted = rasvecInter.outRaster;
						OmsRasterWriter.writeRaster("/Users/administrator/Documents/Adige_project/eudem_cutted.asc", rastercutted);
						
						*/
	
			
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
			extractnet.pThres = 10;
			extractnet.pMode = Variables.TCA;
			//extractnet.pExp = 0.5;
			extractnet.process();
			GridCoverage2D networkcoverage = extractnet.outNet;
		//	OmsRasterWriter.writeRaster("/Users/administrator/Documents/Adige_project/net_A50.asc", networkcoverage);
			
			
			OmsWateroutlet outlet = new OmsWateroutlet();
			outlet.inFlow = markoutcoverage;
			outlet.pEast =  676459.6713;
			outlet.pNorth = 5018091.296;
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
			OmsRasterWriter.writeRaster("/Users/administrator/Documents/Adige_project/net_A10_2.asc", networkcoverage2);
			
			
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
			OmsShapefileFeatureWriter.writeShapefile("/Users/administrator/Documents/Adige_project/net_atribute_A10.shp", netbiuldercoverage, pm);
			
			
			OmsExtractBasin extractbasin = new OmsExtractBasin();
			extractbasin.inFlow = markoutcoverage;
			extractbasin.inNetwork = netbiuldercoverage;
			extractbasin.pEast = 676459.6713;
			extractbasin.pNorth = 5018091.296;
			extractbasin.pSnapbuffer = 1000;
			extractbasin.process();
			
			OmsShapefileFeatureWriter.writeShapefile("/Users/administrator/Documents/Adige_project/maskshape.shp", extractbasin.outVectorBasin, pm);
			
					
			
	       
			OmsNetNumbering netnumber = new OmsNetNumbering();
			netnumber.inFlow = markoutcoverage2;
			netnumber.inNet = networkcoverage2;
			//netnumber.inTca = tcacoverage;
			//netnumber.pThres = 10000;
			netnumber.process();
			
			GridCoverage2D netnumercoverage = netnumber.outNetnum;
			OmsRasterWriter.writeRaster("/Users/administrator/Documents/Adige_project/netnumber_A10.asc", netnumercoverage);
			GridCoverage2D subbasincoverage = netnumber.outBasins;
			OmsRasterWriter.writeRaster("/Users/administrator/Documents/Adige_project/subbasin_A10.asc", subbasincoverage);
			
			
		
			
			OmsVectorizer vectorizer = new OmsVectorizer();
			vectorizer.inRaster = netnumercoverage;
			vectorizer.process();
		    OmsShapefileFeatureWriter.writeShapefile("/Users/administrator/Documents/Adige_project/netnumshapefeature_A10.shp", vectorizer.outVector, pm);
			
			
			OmsBasinShape basinshape = new OmsBasinShape();
			basinshape.inElev = pitcoverage2;
			basinshape.inBasins = subbasincoverage;
			basinshape.process();
			
			SimpleFeatureCollection basinshapefeature = basinshape.outBasins;
			OmsShapefileFeatureWriter.writeShapefile("/Users/administrator/Documents/Adige_project/basinshape_A10.shp", basinshapefeature, pm);
			
	
			OmsCutOut aspectcutout = new OmsCutOut();
			aspectcutout.inMask = mybasin;
			aspectcutout.inRaster = EUdem_projected2_32632;
			aspectcutout.process();
			GridCoverage2D dem4aspectcoverage = aspectcutout.outRaster;
			
			OmsAspect aspect = new OmsAspect();
			aspect.inElev = dem4aspectcoverage;
			aspect.process();
			//OmsRasterWriter.writeRaster("/Users/administrator/Documents/Adige_project/aspect_test.asc", dem4aspectcoverage);
			
			
		

	}

}