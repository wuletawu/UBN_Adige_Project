package org.jgrasstools.hortonmachine.modules.hydrogeomorphology.adige.utils;

import java.text.MessageFormat;
import java.util.*;
import org.geotools.feature.FeatureCollection;
import org.geotools.feature.FeatureIterator;
import org.jgrasstools.gears.libs.monitor.IJGTProgressMonitor;
import org.jgrasstools.hortonmachine.modules.hydrogeomorphology.adige.core.*;
import org.jgrasstools.hortonmachine.modules.network.PfafstetterNumber;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.feature.simple.SimpleFeatureType;

public class AdigeUtilities2 {

    public AdigeUtilities2() {
    }

    public static List generateHillSlopes(FeatureCollection netFeatureCollection, FeatureCollection hillslopeFeatureCollection, String netnumAttr, String pfafAttr, String startelevAttr, String endelevAttr, String baricenterAttr, IJGTProgressMonitor out)
        throws Exception {
        SimpleFeatureType fT = (SimpleFeatureType)netFeatureCollection.getSchema();
        int lAttrIndex = fT.indexOf(netnumAttr);
        if(lAttrIndex == -1) {
            String pattern = "Attribute {0} not found in layer {1}.";
            Object args[] = {
                netnumAttr, fT.getTypeName()
            };
            String newPattern = MessageFormat.format(pattern, args);
            throw new IllegalArgumentException(newPattern);
        }
        int pAttrIndex = fT.indexOf(pfafAttr);
        if(pAttrIndex == -1) {
            String pattern = "Attribute {0} not found in layer {1}.";
            Object args[] = {
                pfafAttr, fT.getTypeName()
            };
            String newPattern = MessageFormat.format(pattern, args);
            throw new IllegalArgumentException(newPattern);
        }
        int startNetElevAttrIndex = -1;
        if(startelevAttr != null) {
            startNetElevAttrIndex = fT.indexOf(startelevAttr);
            if(startNetElevAttrIndex == -1) {
                String pattern = "Attribute {0} not found in layer {1}.";
                Object args[] = {
                    startelevAttr, fT.getTypeName()
                };
                String newPattern = MessageFormat.format(pattern, args);
                throw new IllegalArgumentException(newPattern.getClass().getSimpleName());
            }
        }
        int endNetElevAttrIndex = -1;
        if(endelevAttr != null) {
            endNetElevAttrIndex = fT.indexOf(endelevAttr);
            if(endNetElevAttrIndex == -1) {
                String pattern = "Attribute {0} not found in layer {1}.";
                Object args[] = {
                    endelevAttr, fT.getTypeName()
                };
                String newPattern = MessageFormat.format(pattern, args);
                throw new IllegalArgumentException(newPattern);
            }
        }
        out.message("Analizing the network layer...");
        List netFeaturesList = new ArrayList();
        List netIdsList = new ArrayList();
        ArrayList netPfaffsList = new ArrayList();
        FeatureIterator featureIterator = netFeatureCollection.features();
        PfafstetterNumber mostDownStreamPNumber = null;
        SimpleFeature mostDownStreamNetFeature = null;
        Integer mostDownStreamLinkId = Integer.valueOf(-1);
        PfafstetterNumber current;
        for(; featureIterator.hasNext(); netPfaffsList.add(current)) {
            SimpleFeature f = (SimpleFeature)featureIterator.next();
            String attribute = (String)f.getAttribute(pAttrIndex);
            current = new PfafstetterNumber(attribute);
            Integer tmpId = Integer.valueOf(((Number)f.getAttribute(lAttrIndex)).intValue());
            if(mostDownStreamPNumber == null)
                mostDownStreamPNumber = current;
            else
            if(current.isDownStreamOf(mostDownStreamPNumber)) {
                mostDownStreamLinkId = tmpId;
                mostDownStreamNetFeature = f;
                mostDownStreamPNumber = current;
            }
            netFeaturesList.add(f);
            netIdsList.add(tmpId);
        }

        featureIterator.close();
        out.message("Analyzing the hillslopes layer...");
        SimpleFeatureType ft = (SimpleFeatureType)hillslopeFeatureCollection.getSchema();
        int linkAttrIndexInBasinLayerIndex = ft.indexOf(netnumAttr);
        if(linkAttrIndexInBasinLayerIndex == -1) {
            String pattern = "Attribute {0} not found in layer {1}.";
            Object args[] = {
                netnumAttr, ft.getTypeName()
            };
            pattern = MessageFormat.format(pattern, args);
            throw new IllegalArgumentException(pattern);
        }
        int baricenterAttributeIndex = -1;
        if(baricenterAttr != null) {
            baricenterAttributeIndex = ft.indexOf(baricenterAttr);
            if(baricenterAttributeIndex == -1) {
                String pattern = "Attribute {0} not found in layer {1}.";
                Object args[] = {
                    baricenterAttr, ft.getTypeName()
                };
                pattern = MessageFormat.format(pattern, args);
                throw new IllegalArgumentException(pattern);
            }
        }
        List hillslopeFeaturesList = new ArrayList();
        List hillslopeIdsList = new ArrayList();
        FeatureIterator hillslopeIterator = hillslopeFeatureCollection.features();
        SimpleFeature mostDownstreamHillslopeFeature = null;
        SimpleFeature f;
        for(; hillslopeIterator.hasNext(); hillslopeFeaturesList.add(f)) {
            f = (SimpleFeature)hillslopeIterator.next();
            Integer linkAttribute = Integer.valueOf(((Number)f.getAttribute(linkAttrIndexInBasinLayerIndex)).intValue());
            if(mostDownStreamLinkId == linkAttribute)
                mostDownstreamHillslopeFeature = f;
            hillslopeIdsList.add(linkAttribute);
        }

        out.message("Linking together network and hillslopes layers...");
        ArrayList hillslopeElements = new ArrayList();
        IHillSlope mostDownstreamHillslope = null;
        if(mostDownStreamPNumber.isEndPiece()) {
            Integer basinId = (Integer)hillslopeIdsList.get(0);
            IHillSlope tmpHslp = new HillSlope(mostDownStreamNetFeature, mostDownstreamHillslopeFeature, mostDownStreamPNumber, basinId.intValue());
            hillslopeElements.add(tmpHslp);
            mostDownstreamHillslope = tmpHslp;
        } else {
            ArrayList selectedNetFeatureList = new ArrayList();
            ArrayList selectedNetId = new ArrayList();
            for(int i = 0; i < hillslopeFeaturesList.size(); i++) {
                SimpleFeature basinFeature = (SimpleFeature)hillslopeFeaturesList.get(i);
                Integer link = (Integer)hillslopeIdsList.get(i);
                for(int j = 0; j < netFeaturesList.size(); j++) {
                    Integer netNum = (Integer)netIdsList.get(j);
                    if(!netNum.equals(link))
                        continue;
                    SimpleFeature netFeature = (SimpleFeature)netFeaturesList.get(j);
                    IHillSlope tmpHslp = new HillSlope(netFeature, basinFeature, (PfafstetterNumber)netPfaffsList.get(j), netNum.intValue());
                    hillslopeElements.add(tmpHslp);
                    selectedNetFeatureList.add(netFeature);
                    selectedNetId.add(netNum);
                    break;
                }

            }

            mostDownStreamPNumber = null;
            Integer mostDownStreamNetId = null;
            for(Iterator iterator = selectedNetFeatureList.iterator(); iterator.hasNext();) {
                SimpleFeature feature = (SimpleFeature)iterator.next();
                String attribute = (String)feature.getAttribute(pAttrIndex);
                PfafstetterNumber current2 = new PfafstetterNumber(attribute);
                Integer tmpId = Integer.valueOf(((Number)feature.getAttribute(lAttrIndex)).intValue());
                if(mostDownStreamPNumber == null)
                    mostDownStreamPNumber = current2;
                else
                if(current2.isDownStreamOf(mostDownStreamPNumber)) {
                    mostDownStreamNetId = tmpId;
                    mostDownStreamPNumber = current2;
                }
            }

            for(int i = 0; i < hillslopeElements.size(); i++) {
                Integer hId = (Integer)hillslopeIdsList.get(i);
                if(!hId.equals(mostDownStreamNetId))
                    continue;
                mostDownstreamHillslope = (IHillSlope)hillslopeElements.get(i);
                break;
            }

            if(hillslopeElements.size() == 1)
                mostDownstreamHillslope = (IHillSlope)hillslopeElements.get(0);
        }
        if(mostDownstreamHillslope == null) {
            throw new RuntimeException();
        } else {
            HillSlope.connectElements(hillslopeElements);
            List orderedHillslopes = new ArrayList();
            mostDownstreamHillslope.getAllUpstreamElements(orderedHillslopes, null);
            return orderedHillslopes;
        }
    }

    public static double doRouting(double discharge, IHillSlope hillslope, int routingType) {
        double linkWidth = hillslope.getLinkWidth(8.6600000000000001D, 0.59999999999999998D, 0.0D);
        double linkLength = hillslope.getLinkLength();
        double linkSlope = hillslope.getLinkSlope();
        double chezLawExpon = -0.33333333333333331D;
        double chezLawCoeff = 200D / Math.pow(0.00035791099999999998D, chezLawExpon);
        double linkChezy = hillslope.getLinkChezi(chezLawCoeff, chezLawExpon);
        double K_Q = 0.0D;
        switch(routingType) {
        case 2: // '\002'
            K_Q = 8.7959999999999994D * Math.pow(discharge, 0.33333333333333331D) * Math.pow(linkWidth, -0.33333333333333331D) * Math.pow(linkLength, -1D) * Math.pow(linkSlope, 0.22222222222222221D);
            break;

        case 3: // '\003'
            K_Q = 1.5D * Math.pow(discharge, 0.33333333333333331D) * Math.pow(linkChezy, 0.66666666666666663D) * Math.pow(linkWidth, -0.33333333333333331D) * Math.pow(linkLength, -1D) * Math.pow(linkSlope, 0.33333333333333331D);
            break;

        case 4: // '\004'
            double flowdepth = 0.33333333333333331D * Math.pow(discharge, 0.33333333333333331D);
            double hydrad = (flowdepth * linkWidth) / (2D * flowdepth + linkWidth);
            double mannings_n = 1.0D;
            K_Q = ((Math.pow(hydrad, 0.66666666666666663D) * Math.pow(linkSlope, 0.5D)) / mannings_n) * Math.pow(linkLength, -1D);
            break;
        }
        return K_Q;
    }
}
