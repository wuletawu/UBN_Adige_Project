package org.jgrasstools.hortonmachine.modules.hydrogeomorphology.adige.core;

import java.io.PrintStream;
import java.util.*;

public class PfafstetterNumber
    implements Comparator {

    public PfafstetterNumber(String pfafstetterNumberString) {
        this.pfafstetterNumberString = null;
        pfafstetterUpToLastLeveL = null;
        this.order = -1;
        ordersList = null;
        this.pfafstetterNumberString = pfafstetterNumberString;
        ordersList = new ArrayList();
        int lastDot = pfafstetterNumberString.lastIndexOf('.');
        if(lastDot == -1) {
            this.order = 1;
            ordersList.add(Integer.valueOf(Integer.parseInt(pfafstetterNumberString)));
            pfafstetterUpToLastLeveL = "";
        } else {
            String order[] = pfafstetterNumberString.split("\\.");
            this.order = order.length;
            String as[];
            int j = (as = order).length;
            for(int i = 0; i < j; i++) {
                String string = as[i];
                ordersList.add(Integer.valueOf(Integer.parseInt(string)));
            }

            pfafstetterUpToLastLeveL = pfafstetterNumberString.substring(0, lastDot + 1);
        }
    }

    public int getOrder() {
        return order;
    }

    public List getOrdersList() {
        return ordersList;
    }

    public String toString() {
        return pfafstetterNumberString;
    }

    public String toStringUpToLastLevel() {
        return pfafstetterUpToLastLeveL;
    }

    public boolean isOfOrderOrMinor(int order) {
        return this.order >= order;
    }

    public boolean isDownStreamOf(PfafstetterNumber pfafstetterNumber) {
        int lastDot = pfafstetterNumberString.lastIndexOf('.');
        String pre = pfafstetterNumberString.substring(0, lastDot + 1);
        String lastNum = pfafstetterNumberString.substring(lastDot + 1, pfafstetterNumberString.length());
        int lastNumInt = Integer.parseInt(lastNum);
        if(lastNumInt % 2 == 0)
            return false;
        String pfaff = pfafstetterNumber.toString();
        if(pfaff.startsWith(pre)) {
            String lastPart = pfaff.substring(lastDot + 1, pfaff.length());
            if(Integer.parseInt(lastPart.split("\\.")[0]) >= lastNumInt)
                return true;
        }
        return false;
    }

    public boolean isEndPiece() {
        return ((Integer)ordersList.get(ordersList.size() - 1)).intValue() % 2 == 0;
    }

    public int compare(PfafstetterNumber p1, PfafstetterNumber p2) {
        List p1OrdersList = p1.getOrdersList();
        List p2OrdersList = p2.getOrdersList();
        int levels = p1OrdersList.size();
        if(p2OrdersList.size() < levels)
            levels = p2OrdersList.size();
        for(int i = 0; i < levels; i++) {
            int thisone = ((Integer)p1OrdersList.get(i)).intValue();
            int otherone = ((Integer)p2OrdersList.get(i)).intValue();
            if(thisone > otherone)
                return -1;
            if(thisone < otherone)
                return 1;
        }

        return 0;
    }

    public static synchronized boolean areConnectedUpstream(PfafstetterNumber p1, PfafstetterNumber p2) {
        List p1OrdersList = p1.getOrdersList();
        List p2OrdersList = p2.getOrdersList();
        int levelDiff = p1OrdersList.size() - p2OrdersList.size();
        if(levelDiff == 0) {
            if(p1.toStringUpToLastLevel().equals(p2.toStringUpToLastLevel())) {
                int p1Last = ((Integer)p1OrdersList.get(p1OrdersList.size() - 1)).intValue();
                int p2Last = ((Integer)p2OrdersList.get(p2OrdersList.size() - 1)).intValue();
                if(p2Last == p1Last + 1 || p2Last == p1Last + 2)
                    return p1Last % 2 != 0;
            }
        } else
        if(levelDiff == -1 && p2.toString().startsWith(p1.toStringUpToLastLevel())) {
            int p2Last = ((Integer)p2OrdersList.get(p2OrdersList.size() - 1)).intValue();
            if(p2Last != 1)
                return false;
            int p1Last = ((Integer)p1OrdersList.get(p1OrdersList.size() - 1)).intValue();
            int p2LastMinus1 = ((Integer)p2OrdersList.get(p2OrdersList.size() - 2)).intValue();
            if(p2LastMinus1 == p1Last + 1 || p2Last == p1Last + 2)
                return p1Last % 2 != 0;
        }
        return false;
    }

    public static synchronized boolean areConnectedDownstream(PfafstetterNumber p1, PfafstetterNumber p2) {
        return areConnectedUpstream(p2, p1);
    }

    public static void main(String args[]) {
        PfafstetterNumber n1 = new PfafstetterNumber("2.5");
        PfafstetterNumber n2 = new PfafstetterNumber("2.6.4");
        PfafstetterNumber n3 = new PfafstetterNumber("2.4.3");
        PfafstetterNumber n4 = new PfafstetterNumber("2.7.1");
        PfafstetterNumber n5 = new PfafstetterNumber("2.4.16.45");
        PfafstetterNumber n6 = new PfafstetterNumber("2.6.2.1");
        PfafstetterNumber n7 = new PfafstetterNumber("2.7.6.5.2");
        PfafstetterNumber n8 = new PfafstetterNumber("2.7.6.2.1");
        PfafstetterNumber n9 = new PfafstetterNumber("2.6.2.7");
        List list = new ArrayList();
        list.add(n1);
        list.add(n2);
        list.add(n3);
        list.add(n4);
        list.add(n5);
        list.add(n6);
        list.add(n7);
        list.add(n8);
        System.out.println(n1.isDownStreamOf(n2));
        System.out.println(n1.isDownStreamOf(n4));
        System.out.println(n3.isDownStreamOf(n2));
        System.out.println(n3.isDownStreamOf(n5));
        System.out.println(n4.isDownStreamOf(n1));
        System.out.println(n6.isDownStreamOf(n7));
        System.out.println(n8.isDownStreamOf(n7));
        System.out.println();
        System.out.println(areConnectedUpstream(n1, n2));
        System.out.println(areConnectedUpstream(n1, n4));
        System.out.println(areConnectedUpstream(n3, n2));
        System.out.println(areConnectedUpstream(n3, n5));
        System.out.println(areConnectedUpstream(n4, n1));
        System.out.println(areConnectedUpstream(n4, n6));
        System.out.println(areConnectedUpstream(n8, n7));
        System.out.println(areConnectedUpstream(n6, n9));
        System.out.println();
        PfafstetterNumber array[] = (PfafstetterNumber[])list.toArray(new PfafstetterNumber[list.size()]);
        Arrays.sort(array, n1);
        PfafstetterNumber apfafstetternumber[];
        int j = (apfafstetternumber = array).length;
        for(int i = 0; i < j; i++) {
            PfafstetterNumber pfafstetterNumber = apfafstetternumber[i];
            System.out.println(pfafstetterNumber.toString());
        }

    }

    public int compare(Object obj, Object obj1) {
        return compare((PfafstetterNumber)obj, (PfafstetterNumber)obj1);
    }

    private String pfafstetterNumberString;
    private String pfafstetterUpToLastLeveL;
    private int order;
    private List ordersList;
}
