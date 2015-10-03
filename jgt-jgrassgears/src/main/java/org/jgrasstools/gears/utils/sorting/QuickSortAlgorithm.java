/*
 * JGrass - Free Open Source Java GIS http://www.jgrass.org 
 * (C) HydroloGIS - www.hydrologis.com 
 * 
 * This library is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Library General Public License as published by the Free
 * Software Foundation; either version 2 of the License, or (at your option) any
 * later version.
 * 
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Library General Public License for more
 * details.
 * 
 * You should have received a copy of the GNU Library General Public License
 * along with this library; if not, write to the Free Foundation, Inc., 59
 * Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */
package org.jgrasstools.gears.utils.sorting;

import static org.jgrasstools.gears.libs.modules.JGTConstants.isNovalue;

import org.jgrasstools.gears.libs.monitor.LogProgressMonitor;
import org.jgrasstools.gears.libs.monitor.IJGTProgressMonitor;
public class QuickSortAlgorithm {
    private double[] valuesToSort;
    private int number;
    private double[] valuesToFollow;
    private IJGTProgressMonitor monitor = new LogProgressMonitor();

    public QuickSortAlgorithm( IJGTProgressMonitor monitor ) {
        if (monitor != null)
            this.monitor = monitor;
    }

    /**
     * Sorts an array of values and moves with the sort a second array.
     * 
     * @param values the array to sort.
     * @param valuesToFollow the array that should be sorted following the 
     *              indexes of the first array. Can be null.
     */
    public void sort( double[] values, double[] valuesToFollow ) {
        this.valuesToSort = values;
        this.valuesToFollow = valuesToFollow;

        number = values.length;

        monitor.beginTask("Sorting...", -1);

        monitor.worked(1);
        quicksort(0, number - 1);

        monitor.done();
    }

    private void quicksort( int low, int high ) {
        int i = low, j = high;
        // Get the pivot element from the middle of the list
        double pivot = valuesToSort[(low + high) >>> 1];

        // Divide into two lists
        while( i <= j ) {
            // If the current value from the left list is smaller then the pivot
            // element then get the next element from the left list
            while( valuesToSort[i] < pivot || (isNovalue(valuesToSort[i]) && !isNovalue(pivot)) ) {

                i++;
            }
            // If the current value from the right list is larger then the pivot
            // element then get the next element from the right list
            while( valuesToSort[j] > pivot || (!isNovalue(valuesToSort[j]) && isNovalue(pivot)) ) {
                j--;
            }

            // If we have found a values in the left list which is larger then
            // the pivot element and if we have found a value in the right list
            // which is smaller then the pivot element then we exchange the
            // values.
            // As we are done we can increase i and j
            if (i <= j) {
                exchange(i, j);
                i++;
                j--;
            }
        }
        // Recursion
        if (low < j)
            quicksort(low, j);
        if (i < high)
            quicksort(i, high);
    }

    private void exchange( int i, int j ) {
        double temp = valuesToSort[i];
        valuesToSort[i] = valuesToSort[j];
        valuesToSort[j] = temp;
        if (valuesToFollow != null) {
            double tempFollow = valuesToFollow[i];
            valuesToFollow[i] = valuesToFollow[j];
            valuesToFollow[j] = tempFollow;
        }
    }

}
