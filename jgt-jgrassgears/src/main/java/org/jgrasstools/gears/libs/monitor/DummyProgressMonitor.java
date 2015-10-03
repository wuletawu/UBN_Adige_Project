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
package org.jgrasstools.gears.libs.monitor;

/**
 * As the name says.
 * 
 * @author Andrea Antonello (www.hydrologis.com)
 */
public class DummyProgressMonitor implements IJGTProgressMonitor {

    public void beginTask( String name, int totalWork ) {
    }

    public void beginTask( String name ) {
    }

    public void done() {
    }

    public void internalWorked( double work ) {
    }

    public boolean isCanceled() {
        return false;
    }

    public void setCanceled( boolean value ) {
    }

    public void setTaskName( String name ) {
    }

    public void subTask( String name ) {
    }

    public void worked( int work ) {
    }

    public <T> T adapt( Class<T> adaptee ) {
        return null;
    }

    public void errorMessage( String message ) {
    }

    public void message( String message ) {
    }

    public void exceptionThrown(String message) {
    }

    public void onModuleExit() {
    }

}
