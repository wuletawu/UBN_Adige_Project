package org.osgeo.grass.gis;

import org.jgrasstools.grass.utils.ModuleSupporter;

import oms3.annotations.Author;
import oms3.annotations.Documentation;
import oms3.annotations.Label;
import oms3.annotations.Description;
import oms3.annotations.Execute;
import oms3.annotations.In;
import oms3.annotations.UI;
import oms3.annotations.Keywords;
import oms3.annotations.License;
import oms3.annotations.Name;
import oms3.annotations.Out;
import oms3.annotations.Status;

@Description("GIS manager for GRASS")
@Author(name = "Grass Developers Community", contact = "http://grass.osgeo.org")
@Label("Grass")
@Name("gis__m")
@Status(Status.CERTIFIED)
@License("General Public License Version >=2)")
public class gis__m {

	@Description("Name of GIS manager settings file (.grc) (optional)")
	@In
	public String $$dmrcPARAMETER;

	@Description("Verbose module output")
	@In
	public boolean $$verboseFLAG = false;

	@Description("Quiet module output")
	@In
	public boolean $$quietFLAG = false;


	@Execute
	public void process() throws Exception {
		ModuleSupporter.processModule(this);
	}

}
