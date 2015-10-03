package org.osgeo.grass.r;

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

@Description("GRASS module to produce a raster map layer of gaussian deviates whose mean and standard deviation can be expressed by the user. It uses a gaussian random number generator.")
@Author(name = "Grass Developers Community", contact = "http://grass.osgeo.org")
@Keywords("raster")
@Label("Grass/Raster Modules")
@Name("r__surf__gauss")
@Status(Status.CERTIFIED)
@License("General Public License Version >=2)")
public class r__surf__gauss {

	@UI("outfile,grassfile")
	@Description("Name of the output random surface")
	@In
	public String $$outputPARAMETER;

	@Description("Distribution mean (optional)")
	@In
	public String $$meanPARAMETER = "0.0";

	@Description("Standard deviation (optional)")
	@In
	public String $$sigmaPARAMETER = "1.0";

	@Description("Allow output files to overwrite existing files")
	@In
	public boolean $$overwriteFLAG = false;

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
