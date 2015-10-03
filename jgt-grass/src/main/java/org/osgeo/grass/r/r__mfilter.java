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

@Description("Raster map matrix filter.")
@Author(name = "Grass Developers Community", contact = "http://grass.osgeo.org")
@Keywords("raster, map algebra")
@Label("Grass/Raster Modules")
@Name("r__mfilter")
@Status(Status.CERTIFIED)
@License("General Public License Version >=2)")
public class r__mfilter {

	@UI("infile,grassfile")
	@Description("Name of input raster map")
	@In
	public String $$inputPARAMETER;

	@UI("outfile,grassfile")
	@Description("Name for output raster map")
	@In
	public String $$outputPARAMETER;

	@Description("Name of filter file")
	@In
	public String $$filterPARAMETER;

	@Description("Number of times to repeat the filter (optional)")
	@In
	public String $$repeatPARAMETER = "1";

	@Description("Output raster map title (optional)")
	@In
	public String $$titlePARAMETER;

	@Description("Quiet")
	@In
	public boolean $$qFLAG = false;

	@Description("Apply filter only to zero data values")
	@In
	public boolean $$zFLAG = false;

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
