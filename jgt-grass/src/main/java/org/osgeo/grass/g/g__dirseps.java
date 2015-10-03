package org.osgeo.grass.g;

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

@Description("Converts any directory separator characters in the input string to or from native host format, and writes the changed path to standard output. Useful in scripts for Windows compatibility.")
@Author(name = "Grass Developers Community", contact = "http://grass.osgeo.org")
@Keywords("general")
@Label("Grass/General Modules")
@Name("g__dirseps")
@Status(Status.CERTIFIED)
@License("General Public License Version >=2)")
public class g__dirseps {

	@Description("Path to be converted (read from stdin if not specified) (optional)")
	@In
	public String $$pathPARAMETER;

	@Description("Convert directory separators to native host format")
	@In
	public boolean $$hFLAG = false;

	@Description("Convert directory separators to GRASS internal format")
	@In
	public boolean $$gFLAG = false;

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
