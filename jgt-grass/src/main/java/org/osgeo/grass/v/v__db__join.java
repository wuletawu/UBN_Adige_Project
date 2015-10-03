package org.osgeo.grass.v;

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

@Description("Allows to join a table to a vector map table.")
@Author(name = "Grass Developers Community", contact = "http://grass.osgeo.org")
@Keywords("vector, database, attribute table")
@Label("Grass/Vector Modules")
@Name("v__db__join")
@Status(Status.CERTIFIED)
@License("General Public License Version >=2)")
public class v__db__join {

	@UI("infile,grassfile")
	@Description("Vector map to which to join other table")
	@In
	public String $$mapPARAMETER;

	@Description("Layer where to join (optional)")
	@In
	public String $$layerPARAMETER = "1";

	@Description("Join column in map table")
	@In
	public String $$columnPARAMETER;

	@Description("Other table name")
	@In
	public String $$otablePARAMETER;

	@Description("Join column in other table")
	@In
	public String $$ocolumnPARAMETER;

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
