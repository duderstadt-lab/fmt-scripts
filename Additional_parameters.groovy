@ MoleculeArchive archive
#@ Double (value=2) torqueRecoveryLowerBound
#@ Double (value=6) torqueRecoveryUpperBound
#@ Double (value=1) coilable2p5LowerBound
#@ Double (value=5) coilable2p5UpperBound
#@ Double (value=-5) alphaLowerBound
#@ Double (value=-0.1) alphaUpperBound
#@ Double (value=0.2) chiLowerBound
#@ Double (value=5) chiUpperBound
#@ ImageJ ij
#@ LogService logService

import de.mpg.biochem.mars.molecule.*;
import de.mpg.biochem.mars.table.*;
import de.mpg.biochem.mars.util.*;
import de.mpg.biochem.mars.molecule.commands.*

archive.lock()

//Building log
String titleBlock = LogBuilder.buildTitleBlock("Additional analysis")
logService.info(titleBlock)
archive.addLogMessage(titleBlock)

logger = new LogBuilder()
logger.addParameter("Additional parameters run", 1)
logger.addParameter("torqueRecoveryLowerBound", torqueRecoveryLowerBound)
logger.addParameter("torqueRecoveryUpperBound", torqueRecoveryUpperBound)
logger.addParameter("coilable2p5LowerBound", coilable2p5LowerBound)
logger.addParameter("coilable2p5UpperBound", coilable2p5UpperBound)
logger.addParameter("alphaLowerBound", alphaLowerBound)
logger.addParameter("alphaUpperBound", alphaUpperBound)
logger.addParameter("chiLowerBound", chiLowerBound)
logger.addParameter("chiUpperBound", chiUpperBound)

String parameterList = logger.buildParameterList();
logService.info(parameterList)
archive.addLogMessage(parameterList)

//dont know if checking the metadata is important
if (archive.getNumberOfImageMetadataRecords() > 1) {
	logService.info(logger.buildParameterList())
	logService.info("More than one MarsImageMetadata item found - aborting.")
	logService.info(LogBuilder.endBlock(false))

	archive.addLogMessage(logger.buildParameterList())
	archive.addLogMessage("More than one MarsImageMetadata item found - aborting.")
	archive.addLogMessage(LogBuilder.endBlock(false));

	archive.unlock()

	return;
}

//Get regions
MarsImageMetadata metadata = archive.getImageMetadata(0);

Before_Enzyme = metadata.getRegion("Before Enzyme")
After_Enzyme = metadata.getRegion("After Enzyme")
Before_Torque_Recovery = metadata.getRegion("Before Torque Recovery")
After_Torque_Recovery = metadata.getRegion("After Torque Recovery")
Torque_A = metadata.getRegion("TorqueA")
Torque_B = metadata.getRegion("TorqueB")
Torque_C = metadata.getRegion("TorqueC")
Torque_D = metadata.getRegion("TorqueD")
Torque_E = metadata.getRegion("TorqueE")
Torque_F = metadata.getRegion("TorqueF")
Torque_G = metadata.getRegion("TorqueG")
Torque_H = metadata.getRegion("TorqueH")
Torque_I = metadata.getRegion("TorqueI")
Torque_J = metadata.getRegion("TorqueJ")
Torque_K = metadata.getRegion("TorqueK")

archive.getMoleculeUIDs().parallelStream().forEach({ UID ->
	Molecule molecule = archive.get(UID)
	MarsTable table = molecule.getDataTable()

	//RegionDifferenceCalculatorCommand.calcRegionDifference(molecule, xColumn, yColumn, regionOne, regionTwo, parameterName)
    //coil20
	RegionDifferenceCalculatorCommand.calcRegionDifference(molecule, "slice", "x", After_Torque_Recovery, Before_Torque_Recovery, "recovery")
	RegionDifferenceCalculatorCommand.calcRegionDifference(molecule, "slice", "x", Before_Enzyme, After_Enzyme, "enzymatic")

	if (molecule.getParameter("coil2p5") > coilable2p5LowerBound &&\
		molecule.getParameter("coil2p5") < coilable2p5UpperBound) {
		molecule.addTag("coilable2p5")
	}

	if (molecule.getParameter("recovery") > torqueRecoveryLowerBound &&\
		molecule.getParameter("recovery") < torqueRecoveryUpperBound) {
		molecule.addTag("torqueRecovery")
	}
	
	if (molecule.hasTag("singleTether") &&\
		molecule.hasTag("coilable2p5") &&\
		molecule.getParameter("enzymatic") < alphaLowerBound &&\
		molecule.getParameter("enzymatic") > alphaUpperBound) {
			molecule.addTag("alpha")
	}
	
	if (molecule.hasTag("singleTether") &&\
		molecule.hasTag("coilable2p5") &&\
		molecule.getParameter("enzymatic") < chiLowerBound &&\
		molecule.getParameter("enzymatic") > chiUpperBound) {
			molecule.addTag("chi")
	}


	archive.put(molecule)
})

archive.unlock()
