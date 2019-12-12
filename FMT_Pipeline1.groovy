#@ MoleculeArchive archive
#@ Double (value=0.05) stuckRevThreshold
#@ Double (value=0.2) stuckMSDThreshold
#@ Integer (value=4500) minLength
#@ Double (value=0.0009) doubleCoilingLowerBound
#@ Double (value=10) doubleCoilingUpperBound
#@ Double (value=0.5) coilable20LowerBound
#@ Double (value=5) coilable20UpperBound
#@ Double (value=0.2) enzymaticLowerBound
#@ Double (value=6) enzymaticUpperBound
#@ Double (value=1) turns_per_second
#@ Double (value=2) turns_per_cycle
#@ Integer (value=1) driftZeroRegionStart
#@ Integer (value=100) driftZeroRegionEnd
#@ Double (value=Math.pow(10,-6)*1.56) conversionPixelToMicron
#@ Double (value=296.15) temperature
#@ Double (value=35*Math.pow(10,-9)) persistenceLength
#@ Double (value=6.8*Math.pow(10,-6)) contourLength
#@ ImageJ ij
#@ LogService logService

import de.mpg.biochem.mars.molecule.*;
import de.mpg.biochem.mars.table.*;
import de.mpg.biochem.mars.util.*;
import de.mpg.biochem.mars.molecule.commands.*

archive.lock()

//test
//BUILD LOG
String titleBlock = LogBuilder.buildTitleBlock("FMT Pipeline 1 - Start")
logService.info(titleBlock)
archive.addLogMessage(titleBlock)

logger = new LogBuilder()
logger.addParameter("FMT Pipeline 1 Verion", 2)
logger.addParameter("stuckRevThreshold", stuckRevThreshold)
logger.addParameter("stuckMSDThreshold", stuckMSDThreshold)
logger.addParameter("minLength", minLength)
logger.addParameter("doubleCoilingLowerBound", doubleCoilingLowerBound)
logger.addParameter("doubleCoilingUpperBound", doubleCoilingUpperBound)
logger.addParameter("enzymaticLowerBound", enzymaticLowerBound)
logger.addParameter("enzymaticUpperBound", enzymaticUpperBound)
logger.addParameter("turns_per_second", turns_per_second)
logger.addParameter("turns_per_cycle", turns_per_cycle)
logger.addParameter("driftZeroRegionStart", driftZeroRegionStart)
logger.addParameter("conversionPixelToMicron", conversionPixelToMicron)
logger.addParameter("temperature", temperature)
logger.addParameter("persistenceLength", persistenceLength)
logger.addParameter("contourLength", contourLength)

//ADD PARAMETERS HERE
//logger.addParameter("name", value)

String parameterList = logger.buildParameterList();
logService.info(parameterList)
archive.addLogMessage(parameterList)

//DONE BUILDING LOG

//Check to make sure there is only one metadata item
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

coil20_Positive_Peak = metadata.getRegion("coil20 Positive Peak")
coil20_Negative_Peak = metadata.getRegion("coil20 Negative Peak")

coil2p5_Peak = metadata.getRegion("coil2p5 Peak")
coil2p5_Background = metadata.getRegion("coil2p5 Background")

First_Reversal_RF = metadata.getRegion("First Reversal RF")
First_Reversal_FF = metadata.getRegion("First Reversal FF")

Last_Reversal_RF = metadata.getRegion("Last Reversal RF")
Last_Reversal_FF = metadata.getRegion("Last Reversal FF")

Negative_Coiling_Slope = metadata.getRegion("Negative Coiling Slope")
Positive_Coiling_Slope = metadata.getRegion("Positive Coiling Slope")

Before_Enzyme = metadata.getRegion("Before Enzyme")
After_Enzyme = metadata.getRegion("After Enzyme")

Force2p5 = metadata.getRegion("Force2p5")

archive.getMoleculeUIDs().parallelStream().forEach({ UID ->
	Molecule molecule = archive.get(UID)
	MarsTable table = molecule.getDataTable()

	//RegionDifferenceCalculatorCommand.calcRegionDifference(molecule, xColumn, yColumn, regionOne, regionTwo, parameterName)
    //coil20
	RegionDifferenceCalculatorCommand.calcRegionDifference(molecule, "slice", "x", coil20_Positive_Peak, coil20_Negative_Peak, "coil20")

	//coil2p5
	RegionDifferenceCalculatorCommand.calcRegionDifference(molecule, "slice", "x", coil2p5_Peak, coil2p5_Background, "coil2p5")

	//First Reversal
	RegionDifferenceCalculatorCommand.calcRegionDifference(molecule, "slice", "x", First_Reversal_RF, First_Reversal_FF, "rev_begin_x")
	RegionDifferenceCalculatorCommand.calcRegionDifference(molecule, "slice", "y", First_Reversal_RF, First_Reversal_FF, "rev_begin_y")

	//Last Reversal
	RegionDifferenceCalculatorCommand.calcRegionDifference(molecule, "slice", "x", Last_Reversal_RF, Last_Reversal_FF, "rev_end_x")
	RegionDifferenceCalculatorCommand.calcRegionDifference(molecule, "slice", "y", Last_Reversal_RF, Last_Reversal_FF, "rev_end_y")

	//Enzymatic Activity Detection
	RegionDifferenceCalculatorCommand.calcRegionDifference(molecule, "slice", "x", Before_Enzyme, After_Enzyme, "enzymatic")

	//MSDCalculatorCommand.calcMSD(molecule, column, parameterName)
	//MSD
	MSDCalculatorCommand.calcMSD(molecule, "x", "x_MSD")
	MSDCalculatorCommand.calcMSD(molecule, "y", "y_MSD")

	//Calculate slopes
	output = table.linearRegression("slice","x", Negative_Coiling_Slope.getStart(), Negative_Coiling_Slope.getEnd())
	molecule.setParameter("Slope_Neg_20", output[2])

    //Add Tags
    //stuckRev
    if (Math.abs(molecule.getParameter("rev_begin_y")) < stuckRevThreshold &&\
		Math.abs(molecule.getParameter("rev_end_y")) < stuckRevThreshold &&\
		Math.abs(molecule.getParameter("rev_begin_x")) < stuckRevThreshold &&\
		Math.abs(molecule.getParameter("rev_end_x")) < stuckRevThreshold) {
			molecule.addTag("stuckRev")
		}

	//stuckMSD
	if (molecule.getParameter("x_MSD") < stuckMSDThreshold &&\
		molecule.getParameter("y_MSD") < stuckMSDThreshold) {
			molecule.addTag("stuckMSD")
	}

	//coilable20
	if (molecule.getParameter("coil20") > coilable20LowerBound &&\
		molecule.getParameter("coil20") < coilable20UpperBound) {
			molecule.addTag("coilable20")
	}

	//doubleCoiling
	if (molecule.getParameter("Slope_Neg_20") > doubleCoilingLowerBound &&\
		molecule.getParameter("Slope_Neg_20") < doubleCoilingUpperBound) {
			molecule.addTag("doubleCoiling")
	}

	//shortTrajectory
	int deadSlice = table.getValue("slice",table.getRowCount()-1)
	if (deadSlice < minLength)
		molecule.addTag("shortTrajectory")

	//singleTether
	if(molecule.hasTag("coilable20") &&\
	  !molecule.hasTag("doubleCoiling") &&\
	  !molecule.hasTag("shortTrajectory")) {
	  	molecule.addTag("singleTether")
	}

	//enzymatic
	//if (molecule.hasTag("singleTether") &&\
	//	molecule.getParameter("enzymatic") > enzymaticLowerBound &&\
	//	molecule.getParameter("enzymatic") < enzymaticUpperBound) {
	//		molecule.addTag("chi")
	//}

	archive.put(molecule)
})

//Finish Log
logService.info(LogBuilder.endBlock())
archive.addLogMessage(LogBuilder.endBlock())
archive.unlock()

sleep(1000)

//Drift Calculator
final DriftCalculatorCommand driftCalc = new DriftCalculatorCommand();
driftCalc.setContext(ij.getContext());

//Set all the input parameters
driftCalc.setArchive(archive);
driftCalc.setBackgroundTag("stuckMSD");
driftCalc.setInputX("x");
driftCalc.setInputY("y");
driftCalc.setOutputX("x_drift");
driftCalc.setOutputY("y_drift");
driftCalc.setMode("mean")
driftCalc.setUseIncompleteTraces(false);

driftCalc.run();

sleep(1000)

//Drift Corrector
final DriftCorrectorCommand driftCorr = new DriftCorrectorCommand();
driftCorr.setContext(ij.getContext())

//Set all the input parameters
driftCorr.setArchive(archive)
driftCorr.setFromSlice(driftZeroRegionStart)
driftCorr.setToSlice(driftZeroRegionEnd)
driftCorr.setMetaX("x_drift")
driftCorr.setMetaY("y_drift")
driftCorr.setInputX("x")
driftCorr.setInputY("y")
driftCorr.setOutputX("x_drift_corr")
driftCorr.setOutputY("y_drift_corr")

driftCorr.run()

sleep(1000)

archive.lock()

String titleBlock2 = LogBuilder.buildTitleBlock("FMT Pipeline 1 - End")
logService.info(titleBlock2)
archive.addLogMessage(titleBlock2)

//Force Calculation
archive.getMoleculeUIDs().parallelStream().forEach({ UID ->
      Molecule molecule = archive.get(UID)
      MarsTable table = molecule.getDataTable()

	 double msd = table.msd("y_drift_corr", "slice", Force2p5.getStart(), Force2p5.getEnd());   //inserts different start and stop slice with each loop

	 msd = msd*conversionPixelToMicron*conversionPixelToMicron
	 double[] solution;
	 try {
	 	solution = MarsMath.calculateForceAndLength(persistenceLength, contourLength, temperature, msd);
	 } catch (Exception e) {
	 	return;
	 }

	 double force = solution[0]
	 double length = solution[1]
	 String force1 = "Force_PL35";			//make Force_2.5, Force_5 for different flow rates
	 String length1 = "Length_PL35";			// similarly for length

	 molecule.setParameter(force1, force);			// !!!! which one is correct? how can i insert a string here?
	 molecule.setParameter(length1, length);		// do i have to add something like str(force1)?

      archive.put(molecule)
 })


//Calculate poscycles and negcycles on a molecules-by-molecule basis
archive.getMoleculeUIDs().parallelStream().forEach({ UID ->
      Molecule molecule = archive.get(UID)
      MarsTable table = molecule.getDataTable()

	  double[] pos_coils_per_second = table.linearRegression("Time (s)", "x_drift_corr", Positive_Coiling_Slope.getStart(), Positive_Coiling_Slope.getEnd())

	  molecule.setParameter("pos_coil_slope", pos_coils_per_second[2])

	  	if (!table.hasColumn("poscycles"))
			table.appendColumn("poscycles")

		for (int row = 0; row < table.getRowCount(); row++) {
			double conversion = (-1)*(turns_per_second/pos_coils_per_second[2])/turns_per_cycle
			table.setValue("poscycles", row, conversion*table.getValue("x_drift_corr", row))
		}

	  double[] neg_coils_per_second = table.linearRegression("Time (s)", "x_drift_corr", Negative_Coiling_Slope.getStart(), Negative_Coiling_Slope.getEnd())

	  molecule.setParameter("neg_coil_slope", neg_coils_per_second[2])

	  	if (!table.hasColumn("negcycles"))
			table.appendColumn("negcycles")

		for (int row = 0; row < table.getRowCount(); row++) {
			double conversion = (1)*(turns_per_second/neg_coils_per_second[2])/turns_per_cycle
			table.setValue("negcycles", row, conversion*table.getValue("x_drift_corr", row))
		}

      archive.put(molecule)
})

//Calculate global median values for pos_coil_slope and neg_coil_slope

posCoilTable = new MarsTable("PosCoilRates","poscoil")
negCoilTable = new MarsTable("NegCoilRates","negcoil")

archive.getMoleculeUIDs().stream()\
.filter{ UID -> archive.moleculeHasTag(UID, "singleTether")}\
//.filter{ UID -> archive.moleculeHasTag(UID, "coilable2p5")}\
.forEach{ UID ->
     Molecule molecule = archive.get(UID)

     double pSlope = molecule.getParameter("pos_coil_slope")
     double nSlope = molecule.getParameter("neg_coil_slope")

	if (pSlope != Double.NaN) {
		posCoilTable.appendRow()
		posCoilTable.setValue("poscoil",posCoilTable.getRowCount()-1, pSlope)
	}

	if (nSlope != Double.NaN) {
	    negCoilTable.appendRow()
		negCoilTable.setValue("negcoil",negCoilTable.getRowCount()-1, nSlope)
	}
}

posCoilGlobal = posCoilTable.median("poscoil")
negCoilGlobal = negCoilTable.median("negcoil")

println("median poscoil slope global " + posCoilGlobal)
println("median negcoil slope global " + negCoilGlobal)

//Use global pos_coil_slope and neg_coil_slope to calculate poscycles and negcycles based on global averages.

archive.getMoleculeUIDs().parallelStream().forEach({ UID ->
	      Molecule molecule = archive.get(UID)
	      MarsTable table = molecule.getDataTable()

		  	if (!table.hasColumn("negcyclesG"))
				table.appendColumn("negcyclesG")

			for (int row = 0; row < table.getRowCount(); row++) {
				double conversion = (1)*(turns_per_second/negCoilGlobal)/turns_per_cycle
				table.setValue("negcyclesG", row, conversion*table.getValue("x_drift_corr", row))
			}

			if (!table.hasColumn("poscyclesG"))
				table.appendColumn("poscyclesG")

			for (int row = 0; row < table.getRowCount(); row++) {
				double conversion = (-1)*(turns_per_second/posCoilGlobal)/turns_per_cycle
				table.setValue("poscyclesG", row, conversion*table.getValue("x_drift_corr", row))
			}

	      archive.put(molecule)
	})

logService.info(LogBuilder.endBlock(true))
archive.addLogMessage(LogBuilder.endBlock(true))

archive.unlock()
