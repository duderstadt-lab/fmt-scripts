@ MoleculeArchive archive
#@ Double (value=0.05) stuckRevThreshold
#@ Double (value=0.1) stuckMSDThresholdx
#@ Double (value=0.2) stuckMSDThresholdy
#@ Integer (value=4500) minLength
#@ Double (value=6) reversalLowerBound
#@ Double (value=12) reversalUpperBound
#@ Double (value=-10) reversalLowerBoundNeg
#@ Double (value=-0.5) reversalUpperBoundNeg
#@ Double (value=0.0009) doubleCoilingLowerBound
#@ Double (value=10) doubleCoilingUpperBound
#@ Double (value=0.5) coilable20LowerBound
#@ Double (value=5) coilable20UpperBound
#@ Double (value=1) coilable2p5LowerBound
#@ Double (value=5) coilable2p5UpperBound
#@ Double (value=-5) alphaLowerBound
#@ Double (value=-0.1) alphaUpperBound
#@ Double (value=0.2) chiLowerBound
#@ Double (value=5) chiUpperBound
#@ Double (value=2) torqueRecoveryLowerBound
#@ Double (value=6) torqueRecoveryUpperBound
#@ Double (value=0.25) turns_per_slice
#@ Double (value=2) turns_per_cycle
#@ Integer (value=1600) driftZeroRegionStart
#@ Integer (value=1700) driftZeroRegionEnd
#@ Double (value=1.56E-6) conversionPixelToMicron
#@ Double (value=296.15) temperature
#@ Double (value=3.5E-8) persistenceLength
#@ Double (value=6.8E-6) contourLength
#@ Integer (value=20) slidingForceWindow
#@ Double (value=1.5) stdSlidingForce
#@ ImageJ ij
#@ LogService logService

import de.mpg.biochem.mars.molecule.*;
import de.mpg.biochem.mars.table.*;
import de.mpg.biochem.mars.util.*;
import de.mpg.biochem.mars.molecule.commands.*
import org.scijava.table.*

archive.lock()

//test
//BUILD LOG
String titleBlock = LogBuilder.buildTitleBlock("FMT Pipeline 2 - Start")
logService.info(titleBlock)
archive.addLogMessage(titleBlock)

logger = new LogBuilder()
logger.addParameter("FMT Pipeline 2 Verion", 1)
logger.addParameter("stuckRevThreshold", stuckRevThreshold)
logger.addParameter("stuckMSDThresholdx", stuckMSDThresholdx)
logger.addParameter("stuckMSDThresholdy", stuckMSDThresholdy)
logger.addParameter("minLength", minLength)
logger.addParameter("doubleCoilingLowerBound", doubleCoilingLowerBound)
logger.addParameter("doubleCoilingUpperBound", doubleCoilingUpperBound)
logger.addParameter("torqueRecoveryLowerBound", torqueRecoveryLowerBound)
logger.addParameter("torqueRecoveryUpperBound", torqueRecoveryUpperBound)
logger.addParameter("coilable2p5LowerBound", coilable2p5LowerBound)
logger.addParameter("coilable2p5UpperBound", coilable2p5UpperBound)
logger.addParameter("alphaLowerBound", alphaLowerBound)
logger.addParameter("alphaUpperBound", alphaUpperBound)
logger.addParameter("chiLowerBound", chiLowerBound)
logger.addParameter("chiUpperBound", chiUpperBound)
logger.addParameter("coilable20LowerBound", coilable20LowerBound)
logger.addParameter("coilable20UpperBound", coilable20UpperBound)
logger.addParameter("turns_per_slice", turns_per_slice)
logger.addParameter("turns_per_cycle", turns_per_cycle)
logger.addParameter("driftZeroRegionStart", driftZeroRegionStart)
logger.addParameter("conversionPixelToMicron", conversionPixelToMicron)
logger.addParameter("temperature", temperature)
logger.addParameter("persistenceLength", persistenceLength)
logger.addParameter("contourLength", contourLength)
logger.addParameter("reversalLowerBound", reversalLowerBound)
logger.addParameter("reversalUpperBound", reversalUpperBound)
logger.addParameter("reversalLowerBoundNeg", reversalLowerBoundNeg)
logger.addParameter("reversalUpperBoundNeg", reversalUpperBoundNeg)
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

Slope_Neg_20 = metadata.getRegion("Slope_Neg_20")

Negative_Coiling_Slope = metadata.getRegion("Negative Coiling Slope")
Positive_Coiling_Slope = metadata.getRegion("Positive Coiling Slope")

Before_Enzyme = metadata.getRegion("Before Enzyme")
After_Enzyme = metadata.getRegion("After Enzyme")

Before_Torque_Recovery = metadata.getRegion("Before Torque Recovery")
After_Torque_Recovery = metadata.getRegion("After Torque Recovery")

Force2p5 = metadata.getRegion("Force2p5")

Gyrase_Reaction = metadata.getRegion("Gyrase Reaction")

Magrot_20f = metadata.getRegion("Magrot20f")

Magrot_2p5f = metadata.getRegion("Magrot2p5f")

Region_Torque = metadata.getRegion("regionTorque")

//ADD Time
logService.info("Adding time...")
AddTimeCommand.addTime(archive)

logService.info("Calculating region differences and adding tags...")
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

	//Torque Recovery
	RegionDifferenceCalculatorCommand.calcRegionDifference(molecule, "slice", "x", After_Torque_Recovery, Before_Torque_Recovery, "recovery")

	//MSDCalculatorCommand.calcMSD(molecule, column, parameterName)
	//MSD
	MSDCalculatorCommand.calcMSD(molecule, "x", "x_MSD")
	MSDCalculatorCommand.calcMSD(molecule, "y", "y_MSD")

	//Calculate slopes
	output = table.linearRegression("slice","x", Slope_Neg_20.getStart(), Slope_Neg_20.getEnd())
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
	if (molecule.hasTag("stuckRev") &&\
		molecule.getParameter("x_MSD") < stuckMSDThresholdx &&\
		molecule.getParameter("y_MSD") < stuckMSDThresholdy) {
			molecule.addTag("stuckMSD")
	}

	//coilable20
	if (molecule.getParameter("coil20") > coilable20LowerBound &&\
		molecule.getParameter("coil20") < coilable20UpperBound) {
			molecule.addTag("coilable20")
	}

	//coilable2p5
	if (molecule.getParameter("coil2p5") > coilable2p5LowerBound &&\
		molecule.getParameter("coil2p5") < coilable2p5UpperBound) {
			molecule.addTag("coilable2p5")
	}

	//doubleCoiling
	if (molecule.getParameter("Slope_Neg_20") > doubleCoilingLowerBound &&\
		molecule.getParameter("Slope_Neg_20") < doubleCoilingUpperBound) {
			molecule.addTag("doubleCoiling")
	}

	//alpha
	if (molecule.getParameter("enzymatic") > alphaLowerBound &&\
		molecule.getParameter("enzymatic") < alphaUpperBound) {
			molecule.addTag("alpha")
	}

	//chi
	if (molecule.getParameter("enzymatic") > chiLowerBound &&\
		molecule.getParameter("enzymatic") < chiUpperBound) {
			molecule.addTag("chimol")
	}
	//torqueRecovery
	if (molecule.getParameter("recovery") > torqueRecoveryLowerBound &&\
		molecule.getParameter("recovery") < torqueRecoveryUpperBound) {
			molecule.addTag("torqueRecovery")
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

	//reversal start
	if (molecule.getParameter("rev_begin_x") > reversalLowerBound &&\
		molecule.getParameter("rev_begin_x") < reversalUpperBound) {
		molecule.addTag("mobile_begin")
	}

	//reversal end
	if (molecule.getParameter("rev_end_x") > reversalLowerBound &&\
		molecule.getParameter("rev_end_x") < reversalUpperBound) {
		molecule.addTag("mobile_end")
	}

	//reversal start neg
	if (molecule.getParameter("rev_begin_x") > reversalLowerBoundNeg &&\
		molecule.getParameter("rev_begin_x") < reversalUpperBoundNeg) {
		molecule.addTag("mobile_begin_neg")
	}

	//reversal end neg
	if (molecule.getParameter("rev_end_x") > reversalLowerBoundNeg &&\
		molecule.getParameter("rev_end_x") < reversalUpperBoundNeg) {
		molecule.addTag("mobile_end_neg")
	}

	//taging nicked molecule
	if(molecule.hasTag("mobile_begin") &&\
	  !molecule.hasTag("coilable20") &&\
	  !molecule.hasTag("coilable2p5") &&\
	  !molecule.hasTag("doubleCoiling") &&\
	   molecule.getParameter("enzymatic") < 0.1 &&\
	   molecule.getParameter("enzymatic") > -0.7){
	   molecule.addTag("nicked")
	}
 
	//tagging track loss and bead loss

	if (deadSlice < Magrot_20f.getEnd() && deadSlice > Magrot_20f.getStart()) {
		molecule.addTag("Magrot20f")
	}
	if (deadSlice < Magrot_2p5f.getEnd() && deadSlice > Magrot_2p5f.getStart()) {
		molecule.addTag("Magrot2p5f")
	}

	if (deadSlice < Region_Torque.getEnd() && deadSlice > Region_Torque.getStart()) {
		molecule.addTag("Torque_Break")
	}

	if (deadSlice < Gyrase_Reaction.getEnd() && deadSlice > Gyrase_Reaction.getStart()) {
    molecule.addTag("Gyrase_Break")
	}

	archive.put(molecule)
})

archive.addLogMessage(LogBuilder.endBlock())

//Drift Calculator
logService.info("Calculating drift...")
DriftCalculatorCommand.calcDrift(archive, "stuckMSD", "x", "y", "x_drift", "y_drift", false, "mean", "end")

//Drift Corrector
logService.info("Correcting for drift...")
DriftCorrectorCommand.correctDrift(archive, driftZeroRegionStart, driftZeroRegionEnd, "x_drift", "y_drift", "x", "y", "x_drift_corr", "y_drift_corr", false)

//Force Calculation and tagging
logService.info("Calculating force...")
archive.getMoleculeUIDs().parallelStream().forEach{ UID ->
      Molecule molecule = archive.get(UID)
      MarsTable table = molecule.getDataTable()

	 double msdForce = table.msd("y_drift_corr", "slice", Force2p5.getStart(), Force2p5.getEnd())   //inserts different start and stop slice with each loop

	 msdForce = msdForce*conversionPixelToMicron*conversionPixelToMicron
	 double[] solution;
	 try {
	 	solution = MarsMath.calculateForceAndLength(persistenceLength, contourLength, temperature, msdForce);
	 } catch (Exception e) {
	 	return;
	 }

	 double force = solution[0]
	 double length = solution[1]
	 String force1 = "Force_PL35";			// make Force_2.5, Force_5 for different flow rates
	 String length1 = "Length_PL35";			// similarly for length

	 molecule.setParameter(force1, force);
	 molecule.setParameter(length1, length);

	 DoubleColumn msdColSlice = new DoubleColumn("MSDs")

	  double MaxMSD = 0
	  double MinMSD = 0

      for (int row = 0; row < table.getRowCount() - slidingForceWindow; row++) {
	      	if (table.get("slice",row) >= Force2p5.getStart() && table.get("slice",row) <= Force2p5.getEnd() && table.getRowCount() > row + slidingForceWindow) {
				double msdSlideForce = table.msd("y_drift_corr","slice",table.get("slice",row),table.get("slice",row + slidingForceWindow))

				if (msdSlideForce > MaxMSD)
					MaxMSD = msdSlideForce

				if (msdSlideForce < MinMSD)
					MinMSD = msdSlideForce

				msdColSlice.add(msdSlideForce)
	      	}
	  }

	  MarsTable tempTable = new MarsTable("MSDs table")
	  tempTable.add(msdColSlice)

	  double msdSTD = tempTable.std("MSDs")
	  double meanMSD = tempTable.mean("MSDs")

	  if (msdSTD*stdSlidingForce > Math.abs(MaxMSD - meanMSD)) {
	  	molecule.addTag("stuckForce")
	  } else if (msdSTD*stdSlidingForce > Math.abs(meanMSD - MinMSD)) {
	  	molecule.addTag("stuckForce")
	  }

      archive.put(molecule)
 }

//Calculate poscycles and negcycles on a molecules-by-molecule basis
logService.info("Calculating molecule-by-molecule cycle numbers...")
archive.getMoleculeUIDs().parallelStream().forEach({ UID ->
      Molecule molecule = archive.get(UID)
      MarsTable table = molecule.getDataTable()

	  double[] pos_coils_per_slice = table.linearRegression("slice", "x_drift_corr", Positive_Coiling_Slope.getStart(), Positive_Coiling_Slope.getEnd())

	  molecule.setParameter("pos_coil_slope", pos_coils_per_slice[2])

	  	if (!table.hasColumn("poscycles"))
			table.appendColumn("poscycles")

		for (int row = 0; row < table.getRowCount(); row++) {
			double conversion = (-1)*(turns_per_slice/pos_coils_per_slice[2])/turns_per_cycle
			table.setValue("poscycles", row, conversion*table.getValue("x_drift_corr", row))
		}

	  double[] neg_coils_per_slice = table.linearRegression("slice", "x_drift_corr", Negative_Coiling_Slope.getStart(), Negative_Coiling_Slope.getEnd())

	  molecule.setParameter("neg_coil_slope", neg_coils_per_slice[2])

	  	if (!table.hasColumn("negcycles"))
			table.appendColumn("negcycles")

		for (int row = 0; row < table.getRowCount(); row++) {
			double conversion = (1)*(turns_per_slice/neg_coils_per_slice[2])/turns_per_cycle
			table.setValue("negcycles", row, conversion*table.getValue("x_drift_corr", row))
		}

      archive.put(molecule)
})

//Calculate global median values for pos_coil_slope and neg_coil_slope

posCoilTable = new MarsTable("PosCoilRates","poscoil")
negCoilTable = new MarsTable("NegCoilRates","negcoil")

logService.info("Calculating global median positive and negative coiling slopes...")
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

//Use global pos_coil_slope and neg_coil_slope to calculate poscycles and negcycles based on global averages.
logService.info("Calculating global cycle numbers...")
archive.getMoleculeUIDs().parallelStream().forEach({ UID ->
	      Molecule molecule = archive.get(UID)
	      MarsTable table = molecule.getDataTable()

		  	if (!table.hasColumn("negcyclesG"))
				table.appendColumn("negcyclesG")

			for (int row = 0; row < table.getRowCount(); row++) {
				double conversion = (1)*(turns_per_slice/negCoilGlobal)/turns_per_cycle
				table.setValue("negcyclesG", row, conversion*table.getValue("x_drift_corr", row))
			}

			if (!table.hasColumn("poscyclesG"))
				table.appendColumn("poscyclesG")

			for (int row = 0; row < table.getRowCount(); row++) {
				double conversion = (-1)*(turns_per_slice/posCoilGlobal)/turns_per_cycle
				table.setValue("poscyclesG", row, conversion*table.getValue("x_drift_corr", row))
			}

	      archive.put(molecule)
	})

logService.info("Calculating activity score...")
archive.getMoleculeUIDs().parallelStream()\
.filter{ UID -> archive.moleculeHasTag(UID, "singleTether")}\
.filter{ UID -> archive.moleculeHasTag(UID, "coilable2p5")}\
.forEach({ UID ->
    Molecule molecule = archive.get(UID)
    MarsTable table = molecule.getDataTable()

    if (table == null)
        return

	 double startTime = Gyrase_Reaction.getStart()
     double endTime = Gyrase_Reaction.getEnd()

     if (table.getValue("slice", table.getRowCount()-1) < Gyrase_Reaction.getEnd()) {
       endTime = table.getValue("slice", table.getRowCount()-1) - 10
     }

    double gMax = -1000000
    double gMin =1000000

    for (int row=0; row<table.getRowCount();row++) {
        double start = table.getValue("slice", row)
        double end = table.getValue("slice", row)

        if (start >= startTime && end <= endTime) {
            double cMax = table.getValue("poscyclesG", row)
            if (gMax < cMax)
                gMax = cMax

            double cMin = table.getValue("poscyclesG", row)
            if (gMin > cMin)
            	gMin = cMin
        }
    }

	molecule.setParameter("activity_score", gMax - gMin)

	archive.put(molecule)
})
String titleBlock2 = LogBuilder.buildTitleBlock("FMT Pipeline 1 - End")
logService.info(titleBlock2)
archive.addLogMessage(titleBlock2)

logger2 = new LogBuilder()
String params = logger2.buildParameterList()
logService.info(params)
archive.addLogMessage(params)

logService.info(LogBuilder.endBlock(true))
archive.addLogMessage(LogBuilder.endBlock(true))

archive.unlock()
