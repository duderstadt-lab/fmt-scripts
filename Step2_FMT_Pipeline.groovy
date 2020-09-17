/*******************************************************************************
 * Copyright (C) 2020, Duderstadt Lab
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 ******************************************************************************/
@ MoleculeArchive archive
#@ Double (value=0.05) stuckRevThreshold
#@ Double (value=0.01) stuckVarianceThresholdx
#@ Double (value=0.02) stuckVarianceThresholdy
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
#@ Double (value=0.25) turns_per_T
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

import de.mpg.biochem.mars.molecule.*
import de.mpg.biochem.mars.metadata.*
import de.mpg.biochem.mars.table.*
import de.mpg.biochem.mars.util.*
import de.mpg.biochem.mars.molecule.commands.*
import org.scijava.table.*

archive.lock()

//test
//BUILD LOG
String titleBlock = LogBuilder.buildTitleBlock("FMT Pipeline trunc - Start")
logService.info(titleBlock)
archive.logln(titleBlock)

logger = new LogBuilder()
logger.addParameter("FMT Pipeline truncated version", 2)
logger.addParameter("stuckRevThreshold", stuckRevThreshold)
logger.addParameter("stuckVarianceThresholdx", stuckVarianceThresholdx)
logger.addParameter("stuckVarianceThresholdy", stuckVarianceThresholdy)
logger.addParameter("minLength", minLength)
logger.addParameter("doubleCoilingLowerBound", doubleCoilingLowerBound)
logger.addParameter("doubleCoilingUpperBound", doubleCoilingUpperBound)
logger.addParameter("coilable2p5LowerBound", coilable2p5LowerBound)
logger.addParameter("coilable2p5UpperBound", coilable2p5UpperBound)
logger.addParameter("alphaLowerBound", alphaLowerBound)
logger.addParameter("alphaUpperBound", alphaUpperBound)
logger.addParameter("chiLowerBound", chiLowerBound)
logger.addParameter("chiUpperBound", chiUpperBound)
logger.addParameter("coilable20LowerBound", coilable20LowerBound)
logger.addParameter("coilable20UpperBound", coilable20UpperBound)
logger.addParameter("turns_per_T", turns_per_T)
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
logger.addParameter("slidingForceWindow", slidingForceWindow)
logger.addParameter("stdSlidingForce", stdSlidingForce)
//ADD PARAMETERS HERE
//logger.addParameter("name", value)

String parameterList = logger.buildParameterList();
logService.info(parameterList)
archive.addLogMessage(parameterList)

//DONE BUILDING LOG

//Get regions
MarsMetadata metadata = archive.getMetadata(0);

coil20_Positive_Peak = metadata.getRegion("coil20 Positive Peak")
coil20_Negative_Peak = metadata.getRegion("coil20 Negative Peak")

coil2p5_Peak = metadata.getRegion("coil2p5 Peak")
coil2p5_Background = metadata.getRegion("coil2p5 Background")

First_Reversal_RF = metadata.getRegion("First Reversal RF")
First_Reversal_FF = metadata.getRegion("First Reversal FF")

Slope_Neg_20 = metadata.getRegion("Slope_Neg_20")

Negative_Coiling_Slope = metadata.getRegion("Negative Coiling Slope")
Positive_Coiling_Slope = metadata.getRegion("Positive Coiling Slope")

Before_Enzyme = metadata.getRegion("Before Enzyme")
After_Enzyme = metadata.getRegion("After Enzyme")

Force2p5 = metadata.getRegion("Force2p5")

Gyrase_Reaction = metadata.getRegion("Gyrase Reaction")

Magrot_20f = metadata.getRegion("Magrot20f")

Magrot_2p5f = metadata.getRegion("Magrot2p5f")

//ADD Time
logService.info("Adding time...")
AddTimeCommand.addTime(archive)

logService.info("Calculating region differences and adding tags...")
archive.getMoleculeUIDs().parallelStream().forEach({ UID ->
	Molecule molecule = archive.get(UID)
	MarsTable table = molecule.getTable()

	//RegionDifferenceCalculatorCommand.calcRegionDifference(molecule, xColumn, yColumn, regionOne, regionTwo, parameterName)
    //coil20
	RegionDifferenceCalculatorCommand.calcRegionDifference(molecule, "T", "x", coil20_Positive_Peak, coil20_Negative_Peak, "coil20")

	//coil2p5
	RegionDifferenceCalculatorCommand.calcRegionDifference(molecule, "T", "x", coil2p5_Peak, coil2p5_Background, "coil2p5")

	//First Reversal
	RegionDifferenceCalculatorCommand.calcRegionDifference(molecule, "T", "x", First_Reversal_RF, First_Reversal_FF, "rev_begin_x")
	RegionDifferenceCalculatorCommand.calcRegionDifference(molecule, "T", "y", First_Reversal_RF, First_Reversal_FF, "rev_begin_y")

	//Enzymatic Activity Detection
	RegionDifferenceCalculatorCommand.calcRegionDifference(molecule, "T", "x", Before_Enzyme, After_Enzyme, "enzymatic")

	//VarianceCalculatorCommand.calcVariance(molecule, column, parameterName)
	//Variance
	VarianceCalculatorCommand.calcVariance(molecule, "x", "x_Variance")
	VarianceCalculatorCommand.calcVariance(molecule, "y", "y_Variance")

	//Calculate slopes
	output = table.linearRegression("T","x", Slope_Neg_20.getStart(), Slope_Neg_20.getEnd())
	molecule.setParameter("Slope_Neg_20", output[2])

    //Add Tags
    //stuckRev
    if (Math.abs(molecule.getParameter("rev_begin_x")) < stuckRevThreshold &&\
		Math.abs(molecule.getParameter("rev_begin_y")) < stuckRevThreshold) {
			molecule.addTag("stuckRev")
		}

	//stuckVariance
	if (molecule.hasTag("stuckRev") &&\
		molecule.getParameter("x_Variance") < stuckVarianceThresholdx &&\
		molecule.getParameter("y_Variance") < stuckVarianceThresholdy) {
			molecule.addTag("stuckVariance")
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

	//shortTrajectory
	int deadT = table.getValue("T",table.getRowCount()-1)
	if (deadT < minLength)
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

	//reversal start neg
	if (molecule.getParameter("rev_begin_x") > reversalLowerBoundNeg &&\
		molecule.getParameter("rev_begin_x") < reversalUpperBoundNeg) {
		molecule.addTag("mobile_begin_neg")
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

	if (deadT < Magrot_20f.getEnd() && deadT > Magrot_20f.getStart()) {
		molecule.addTag("Magrot20f")
	}
	if (deadT < Magrot_2p5f.getEnd() && deadT > Magrot_2p5f.getStart()) {
		molecule.addTag("Magrot2p5f")
	}

	if (deadT < Gyrase_Reaction.getEnd() && deadT > Gyrase_Reaction.getStart()) {
    molecule.addTag("Gyrase_Break")
	}

	archive.put(molecule)
})

archive.addLogMessage(LogBuilder.endBlock())

//Drift Calculator
logService.info("Calculating drift...")
DriftCalculatorCommand.calcDrift(archive, "stuckVariance", "x", "y", "x_drift", "y_drift", false, "mean", "end")

//Drift Corrector
logService.info("Correcting for drift...")
DriftCorrectorCommand.correctDrift(archive, driftZeroRegionStart, driftZeroRegionEnd, "x_drift", "y_drift", "x", "y", "x_drift_corr", "y_drift_corr", false)

//Force Calculation and tagging
logService.info("Calculating force...")
archive.getMoleculeUIDs().parallelStream().forEach{ UID ->
      Molecule molecule = archive.get(UID)
      MarsTable table = molecule.getTable()

	 double varianceForce = table.variance("y_drift_corr", "T", Force2p5.getStart(), Force2p5.getEnd())   //inserts different start and stop T with each loop

	 varianceForce = varianceForce*conversionPixelToMicron*conversionPixelToMicron
	 double[] solution;
	 try {
	 	solution = MarsMath.calculateForceAndLength(persistenceLength, contourLength, temperature, varianceForce);
	 } catch (Exception e) {
	 	return;
	 }

	 double force = solution[0]
	 double length = solution[1]
	 String force1 = "Force_PL35";			// make Force_2.5, Force_5 for different flow rates
	 String length1 = "Length_PL35";			// similarly for length

	 molecule.setParameter(force1, force);
	 molecule.setParameter(length1, length);

	 DoubleColumn varianceColSlide = new DoubleColumn("Variances")

	  double MaxVariance = 0
	  double MinVariance = 0

      for (int row = 0; row < table.getRowCount() - slidingForceWindow; row++) {
	      	if (table.get("T",row) >= Force2p5.getStart() && table.get("T",row) <= Force2p5.getEnd() && table.getRowCount() > row + slidingForceWindow) {
				double varianceSlideForce = table.variance("y_drift_corr","T",table.get("T",row),table.get("T",row + slidingForceWindow))

				if (varianceSlideForce > MaxVariance)
					MaxVariance = varianceSlideForce

				if (varianceSlideForce < MinVariance)
					MinVariance = varianceSlideForce

				varianceColSlide.add(varianceSlideForce)
	      	}
	  }

	  MarsTable tempTable = new MarsTable("Variances table")
	  tempTable.add(varianceColSlide)

	  double varianceSTD = tempTable.std("Variances")
	  double meanVariance = tempTable.mean("Variances")

	  if (varianceSTD*stdSlidingForce > Math.abs(MaxVariance - meanVariance)) {
	  	molecule.addTag("stuckForce")
	  } else if (varianceSTD*stdSlidingForce > Math.abs(meanVariance - MinVariance)) {
	  	molecule.addTag("stuckForce")
	  }

      archive.put(molecule)
 }

//Calculate poscycles and negcycles on a molecules-by-molecule basis
logService.info("Calculating molecule-by-molecule cycle numbers...")
archive.getMoleculeUIDs().parallelStream().forEach({ UID ->
      Molecule molecule = archive.get(UID)
      MarsTable table = molecule.getTable()

	  double[] pos_coils_per_T = table.linearRegression("T", "x_drift_corr", Positive_Coiling_Slope.getStart(), Positive_Coiling_Slope.getEnd())

	  molecule.setParameter("pos_coil_slope", pos_coils_per_T[2])

	  	if (!table.hasColumn("poscycles"))
			table.appendColumn("poscycles")

		for (int row = 0; row < table.getRowCount(); row++) {
			double conversion = (-1)*(turns_per_T/pos_coils_per_T[2])/turns_per_cycle
			table.setValue("poscycles", row, conversion*table.getValue("x_drift_corr", row))
		}

	  double[] neg_coils_per_T = table.linearRegression("T", "x_drift_corr", Negative_Coiling_Slope.getStart(), Negative_Coiling_Slope.getEnd())

	  molecule.setParameter("neg_coil_slope", neg_coils_per_T[2])

	  	if (!table.hasColumn("negcycles"))
			table.appendColumn("negcycles")

		for (int row = 0; row < table.getRowCount(); row++) {
			double conversion = (1)*(turns_per_T/neg_coils_per_T[2])/turns_per_cycle
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
	      MarsTable table = molecule.getTable()

		  	if (!table.hasColumn("negcyclesG"))
				table.appendColumn("negcyclesG")

			for (int row = 0; row < table.getRowCount(); row++) {
				double conversion = (1)*(turns_per_T/negCoilGlobal)/turns_per_cycle
				table.setValue("negcyclesG", row, conversion*table.getValue("x_drift_corr", row))
			}

			if (!table.hasColumn("poscyclesG"))
				table.appendColumn("poscyclesG")

			for (int row = 0; row < table.getRowCount(); row++) {
				double conversion = (-1)*(turns_per_T/posCoilGlobal)/turns_per_cycle
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
    MarsTable table = molecule.getTable()

    if (table == null)
        return

	 double startTime = Gyrase_Reaction.getStart()
     double endTime = Gyrase_Reaction.getEnd()

     if (table.getValue("T", table.getRowCount()-1) < Gyrase_Reaction.getEnd()) {
       endTime = table.getValue("T", table.getRowCount()-1) - 10
     }

    double gMax = -1000000
    double gMin =1000000

    for (int row=0; row<table.getRowCount();row++) {
        double start = table.getValue("T", row)
        double end = table.getValue("T", row)

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
String titleBlock2 = LogBuilder.buildTitleBlock("FMT Pipeline trunc - End")
logService.info(titleBlock2)
archive.addLogMessage(titleBlock2)

logger2 = new LogBuilder()
String params = logger2.buildParameterList()
logService.info(params)
archive.addLogMessage(params)

logService.info(LogBuilder.endBlock(true))
archive.addLogMessage(LogBuilder.endBlock(true))

archive.unlock()
