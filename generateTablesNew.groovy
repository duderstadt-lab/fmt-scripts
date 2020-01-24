#@ MoleculeArchive archive
#@ ImageJ ij

import de.mpg.biochem.mars.molecule.*;
import de.mpg.biochem.mars.table.*;
import de.mpg.biochem.mars.molecule.commands.*
import de.mpg.biochem.mars.util.*
import org.scijava.table.*
import java.util.concurrent.ConcurrentHashMap
import java.util.concurrent.ConcurrentMap
import java.io.File


class RatesObs {
    double posburstRate
    double posburstPosition
    double negburstRate
    double negburstPosition
    double recoveryBurstRate
    double recoveryBurstPosition
    double Force
    String UIDCol
    String tag
}

Gyrase_Reaction = metadata.getRegion("Gyrase Reaction")
Before_Torque_Recovery = metadata.getRegion("Before Torque Recovery")
After_Torque_Recovery = metadata.getRegion("After Torque Recovery")


def generateTables(MoleculeArchive archive, String filePath) {

	//PosCycles calculation
	def observations = new ConcurrentHashMap<String, RatesObs>()

	archive.lock()

	archive.getMoleculeUIDs().parallelStream()\
	.filter{ UID -> archive.moleculeHasTag(UID, "singleTether")}\
	.filter{ UID -> archive.moleculeHasTag(UID, "coilable2p5")}\
	.filter{ UID -> !archive.moleculeHasTag(UID, "stuckForce")}\
	.forEach{ UID ->
	    Molecule molecule = archive.get(UID)
	    MarsTable table = molecule.getDataTable()

	    if (table == null)
	        return

	    if (molecule.getParameter("activity_score") < 12)
	    	return

      RatesObs obs = new RatesObs()
      obs.UIDCol = UID

      if (molecule.hasTag("bead_loss_manual"))
          obs.tag = "reactionBreak"
      //else if (molecule.hasTag("manual_torque"))
          //obs.tag = "torqueBreak"
      else
          obs.tag = "noBreak"
          
	   //if (molecule.getParameter("Force_PL35") == Double.NaN)
	    	//return

	   obs.Force = molecule.getParameter = ("Force_PL35")
	   
      //set defaults
      obs.posburstRate = Double.NaN
      obs.posburstPosition = Double.NaN
      obs.negburstRate = Double.NaN
      obs.negburstPosition = Double.NaN
      obs.recoveryBurstRate = Double.NaN
      obs.recoveryBurstPosition = Double.NaN

      //positive coil slope burst
     if (!Double.isNaN(molecule.getParameter("pos_coil_slope"))) {
         double timeWindowSize = 12.5
         double startTime = archive.getImageMetadata(0).getDataTable().rowStream().filter{row -> row.getValue("slice") == Gyrase_Reaction.getStart()}.findFirst().get().getValue("Time (s)")
         double endTime = archive.getImageMetadata(0).getDataTable().rowStream().filter{row -> row.getValue("slice") == Gyrase_Reaction.getEnd()}.findFirst().get().getValue("Time (s)")

         if (table.getValue("Time (s)", table.getRowCount()-1) < endTime) {
           endTime = table.getValue("Time (s)", table.getRowCount()-1) - 4
         }

    	    //First find max value position and set as new end point
    	    double newEndPoint = endTime
    	    double gMax = 0
    	    for (int row=0; row<table.getRowCount();row++) {
    	        double start = table.getValue("Time (s)", row)
    	        double end = table.getValue("Time (s)", row) + timeWindowSize

    	        if (start >= startTime && end <= endTime) {
    	            double cMax = table.getValue("poscyclesG", row)
    	            if (gMax < cMax) {
    	                gMax = cMax
    	                newEndPoint = table.getValue("Time (s)", row)
    	            }
    	        }
    	    }
          endTime = newEndPoint

          double maxRate = 0
        	double maxPos = 0
        	double[] linearFit = new double[4]

        	for (int row=0; row<table.getRowCount();row++) {
        	    double start = table.getValue("Time (s)", row)
        	    double end = table.getValue("Time (s)", row) + timeWindowSize

        	    if (start >= startTime && end <= endTime) {
        	        linearFit = table.linearRegression("Time (s)", "poscyclesG", start, end)
        	        if (maxRate < linearFit[2]) {
        	            maxRate = linearFit[2]
        	            maxPos = table.getValue("Time (s)", row) + timeWindowSize/2
        	        }
        	    }
        	}

          obs.posburstRate = maxRate
          obs.posburstPosition = maxPos
    	}

      //negative coil slope burst
      if (!molecule.hasTag("chi") && !Double.isNaN(molecule.getParameter("neg_coil_slope"))) {
          double timeWindowSize = 25
          double startTime = archive.getImageMetadata(0).getDataTable().rowStream().filter{row -> row.getValue("slice") == Gyrase_Reaction.getStart()}.findFirst().get().getValue("Time (s)")
          double endTime = archive.getImageMetadata(0).getDataTable().rowStream().filter{row -> row.getValue("slice") == Gyrase_Reaction.getEnd()}.findFirst().get().getValue("Time (s)")

          if (table.getValue("Time (s)", table.getRowCount()-1) < endTime) {
            endTime = table.getValue("Time (s)", table.getRowCount()-1) - 4
          }

    	    //First find min value position and set as new start point
    	    double newStartPoint = startTime
    	    double gMin = 100000
    	    for (int row=0; row<table.getRowCount();row++) {
    	        double start = table.getValue("Time (s)", row)
    	        double end = table.getValue("Time (s)", row) + timeWindowSize

    	        if (start >= startTime && end <= endTime) {
    	            double cMin = table.getValue("negcyclesG", row)
    	            if (gMin > cMin) {
    	                gMin = cMin
    	                newStartPoint = table.getValue("Time (s)", row)
    	            }
    	        }
    	    }
          startTime = newStartPoint

          double maxRate = 0
        	double maxPos = 0
        	double[] linearFit = new double[4]

        	for (int row=0; row<table.getRowCount();row++) {
        	    double start = table.getValue("Time (s)", row)
        	    double end = table.getValue("Time (s)", row) + timeWindowSize

        	    if (start >= startTime && end <= endTime) {
        	        linearFit = table.linearRegression("Time (s)", "negcyclesG", start, end)
        	        if (maxRate < linearFit[2]) {
        	            maxRate = linearFit[2]
        	            maxPos = table.getValue("Time (s)", row) + timeWindowSize/2
        	        }
        	    }
        	}

          obs.negburstRate = maxRate*(-1)
          obs.negburstPosition = maxPos
      }

      if (!Double.isNaN(molecule.getParameter("neg_coil_slope"))) {
          double timeWindowSize = 25
          double startTime = archive.getImageMetadata(0).getDataTable().rowStream().filter{row -> row.getValue("slice") == Before_Torque_Recovery.getEnd()}.findFirst().get().getValue("Time (s)")
          double endTime = archive.getImageMetadata(0).getDataTable().rowStream().filter{row -> row.getValue("slice") == After_Torque_Recovery.getStart()}.findFirst().get().getValue("Time (s)")


          if (table.getValue("Time (s)", table.getRowCount()-1) < endTime) {
            endTime = table.getValue("Time (s)", table.getRowCount()-1) - 4
          }

          double maxRate = 0
        	double maxPos = 0
        	double[] linearFit = new double[4]

        	for (int row=0; row<table.getRowCount();row++) {
        	    double start = table.getValue("Time (s)", row)
        	    double end = table.getValue("Time (s)", row) + timeWindowSize

        	    if (start >= startTime && end <= endTime) {
        	        linearFit = table.linearRegression("Time (s)", "negcyclesG", start, end)
        	        if (maxRate < linearFit[2]) {
        	            maxRate = linearFit[2]
        	            maxPos = table.getValue("Time (s)", row) + timeWindowSize/2
        	        }
        	    }
        	}

          obs.recoveryBurstRate = maxRate*(-1)
          obs.recoveryBurstPosition = maxPos
      }

      observations.put(UID, obs)
	}

    MarsTable ratesTable = new MarsTable()
	ratesTable.add(new DoubleColumn("PosBurstRate"))
	ratesTable.add(new DoubleColumn("PosBurstPosition"))
    ratesTable.add(new DoubleColumn("NegBurstRate"))
	ratesTable.add(new DoubleColumn("NegBurstPosition"))
    //ratesTable.add(new DoubleColumn("TorqueMeanRate"))
    ratesTable.add(new DoubleColumn("RecoveryBurstRate"))
    ratesTable.add(new DoubleColumn("RecoveryBurstPosition"))
    ratesTable.add(new DoubleColumn("Force"))
	ratesTable.add(new GenericColumn("UID"))
    ratesTable.add(new GenericColumn("Tag"))

	int row = 0
	for (RatesObs obs : observations.values()) {
	    ratesTable.appendRow()
	    ratesTable.setValue("PosBurstRate", row, obs.posburstRate)
        ratesTable.setValue("PosBurstPosition", row, obs.posburstPosition)
        ratesTable.setValue("NegBurstRate", row, obs.negburstRate)
        ratesTable.setValue("NegBurstPosition", row, obs.negburstPosition)
        ratesTable.setValue("Force", row, obs.Force)
        ratesTable.setValue("RecoveryBurstRate", row, obs.recoveryBurstRate)
        ratesTable.setValue("RecoveryBurstPosition", row, obs.recoveryBurstPosition)
	    ratesTable.setValue("UID", row, obs.UIDCol)
        ratesTable.setValue("Tag", row, obs.tag)
	    row++
	}

	ratesTable.saveAsCSV(filePath + "RatesTable.csv")

	println("Generated rates table for " + filePath)

	archive.unlock()

	println("Done!!!")
}

//ratesTable.saveAsCSV(filePath + "RatesTable.csv")String dirBase = "/fs/pool/pool-duderstadt/Karl/Projects/FMT/BigBubblePlots/CiproDatasets/"

//generateTables(new SingleMoleculeArchive(new File(dirBase + "190822gyrcipro0p2nM1p25uMposcoil.yama.store")), dirBase + "190822gyrcipro0p2nM1p25uMposcoil_ActivityFilter")
