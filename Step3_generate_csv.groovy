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
#@ MoleculeArchive archive
#@OUTPUT MarsTable table
#@ ImageJ ij

import de.mpg.biochem.mars.molecule.*;
import de.mpg.biochem.mars.table.*;
import de.mpg.biochem.mars.molecule.commands.*
import de.mpg.biochem.mars.util.*
import org.scijava.table.*
import java.util.concurrent.ConcurrentHashMap
import java.util.concurrent.ConcurrentMap
import java.io.File

// Defining columns
class RatesObs {
    double posburstRate
    double posburstPosition
    double negburstRate
    double negburstPosition
    double Force
    String UIDCol
    String tag
}

Gyrase_Reaction = archive.getMetadata(0).getRegion("Gyrase Reaction")

//PosCycles calculation
def observations = new ConcurrentHashMap<String, RatesObs>()

archive.lock()

archive.getMoleculeUIDs().parallelStream()\
.filter{ UID -> archive.moleculeHasTag(UID, "singleTether")}\
.filter{ UID -> archive.moleculeHasTag(UID, "coilable2p5")}\
//.filter{ UID -> !archive.moleculeHasTag(UID, "stuckForce")}\
.forEach{ UID ->
    Molecule molecule = archive.get(UID)
    MarsTable table = molecule.getTable()

    if (table == null)
        return

//    if (molecule.getParameter("activity_score") < 12)
//    	return

    RatesObs obs = new RatesObs()
    obs.UIDCol = UID

    if (molecule.hasTag("script_bead_loss_gyrase") &&\
        molecule.hasTag("bead_loss_manual"))
        obs.tag = "reactionBreak"
    //else if (molecule.hasTag("manual_torque"))
        //obs.tag = "torqueBreak"
    else
        obs.tag = "noBreak"


   obs.Force = molecule.getParameter("Force_PL35")

    //set defaults
    obs.posburstRate = Double.NaN
    obs.posburstPosition = Double.NaN
    obs.negburstRate = Double.NaN
    obs.negburstPosition = Double.NaN

    //positive coil slope burst
   if (!Double.isNaN(molecule.getParameter("pos_coil_slope"))) {
       double timeWindowSize = 12.5
       double startTime = archive.getMetadata(0).getTable().rowStream().filter{row -> row.getValue("T") == Gyrase_Reaction.getStart()}.findFirst().get().getValue("Time (s)")
       double endTime = archive.getMetadata(0).getTable().rowStream().filter{row -> row.getValue("T") == Gyrase_Reaction.getEnd()}.findFirst().get().getValue("Time (s)")

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
        double startTime = archive.getMetadata(0).getTable().rowStream().filter{row -> row.getValue("T") == Gyrase_Reaction.getStart()}.findFirst().get().getValue("Time (s)")
        double endTime = archive.getMetadata(0).getTable().rowStream().filter{row -> row.getValue("T") == Gyrase_Reaction.getEnd()}.findFirst().get().getValue("Time (s)")

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
    observations.put(UID, obs)

}

  	ratesTable = new MarsTable()
	ratesTable.add(new DoubleColumn("PosBurstRate"))
	ratesTable.add(new DoubleColumn("PosBurstPosition"))
  	ratesTable.add(new DoubleColumn("NegBurstRate"))
	ratesTable.add(new DoubleColumn("NegBurstPosition"))
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
      ratesTable.setValue("UID", row, obs.UIDCol)
      ratesTable.setValue("Tag", row, obs.tag)
    row++
}

ratesTable.saveAsCSV(archive.getFile().getParent() + "/Gyrase_Scatter.csv")

println("Generated rates table for " + archive.getFile().getName())

archive.unlock()

println("Done!!!")
