@ MoleculeArchive archive

#@ ImageJ ij
#@ LogService logService

import de.mpg.biochem.mars.molecule.*;
import de.mpg.biochem.mars.table.*;
import de.mpg.biochem.mars.util.*;
import de.mpg.biochem.mars.molecule.commands.*

archive.lock()

//Get regions
MarsImageMetadata metadata = archive.getImageMetadata(0);


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

archive.getMoleculeUIDs().parallelStream()\

.filter({ UID -> archive.moleculeHasTag(UID,"singleTether")})\
.filter({ UID -> archive.moleculeHasTag(UID,"coilable2p5")})\

.forEach({ UID ->
      Molecule molecule = archive.get(UID)
      MarsTable table = molecule.getDataTable()

int deadSlice = table.getValue("slice",table.getRowCount()-1)

if (deadSlice < Torque_A.getEnd() && deadSlice > Torque_A.getStart()) {
    molecule.addTag("TorqueA")
}
if (deadSlice < Torque_B.getEnd() && deadSlice > Torque_B.getStart()) {
    molecule.addTag("TorqueB")
}
if (deadSlice < Torque_C.getEnd() && deadSlice > Torque_C.getStart()) {
    molecule.addTag("TorqueC")
}
if (deadSlice < Torque_D.getEnd() && deadSlice > Torque_D.getStart()) {
    molecule.addTag("TorqueD")
}
if (deadSlice < Torque_E.getEnd() && deadSlice > Torque_E.getStart()) {
    molecule.addTag("TorqueE")	
}
if (deadSlice < Torque_F.getEnd() && deadSlice > Torque_F.getStart()) {
    molecule.addTag("TorqueF")
}
if (deadSlice < Torque_G.getEnd() && deadSlice > Torque_G.getStart()) {
    molecule.addTag("TorqueG")
}
if (deadSlice < Torque_H.getEnd() && deadSlice > Torque_H.getStart()) {
    molecule.addTag("TorqueH")
}
if (deadSlice < Torque_I.getEnd() && deadSlice > Torque_I.getStart()) {
    molecule.addTag("TorqueI")
}
if (deadSlice < Torque_J.getEnd() && deadSlice > Torque_J.getStart()) {
    molecule.addTag("TorqueJ")
}
if (deadSlice < Torque_K.getEnd() && deadSlice > Torque_K.getStart()) {
    molecule.addTag("TorqueK")
}

        
	archive.put(molecule)
})

archive.unlock()
