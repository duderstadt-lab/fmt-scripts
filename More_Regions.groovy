#@ MoleculeArchive archive

import de.mpg.biochem.mars.molecule.*;
import de.mpg.biochem.mars.table.*;
import de.mpg.biochem.mars.util.*;

MarsImageMetadata metadata = archive.getImageMetadata(0);

archive.lock()
//MarsRegion(name, column, start, end, hex color, opacity (0-1))

metadata.putRegion(new MarsRegion("Before Torque Recovery", "slice", 19500, 19520, "#66BB6A", 0.2))
metadata.putRegion(new MarsRegion("After Torque Recovery", "slice", 20060, 20080, "#66BB6A", 0.2))

metadata.putRegion(new MarsRegion("Magrot20f", "slice", 238, 1438, "#F44336", 0.2))

metadata.putRegion(new MarsRegion("Magrot2p5f", "slice", 2770, 4130, "#F44336", 0.2))

metadata.putRegion(new MarsRegion("Gyrase Reaction", "slice", 5000, 7957, "#42A5F5", 0.2))

metadata.putRegion(new MarsRegion("TorqueA", "slice", 7958, 9000, "#42A5F5", 0.2))
metadata.putRegion(new MarsRegion("TorqueB", "slice", 9001, 10000, "#42A5F5", 0.2))
metadata.putRegion(new MarsRegion("TorqueC", "slice", 10001, 11000, "#42A5F5", 0.2))
metadata.putRegion(new MarsRegion("TorqueD", "slice", 11001, 12000, "#42A5F5", 0.2))
metadata.putRegion(new MarsRegion("TorqueE", "slice", 12001, 13000, "#42A5F5", 0.2))
metadata.putRegion(new MarsRegion("TorqueF", "slice", 13001, 14000, "#42A5F5", 0.2))
metadata.putRegion(new MarsRegion("TorqueG", "slice", 14001, 15000, "#42A5F5", 0.2))
metadata.putRegion(new MarsRegion("TorqueH", "slice", 15001, 16000, "#42A5F5", 0.2))
metadata.putRegion(new MarsRegion("TorqueI", "slice", 16001, 17000, "#42A5F5", 0.2))
metadata.putRegion(new MarsRegion("TorqueJ", "slice", 17001, 18000, "#42A5F5", 0.2))
metadata.putRegion(new MarsRegion("TorqueK", "slice", 18001, 19718, "#42A5F5", 0.2))


archive.putImageMetadata(metadata)

archive.unlock()
