#@ MoleculeArchive archive

import de.mpg.biochem.mars.molecule.*
import de.mpg.biochem.mars.table.*
import de.mpg.biochem.mars.util.*

MarsMetadata metadata = archive.getMetadata(0)

archive.lock()
//MarsRegion(name, column, start, end, hex color, opacity (0-1))

metadata.putRegion(new MarsRegion("coil20 Positive Peak", "slice", 1135, 1145, "#FFCA28", 0.2))
metadata.putRegion(new MarsRegion("coil20 Negative Peak", "slice", 535, 545, "#42A5F5", 0.2))

metadata.putRegion(new MarsRegion("coil2p5 Peak", "slice", 3900, 3910, "#FFCA28", 0.2))
metadata.putRegion(new MarsRegion("coil2p5 Background", "slice", 2700, 2710, "#42A5F5", 0.2))

metadata.putRegion(new MarsRegion("First Reversal RF", "slice", 100, 120, "#FFCA28", 0.2))
metadata.putRegion(new MarsRegion("First Reversal FF", "slice", 200, 220, "#42A5F5", 0.2))

metadata.putRegion(new MarsRegion("Last Reversal RF", "slice", 31115, 31135, "#42A5F5", 0.2))
metadata.putRegion(new MarsRegion("Last Reversal FF", "slice", 31215, 31235, "#FFCA28", 0.2))

metadata.putRegion(new MarsRegion("Slope_Neg_20", "slice", 350, 450, "#66BB6A", 0.2))

metadata.putRegion(new MarsRegion("Negative Coiling Slope", "slice", 2950, 3050, "#66BB6A", 0.2))
metadata.putRegion(new MarsRegion("Positive Coiling Slope", "slice", 3750, 3850, "#66BB6A", 0.2))

metadata.putRegion(new MarsRegion("Before Enzyme", "slice", 4500, 4520, "#42A5F5", 0.2))
metadata.putRegion(new MarsRegion("After Enzyme", "slice", 7800, 7820, "#FFCA28", 0.2))

metadata.putRegion(new MarsRegion("Before Torque Recovery", "slice", 15000, 15020, "#66BB6A", 0.2))
metadata.putRegion(new MarsRegion("After Torque Recovery", "slice", 15580, 15600, "#66BB6A", 0.2))

metadata.putRegion(new MarsRegion("Force2p5", "slice", 1660, 2060, "#F44336", 0.2))

metadata.putRegion(new MarsRegion("Gyrase Reaction", "slice", 4600, 7957, "#42A5F5", 0.2))

metadata.putRegion(new MarsRegion("Magrot20f", "slice", 238, 1438, "#F44336", 0.2))

metadata.putRegion(new MarsRegion("Magrot2p5f", "slice", 2770, 4130, "#F44336", 0.2))

metadata.putRegion(new MarsRegion("regionTorque", "slice", 7960, 15155, "#42A5F5", 0.2))

archive.putMetadata(metadata)

archive.unlock()
