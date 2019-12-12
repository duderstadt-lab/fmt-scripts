#@ MoleculeArchive archive

import de.mpg.biochem.mars.molecule.*;
import de.mpg.biochem.mars.table.*;
import de.mpg.biochem.mars.util.*;

MarsImageMetadata metadata = archive.getImageMetadata(0);

//MarsRegion(name, column, start, end, hex color, opacity (0-1))

metadata.putRegion(new MarsRegion("coil20 Positive Peak", "slice", 1135, 1145, "#FFCA28", 0.2))
metadata.putRegion(new MarsRegion("coil20 Negative Peak", "slice", 535, 545, "#42A5F5", 0.2))

metadata.putRegion(new MarsRegion("coil2p5 Peak", "slice", 4010, 4020, "#FFCA28", 0.2))
metadata.putRegion(new MarsRegion("coil2p5 Background", "slice", 2700, 2710, "#42A5F5", 0.2))

metadata.putRegion(new MarsRegion("First Reversal RF", "slice", 100, 120, "#FFCA28", 0.2))
metadata.putRegion(new MarsRegion("First Reversal FF", "slice", 200, 220, "#42A5F5", 0.2))

//metadata.putRegion(new MarsRegion("Last Reversal RF", "slice", 14100, 14120, "#42A5F5", 0.2))
//metadata.putRegion(new MarsRegion("Last Reversal FF", "slice", 13960, 13980, "#FFCA28", 0.2))

metadata.putRegion(new MarsRegion("Last Reversal RF", "slice", 17800, 17820, "#42A5F5", 0.2))
metadata.putRegion(new MarsRegion("Last Reversal FF", "slice", 17900, 17920, "#FFCA28", 0.2))

metadata.putRegion(new MarsRegion("Negative Coiling Slope", "slice", 350, 450, "#66BB6A", 0.2))

metadata.putRegion(new MarsRegion("Before Enzyme", "slice", 4500, 4520, "#42A5F5", 0.2))
metadata.putRegion(new MarsRegion("After Enzyme", "slice", 7800, 7820, "#FFCA28", 0.2))

metadata.putRegion(new MarsRegion("Force2p5", "slice", 1710, 2110, "#F44336", 0.2))

archive.putImageMetadata(metadata)
