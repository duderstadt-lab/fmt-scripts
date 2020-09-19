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

import de.mpg.biochem.mars.molecule.*
import de.mpg.biochem.mars.metadata.*
import de.mpg.biochem.mars.table.*
import de.mpg.biochem.mars.util.*

MarsMetadata metadata = archive.getMetadata(0)
//MarsRegion(name, column, start, end, hex color, opacity (0-1))
metadata.putRegion(new MarsRegion("coil20 Positive Peak", "T", 1135, 1145, "#FFCA28", 0.2))
metadata.putRegion(new MarsRegion("coil20 Negative Peak", "T", 535, 545, "#42A5F5", 0.2))
metadata.putRegion(new MarsRegion("coil2p5 Peak", "T", 3940, 3950, "#FFCA28", 0.2))
metadata.putRegion(new MarsRegion("coil2p5 Background", "T", 2700, 2710, "#42A5F5", 0.2))
metadata.putRegion(new MarsRegion("First Reversal RF", "T", 100, 120, "#FFCA28", 0.2))
metadata.putRegion(new MarsRegion("First Reversal FF", "T", 200, 220, "#42A5F5", 0.2))
metadata.putRegion(new MarsRegion("Slope_Neg_20", "T", 350, 450, "#66BB6A", 0.2))
metadata.putRegion(new MarsRegion("Negative Coiling Slope", "T", 2950, 3050, "#66BB6A", 0.2))
metadata.putRegion(new MarsRegion("Positive Coiling Slope", "T", 3750, 3850, "#66BB6A", 0.2))
metadata.putRegion(new MarsRegion("Before Enzyme", "T", 4500, 4520, "#42A5F5", 0.2))
metadata.putRegion(new MarsRegion("After Enzyme", "T", 6900, 6920, "#FFCA28", 0.2))
metadata.putRegion(new MarsRegion("Force2p5", "T", 1660, 2060, "#F44336", 0.2))
metadata.putRegion(new MarsRegion("Gyrase Reaction", "T", 4600, 6990, "#42A5F5", 0.2))
metadata.putRegion(new MarsRegion("Magrot20f", "T", 238, 1438, "#F44336", 0.2))
metadata.putRegion(new MarsRegion("Magrot2p5f", "T", 2770, 4130, "#F44336", 0.2))

archive.putMetadata(metadata)
