This program is designed as a complementary tool to CREST, enabling a fully automated exploration of the conformational degrees of freedom within reaction profiles at DFT level (with Gaussian) and with a reasonable computational effort.
Starting from a Gaussian output file (opt+freq) and a user-defined parameters list, FASTCAR proceeds accordingly:
- run a CREST calculation, producing a conformer ensemble within a 6 kcal/mol energy window (DFTB);
- reject rotamers by the use of a internal tool from CREST (cregen);
- reject structures by the computation of a similarity index (invariant RMSD, threshold given by user);
- reoptimise all structures at a first DFT level using Gaussian (+freq if transition state);
- reject duplicate, unconverged structures ; in the case of transition states, also reject structures with a first vibration modulus less than 0.5 times the reference one;
- if required, recompute at a higher level of theory, and/or perform Intrinsic Reaction Profile calculations or charge (Hirshfeld, NPA) or NBO analyses.

See our publication in PCCP: https://doi.org/10.1039/D4CP01721H

Python 3 program, relying on the following dependencies: SpyRMSD, CREST (including XTB), Gaussian, Slurm. 
Required Python libraries: glob, math, numpy, os, re, shutil, subprocess, sys, unicodedata.
