This package contains modules and configuration used to reduce the MC truth payload by removing
objects that don't result in observable signals.  Because many MC truth objects contain Ptrs to
other objects, the entire network of objects must be re-written to maintain coherence.  There
are 2 types of compression: StepCompression compresses truth based on aggregate MC energy deposits
in sensitive volumes (DetectorSteps), while DigiCompression works on the output of
simulated digitization, after noise, electronics effects, etc.
