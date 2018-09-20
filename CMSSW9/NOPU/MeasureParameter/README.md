Parameter measure
=================

The test-DeltaZtotalRegion.C code measures deltaZ parameter in all SW eta range <br>
The test-DeltaAngleTotalRegion.C code measures deltaEta(pixel, pixel), deltaEta(PV, pixel), deltaPhi-difference overall region based on Junho's eta range. <br>
<br>

## How to use deltaZ code
"test-DzValid.C" measures deltaZ as a distribution by calculating the distance between gen-level vertex and reconstrcted vertex on z-axis. <br>
**Caution!!** You have to change file name "test-DeltaZTotalRegion.h" or "test-DeltaZHakseongSW.h"to "test.h". <br>
To run this file, open "test-DeltaZtotalRegion.C" or "test-DeltaZHakseongSW.C". <br>
Turn on and off eta\_region flag depending on your interest region. <br>

**test-DeltaZtotalRegion.C**
```
//if( fabs(propGenEta) < 0.8  ) eta_region = 1;
//if( fabs(propGenEta) > 0.8 && fabs(propGenEta) < 1.4 ) eta_region = 2;
//if( fabs(propGenEta) > 1.4 && fabs(propGenEta) < 1.8 ) eta_region = 3;
//if( fabs(propGenEta) > 1.8 && fabs(propGenEta) < 2.7 ) eta_region = 4;
//if( fabs(propGenEta) > 2.7 && fabs(propGenEta) < 2.9 ) eta_region = 5;
//if( fabs(propGenEta) > 2.9 && fabs(propGenEta) < 3.0 ) eta_region = 6;
``` 

**test-DeltaZHakseongSW.C**
```
//if( fabs(propGenEta) < 0.8  ) eta_region = 1;                              // L1234
//if( fabs(propGenEta) > 0.8 && fabs(propGenEta) < 1.15 ) eta_region = 2;    // L1234 
//if( fabs(propGenEta) > 1.15 && fabs(propGenEta) < 1.4 ) eta_region = 3;    // L123D1
//if( fabs(propGenEta) > 1.4 && fabs(propGenEta) < 1.7 ) eta_region = 4;     // L12D12
//if( fabs(propGenEta) > 1.7 && fabs(propGenEta) < 2.25 ) eta_region = 5;    // L1D123
//if( fabs(propGenEta) > 2.25 && fabs(propGenEta) < 2.7 ) eta_region = 6;    // D2345
//if( fabs(propGenEta) > 2.7 && fabs(propGenEta) < 3.0 ) eta_region = 7;     // D3456
```


## How to run test-DeltaAngleTotalRegion code
**Caution!!** You have to change file name "test-DeltaAngleTotalRegion.h" or "test-DeltaAngleHakseongSW.h" to "test.h". <br>
To run this file, open "test-DeltaAngleTotalRegion.C" or "test-DeltaAngleHakseongSW.C". <br>
Turn on and off eta\_region flag depending on your interest region.

**test-DeltaAngleTotalRegion.C**
```
//if( fabs(propGenEta) < 0.8  ) eta_region = 1;
//if( fabs(propGenEta) > 0.8 && fabs(propGenEta) < 1.4 ) eta_region = 2;
//if( fabs(propGenEta) > 1.4 && fabs(propGenEta) < 1.8 ) eta_region = 3;
//if( fabs(propGenEta) > 1.8 && fabs(propGenEta) < 2.7 ) eta_region = 4;
//if( fabs(propGenEta) > 2.7 && fabs(propGenEta) < 2.9 ) eta_region = 5;
//if( fabs(propGenEta) > 2.9 && fabs(propGenEta) < 3.0 ) eta_region = 6;
```

**test-DeltaAngleHakseongSW.C**
```
//if( fabs(propGenEta) < 0.8  ) eta_region = 1;                              // L1234
//if( fabs(propGenEta) > 0.8 && fabs(propGenEta) < 1.15 ) eta_region = 2;    // L1234 
//if( fabs(propGenEta) > 1.15 && fabs(propGenEta) < 1.4 ) eta_region = 3;    // L123D1
//if( fabs(propGenEta) > 1.4 && fabs(propGenEta) < 1.7 ) eta_region = 4;     // L12D12
//if( fabs(propGenEta) > 1.7 && fabs(propGenEta) < 2.25 ) eta_region = 5;    // L1D123
//if( fabs(propGenEta) > 2.25 && fabs(propGenEta) < 2.7 ) eta_region = 6;    // D2345
//if( fabs(propGenEta) > 2.7 && fabs(propGenEta) < 3.0 ) eta_region = 7;     // D3456
```
