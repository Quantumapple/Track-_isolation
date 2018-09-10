Parameter measure
=================

The test.C code measures deltaZ parameter. <br>
The test-DeltaAngleTotalRegion.C code measures deltaEta(pixel, pixel), deltaEta(PV, pixel), deltaPhi-difference overall region based on Junho's eta range. <br>
<br>

## How to use deltaZ code
"test-DzValid.C" measures deltaZ as a distribution by calculating the distance between gen-level vertex and reconstrcted vertex on z-axis.


## How to run test-DeltaAngleTotalRegion code
**Caution!!** You have to change file name "test-DeltaAngleTotalRegion.h" to "test.h". <br>
To run this file, open "test-DeltaAngleTotalRegion.C". <br>
Turn on and off eta\_region flag depending on your interest region.

```
//if( fabs(propGenEta) < 0.8  ) eta_region = 1;
//if( fabs(propGenEta) > 0.8 && fabs(propGenEta) < 1.4 ) eta_region = 2;
//if( fabs(propGenEta) > 1.4 && fabs(propGenEta) < 1.8 ) eta_region = 3;
//if( fabs(propGenEta) > 1.8 && fabs(propGenEta) < 2.7 ) eta_region = 4;
//if( fabs(propGenEta) > 2.7 && fabs(propGenEta) < 2.9 ) eta_region = 5;
//if( fabs(propGenEta) > 2.9 && fabs(propGenEta) < 3.0 ) eta_region = 6;
``` 
