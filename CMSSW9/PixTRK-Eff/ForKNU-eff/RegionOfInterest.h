#ifndef RegionOfInterest_h
#define RegionOfInterest_h

double ROI_func(int region, double eget){
  double p[5];

  if(region == 1){
  // fit for median
p[0] = -0.00296036;
p[1] = -1.30531;
p[2] = -0.366932;
p[3] = 0.296419;

  }

  if(region == 2){
p[0] = 0.000775971;
p[1] = -0.785091;
p[2] = -0.989296;
p[3] = -26.4964;


  // fit for median
/*  p[0] = -0.00154114;
  p[1] = -0.401004;
  p[2] = -0.43574;
  p[3] = 0.243757;
  p[4] = 0.828763;
*/
/*p[0] = -0.002281;
p[1] = -1.29856;
p[2] = -0.424359;
p[3] = 0.277163;
*/

  return p[0]*pow(eget,0) + p[1]*pow(eget,p[2])*exp(-pow(eget,p[3]));
  }

  if(region == 3){
  // fit for median
p[0] = -0.00046898;
p[1] = -0.675064;
p[2] = -0.989122;
p[3] = -0.661899;
  }

  if(region == 4){
  // fit for median
p[0] = -0.00269004;
p[1] = -0.262797;
p[2] = 0.0424555;
p[3] = 0.35287;

  }

  if(region == 5){
  // fit for median
p[0] = -0.000525564;
p[1] = -0.340953;
p[2] = -0.56254;
p[3] = 0.206166;
  }



  return p[0]*pow(eget,0) + p[1]*pow(eget,p[2])*exp(-pow(eget,p[3]));
}

#endif 
