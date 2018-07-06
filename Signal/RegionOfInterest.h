#ifndef RegionOfInterest_h
#define RegionOfInterest_h

double ROI_func(double eget, int region, int i, int up_down){
double p[5];

if(region == 1){
     if( i == 0 && up_down == 0){
p[0] = -0.00616578;
p[1] = 0.0130247;
p[2] = 1.06586;
p[3] = 0.291827;
p[4] = 0.418001;
      return p[0]*pow(eget,0) + p[1]*pow(eget,p[2])*exp(-pow(eget,p[3])+p[4]);
     }
     if(  i == 0 && up_down == 1){
p[0] = -0.0476494;
p[1] = -0.51066;
p[2] = -0.175453;
p[3] = 0.305778;
p[4] = 0.798656;
      return p[0]*pow(eget,0) + p[1]*pow(eget,p[2])*exp(-pow(eget,p[3])+p[4]);
     }

     if( i == 1 && up_down == 0){
p[0] = -0.0132098;
p[1] = 0.0143773;
p[2] = 0.779281;
p[3] = 0.238622;
p[4] = 0.940267;
      return p[0]*pow(eget,0) + p[1]*pow(eget,p[2])*exp(-pow(eget,p[3])+p[4]);
     }
     if( i == 1 && up_down == 1){
p[0] = -0.0417477;
p[1] = -0.670655;
p[2] = -0.311782;
p[3] = 0.24879;
p[4] = 0.62779;
      return p[0]*pow(eget,0) + p[1]*pow(eget,p[2])*exp(-pow(eget,p[3])+p[4]);
     }

     if( i == 2 && up_down == 0){
p[0] = -0.0314331;
p[1] = 0.0211268;
p[2] = 0.747704;
p[3] = 0.239927;
p[4] = 0.956927;
      return  p[0]*pow(eget,0) + p[1]*pow(eget,p[2])*exp(-pow(eget,p[3])+p[4]);
     }
     if( i == 2 && up_down == 1){
p[0] = -0.0452316;
p[1] = -0.432685;
p[2] = -0.132724;
p[3] = 0.304626;
p[4] = 0.832704;
      return p[0]*pow(eget,0) + p[1]*pow(eget,p[2])*exp(-pow(eget,p[3])+p[4]);
     }

     if( i == 3 && up_down == 0){
p[0] = -0.0135549;
p[1] = 0.0179451;
p[2] = 0.890612;
p[3] = 0.265159;
p[4] = 0.604681;
      return p[0]*pow(eget,0) + p[1]*pow(eget,p[2])*exp(-pow(eget,p[3])+p[4]);
     }
     if( i == 3 && up_down == 1){
p[0] = -0.0434558;
p[1] = -0.355679;
p[2] = -0.150349;
p[3] = 0.288832;
p[4] = 0.974717;
      return p[0]*pow(eget,0) + p[1]*pow(eget,p[2])*exp(-pow(eget,p[3])+p[4]);
     }
}

if(region == 2){
     if( i == 0 && up_down == 0){
p[0] = 0.0368682;
p[1] = 0.00862309;
p[2] = 0.552664;
p[3] = 0.0998663;
p[4] = -0.117296;
      return p[0]*pow(eget,0) + p[1]*pow(eget,p[2])*exp(-pow(eget,p[3])+p[4]);
     }
     if(  i == 0 && up_down == 1){
p[0] = -0.0411617;
p[1] = -0.431278;
p[2] = -0.667796;
p[3] = 0.100234;
p[4] = 1.41092;
      return p[0]*pow(eget,0) + p[1]*pow(eget,p[2])*exp(-pow(eget,p[3])+p[4]);
     }

     if( i == 1 && up_down == 0){
p[0] = 0.0335665;
p[1] = 0.00734104;
p[2] = 0.584482;
p[3] = 0.119388;
p[4] = -0.183816;
      return p[0]*pow(eget,0) + p[1]*pow(eget,p[2])*exp(-pow(eget,p[3])+p[4]);
     }
     if( i == 1 && up_down == 1){
p[0] = -0.0465157;
p[1] = -0.394909;
p[2] = -0.259955;
p[3] = 0.271155;
p[4] = 1.07077;
      return p[0]*pow(eget,0) + p[1]*pow(eget,p[2])*exp(-pow(eget,p[3])+p[4]);
     }

     if( i == 2 && up_down == 0){
p[0] = 0.0283885;
p[1] = 0.00881962;
p[2] = 0.463245;
p[3] = 0.0948981;
p[4] = -0.115962;
      return  p[0]*pow(eget,0) + p[1]*pow(eget,p[2])*exp(-pow(eget,p[3])+p[4]);
     }
     if( i == 2 && up_down == 1){
p[0] = -0.0442438;
p[1] = -0.400488;
p[2] = -0.452698;
p[3] = 0.221728;
p[4] = 1.29214;
      return p[0]*pow(eget,0) + p[1]*pow(eget,p[2])*exp(-pow(eget,p[3])+p[4]);
     }

     if( i == 3 && up_down == 0){
p[0] = 0.0835003;
p[1] = 0.00531508;
p[2] = 1.03016;
p[3] = 0.281926;
p[4] = 0.301902;
      return p[0]*pow(eget,0) + p[1]*pow(eget,p[2])*exp(-pow(eget,p[3])+p[4]);
     }
     if( i == 3 && up_down == 1){
p[0] = -0.0563229;
p[1] = -0.462875;
p[2] = -0.134261;
p[3] = 0.317652;
p[4] = 0.845493;
      return p[0]*pow(eget,0) + p[1]*pow(eget,p[2])*exp(-pow(eget,p[3])+p[4]);
     }
}

if(region == 3){
     if( i == 0 && up_down == 0){
p[0] = 0.0470136;
p[1] = 0.0103264;
p[2] = 0.494667;
p[3] = 0.0967326;
p[4] = 0.0195524;
      return p[0]*pow(eget,0) + p[1]*pow(eget,p[2])*exp(-pow(eget,p[3])+p[4]);
     }
     if(  i == 0 && up_down == 1){
p[0] = -0.0561412;
p[1] = -0.384896;
p[2] = -0.0704588;
p[3] = 0.327904;
p[4] = 0.763253;
      return p[0]*pow(eget,0) + p[1]*pow(eget,p[2])*exp(-pow(eget,p[3])+p[4]);
     }

     if( i == 1 && up_down == 0){
p[0] = 0.0406149;
p[1] = 0.0080562;
p[2] = 0.547599;
p[3] = 0.130338;
p[4] = -0.151158;
      return p[0]*pow(eget,0) + p[1]*pow(eget,p[2])*exp(-pow(eget,p[3])+p[4]);
     }
     if( i == 1 && up_down == 1){
p[0] = -0.0529226;
p[1] = -0.392632;
p[2] = -0.318888;
p[3] = 0.262534;
p[4] = 1.03539;
      return p[0]*pow(eget,0) + p[1]*pow(eget,p[2])*exp(-pow(eget,p[3])+p[4]);
     }

     if( i == 2 && up_down == 0){
p[0] = 0.044638;
p[1] = 0.00767674;
p[2] = 0.535735;
p[3] = 0.158043;
p[4] = -0.174106;
      return  p[0]*pow(eget,0) + p[1]*pow(eget,p[2])*exp(-pow(eget,p[3])+p[4]);
     }
     if( i == 2 && up_down == 1){
p[0] = -0.0536385;
p[1] = -0.345268;
p[2] = -0.0993922;
p[3] = 0.305879;
p[4] = 0.746996;
      return p[0]*pow(eget,0) + p[1]*pow(eget,p[2])*exp(-pow(eget,p[3])+p[4]);
     }

     if( i == 3 && up_down == 0){
p[0] = 0.0840861;
p[1] = 0.00759447;
p[2] = 1.52426;
p[3] = 0.354157;
p[4] = -0.795288;
      return p[0]*pow(eget,0) + p[1]*pow(eget,p[2])*exp(-pow(eget,p[3])+p[4]);
     }
     if( i == 3 && up_down == 1){
p[0] = -0.0602759;
p[1] = -0.400662;
p[2] = -0.222637;
p[3] = 0.319729;
p[4] = 1.11899;
      return p[0]*pow(eget,0) + p[1]*pow(eget,p[2])*exp(-pow(eget,p[3])+p[4]);
     }
}

if(region == 4){
     if( i == 0 && up_down == 0){
p[0] = 0.0415527;
p[1] = 0.0150208;
p[2] = 0.951691;
p[3] = 0.276098;
p[4] = 0.125928;
      return p[0]*pow(eget,0) + p[1]*pow(eget,p[2])*exp(-pow(eget,p[3])+p[4]);
     }
     if(  i == 0 && up_down == 1){
p[0] = -0.056355;
p[1] = -0.398047;
p[2] = -0.880233;
p[3] = -0.1076;
p[4] = 1.03668;
      return p[0]*pow(eget,0) + p[1]*pow(eget,p[2])*exp(-pow(eget,p[3])+p[4]);
     }

     if( i == 1 && up_down == 0){
p[0] = 0.0513878;
p[1] = 0.00945545;
p[2] = 0.490023;
p[3] = 0.151189;
p[4] = -0.083779;
      return p[0]*pow(eget,0) + p[1]*pow(eget,p[2])*exp(-pow(eget,p[3])+p[4]);
     }
     if( i == 1 && up_down == 1){
p[0] = -0.057823;
p[1] = -0.363936;
p[2] = -0.854073;
p[3] = -0.112329;
p[4] = 0.905337;
      return p[0]*pow(eget,0) + p[1]*pow(eget,p[2])*exp(-pow(eget,p[3])+p[4]);
     }

     if( i == 2 && up_down == 0){
p[0] = 0.0574148;
p[1] = 0.0105657;
p[2] = 0.458329;
p[3] = 0.170437;
p[4] = 0.00875068;
      return  p[0]*pow(eget,0) + p[1]*pow(eget,p[2])*exp(-pow(eget,p[3])+p[4]);
     }
     if( i == 2 && up_down == 1){
p[0] = -0.0587855;
p[1] = -0.391413;
p[2] = -0.932153;
p[3] = -0.139259;
p[4] = 1.02293;
      return p[0]*pow(eget,0) + p[1]*pow(eget,p[2])*exp(-pow(eget,p[3])+p[4]);
     }

     if( i == 3 && up_down == 0){
p[0] = 0.0898407;
p[1] = 0.00472105;
p[2] = 0.725232;
p[3] = 0.212939;
p[4] = -0.256836;
      return p[0]*pow(eget,0) + p[1]*pow(eget,p[2])*exp(-pow(eget,p[3])+p[4]);
     }
     if( i == 3 && up_down == 1){
p[0] = -0.0735164;
p[1] = -0.836734;
p[2] = 0.170064;
p[3] = 0.373057;
p[4] = -0.682527;
      return p[0]*pow(eget,0) + p[1]*pow(eget,p[2])*exp(-pow(eget,p[3])+p[4]);
     }
}

return 0.;
}

#endif 
