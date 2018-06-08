// Include libraries
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

#include <time.h>
clock_t start, end;
double cpu_time_used;

int main(void)
{
  char TmpStr[32];
  double Tmp;

  int Yr,Mth,Day,JDN;
  double JC,GMSTS,GMHA;
  double EstLat,EstLon;
  double RA1,Dec1,Alt1,UT1,LHA1;
  double RA2,Dec2,Alt2,UT2,LHA2;
  double P1,U1,V1;
  double P2,U2,V2;
  double X0,Y0;


  printf("Star Find\n");
  // printf("Year  ");fgets(TmpStr,32,stdin);sscanf(TmpStr,"%d",&Yr);
  // printf("Month ");fgets(TmpStr,32,stdin);sscanf(TmpStr,"%d",&Mth);
  // printf("Day   ");fgets(TmpStr,32,stdin);sscanf(TmpStr,"%d",&Day);
  Yr=1984;
  Mth=5;
  Day=20;

  JDN=367*Yr-(7*(Yr+(Mth+9)/12))/4+275*Mth/9+Day;                         // Julian Day Number
  JC=(JDN-730531.5)/36525.0;                                              // Julian Century
  GMSTS=((((-0.0000062*JC)+0.093104)*JC)+8640184.812866)*JC+24110.54841;  // Greenwich Mean Sidereal Time Seconds
  GMHA=(fmod(GMSTS,86400.0)+86400.0)/240.0;                               // Greenwich Mean Hour Angle

  printf("JDN   %d\n",JDN);
  printf("JC    %lf\n",JC);
  printf("GMSTS %lf\n",GMSTS);
  printf("GMHA  %lf\n",GMHA);

  // printf("Estimated Latitude and Longitude\n");
  // printf("Est Lat ");fgets(TmpStr,32,stdin);sscanf(TmpStr,"%lf",&EstLat);
  // printf("Est Lon ");fgets(TmpStr,32,stdin);sscanf(TmpStr,"%lf",&EstLon);
  EstLat=-30.0;
  EstLon=120.0;

  // printf("Star 1 Observation Data\n");
  // printf("  R.A ");fgets(TmpStr,32,stdin);sscanf(TmpStr,"%lf",&RA1);
  // printf("  Dec ");fgets(TmpStr,32,stdin);sscanf(TmpStr,"%lf",&Dec1);
  // printf("  Alt ");fgets(TmpStr,32,stdin);sscanf(TmpStr,"%lf",&Alt1);
  // printf("  UT  ");fgets(TmpStr,32,stdin);sscanf(TmpStr,"%lf",&UT1);
  RA1=6.0+44.0/60+19.0/3600;
  Dec1=-16.0-37.0/60-10.0/3600;
  Alt1=16.0+40.0/60;
  UT1=12.0+4.0/60+45.0/3600;

  RA1=RA1*15.0;
  if (Alt1>0.0) Alt1=Alt1-58.0/3600.0/tan(Alt1*M_PI/180.0); // Correction for refraction
  LHA1=fmod(GMHA+UT1*15.0*1.00273791,360.0);

  // printf("Star 2 Observation Data\n");
  // printf("R.A ");fgets(TmpStr,32,stdin);sscanf(TmpStr,"%lf",&RA2);
  // printf("Dec ");fgets(TmpStr,32,stdin);sscanf(TmpStr,"%lf",&Dec2);
  // printf("Alt ");fgets(TmpStr,32,stdin);sscanf(TmpStr,"%lf",&Alt2);
  // printf("UT  ");fgets(TmpStr,32,stdin);sscanf(TmpStr,"%lf",&UT2);

  RA2=14.0+14.0/60+43.0/3600;
  Dec2=19.0+17.0/60+24.0/3600;
  Alt2=30.0+30.0/60;
  UT2=12.0+5.0/60+30.0/3600;

  RA2=RA2*15.0;
  if (Alt2>0.0) Alt2=Alt2-58.0/3600.0/tan(Alt2*M_PI/180.0); // Correction for refraction
  LHA2=fmod(GMHA+UT2*15.0*1.00273791,360.0);

  printf("RA1   %lf\n",RA1);
  printf("Dec1  %lf\n",Dec1);
  printf("Alt1  %lf\n",Alt1);
  printf("LHA1  %lf\n",LHA1);
  printf("RA2   %lf\n",RA2);
  printf("Dec2  %lf\n",Dec2);
  printf("ALt2  %lf\n",Alt2);
  printf("LHA2  %lf\n",LHA2);
  do {
    P1=asin(sin(Dec1*M_PI/180)*sin(EstLat*M_PI/180)+cos(Dec1*M_PI/180)*cos(EstLat*M_PI/180)*cos((LHA1-RA1+EstLon)*M_PI/180))*180/M_PI;
    P2=asin(sin(Dec2*M_PI/180)*sin(EstLat*M_PI/180)+cos(Dec2*M_PI/180)*cos(EstLat*M_PI/180)*cos((LHA2-RA2+EstLon)*M_PI/180))*180/M_PI;
    U1=(P1-Alt1)*cos(Dec1*M_PI/180)*sin((LHA1-RA1+EstLon)*M_PI/180)/cos(P1*M_PI/180);
    U2=(P2-Alt2)*cos(Dec2*M_PI/180)*sin((LHA2-RA2+EstLon)*M_PI/180)/cos(P2*M_PI/180);
    V1=(Alt1-P1)*(sin(Dec1*M_PI/180)-sin(EstLat*M_PI/180)*sin(P1*M_PI/180))/cos(EstLat*M_PI/180)/cos(P1*M_PI/180);
    V2=(Alt2-P2)*(sin(Dec2*M_PI/180)-sin(EstLat*M_PI/180)*sin(P2*M_PI/180))/cos(EstLat*M_PI/180)/cos(P2*M_PI/180);
    X0=0.0;
    Y0=0.0;
    if (V1*U2!=U1*V2) {
      X0=((V1*V1+U1*U1)*U2-(V2*V2+U2*U2)*U1)/(V1*U2-U1*V2);
      Y0=((V2*V2+U2*U2)*V1-(V1*V1+U1*U1)*V2)/(V1*U2-U1*V2);
    }
    EstLat=EstLat+X0;
    EstLon=EstLon+Y0;

    printf("P1  %lf\n",P1);
    printf("P2  %lf\n",P2);
    printf("U1  %lf\n",U1);
    printf("U2  %lf\n",U2);
    printf("V1  %lf\n",V1);
    printf("V2  %lf\n",V2);
    printf("X0  %lf\n",X0);
    printf("Y0  %lf\n",Y0);
    printf("Lat %lf\n",EstLat);
    printf("Lon %lf\n\n",EstLon);
  } while ((fabs(X0)>0.001)||(fabs(Y0)>0.001));

  return(0);
}
