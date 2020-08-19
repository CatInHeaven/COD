#ifndef _REFSYS_H
#define _REFSYS_H


#include "GPSTime.h"
#define TT_TAI   32.184
#define BDT_TAI -33.0
#define Omega_BDS 7.2921150e-5      /* [rad/s], the earth rotation rate */
#define EOPDataFileName "../../data/input/eop_multi.txt"

#define pi 3.1415926535897932384626433832795
#define Rad (pi/180.0)                  /* Radians per degree */
#define Arcs (3600.0*180.0/pi)          /* Arcseconds per radian */

/* 地球自转参数结构体定义 */
typedef struct EOPPARA
{
	double Mjd;          /* 简化儒略日 */
	double x;            /* 极移参数x[角秒] */
	double y;            /* 极移参数y[角秒] */
	double dUT1;         /* UT1-UTC [s]     */
	char   Status;       /* 初始化标志 */

}EOPPARA;


int InitEOPPara(const MJDTIME* CurrTime);
void ICRF_ITRF_GPST(const double Mjd1, const GPSTIME* GT, int flag, double ICRF[6], double ITRF[6]);
void ICRF_ITRF_MJD(const double Mjd1, const MJDTIME* Mjd2, int flag, double ICRF[6], double ITRF[6]);
void InterposeEOP(const MJDTIME* time, EOPPARA* CurrEop);

void PrecMatrix(double Mjd_1, const MJDTIME* Mjd_2, double Mat[]);
void NutMatrix(const MJDTIME* Mjd_TT, double* dpsi, double Mat[]);
double MeanObliquity(const MJDTIME* Mjd_TT);
void NutAngles(const MJDTIME* Mjd_TT, double* dpsi, double* deps);
void GHAMatrix(const MJDTIME* Mjd_UT1, double* dpsi, double Mat[]);
void PoleMatrix(double x, double y, double Mat[]);

void ICRF_ITRF_Matrix(const double Mjd1, const MJDTIME* GT, int flag, double Mat[9]);

double EqnEquinox(const MJDTIME* Mjd_TT, double* dpsi);
double GMST(const MJDTIME* Mjd_UT1);
double GAST(const MJDTIME* Mjd_UT1, double* dpsi);

double TT_UTC();
#endif