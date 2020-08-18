#pragma once

#define JAN12006  53736                 /* MJD of 2006.1.1 */
#define SECPERHOUR 3600.0                /* Seconds per hour */
#define SECPERDAY  86400.0                /* Seconds per day */
#define SECPERWEEK 604800.0               /* Seconds per week */
#define MJD_J2000  51544.5       /* Modif. Julian Date of J2000.0 */

typedef struct COMMONTIME {
	unsigned short Year;
	unsigned short Month;
	unsigned short Day;
	unsigned short Hour;
	unsigned short Minute;
	double Second;
}COMMONTIME;

typedef struct GPSTIME {
	unsigned short Week;
	double         SecOfWeek;
}GPSTIME;

typedef struct MJDTIME {
	unsigned int Days;
	double FracDay;
}MJDTIME;


void CheckGPSTime(GPSTIME* GT);
void CheckMjdTime(MJDTIME* Mjd);
void CommonTimeToMJDTime(const COMMONTIME* CT, MJDTIME* MJDT);
void MJDTimeToCommonTime(const MJDTIME* MJDT, COMMONTIME* CT);
void GPSTimeToMJDTime(const GPSTIME* GT, MJDTIME* MJDT);
void MJDTimeToGPSTime(const MJDTIME* MJDT, GPSTIME* GT);
void CommonTimeToGPSTime(const COMMONTIME* CT, GPSTIME* GT);
void GPSTimeToCommonTime(const GPSTIME* GT, COMMONTIME* CT);
double GetDifGPSTime(const GPSTIME* GT1, const GPSTIME* GT2);
