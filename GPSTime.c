#include <math.h>
#include "GPSTime.h"

void CheckGPSTime(GPSTIME* GT) {
	int dw;

	dw = (int)floor(GT->SecOfWeek / SECPERWEEK);
	GT->Week += dw;
	GT->SecOfWeek -= dw * SECPERWEEK;

	if (GT->SecOfWeek < 0.0)
	{
		GT->Week = GT->Week - 1;
		GT->SecOfWeek = GT->SecOfWeek + SECPERWEEK;
	}
	else if (GT->SecOfWeek + 1E-10 > SECPERWEEK)
	{
		GT->Week = GT->Week + 1;
		GT->SecOfWeek = GT->SecOfWeek - SECPERWEEK;
	}
}

void CheckMjdTime(MJDTIME* Mjd)
{
	int d;

	d = (int)floor(Mjd->FracDay);
	if (d == 0)  return;

	Mjd->Days += d;
	Mjd->FracDay -= d;
}

void CommonTimeToMJDTime(const COMMONTIME* CT, MJDTIME* MJDT)
{
	int y, m, temp;

	y = CT->Year;

	if (CT->Year < 100)
	{
		if (CT->Year < 80)
		{
			y = y + 2000;
		}
		else
		{
			y = y + 1900;
		}
	}

	if (CT->Month <= 2)
	{
		y = y - 1;
		m = CT->Month + 12;
	}
	else
	{
		m = CT->Month;
	}

	temp = (int)(365.25 * y);
	temp += (int)(30.6001 * (m + 1));
	temp += CT->Day;
	temp += -679019;

	MJDT->Days = temp;

	MJDT->FracDay = CT->Hour + CT->Minute / 60.0 + CT->Second / SECPERHOUR;
	MJDT->FracDay = MJDT->FracDay / 24.0;

}

void MJDTimeToCommonTime(const MJDTIME* MJDT, COMMONTIME* CT)
{
	int a, b, c, d, e;

	a = (int)(MJDT->Days + MJDT->FracDay + 2400000.5 + 0.5);
	b = a + 1537;
	c = (int)((b - 122.1) / 365.25);
	d = (int)(365.25 * c);
	e = (int)((b - d) / 30.6001);

	CT->Day = b - d - (int)(30.6001 * e);
	CT->Month = e - 1 - 12 * (int)(e / 14);
	CT->Year = c - 4715 - (int)((7 + CT->Month) / 10);

	CT->Hour = (int)(MJDT->FracDay * 24);
	CT->Minute = (int)((MJDT->FracDay * 24 - CT->Hour) * 60);
	CT->Second = ((MJDT->FracDay * 24 - CT->Hour) * 60 - CT->Minute) * 60;
}

void GPSTimeToMJDTime(const GPSTIME* GT, MJDTIME* MJDT)
{
	int day;

	day = (int)(GT->SecOfWeek / SECPERDAY);
	MJDT->FracDay = GT->SecOfWeek / SECPERDAY - day;

	MJDT->Days = JAN12006 + GT->Week * 7 + day;
}

void MJDTimeToGPSTime(const MJDTIME* MJDT, GPSTIME* GT)
{
	int RemainDay;

	GT->Week = (int)((MJDT->Days - JAN12006) / 7);

	RemainDay = MJDT->Days - GT->Week * 7 - JAN12006;

	GT->SecOfWeek = (RemainDay + MJDT->FracDay) * SECPERDAY;
}

void CommonTimeToGPSTime(const COMMONTIME* CT, GPSTIME* GT)
{
	MJDTIME mjd;

	CommonTimeToMJDTime(CT, &mjd);
	MJDTimeToGPSTime(&mjd, GT);
}

void GPSTimeToCommonTime(const GPSTIME* GT, COMMONTIME* CT)
{
	MJDTIME mjd;

	GPSTimeToMJDTime(GT, &mjd);
	MJDTimeToCommonTime(&mjd, CT);
}

double GetDifGPSTime(const GPSTIME* GT1, const GPSTIME* GT2)
{
	double dT;

	dT = (GT1->Week - GT2->Week) * SECPERWEEK;
	dT = dT + GT1->SecOfWeek - GT2->SecOfWeek;

	return dT;
}