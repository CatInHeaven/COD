#include "RefSys.h"

#include "matrix.h"
#include <math.h>
#include <stdio.h>

static double BDT_UTC = 3;        /* 2006年1月1日开始, 跳秒数[s],在上星前要设置为正确的跳秒 */
EOPPARA EOP[3];
int InitEOPPara(const MJDTIME* CurrTime)
{
	char    line[256] = { 0 };
	//int     index=0;
	double  CurrMjd = 0, TAI_UT1 = 0, TAI_UTC = 0;
	FILE* in = NULL;
	COMMONTIME ct;
	MJDTIME Mjd;
	EOPPARA TempEOP;

	Mjd = *CurrTime;
	CheckMjdTime(&Mjd);
	MJDTimeToCommonTime(&Mjd, &ct);


	if ((in = fopen(EOPDataFileName, "rt")) == NULL)
	{
		printf("The file %s was not opened\n", EOPDataFileName);
		return 0;
	}

	//index = 0;
	CurrMjd = CurrTime->Days + CurrTime->FracDay;
	while (!feof(in))
	{
		if (fgets(line, 256, in) != NULL)
		{
			if (sscanf(line, "%*s %hd %hd %hd %hd %hd %lf %lf %lf %lf", &ct.Year, &ct.Month, &ct.Day, &ct.Hour, &ct.Minute,
				&TempEOP.x, &TempEOP.y, &TAI_UT1, &TAI_UTC) != 9)  continue;
			ct.Second = 0.0;
			CommonTimeToMJDTime(&ct, &Mjd);
			TempEOP.Mjd = Mjd.Days + Mjd.FracDay;
			TempEOP.dUT1 = TAI_UTC - TAI_UT1;
			BDT_UTC = BDT_TAI + TAI_UTC;
			TempEOP.Status = 1;

			if (((CurrMjd - TempEOP.Mjd) > 1.0E-5) && ((CurrMjd - TempEOP.Mjd) < 1.1))
			{
				EOP[0] = TempEOP;   /* Set the first EOP struct */

				if (fgets(line, 256, in) != NULL)    /* Set the second EOP struct */
				{
					if (sscanf(line, "%*s %hd %hd %hd %hd %hd %lf %lf %lf %lf", &ct.Year, &ct.Month, &ct.Day, &ct.Hour, &ct.Minute,
						&TempEOP.x, &TempEOP.y, &TAI_UT1, &TAI_UTC) != 9)  continue;

					CommonTimeToMJDTime(&ct, &Mjd);
					TempEOP.Mjd = Mjd.Days + Mjd.FracDay;
					TempEOP.dUT1 = TAI_UTC - TAI_UT1;
					BDT_UTC = BDT_TAI + TAI_UTC;
					EOP[1] = TempEOP;
					EOP[1].Status = 1;
				}
				else
				{
					printf("No new EOP data.\n");
					fclose(in);
					return 0;
				}

				if (fgets(line, 256, in) != NULL)   /* Set the third EOP struct */
				{
					if (sscanf(line, "%*s %hd %hd %hd %hd %hd %lf %lf %lf %lf", &ct.Year, &ct.Month, &ct.Day, &ct.Hour, &ct.Minute,
						&TempEOP.x, &TempEOP.y, &TAI_UT1, &TAI_UTC) != 9)  continue;

					CommonTimeToMJDTime(&ct, &Mjd);
					TempEOP.Mjd = Mjd.Days + Mjd.FracDay;
					TempEOP.dUT1 = TAI_UTC - TAI_UT1;
					BDT_UTC = BDT_TAI + TAI_UTC;
					EOP[2] = TempEOP;
					EOP[2].Status = 1;
					fclose(in);
					return 1;
				}
				else
				{
					printf("No new EOP data.\n");
					fclose(in);
					return 0;
				}
			}
		}
		else
		{
			EOP[0].Status = 0;
			printf("No new EOP data.\n");
			fclose(in);
			return 0;
		}
	}

	fclose(in);
	return 0;
}

void ICRF_ITRF_GPST(const double Mjd1, const GPSTIME* GT,
	int flag, double ICRF[6], double ITRF[6])
{
	MJDTIME Mjd_TT;

	GPSTimeToMJDTime(GT, &Mjd_TT);    /* 此处Mjd_TT为BDT时间的MJD表示 */
	Mjd_TT.FracDay = Mjd_TT.FracDay - (BDT_TAI - TT_TAI) / SECPERDAY;

	ICRF_ITRF_MJD(Mjd1, &Mjd_TT, flag, ICRF, ITRF);
}

void ICRF_ITRF_MJD(const double Mjd1, const MJDTIME* Mjd2,
	int flag, double ICRF[6], double ITRF[6])

{
	int i;
	MJDTIME Mjd_UT1, Mjd_UTC;
	double dpsi;
	double  U[9], U_dot[9], UT[9], U_dotT[9];
	double Prec[9], Nut[9], GH[9], Pole[9];
	double GH_dot[9];
	double tmp1[9], tmp2[9];
	double v1[3], v2[3];

	static EOPPARA CurrEop;

	Mjd_UTC.Days = Mjd2->Days;
	Mjd_UTC.FracDay = Mjd2->FracDay - TT_UTC() / SECPERDAY;
	InterposeEOP(&Mjd_UTC, &CurrEop);

	Mjd_UT1.Days = Mjd_UTC.Days;
	Mjd_UT1.FracDay = Mjd_UTC.FracDay + CurrEop.dUT1 / SECPERDAY;

	PrecMatrix(Mjd1, Mjd2, Prec);
	NutMatrix(Mjd2, &dpsi, Nut);
	GHAMatrix(&Mjd_UT1, &dpsi, GH);
	PoleMatrix(CurrEop.x, CurrEop.y, Pole);

	for (i = 0; i < 3; i++)
	{
		GH_dot[i] = Omega_BDS * GH[3 + i];
		GH_dot[3 + i] = -Omega_BDS * GH[i];
		GH_dot[6 + i] = 0.0;
	}

	MatrixMultiply(3, 3, 3, 3, Nut, Prec, tmp1);
	MatrixMultiply(3, 3, 3, 3, Pole, GH, tmp2);
	MatrixMultiply(3, 3, 3, 3, tmp2, tmp1, U);   /*坐标转换矩阵*/

	MatrixMultiply(3, 3, 3, 3, Pole, GH_dot, tmp2);
	MatrixMultiply(3, 3, 3, 3, tmp2, tmp1, U_dot);  /*速度转换矩阵 */

	if (flag == 1)     /*天球系——>地固系 */
	{
		MatrixMultiply(3, 3, 3, 1, U, ICRF, ITRF);

		MatrixMultiply(3, 3, 3, 1, U, &ICRF[3], v1);
		MatrixMultiply(3, 3, 3, 1, U_dot, ICRF, v2);

		MatrixAddition(3, 1, v1, v2, &ITRF[3]);

	}
	else  /*地固系——>天球系 */
	{
		MatrixTranspose(3, 3, U, UT);
		MatrixTranspose(3, 3, U_dot, U_dotT);

		MatrixMultiply(3, 3, 3, 1, UT, ITRF, ICRF);

		MatrixMultiply(3, 3, 3, 1, UT, &ITRF[3], v1);
		MatrixMultiply(3, 3, 3, 1, U_dotT, ITRF, v2);

		MatrixAddition(3, 1, v1, v2, &ICRF[3]);
	}
}

void InterposeEOP(const MJDTIME* time, EOPPARA* CurrEop)
{
	int i = 0, j = 0;
	double t = 0;
	double temp = 0;

	t = time->Days + time->FracDay;

	if (fabs(t - EOP[1].Mjd) > 0.8)  /* time不在EOP[2]的时间段内, 在线运行时要考虑EOP文件实际推送时间，避免无法打开文件 */
	{
		if (InitEOPPara(time) == 0 && EOP[0].Mjd < 1.0)  // 无法初始化
		{
			printf("No new EOP data in MJD %d.\n", time->Days); // 实时程序不能退出运行
			return;
		}
	}

	CurrEop->x = 0.0;
	CurrEop->y = 0.0;
	CurrEop->dUT1 = 0.0;

	if ((EOP[0].Status == 1) && (EOP[1].Status == 1) && (EOP[2].Status == 1))
	{
		for (i = 0; i < 3; i++)
		{
			temp = 1.0;

			for (j = 0; j < 3; j++)
			{
				if (i == j)
				{
					continue;
				}

				temp = temp * (t - EOP[j].Mjd) / (EOP[i].Mjd - EOP[j].Mjd);
			}

			CurrEop->x = CurrEop->x + temp * EOP[i].x;
			CurrEop->y = CurrEop->y + temp * EOP[i].y;
			CurrEop->dUT1 = CurrEop->dUT1 + temp * EOP[i].dUT1;
		}
		CurrEop->Mjd = t;
		CurrEop->Status = 1;
	}
	else
	{
		printf("No new EOP data.\n");
		exit(1);
	}

}

void PrecMatrix(double Mjd_1, const MJDTIME* Mjd_2, double Mat[])
{
	double T = 0, dT = 0;
	double zeta = 0, z = 0, theta = 0;
	double M1[9] = { 0 }, M2[9] = { 0 }, M3[9] = { 0 };

	T = (Mjd_1 - MJD_J2000) / 36525.0;
	dT = (Mjd_2->Days - Mjd_1 + Mjd_2->FracDay) / 36525.0;

	/* Precession angles */

	zeta = ((2306.2181 + (1.39656 - 0.000139 * T) * T) +
		((0.30188 - 0.000344 * T) + 0.017998 * dT) * dT) * dT / Arcs;
	z = zeta + ((0.79280 + 0.000411 * T) + 0.000205 * dT) * dT * dT / Arcs;
	theta = ((2004.3109 - (0.85330 + 0.000217 * T) * T) -
		((0.42665 + 0.000217 * T) + 0.041833 * dT) * dT) * dT / Arcs;


	Rotation_z(-z, M1);
	Rotation_y(theta, M2);
	MatrixMultiply(3, 3, 3, 3, M1, M2, M3);

	Rotation_z(-zeta, M1);
	MatrixMultiply(3, 3, 3, 3, M3, M1, Mat);

}

void NutMatrix(const MJDTIME* Mjd_TT, double* dpsi, double Mat[])
{
	double deps = 0, eps = 0;
	double M1[9] = { 0 }, M2[9] = { 0 }, M3[9] = { 0 };

	eps = MeanObliquity(Mjd_TT);  /* Mean obliquity of the ecliptic */

	NutAngles(Mjd_TT, dpsi, &deps); /*Nutation in longitude and obliquity*/

	Rotation_x(-eps - deps, M1);
	Rotation_z(-1.0 * *dpsi, M2);
	MatrixMultiply(3, 3, 3, 3, M1, M2, M3);

	Rotation_x(eps, M1);
	MatrixMultiply(3, 3, 3, 3, M3, M1, Mat);
}

double MeanObliquity(const MJDTIME* Mjd_TT)
{
	double T = 0, Val = 0;

	T = (Mjd_TT->Days + Mjd_TT->FracDay - MJD_J2000) / 36525.0;

	Val = Rad * (23.43929111 -
		(46.8150 + (0.00059 - 0.001813 * T) * T) * T / 3600.0);

	return (Val);

}

void NutAngles(const MJDTIME* Mjd_TT, double* dpsi, double* deps)
{
	int i = 0;
	double T = 0, T2 = 0, T3 = 0, T4 = 0, rev = 0;

	double  l = 0, lp = 0, F = 0, D = 0, Om = 0;
	double  arg = 0;

	static int C[106][9] = {

		/* l  l' F  D Om    dpsi    *T     deps     *T     */

		 {  0, 0, 0, 0, 1,-1719960,-1742,  920250,   89 },   //   1
		 {  0, 0, 0, 0, 2,   20620,    2,   -8950,    5 },   //   2
		 { -2, 0, 2, 0, 1,     460,    0,    -240,    0 },   //   3
		 {  2, 0,-2, 0, 0,     110,    0,       0,    0 },   //   4
		 { -2, 0, 2, 0, 2,     -30,    0,      10,    0 },   //   5
		 {  1,-1, 0,-1, 0,     -30,    0,       0,    0 },   //   6
		 {  0,-2, 2,-2, 1,     -20,    0,      10,    0 },   //   7
		 {  2, 0,-2, 0, 1,      10,    0,       0,    0 },   //   8
		 {  0, 0, 2,-2, 2, -131870,  -16,   57360,  -31 },   //   9
		 {  0, 1, 0, 0, 0,   14260,  -34,     540,   -1 },   //  10
		 {  0, 1, 2,-2, 2,   -5170,   12,    2240,   -6 },   //  11
		 {  0,-1, 2,-2, 2,    2170,   -5,    -950,    3 },   //  12
		 {  0, 0, 2,-2, 1,    1290,    1,    -700,    0 },   //  13
		 {  2, 0, 0,-2, 0,     480,    0,      10,    0 },   //  14
		 {  0, 0, 2,-2, 0,    -220,    0,       0,    0 },   //  15
		 {  0, 2, 0, 0, 0,     170,   -1,       0,    0 },   //  16
		 {  0, 1, 0, 0, 1,    -150,    0,      90,    0 },   //  17
		 {  0, 2, 2,-2, 2,    -160,    1,      70,    0 },   //  18
		 {  0,-1, 0, 0, 1,    -120,    0,      60,    0 },   //  19
		 { -2, 0, 0, 2, 1,     -60,    0,      30,    0 },   //  20
		 {  0,-1, 2,-2, 1,     -50,    0,      30,    0 },   //  21
		 {  2, 0, 0,-2, 1,      40,    0,     -20,    0 },   //  22
		 {  0, 1, 2,-2, 1,      40,    0,     -20,    0 },   //  23
		 {  1, 0, 0,-1, 0,     -40,    0,       0,    0 },   //  24
		 {  2, 1, 0,-2, 0,      10,    0,       0,    0 },   //  25
		 {  0, 0,-2, 2, 1,      10,    0,       0,    0 },   //  26
		 {  0, 1,-2, 2, 0,     -10,    0,       0,    0 },   //  27
		 {  0, 1, 0, 0, 2,      10,    0,       0,    0 },   //  28
		 { -1, 0, 0, 1, 1,      10,    0,       0,    0 },   //  29
		 {  0, 1, 2,-2, 0,     -10,    0,       0,    0 },   //  30
		 {  0, 0, 2, 0, 2,  -22740,   -2,    9770,   -5 },   //  31
		 {  1, 0, 0, 0, 0,    7120,    1,     -70,    0 },   //  32
		 {  0, 0, 2, 0, 1,   -3860,   -4,    2000,    0 },   //  33
		 {  1, 0, 2, 0, 2,   -3010,    0,    1290,   -1 },   //  34
		 {  1, 0, 0,-2, 0,   -1580,    0,     -10,    0 },   //  35
		 { -1, 0, 2, 0, 2,    1230,    0,    -530,    0 },   //  36
		 {  0, 0, 0, 2, 0,     630,    0,     -20,    0 },   //  37
		 {  1, 0, 0, 0, 1,     630,    1,    -330,    0 },   //  38
		 { -1, 0, 0, 0, 1,    -580,   -1,     320,    0 },   //  39
		 { -1, 0, 2, 2, 2,    -590,    0,     260,    0 },   //  40
		 {  1, 0, 2, 0, 1,    -510,    0,     270,    0 },   //  41
		 {  0, 0, 2, 2, 2,    -380,    0,     160,    0 },   //  42
		 {  2, 0, 0, 0, 0,     290,    0,     -10,    0 },   //  43
		 {  1, 0, 2,-2, 2,     290,    0,    -120,    0 },   //  44
		 {  2, 0, 2, 0, 2,    -310,    0,     130,    0 },   //  45
		 {  0, 0, 2, 0, 0,     260,    0,     -10,    0 },   //  46
		 { -1, 0, 2, 0, 1,     210,    0,    -100,    0 },   //  47
		 { -1, 0, 0, 2, 1,     160,    0,     -80,    0 },   //  48
		 {  1, 0, 0,-2, 1,    -130,    0,      70,    0 },   //  49
		 { -1, 0, 2, 2, 1,    -100,    0,      50,    0 },   //  50
		 {  1, 1, 0,-2, 0,     -70,    0,       0,    0 },   //  51
		 {  0, 1, 2, 0, 2,      70,    0,     -30,    0 },   //  52
		 {  0,-1, 2, 0, 2,     -70,    0,      30,    0 },   //  53
		 {  1, 0, 2, 2, 2,     -80,    0,      30,    0 },   //  54
		 {  1, 0, 0, 2, 0,      60,    0,       0,    0 },   //  55
		 {  2, 0, 2,-2, 2,      60,    0,     -30,    0 },   //  56
		 {  0, 0, 0, 2, 1,     -60,    0,      30,    0 },   //  57
		 {  0, 0, 2, 2, 1,     -70,    0,      30,    0 },   //  58
		 {  1, 0, 2,-2, 1,      60,    0,     -30,    0 },   //  59
		 {  0, 0, 0,-2, 1,     -50,    0,      30,    0 },   //  60
		 {  1,-1, 0, 0, 0,      50,    0,       0,    0 },   //  61
		 {  2, 0, 2, 0, 1,     -50,    0,      30,    0 },   //  62
		 {  0, 1, 0,-2, 0,     -40,    0,       0,    0 },   //  63
		 {  1, 0,-2, 0, 0,      40,    0,       0,    0 },   //  64
		 {  0, 0, 0, 1, 0,     -40,    0,       0,    0 },   //  65
		 {  1, 1, 0, 0, 0,     -30,    0,       0,    0 },   //  66
		 {  1, 0, 2, 0, 0,      30,    0,       0,    0 },   //  67
		 {  1,-1, 2, 0, 2,     -30,    0,      10,    0 },   //  68
		 { -1,-1, 2, 2, 2,     -30,    0,      10,    0 },   //  69
		 { -2, 0, 0, 0, 1,     -20,    0,      10,    0 },   //  70
		 {  3, 0, 2, 0, 2,     -30,    0,      10,    0 },   //  71
		 {  0,-1, 2, 2, 2,     -30,    0,      10,    0 },   //  72
		 {  1, 1, 2, 0, 2,      20,    0,     -10,    0 },   //  73
		 { -1, 0, 2,-2, 1,     -20,    0,      10,    0 },   //  74
		 {  2, 0, 0, 0, 1,      20,    0,     -10,    0 },   //  75
		 {  1, 0, 0, 0, 2,     -20,    0,      10,    0 },   //  76
		 {  3, 0, 0, 0, 0,      20,    0,       0,    0 },   //  77
		 {  0, 0, 2, 1, 2,      20,    0,     -10,    0 },   //  78
		 { -1, 0, 0, 0, 2,      10,    0,     -10,    0 },   //  79
		 {  1, 0, 0,-4, 0,     -10,    0,       0,    0 },   //  80
		 { -2, 0, 2, 2, 2,      10,    0,     -10,    0 },   //  81
		 { -1, 0, 2, 4, 2,     -20,    0,      10,    0 },   //  82
		 {  2, 0, 0,-4, 0,     -10,    0,       0,    0 },   //  83
		 {  1, 1, 2,-2, 2,      10,    0,     -10,    0 },   //  84
		 {  1, 0, 2, 2, 1,     -10,    0,      10,    0 },   //  85
		 { -2, 0, 2, 4, 2,     -10,    0,      10,    0 },   //  86
		 { -1, 0, 4, 0, 2,      10,    0,       0,    0 },   //  87
		 {  1,-1, 0,-2, 0,      10,    0,       0,    0 },   //  88
		 {  2, 0, 2,-2, 1,      10,    0,     -10,    0 },   //  89
		 {  2, 0, 2, 2, 2,     -10,    0,       0,    0 },   //  90
		 {  1, 0, 0, 2, 1,     -10,    0,       0,    0 },   //  91
		 {  0, 0, 4,-2, 2,      10,    0,       0,    0 },   //  92
		 {  3, 0, 2,-2, 2,      10,    0,       0,    0 },   //  93
		 {  1, 0, 2,-2, 0,     -10,    0,       0,    0 },   //  94
		 {  0, 1, 2, 0, 1,      10,    0,       0,    0 },   //  95
		 { -1,-1, 0, 2, 1,      10,    0,       0,    0 },   //  96
		 {  0, 0,-2, 0, 1,     -10,    0,       0,    0 },   //  97
		 {  0, 0, 2,-1, 2,     -10,    0,       0,    0 },   //  98
		 {  0, 1, 0, 2, 0,     -10,    0,       0,    0 },   //  99
		 {  1, 0,-2,-2, 0,     -10,    0,       0,    0 },   // 100
		 {  0,-1, 2, 0, 1,     -10,    0,       0,    0 },   // 101
		 {  1, 1, 0,-2, 1,     -10,    0,       0,    0 },   // 102
		 {  1, 0,-2, 2, 0,     -10,    0,       0,    0 },   // 103
		 {  2, 0, 0, 2, 0,      10,    0,       0,    0 },   // 104
		 {  0, 0, 2, 4, 2,     -10,    0,       0,    0 },   // 105
		 {  0, 1, 0, 1, 0,      10,    0,       0,    0 }    // 106
	};

	T = (Mjd_TT->Days + Mjd_TT->FracDay - MJD_J2000) / 36525.0;
	T2 = T * T;
	T3 = T2 * T;
	T4 = T3 * T;
	rev = 360.0 * 3600.0;


	l = fmod(485868.249036 + 1717915923.2178 * T + 31.8792 * T2 + 0.051635 * T3 - 0.00024470 * T4, rev) / Arcs;
	lp = fmod(1287104.79305 + 129596581.0481 * T - 0.5532 * T2 + 0.000136 * T3 - 0.00001149 * T4, rev) / Arcs;
	F = fmod(335779.526232 + 1739527262.8478 * T - 12.7512 * T2 - 0.001037 * T3 + 0.00000417 * T4, rev) / Arcs;
	D = fmod(1072260.70369 + 1602961601.2090 * T - 6.3706 * T2 + 0.006593 * T3 - 0.00003169 * T4, rev) / Arcs;
	Om = fmod(450160.398036 - 6962890.5431 * T + 7.4722 * T2 + 0.007702 * T3 - 0.00005939 * T4, rev) / Arcs;

	/* Nutation in longitude and obliquity [rad] */

	*deps = *dpsi = 0.0;

	for (i = 0; i < 106; i++)
	{
		arg = C[i][0] * l + C[i][1] * lp + C[i][2] * F + C[i][3] * D + C[i][4] * Om;
		*dpsi += (C[i][5] + C[i][6] * T) * sin(arg);
		*deps += (C[i][7] + C[i][8] * T) * cos(arg);
	}

	*dpsi = (1.0E-5 * *dpsi) / Arcs;
	*deps = (1.0E-5 * *deps) / Arcs;
}

void GHAMatrix(const MJDTIME* Mjd_UT1, double* dpsi, double Mat[])
{
	double Angle = 0;

	Angle = GAST(Mjd_UT1, dpsi);

	Rotation_z(Angle, Mat);
}

void PoleMatrix(double x, double y, double Mat[])
{
	double M1[9] = { 0 }, M2[9] = { 0 };

	Rotation_y(-x / Arcs, M1);
	Rotation_x(-y / Arcs, M2);

	MatrixMultiply(3, 3, 3, 3, M1, M2, Mat);
}

void ICRF_ITRF_Matrix(const double Mjd1, const MJDTIME* GT,
	int flag, double Mat[9])
{
	MJDTIME Mjd_UT1, Mjd_UTC, Mjd_TT;
	double dpsi;
	double Prec[9], Nut[9], GH[9], Pole[9];
	double tmp1[9], tmp2[9];

	static EOPPARA CurrEop;

	Mjd_TT.Days = GT->Days;
	Mjd_TT.FracDay = GT->FracDay - (BDT_TAI - TT_TAI) / SECPERDAY;

	Mjd_UTC.Days = Mjd_TT.Days;
	Mjd_UTC.FracDay = Mjd_TT.FracDay - TT_UTC() / SECPERDAY;
	InterposeEOP(&Mjd_UTC, &CurrEop);

	Mjd_UT1.Days = Mjd_UTC.Days;
	Mjd_UT1.FracDay = Mjd_UTC.FracDay + CurrEop.dUT1 / SECPERDAY;

	PrecMatrix(Mjd1, &Mjd_TT, Prec);
	NutMatrix(&Mjd_TT, &dpsi, Nut);
	GHAMatrix(&Mjd_UT1, &dpsi, GH);
	PoleMatrix(CurrEop.x, CurrEop.y, Pole);

	MatrixMultiply(3, 3, 3, 3, Nut, Prec, tmp1);
	MatrixMultiply(3, 3, 3, 3, Pole, GH, tmp2);
	MatrixMultiply(3, 3, 3, 3, tmp2, tmp1, Mat);   /*坐标转换矩阵*/

	if (flag == 0)     /*天球系——>地固系 */
	{
		MatrixTranspose(3, 3, Mat, tmp1);
		CopyArray(9, Mat, tmp1);
	}
}

double EqnEquinox(const MJDTIME* Mjd_TT, double* dpsi)
{
	double T = 0, T2 = 0, T3 = 0, rev = 0;
	double Om = 0, Angle = 0;

	T = (Mjd_TT->Days + Mjd_TT->FracDay - MJD_J2000) / 36525.0;
	T2 = T * T;
	T3 = T2 * T;
	rev = 360.0 * 3600.0;      /* arcsec/revolution */

	Om = (450160.280 - (5.0 * rev + 482890.539) * T +
		7.455 * T2 + 0.008 * T3) / Arcs;

	Angle = (0.00264 * sin(Om) + 0.000063 * sin(2 * Om)) / Arcs;

	return  *dpsi * cos(MeanObliquity(Mjd_TT)) + Angle;
}


double GMST(const MJDTIME* Mjd_UT1)
{
	double T0 = 0, T = 0, T2 = 0, T3 = 0, RetVal = 0;

	T0 = (Mjd_UT1->Days - MJD_J2000) / 36525.0;
	T = (Mjd_UT1->Days + Mjd_UT1->FracDay - MJD_J2000) / 36525;
	T2 = T * T;
	T3 = T * T2;

	RetVal = 24110.54841 + 8640184.812866 * T0 + 0.093104 * T2 - 6.2E-6 * T3;
	RetVal = RetVal / SECPERDAY + 1.002737909350795 * Mjd_UT1->FracDay;

	RetVal = fmod(RetVal, 1.0);

	if (RetVal < 0.0)
	{
		RetVal = RetVal + 1.0;
	}

	return (RetVal * pi * 2);
}

double GAST(const MJDTIME* Mjd_UT1, double* dpsi)
{
	double Res = 0;

	Res = GMST(Mjd_UT1) + EqnEquinox(Mjd_UT1, dpsi);

	Res = fmod(Res, pi * 2);

	return Res;
}

double TT_UTC()
{
	return (TT_TAI - BDT_TAI + BDT_UTC);
}