#pragma once

#include "util.h"
#include "GPSTime.h"
#include "RefSys.h"
#include "matrix.h"

#define false 0
#define true 1

#define MAXSATNUM 40     // �������ܵ�������
#define SRC_DATA_PATH "../"

#define DIM 6
#define Dim_STM  (6*(DIM+1))            // ����״̬����״̬ת�ƾ����ά��

// constance
#define C_Light 299792458.0      /* Speed of light  [m/s]; IAU 1976  */

// structure definition 
typedef struct CONSTSTATE{
	double	Clk[MAXSATNUM * 2], ClkCov[MAXSATNUM * MAXSATNUM * 4];
	double	Orb[MAXSATNUM * DIM], OrbCov[MAXSATNUM * MAXSATNUM * DIM * DIM];
}CONSTSTATE;

#define MAXCLKSER 40                    // ����Ӳ���ϵ�����Ԫ��
#define MAXANCHORNUM 5	 // ���ê��վ����
#define MAXANTENNA  6    // ê��վ�������������
typedef struct TIMESYC
{
	int     CurNum;                          // �Ӳ����еĵ�ǰֵ
	int     TotalNum;                        // �Ӳ����е�����
	double  ClkSeq[MAXCLKSER];               // �����Ӳ����У������Ӳ����Ӳ�仯�ʵĲ������
	GPSTIME Time[MAXCLKSER];                 // �Ӳ����ж�Ӧ��ʱ��
	double  Residual[MAXCLKSER];             // �����Ӳ���ϲв�
}TIMESYC;
typedef struct SCSTM
{
	GPSTIME   Time;
	double    Phi[Dim_STM];  // ��������λ�á��ٶȺ�״̬ת�ƾ���
	double    dPhi[Dim_STM];
	int      Valid;
}SCSTM;
typedef enum ANSTATEID { UNKNST, NOINIT, MANEUVER, BREAKING, ANSOK } ANSTATEID;
typedef struct SATINFO{
	POINT_STRUCT;
	unsigned char index;
	unsigned char id;
	double Tgd[4];              // ���ߵ�Ӳ���ӳ�
	GPSTIME TOE;             // ���������Ĳο�ʱ�䣬Ϊ�����ʼ֡ʱ�̵�GPSʱ
	double X[DIM];          // ����λ��3ά���ٶ�3ά
	double  dX[DIM];         // �˲�����ʱ�Ĺ���Ľ���
	double  Clk[2];          // �Ӳ��仯��
	double  CovX[DIM * DIM];   // ����״̬Э�������
	double  CovC[4];         // �����Ӳ�Э�������
	short   Health;          // ���ǽ���״��, �Ƿ������������ı��
	double  GapTime;         // û�н��в������µ�ʱ����
	double MeasStep;                        // �������µ�ʱ����
	double X_LSQ[3], Clk_LSQ, PDOP;         // ��С���˼��������
	double MeanClk_apr, StdClk_apr;         // ���ڳ�ʱ���жϺ�, ���������˲��������жϱ�־
	double MeanOrb_apr, StdOrb_apr;
	double MeanClk_pst, StdClk_pst;         // ���ڳ�ʱ���жϺ�, ���
	double MeanOrb_pst, StdOrb_pst;
	float   URA;             // ��������ľ���ָ��
	ANSTATEID  Valid;                       // �������,û�г�ʼ��Ϊ0,�ɹ���ʼ��Ϊtrue
	SCSTM    STM[2];                        // ���֡��ʼʱ��0-120s�Ķ˵�Ļ���״̬
	TIMESYC SatClk;                         // ʱ��ͬ���ṹ��
}SATINFO;

typedef struct ANCHSTN {
	POINT_STRUCT;
	short StnId;
	//short AtnId;
	double clk;
	//double Elv;
	double Pos[6];
	double Tgd[2];
	short RefClkFlag;
	short Valid;
}ANCHSTN;

typedef struct DERISLOMC
{
	double ApriClkResid;        // ��ǰ�Ӳ�в�
	double PostClkResid;        // ����Ӳ�в�
	double ApriOrbResid;        // ��ǰ��������в�
	double PostOrbResid;        // �����������в�
	double dt[2];               // �Ӳ���ص�ʱ������t2+t1-2*t0��dt[0]��ӦPrn1���ǵ�ֵ��dt[1]��ӦPrn2���ǵ�ֵ
}DERISLOMC;

typedef struct ISLPROBS {
	LIST_STRUCT;
	short   TrScid;          // �����źŵ�����SCID
	short   RvScid;          // �����źŵ�����SCID
	short	TrAnt;
	short	RvAnt;
	GPSTIME RvLocTime;       // ����ʱ�̵������ӵı���ʱ
	double  PRObs;           // �Ǽ����۲�ֵ
	//unsigned short Quality;         // �۲����ݵ�����, δ֪Ϊ0
	int    Valid;                  // �Ǽ�α��۲�ֵ����Ч��

	double Corr[10];     // 0: 14����ģ�͸����ĵ�����ӳ�, 1: ����������ο�ֵ, 2: ���������������ο�ֵ����
	//// 3: ���峱�����ο�ֵ����, 4-5:����������λ���ĸ����ο�ֵ,��X73B [0]��ʾ�ź�Դ��[1]:��ʾ�ź��ն� ,�������Ƿ���[0],[1]��0;
	//// 6-7: ��������ʱ��ֵ,����������һ��, 8-9:����ֵ1
	////struct ISLPROBS* next;
}ISLPROBS;

typedef struct DEROBS
{
	EDGE_STRUCT;
	short Scid1;       // ������������SCID��PRN��
	short Scid2;       // ��Ϊ��׼�Ĳο�����SCID��PRN��
	short TrAnt, RvAnt;      // ����ͽ����źŵ����ߺ�
	//short I1, I2;            // I1Ϊ���ǽ����źŵĹ۲�ֵ���ڵ������ţ�I2Ϊ�ο��ǵ�������
	ISLPROBS *isl1, *isl2;
	GPSTIME T1, T2;          // T1Ϊ���ǲ���ʱ�̣�T2Ϊ�ο��ǵĲ���ʱ�̣���Ϊ�����ӱ���ʱ
	double DerCObs;          // �������Ӳ�۲�ֵ��m��
	double CCObs;            // �������Ӳ�۲�ֵ�ļ���ֵ[m]
	double CConCorr;         // ���������ЧӦ���Ӳ�仯�����Ӳ��������[m]
	double DerANObs;         // ��������������۲�ֵ[m]
	double CANObs;           // ��������������۲�ֵ�ļ���ֵ[m]
	double AConCorr;         // ������������Ӳ�仯���ĳ������������[m]
	double P1RvState[Dim_STM];    // ���������źŽ���ʱ�̵�״̬��״̬ת�ƾ���
	double P1TrState[Dim_STM];    // ���������źŷ���ʱ�̵�״̬��״̬ת�ƾ���
	double P2TrState[Dim_STM];    // �ο������źŷ���ʱ�̵�״̬
	double P2RvState[Dim_STM];    // �ο������źŽ���ʱ�̵�״̬

	DERISLOMC ObsQua;            // �۲����ݵ�����ָ��
	//unsigned short  Quality;      // �����۲�ֵ���������,��ԭʼ�۲����ݻ�ȡ
	short  Valid;                 // ������Ϊ-1�����ɵ����۲�ֵ�������ֵΪ0
								  // ����̽������Ϊ1���ֲ�Ϊ-1, �ظ�����Ϊ3
	short ClkBlunder[2], OrbBlunder[2];        // [0]Prn1,[1]Prn2��������ֲ�Ϊ2������Ϊ1��δ���Ϊ0
}DEROBS;

typedef struct SATNET {
	GRAPH_STRUCT;
	int ObsNum; // number of links
	CONSTSTATE AllSatCov;
}SATNET;

typedef struct SATINDEX {
	int SID;
	SATINFO* satinfo;
}SATINDEX;


#define max(a,b) a>b?a:b

int init();
void run();

int SearchSatIndex(const int SID);

int ReadSimObsData(GPSTIME* Time, ISLPROBS* islist);
double GetPseudoRange(int tid, int w, double sec, double rpv[6], double tpv[6]);
int OpenANSResFileDaily(GPSTIME* Time);
int AssignEpkISLObs(GPSTIME* Time, ISLPROBS* il, ISLPROBS* EpkObs);
void TimeUpdate(const GPSTIME* Time, SATNET* SatNet);
void InitStateTranMatrix(int row, int col, double STM[]);
void CompStateNoiseCov(const double Step, const ANSTATEID Valid, double Q[]);
int GenDerPrObs(ISLPROBS* EpkObs, SATNET* SatNet);
int GenDerObsPredict(DEROBS* derobs);
void DectectDerObsOutlier(SATNET* SatNet);
int ScalarTimeMeasUpdate(double O_C, double sigma2, double H[], int Scale, CONSTSTATE* AllSatCov);
int CompVectStat(const int n, const double Dat[], double* Mean, double* Std, SATINFO* SatAtod);
void ClkMeasUpdate(SATNET* SatNet);
ScalarOrbitMeasUpdate(double O_C, double sigma2, double H[], int Scale, CONSTSTATE* AllSatCov);
void ANSMeasUpdate(SATNET* SatNet);
void OutputSatOrbit(SATINFO* SatAtod);

void OrbitIntegToGivenTime(const GPSTIME* BegGPS, const GPSTIME* GivenTime, const double Pace, double Y0[]);