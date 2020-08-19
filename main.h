#pragma once

#include "util.h"
#include "GPSTime.h"
#include "RefSys.h"
#include "matrix.h"

#define false 0
#define true 1

#define MAXSATNUM 40     // 星座中总的卫星数
#define SRC_DATA_PATH "../../"

#define DIM 6
#define Dim_STM  (6*(DIM+1))            // 卫星状态及其状态转移矩阵的维数

// constance
#define C_Light 299792458.0      /* Speed of light  [m/s]; IAU 1976  */

// structure definition 
typedef struct CONSTSTATE{
	double	Clk[MAXSATNUM * 2], ClkCov[MAXSATNUM * MAXSATNUM * 4];
	double	Orb[MAXSATNUM * DIM], OrbCov[MAXSATNUM * MAXSATNUM * DIM * DIM];
}CONSTSTATE;

#define MAXCLKSER 40                    // 最大钟差拟合的序列元素
#define MAXANCHORNUM 5	 // 最大锚固站数量
#define MAXANTENNA  6    // 锚固站的最大天线数量
typedef struct TIMESYC
{
	int     CurNum;                          // 钟差序列的当前值
	int     TotalNum;                        // 钟差序列的总数
	double  ClkSeq[MAXCLKSER];               // 保存钟差序列，用于钟差与钟差变化率的参数拟合
	GPSTIME Time[MAXCLKSER];                 // 钟差序列对应的时间
	double  Residual[MAXCLKSER];             // 保存钟差拟合残差
}TIMESYC;
typedef struct SCSTM
{
	GPSTIME   Time;
	double    Phi[Dim_STM];  // 包含卫星位置、速度和状态转移矩阵
	double    dPhi[Dim_STM];
	int      Valid;
}SCSTM;
typedef enum ANSTATEID { UNKNST, NOINIT, MANEUVER, BREAKING, ANSOK } ANSTATEID;
typedef struct SATINFO{
	POINT_STRUCT;
	unsigned char index;
	unsigned char id;
	double Tgd[4];              // 天线的硬件延迟
	GPSTIME TOE;             // 卫星星历的参考时间，为测距起始帧时刻的GPS时
	double X[DIM];          // 卫星位置3维、速度3维
	double  dX[DIM];         // 滤波更新时的轨道改进量
	double  Clk[2];          // 钟差及其变化率
	double  CovX[DIM * DIM];   // 卫星状态协方差矩阵
	double  CovC[4];         // 卫星钟差协方差矩阵
	short   Health;          // 卫星健康状况, 是否参与自主定轨的标记
	double  GapTime;         // 没有进行测量更新的时间间隔
	double MeasStep;                        // 测量更新的时间间隔
	double X_LSQ[3], Clk_LSQ, PDOP;         // 最小二乘计算的坐标
	double MeanClk_apr, StdClk_apr;         // 用于长时间中断后, 自主定轨滤波收敛的判断标志
	double MeanOrb_apr, StdOrb_apr;
	double MeanClk_pst, StdClk_pst;         // 用于长时间中断后, 验后
	double MeanOrb_pst, StdOrb_pst;
	float   URA;             // 自主定轨的精度指标
	ANSTATEID  Valid;                       // 可用情况,没有初始化为0,成功初始化为true
	SCSTM    STM[2];                        // 测距帧起始时刻0-120s的端点的积分状态
	TIMESYC SatClk;                         // 时间同步结构体
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
	double ApriClkResid;        // 验前钟差残差
	double PostClkResid;        // 验后钟差残差
	double ApriOrbResid;        // 验前自主定轨残差
	double PostOrbResid;        // 验后自主定轨残差
	double dt[2];               // 钟差相关的时间间隔，t2+t1-2*t0，dt[0]对应Prn1卫星的值；dt[1]对应Prn2卫星的值
}DERISLOMC;

typedef struct ISLPROBS {
	LIST_STRUCT;
	short   TrScid;          // 发射信号的卫星SCID
	short   RvScid;          // 接收信号的卫星SCID
	short	TrAnt;
	short	RvAnt;
	GPSTIME RvLocTime;       // 测量时刻的卫星钟的表面时
	double  PRObs;           // 星间距离观测值
	//unsigned short Quality;         // 观测数据的质量, 未知为0
	int    Valid;                  // 星间伪距观测值的有效性

	double Corr[10];     // 0: 14参数模型改正的电离层延迟, 1: 对流层改正参考值, 2: 相对论周期项改正参考值，米
	//// 3: 固体潮改正参考值，米, 4-5:卫星天线相位中心改正参考值,对X73B [0]表示信号源，[1]:表示信号终端 ,其它均是放在[0],[1]填0;
	//// 6-7: 卫星天线时延值,与卫星天线一致, 8-9:理论值1
	////struct ISLPROBS* next;
}ISLPROBS;

typedef struct DEROBS {
	EDGE_STRUCT;
	short Scid1;       // 自主导航卫星SCID与PRN号
	short Scid2;       // 作为基准的参考卫星SCID与PRN号
	short TrAnt, RvAnt;      // 发射和接收信号的天线号
	//short I1, I2;            // I1为本星接收信号的观测值所在的索引号，I2为参考星的索引号
	ISLPROBS *isl1, *isl2;
	GPSTIME T1, T2;          // T1为本星测量时刻，T2为参考星的测量时刻，均为卫星钟表面时
	double DerCObs;          // 导出的钟差观测值【m】
	double CCObs;            // 导出的钟差观测值的计算值[m]
	double CConCorr;         // 包括相对论效应与钟差变化量的钟差常量改正项[m]
	double DerANObs;         // 导出的自主定轨观测值[m]
	double CANObs;           // 导出的自主定轨观测值的计算值[m]
	double AConCorr;         // 包括相对论与钟差变化量的常量轨道改正项[m]
	double P1RvState[Dim_STM];    // 自主星在信号接收时刻的状态和状态转移矩阵
	double P1TrState[Dim_STM];    // 自主星在信号发射时刻的状态和状态转移矩阵
	double P2TrState[Dim_STM];    // 参考星在信号发射时刻的状态
	double P2RvState[Dim_STM];    // 参考星在信号接收时刻的状态

	DERISLOMC ObsQua;            // 观测数据的质量指标
	//unsigned short  Quality;      // 导出观测值的质量标记,从原始观测数据获取
	short  Valid;                 // 不可用为-1，生成导出观测值及其计算值为0
								  // 数据探测正常为1，粗差为-1, 重复数据为3
	short ClkBlunder[2], OrbBlunder[2];        // [0]Prn1,[1]Prn2检测结果，粗差为2，正常为1，未检查为0
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
int ScalarOrbitMeasUpdate(double O_C, double sigma2, double H[], int Scale, CONSTSTATE* AllSatCov);
void ANSMeasUpdate(SATNET* SatNet);
void OutputSatOrbit(SATINFO* SatAtod);

void OrbitIntegToGivenTime(const GPSTIME* BegGPS, const GPSTIME* GivenTime, const double Pace, double Y0[]);