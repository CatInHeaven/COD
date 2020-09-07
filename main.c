#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "main.h"
#include <sys/time.h>
#include <mpi.h>

FILE* FDerOrb, * FDerClk, * FPProc, * FPLOG, * FDerMAO;
FILE* Fres;
GPSTIME    FrameTime;
// configural iformation
int        SatNum;
int		AncNum = 1;

//SATINDEX        SatID[MAXSATNUM][2] = { {19,},{20,},{21,},{22,},{23,},{24,},{25,},{26,},{27,},{28,},{29,},{30,},{36,},{37,},{38,},{39,},{40,},{41,},{42,},{43,},{44,},{45,},{46,},{47,} };
int   IDList[MAXSATNUM][2] = { {19,0},{20,1},{21,2},{22,0},{23,1},{24,2},{25,0},{26,1},{27,2},{28,0},{29,1},{30,2},{36,0},{37,1},{38,2},{39,0},{40,1},{41,2},{42,0},{43,1},{44,2},{45,0},{46,1},{47,2},{76,0} };
double     PosAccu, VelAccu;         

double StepOfAutoNav = 300;
double NoiseOfISL = 5.0;
double OrbDiff = 100;
double ClkDiff = 50;
int RunFlag = 0;
int TotalSatNum = 24;
// satellite network
SATNET* SatNet;

// nodes in graph
SATINFO* SatNodes;
ANCHSTN* Anchor;
// edges in graph
DEROBS* EpkDerObs;

// simulate link data
ISLPROBS* SimData;
//ISLPROBS* EpkISLObs;

// MPI Environment
int world_size;
int world_rank;

int main() {
	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	//if (init() < 0)	return -1;
	init();
	int n = 0;
	struct timeval tv_begin;
	struct timeval tv_end;
	double dDuration;
	gettimeofday(&tv_begin, NULL);
	while (n < 10) {
		run();
		printf("%4d %10.1f\n", FrameTime.Week, FrameTime.SecOfWeek);
		n++;
	}
	gettimeofday(&tv_end, NULL);
	dDuration = 1000 * (tv_end.tv_sec - tv_begin.tv_sec) + ((tv_end.tv_usec - tv_begin.tv_usec) / 1000.0);
	printf("Process %d Duration: %f ms\n",world_rank, dDuration);
	// TODO memory free!
	DestroySpraseMatrix(SatNet->AllSatCov.ClkCov);
	DestroySpraseMatrix(SatNet->AllSatCov.OrbCov);
	graph_destroy((Graph*)SatNet);
	if (world_rank == 0) list_destroy((list*)SimData);
	fclose(FDerOrb);
	fclose(FDerClk);
	fclose(FPProc);
	fclose(FPLOG);
	fclose(FDerMAO);
	MPI_Finalize();
	printf("I'm %d, I finished here.\n",world_rank);
	return 0;
}

int init() {
	int send_rank;
	int receive_rank;
	receive_rank = 0; // main process

	int i;

	// start time
	COMMONTIME StartDate;
	StartDate.Year = 2018;
	StartDate.Month = 1;
	StartDate.Day = 21;
	StartDate.Hour = StartDate.Minute = 0;
	StartDate.Second = 0.0;
	CommonTimeToGPSTime(&StartDate, &FrameTime);


	// configuration
	//SatNum = 24;
	PosAccu = 2.0;
	VelAccu = 0.05;
	StepOfAutoNav = 300.0;

	// initialize satellite network
	if ((SatNet = (SATNET*)malloc(sizeof(SATNET))) == NULL) return -1;
	if ((SatNodes = (SATINFO*)malloc(sizeof(SATINFO))) == NULL) return -1;
	SatNet->points = (Point*)SatNodes;
	if ((EpkDerObs = (DEROBS*)malloc(sizeof(DEROBS))) == NULL) return -1;
	SatNet->edges = (Edge*)EpkDerObs;
	if ((SatNet->border_points = (Point*)malloc(sizeof(SATINFO))) == NULL) return -1;
	graph_init((Graph*)SatNet);

	SATINFO* Info;
	// read satllite information
	if (world_rank == 0) {
		char FName[256] = { 0 };
		FILE* Frin;
		char Line[512];
		int SCID;
		unsigned char id;
		int index;

		sprintf(FName, "../../data/input/InitState0629000000.txt");
		if ((Frin = fopen(FName, "rt")) == NULL)	return -1;

		double x[DIM];
		double clk[2];
		short health;
		GPSTIME toe;
		double SigmaPos, SigmaVel, SigmaClkoff, SigmaClkSht;

		SigmaPos = PosAccu;
		SigmaVel = VelAccu;
		SigmaClkoff = PosAccu / C_Light;
		SigmaClkSht = VelAccu / C_Light;

		while (!feof(Frin)){
			fgets(Line, 512, Frin);
			if (strncmp(Line, "EOF", 3) == 0)   break;

			if (sscanf(Line, "%d %*d %hd %lf %lf %lf %lf %lf %lf %lf %lf %lf %hd", &SCID, &(toe.Week), &(toe.SecOfWeek), x, x + 1, x + 2, x + 3, x + 4, x + 5, clk, clk + 1, &health) != 12)
				continue;

			id = (unsigned char)SCID;
			if ((send_rank = GetPointRank(id)) < 0)	continue;
			//if ((index = SearchSatIndex(id)) < 0) continue;
			if ((Info = (SATINFO*)malloc(sizeof(SATINFO))) == NULL) return -1;
			memset(Info, 0, sizeof(SATINFO));
			Info->id = id;
			//TODO: other information initialization
			//SatID[index]->satinfo = Info;
			Info->URA = 2;
			Info->type = 1;
			Info->Valid = ANSOK;
			Info->MeasStep = 300;
			Info->Health = health;
			for (i = 0; i < DIM; i++){
				Info->X[i] = x[i];
			}
			Info->Clk[0] = clk[0];
			Info->Clk[1] = clk[1];
			Info->TOE = toe;

			//initialize the covariance of clk and X
			//diagonal matrix
			Info->CovC[0] = pow(SigmaClkoff, 2.0);
			Info->CovC[3] = pow(SigmaClkSht, 2.0);
			for (i = 0; i < 3; i++)      Info->CovX[i * (DIM + 1)] = pow(SigmaPos, 2.0);
			for (i = 3; i < 6; i++)      Info->CovX[i * (DIM + 1)] = pow(SigmaVel, 2.0);

			MPI_Send((char*)Info, sizeof(SATINFO), MPI_BYTE, 1, 0, MPI_COMM_WORLD);
			MPI_Send((char*)Info, sizeof(SATINFO), MPI_BYTE, 2, 0, MPI_COMM_WORLD);
			if (send_rank != 0) {
				border_point_append((Graph*)SatNet, (Point*)Info);
				//free(Info);
			}
			else {
				Info->rank = world_rank;
				point_append((Graph*)SatNet, (Point*)Info);
			}
		}
		fclose(Frin);
		MPI_Send(Info, 0, MPI_BYTE, 1, 0, MPI_COMM_WORLD);
		MPI_Send(Info, 0, MPI_BYTE, 2, 0, MPI_COMM_WORLD);
	}
	else {
		MPI_Status status;
		int count;
		while (1) {
			if ((Info = (SATINFO*)malloc(sizeof(SATINFO))) == NULL) return -1;
			MPI_Recv((char*)Info, sizeof(SATINFO), MPI_BYTE, 0, 0, MPI_COMM_WORLD, &status);
			MPI_Get_count(&status, MPI_BYTE, &count);
			if (count <= 0) break; 

			if ((receive_rank = GetPointRank(Info->id)) < 0)	continue;

			if (receive_rank == world_rank) {
				Info->rank = world_rank;
				point_append((Graph*)SatNet, (Point*)Info);
			}
			else {
				border_point_append((Graph*)SatNet, (Point*)Info);
			}
		}
		free(Info);
	}

	SatNum = SatNet->points_num;
	CONSTSTATE* AllSatCov;
	AllSatCov = &SatNet->AllSatCov;
	memset(AllSatCov->Clk, 0, sizeof(double) * MAXSATNUM * 2);
	memset(AllSatCov->Orb, 0, sizeof(double) * MAXSATNUM * DIM);
	AllSatCov->ClkCov = (SpraseMatrix*)malloc(sizeof(SpraseMatrix));
	AllSatCov->OrbCov = (SpraseMatrix*)malloc(sizeof(SpraseMatrix));
	CreateSpraseMatrix(AllSatCov->ClkCov, SatNum, SatNum, 2, 2);
	CreateSpraseMatrix(AllSatCov->OrbCov, SatNum, SatNum, DIM, DIM);
	SATINFO* sat;
	sat = (SATINFO*)SatNet->points;
	for (i = 0; i < SatNum; i++) {
		sat = (SATINFO*)sat->next;
		CopyArray(2, AllSatCov->Clk + 2 * i, sat->Clk);
		CopySubSpraseMatrix(AllSatCov->OrbCov, i, i, DIM, DIM, sat->CovX);
		CopySubSpraseMatrix(AllSatCov->ClkCov, i, i, 2, 2, sat->CovC);
	}

	if (world_rank == 0) {
		// initialize anchor station information
		if ((Anchor = (ANCHSTN*)malloc(sizeof(ANCHSTN))) == NULL) return -1;
		point_append((Graph*)SatNet, (Point*)Anchor);
		Anchor[0].StnId = 76;
		Anchor[0].type = 0;
		Anchor[0].rank = world_rank;
		Anchor[0].clk = 0.0;
		Anchor[0].Pos[0] = -2003066.221396;
		Anchor[0].Pos[1] = 5716340.837106;
		Anchor[0].Pos[2] = 1991308.673030;
		Anchor[0].Tgd[0] = 0.0;
		Anchor[0].Tgd[1] = 0.0;
		Anchor[0].RefClkFlag = 1;
		Anchor[0].Valid = 1;

		// read simulate data (small dataset)
		if ((SimData = (ISLPROBS*)malloc(sizeof(ISLPROBS))) == NULL) return -1;
		//if ((EpkISLObs = (ISLPROBS*)malloc(sizeof(ISLPROBS))) == NULL) return -1;
		list_init((list*)SimData);
		//list_init((list*)EpkISLObs);
		if (ReadSimObsData(&FrameTime, SimData) < 0) return -1;
	}
	return 0;
}

int SearchSatIndex(int id) {
	int i;
	for (i = 0; i < SatNum; i++) {
		if (id == IDList[i][0])
			return i;
	}
	return -1;
}

char* GetPoint(int id) {
	int index;
	char* node = NULL;
	Point* iter;
	for (iter = (Point*)SatNet->points->next; iter != SatNet->points; iter = (Point*)iter->next){
		if (iter->type == 1) {
			//printf("%d %d %d\n",id, ((SATINFO*)iter)->id, world_rank);
			if (((SATINFO*)iter)->id == id) {
				node = (char*)iter;
				break;
			}
			else continue;
			
		}
		if (iter->type == 0) {
			if (((ANCHSTN*)iter)->StnId == id) {
				node = (char*)iter;
				break;
			}
		}
		
	}
	return node;
}

char* GetBorderPoint(int id) {
	int index;
	char* node = NULL;
	Point* iter;
	for (iter = (Point*)SatNet->border_points->next; iter != SatNet->border_points; iter = (Point*)iter->next){
		if (iter->type == 1) {
			//printf("%d %d %d\n",id, ((SATINFO*)iter)->id, world_rank);
			if (((SATINFO*)iter)->id == id) {
				node = (char*)iter;
				break;
			}
			else continue;
			
		}
		if (iter->type == 0) {
			if (((ANCHSTN*)iter)->StnId == id) {
				node = (char*)iter;
				break;
			}
		}
		
	}
	return node;
}

int GetPointRank(int id) {
	int i;
	for (i = 0; i < TotalSatNum + AncNum; i++) {
		if (id == IDList[i][0])
			return IDList[i][1];
	}
	printf("can't find id = %d\n",id);
	return -1;
}

void run() {
	int i;
	GPSTIME NextFramTime;
	NextFramTime.Week = FrameTime.Week;
	NextFramTime.SecOfWeek = FrameTime.SecOfWeek + StepOfAutoNav;
	CheckGPSTime(&NextFramTime);

	// output files
	if (fmod(NextFramTime.SecOfWeek + 0.001, SECPERDAY) < 0.5 || RunFlag == 0) {

		OpenANSResFileDaily(&NextFramTime);
	}

	TimeUpdate(&NextFramTime, SatNet);

	// add link to Satellite Network
	AssignEpkISLObs(&NextFramTime, SimData, SatNet);
	GenDerPrObs(SatNet);

	printf("开始进行数据预处理\n");
	DectectDerObsOutlier(SatNet);
	printf("数据预处理结束\n");

	printf("开始进行钟差测量更新\n");
	ClkMeasUpdate(SatNet);
	printf("完成钟差测量更新\n");

	printf("开始进行钟轨道测量更新\n");
	ANSMeasUpdate(SatNet);
	printf("完成轨道测量更新\n");

	SATINFO* sat;
	sat = (SATINFO*)SatNet->points;
	for (i = 0; i < SatNum; i++) {
		sat = (SATINFO*)sat->next;
		if (sat->Valid <= NOINIT)   continue;
		OutputSatOrbit(sat);
	}

	FrameTime.Week = NextFramTime.Week;
	FrameTime.SecOfWeek = NextFramTime.SecOfWeek;
	RunFlag = 1;
}

int ReadSimObsData(GPSTIME* Time, ISLPROBS* islist) {
	char Line[1024];
	char DatFileName[256];
	FILE* in;
	GPSTIME CurTime;
	ISLPROBS* obs;

	CurTime.Week = Time->Week;
	CurTime.SecOfWeek = Time->SecOfWeek;

	memset(DatFileName, 0, sizeof(DatFileName));
	//sprintf(DatFileName, "%sSimu_%d.txt", SRC_DATA_PATH, (int)((CurTime.Week - 629) * 7 + CurTime.SecOfWeek / SECPERDAY));

	sprintf(DatFileName, "%ssimu_0.txt", SRC_DATA_PATH);

	if ((in = fopen(DatFileName, "rt")) == NULL)    return -1;

	memset(Line, 0, sizeof(Line));
	while (fgets(Line, 1024, in) != NULL) {
		obs = (ISLPROBS*)malloc(sizeof(ISLPROBS));
		if (sscanf(Line, "%hd %lf %hd %hd %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
			&obs->RvLocTime.Week, &obs->RvLocTime.SecOfWeek, &obs->TrScid, &obs->RvScid,
			&obs->tpv[0], &obs->tpv[1], &obs->tpv[2], &obs->tpv[3], &obs->tpv[4], &obs->tpv[5],
			&obs->rpv[0], &obs->rpv[1], &obs->rpv[2], &obs->rpv[3], &obs->rpv[4], &obs->rpv[5]) < 16)  continue;

		obs->TrAnt = obs->RvAnt = 1;
		if (obs->TrScid == 76) {
			obs->TrAnt = 6;
			ICRF_ITRF_GPST(MJD_J2000, &obs->RvLocTime, 0, obs->tpv, Anchor[0].Pos);
		}
		if (obs->RvScid == 76) {
			obs->RvAnt = 6;
			ICRF_ITRF_GPST(MJD_J2000, &obs->RvLocTime, 0, obs->rpv, Anchor[0].Pos);
		}

		//obs->PRObs = GetPseudoRange(obs->TrScid, obs->RvLocTime.Week, obs->RvLocTime.SecOfWeek, rpv, tpv);

		//obs->Valid = true;
		for (int i = 0; i < 10; i++)  *(obs->Corr + i) = 0.0;
		//if (obs->TrScid != 0 && obs->RvScid != 0 && fabs(obs->PRObs) > 0.1)
		list_append((list*)islist, (node_t*)obs);
	}
	return 0;
}


#define GM_D 4
#define GM_Earth   398600.4415e+9     /* [m^3/s^2]; JGM3  */
#define R_Earth    6378136.3          /* Radius Earth [m]; JGM3  */
static double EGM2008[12 + 1][12 + 1] = {   //EGM2008 12*12    non-normalization coefficients
{  1.00000000000000e+00,  0.00000000000000e+00,  1.78727064852404e-09,  2.68087089400898e-07,
  -4.49265432143808e-07, -8.08134749120456e-08,  2.08597340829085e-08,  6.96250561212826e-08,
   4.04734040507365e-08,  1.39057297289347e-08, -8.10038372974958e-08, -1.60117307756400e-08,
  -2.44377108803877e-08 },
{  0.00000000000000e+00,  0.00000000000000e+00, -9.03872789196567e-07, -2.11430620933483e-07,
   1.48135037248860e-07, -5.23297729692969e-08, -4.65006396354240e-08,  9.26271584963605e-09,
   5.36176396348744e-09, -2.19566872878597e-09, -3.04873277966174e-09, -5.12574744185821e-09,
   1.41851840209032e-09 },
{ -1.08262617385222e-03, -2.66739475237484e-10,  1.57461532572292e-06,  1.97221581835718e-07,
  -1.20094612622961e-08, -7.10091727223769e-09,  1.85610014415166e-10, -3.05830741090084e-09,
  -8.69097506975657e-10, -5.61279122373972e-10, -8.98700350298437e-10, -6.86500851352006e-10,
   9.33548734028760e-11 },
{  2.53241051856772e-06,  2.19314963131333e-06,  3.09043900391649e-07,  1.00583513408823e-07,
   6.52467287187687e-09,  3.87811503746114e-10, -1.78456878985903e-09, -2.63441737188836e-10,
   9.11123256683835e-11,  1.70315311872616e-11, -4.65429878797259e-11, -2.68497529620093e-11,
   1.19142578228308e-12 },
{  1.61989759991697e-06, -5.08643560439584e-07,  7.83745457404552e-08,  5.92150177639666e-08,
  -3.98320424873187e-09, -1.64817193269055e-09, -4.32985145712914e-10,  6.34517062883340e-12,
   1.61456958686981e-11, -5.52719054744610e-12, -3.14219976112891e-12,  1.97304332557368e-12,
   2.02043921993195e-13 },
{  2.27753590730836e-07, -5.38824899516695e-08,  1.05528866715599e-07, -1.49264876546633e-08,
  -2.29951140350422e-09,  4.30428041979182e-10, -5.53052835075614e-11,  1.05362788327707e-11,
   8.62847180421363e-12,  2.94407727172148e-12, -5.53673660391750e-13,  1.35018682098682e-13,
   9.24324908420780e-14 },
{ -5.40666576283813e-07, -5.97343298033422e-08,  6.05208460634626e-09,  1.18691485343662e-09,
  -3.25640750330606e-10, -2.15620819339786e-10,  2.20647849206413e-12,  4.47195472935274e-13,
   3.81759103462575e-13, -1.84722722000921e-13, -2.56642580235160e-15, -3.73002990453661e-14,
   7.93456619902649e-15 },
{  3.50551795713742e-07,  2.05588639629593e-07,  3.29094239021408e-08,  3.52793316982768e-09,
  -5.83957862236059e-10,  5.83168243827302e-13, -2.49041082423985e-11,  2.79642779894244e-14,
   1.53675163200403e-13, -9.82332112412755e-16, -1.05051949374750e-14,  1.16913366043005e-15,
   3.76136664795241e-16 },
{  2.03993125929884e-07,  1.59157368608999e-08,  6.57191789050830e-09, -1.95877062127080e-10,
  -3.18938954974722e-10, -4.65186864900789e-12, -1.84231113950009e-12,  3.42946981235784e-13,
  -1.58099717990021e-13,  7.46375383846391e-15, -7.05044513474636e-16,  2.58693871035618e-16,
   6.04892021131941e-17 },
{  1.22127958919496e-07,  9.23680159834313e-08,  1.48332354779638e-09, -1.21385977402798e-09,
  -8.01426596948314e-12, -1.66854476583190e-12,  8.29075183256221e-13, -2.24863950804731e-13,
   6.14935949868439e-14, -3.66382151027340e-15, -9.91347388867048e-17, -1.74797497730406e-17,
   9.22857076268619e-18 },
{  2.44390769772693e-07,  5.17579366855752e-08, -5.58850515750517e-09, -4.08543277838639e-11,
  -4.97504625744296e-11, -3.05998317245760e-12, -2.60875859956474e-13,  6.95434264606992e-15,
   4.65024152588070e-15,  2.32966745332555e-15,  4.17303195090727e-16, -1.40944600782318e-17,
  -2.80823506805120e-19 },
{ -2.43476590998147e-07,  9.21662362461700e-09,  1.04137815696613e-09, -1.41037594062331e-10,
  -1.59791999662987e-11,  1.48877935544880e-12, -6.16244289111870e-15,  1.93284906670200e-15,
  -3.00170318252836e-16, -1.91078571621793e-16, -4.95732591829582e-17,  9.35314409371331e-18,
  -9.96392594151775e-20 },
{  1.82180961307286e-07, -3.03368874316361e-08,  6.50852586593136e-10,  1.47585448435234e-10,
  -2.10235803667400e-11,  8.21880584565280e-13,  7.43205597867172e-15, -4.23120659997131e-15,
  -5.74916521515487e-16,  1.01567724286930e-16, -1.84917037330171e-18,  4.99790432119338e-19,
  -2.17582516336317e-20 }
};
void AccelHarmonic(const double Pos[], const double E[9], double Acc[3],
	int HasPDev, double dadr[9])
{

	int     n, m;                           /* Loop counters */
	double  r_sqr, rho, Fac;               /* Auxiliary quantities  */
	double  x0, y0, z0;                      /* Normalized coordinates */
	double  ax, ay, az;                      /* Acceleration vector  */
	double  C, S;                           /* Gravitational coefficients */
	double  r_bf[3];                       /* Body-fixed position */
	double  a_bf[3];                       /* Body-fixed acceleration */
	double  ET[9];
	double  V[GM_D + 2][GM_D + 2] = { 0.0 };            /* Harmonic functions */
	double  W[GM_D + 2][GM_D + 2] = { 0.0 };            /* work array (0..n_max+1,0..n_max+1) */
	double  GM_R2 = GM_Earth / R_Earth / R_Earth;

	/* Body-fixed position  */

	MatrixMultiply(3, 3, 3, 1, E, Pos, r_bf);

	/* Auxiliary quantities  */

	r_sqr = VectDot(3, 3, r_bf, r_bf);        /* Square of distance  */
	rho = R_Earth * R_Earth / r_sqr;

	x0 = R_Earth * r_bf[0] / r_sqr;          /* Normalized */
	y0 = R_Earth * r_bf[1] / r_sqr;          /* coordinates */
	z0 = R_Earth * r_bf[2] / r_sqr;

	/* Calculate zonal terms V(n,0); set W(n,0)=0.0  */

	V[0][0] = R_Earth / sqrt(r_sqr);
	W[0][0] = 0.0;

	V[1][0] = z0 * V[0][0];
	W[1][0] = 0.0;

	for (n = 2; n <= GM_D + 1; n++)
	{
		V[n][0] = ((2 * n - 1) * z0 * V[n - 1][0] - (n - 1) * rho * V[n - 2][0]) / n;
		W[n][0] = 0.0;
	}

	/* Calculate tesseral and sectorial terms  */

	for (m = 1; m <= GM_D + 1; m++)
	{
		/* Calculate V(m,m) .. V(n_max+1,m)  */

		V[m][m] = (2 * m - 1) * (x0 * V[m - 1][m - 1] - y0 * W[m - 1][m - 1]);
		W[m][m] = (2 * m - 1) * (x0 * W[m - 1][m - 1] + y0 * V[m - 1][m - 1]);

		if (m <= GM_D)
		{
			V[m + 1][m] = (2 * m + 1) * z0 * V[m][m];
			W[m + 1][m] = (2 * m + 1) * z0 * W[m][m];
		}

		for (n = m + 2; n <= GM_D + 1; n++)
		{
			V[n][m] = ((2 * n - 1) * z0 * V[n - 1][m] - (n + m - 1) * rho * V[n - 2][m]) / (n - m);
			W[n][m] = ((2 * n - 1) * z0 * W[n - 1][m] - (n + m - 1) * rho * W[n - 2][m]) / (n - m);
		}
	}

	/* Calculate accelerations ax,ay,az  */

	ax = ay = az = 0.0;

	for (m = 0; m <= GM_D; m++)
	{
		for (n = m; n <= GM_D; n++)
		{
			if (m == 0)
			{
				C = EGM2008[n][0];
				ax = ax - C * V[n + 1][1];
				ay = ay - C * W[n + 1][1];
				az = az - (n + 1) * C * V[n + 1][0];
			}
			else
			{
				C = EGM2008[n][m];     /* = C_n,m  */
				S = EGM2008[m - 1][n];   /* = S_n,m  */

				Fac = 0.5 * (n - m + 1) * (n - m + 2);
				ax = ax + 0.5 * (-C * V[n + 1][m + 1] - S * W[n + 1][m + 1])
					+ Fac * (C * V[n + 1][m - 1] + S * W[n + 1][m - 1]);
				ay = ay + 0.5 * (-C * W[n + 1][m + 1] + S * V[n + 1][m + 1])
					+ Fac * (-C * W[n + 1][m - 1] + S * V[n + 1][m - 1]);
				az = az + (n - m + 1) * (-C * V[n + 1][m] - S * W[n + 1][m]);
			}
		}
	}

	/* Body-fixed acceleration  */

	a_bf[0] = GM_R2 * ax;
	a_bf[1] = GM_R2 * ay;
	a_bf[2] = GM_R2 * az;

	/* Inertial acceleration */
	MatrixTranspose(3, 3, E, ET);
	MatrixMultiply(3, 3, 3, 1, ET, a_bf, Acc);

	if (HasPDev == true) 
	{
		double  GM_R3, GM_r5;
		double  dadr_ef[9], Temp[9];

		C = EGM2008[2][0];                 // only J2
		GM_r5 = GM_Earth / pow(r_sqr, 2.5);
		GM_R3 = GM_R2 / R_Earth;

		dadr_ef[0] = GM_r5 * (3.0 * r_bf[0] * r_bf[0] - r_sqr) + 0.5 * GM_R3 * C * (V[4][2] - 12.0 * V[4][0]);
		dadr_ef[1] = GM_r5 * 3.0 * r_bf[0] * r_bf[1] + 0.5 * GM_R3 * C * W[4][2];
		dadr_ef[2] = GM_r5 * 3.0 * r_bf[0] * r_bf[2] + 3.0 * GM_R3 * C * V[4][1];
		dadr_ef[3] = dadr_ef[1];
		dadr_ef[5] = GM_r5 * 3.0 * r_bf[1] * r_bf[2] + GM_R3 * 3.0 * C * W[4][1];
		dadr_ef[6] = dadr_ef[2];
		dadr_ef[7] = dadr_ef[5];
		dadr_ef[8] = GM_r5 * (3.0 * r_bf[2] * r_bf[2] - r_sqr) + GM_R3 * 12.0 * C * V[4][0];
		dadr_ef[4] = -(dadr_ef[0] + dadr_ef[8]);

		MatrixMultiply(3, 3, 3, 3, ET, dadr_ef, Temp);
		MatrixMultiply(3, 3, 3, 3, Temp, E, dadr);
	}
}

void Accel(const MJDTIME* Mjd_GPS, const double Y[6], double dY[6])
{
	int     i;
	double  Acc[3];
	double  E[9], dadr[9];                       /* transform matrix   */
	MJDTIME Mjd_TT;

	Mjd_TT.Days = Mjd_GPS->Days;
	Mjd_TT.FracDay = Mjd_GPS->FracDay - (BDT_TAI - TT_TAI) / SECPERDAY;

	ICRF_ITRF_Matrix(MJD_J2000, Mjd_GPS, 1, E);

	//AccelMain(&Mjd_TT, Y, &Y[3], E, SBInfo, Acc, 0, dadr);
	AccelHarmonic(Y, E, Acc, false, dadr);

	for (i = 0; i < 3; i++)
	{
		dY[i] = Y[3 + i];
		dY[3 + i] = Acc[i];
	}
}

void RK4Step(MJDTIME* Mjd_GPS, double step, double Y0[6])
{
	int i;
	double Y[6], dY[4][6];
	double h = step / SECPERDAY; 

	CopyArray(6, Y, Y0);
	Accel(Mjd_GPS, Y, dY[0]);

	Mjd_GPS->FracDay = Mjd_GPS->FracDay + h / 2.0;
	for (i = 0; i < 6; i++)
	{
		Y[i] = Y0[i] + dY[0][i] * step / 2.0;
	}
	Accel(Mjd_GPS, Y, dY[1]);     /* 2nd time */

	Mjd_GPS->FracDay = Mjd_GPS->FracDay + h / 2.0;
	for (i = 0; i < 6; i++)
	{
		Y[i] = Y0[i] + dY[1][i] * step / 2.0;
	}
	Accel(Mjd_GPS, Y, dY[2]);     /* 3rd time */

	Mjd_GPS->FracDay = Mjd_GPS->FracDay + h / 2.0;
	for (i = 0; i < 6; i++)
	{
		Y[i] = Y0[i] + dY[2][i] * step;
	}
	Accel(Mjd_GPS, Y, dY[3]);    /* 4th time */


	if (Mjd_GPS->FracDay >= 1.0)
	{
		Mjd_GPS->FracDay -= 1.0;
		Mjd_GPS->Days += 1;
	}

	for (i = 0; i < 6; i++)
	{
		Y0[i] += (dY[0][i] + 2.0 * dY[1][i] + 2.0 * dY[2][i] + dY[3][i]) * step / 6.0;
	}
}

void OrbitIntegToGivenTime(const GPSTIME* BegGPS, const GPSTIME* GivenTime, const double Pace,
	double Y0[])
{
	short   sign;
	double  Step;              
	double  dT;                

	MJDTIME Mjd_GPS0, Mjd_GivenTime;
	GPSTimeToMJDTime(BegGPS, &Mjd_GPS0);
	GPSTimeToMJDTime(GivenTime, &Mjd_GivenTime);

	dT = GetDifGPSTime(GivenTime, BegGPS);

	if (dT >= 0.0)
		sign = 1;       
	else
		sign = -1;

	do {
		if (fabs(dT) > Pace)
			Step = Pace * sign;   
		else
			Step = dT;

		RK4Step(&Mjd_GPS0, Step, Y0);
		dT = dT - Step;

	} while (fabs(dT) > 1E-8);
}

double GetPseudoRange(int tid, int w, double sec, double rpv[6], double tpv[6])
{

	GPSTIME   T, T0, T1;
	double dPos[3];
	double prePos[6];
	double range;
	double dt = 0, dt0 = 0;

	T0.Week = w;
	T0.SecOfWeek = sec;

	for (int j = 0; j < 3; j++)     dPos[j] = rpv[j] - tpv[j];
	range = sqrt(VectDot(3, 3, dPos, dPos));

	dt = range / C_Light;

	for (int i = 0; i < 10; i++)
	{
		T1 = T0;
		T1.SecOfWeek -= dt;
		CheckGPSTime(&T1);
		memcpy(prePos, tpv, sizeof(prePos));
		if (tid != 76)
			OrbitIntegToGivenTime(&T0, &T1, dt, prePos);
		else {
			ICRF_ITRF_GPST(MJD_J2000, &T1, 0, prePos, Anchor[0].Pos);
		}
		for (int j = 0; j < 3; j++)     dPos[j] = rpv[j] - prePos[j];

		range = sqrt(VectDot(3, 3, dPos, dPos));
		dt0 = range / C_Light;
		if (fabs(dt0 - dt) < 1E-10)
		{
			for (size_t i = 0; i < 6; i++)  tpv[i] = prePos[i];
			break;
		}

		dt = dt0;
	}

	return range;
}

#define MID_RES_PATH "../res/"
int OpenANSResFileDaily(GPSTIME* Time) {
	int RetVal;
	int i = 0, len = 0;
	COMMONTIME ct;
	char FileName[512] = { 0 };

	RetVal = 0;
	GPSTimeToCommonTime(Time, &ct);

	memset(FileName, 0, sizeof(FileName));
	sprintf(FileName, "%s%4d_%02d_%02d/DerOrb_%4d_%02d_%02d.txt",
		MID_RES_PATH, ct.Year, ct.Month, ct.Day, ct.Year, ct.Month, ct.Day);
	if (FDerOrb != NULL)		fclose(FDerOrb);
	if ((FDerOrb = fopen(FileName, "w")) == NULL)
	{
		RetVal = -1;
	}

	memset(FileName, 0, sizeof(FileName));
	sprintf(FileName, "%s%4d_%02d_%02d/DerClk_%4d_%02d_%02d.txt",
		MID_RES_PATH, ct.Year, ct.Month, ct.Day, ct.Year, ct.Month, ct.Day);
	if (FDerClk != NULL)	fclose(FDerClk);
	if ((FDerClk = fopen(FileName, "w")) == NULL)
	{
		RetVal = -1;
	}

	memset(FileName, 0, sizeof(FileName));
	sprintf(FileName, "%s%4d_%02d_%02d/MidRes_%4d_%02d_%02d.txt",
		MID_RES_PATH, ct.Year, ct.Month, ct.Day, ct.Year, ct.Month, ct.Day);
	if (FPProc != NULL)		fclose(FPProc);
	if ((FPProc = fopen(FileName, "w")) == NULL)
	{
		RetVal = -1;
	}

	memset(FileName, 0, sizeof(FileName));
	sprintf(FileName, "%s%4d_%02d_%02d/OutlierLog_%4d_%02d_%02d.txt",
		MID_RES_PATH, ct.Year, ct.Month, ct.Day, ct.Year, ct.Month, ct.Day);
	if (FPLOG != NULL)		fclose(FPLOG);
	if ((FPLOG = fopen(FileName, "w")) == NULL)
	{
		RetVal = -1;
	}

	memset(FileName, 0, sizeof(FileName));
	sprintf(FileName, "%s%4d_%02d_%02d/ANSRes_%4d_%02d_%02d.txt", MID_RES_PATH, ct.Year, ct.Month, ct.Day, ct.Year, ct.Month, ct.Day);
	if (Fres != NULL)		fclose(Fres);
	if ((Fres = fopen(FileName, "w")) == NULL)
	{
		RetVal = -1;
	}

	memset(FileName, 0, sizeof(FileName));
	sprintf(FileName, "%s%4d_%02d_%02d/MAO_%4d_%02d_%02d.txt",
		MID_RES_PATH, ct.Year, ct.Month, ct.Day, ct.Year, ct.Month, ct.Day);
	//if (FDerMAO != NULL)	fclose(FDerMAO);
	if ((FDerMAO = fopen(FileName, "wb")) == NULL)
	{
		RetVal = -1;
	}

	return RetVal;
}

int AssignEpkISLObs(GPSTIME* Time, ISLPROBS* il, SATNET* SatNet) {
	int n, id;
	double dt = 0.0;
	Point *p1, *p2;
	node_t *iter;
	ISLPROBS* isl;
	EdgeList *edge_in, *edge_out;

	graph_refresh((Graph*)SatNet);
	n = 0;

	int Rvrank, Trrank;

	if (world_rank == 0) {
		for (iter = il->next; iter != (node_t*)il;) {
			isl = (ISLPROBS*)iter;
			dt = GetDifGPSTime(&isl->RvLocTime, Time);
			iter = iter->next;
			if (dt < StepOfAutoNav) {	
				if (dt >= 0.0) {
					list_delete((list*)il, iter->prev);
					isl->PRObs = GetPseudoRange(isl->TrScid, isl->RvLocTime.Week, isl->RvLocTime.SecOfWeek, isl->rpv, isl->tpv);
					if (isl->TrScid < 0 || isl->RvScid != 0 && fabs(isl->PRObs) < 0.1) {
						free(isl);
						continue;
					}
					Rvrank = GetPointRank(isl->RvScid);
					Trrank = GetPointRank(isl->TrScid);
					p1 = (Point*)GetPoint(isl->RvScid);
					p2 = (Point*)GetPoint(isl->TrScid);
					isl->endpoints[0] = p1;
					isl->endpoints[1] = p2;
					if (Rvrank == Trrank) {
						// the same process
						
						if (Rvrank == 0) {
							// both endpoints in the main process
							if (p1 == NULL || p2 == NULL) continue;
							edge_in = (EdgeList*)malloc(sizeof(EdgeList));
							edge_in->edge = (Edge*)isl;
							edge_in->flag = 0;
							list_append((list*)&p1->in_edges, (node_t*)edge_in);
							p1->in_num++;
							edge_out = (EdgeList*)malloc(sizeof(EdgeList));
							edge_out->edge = (Edge*)isl;
							edge_out->flag = 0;
							list_append((list*)&p2->out_edges, (node_t*)edge_out);
							p2->out_num++;
							edge_in->reverse = edge_out;
							edge_out->reverse = edge_in;
							n++;
						}
						else {
							// only need to send to one other process
							MPI_Send((char*)isl, sizeof(ISLPROBS), MPI_BYTE, Rvrank, 0, MPI_COMM_WORLD);
							free(isl);
						}
					} else {
						// different processes
						if (Rvrank == 0) {
							// p1 not NULL
							if (p1 == NULL) continue;
							edge_in = (EdgeList*)malloc(sizeof(EdgeList));
							edge_in->edge = (Edge*)isl;
							edge_in->flag = 0;
							list_append((list*)&p1->in_edges, (node_t*)edge_in);
							p1->in_num++;
							edge_in->reverse = NULL;
							MPI_Send((char*)isl, sizeof(ISLPROBS), MPI_BYTE, Trrank, 0, MPI_COMM_WORLD);
							p2 = (Point*)GetBorderPoint(isl->TrScid);
							if (p2 == NULL) continue;
							isl->endpoints[1] = p2;
							n++;
						}
						else if (Trrank == 0) {
							// p2 not NULL
							if (p2 == NULL) continue;
							edge_out = (EdgeList*)malloc(sizeof(EdgeList));
							edge_out->edge = (Edge*)isl;
							edge_out->flag = 0;
							list_append((list*)&p2->out_edges, (node_t*)edge_out);
							p2->out_num++;
							edge_out->reverse = NULL;
							MPI_Send((char*)isl, sizeof(ISLPROBS), MPI_BYTE, Rvrank, 0, MPI_COMM_WORLD);
							p1 = (Point*)GetBorderPoint(isl->RvScid);
							if (p1 == NULL) continue;
							isl->endpoints[0] = p1;
							n++;
						} else {
							// send to p1 and p2 processes
							MPI_Send((char*)isl, sizeof(ISLPROBS), MPI_BYTE, Rvrank, 0, MPI_COMM_WORLD);
							MPI_Send((char*)isl, sizeof(ISLPROBS), MPI_BYTE, Trrank, 0, MPI_COMM_WORLD);
							free(isl);
						}
					}
				}
				else {
					list_erase((list*)il, iter->prev);
				}
			}
			else {
				break;
			}
			//TODO find repeat links
		}
		MPI_Send(isl, 0, MPI_BYTE, 1, 0, MPI_COMM_WORLD);
		MPI_Send(isl, 0, MPI_BYTE, 2, 0, MPI_COMM_WORLD);
	}
	else {
		MPI_Status status;
		int count;
		while (1) {
			isl = (ISLPROBS*)malloc(sizeof(ISLPROBS));
			MPI_Recv((char*)isl, sizeof(ISLPROBS), MPI_BYTE, 0, 0, MPI_COMM_WORLD, &status);
			MPI_Get_count(&status, MPI_BYTE, &count);
			if (count <= 0) break; 
			// isl->PRObs = GetPseudoRange(isl->TrScid, isl->RvLocTime.Week, isl->RvLocTime.SecOfWeek, isl->rpv, isl->tpv);
			// if (isl->TrScid < 0 || isl->RvScid != 0 && fabs(isl->PRObs) < 0.1) {
			// 	free(isl);
			// 	continue;
			// }
			p1 = (Point*)GetPoint(isl->RvScid);
			p2 = (Point*)GetPoint(isl->TrScid);
			isl->endpoints[0] = p1;
			isl->endpoints[1] = p2;
			if (p1 != NULL) {
				edge_in = (EdgeList*)malloc(sizeof(EdgeList));
				edge_in->edge = (Edge*)isl;
				edge_in->flag = 0;
				list_append((list*)&p1->in_edges, (node_t*)edge_in);
				p1->in_num++;
				edge_in->reverse = NULL;
				if (p2 != NULL) {
					edge_out = (EdgeList*)malloc(sizeof(EdgeList));
					edge_out->edge = (Edge*)isl;
					edge_out->flag = 0;
					list_append((list*)&p2->out_edges, (node_t*)edge_out);
					p2->out_num++;
					edge_in->reverse = edge_out;
					edge_out->reverse = edge_in;
				} else {
					p2 = (Point*)GetBorderPoint(isl->TrScid);
					if (p2 == NULL) continue;
					isl->endpoints[1] = p2;
				}
			} else {
				edge_out = (EdgeList*)malloc(sizeof(EdgeList));
				edge_out->edge = (Edge*)isl;
				edge_out->flag = 0;
				list_append((list*)&p2->out_edges, (node_t*)edge_out);
				p2->out_num++;
				edge_out->reverse = NULL;
				p1 = (Point*)GetBorderPoint(isl->RvScid);
				if (p1 == NULL) continue;
				isl->endpoints[0] = p1;
			}			
			n++;
		}
		free(isl);
		
	}
	
	printf("读取数据成功，共读取%d条星间/星地单向链路数据\n", n);
	return 0;
}

void VarEquation(const MJDTIME* Mjd_GPS, const double Phi[],
	double dPhi[])
{
	int       i, j;
	double    Pos[3], Vel[3], Acc[3];      /* Position, velocity of spacecraft */
	double    dfdy[36];                    /*  df/dy  ( y = {r,v } )  */
	double    E[9], dadr[9];               /* transform matrix   */
	MJDTIME   Mjd_TT;

	Mjd_TT.Days = Mjd_GPS->Days;
	Mjd_TT.FracDay = Mjd_GPS->FracDay - (BDT_TAI - TT_TAI) / SECPERDAY;

	ICRF_ITRF_Matrix(MJD_J2000, Mjd_GPS, 1, E);

	for (i = 0; i < 3; i++)    /* initialize spacecraft position and velocity */
	{
		Pos[i] = Phi[i];
		Vel[i] = Phi[3 + i];
	}

	AccelHarmonic(Pos, E, Acc, true, dadr);

	/*  set df/dy matrix
		|  0      I   |
		|  da/dr  0   |
	*/

	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			*(dfdy + i * 6 + j) = 0.0;     /*  dv/dr  */
			if (i == j)               /*  dv/dv  */
			{
				*(dfdy + i * 6 + 3 + j) = 1.0;
			}
			else
			{
				*(dfdy + i * 6 + 3 + j) = 0.0;
			}

			*(dfdy + (3 + i) * 6 + j) = *(dadr + i * 3 + j);
			*(dfdy + (3 + i) * 6 + 3 + j) = 0.0;
		}
	}

	MatrixMultiply(6, 6, 6, DIM, dfdy, Phi + 6, dPhi + 6);

	for (i = 0; i < 3; i++)
	{
		*(dPhi + i) = *(Phi + 3 + i);     /* velocity set */
		*(dPhi + 3 + i) = Acc[i];         /* acceleration set */
	}
}

void RKF4OrbitSTM(MJDTIME* Mjd_GPS, double step, double Y0[])
{
	int i;
	double Y[Dim_STM], dY[5][Dim_STM];
	double h = step / SECPERDAY;

	memset(Y, 0, Dim_STM * sizeof(double));
	memset(dY, 0, 5 * Dim_STM * sizeof(double));

	CopyArray(Dim_STM, Y, Y0);
	VarEquation(Mjd_GPS, Y, dY[0]);     /* first time */

	Mjd_GPS->FracDay += h / 4.0;
	for (i = 0; i < Dim_STM; i++)
	{
		Y[i] = Y0[i] + dY[0][i] * step / 4.0;
	}
	VarEquation(Mjd_GPS, Y, dY[1]);     /* 2nd time */

	Mjd_GPS->FracDay += h / 8.0;
	for (i = 0; i < Dim_STM; i++)
	{
		Y[i] = Y0[i] + (dY[0][i] * 3 / 32 + dY[1][i] * 9 / 32.0) * step;
	}
	VarEquation(Mjd_GPS, Y, dY[2]);     /* 3rd time */

	Mjd_GPS->FracDay += h * 57.0 / 104.0;
	for (i = 0; i < Dim_STM; i++)
	{
		Y[i] = Y0[i] + (dY[0][i] * 1932 / 2197.0 - dY[1][i] * 7200 / 2197.0
			+ dY[2][i] * 7296 / 2197.0) * step;
	}
	VarEquation(Mjd_GPS, Y, dY[3]);     /* 4th time */

	Mjd_GPS->FracDay += h / 13.0;
	for (i = 0; i < Dim_STM; i++)
	{
		Y[i] = Y0[i] + (dY[0][i] * 439 / 216.0 - dY[1][i] * 8.0
			+ dY[2][i] * 3680 / 513.0 - dY[3][i] * 845 / 4104.0) * step;
	}
	VarEquation(Mjd_GPS, Y, dY[4]);     /* 5th time */

	if (Mjd_GPS->FracDay >= 1.0)
	{
		Mjd_GPS->FracDay -= 1.0;
		Mjd_GPS->Days += 1;
	}

	for (i = 0; i < Dim_STM; i++)
	{
		Y0[i] += (dY[0][i] * 25 / 216 + dY[2][i] * 1408 / 2565
			+ dY[3][i] * 2197 / 4104 - dY[4][i] / 5) * step;
	}
}

void OrbitSTMIntegToGivenTime(const GPSTIME* BegGPS, const GPSTIME* GivenTime, const double Pace,
	double Y0[])
{
	short   sign;
	double  Step;
	double  dT;

	MJDTIME Mjd_GPS0, Mjd_GivenTime;
	GPSTimeToMJDTime(BegGPS, &Mjd_GPS0);
	GPSTimeToMJDTime(GivenTime, &Mjd_GivenTime);

	dT = GetDifGPSTime(GivenTime, BegGPS);

	if (dT >= 0.0)
	{
		sign = 1;
	}
	else
	{
		sign = -1;
	}

	do {
		if (fabs(dT) > Pace)
		{
			Step = Pace * sign;
		}
		else
		{
			Step = dT;
		}

		RKF4OrbitSTM(&Mjd_GPS0, Step, Y0);
		dT = dT - Step;

	} while (fabs(dT) > 1E-8);
}

void RKF4OrbitSTMForInterp(const GPSTIME* GT, const double step, double Y0[Dim_STM],
	SCSTM Interp[2])
{
	int i;
	double Y[Dim_STM], dY[5][Dim_STM];
	double h = step / SECPERDAY;
	MJDTIME Mjd_GPS;

	GPSTimeToMJDTime(GT, &Mjd_GPS);

	CopyArray(Dim_STM, Y, Y0);
	VarEquation(&Mjd_GPS, Y, dY[0]);     /* first time */ 

	Interp[0].Time.Week = GT->Week;
	Interp[0].Time.SecOfWeek = GT->SecOfWeek;
	CopyArray(Dim_STM, Interp[0].Phi, Y0);
	CopyArray(Dim_STM, Interp[0].dPhi, dY[0]);

	Mjd_GPS.FracDay += h / 4.0;
	for (i = 0; i < Dim_STM; i++)
	{
		Y[i] = Y0[i] + dY[0][i] * step / 4.0;
	}
	VarEquation(&Mjd_GPS, Y, dY[1]);     /* 2nd time */

	Mjd_GPS.FracDay += h / 8.0;
	for (i = 0; i < Dim_STM; i++)
	{
		Y[i] = Y0[i] + (dY[0][i] * 3 / 32 + dY[1][i] * 9 / 32.0) * step;
	}
	VarEquation(&Mjd_GPS, Y, dY[2]);     /* 3rd time */

	Mjd_GPS.FracDay += h * 57.0 / 104.0;
	for (i = 0; i < Dim_STM; i++)
	{
		Y[i] = Y0[i] + (dY[0][i] * 1932 / 2197.0 - dY[1][i] * 7200 / 2197.0
			+ dY[2][i] * 7296 / 2197.0) * step;
	}
	VarEquation(&Mjd_GPS, Y, dY[3]);     /* 4th time */

	Mjd_GPS.FracDay += h / 13.0;
	for (i = 0; i < Dim_STM; i++)
	{
		Y[i] = Y0[i] + (dY[0][i] * 439 / 216.0 - dY[1][i] * 8.0
			+ dY[2][i] * 3680 / 513.0 - dY[3][i] * 845 / 4104.0) * step;
	}
	VarEquation(&Mjd_GPS, Y, dY[4]);     /* 5th time */

	Interp[1].Time.Week = GT->Week;
	Interp[1].Time.SecOfWeek = GT->SecOfWeek + step;

	for (i = 0; i < Dim_STM; i++)
	{
		Interp[1].Phi[i] = Y0[i] + (dY[0][i] * 25 / 216 + dY[2][i] * 1408 / 2565
			+ dY[3][i] * 2197 / 4104 - dY[4][i] / 5) * step;
	}

	VarEquation(&Mjd_GPS, Interp[1].Phi, Interp[1].dPhi);
}

void TimeUpdate(const GPSTIME* Time, SATNET* SatNet) {
	int	i;
	SATINFO* sat;
	CONSTSTATE* AllSatCov;
	double dt, Stm[Dim_STM], q;
	SpraseMatrix *P_Clk, *P_Orb;
	SpraseMatrix *Q_Clk, *Q_Orb;
	SpraseMatrix *Cov_Clk, *Cov_Orb;
	double P1_Clk[4];
	double Q1_Clk[4];
	double Q1_Orb[DIM * DIM];
	double Pace = 60.0, PredPace = 90.0;
	MJDTIME Mjd;
	double h[3] = { 1.43E-11, 5.0E-14, 7.61E-17 };

	int num = SatNum;
	AllSatCov = &SatNet->AllSatCov;
	P_Clk = (SpraseMatrix*)malloc(sizeof(SpraseMatrix));
	Q_Clk = (SpraseMatrix*)malloc(sizeof(SpraseMatrix));
	Cov_Clk = (SpraseMatrix*)malloc(sizeof(SpraseMatrix));
	P_Orb = (SpraseMatrix*)malloc(sizeof(SpraseMatrix));
	Q_Orb = (SpraseMatrix*)malloc(sizeof(SpraseMatrix));
	Cov_Orb = (SpraseMatrix*)malloc(sizeof(SpraseMatrix));

	CreateSpraseMatrix(P_Clk, num, num, 2, 2);
	CreateSpraseMatrix(Q_Clk, num, num, 2, 2);
	CreateSpraseMatrix(Cov_Clk, num, num, 2, 2);
	CreateSpraseMatrix(P_Orb, num, num, DIM, DIM);
	CreateSpraseMatrix(Q_Orb, num, num, DIM, DIM);
	CreateSpraseMatrix(Cov_Orb, num, num, DIM, DIM);

	sat = (SATINFO*)SatNet->points;
	for (i = 0; i < num; i++) {
		sat = (SATINFO*)sat->next;

		dt = GetDifGPSTime(Time, &(sat->TOE));
		q = h[0] * h[0] / dt + h[1] * h[1] + h[2] * h[2] * dt;
		/*
		Q1_Clk = [dt^3/3*q    dt^2/2*q]
				 [dt^2/2*q    dt*q ]
		*/
		Q1_Clk[0] = pow(dt, 3.0) * q / 3.0;
		Q1_Clk[1] = pow(dt, 2.0) * q / 2.0;
		Q1_Clk[2] = pow(dt, 2.0) * q / 2.0;
		Q1_Clk[3] = dt * q;
		/*
		 P1_Clk= [1  dt
				   0  1 ]
		*/
		P1_Clk[0] = P1_Clk[3] = 1.0;
		P1_Clk[1] = dt;
		P1_Clk[2] = 0.0;

		CopySubSpraseMatrix(P_Clk, i, i, 2, 2, P1_Clk);
		CopySubSpraseMatrix(Q_Clk, i, i, 2, 2, Q1_Clk);


		sat->Clk[0] += sat->Clk[1] * dt;
		memcpy(AllSatCov->Clk + 2 * i, sat->Clk, sizeof(double) * 2);


		memcpy(Stm, sat->X, sizeof(double) * 6);
		InitStateTranMatrix(6, DIM, Stm);
		OrbitSTMIntegToGivenTime(&(sat->TOE), Time, Pace, Stm);  
		CopyArray(6, sat->X, Stm);
		
		CompStateNoiseCov(dt, sat->Valid, Q1_Orb);

		CopySubSpraseMatrix(P_Orb, i, i, DIM, DIM, Stm + 6);
		CopySubSpraseMatrix(Q_Orb, i, i, DIM, DIM, Q1_Orb);

		InitStateTranMatrix(6, DIM, Stm);
		RKF4OrbitSTMForInterp(Time, PredPace, Stm, sat->STM);

		memcpy(&(sat->TOE), Time, sizeof(GPSTIME));
		sat->GapTime = sat->GapTime + dt;

		if (sat->GapTime > 2.0 * SECPERHOUR && sat->Valid != BREAKING) {
																				
			sat->Valid = BREAKING;
			sat->Health = 2;

			GPSTimeToMJDTime(Time, &Mjd);
			fprintf(FPLOG, "OBS breaking interval: %2d %12.5f %6d %10.1f  dt %10.1f\n",
				sat->id, Mjd.Days + Mjd.FracDay,
				Time->Week, Time->SecOfWeek,
				sat->GapTime);

		}
	}

	SpraseMatrixMultiply_MPMT(P_Clk, AllSatCov->ClkCov, Cov_Clk);  // M*P*MT+Q
	SpraseMatrixAddition(Cov_Clk, Q_Clk, AllSatCov->ClkCov);
	SpraseMatrixMultiply_MPMT(P_Orb, AllSatCov->OrbCov, Cov_Orb);
	SpraseMatrixAddition(Cov_Orb, Q_Orb, AllSatCov->OrbCov);

	sat = (SATINFO*)SatNet->points;
	for (i = 0; i < num; i++) {
		sat = (SATINFO*)sat->next;
		GetSubSpraseMatrix(AllSatCov->ClkCov, i, i, 2, 2, sat->CovC);
		GetSubSpraseMatrix(AllSatCov->OrbCov, i, i, DIM, DIM, sat->CovX);
	}

	DestroySpraseMatrix(P_Clk);
	DestroySpraseMatrix(P_Orb);
	DestroySpraseMatrix(Q_Clk);
	DestroySpraseMatrix(Q_Orb);
	DestroySpraseMatrix(Cov_Clk);
	DestroySpraseMatrix(Cov_Orb);


	// update all node
	int n, type, rank;
	SATINFO* p, *p_old;
	for (rank = 0; rank < world_size; rank++) {
		if (world_rank == rank) {
			p = (SATINFO*)SatNet->points;
			n = SatNet->points_num;
		}
		MPI_Bcast(&n, 1, MPI_INT, rank, MPI_COMM_WORLD);
		for (i = 0; i < n; i++) {
			if (world_rank == rank) {
				p = (SATINFO*)p->next;
				type = p->type;
			} 
			MPI_Bcast(&type, 1, MPI_INT, rank, MPI_COMM_WORLD);
			if (type == 1) {
				if (world_rank != rank) p = (SATINFO*)malloc(sizeof(SATINFO));
				MPI_Bcast(p, sizeof(SATINFO), MPI_BYTE, rank, MPI_COMM_WORLD);
				if (world_rank != rank) {
					p_old = (SATINFO*)GetBorderPoint(p->id);
					memcpy(p_old->Clk, p->Clk, sizeof(double) * 2);
					memcpy(p_old->X, p->X, sizeof(double) * DIM);
					memcpy(p_old->CovC, p->CovC, sizeof(double) * 4);
					memcpy(p_old->CovX, p->CovX, sizeof(double) * DIM * DIM);
					p_old->STM[0] = p->STM[0];
					p_old->STM[1] = p->STM[1];
					p_old->Valid = p->Valid;
					p_old->Health = p->Health;
					p_old->GapTime = p->GapTime;
					free(p);
				}

			}
		}
	}

	printf("完成时间更新\n");
}

void InitStateTranMatrix(int row, int col, double STM[])
{
	int j, k;

	for (j = 0; j < row; j++)
	{
		for (k = 0; k < col; k++)
		{
			if (j == k)
			{
				*(STM + 6 + j * DIM + k) = 1.0;
			}
			else
			{
				*(STM + 6 + j * DIM + k) = 0.0;
			}
		}
	}
}

double PSDofBreak = 1E-5;
double PSDofOK = 1E-5;
void CompStateNoiseCov(const double Step, const ANSTATEID Valid, double Q[])
{
	int i;
	double q, q2;

	switch (Valid)
	{
	case MANEUVER:
		q = 1.0E-2;
		break;
	case BREAKING:
		q = PSDofBreak;
		break;
	case ANSOK:
		q = PSDofOK;
		break;
	default:
		q = 1.0;
	}

	q2 = pow(q, 2.0);

	for (i = 0; i < DIM * DIM; i++)      Q[i] = 0.0;

	for (i = 0; i < 3; i++)
	{
		Q[i * DIM + i] = pow(Step, 3.0) / 3.0 * q2;
		Q[i * DIM + i + 3] = pow(Step, 2.0) / 2.0 * q2;
		Q[(i + 3) * DIM + i] = pow(Step, 2.0) / 2.0 * q2;
		Q[(i + 3) * (DIM + 1)] = Step * q2;
	}
}

int GenDerPrObs(SATNET* SatNet) {
	int n, anc_id;
	DEROBS* obs;
	ISLPROBS* isl, * reisl;
	DEROBS* EpkDerObs = (DEROBS*)SatNet->edges;
	Point* p1, * p2;
	
	//SatNet->ObsNum = 0;
	n = 0;

	Point* node;
	EdgeList *e_in, *e_in_t, *e_out, *e_out_t, *obs_el, *el;
	ISLPROBS *in_edge, *out_edge;
	int flag, req_flag;
	int req_rank;
	MPI_Status status;
	int count;
	int req_buff[2];

	// if (world_rank == 0){
	// 	p1 = (Point*)SatNet->points->next;
	// 	for (e_in_t = (EdgeList*)p1->in_edges.next; e_in_t != &p1->in_edges; e_in_t = (EdgeList*)e_in_t->next){
	// 		in_edge = (ISLPROBS*)e_in_t->edge;
	// 		printf("process %d edge %d %d\n", world_rank, in_edge->RvScid, in_edge->TrScid);
	// 	}
	// }
	for (node = (Point*)SatNet->points->next; node != SatNet->points; node = (Point*)node->next) {
		
		for (e_in_t = (EdgeList*)node->in_edges.next; e_in_t != &node->in_edges; ) {
			e_in = e_in_t;
			e_in_t = (EdgeList*)e_in_t->next;
			if (e_in->flag == 1){
				continue;
			} 

			in_edge = (ISLPROBS*)e_in->edge;
			//printf("process %d in edge %d %d\n", world_rank, in_edge->RvScid, in_edge->TrScid);
			flag = 0;
			for (e_out_t = (EdgeList*)node->out_edges.next; e_out_t != &node->out_edges; ) {
				e_out = e_out_t;
				e_out_t = (EdgeList*)e_out_t->next;
				if (e_out->flag == 1) {
					continue;
				}
				
				out_edge = (ISLPROBS*)e_out->edge;
				//printf("process %d out edge %d %d\n", world_rank, out_edge->RvScid, out_edge->TrScid);
				if (out_edge->RvScid == in_edge->TrScid) {
					if (e_in->reverse != NULL && e_out->reverse != NULL) {
						if (e_in->reverse->flag == 1 && e_out->flag == 2) {
							free(e_in->edge);
							e_in->edge = e_in->reverse->edge;
							e_in->flag = 1;
							free(e_out->edge);
							list_erase((list*)&node->out_edges, (node_t*)e_out);
							node->out_num--;
							flag = 1;
							break;
						}
						if (e_in->flag == 2 && e_out->reverse->flag == 1) {
							free(e_in->edge);
							list_erase((list*)&node->in_edges, (node_t*)e_in);
							node->in_num--;
							free(e_out->edge);
							e_out->edge = e_out->reverse->edge;
							e_out->flag = 1;
							flag = 1;
							break;
						}

						if (GetDifGPSTime(&((ISLPROBS*)in_edge)->RvLocTime, &((ISLPROBS*)out_edge)->RvLocTime) < 0) {
							isl = (ISLPROBS*)in_edge;
							reisl = (ISLPROBS*)out_edge;
							obs_el = (EdgeList*)e_in;
							e_in->flag = 1;
							e_out->reverse->flag = 2;
							list_erase((list*)&node->out_edges, (node_t*)e_out);
							node->out_num--;
						}
						else {
							reisl = (ISLPROBS*)in_edge;
							isl = (ISLPROBS*)out_edge;
							obs_el = (EdgeList*)e_out;
							e_out->flag = 1;
							e_in->reverse->flag = 2;
							list_erase((list*)&node->in_edges, (node_t*)e_in);
							node->in_num--;
						}
						
					} 
					else {
						if (GetDifGPSTime(&((ISLPROBS*)in_edge)->RvLocTime, &((ISLPROBS*)out_edge)->RvLocTime) < 0) {
							isl = in_edge;
							reisl = out_edge;
							obs_el = e_in;
							e_in->flag = 1;
							list_erase((list*)&node->out_edges, (node_t*)e_out);
							node->out_num--;
						}
						else {
							reisl = in_edge;
							isl = out_edge;
							obs_el = e_out;
							e_out->flag = 1;
							list_erase((list*)&node->in_edges, (node_t*)e_in);
							node->in_num--;
						}
					}
					
					obs = (DEROBS*)malloc(sizeof(DEROBS));
					obs->Scid1 = isl->RvScid;
					obs->Scid2 = isl->TrScid;
					obs->RvAnt = isl->RvAnt;
					obs->TrAnt = isl->TrAnt;
					obs->T1 = isl->RvLocTime;
					obs->T2 = reisl->RvLocTime;
					obs->DerCObs = isl->PRObs - reisl->PRObs;
					obs->DerANObs = isl->PRObs + reisl->PRObs;
					obs->ClkBlunder[0] = 0;
					obs->ClkBlunder[1] = 0;
					obs->OrbBlunder[0] = 0;
					obs->OrbBlunder[1] = 0;
					obs->Valid = 0;
					obs->endpoints[0] = isl->endpoints[0];
					obs->endpoints[1] = isl->endpoints[1];
					obs_el->edge = (Edge*)obs;
					if (GenDerObsPredict(obs, isl, reisl) == 0) {
						free(obs);
						continue;
					}
					flag = 1;
					// list_append((list*)EpkDerObs, (node_t*)obs);
					
					n++;
					break;
				}
			}
			if (flag == 0) {
				e_in->reverse->flag = 2;
				if (world_rank==0)printf("here in %d %d\n",e_in->edge->endpoints[0]->index,e_in->edge->endpoints[1]->index);
				free(e_in->edge);
				e_in = (EdgeList*)e_in->next;
				list_erase((list*)&node->in_edges, ((node_t*)e_in)->prev);
				node->in_num--;
			}
		}

		for (e_out = (EdgeList*)node->out_edges.next; e_out != &node->out_edges;) {
			if (e_out->flag == 0) {
				e_out->reverse->flag = 2;
				free(e_out->edge);
				e_out = (EdgeList*)e_out->next;
				list_erase((list*)&node->out_edges,(node_t*)e_out->prev);
				printf("here out\n");
				node->out_num--;
			}
			else e_out = (EdgeList*)e_out->next;
		}
	}
	//SatNet->ObsNum = n;
	printf("生成导出观测值完成，共生成%d个导出观测值\n", n);
	return true;
}

int Hermite3ForSTM(const SCSTM S[2], const GPSTIME* GT, double STM[Dim_STM])
{
	int i;
	double delta, step;   /* delta is coefficent, step is time space */
	double delta2, delta3;
	double d[4];  /* d is coefficent for pos and vel interpolation */

	step = GetDifGPSTime(&S[1].Time, &S[0].Time);

	if (step < 0.0)
	{
		return(0);
	}

	/*  0<delta<1   */
	delta = GetDifGPSTime(GT, &S[0].Time) / step;
	delta2 = delta * delta;
	delta3 = delta * delta2;

	if (delta > -0.1 && delta < 1.0)
	{
		d[0] = 2 * delta3 - 3 * delta2 + 1.0;
		d[1] = delta * (delta2 - 2.0 * delta + 1.0);
		d[2] = delta2 * (3.0 - 2.0 * delta);
		d[3] = delta2 * (delta - 1.0);

		for (i = 0; i < Dim_STM; i++)
		{
			STM[i] = d[0] * S[0].Phi[i] + d[1] * step * S[0].dPhi[i]
				+ d[2] * S[1].Phi[i] + d[3] * step * S[1].dPhi[i];
		}
	}
	else
	{
		return (0);
	}

	return (1);
}

double GetRelCorr(const double Pos[])
{
	double Corr;

	Corr = -2.0 * VectDot(3, 3, Pos, Pos + 3) / C_Light;

	return Corr;
}

int GenDerObsPredict(DEROBS* derobs, ISLPROBS* isl1, ISLPROBS* isl2) {
	int		j;
	double	dt[4], RelCorr[4], dT_Prn1[2], dT_Prn2[2];
	double  trop[2], Tgd[4];
	double	dPos[3], Range1, Range2;
	GPSTIME GT;
	double  m_anchor_X_I[6];
	SATINFO* sat1, * sat2;
	ANCHSTN* anc;

	if (derobs->endpoints[0]->type == 1 && derobs->endpoints[1]->type == 1) {
		sat1 = (SATINFO*)derobs->endpoints[0];
		sat2 = (SATINFO*)derobs->endpoints[1];
		if (sat1->Valid <= NOINIT || sat2->Valid <= NOINIT)   return 0;

		dt[0] = GetDifGPSTime(&derobs->T1, &sat1->TOE);
		dT_Prn1[0] = sat1->Clk[0] + sat1->Clk[1] * dt[0];
		GT.Week = derobs->T1.Week;
		GT.SecOfWeek = derobs->T1.SecOfWeek - dT_Prn1[0];
		if (Hermite3ForSTM(sat1->STM, &GT, derobs->P1RvState) == 0) 	return 0;

		RelCorr[0] = 0;

		if (sat2->Health != 0) {
			dt[2] = GetDifGPSTime(&derobs->T2, &sat2->TOE);
			dT_Prn2[0] = sat2->Clk[0] + sat2->Clk[1] * dt[2];
			GT.Week = derobs->T2.Week;
			GT.SecOfWeek = derobs->T2.SecOfWeek - dT_Prn2[0];
			if (Hermite3ForSTM(sat2->STM, &GT, derobs->P2RvState) == 0)   return 0;
			RelCorr[2] = 0;

			dt[3] = GetDifGPSTime(&derobs->T1, &sat2->TOE);
			dT_Prn2[1] = sat2->Clk[0] + sat2->Clk[1] * dt[3];
			GT.Week = derobs->T1.Week;
			GT.SecOfWeek = derobs->T1.SecOfWeek - isl1->PRObs / C_Light - dT_Prn2[1];
			if (Hermite3ForSTM(sat2->STM, &GT, derobs->P2TrState) == 0)   return 0; 
			//RelCorr[3] = GetRelCorr( EpkDerObs.DerObsList[i].P2TrState );
			//Tgd[3] = SatAtod[id2].SatInfo.Tgd[EpkDerObs.DerObsList[i].TrAnt-1];
			RelCorr[3] = 0;

			dt[1] = GetDifGPSTime(&derobs->T2, &sat1->TOE);
			dT_Prn1[1] = sat1->Clk[0] + sat1->Clk[1] * dt[1];
			GT.Week = derobs->T2.Week;
			GT.SecOfWeek = derobs->T2.SecOfWeek - isl2->PRObs / C_Light - dT_Prn1[1];
			if (Hermite3ForSTM(sat1->STM, &GT, derobs->P1TrState) == 0)   return 0;
			//RelCorr[1] = GetRelCorr( EpkDerObs->DerObsList[i].P1TrState );
			//Tgd[1] = SatAtod[id1].SatInfo.Tgd[EpkDerObs->DerObsList[i].RvAnt-1];
			RelCorr[1] = 0;

			for (j = 0; j < 3; j++)     dPos[j] = derobs->P1RvState[j] - derobs->P2TrState[j];
			Range1 = sqrt(VectDot(3, 3, dPos, dPos));

			for (j = 0; j < 3; j++)     dPos[j] = derobs->P2RvState[j] - derobs->P1TrState[j];
			Range2 = sqrt(VectDot(3, 3, dPos, dPos));

			derobs->CConCorr = RelCorr[0] + RelCorr[1] - RelCorr[2] - RelCorr[3]
				- sat1->Tgd[0] + sat2->Tgd[0];
			derobs->ObsQua.dt[0] = dt[0] + dt[1];
			derobs->ObsQua.dt[1] = dt[2] + dt[3];
			derobs->CCObs = Range1 - Range2;
			derobs->ObsQua.ApriClkResid = derobs->DerCObs - derobs->CCObs - derobs->CConCorr
				- (dT_Prn1[0] + dT_Prn1[1] - dT_Prn2[0] - dT_Prn2[1]) * C_Light;

			derobs->CANObs = Range1 + Range2;
			derobs->AConCorr = RelCorr[0] - RelCorr[1] + RelCorr[2] - RelCorr[3]  
				+ (dT_Prn1[0] - dT_Prn1[1] + dT_Prn2[0] - dT_Prn2[1]) * C_Light
				+ 2 * (isl1->Corr[4] + isl2->Corr[4])
				+ (sat1->Tgd[1] + sat2->Tgd[1]) * 2;

			derobs->ObsQua.ApriOrbResid = derobs->DerANObs - derobs->CANObs - derobs->AConCorr;
			derobs->Valid = 0;

		}
	}

	else if (derobs->endpoints[0]->type == 1 && derobs->endpoints[1]->type == 0) {
		sat1 = (SATINFO*)derobs->endpoints[0];
		anc = (ANCHSTN*)derobs->endpoints[1];

		if (sat1->Valid <= NOINIT)  return 0;

		dt[0] = GetDifGPSTime(&derobs->T1, &sat1->TOE);
		dT_Prn1[0] = sat1->Clk[0] + sat1->Clk[1] * dt[0];
		GT.Week = derobs->T1.Week;
		GT.SecOfWeek = derobs->T1.SecOfWeek - dT_Prn1[0];
		if (Hermite3ForSTM(sat1->STM, &GT, derobs->P1RvState) == 0) 	return 0;
		RelCorr[0] = GetRelCorr(derobs->P1RvState);
		Tgd[0] = sat1->Tgd[derobs->RvAnt - 1];

		dt[1] = GetDifGPSTime(&derobs->T2, &sat1->TOE);
		dT_Prn1[1] = sat1->Clk[0] + sat1->Clk[1] * dt[1];
		GT.Week = derobs->T2.Week;
		GT.SecOfWeek = derobs->T2.SecOfWeek - isl2->PRObs / C_Light - dT_Prn1[1];
		if (Hermite3ForSTM(sat1->STM, &GT, derobs->P1TrState) == 0)   return 0;
		RelCorr[1] = GetRelCorr(derobs->P1TrState);
		Tgd[1] = sat1->Tgd[derobs->RvAnt - 1];

		dT_Prn2[1] = anc->clk;
		GT.Week = derobs->T1.Week;
		GT.SecOfWeek = derobs->T1.SecOfWeek - isl1->PRObs / C_Light - dT_Prn2[1];
		ICRF_ITRF_GPST(MJD_J2000, &GT, 0, m_anchor_X_I, anc->Pos);

		//trop[0] = hopfield(EpkDerObs->DerObsList[i].P1RvState,m_anchor_X_I);
		trop[0] = 0;
		Tgd[3] = anc->Tgd[0];

		for (j = 0; j < 3; j++)
		{
			dPos[j] = derobs->P1RvState[j] - m_anchor_X_I[j];
			derobs->P2TrState[j] = m_anchor_X_I[j];
		}
		Range1 = sqrt(VectDot(3, 3, dPos, dPos));

		dT_Prn2[0] = anc->clk;
		GT.Week = derobs->T2.Week;
		GT.SecOfWeek = derobs->T2.SecOfWeek - dT_Prn2[0];
		ICRF_ITRF_GPST(MJD_J2000, &GT, 0, m_anchor_X_I, anc->Pos);
		//trop[1] = hopfield(EpkDerObs->DerObsList[i].P1TrState,m_anchor_X_I);
		trop[1] = 0;
		Tgd[2] = anc->Tgd[0];
		for (j = 0; j < 3; j++)
		{
			dPos[j] = m_anchor_X_I[j] - derobs->P1TrState[j];
			derobs->P2RvState[j] = m_anchor_X_I[j];
		}
		Range2 = sqrt(VectDot(3, 3, dPos, dPos));
	
		derobs->CConCorr = RelCorr[0] + RelCorr[1] - sat1->Tgd[0] + anc->Tgd[0];
		derobs->ObsQua.dt[0] = dt[0] + dt[1];
		derobs->ObsQua.dt[1] = 0.0;
		derobs->CCObs = Range1 - Range2;
		derobs->ObsQua.ApriClkResid = derobs->DerCObs - derobs->CCObs - derobs->CConCorr - (dT_Prn1[0] + dT_Prn1[1]) * C_Light;

		derobs->CANObs = Range1 + Range2;
		derobs->AConCorr = RelCorr[0] - RelCorr[1] + (dT_Prn1[0] - dT_Prn1[1]) * C_Light
			+ 2 * (isl1->Corr[0] + isl1->Corr[1] + isl1->Corr[3] + isl1->Corr[4] + isl2->Corr[4])
			+ (sat1->Tgd[1] + anc->Tgd[1]) * 2;
		derobs->ObsQua.ApriOrbResid = derobs->DerANObs - derobs->CANObs - derobs->AConCorr;
		derobs->Valid = 0;

	}

	else if (derobs->endpoints[0]->type == 0 && derobs->endpoints[1]->type == 1) {
		anc = (ANCHSTN*)derobs->endpoints[0];
		sat2 = (SATINFO*)derobs->endpoints[1];

		if (sat2->Valid <= NOINIT)  return 0;

		dt[2] = GetDifGPSTime(&derobs->T2, &sat2->TOE);
		dT_Prn2[0] = sat2->Clk[0] + sat2->Clk[1] * dt[2];
		GT.Week = derobs->T2.Week;
		GT.SecOfWeek = derobs->T2.SecOfWeek - dT_Prn2[0];
		if (Hermite3ForSTM(sat2->STM, &GT, derobs->P2RvState) == 0)   return 0;
		RelCorr[2] = GetRelCorr(derobs->P2RvState);
		Tgd[2] = sat2->Tgd[derobs->TrAnt - 1];

		dt[3] = GetDifGPSTime(&derobs->T1, &sat2->TOE);
		dT_Prn2[1] = sat2->Clk[0] + sat2->Clk[1] * dt[3];
		GT.Week = derobs->T1.Week;
		GT.SecOfWeek = derobs->T1.SecOfWeek - isl1->PRObs / C_Light - dT_Prn2[1];
		if (Hermite3ForSTM(sat2->STM, &GT, derobs->P2TrState) == 0)  return 0;
		RelCorr[3] = GetRelCorr(derobs->P2TrState);
		Tgd[3] = sat2->Tgd[derobs->TrAnt - 1];

		dT_Prn1[0] = anc->clk;
		GT.Week = derobs->T1.Week;
		GT.SecOfWeek = derobs->T1.SecOfWeek - dT_Prn1[0];
		ICRF_ITRF_GPST(MJD_J2000, &GT, 0, m_anchor_X_I, anc->Pos);

		trop[0] = 0;
		Tgd[0] = anc->Tgd[0];

		for (j = 0; j < 3; j++)
		{
			dPos[j] = m_anchor_X_I[j] - derobs->P2TrState[j];
			derobs->P1RvState[j] = m_anchor_X_I[j];
		}

		Range1 = sqrt(VectDot(3, 3, dPos, dPos));

		dT_Prn1[1] = anc->clk;
		GT.Week = derobs->T2.Week;
		GT.SecOfWeek = derobs->T2.SecOfWeek - isl2->PRObs / C_Light - dT_Prn1[1];
		ICRF_ITRF_GPST(MJD_J2000, &GT, 0, m_anchor_X_I, anc->Pos);

		trop[1] = 0;
		Tgd[1] = anc->Tgd[0];
		for (j = 0; j < 3; j++)
		{
			dPos[j] = derobs->P2RvState[j] - m_anchor_X_I[j];
			derobs->P1TrState[j] = m_anchor_X_I[j];

		}
		Range2 = sqrt(VectDot(3, 3, dPos, dPos));

		derobs->CConCorr = (-1.0) * (RelCorr[2] + RelCorr[3]) - anc->Tgd[0] + sat2->Tgd[0];
		derobs->ObsQua.dt[0] = 0.0;
		derobs->ObsQua.dt[1] = dt[2] + dt[3];
		derobs->CCObs = Range1 - Range2;
		derobs->ObsQua.ApriClkResid = derobs->DerCObs - derobs->CCObs - derobs->CConCorr + (dT_Prn2[0] + dT_Prn2[1]) * C_Light;

		derobs->CANObs = Range1 + Range2;
		derobs->AConCorr = RelCorr[2] - RelCorr[3] + (dT_Prn2[0] - dT_Prn2[1]) * C_Light
			+ 2 * (isl1->Corr[0] + isl1->Corr[1] + isl1->Corr[3] + isl1->Corr[4] + isl2->Corr[4])
			+ (sat2->Tgd[1] + anc->Tgd[1]) * 2;
		derobs->ObsQua.ApriOrbResid = derobs->DerANObs - derobs->CANObs - derobs->AConCorr;
		derobs->Valid = 0;

	}
	return 1;
}

void PrintOutlierErr(DEROBS* obs) {
	MJDTIME Mjd;

	GPSTimeToMJDTime(&FrameTime, &Mjd);
	printf("DerObs outlier: %12.5f %6d %10.1lf %3d %3d C %10.3lf %2d %2d O %10.3lf %2d %2d\n",
		Mjd.Days + Mjd.FracDay, FrameTime.Week, FrameTime.SecOfWeek,
		obs->Scid1, obs->Scid2,
		obs->ObsQua.ApriClkResid, obs->ClkBlunder[0],
		obs->ClkBlunder[1], obs->ObsQua.ApriOrbResid,
		obs->OrbBlunder[0], obs->OrbBlunder[1]);

	fprintf(FPLOG, "DerObs outlier: %12.5f %6d %10.1lf %3d %3d C %10.3lf %2d %2d O %10.3lf %2d %2d\n",
		Mjd.Days + Mjd.FracDay, FrameTime.Week, FrameTime.SecOfWeek,
		obs->Scid1, obs->Scid2,
		obs->ObsQua.ApriClkResid, obs->ClkBlunder[0],
		obs->ClkBlunder[1], obs->ObsQua.ApriOrbResid,
		obs->OrbBlunder[0], obs->OrbBlunder[1]);
} 

void CheckOutlierObs_LG3(SATINFO* SatAtod) {
	int i, num;
	double MeanClk, MeanOrb, StdClk, StdOrb, StdClk1, StdOrb1;
	double* CO_C, * ANO_C;
	double ChkClk, ChkOrb;
	DEROBS* EpkDerObs;
	EdgeList* el;

	num = SatAtod->in_num + SatAtod->out_num;
	CO_C = (double*)malloc(num * sizeof(double));
	ANO_C = (double*)malloc(num * sizeof(double));

	memset(CO_C, 0, sizeof(double) * num);
	memset(ANO_C, 0, sizeof(double) * num);

	i = 0;
	for (el = (EdgeList*)SatAtod->in_edges.next; el != &SatAtod->in_edges; el = (EdgeList*)el->next) {
		EpkDerObs = (DEROBS*)el->edge;
		CO_C[i] = EpkDerObs->ObsQua.ApriClkResid;
		ANO_C[i] = EpkDerObs->ObsQua.ApriOrbResid;
		i++;
	}
	for (el = (EdgeList*)SatAtod->out_edges.next; el != &SatAtod->out_edges; el = (EdgeList*)el->next) {
		EpkDerObs = (DEROBS*)el->edge;
		CO_C[i] = EpkDerObs->ObsQua.ApriClkResid;
		ANO_C[i] = EpkDerObs->ObsQua.ApriOrbResid;
		i++;
	}

	CompVectStat(num, CO_C, &MeanClk, &StdClk, SatAtod);
	CompVectStat(num, ANO_C, &MeanOrb, &StdOrb, SatAtod);

	StdClk1 = max(StdClk, ClkDiff);
	StdOrb1 = max(StdOrb, OrbDiff);

	for (el = (EdgeList*)SatAtod->in_edges.next; el != &SatAtod->in_edges; el = (EdgeList*)el->next) {
		EpkDerObs = (DEROBS*)el->edge;
		ChkClk = EpkDerObs->ObsQua.ApriClkResid - MeanClk;
		ChkOrb = EpkDerObs->ObsQua.ApriOrbResid - MeanOrb;

		if (fabs(ChkClk) < 3.0 * StdClk1)	EpkDerObs->ClkBlunder[0] = 1;
		else								EpkDerObs->ClkBlunder[0] = 2;
		if (fabs(ChkOrb) < 3.0 * StdOrb1)	EpkDerObs->OrbBlunder[0] = 1;
		else								EpkDerObs->OrbBlunder[0] = 2;
		if (EpkDerObs->ClkBlunder[0] == 1 && EpkDerObs->OrbBlunder[0] == 1) {
			if (EpkDerObs->Valid != 3)
				EpkDerObs->Valid = 1;
		}
		else if (EpkDerObs->Valid != 1){
			EpkDerObs->Valid = -1;
		}
	}

	for (el = (EdgeList*)SatAtod->out_edges.next; el != &SatAtod->out_edges; el = (EdgeList*)el->next) {
		EpkDerObs = (DEROBS*)el->edge;
		ChkClk = EpkDerObs->ObsQua.ApriClkResid - MeanClk;
		ChkOrb = EpkDerObs->ObsQua.ApriOrbResid - MeanOrb;

		if (fabs(ChkClk) < 3.0 * StdClk1)	EpkDerObs->ClkBlunder[1] = 1;
		else								EpkDerObs->ClkBlunder[1] = 2;
		if (fabs(ChkOrb) < 3.0 * StdOrb1)	EpkDerObs->OrbBlunder[1] = 1;
		else								EpkDerObs->OrbBlunder[1] = 2;
		if (EpkDerObs->ClkBlunder[1] == 1 && EpkDerObs->OrbBlunder[1] == 1) {
			if (EpkDerObs->Valid != 3)
				EpkDerObs->Valid = 1;
		}
		else if (EpkDerObs->Valid != 1){
			EpkDerObs->Valid = -1;
		}
	}
	SatAtod->MeanClk_apr = MeanClk;
	SatAtod->StdClk_apr = StdClk;
	SatAtod->MeanOrb_apr = MeanOrb;
	SatAtod->StdOrb_apr = StdOrb;
	free(CO_C);
	free(ANO_C);
}

void CheckOutlierObs_LE3(SATINFO* SatAtod) {
	int n, num;
	double MeanClk, MeanOrb, StdClk, StdOrb;
	double ChkClk, ChkOrb;
	DEROBS* EpkDerObs;
	EdgeList* el;

	num = SatAtod->in_num + SatAtod->out_num;

	n = 0;
	MeanClk = MeanOrb = StdClk = StdOrb = 0.0;

	for (el = (EdgeList*)SatAtod->in_edges.next; el != &SatAtod->in_edges; el = (EdgeList*)el->next) {
		EpkDerObs = (DEROBS*)el->edge;
		if ((EpkDerObs->ClkBlunder[0] == 1 && EpkDerObs->OrbBlunder[0] == 1) ||
			(EpkDerObs->ClkBlunder[1] == 1 && EpkDerObs->OrbBlunder[1] == 1))
		{
			MeanClk += EpkDerObs->ObsQua.ApriClkResid;
			MeanOrb += EpkDerObs->ObsQua.ApriOrbResid;

			StdClk += EpkDerObs->ObsQua.ApriClkResid * EpkDerObs->ObsQua.ApriClkResid;
			StdOrb += EpkDerObs->ObsQua.ApriOrbResid * EpkDerObs->ObsQua.ApriOrbResid;
			n++;
		}
	}

	for (el = (EdgeList*)SatAtod->in_edges.next; el != &SatAtod->in_edges; el = (EdgeList*)el->next) {
		EpkDerObs = (DEROBS*)el->edge;
		if ((EpkDerObs->ClkBlunder[0] == 1 && EpkDerObs->OrbBlunder[0] == 1) ||
			(EpkDerObs->ClkBlunder[1] == 1 && EpkDerObs->OrbBlunder[1] == 1))
		{
			MeanClk += EpkDerObs->ObsQua.ApriClkResid;
			MeanOrb += EpkDerObs->ObsQua.ApriOrbResid;

			StdClk += EpkDerObs->ObsQua.ApriClkResid * EpkDerObs->ObsQua.ApriClkResid;
			StdOrb += EpkDerObs->ObsQua.ApriOrbResid * EpkDerObs->ObsQua.ApriOrbResid;
			n++;
		}
	}


	if (n > 1)
	{
		MeanOrb /= n;
		MeanClk /= n;
		StdClk = max(sqrt(StdClk / n - MeanClk * MeanClk), ClkDiff);
		StdOrb = max(sqrt(StdOrb / n - MeanOrb * MeanOrb), OrbDiff);
	}
	else if (n == 1) {
		MeanClk = MeanOrb = 0.0;
		StdOrb = 2.5;
		StdClk = 2.1;
	}
	else {
		printf("n=0\n");
		return;
	}

	for (el = (EdgeList*)SatAtod->in_edges.next; el != &SatAtod->in_edges; el = (EdgeList*)el->next) {
		EpkDerObs = (DEROBS*)el->edge;
		ChkClk = EpkDerObs->ObsQua.ApriClkResid - MeanClk;
		ChkOrb = EpkDerObs->ObsQua.ApriOrbResid - MeanOrb;

		if (fabs(ChkClk) < 3.0 * StdClk)	EpkDerObs->ClkBlunder[0] = 1;
		else								EpkDerObs->ClkBlunder[0] = 2;
		if (fabs(ChkOrb) < 3.0 * StdOrb)	EpkDerObs->OrbBlunder[0] = 1;
		else								EpkDerObs->OrbBlunder[0] = 2;
		if (EpkDerObs->ClkBlunder[0] == 1 && EpkDerObs->OrbBlunder[0] == 1) {
			if (EpkDerObs->Valid != 3)
				EpkDerObs->Valid = 1;
		}
		else if (EpkDerObs->Valid != 1){
			EpkDerObs->Valid = -1;
		}
	}

	for (el = (EdgeList*)SatAtod->out_edges.next; el != &SatAtod->out_edges; el = (EdgeList*)el->next) {
		EpkDerObs = (DEROBS*)el->edge;
		ChkClk = EpkDerObs->ObsQua.ApriClkResid - MeanClk;
		ChkOrb = EpkDerObs->ObsQua.ApriOrbResid - MeanOrb;

		if (fabs(ChkClk) < 3.0 * StdClk)	EpkDerObs->ClkBlunder[1] = 1;
		else								EpkDerObs->ClkBlunder[1] = 2;
		if (fabs(ChkOrb) < 3.0 * StdOrb)	EpkDerObs->OrbBlunder[1] = 1;
		else								EpkDerObs->OrbBlunder[1] = 2;
		if (EpkDerObs->ClkBlunder[1] == 1 && EpkDerObs->OrbBlunder[1] == 1) {
			if (EpkDerObs->Valid != 3)
				EpkDerObs->Valid = 1;
		}
		else if (EpkDerObs->Valid != 1){
			EpkDerObs->Valid = -1;
		}
	}
}

void WriteDerObsResidual(SATINFO* SatAtod)
{
	int j;
	MJDTIME mjd;
	GPSTIME  Time;
	SATINFO* sat;
	DEROBS* EpkDerObs;
	EdgeList* el;

	Time.Week = FrameTime.Week;
	Time.SecOfWeek = FrameTime.SecOfWeek + StepOfAutoNav;
	CheckGPSTime(&Time);

	sat = SatAtod;

	for (j = 0; j < SatNum; j++) {
		sat = (SATINFO*)sat->next;
		GPSTimeToMJDTime(&sat->TOE, &mjd);

		fprintf(FPProc, "%12.5lf %2d C %10.1lf %9.3lf %9.3lf %3d %8.1lf", mjd.Days + mjd.FracDay, sat->id, sqrt(sat->CovC[0]) * C_Light, sat->MeanClk_apr, sat->StdClk_apr, sat->Valid, sat->GapTime);

		//fprintf(FDerClk,"%12.5f %2d ",mjd.Days+mjd.FracDay,SatAtod[j].SCID);
		fprintf(FDerClk, "%4d %6.0lf %4d ", Time.Week, Time.SecOfWeek, sat->id);

		for (el = (EdgeList*)sat->in_edges.next; el != &sat->in_edges; el = (EdgeList*)el->next) {
			EpkDerObs = (DEROBS*)el->edge;
			//n = (EpkDerObs->DerObsList[i].Scid1 == SatAtod[j].SCID)? 0 : 1;
			fprintf(FPProc, "%10.3lf %2d %2d %3d %3d",
				EpkDerObs->ObsQua.ApriClkResid, EpkDerObs->ClkBlunder[0],
				EpkDerObs->ClkBlunder[1], EpkDerObs->Valid,
				((EpkDerObs->Scid1 == sat->id) ? EpkDerObs->Scid2 : EpkDerObs->Scid1));

			fprintf(FDerClk, "%4d %14.3lf %14.3lf %14.3lf %4d %4d %4d", ((EpkDerObs->Scid1 == sat->id) ? EpkDerObs->Scid2 : EpkDerObs->Scid1),
				EpkDerObs->DerCObs, EpkDerObs->CCObs,
				EpkDerObs->CConCorr, EpkDerObs->ClkBlunder[0], EpkDerObs->ClkBlunder[1],
				EpkDerObs->Valid);

		}
		for (el = (EdgeList*)sat->out_edges.next; el != &sat->out_edges; el = (EdgeList*)el->next) {
			EpkDerObs = (DEROBS*)el->edge;
			//n = (EpkDerObs->DerObsList[i].Scid1 == SatAtod[j].SCID)? 0 : 1;
			fprintf(FPProc, "%10.3lf %2d %2d %3d %3d",
				EpkDerObs->ObsQua.ApriClkResid, EpkDerObs->ClkBlunder[0],
				EpkDerObs->ClkBlunder[1], EpkDerObs->Valid,
				((EpkDerObs->Scid1 == sat->id) ? EpkDerObs->Scid2 : EpkDerObs->Scid1));

			fprintf(FDerClk, "%4d %14.3lf %14.3lf %14.3lf %4d %4d %4d", ((EpkDerObs->Scid1 == sat->id) ? EpkDerObs->Scid2 : EpkDerObs->Scid1),
				EpkDerObs->DerCObs, EpkDerObs->CCObs,
				EpkDerObs->CConCorr, EpkDerObs->ClkBlunder[0], EpkDerObs->ClkBlunder[1],
				EpkDerObs->Valid);

		}
		fprintf(FPProc, "\n");
		fprintf(FDerClk, "\n");

		fprintf(FPProc, "%12.5lf %2d O %10.1lf %9.3lf %9.3lf %3d %8.1lf",
			mjd.Days + mjd.FracDay, sat->id,
			sqrt(sat->CovX[0] + sat->CovX[7] + sat->CovX[14]), sat->MeanOrb_apr,
			sat->StdOrb_apr, sat->Valid,
			sat->GapTime);

		//fprintf(FDerOrb,"%10.3f %2d ",mjd.Days+mjd.FracDay,SatAtod[j].SCID);
		fprintf(FDerOrb, "%4d %6.0lf %4d ", Time.Week, Time.SecOfWeek, sat->id);

		for (el = (EdgeList*)sat->in_edges.next; el != &sat->in_edges; el = (EdgeList*)el->next) {
			EpkDerObs = (DEROBS*)el->edge;
			//n = (EpkDerObs->DerObsList[i].Scid1 == SatAtod[j].SCID)? 0 : 1;
			fprintf(FPProc, "%10.3lf %2d %2d %3d %3d",
				EpkDerObs->ObsQua.ApriOrbResid, EpkDerObs->OrbBlunder[0],
				EpkDerObs->OrbBlunder[1], EpkDerObs->Valid,
				((EpkDerObs->Scid1 == sat->id) ? EpkDerObs->Scid2 : EpkDerObs->Scid1));

			fprintf(FDerOrb, "%4d %14.3lf %14.3lf %14.3lf %4d %4d %4d",
				((EpkDerObs->Scid1 == sat->id) ? EpkDerObs->Scid2 : EpkDerObs->Scid1),
				EpkDerObs->DerANObs, EpkDerObs->CANObs,
				EpkDerObs->AConCorr, EpkDerObs->OrbBlunder[0],
				EpkDerObs->OrbBlunder[1], EpkDerObs->Valid);

		}
		for (el = (EdgeList*)sat->out_edges.next; el != &sat->out_edges; el = (EdgeList*)el->next) {
			EpkDerObs = (DEROBS*)el->edge;
			//n = (EpkDerObs->DerObsList[i].Scid1 == SatAtod[j].SCID)? 0 : 1;
			fprintf(FPProc, "%10.3lf %2d %2d %3d %3d",
				EpkDerObs->ObsQua.ApriOrbResid, EpkDerObs->OrbBlunder[0],
				EpkDerObs->OrbBlunder[1], EpkDerObs->Valid,
				((EpkDerObs->Scid1 == sat->id) ? EpkDerObs->Scid2 : EpkDerObs->Scid1));

			fprintf(FDerOrb, "%4d %14.3lf %14.3lf %14.3lf %4d %4d %4d",
				((EpkDerObs->Scid1 == sat->id) ? EpkDerObs->Scid2 : EpkDerObs->Scid1),
				EpkDerObs->DerANObs, EpkDerObs->CANObs,
				EpkDerObs->AConCorr, EpkDerObs->OrbBlunder[0],
				EpkDerObs->OrbBlunder[1], EpkDerObs->Valid);

		}
		fprintf(FPProc, "\n");
		fprintf(FDerOrb, "\n");
	}

	// for (j = 0; j < AncNum; j++)
	// {
	// 	if (Anchor[j].Valid == false)  continue;

	// 	fprintf(FPProc, "%12.5lf %2d C        0.0     0.0     0.0    1    0.0",
	// 		mjd.Days + mjd.FracDay, Anchor[j].StnId);

	// 	//fprintf(FDerClk,"%12.5lf %2d C ",mjd.Days+mjd.FracDay,Anchor[j].StnId);
	// 	fprintf(FDerClk, "%4d %6.0lf %4d ", Time.Week, Time.SecOfWeek, Anchor[j].StnId);
	// 	for (el = (EdgeList*)Anchor[j].in_edges.next; el != &Anchor[j].in_edges; el = (EdgeList*)el->next) {
	// 		EpkDerObs = (DEROBS*)el->edge;
	// 		//n = (EpkDerObs->DerObsList[i].Scid1 == Anchor[j].StnId)? 0 : 1;

	// 		fprintf(FPProc, "%10.3lf %2d %2d %3d %3d",
	// 			EpkDerObs->ObsQua.ApriClkResid, EpkDerObs->ClkBlunder[0],
	// 			EpkDerObs->ClkBlunder[1], EpkDerObs->Valid,
	// 			((EpkDerObs->Scid1 == Anchor[j].StnId) ? EpkDerObs->Scid2 : EpkDerObs->Scid1));

	// 		fprintf(FDerClk, "%4d %14.3lf %14.3lf %14.3lf %4d %4d %4d", ((EpkDerObs->Scid1 == Anchor[j].StnId) ? EpkDerObs->Scid2 : EpkDerObs->Scid1),
	// 			EpkDerObs->DerCObs, EpkDerObs->CCObs,
	// 			EpkDerObs->CConCorr, EpkDerObs->ClkBlunder[0], EpkDerObs->ClkBlunder[1],
	// 			EpkDerObs->Valid);
	// 	}
	// 	for (el = (EdgeList*)Anchor[j].out_edges.next; el != &Anchor[j].out_edges; el = (EdgeList*)el->next) {
	// 		EpkDerObs = (DEROBS*)el->edge;
	// 		//n = (EpkDerObs->DerObsList[i].Scid1 == Anchor[j].StnId)? 0 : 1;

	// 		fprintf(FPProc, "%10.3lf %2d %2d %3d %3d",
	// 			EpkDerObs->ObsQua.ApriClkResid, EpkDerObs->ClkBlunder[0],
	// 			EpkDerObs->ClkBlunder[1], EpkDerObs->Valid,
	// 			((EpkDerObs->Scid1 == Anchor[j].StnId) ? EpkDerObs->Scid2 : EpkDerObs->Scid1));

	// 		fprintf(FDerClk, "%4d %14.3lf %14.3lf %14.3lf %4d %4d %4d", ((EpkDerObs->Scid1 == Anchor[j].StnId) ? EpkDerObs->Scid2 : EpkDerObs->Scid1),
	// 			EpkDerObs->DerCObs, EpkDerObs->CCObs,
	// 			EpkDerObs->CConCorr, EpkDerObs->ClkBlunder[0], EpkDerObs->ClkBlunder[1],
	// 			EpkDerObs->Valid);
	// 	}
	// 	fprintf(FPProc, "\n");
	// 	fprintf(FDerClk, "\n");

	// 	fprintf(FPProc, "%12.5lf %2d O        0.0     0.0     0.0     1    0.0", mjd.Days + mjd.FracDay, Anchor[j].StnId);

	// 	//fprintf(FDerOrb,"%10.3f %2d ",mjd.Days+mjd.FracDay,Anchor[j].StnId);
	// 	fprintf(FDerOrb, "%4d %6.0lf %4d ", Time.Week, Time.SecOfWeek, Anchor[j].StnId);
	// 	for (el = (EdgeList*)Anchor[j].in_edges.next; el != &Anchor[j].in_edges; el = (EdgeList*)el->next) {
	// 		EpkDerObs = (DEROBS*)el->edge;
	// 		//n = (EpkDerObs->DerObsList[i].Scid1 == Anchor[j].StnId)? 0 : 1;
	// 		fprintf(FPProc, "%10.3lf %2d %2d %3d %3d",
	// 			EpkDerObs->ObsQua.ApriOrbResid, EpkDerObs->OrbBlunder[0],
	// 			EpkDerObs->OrbBlunder[1], EpkDerObs->Valid,
	// 			((EpkDerObs->Scid1 == Anchor[j].StnId) ? EpkDerObs->Scid2 : EpkDerObs->Scid1));
	// 		fprintf(FDerOrb, "%4d %14.3lf %14.3lf %14.3lf %4d %4d %4d",
	// 			((EpkDerObs->Scid1 == Anchor[j].StnId) ? EpkDerObs->Scid2 : EpkDerObs->Scid1),
	// 			EpkDerObs->DerANObs, EpkDerObs->CANObs,
	// 			EpkDerObs->AConCorr, EpkDerObs->OrbBlunder[0],
	// 			EpkDerObs->OrbBlunder[1], EpkDerObs->Valid);
	// 	}
	// 	for (el = (EdgeList*)Anchor[j].out_edges.next; el != &Anchor[j].out_edges; el = (EdgeList*)el->next) {
	// 		EpkDerObs = (DEROBS*)el->edge;
	// 		//n = (EpkDerObs->DerObsList[i].Scid1 == Anchor[j].StnId)? 0 : 1;
	// 		fprintf(FPProc, "%10.3lf %2d %2d %3d %3d",
	// 			EpkDerObs->ObsQua.ApriOrbResid, EpkDerObs->OrbBlunder[0],
	// 			EpkDerObs->OrbBlunder[1], EpkDerObs->Valid,
	// 			((EpkDerObs->Scid1 == Anchor[j].StnId) ? EpkDerObs->Scid2 : EpkDerObs->Scid1));
	// 		fprintf(FDerOrb, "%4d %14.3lf %14.3lf %14.3lf %4d %4d %4d",
	// 			((EpkDerObs->Scid1 == Anchor[j].StnId) ? EpkDerObs->Scid2 : EpkDerObs->Scid1),
	// 			EpkDerObs->DerANObs, EpkDerObs->CANObs,
	// 			EpkDerObs->AConCorr, EpkDerObs->OrbBlunder[0],
	// 			EpkDerObs->OrbBlunder[1], EpkDerObs->Valid);
	// 	}
	// 	fprintf(FPProc, "\n");
	// 	fprintf(FDerOrb, "\n");
	// }
	fprintf(FDerOrb, "--------------------------------------------------\n");
	fprintf(FDerClk, "--------------------------------------------------\n");

	fflush(FDerOrb);
	fflush(FDerClk);
}

void CheckSatOutlier(SATINFO* sat) {

	sat->MeanClk_apr = sat->StdClk_apr = sat->MeanOrb_apr = sat->StdOrb_apr = 999.99;
	if (sat->in_num + sat->out_num > 3)	CheckOutlierObs_LG3(sat);
	else								CheckOutlierObs_LE3(sat);

	if (sat->Valid >= BREAKING && fabs(sat->MeanClk_apr) > 3.5 && sat->StdClk_apr < 2.0) {
		sat->CovC[0] += fabs(sat->MeanClk_apr) * 1.0E-10;
	}
}

void CheckAnchorOutlier(ANCHSTN* anc)
{
	int i, num;
	double MeanClk, MeanOrb, StdClk, StdOrb, StdClk1, StdOrb1;
	double* CO_C, * ANO_C;
	double ChkClk, ChkOrb;
	DEROBS* EpkDerObs;
	EdgeList* el;

	MeanClk = MeanOrb = StdClk = StdOrb = 0.0;
	num = anc->in_num + anc->out_num;

	CO_C = (double*)malloc(num * sizeof(double));
	ANO_C = (double*)malloc(num * sizeof(double));

	memset(CO_C, 0, sizeof(double) * num);
	memset(ANO_C, 0, sizeof(double) * num);

	i = 0;
	for (el = (EdgeList*)anc->in_edges.next ; el != &anc->in_edges; el = (EdgeList*)el->next) {
		EpkDerObs = (DEROBS*)el->edge;
		CO_C[i] = EpkDerObs->ObsQua.ApriClkResid;
		ANO_C[i] = EpkDerObs->ObsQua.ApriOrbResid;
		i++;
	}
	for (el = (EdgeList*)anc->out_edges.next ; el != &anc->out_edges; el = (EdgeList*)el->next) {
		EpkDerObs = (DEROBS*)el->edge;
		CO_C[i] = EpkDerObs->ObsQua.ApriClkResid;
		ANO_C[i] = EpkDerObs->ObsQua.ApriOrbResid;
		i++;
	}

	if (num > 3)
	{
		//ClkAccu = 2.1;
		//OrbAccu = 2.5;
		CompVectStat(num, CO_C, &MeanClk, &StdClk, NULL);
		CompVectStat(num, ANO_C, &MeanOrb, &StdOrb, NULL);

		StdClk1 = max(StdClk, ClkDiff);
		StdOrb1 = max(StdOrb, OrbDiff);
		for (el = (EdgeList*)anc->in_edges.next ; el != &anc->in_edges; el = (EdgeList*)el->next) {
			EpkDerObs = (DEROBS*)el->edge;
			ChkClk = EpkDerObs->ObsQua.ApriClkResid - MeanClk;
			ChkOrb = EpkDerObs->ObsQua.ApriOrbResid - MeanOrb;

			if (fabs(ChkClk) < 3.0 * StdClk1)   EpkDerObs->ClkBlunder[0] = 1;
			else                           		EpkDerObs->ClkBlunder[0] = 2;
			if (fabs(ChkOrb) < 3.0 * StdOrb1)   EpkDerObs->OrbBlunder[0] = 1;
			else                           		EpkDerObs->OrbBlunder[0] = 2;
			if (EpkDerObs->ClkBlunder[0] == 1 && EpkDerObs->OrbBlunder[0] == 1) {
				if (EpkDerObs->Valid != 3)
					EpkDerObs->Valid = 1;
			}
			else if (EpkDerObs->Valid != 1){
				EpkDerObs->Valid = -1;
			}
		}
		for (el = (EdgeList*)anc->out_edges.next ; el != &anc->out_edges; el = (EdgeList*)el->next) {
			EpkDerObs = (DEROBS*)el->edge;
			ChkClk = EpkDerObs->ObsQua.ApriClkResid - MeanClk;
			ChkOrb = EpkDerObs->ObsQua.ApriOrbResid - MeanOrb;

			if (fabs(ChkClk) < 3.0 * StdClk1)   EpkDerObs->ClkBlunder[1] = 1;
			else                           		EpkDerObs->ClkBlunder[1] = 2;
			if (fabs(ChkOrb) < 3.0 * StdOrb1)   EpkDerObs->OrbBlunder[1] = 1;
			else                           		EpkDerObs->OrbBlunder[1] = 2;
			if (EpkDerObs->ClkBlunder[0] == 1 && EpkDerObs->OrbBlunder[0] == 1) {
				if (EpkDerObs->Valid != 3)
					EpkDerObs->Valid = 1;
			}
			else if (EpkDerObs->Valid != 1){
				EpkDerObs->Valid = -1;
			}
		}
	}
	free(CO_C);
	free(ANO_C);
}

void DectectDerObsOutlier(SATNET* SatNet)
{
	int i, id1, id2, k;
	Point* p;
	DEROBS* EpkDerObs;
	DEROBS* obs, * obs1;

	EpkDerObs = (DEROBS*)SatNet->edges;

	for (p = (Point*)SatNet->points->next; p != SatNet->points; p = (Point*)p->next) {
		if (p->type == 1) 						CheckSatOutlier((SATINFO*)p);
		else if (((ANCHSTN*)p)->Valid == 1)		CheckAnchorOutlier((ANCHSTN*)p);
	}

	// for (obs = (DEROBS*)EpkDerObs->next; obs != EpkDerObs; obs = (DEROBS*)obs->next) {
	// 	if (obs->Valid != 1)   continue;

	// 	id1 = obs->Scid1;
	// 	id2 = obs->Scid2;
	// 	for (obs1 = (DEROBS*)obs->next; obs1 != EpkDerObs; obs1 = (DEROBS*)obs1->next) {
	// 		if (obs1->Valid != 1)  continue;

	// 		if ((obs1->Scid1 == id1 || obs1->Scid2 == id1) && (obs1->Scid1 == id2 || obs1->Scid2 == id2)) {
	// 			obs->Valid = 3;
	// 		}
	// 	}
	// }

	WriteDerObsResidual((SATINFO*)SatNet->points);

}

void Dyadic(int m, int n, const double A[], const double B[], double Mat[])
{
	int i, j;

	if (m < 1 || n < 1)
	{
		printf("Dyadic fail: n<1 or m<1. \n");
		return;
	}

	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			*(Mat + i * n + j) = *(A + i) * *(B + j);
		}
	}
}

int ScalarMeasUpdate(double O_C, double sigma2, SpraseMatrix* H, int Scale, int size, double X[], SpraseMatrix* AllCov) {
	int	i, j;
	double	Error;

	SpraseMatrix* K, * Mat, * Cov;
	K = (SpraseMatrix*)malloc(sizeof(SpraseMatrix));
	Mat = (SpraseMatrix*)malloc(sizeof(SpraseMatrix));
	Cov = (SpraseMatrix*)malloc(sizeof(SpraseMatrix));

	CreateSpraseMatrix(K, SatNum, 1, size, 1);
	CreateSpraseMatrix(Mat, SatNum, SatNum, size, size);
	CreateSpraseMatrix(Cov, SatNum, SatNum, size, size);

	SpraseMatrixMultiply(AllCov, H, K);
	Error = sigma2 + SpraseVectDot(H, K);

	if ((Scale == 1 || Scale == 3) && fabs(O_C) / Error > 3.0) {
		return false;
	}

	if (Scale == 3)  return true;

	SpraseMatrixNode* p;
	p = K->chead[0].down;
	while (p != NULL) {
		for (i = 0; i < size; i++)
			*(p->submatrix + i) = *(p->submatrix + i) / Error;
		for (i = 0; i < size; i++)
			X[p->i * size + i] += *(p->submatrix + i) * O_C;
		p = p->down;
	}

	double* buff_;
	SpraseMatrixNode* pa, * pb;

	buff_ = (double*)malloc(sizeof(double) * size * size);
	pa = K->chead[0].down;
	while (pa != NULL) {
		pb = H->chead[0].down;
		while (pb != NULL) {
			Dyadic(size, size, pa->submatrix, pb->submatrix, buff_);
			for (i = 0; i < size; i++) {
				for (j = 0; j < size; j++) {
					*(buff_ + i * size + j) = 0 - *(buff_ + i * size + j);
				}
			}
			CopySubSpraseMatrix(Mat, pa->i, pb->i, size, size, buff_);
			pb = pb->down;
		}
		pa = pa->down;
	}
	free(buff_);
	SpraseMatrixMultiply(Mat, AllCov, Cov);
	SpraseMatrixAddition1(Cov, AllCov);
	DestroySpraseMatrix(Cov);
	DestroySpraseMatrix(K);
	DestroySpraseMatrix(Mat);

	return true;
}

void ClkMeasUpdate(SATNET* SatNet) {
	int i;
	double	O_C = 0, R = 0, ClkRate[MAXSATNUM] = { 0 }, rms0 = 0, rms1 = 0;
	DEROBS* derobs;
	CONSTSTATE* AllSatCov;
	SATINFO* sat;
	ANCHSTN* anc;
	double  Clk[2];
	double  CovC[4];
	SATINFO* sat1;
	SATINFO* sat2;

	sat = (SATINFO*)SatNet->points;
	AllSatCov = &SatNet->AllSatCov;
	derobs = (DEROBS*)SatNet->edges;
	for (i = 0; i < SatNum; i++) {
		sat = (SATINFO*)sat->next;
		if (sat->Valid > NOINIT) 	ClkRate[i] = sat->Clk[1];
	}

	Point* p;
	EdgeList* el;

	SpraseMatrix* H;
	H = (SpraseMatrix*)malloc(sizeof(SpraseMatrix));
	CreateSpraseMatrix(H, SatNum, 1, 2, 1);

	double H_sub[2];

	for (p = (Point*)SatNet->points->next; p != SatNet->points; p = (Point*)p->next) {
		// anc
		if (p->type == 0) {
			anc = (ANCHSTN*)p;
			if (anc->RefClkFlag != 1)  continue;
			for (el = (EdgeList*)anc->in_edges.next; el != &anc->in_edges; el = (EdgeList*)el->next) {
				derobs = (DEROBS*)el->edge;
				if (derobs->Valid < 1)  continue;
				//sat = (SATINFO*)derobs->endpoints[1];
				sat = (SATINFO*)GetPoint(derobs->Scid2);
				if (sat == NULL) continue;
				if (sat->type != 1 || sat->Valid <= NOINIT || sat->Health == 0) continue;
				EmptySpraseMatrix(H);
				O_C = derobs->DerCObs - derobs->CCObs - derobs->CConCorr
					- (2.0 * AllSatCov->Clk[sat->index * 2] + AllSatCov->Clk[sat->index * 2 + 1] * derobs->ObsQua.dt[0]) * C_Light;

				H_sub[0] = 2.0 * C_Light;
				H_sub[1] = derobs->ObsQua.dt[0] * C_Light;
				CopySubSpraseMatrix(H, sat->index, 0, 2, 1, H_sub);

				R = NoiseOfISL * NoiseOfISL;

				if (ScalarMeasUpdate(O_C, R, H, max((int)derobs->Valid, 2), 2, AllSatCov->Clk, AllSatCov->ClkCov) == false) {
					derobs->Valid = -1;

					printf("Clock MeasUpdate fail: SCID %2d %6d %10.1f %10.2f %8.3lf  Ref SCID:%2d\n",
						sat->id, sat->TOE.Week,
						sat->TOE.SecOfWeek, O_C,
						R, anc->StnId);

					fprintf(FPLOG, "Clock MeasUpdate fail: SCID %2d %6d %10.1f %10.2f %8.3lf  Ref SCID:%2d\n",
						sat->id, sat->TOE.Week,
						sat->TOE.SecOfWeek, O_C,
						R, anc->StnId);
				}
			}
		}
		else {
			sat1 = (SATINFO*)p;
			if (sat1->Valid <= NOINIT || sat1->Health == 0)    continue;
			for (el = (EdgeList*)sat1->in_edges.next; el != &sat1->in_edges; el = (EdgeList*)el->next) {
				derobs = (DEROBS*)el->edge;
				if (derobs->Valid < 1)  continue;
				if (derobs->endpoints[1]->type == 1) {
					EmptySpraseMatrix(H);
					//sat2 = (SATINFO*)derobs->endpoints[1];
					sat2 = (SATINFO*)GetPoint(derobs->Scid2);
					if (sat2 == NULL) continue;
					if (sat2->Valid <= NOINIT || sat2->Health == 0)    continue;
					O_C = derobs->DerCObs - derobs->CCObs - derobs->CConCorr
						- (2.0 * AllSatCov->Clk[sat1->index * 2] + AllSatCov->Clk[sat1->index * 2 + 1] * derobs->ObsQua.dt[0]) * C_Light
						+ (2.0 * AllSatCov->Clk[sat2->index * 2] + AllSatCov->Clk[sat2->index * 2 + 1] * derobs->ObsQua.dt[1]) * C_Light;

					H_sub[0] = 2.0 * C_Light;
					H_sub[1] = derobs->ObsQua.dt[0] * C_Light;
					CopySubSpraseMatrix(H, sat1->index, 0, 2, 1, H_sub);
					//H[id2*2+0]	= -2.0*C_Light;
					//H[id2*2+1]	= -1.0*EpkDerObs->DerObsList[i].ObsQua.dt[1]*C_Light;

					R = NoiseOfISL * NoiseOfISL;

					if (ScalarMeasUpdate(O_C, R, H, derobs->Valid, 2, AllSatCov->Clk, AllSatCov->ClkCov) == false) {
						derobs->Valid = -1;

						printf("Clock MeasUpdate fail: %2d %6d %10.1f %10.2f %8.2f  Ref SCID: %2d\n",
							derobs->Scid1, sat1->TOE.Week,
							sat1->TOE.SecOfWeek, O_C,
							R, derobs->Scid2);

						fprintf(FPLOG, "Clock MeasUpdate fail: %2d %6d %10.1f %10.2f %8.2f  Ref SCID: %2d\n",
							derobs->Scid1, sat1->TOE.Week,
							sat1->TOE.SecOfWeek, O_C,
							R, derobs->Scid2);
					}
				}
				else continue;
			}
		}
	}

	DestroySpraseMatrix(H);

	int j, k;
	double ANO_C[(MAXSATNUM + MAXANCHORNUM) * MAXANTENNA];

	sat = (SATINFO*)SatNet->points;
	for (i = 0; i < SatNum; i++) {
		sat = (SATINFO*)sat->next;
		memcpy(Clk, sat->Clk, 2 * sizeof(double));
		memcpy(CovC, sat->CovC, 4 * sizeof(double));
		CopyArray(2, sat->Clk, AllSatCov->Clk + 2 * i);
		GetSubSpraseMatrix(AllSatCov->ClkCov, i, i, 2, 2, sat->CovC);

		if (sat->Valid <= NOINIT)	continue;
		{
			sat->Clk[1] = 0.99999 * ClkRate[sat->index] + 0.00001 * sat->Clk[1];

			if (fmod(sat->TOE.SecOfWeek + 0.01, 300.0) < 0.5)
			{
				sat->SatClk.ClkSeq[sat->SatClk.CurNum] = sat->Clk[0];
				sat->SatClk.Time->Week = sat->TOE.Week;
				sat->SatClk.Time->SecOfWeek = sat->TOE.SecOfWeek;
				sat->SatClk.CurNum = (sat->SatClk.CurNum + 1) % MAXCLKSER;
				sat->SatClk.TotalNum++;
			}
		}

		if (sat->Health == 0)	continue;
		j = 0;
		for (el = (EdgeList*)sat->in_edges.next; el != &sat->in_edges; el = (EdgeList*)el->next) {
			derobs = (DEROBS*)el->edge;
			sat1 = (SATINFO*)derobs->endpoints[0];
			sat2 = (SATINFO*)derobs->endpoints[1];

			if (derobs->Valid < 1)
				continue;

			if (sat2->type == 1) {
				if (sat2->Valid <= NOINIT || sat2->Health == 0)		continue;

				ANO_C[j] = derobs->DerCObs - derobs->CCObs - derobs->CConCorr
					- (2.0 * sat->Clk[0] + sat->Clk[1] * derobs->ObsQua.dt[0]) * C_Light
					+ (2.0 * sat2->Clk[0] + sat2->Clk[1] * derobs->ObsQua.dt[1]) * C_Light;
				j++;
			}

		}
		for (el = (EdgeList*)sat->out_edges.next; el != &sat->out_edges; el = (EdgeList*)el->next) {
			derobs = (DEROBS*)el->edge;
			sat1 = (SATINFO*)derobs->endpoints[0];
			sat2 = (SATINFO*)derobs->endpoints[1];

			if (derobs->Valid < 1)
				continue;

			if (sat1->type == 1) {
				if (sat1->Valid <= NOINIT || sat1->Health == 0)		continue;

				ANO_C[j] = derobs->DerCObs - derobs->CCObs - derobs->CConCorr
					- (2.0 * sat1->Clk[0] + sat1->Clk[1] * derobs->ObsQua.dt[0]) * C_Light
					+ (2.0 * sat->Clk[0] + sat->Clk[1] * derobs->ObsQua.dt[1]) * C_Light;
				j++;
			}
		}
		CompVectStat(j, ANO_C, &(sat->MeanClk_pst), &(sat->StdClk_pst), sat);
		rms0 = sqrt(sat->MeanClk_apr * sat->MeanClk_apr + sat->StdClk_apr * sat->StdClk_apr);
		rms1 = sqrt(sat->MeanClk_pst * sat->MeanClk_pst + sat->StdClk_pst * sat->StdClk_pst);

		if ((rms1 - rms0) > NoiseOfISL)
		{
			CopyArray(2, sat->Clk, Clk);
			CopyArray(4, sat->CovC, CovC);
			sat->SatClk.TotalNum--;
			if (sat->SatClk.CurNum == 0)   sat->SatClk.CurNum = MAXCLKSER - 1;
			else                           sat->SatClk.CurNum--;
		}
	}
}

void mbbub(int n, double p[])
{
	int m, k, i, j;
	double d;

	if (n <= 0)
	{
		printf("mbbub fail: m<1. \n");
		return;
	}

	k = 0;
	m = n - 1;

	while (k < m)
	{
		j = m - 1;
		m = 0;
		for (i = k; i <= j; i++)
		{
			if (p[i] > p[i + 1])
			{
				d = p[i];
				p[i] = p[i + 1];
				p[i + 1] = d;
				m = i;
			}
		}

		j = k + 1;
		k = 0;
		for (i = m; i >= j; i--)
		{
			if (p[i - 1] > p[i])
			{
				d = p[i];
				p[i] = p[i - 1];
				p[i - 1] = d;
				k = i;
			}
		}
	}
}

int CompVectStat(const int n, const double Dat[], double* Mean, double* Std, SATINFO* SatAtod)
{
	int  i, j, k, m;
	double Limit, Mean1;
	double Val[(MAXSATNUM + MAXANCHORNUM) * MAXANTENNA];

	if (n <= 3)   return false;

	if (SatAtod == NULL)   Limit = 300.0;
	else
	{
		switch (SatAtod->Valid)
		{
		case BREAKING:
			Limit = 300.0 * exp(SatAtod->GapTime / SatAtod->MeasStep);
			break;
		case MANEUVER:
			Limit = 1E20;
			break;
		case ANSOK:
			Limit = 300.0;
			break;
		default:
			Limit = 1E20;
			break;
		}
	}

	CopyArray(n, Val, Dat);
	mbbub(n, Val);

	m = n / 2;
	Mean1 = Val[n / 2];

	j = 1;
	while (m - j > 0 && m + j < n)
	{
		k = 0;
		*Mean = *Std = 0.0;
		for (i = m - j; i <= m + j; i++)
		{
			if (i == m)  continue;
			if (fabs(Val[i] - Mean1) < Limit)
			{
				*Mean += Val[i];
				*Std += (Val[i] - Mean1) * (Val[i] - Mean1);
				k++;
			}
		}
		if (k > 1)
		{
			Mean1 = *Mean / k;
			*Std = sqrt(*Std / (k - 1));
		}
		j++;
	}

	*Mean = Mean1;
	return true;
}

void ANSMeasUpdate(SATNET* SatNet) {

	int	   i = 0, j = 0, k = 0, n = 0;
	double dPos[DIM], O_C = 0, R = 0, H1[DIM], H2[DIM], H3[DIM], H4[DIM], H_[DIM];
	double X1[DIM], X2[DIM], X3[DIM], X4[DIM], Range1 = 0, Range2 = 0;
	double ANO_C[MAXSATNUM * MAXSATNUM];
	DEROBS* derobs;
	CONSTSTATE* AllSatCov;
	SATINFO* sat1;
	SATINFO* sat2;
	SATINFO* sat;
	ANCHSTN* anc;

	derobs = (DEROBS*)SatNet->edges;
	AllSatCov = &SatNet->AllSatCov;
	memset(dPos, 0, DIM * sizeof(double));
	sat = (SATINFO*)SatNet->points;

	for (i = 0; i < SatNum; i++) {
		sat = (SATINFO*)sat->next;
		if (sat->Valid <= NOINIT)	continue;

		//sat->TotalSatNum = 0;
		//sat->ValidSatNum = 0;
		memset(sat->dX, 0, DIM * sizeof(double));
	}

	memset(AllSatCov->Orb, 0, MAXSATNUM * DIM * sizeof(double));

	Point* p;
	EdgeList* el;

	SpraseMatrix* H;
	H = (SpraseMatrix*)malloc(sizeof(SpraseMatrix));
	CreateSpraseMatrix(H, SatNum, 1, DIM, 1);

	for (p = (Point*)SatNet->points->next; p != SatNet->points; p = (Point*)p->next) {
		if (p->type == 1) {
			sat1 = (SATINFO*)p;
			if (sat1->Valid <= NOINIT || sat1->Health == 0 || sat1->in_num + sat1->out_num < 2)    continue;
			for (el = (EdgeList*)sat1->in_edges.next; el != &sat1->in_edges; el = (EdgeList*)el->next) {
				derobs = (DEROBS*)el->edge;
				if (derobs->Valid < 1)  continue;
				if (derobs->endpoints[1]->type == 1) {
					//sat2 = (SATINFO*)derobs->endpoints[1];
					sat2 = (SATINFO*)GetPoint(derobs->Scid2);
					if (sat2 == NULL) {
						sat2 = (SATINFO*)derobs->endpoints[1];
						if (sat2->Valid <= NOINIT || sat2->Health == 0 || sat2->in_num + sat2->out_num < 2)    continue;
						memset(dPos, 0, DIM * sizeof(double));
						EmptySpraseMatrix(H);

						MatrixMultiply(DIM, DIM, DIM, 1, derobs->P1RvState + 6, AllSatCov->Orb + sat1->index * DIM, X1);
						MatrixMultiply(DIM, DIM, DIM, 1, derobs->P1TrState + 6, AllSatCov->Orb + sat1->index * DIM, X2);
						MatrixMultiply(DIM, DIM, DIM, 1, derobs->P2RvState + 6, AllSatCov->Orb + sat2->index * DIM, X3);
						MatrixMultiply(DIM, DIM, DIM, 1, derobs->P2TrState + 6, AllSatCov->Orb + sat2->index * DIM, X4);

						MatrixAddition2(1, DIM, derobs->P1RvState, X1);    // X1=X1+Phi*dX
						MatrixAddition2(1, DIM, derobs->P1TrState, X2);
						MatrixAddition2(1, DIM, derobs->P2RvState, X3);    // X1=X1+Phi*dX
						MatrixAddition2(1, DIM, derobs->P2TrState, X4);
						for (j = 0; j < 3; j++)      dPos[j] = X1[j] - X4[j];
						Range1 = sqrt(VectDot(3, 3, dPos, dPos));
						for (j = 0; j < 3; j++)      dPos[j] = dPos[j] / Range1;
						MatrixMultiply(1, DIM, DIM, DIM, dPos, derobs->P1RvState + 6, H1);

						for (j = 0; j < 3; j++)      dPos[j] = X2[j] - X3[j];
						Range2 = sqrt(VectDot(3, 3, dPos, dPos));
						for (j = 0; j < 3; j++)      dPos[j] = dPos[j] / Range2;
						MatrixMultiply(1, DIM, DIM, DIM, dPos, derobs->P1TrState + 6, H2);

						for (j = 0; j < DIM; j++) {
							H_[j] = H1[j] + H2[j];
						}
						CopySubSpraseMatrix(H, sat1->index, 0, DIM, 1, H_);

						O_C = derobs->DerANObs - Range1 - Range2 - derobs->AConCorr;
						R = NoiseOfISL * NoiseOfISL;

						if (ScalarMeasUpdate(O_C, R, H, derobs->Valid, DIM, AllSatCov->Orb, AllSatCov->OrbCov)) {
							//sat1->TotalSatNum++;
							//sat2->TotalSatNum++;
							//if (derobs->Valid != 3)
							//{
							//	sat1->ValidSatNum++;
							//	sat1->ValidSatNum++;
							//}
						}
						else {
							derobs->Valid = -1;

							printf("Orbit MeasUpdate fail: SCID %2d %6d %10.1f %10.2f %8.2f  Ref SCID: %2d\n",
								sat1->id, sat1->TOE.Week,
								sat1->TOE.SecOfWeek, O_C,
								R, sat2->id);

							fprintf(FPLOG, "Orbit MeasUpdate fail: SCID %2d %6d %10.1f %10.2f %8.2f  Ref SCID: %2d\n",
								sat1->id, sat1->TOE.Week,
								sat1->TOE.SecOfWeek, O_C,
								R, sat2->id);
						}

					}
					else {
						if (sat2->Valid <= NOINIT || sat2->Health == 0 || sat2->in_num + sat2->out_num < 2)    continue;
						memset(dPos, 0, DIM * sizeof(double));
						
						EmptySpraseMatrix(H);

						MatrixMultiply(DIM, DIM, DIM, 1, derobs->P1RvState + 6, AllSatCov->Orb + sat1->index * DIM, X1);
						MatrixMultiply(DIM, DIM, DIM, 1, derobs->P1TrState + 6, AllSatCov->Orb + sat1->index * DIM, X2);
						MatrixMultiply(DIM, DIM, DIM, 1, derobs->P2RvState + 6, AllSatCov->Orb + sat2->index * DIM, X3);
						MatrixMultiply(DIM, DIM, DIM, 1, derobs->P2TrState + 6, AllSatCov->Orb + sat2->index * DIM, X4);

						MatrixAddition2(1, DIM, derobs->P1RvState, X1);    // X1=X1+Phi*dX
						MatrixAddition2(1, DIM, derobs->P1TrState, X2);
						MatrixAddition2(1, DIM, derobs->P2RvState, X3);    // X1=X1+Phi*dX
						MatrixAddition2(1, DIM, derobs->P2TrState, X4);
						for (j = 0; j < 3; j++)      dPos[j] = X1[j] - X4[j];
						Range1 = sqrt(VectDot(3, 3, dPos, dPos));
						for (j = 0; j < 3; j++)      dPos[j] = dPos[j] / Range1;
						MatrixMultiply(1, DIM, DIM, DIM, dPos, derobs->P1RvState + 6, H1);
						for (j = 0; j < 3; j++)      dPos[j] = -1.0 * dPos[j];
						MatrixMultiply(1, DIM, DIM, DIM, dPos, derobs->P2TrState + 6, H4);

						for (j = 0; j < 3; j++)      dPos[j] = X2[j] - X3[j];
						Range2 = sqrt(VectDot(3, 3, dPos, dPos));
						for (j = 0; j < 3; j++)      dPos[j] = dPos[j] / Range2;
						MatrixMultiply(1, DIM, DIM, DIM, dPos, derobs->P1TrState + 6, H2);
						for (j = 0; j < 3; j++)      dPos[j] = -1.0 * dPos[j];
						MatrixMultiply(1, DIM, DIM, DIM, dPos, derobs->P2RvState + 6, H3);

						for (j = 0; j < DIM; j++) {
							H_[j] = H1[j] + H2[j];
						}
						CopySubSpraseMatrix(H, sat1->index, 0, DIM, 1, H_);
						for (j = 0; j < DIM; j++) {
							H_[j] = H3[j] + H4[j];
						}
						CopySubSpraseMatrix(H, sat2->index, 0, DIM, 1, H_);

						O_C = derobs->DerANObs - Range1 - Range2 - derobs->AConCorr;
						R = NoiseOfISL * NoiseOfISL;

						if (ScalarMeasUpdate(O_C, R, H, derobs->Valid, DIM, AllSatCov->Orb, AllSatCov->OrbCov)) {
							//sat1->TotalSatNum++;
							//sat2->TotalSatNum++;
							//if (derobs->Valid != 3)
							//{
							//	sat1->ValidSatNum++;
							//	sat1->ValidSatNum++;
							//}
						}
						else {
							derobs->Valid = -1;

							printf("Orbit MeasUpdate fail: SCID %2d %6d %10.1f %10.2f %8.2f  Ref SCID: %2d\n",
								sat1->id, sat1->TOE.Week,
								sat1->TOE.SecOfWeek, O_C,
								R, sat2->id);

							fprintf(FPLOG, "Orbit MeasUpdate fail: SCID %2d %6d %10.1f %10.2f %8.2f  Ref SCID: %2d\n",
								sat1->id, sat1->TOE.Week,
								sat1->TOE.SecOfWeek, O_C,
								R, sat2->id);
						}
					}
				}
				else continue;
			}
		}
		else {
			anc = (ANCHSTN*)p;
			if (anc->RefClkFlag != 1)  continue;
			for (el = (EdgeList*)anc->in_edges.next; el != &anc->in_edges; el = (EdgeList*)el->next) {
				derobs = (DEROBS*)el->edge;
				//if (derobs->Valid < 1)  continue;
				sat = (SATINFO*)derobs->endpoints[1];
				sat = (SATINFO*)GetPoint(derobs->Scid2);
				if (sat == NULL) continue;
				if (sat->type != 1 || sat->Valid <= NOINIT || sat->Health == 0 || sat->in_num + sat->out_num < 2) continue;
				memset(dPos, 0, DIM * sizeof(double));

				EmptySpraseMatrix(H);

				MatrixMultiply(DIM, DIM, DIM, 1, derobs->P2RvState + 6, AllSatCov->Orb + sat->index * DIM, X3);
				MatrixMultiply(DIM, DIM, DIM, 1, derobs->P2TrState + 6, AllSatCov->Orb + sat->index * DIM, X4);
				MatrixAddition2(1, DIM, derobs->P2RvState, X3);    // X1=X1+Phi*dX
				MatrixAddition2(1, DIM, derobs->P2TrState, X4);

				for (j = 0; j < 3; j++)      dPos[j] = X4[j] - derobs->P1RvState[j];
				Range1 = sqrt(VectDot(3, 3, dPos, dPos));
				for (j = 0; j < 3; j++)      dPos[j] = dPos[j] / Range1;
				MatrixMultiply(1, DIM, DIM, DIM, dPos, derobs->P2TrState + 6, H2);

				for (j = 0; j < 3; j++)      dPos[j] = X3[j] - derobs->P1TrState[j];
				Range2 = sqrt(VectDot(3, 3, dPos, dPos));
				for (j = 0; j < 3; j++)      dPos[j] = dPos[j] / Range2;
				MatrixMultiply(1, DIM, DIM, DIM, dPos, derobs->P2RvState + 6, H1);

				for (j = 0; j < DIM; j++)   H_[j] = H1[j] + H2[j];
				CopySubSpraseMatrix(H, sat->index, 0, DIM, 1, H_);
				O_C = derobs->DerANObs - Range1 - Range2 - derobs->AConCorr;
				R = NoiseOfISL * NoiseOfISL;

				if (ScalarMeasUpdate(O_C, R, H, max((int)derobs->Valid, 2), DIM, AllSatCov->Orb, AllSatCov->OrbCov)) {
					//sat->TotalSatNum++;
					//if (derobs->Valid != 3)  sat->ValidSatNum++;
				}
				else {
					derobs->Valid = -1;

					printf("Orbit MeasUpdate fail: SCID %2d %6d %10.1f %10.2f %8.2f  Ref SCID: %2d\n",
						sat->id, sat->TOE.Week,
						sat->TOE.SecOfWeek, O_C,
						R, anc->StnId);

					fprintf(FPLOG, "Orbit MeasUpdate fail: SCID %2d %6d %10.1f %10.2f %8.2f  Ref SCID: %2d\n",
						sat->id, sat->TOE.Week,
						sat->TOE.SecOfWeek, O_C,
						R, anc->StnId);
				}
			}
		}
	}
	DestroySpraseMatrix(H);

	sat = (SATINFO*)SatNet->points;
	for (i = 0; i < SatNum; i++) {
		sat = (SATINFO*)sat->next;

		if (sat->Valid <= NOINIT || sat->Health == 0)        continue;

		n = 0;

		for (el = (EdgeList*)sat->in_edges.next; el != &sat->in_edges; el = (EdgeList*)el->next) {
			derobs = (DEROBS*)el->edge;
			if (derobs->Valid < 1)   continue;

			sat1 = (SATINFO*)derobs->endpoints[0];
			sat2 = (SATINFO*)derobs->endpoints[1];

			if (sat2->type == 1) {
				if (sat2->Valid <= NOINIT || sat2->Health == 0)
					continue;

				MatrixMultiply(DIM, DIM, DIM, 1, derobs->P1RvState + 6, AllSatCov->Orb + sat1->index * DIM, X1);
				MatrixMultiply(DIM, DIM, DIM, 1, derobs->P1TrState + 6, AllSatCov->Orb + sat1->index * DIM, X2);
				MatrixMultiply(DIM, DIM, DIM, 1, derobs->P2RvState + 6, AllSatCov->Orb + sat2->index * DIM, X3);
				MatrixMultiply(DIM, DIM, DIM, 1, derobs->P2TrState + 6, AllSatCov->Orb + sat2->index * DIM, X4);
				MatrixAddition2(1, DIM, derobs->P1RvState, X1);    // X1=X1+Phi*dX
				MatrixAddition2(1, DIM, derobs->P1TrState, X2);
				MatrixAddition2(1, DIM, derobs->P2RvState, X3);    // X1=X1+Phi*dX
				MatrixAddition2(1, DIM, derobs->P2TrState, X4);

				for (k = 0; k < 3; k++)      dPos[k] = X1[k] - X4[k];
				Range1 = sqrt(VectDot(3, 3, dPos, dPos));
				for (k = 0; k < 3; k++)      dPos[k] = X3[k] - X2[k];
				Range2 = sqrt(VectDot(3, 3, dPos, dPos));
				ANO_C[n++] = derobs->DerANObs - Range1 - Range2 - derobs->AConCorr;
			}
			else {

				MatrixMultiply(DIM, DIM, DIM, 1, derobs->P1RvState + 6, AllSatCov->Orb + sat1->index * DIM, X1);
				MatrixMultiply(DIM, DIM, DIM, 1, derobs->P1TrState + 6, AllSatCov->Orb + sat1->index * DIM, X2);
				MatrixAddition2(1, DIM, derobs->P1RvState, X1);    // X1=X1+Phi*dX
				MatrixAddition2(1, DIM, derobs->P1TrState, X2);

				for (k = 0; k < 3; k++)      dPos[k] = X1[k] - derobs->P2TrState[k];
				Range1 = sqrt(VectDot(3, 3, dPos, dPos));
				for (k = 0; k < 3; k++)      dPos[k] = X2[k] - derobs->P2RvState[k];
				Range2 = sqrt(VectDot(3, 3, dPos, dPos));
				ANO_C[n++] = derobs->DerANObs - Range1 - Range2 - derobs->AConCorr;
			}
		}
		for (el = (EdgeList*)sat->out_edges.next; el != &sat->out_edges; el = (EdgeList*)el->next) {
			derobs = (DEROBS*)el->edge;
			if (derobs->Valid < 1)   continue;

			sat1 = (SATINFO*)derobs->endpoints[0];
			sat2 = (SATINFO*)derobs->endpoints[1];

			if (sat1->type == 1) {
				if (sat1->Valid <= NOINIT || sat1->Health == 0)
					continue;

				MatrixMultiply(DIM, DIM, DIM, 1, derobs->P1RvState + 6, AllSatCov->Orb + sat1->index * DIM, X1);
				MatrixMultiply(DIM, DIM, DIM, 1, derobs->P1TrState + 6, AllSatCov->Orb + sat1->index * DIM, X2);
				MatrixMultiply(DIM, DIM, DIM, 1, derobs->P2RvState + 6, AllSatCov->Orb + sat2->index * DIM, X3);
				MatrixMultiply(DIM, DIM, DIM, 1, derobs->P2TrState + 6, AllSatCov->Orb + sat2->index * DIM, X4);
				MatrixAddition2(1, DIM, derobs->P1RvState, X1);    // X1=X1+Phi*dX
				MatrixAddition2(1, DIM, derobs->P1TrState, X2);
				MatrixAddition2(1, DIM, derobs->P2RvState, X3);    // X1=X1+Phi*dX
				MatrixAddition2(1, DIM, derobs->P2TrState, X4);

				for (k = 0; k < 3; k++)      dPos[k] = X1[k] - X4[k];
				Range1 = sqrt(VectDot(3, 3, dPos, dPos));
				for (k = 0; k < 3; k++)      dPos[k] = X3[k] - X2[k];
				Range2 = sqrt(VectDot(3, 3, dPos, dPos));
				ANO_C[n++] = derobs->DerANObs - Range1 - Range2 - derobs->AConCorr;
			}
			else {

				MatrixMultiply(DIM, DIM, DIM, 1, derobs->P2RvState + 6, AllSatCov->Orb + sat2->index * DIM, X3);
				MatrixMultiply(DIM, DIM, DIM, 1, derobs->P2TrState + 6, AllSatCov->Orb + sat2->index * DIM, X4);
				MatrixAddition2(1, DIM, derobs->P2RvState, X3);    // X1=X1+Phi*dX
				MatrixAddition2(1, DIM, derobs->P2TrState, X4);

				for (k = 0; k < 3; k++)      dPos[k] = X4[k] - derobs->P1RvState[k];
				Range1 = sqrt(VectDot(3, 3, dPos, dPos));
				for (k = 0; k < 3; k++)      dPos[k] = X3[k] - derobs->P1TrState[k];
				Range2 = sqrt(VectDot(3, 3, dPos, dPos));
				ANO_C[n++] = derobs->DerANObs - Range1 - Range2 - derobs->AConCorr;
			}
		}

		if (n > 1) {
			CopyArray(DIM, sat->dX, AllSatCov->Orb + i * DIM);
			MatrixAddition2(1, DIM, sat->dX, sat->X);
			GetSubSpraseMatrix(AllSatCov->OrbCov, i, i, DIM, DIM, sat->CovX);
		}

		if (n > 1) {
			sat->GapTime = 0.0;

			R = sqrt(sat->CovX[0] + sat->CovX[7] + sat->CovX[14]);
			CompVectStat(n, ANO_C, &(sat->MeanOrb_pst), &(sat->StdOrb_pst), sat);
			Range1 = sqrt(sat->MeanOrb_pst * sat->MeanOrb_pst + sat->StdOrb_pst * sat->StdOrb_pst);
			sat->URA = (float)(max(sat->PDOP, 5.0) * Range1);

			if (sat->URA < 2.0)
			{
				sat->Health = 1;
				sat->URA = 2.0;
				if (sat->Valid == BREAKING)    sat->Valid = ANSOK;
			}
			if (sat->URA > 2.0 && sat->URA <= 5.0)           sat->Health = 2;
			if (sat->PDOP > 5.0 || Range1 > 5.0)    sat->Health = 3;
		}
	}
}

int ClockOffsetFitting(TIMESYC* SatClk, GPSTIME* TOC, double Clk[3])
{
	int i, sum;
	double A[MAXCLKSER * 3], B[MAXCLKSER];
	double AT[MAXCLKSER * 3], ATA[9], ATB[3], ATA_[9];
	double dt, sigma;

	if (SatClk->TotalNum >= MAXCLKSER)
		sum = MAXCLKSER;
	else
		sum = SatClk->TotalNum;


	if (sum < 10)
		return false;


	for (i = 0; i < sum; i++) {
		dt = GetDifGPSTime(SatClk->Time + i, TOC) / SECPERDAY;
		A[3 * i] = AT[i] = 1.0;
		A[3 * i + 1] = AT[sum + i] = dt;
		A[3 * i + 2] = AT[2 * sum + i] = dt * dt;
		B[i] = SatClk->ClkSeq[i];
	}

	MatrixMultiply(3, sum, sum, 3, AT, A, ATA);
	MatrixMultiply(3, sum, sum, 1, AT, B, ATB);
	MatrixInv(3, ATA, ATA_);
	MatrixMultiply(3, 3, 3, 1, ATA_, ATB, Clk);

	sigma = 0.0;
	for (i = 0; i < sum; i++) {
		B[i] = B[i] - Clk[0] - A[3 * i + 1] * Clk[1] - A[3 * i + 2] * Clk[2];
		sigma = sigma + B[i] * B[i];
	}

	sigma = sigma / (sum - 1);

	return true;
}

double EphPredictTime = 0;
void OutputSatOrbit(SATINFO* SatAtod) {
	MJDTIME Mjd;
	double dt, XECF[6], XECI[6];      // output ECF, to test EOP predictions
	double Ceof[3];
	GPSTIME T;
	double Clk;

	dt = SECPERHOUR * EphPredictTime;
	T = SatAtod->TOE;
	T.SecOfWeek += dt;
	CheckGPSTime(&T);

	if (SatAtod->Valid < BREAKING)    return;

	GPSTimeToMJDTime(&T, &Mjd);

	if (!Fres)
	{
		printf("AutoNave result file cannot be opened!\n");
		return;
	}
	else
	{
		if (fabs(dt) < 1.0E-8)
		{
			Clk = SatAtod->Clk[0];
			CopyArray(6, XECI, SatAtod->X);
		}
		else
		{
			if (ClockOffsetFitting(&SatAtod->SatClk, &SatAtod->TOE, Ceof) == false)
			{
				Ceof[0] = SatAtod->Clk[0];
				Ceof[1] = SatAtod->Clk[1];
				Clk = Ceof[0] + Ceof[1] * dt;
			}
			else
			{
				Clk = Ceof[0] + Ceof[1] * dt / SECPERDAY + Ceof[2] * pow(dt / SECPERDAY, 2.0);
			}

			CopyArray(6, XECI, SatAtod->X);
			OrbitIntegToGivenTime(&SatAtod->TOE, &T, 120.0, XECI);
		}

		ICRF_ITRF_GPST(MJD_J2000, &T, 1, XECI, XECF);

		fprintf(Fres, "%12.5lf %3d %2d %14.3lf %14.3lf %14.3lf %14.5lf %14.5lf %14.5lf %14.3lf %14.3lf %14.3lf %14.5lf %14.5lf %14.5lf %14.6lf %17.8lf %4d %4d %10.3lf %4d %10.1lf %6.1lf %6.1lf %6.1lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf\n",
			Mjd.Days + Mjd.FracDay, SatAtod->id, SatAtod->Valid, XECF[0], XECF[1], XECF[2], XECF[3], XECF[4], XECF[5],
			XECI[0], XECI[1], XECI[2], XECI[3], XECI[4], XECI[5], Clk * 1E6, SatAtod->Clk[1] * 1E9,
			SatAtod->in_num + SatAtod->out_num, SatAtod->in_num + SatAtod->out_num, SatAtod->PDOP, T.Week, T.SecOfWeek,
			sqrt(SatAtod->CovX[0]), sqrt(SatAtod->CovX[7]), sqrt(SatAtod->CovX[14]), SatAtod->URA,
			SatAtod->MeanClk_apr, SatAtod->StdClk_apr, SatAtod->MeanClk_pst, SatAtod->StdClk_pst,
			SatAtod->MeanOrb_apr, SatAtod->StdOrb_apr, SatAtod->MeanOrb_pst, SatAtod->StdOrb_pst);
	}
}