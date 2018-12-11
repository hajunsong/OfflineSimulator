#pragma once

#include <memory.h>

#define DynamicCPUTimeCheck 1

#pragma warning(disable:4996)

/************************************************************************/
/*                              참고                                    */
/************************************************************************/
// 각 라이브러리 별로, 설명, lib 링크, include, namespace 지정 순으로 나열
// #include ""의 주소는 소스코드위치 기준
// #pragma comment의 주소는 프로젝트위치 기준

/************************************************************************/
/*                              std                                     */
/************************************************************************/
//#include <stdio.h>

typedef unsigned int uint;

#include <string>		
// using std::string;

#include <iostream>	
// using std::cout; using std::endl;

#include <vector>		
// using std::vector; 

#include <deque>		
// using std::deque;

namespace L
{
	/************************************************************************/
	/*                           전역변수                                   */
	/************************************************************************/

	// 	const uword NULL_uword = uword(-1); //matlab에서는 0으로 초기화했는데, c++에서는 0도 index로 의미가 있어서 이렇게 처리, uword(-1)은 자료형의 최대치
	// 	const double PSV=1e-6; // properly small value적당히 작은값
	// 	const double pi = datum::pi;
	// 	const double eps = datum::eps;


	/************************************************************************/
	/*                         템플릿 함수                                  */
	/************************************************************************/

	// 	template <typename T>
	// 	T fn_sign(T x)
	// 	{
	// 		T sol;
	// 		if		(x>0)	{sol=+1;}
	// 		else if (x<0)	{sol=-1;}
	// 		else			{sol=0;}
	// 		return sol; 
	// 	}


	/************************************************************************/
	/*                         인라인 함수                                  */
	/************************************************************************/

	inline void fn_print2DArray(double **matrix, uint rows, uint cols)
	{
		uint x = 0, y = 0;
		puts("");
		for (; x < rows; x++) {
			for (y = 0; y < cols; y++) {
				printf("%lf\t", matrix[x][y]);
			}
			puts("");
		}
	}

	inline double** fn_get2DArray(uint rows, uint cols)
	{
		uint index = 0;
		double **array = (double**)calloc(rows, sizeof(double*));

		for (; index < rows; index++) {
			array[index] = (double*)calloc(cols, sizeof(double));
		}

		return array;
	}

	inline void fn_freeMatrix(double **matrix, uint rows)
	{
		rows--;
		for (; rows > -1; rows--) {
			free(matrix[rows]);
		}

		free(matrix);
		//matrix = NULL;
	}

	inline void fn_file2DArray(double** matrix, uint rows, uint cols, char *file_name)
	{
		FILE *file;
		int err = 0;
		uint x = 0, y = 0;

		if ((err = fopen_s(&file, file_name, "rt")) != 0)
		{
			printf_s("\n File open error! \n");
		}

		for (x = 0; x < rows; x++) {
			for (y = 0; y < cols; y++) {
				fscanf_s(file, "%lf", &matrix[x][y]);
			}
		}

		fclose(file);
	}

	inline void fn_Array2Dfile(double** matrix, uint rows, uint cols, char *file_name)
	{
		FILE *file;
		int err = 0;
		uint x = 0, y = 0;

		if ((err = fopen_s(&file, file_name, "wt")) != 0)
		{
			printf_s("\n File open error! \n");
		}

		for (x = 0; x < rows; x++) {
			for (y = 0; y < cols; y++) {
				fprintf_s(file, "%lf\t", matrix[x][y]);
			}
			fprintf_s(file, "\n");
		}

		fclose(file);
	}
}

class Obj_Map
{

private:
	static const double*const*_map_data;
	static uint _num_x; //포인터로 수정
	static uint _num_y;
	static double _gab;
	static double _x_c;
	static double _y_c;
	static double *_x_data; //갯수가 런타임 결정이기 때문에 배열 사용곤란
	static double *_y_data;

public:
	~Obj_Map();
	static void PNU_set_map(double **map_data, uint rows, uint cols, double gab, double x_c, double y_c);
	static void PNU_fn_map(double x_p, double y_p, double *z);

private:
	static void fn_linspace(double *data, double L, unsigned int num_data, double c);
	static void fn_find_patch_uv(uint* patch_r, uint* patch_c, double* u, double* v, double x_p, double y_p);
	static void fn_MBMt_uv(double(*MBMt)[4], double* u, double* v, double x_p, double y_p);
	static void fn_patch_idxs(int row_idx[4], int col_idx[4], uint patch_r, uint patch_c);
	static void fn_B_mat(double(*B)[4], uint patch_r, uint patch_c);
	static void fn_mat444_multi(double(*C)[4], double(*A)[4], double(*B)[4]);
	static void fn_mat1441_multi(double *A, double(*B)[4], double *C, double *output);
};

class Obj_Nan
{
private:
	static double z_min; //포인터로 수정

public:
	static void fn_nan_interp(double **Map_data, uint n_rows, uint n_cols);

private:
	static double fn_nan_patch(double **Map_data, uint n_rows, uint n_cols, uint i, uint j);
	static int is_nan(double z);
	static void fn_ij_p(int* i_p, int* j_p, uint p, uint k, uint i, uint j);
};

class Obj_Tire
{
private:
	//Obj_Map* _C_Map;

	double _c_t, _k_t, _Iyy, _mus, _mud, _w, _rr, _Ca, _CSLIP, _R_u, _t_s, _slope_rr;
	double _eps, _pi;
	double _c_x, _c_y;
	double _slip_x, _slip_y;

	double _Fx, _Fy, _Fz, _R_e, _R_d, _My; // 출력을 위한 저장변수
	double _F_tire[3], _M_tire[3], _F_global[3], _M_global[3];
	double _slip, _angle;
	double _pen, _road_h; // 디버깅을 위한 임시 변수

	double _Road_h_i[6], _Fric_i[6], _N_vec[3];

	double _P[8]; //타이어 파라미터

	double _mu;

public:
	Obj_Tire(); //생성자
	void PNU_Pre_road(double R_c_RF[3], double R_c_RM[3], double R_c_RR[3], double R_c_LF[3], double R_c_LM[3], double R_c_LR[3]);
	void PNU_Tire_force(double R_c[3], double dR_c[3], double u_vec[3], double omega, int index);
	void PNU_get_data(double F[3], double M[3], double *slip, double *angle, double *road_h, double *pen, double *R_d, double *Fx, double *Fy, double *Fz, double *My);
	double PNU_get_unloaded_radius() { return _R_u; };
	void set_Road_h(double road[6]) { memcpy(_Road_h_i, road, sizeof(double) * 6); }
	void set_N_vec(double vec[3]) { memcpy(_N_vec, vec, sizeof(double) * 3); }
	double PNU_get_mu() { return _mu; }
	void set_mu(double road_mu) { _mu = road_mu; };

private:
	void fn_R_eff(double *R_e, double *a, double R_d);
	void fn_road_tire_A(double road_A[3][3], double tire_A[3][3], double road_z_Vec[3], double tire_y_Vec[3]);
	//void fn_slip_ratio(double *slip, double *angle, double V_x, double omega, double R_e, double V_sy, double Fz);
	void fn_slip_ratio(double *slip_x, double *slip_y, double V_x, double omega, double R_e, double V_sy, double Fz, double a);
	void fn_Fiala_tire(double F[3], double M[3], double Fz, double slip, double alpha, double omega);
	void fn_Brush_tire(double F[3], double M[3], double Fz, double slip_x, double slip_y, double mu, double a, double pen, double omega);
	void fn_step_s(double x, double slope, double *output);
	void fn_cross(double c[3], double a[3], double b[3]);
	void fn_mat3331(double c[3], double a[3][3], double b[3]);
	void fn_mat33T31(double c[3], double a[3][3], double b[3]);
	void fn_sign(double a, double *output);
	void fn_max(double a, double b, double *output);
	void fn_min(double a, double b, double *output);
	void fn_normalize(double a[], uint num);
};

// suspension index
#define RF 0
#define RM 1
#define RR 2
#define LF 3
#define LM 4
#define LR 5

typedef struct RTTInputData {
	// Map Info
	int MapInfo_Row, MapInfo_Col, MapInfo_Robot_Row, MapInfo_Robot_Col;
	double MapInfo_Resolution_Row, MapInfo_Resolution_Col, MapInfo_Resolution_Down, MapInfo_ReferenceNorth, MapInfo_ReferenceEast, **ElevationData, Map_Time;
	// Navigation Info
	double Position_East, Position_North, Velocity_Longitude, Velocity_North, Velocity_East, Velocity_Down, Velocity, Attitude_Yaw, Navi_Time;
	// Local Path Info
	int WaypointSize;
	double **LocalPath, Path_Time;
	// 고도맵과 경로 데이터가 맞지 않은 경우 에러 메시지
	char input_data_err[255];
	int err_flag;
	// 고도맵 마찰 계수
	double road_mu;

	int sim_count;
}RINPUTDATA;

typedef struct RealtimeDynamicInputData {
	int NumberofThreads, WaypointSize;
	double **LocalPath, **VehicleVelocity, CurrentVelocity, step_size;
	int **StabilityIndex;
	int err_flag;
	double LocalPathSlopeAngle;
	double vd_mission;
	double **RSM, **PSM, **LSM, **VSM;
}DynamicInputData;

class realtime_dynamics_analysis {
public:
	realtime_dynamics_analysis(int WaypointSize);
	~realtime_dynamics_analysis();

	void init(RINPUTDATA *RInputData);
	void LPP(RINPUTDATA *RTT_input);
	void run(RINPUTDATA *RInputData);	

	DynamicInputData *DID;

	double get_read_time() { return read_time; }
	double get_pre_time() { return equilibrium_time; }
	double get_parallel_time() { return parallel_time; }
	void set_indx(int simulation_indx) { sim_indx = simulation_indx; }

private:
	typedef struct suspension_variable {
		int id;	// suspension id
		// read suspension
		double q1, dq1, rho1p[3], C11[3][3], m1, J1[3][3], s01p[3], s12p[3], C01[3][3], C12[3][3], s0sp[3], s1sp[3];
		// orientation suspension
		double A1[3][3], H1[3];
		// position suspension
		double r1[3], s12[3], rw[3], rho1[3], r1c[3];
		// velocity state suspension
		double r1t[3][3], B1[6], Yh1[6];
		// cartesian velocity suspension
		double T1[6][6], Yb1[6], dr1[3], w1[3], w1t[3][3], dr1c[3], drwc[3];
		// mass force state suspension
		double J1c[3][3], r1ct[3][3], dr1ct[3][3], Mh1[6][6], F1c[3], s12t[3][3], rho1t[3][3], T1c[3], Qh1[6], F_tire[3], M_tire[3];
		// CNU TSDA
		double d01[3], L_spring, defo, T_spring, T_damper, Qh0_TSDA[6], Qh1_TSDA[6], f;
		// velocity coupling
		double dr1t[3][3], D1[6], dH1[3];
		// effective mass force
		double Myq[6], Pq, inv_Mqq, Mhc[6][6], Phc[6];
		// acceleration state suspension
		double dYh1[6], ddq1;
		// cartesian acceleration suspension
		double ddr1[3], dw1[3], dw1t[3][3], ddr1c[3];
		// wheel & tire
		double theta_wh, w_wh, dw_wh, T_in;
		// Tire 계산 후 결과 저장 변수 - slip ratio, slip angle, Tire longitudinal force(Local), Tire lateral force(Local), Tire vertical force(Local), deformed radius, penetration, height of road, Moment y axis(타이어 회전 축, local)
		double slip, angle, Fx, Fy, Fz, R_d, pen, road_h, My;

		double pen_old;	// 이전 time step 타이어 침투량
		double pen_old2; // 두 스텝 전 타이어 침투량
	}SUS;

	typedef struct chassis_variable {
		// read chassis
		double r0[3], e0[4], dr0[3], w0[3], dr0p[3], rho0p[3], C00[3][3], m0, J0[3][3];
		// orientation chassis
		double A0[3][3], E0[3][4], G0[3][4], roll_ang, pitch_ang, yaw_ang;
		// position chassis
		double rho0[3], r0c[3];
		// velocity state chassis
		double r0t[3][3], Yh0[6];
		// cartesian velocity chassis
		double de0[4], T0[6][6], w0t[3][3], dr0c[3], dr0t[3][3];
		// mass force state chassis
		double J0c[3][3], r0ct[3][3], dr0ct[3][3], Mh0[6][6], F0c[3], T0c[3], Qh0[6];
		// acceleration state chassis
		double dYh0[6];
		// cartesian acceleration chassis
		double ddr0[3], dw0[3], dw0t[3][3], ddr0c[3];
		// local variable
		double dr0cp[3], ddr0cp[3], w0p[3];
	}CHASSIS;

	typedef struct simulation_variable {
		// state vector
		double Y[31], Yp[31];
		// absh3 - explicit integrator variables
		int intcount;
		double AW[31][2], AW1[31][2], Y_next[31], t_next;
		// system variable
		double h, g, t_current;
		// simulation control variables
		int simulation_flag, sim_indx, thread_indx, singular_flag, singular_flag2;
		// heading direction, current velocity
		double heading, vx;
		// equilibrium state 해석 후 결과 저장을 위한 state vector
		double Y_equil[13];

		// 한 스텝 전 노면, 두 스텝 전 노면
		double pre_road_h[6], old_road_h[6];
		// 노면 에러 알고리즘을 적용한 횟수, 기준 값 이상이면 해석을 종료함
		int peak_error_count, peak_error_count_max;
		// 차량의 전복 상태를 저장하는 flag 변수
		int overturn_flag;
		int main_count;

		FILE *fp[50];
		char buf[50][256];
	}SIM;

	typedef struct controller_variable {
		int WP_indx, WP_size, WP_indx_update_flag;	// parallel simulation 코드와 optimal velocity control 코드에서 사용	
		double yaw_d, e_v_sum, M_d_hat, yaw_hat;			// read control에서 사용				   
		double v_di, v_d[7], v_x;					// LP stability metric, save data, equilibrium, parallel process 등에서 사용
		double e_l, e_psi;						// LP stability metric에서 사용 (lateral position error, 지향각 오차)
		double F_xd_total, e_v;						// save data에서 사용
		double motor_torque[6];	//output		
		double ddr0c_p[3];      // 차체의 로컬 가속도	, LP stability metric에서 사용		
		int RSM_indx, PSM_indx, LSM_indx, VSM_indx, LPE_indx;	// save data에서 사용
		double RSM, PSM, LSM, VSM, LPE;			// save data에서 사용
	}CTRL;

	int sim_indx; // simulation index(6000개 데이터 중 몇번째 데이터를 시뮬레이션 하는지 구분하기 위한 변수)
	double heading, vx; // heading direction, current velocity
	int WP_size; // way point size
	int NumberofThreads; // Number of Threads
	double road_mu; // 고도맵 마찰계수

	// CPU time 측정 변수, 실제 차량 적용시 불필요
	double read_start1, read_end1, read_start2, read_end2, read_time;
	double equilibrium_start, equilibrium_end, equilibrium_time;
	double parallel_start, parallel_end, parallel_time;

	void tilde(double a[3], double b[3][3]);
	void mat33T31(double a[3][3], double b[3], double c[3]);
	void mat3333(double a[3][3], double b[3][3], double c[3][3]);
	void mat3331(double a[3][3], double b[3], double c[3]);
	void mat333333(double a[3][3], double b[3][3], double c[3][3], double d[3][3]);
	void mat34T31(double a[3][4], double b[3], double c[4]);
	void mat3333T(double a[3][3], double b[3][3], double c[3][3]);
	void mat333331(double a[3][3], double b[3][3], double c[3], double d[3]);
	void mat6661(double a[6][6], double b[6], double c[6]);
	int ludcmp6(double a[6][6], int n, int indx[6], double d, double a_fac[6][6]);
	void lubksb6(double a_fac[6][6], int n, int indx[6], double b[6], double x[6]);

	//	void fMulti(double(*M)[6], double(*A)[6], double(*B)[6], int row1, int column1, int row2, int column2);
	//	void fSum(double(*M)[6], double(*A)[6], double(*B)[6], int row1, int column1, int row2, int column2);
	//	void fInverse22(double(*M)[6], double(*A)[6]);
	//	void fTranspose(double(*Result)[6], double(*A)[6]);
	//	void fInverse44(double(*M)[6], double(*A)[6], double t_current, int thread_indx);
	double Find_Max(double Array[], int length);
	double Find_min(double data, double min_val);
	//	void solve_tridiagonal_in_place_destructive(double *x, int X, double * a, double * b, double * c);
	double fsign(double data);
	int ludcmp4(double a[4][4], int n, int indx[4], double d, double fac[4][4]);
	void lusolve4(double a_fac[4][4], int n, int indx[4], double b[4], double x[4]);

	void read_system(SIM *sim);
	void read_chassis(SIM *sim, CHASSIS *chassis, SUS sus[6], double R_u);
	void read_RF(SUS *sus, CHASSIS *chassis);
	void read_RM(SUS *sus, CHASSIS *chassis);
	void read_RR(SUS *sus, CHASSIS *chassis);
	void read_LF(SUS *sus, CHASSIS *chassis);
	void read_LM(SUS *sus, CHASSIS *chassis);
	void read_LR(SUS *sus, CHASSIS *chassis);
	void read_control(CTRL *ctrl);

	void define_Y_vector(SIM *sim, CHASSIS *chassis, SUS sus[6]);

	void equilibrium_process(SIM *sim, CHASSIS *chassis, SUS sus[6], CTRL *ctrl, Obj_Tire *tire);
		void velocity_candidate(CTRL *ctrl, SIM *sim, double mu);

	void parallel_process(SIM *sim, CHASSIS *chassis, SUS sus[6], CTRL *ctrl, Obj_Tire *tire);

	void analysis(SIM *sim, CHASSIS *chassis, SUS sus[6], CTRL *ctrl, Obj_Tire *tire);

	void Y2qdq(SIM *sim, CHASSIS *chassis, SUS sus[6]);

	void orientation_chassis(CHASSIS *chassis);
	void position_chassis(CHASSIS *chassis);
	void velocity_state_chassis(CHASSIS *chassis);
	void cartesian_velocity_chassis(CHASSIS *chassis);
	void mass_force_state_chassis(CHASSIS *chassis, SIM *sim);

	void orientation_suspension(SUS *sus, CHASSIS *chassis);
	void position_suspension(SUS *sus, CHASSIS *chassis);
	void velocity_state_suspension(SUS *sus, CHASSIS *chassis);
	void cartesian_velocity_suspension(SUS *sus);

	void pre_road(SUS sus[6], SIM *sim, Obj_Tire *tire);
		int peak_error_detect(double current_road[6], double previous_road[6], double old_road[6]);

	void mass_force_state_suspension(SUS *sus, SIM *sim, Obj_Tire *tire, CHASSIS *chassis);
		void CNU_TSDA(SUS *sus, CHASSIS *chassis, SIM *sim);
	void velocity_coupling(CHASSIS *chassis, SUS *sus);
	void effective_mass_force(SUS *sus);

	int overturn_detect(CHASSIS *chassis, SUS sus[6]);

	void acceleration_state_chassis(CHASSIS *chassis, SUS sus[6], SIM *sim);
	void cartesian_acceleration_chassis(CHASSIS *chassis);

	void acceleration_state_suspension(CHASSIS *chassis, SUS *sus);
	void cartesian_acceleration_suspension(SUS *sus);

	void wheel_spin_dyn(SUS *sus, CTRL *ctrl);

	void LP_control(CHASSIS *chassis, SUS sus[6], CTRL *ctrl, SIM *sim);
		void LP_stability_metric(CHASSIS *chassis, SUS sus[6], SIM *sim, CTRL *ctrl, int thread_indx);

	void dqddq2Yp(SIM *sim, CHASSIS *chassis, SUS sus[6]);

	void save_data(SIM *sim, CHASSIS *chassis, SUS sus[6], CTRL *ctrl);

	void absh3(double step_size, const int n, int *intcount, double AW1[31][2], double AW[31][2], double *t_current, double Y[31], double Yp[31], double Y_next[31], double *t_next);
};