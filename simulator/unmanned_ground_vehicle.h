#pragma once

#include "pnu_map_tire.h"
#include <vector>

using namespace std;

#define CPUTimeCheckMain 0

// suspension index
#define RF 0
#define RM 1
#define RR 2
#define LF 3
#define LM 4
#define LR 5

typedef struct GlobalInput {
	// Map Info
	int MapInfo_Row, MapInfo_Col, MapInfo_Robot_Row, MapInfo_Robot_Col;
	double MapInfo_Resolution_Row, MapInfo_Resolution_Col, MapInfo_Resolution_Down, MapInfo_ReferenceNorth, MapInfo_ReferenceEast, **ElevationData;
	// Navigation Info
	double Velocity_North, Velocity_East, Velocity_Down, Velocity, Attitude_Yaw;
	// Local Path Info
	int WaypointSize;
	double **LocalPath;

	vector <double> distance_by_time, velocity_by_time;
	double *velocity_command;
	int DataLength;
	int flag;
}GINPUT;

class unmanned_ground_vehicle {
public :
	unmanned_ground_vehicle();
	~unmanned_ground_vehicle();

	bool HilbertMap;
	double LocalPathSlopeAngle;
	int sim_count;
private :
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
		double r0cp[3], dr0cp[3], ddr0cp[3], w0p[3];
	}CHASSIS;

	typedef struct suspension_variable {
		// suspension id
		int id;
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
	}SUS;

	typedef struct simulation_variable {
		// state vector
		double Y[31], Yp[31];
		// absh3 - explicit integrator variables
		int intcount;
		double AW[31][2], AW1[31][2], Y_next[31], t_next;
		// system variable
		double h, g, t_current;
		// simulation control variables
		int simulation_flag, sim_indx, singular_flag, singular_flag2;
		// heading direction, current velocity
		double heading, vx;
		// equilibrium state 해석 후 결과 저장을 위한 state vector
		double Y_equil[13];
	}SIM;

	typedef struct controller_variable {
		// controller variable
		int WP_indx, WP_size, WP_indx_update_flag;	// parallel simulation 코드와 optimal velocity control 코드에서 사용	
		double yaw_d, e_v_sum, M_d_hat, yaw_hat, eq_prev;			// read control에서 사용				   
		double v_di, v_d, v_x, v_d_init;					// LP stability metric, save data, equilibrium, parallel process 등에서 사용
		double e_l, e_psi;						// LP stability metric에서 사용 (lateral position error, 지향각 오차)
		double F_xd_total, e_v;						// save data에서 사용
		double F_x_error[10][6], F_x_error_sum[6];		// F_xi(desired - actual)
		double motor_torque[6];	//output		
		double ddr0c_p[3];      // 차체의 로컬 가속도	, LP stability metric에서 사용		
		int RSM_indx, PSM_indx, LSM_indx, VSM_indx, LPE_indx;	// save data에서 사용
		double RSM, PSM, LSM, VSM, LPE;			// save data에서 사용

		double u_LQ_command[300];	// 경로추종제어기의 종방향 힘 입력
		double v_LQ_command[300];	// 경로추종제어기의 desired velocity 명령 입력
		double v_d_command[300];		// 경로추종제어기의 desired velocity 명령 입력
		double v_d_LQ;
		int sim_5_count_for_control;

		// v_d Low pass filtering parameters
		double v_d_before, tau_LPF, a_LPF;
	}CTRL;

	typedef struct _lq {	// LQ 속도 프로파일 처리를 위한 적분기(absh3_for_control 등)에 사용 (2016.11.11 홍효성)
		double start_time, end_time, h, g;
		double t_current;
		// define Y vector
		double Y[2], Yp[2];
		int velocity_index, vel_by_time_size;

#define LQ_data_point 10
#define time_indx_MAX_global 5000

		//double h_i[LQ_data_point - 1];
		double element_1[LQ_data_point - 2], element_2[LQ_data_point - 2], element_3[LQ_data_point - 2], element_4[LQ_data_point - 2], S_i[LQ_data_point];
		double a_i[LQ_data_point - 1], b_i[LQ_data_point - 1], c_i[LQ_data_point - 1], d_i[LQ_data_point - 1];

		//double input_distance[LQ_data_point];
		//double input_velocity[LQ_data_point];
		double distance_by_time[time_indx_MAX_global];	// 시간에 대한 거리
		//double velocity_by_time[time_indx_MAX_global];	// 시간에 대한 속도
		double *velocity_by_time;
		double LocalPathSlopeAngle;

		// absh3
		int intcount, num_state;
		double AW[2][2], AW1[2][2], Y_next[2], t_next;

		double *h_in, *input_distance, *v_out;
	}LQ;

	SUS *sus;
	CHASSIS *chassis;
	SIM *sim;
	CTRL *ctrl;
	Obj_Tire_Global *tire;
	LQ *lq;

	FILE *fp;
	char buf[256];

	int WP_size, RTT_WP_size;
	bool RTT_run_flag, vel_command_flag;
	double heading;
	double **LocalPath, *VehicleVelocity, MaximumVelocity, CurrentVelocity;
	int *StabilityIndex;
	double *VelocityCommand;

public :
	void init(GINPUT *GInput);
	void run();

	int run_RTT(int sim_indx);
	
private :
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

	double Find_Max(double Array[], int length);
	double Find_min(double data, double min_val);

	double fsign(double data);
	int ludcmp4(double a[4][4], int n, int indx[4], double d, double fac[4][4]);
	void lubksb4(double a_fac[4][4], int n, int indx[4], double b[4], double x[4]);

	void equilibrium_analysis();

	void read_system();
	void read_chassis();
	void read_RF(SUS *sus, CHASSIS *chassis);
	void read_RM(SUS *sus, CHASSIS *chassis);
	void read_RR(SUS *sus, CHASSIS *chassis);
	void read_LF(SUS *sus, CHASSIS *chassis);
	void read_LM(SUS *sus, CHASSIS *chassis);
	void read_LR(SUS *sus, CHASSIS *chassis);
	void read_control();

	void define_Y_vector();

	void analysis();

	void Y2qdq();

	void orientation_chassis();
	void position_chassis();
	void velocity_state_chassis();
	void cartesian_velocity_chassis();
	void mass_force_state_chassis();

	void orientation_suspension(SUS *sus);
	void position_suspension(SUS *sus);
	void velocity_state_suspension(SUS *sus);
	void cartesian_velocity_suspension(SUS *sus);

	void pre_road(SUS sus[6]);

	void mass_force_state_suspension(SUS *sus);
		void CNU_TSDA(SUS *sus);
	void velocity_coupling(SUS *sus);
	void effective_mass_force(SUS *sus);

	void acceleration_state_chassis();
	void cartesian_acceleration_chassis();

	void acceleration_state_suspension(SUS *sus);
	void cartesian_acceleration_suspension(SUS *sus);

	void wheel_spin_dyn();

	void LP_control();
	void LQ_velocity_control(int data_indx, int sim_count);
	void LP_stability_metric();
// 	void fuzzification(double fuzzy[][2], double e, double de);
// 	void fuzzy_inference(vector<vector<double>> &u_indx, double fuzzy[][2]);
// 	double defuzzification(vector<vector<double>> &u_indx);
	double Find_max2(double data, double max_val);

	void dqddq2Yp();

	void save_data(int sim_indx);

	void absh3(double step_size, const int n, int *intcount, double AW1[31][2], double AW[31][2], double *t_current, double Y[31], double Yp[31], double Y_next[31], double *t_next);

	void solve_tridiagonal_in_place_destructive(double *x, int X, double * a, double * b, double * c);
	void absh3_for_control(double step_size, const int n, int *intcount, double AW1[2][2], double AW[2][2], double *t_current, double Y[2], double Yp[2], double Y_next[2], double *t_next);
	void define_Y_vector_for_control(LQ *lq, double* input_distance, double* input_velocity);
	void dqddq2Yp_for_control(LQ *lq, double* input_distance);

	double gaussianRandom(double m, double sv);
};