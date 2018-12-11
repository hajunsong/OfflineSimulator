#pragma once

#include <vector>

using namespace std;

#define OptimalCPUTimeCheck 1
#define time_indx_MAX 10000

#pragma warning(disable:4996)

typedef struct RTTOutputData {
	vector <double> distance_by_time, velocity_by_time;
	double *velocity_command;
	int DataLength;
	int sim_count;
}ROUTPUTDATA;

typedef struct OptimalVelocityInputData {
	int NumberofThreads, WaypointSize;
	double **LocalPath, **VehicleVelocity, CurrentVelocity, step_size;	
	int **StabilityIndex;
}OptimalInputData;

class optimal_velocity_planning {
public:
	optimal_velocity_planning(int WaypointSize, int NumberofThreads);
	~optimal_velocity_planning();
	
	void init();
	void run(ROUTPUTDATA *ROutputData);

	OptimalInputData *OID;

	void set_indx(int indx) { data_indx = indx; }
	double get_opt_time() { return (opt_end - opt_start); }


public:

	double v_x;				// Waypoint 배열, 주행 속도
	int data_indx;					// map data의 순번 (시뮬레이션 순번)
	int WP_size;					// waypoint 개수
	double *input_distance;			// 누적거리
	double *h_in;					// 각 waypoint 사이의 거리(간격)
	double *v_in;					// 속도 명령 입력 (가감속 적용 되지 않은 raw data)	
	double *v_out;					// 속도 명령 출력 (가감속이 적용된 속도 명령)
	double opt_start, opt_end;		// CPU 계산 시간 측정을 위한 변수
	double start_time, end_time;	// 시뮬레이션 시작 및 종료 시간
	double h, g;					// step time 및 중력가속도
	double t_current;				// 시뮬레이션 현재 시간
									// define Y vector
	double Y[2], Yp[2];				// 적분기 사용을 위한 상태 변수
	int velocity_indx;				// v_in 배열의 순번
	double vel_by_time_size;		// 적분기를 사용하여 시간에 대한 속도 변수로 나타낼 때 사용되는 배열의 전체 크기

// 	double element_1[WP_data_point - 2], element_2[WP_data_point - 2], element_3[WP_data_point - 2], element_4[WP_data_point - 2], S_i[WP_data_point];
// 	double a_i[WP_data_point - 1], b_i[WP_data_point - 1], c_i[WP_data_point - 1], d_i[WP_data_point - 1];
	double *element_1, *element_2, *element_3, *element_4, *S_i, *a_i, *b_i, *c_i, *d_i;	// Cubic spline 계산에 사용되는 변수


// 	double distance_by_time[time_indx_MAX];	// 시간에 대한 거리
// 	double velocity_by_time[time_indx_MAX];	// 시간에 대한 속도

	// absh3
	int intcount, num_state;						// absh3 적분기에 사용되는 변수
	double AW[2][2], AW1[2][2], Y_next[2], t_next;	// absh3 적분기에 사용되는 변수

	void absh3_for_control(double step_size, const int n, int *intcount, double AW1[2][2], double AW[2][2], double *t_current, double Y[2], double Yp[2], double Y_next[2], double *t_next);
	void define_Y_vector_for_control(double* input_distance, double* input_velocity);
	void dqddq2Yp_for_control(double *input_distance);
	void solve_tridiagonal_in_place_destructive(double *x, int X, double * a, double * b, double * c);
};