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

	double v_x;				// Waypoint �迭, ���� �ӵ�
	int data_indx;					// map data�� ���� (�ùķ��̼� ����)
	int WP_size;					// waypoint ����
	double *input_distance;			// �����Ÿ�
	double *h_in;					// �� waypoint ������ �Ÿ�(����)
	double *v_in;					// �ӵ� ��� �Է� (������ ���� ���� ���� raw data)	
	double *v_out;					// �ӵ� ��� ��� (�������� ����� �ӵ� ���)
	double opt_start, opt_end;		// CPU ��� �ð� ������ ���� ����
	double start_time, end_time;	// �ùķ��̼� ���� �� ���� �ð�
	double h, g;					// step time �� �߷°��ӵ�
	double t_current;				// �ùķ��̼� ���� �ð�
									// define Y vector
	double Y[2], Yp[2];				// ���б� ����� ���� ���� ����
	int velocity_indx;				// v_in �迭�� ����
	double vel_by_time_size;		// ���б⸦ ����Ͽ� �ð��� ���� �ӵ� ������ ��Ÿ�� �� ���Ǵ� �迭�� ��ü ũ��

// 	double element_1[WP_data_point - 2], element_2[WP_data_point - 2], element_3[WP_data_point - 2], element_4[WP_data_point - 2], S_i[WP_data_point];
// 	double a_i[WP_data_point - 1], b_i[WP_data_point - 1], c_i[WP_data_point - 1], d_i[WP_data_point - 1];
	double *element_1, *element_2, *element_3, *element_4, *S_i, *a_i, *b_i, *c_i, *d_i;	// Cubic spline ��꿡 ���Ǵ� ����


// 	double distance_by_time[time_indx_MAX];	// �ð��� ���� �Ÿ�
// 	double velocity_by_time[time_indx_MAX];	// �ð��� ���� �ӵ�

	// absh3
	int intcount, num_state;						// absh3 ���б⿡ ���Ǵ� ����
	double AW[2][2], AW1[2][2], Y_next[2], t_next;	// absh3 ���б⿡ ���Ǵ� ����

	void absh3_for_control(double step_size, const int n, int *intcount, double AW1[2][2], double AW[2][2], double *t_current, double Y[2], double Yp[2], double Y_next[2], double *t_next);
	void define_Y_vector_for_control(double* input_distance, double* input_velocity);
	void dqddq2Yp_for_control(double *input_distance);
	void solve_tridiagonal_in_place_destructive(double *x, int X, double * a, double * b, double * c);
};