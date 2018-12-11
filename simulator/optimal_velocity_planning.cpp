#include "optimal_velocity_planning.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <memory.h>
#include <omp.h>

optimal_velocity_planning::optimal_velocity_planning(int WaypointSize, int NumberofThreads) {
	OID = new OptimalInputData;
	OID->NumberofThreads = NumberofThreads;
	OID->WaypointSize = WaypointSize;
	OID->LocalPath = new double*[WaypointSize];
	OID->StabilityIndex = new int*[WaypointSize];
	OID->VehicleVelocity = new double*[WaypointSize];
	for (int i = 0; i < WaypointSize; i++) {
		OID->LocalPath[i] = new double[3];
		OID->StabilityIndex[i] = new int[NumberofThreads];
		OID->VehicleVelocity[i] = new double[NumberofThreads];
	}
	
	v_x = 0;
	data_indx = 0;
	opt_start = 0;
	opt_end = 0;
}

optimal_velocity_planning::~optimal_velocity_planning() {
	delete[] input_distance, v_in, v_out, h_in;
	delete[] element_1, element_2, element_3, element_4, S_i, a_i, b_i, c_i, d_i;

	for (int i = 0; i < WP_size; i++) {
		delete[] OID->LocalPath[i], OID->StabilityIndex[i], OID->VehicleVelocity[i];
	}
	delete[] OID->LocalPath, OID->StabilityIndex, OID->VehicleVelocity;
	delete OID;
}

void optimal_velocity_planning::init() {
	WP_size = OID->WaypointSize;		// waypoint ���� �ڵ� �Է�
	
	input_distance = new double[WP_size];	// input_distance �迭 �޸� �Ҵ� input_distance[WP_size]
	v_in = new double[WP_size];				// v_in �迭 �޸� �Ҵ� v_in[WP_size]
	v_out = new double[WP_size];			// v_out �迭 �޸� �Ҵ� v_in[WP_size]
	h_in = new double[WP_size - 1];			// h_in �迭 �޸� �Ҵ� v_in[WP_size - 1]	

	v_x = OID->CurrentVelocity;			// ������ ���� ���� �ӵ� �Է�

	element_1 = new double[WP_size - 2];	// Cubic spline Ǯ�̿� ����ϱ� ���� �迭 �޸� �Ҵ� (���� �밢 ��Ŀ����� ���� �밢 ���)
	element_2 = new double[WP_size - 2];	// Cubic spline Ǯ�̿� ����ϱ� ���� �迭 �޸� �Ҵ� (���� �밢 ��Ŀ����� ���� �밢 ���)
	element_3 = new double[WP_size - 2];	// Cubic spline Ǯ�̿� ����ϱ� ���� �迭 �޸� �Ҵ� (���� �밢 ��Ŀ����� ���� �밢 ���)
	element_4 = new double[WP_size - 2];	// Cubic spline Ǯ�̿� ����ϱ� ���� �迭 �޸� �Ҵ� (�캯 ����)
	S_i = new double[WP_size];				// Cubic spline Ǯ�̿� ����ϱ� ���� �迭 �޸� �Ҵ� (���ӵ�)
	a_i = new double[WP_size - 1];			// Cubic spline Ǯ�̿� ����ϱ� ���� �迭 �޸� �Ҵ� (���׽��� ���)
	b_i = new double[WP_size - 1];			// Cubic spline Ǯ�̿� ����ϱ� ���� �迭 �޸� �Ҵ� (���׽��� ���)
	c_i = new double[WP_size - 1];			// Cubic spline Ǯ�̿� ����ϱ� ���� �迭 �޸� �Ҵ� (���׽��� ���)
	d_i = new double[WP_size - 1];			// Cubic spline Ǯ�̿� ����ϱ� ���� �迭 �޸� �Ҵ� (���׽��� ���)

	h = OID->step_size;	// ���� ������. ������ �𵨿��� �Է����� �޾ƿ;� ������ �ӽ÷� �׳� ���, ���� ���� �� ��(2017-02-28 song hajun)
}

void optimal_velocity_planning::run(ROUTPUTDATA *RTT_output) {

	////////// �Է� �޾ƾ��� ���� �� �Ķ���� //////////
	double a_acc_max = 4;			// �ִ� ���ӵ�
	double a_dec_max = -2;			// �ִ� ���ӵ�
	double lowest_velocity_command = 1.0;	// ��� thread�� ���� �������� fail �Ǿ��� �� ������ ���� �ӵ�

											////////// Sectoin 0. ���� ���� ���� //////////
	int WP_thread, WP_indx;	// CPU ������ ��ȣ, waypoint�� ����(�ε���)
	int i, k;					// for �ݺ��� ����� ���� ����
	int time_indx = 1;			// �ð��� ���� �ӵ��� ��� ������ ������ �������� �������Ƿ� ������ �ε����� ����Ѵ�. time_indx = 1���� ���� ���е� ������ ���� (2016.11.18 ȫȿ��)
	int check_flag = 1;			// ������ �����Ͽ� �ӵ� ����ȭ �� �� �ӵ��� ������ �޶������� ���θ� Ȯ���ϱ� ���� flag
	double v_acc, v_dec;		// �ִ� ���ӵ��� ����� �ӵ�, �ִ� ���ӵ��� ����� �ӵ�	
	double WP_previous[4], WP_next[4];		// ������ ���� �Ҽӵ� waypoint �� �� ���� waypoint���� ��ġ, �ӵ�, ���༺��ǥ��� ������ ���� �迭

	double *v_in_temp = new double[WP_size];

	opt_start = omp_get_wtime();		// CPU time ���� ����

	////////// Section 1. ���� ������ ��ǥ�� �������� waypoint �� �ִ� ���� ��� �ӵ� ��� //////////

// 	// 7���� thread�鿡 �ִ� �ӵ� �߿��� stability indx�� pass ó���� ���� ���� �ӵ��� �ش� WP�� �ӵ��� ����
 	for (WP_indx = 0; WP_indx < WP_size; WP_indx++) {
 		for (WP_thread = 0; WP_thread < 7; WP_thread++) {	// thread ������ �ӵ��� ���� ��(0�� thread)���� ���� ��(6�� thread)����
 			//if (WP[WP_indx][3][WP_thread] == 0) {			// ���� ���� �ӵ� �ĺ������� stability indx�� 0 (fail)�� �κ� ã��
 			if(OID->StabilityIndex[WP_indx][WP_thread] == 0) {
 				if (WP_thread == 0) {
 					//WP[WP_indx][2][7] = lowest_velocity_command;		// ���� �ӵ� thread���� ���� ������ ��ǥ fail�� ��� ���� �ӵ� ���� �Է� (2016.12.22)
 					v_in_temp[WP_indx] = lowest_velocity_command;
 				}
 				else {
 					//WP[WP_indx][2][7] = WP[WP_indx][2][WP_thread - 1];	// ��ǥ�� fail �Ǳ� ������ thread�� �ӵ� �ĺ����� �ش� WP������ �ִ� ��� �ӵ��� �Ҵ�
 					v_in_temp[WP_indx] = OID->VehicleVelocity[WP_indx][WP_thread - 1];
 				}
 				break;
 			}
 			//WP[WP_indx][2][7] = WP[WP_indx][2][WP_thread];	// ��� thread�� ��ǥ�� ������ ��� �ִ� �ӵ��� 6�� thread�� �ӵ��� �Ҵ�
 			v_in_temp[WP_indx] = OID->VehicleVelocity[WP_indx][WP_thread];
 		}
 	}
// 7���� thread�鿡 �ִ� �ӵ� �߿��� stability index�� pass ó���� ���� ���� �ӵ��� �ش� WP�� �ӵ��� ����
//	for (WP_indx = 0; WP_indx < WP_size; WP_indx++) {
//		for (WP_thread = 6; WP_thread >= 0; WP_thread--) {	// thread ������ �ӵ��� ���� ��(6�� thread)���� ���� ��(0�� thread)����
//// 			if (WP[WP_indx][3][WP_thread] == 1) {			// ���� ���� �ӵ� �ĺ������� stability index�� 1 (pass)�� �κ� ã��
//// 				WP[WP_indx][2][7] = WP[WP_indx][2][WP_thread];	// ��ǥ�� fail �Ǳ� ������ thread�� �ӵ� �ĺ����� �ش� WP������ �ִ� ��� �ӵ��� �Ҵ�
//// 				break;
//// 			}
//// 			WP[WP_indx][2][7] = lowest_velocity_command;	// ���� �ӵ� thread���� ���� ������ ��ǥ fail�� ��� ���� �ӵ� 1.5 ���� �Է� (2016.12.22)
//			if (OID->StabilityIndex[WP_indx][WP_thread] == 1) {	// ���� ���� �ӵ� �ĺ������� stability index�� 1 (pass)�� �κ� ã��
//				v_in_temp[WP_indx] = OID->VehicleVelocity[WP_indx][WP_thread];	// ��ǥ�� fail �Ǳ� ������ thread�� �ӵ� �ĺ����� �ش� WP������ �ִ� ��� �ӵ��� �Ҵ�
//				break;
//			}
//			v_in_temp[WP_indx] = lowest_velocity_command;	// ���� �ӵ� thread���� ���� ������ ��ǥ fail�� ��� ���� �ӵ� 1.5 ���� �Է� (2016.12.22)
//		}
//	}

	//WP[WP_size - 1][2][7] = WP[WP_size - 2][2][7];			// ������ waypoint�� �ӵ� �ĺ����� ���� waypoint�� �ӵ� �ĺ����� ������ ����
	v_in_temp[WP_size - 1] = v_in_temp[WP_size - 2];


	////////// Section 2. ������ �Ѱ踦 ����� �ӵ� ��� ��� //////////

	// ��ġ(x��)�� �ӵ�(y��) �Է�
	for (WP_indx = 0; WP_indx < WP_size; WP_indx++) {
		if (WP_indx == 0) {
			input_distance[0] = 0;		// �ʱ� ��ġ 0
			if (v_x < OID->VehicleVelocity[0][0])
				v_in[0] = OID->VehicleVelocity[0][0];
			//else if (v_x > OID->VehicleVelocity[0][6])
			//	v_in[0] = OID->VehicleVelocity[0][6];
			else
				v_in[0] = v_x;			
		}
		else {
// 			for (i = 0; i < 3; i++) {
// 				WP_previous[i] = WP[WP_indx - 1][i][7];				// waypoint 1 (x, y, velocity)
// 				WP_next[i] = WP[WP_indx][i][7];					// waypoint 2 (x, y, velocity)
// 			}
			WP_previous[0] = OID->LocalPath[WP_indx - 1][0];
			WP_previous[1] = OID->LocalPath[WP_indx - 1][1];
			WP_previous[2] = v_in_temp[WP_indx - 1];
			WP_next[0] = OID->LocalPath[WP_indx][0];
			WP_next[1] = OID->LocalPath[WP_indx][1];
			WP_next[2] = v_in_temp[WP_indx];

			input_distance[WP_indx] = input_distance[WP_indx - 1] + sqrt((WP_next[0] - WP_previous[0])*(WP_next[0] - WP_previous[0]) + (WP_next[1] - WP_previous[1])*(WP_next[1] - WP_previous[1]));	// ���� �Ÿ�
			v_in[WP_indx] = WP_next[2];	// �ش� waypoint������ �ִ� ���� ��� �ӵ� �Է�
		}
	}

	// �� waypoint ������ �Ÿ� ���� ���
	for (i = 0; i < WP_size - 1; i++) {
		h_in[i] = input_distance[i + 1] - input_distance[i];		// h_i = x(i + 1) - x(i); x��(���⼭�� �Ÿ�) ����
	}

	//v_in[0] = v_in_temp[0];		// �ʱ� �ӵ� ����� ����ӵ��� �����ϴ� �ڵ� ����

	for (i = 0; i < WP_size; i++) {
		v_out[i] = v_in[i];          // ���� v_out�� v_in�� ����ȭ ��Ű��, ���� �ڵ忡�� v_out�� ������ �Ѱ迡 ���߾� �����Ѵ�.
	}

	// ������ �Ѱ踦 ����� �ӵ� ��� ���
	while (check_flag == 1) {			//check_flag == 1�� ��� v_out �迭�� �ӵ� ��� �� ��� �ϳ��� ���� �Ǿ����� �ǹ��Ѵ�.
		check_flag = 0;
		for (k = WP_size - 1; k > 1; k--) {
			if (v_out[k] < v_out[k - 1]) {	// ������ �Ǵ� ���� �˻�
				v_dec = sqrt(v_out[k] * v_out[k] - 2 * a_dec_max*h_in[k - 1]);	// ���� �������� k��° ��ġ�� �ӵ� ����� ������Ű�� ���� k-1��°���� ���� �� �ִ� �ִ� �Ѱ� �ӵ�
				if (v_dec < v_out[k - 1]) { // �ִ� �Ѱ� �ӵ��� k-1��°������ ���� �ӵ� ��ɺ��� ���� ���
					v_out[k - 1] = v_dec;	// k-1��°������ �ӵ� ����� �ִ� �Ѱ� �ӵ��� �����Ѵ�.
					check_flag = 1;			// �ӵ��� ���� �Ǿ���
				}
			}
		}

		for (k = 1; k < WP_size; k++) {
			if (v_out[k] > v_out[k - 1]) {	// ������ �Ǵ� ���� �˻�
				v_acc = sqrt(v_out[k - 1] * v_out[k - 1] + 2 * a_acc_max*h_in[k - 1]);	// k-1��°���� �ִ� �����Ͽ� k��°���� ������ �� �ִ� �ִ� �Ѱ� �ӵ�
				if (v_acc < v_out[k]) {		// �ִ� �Ѱ� �ӵ��� k��°������ ���� �ӵ� ��ɺ��� ���� ���
					v_out[k] = v_acc;		// k��°������ �ӵ� ����� �ִ� �Ѱ� �ӵ��� �����Ѵ�.
					check_flag = 1;			// �ӵ��� ���� �Ǿ���
				}
			}
		}
	}

	// �ʱ� �ӵ��� ����� ���� ���� �˻�
	for (k = 1; k < WP_size; k++) {
		if (v_out[k] < v_out[k - 1]) {		// ������ �Ǵ� ���� �˻�
			v_dec = sqrt(v_out[k - 1] * v_out[k - 1] + 2 * a_dec_max*h_in[k - 1]);	// k-1���� k���� ���� ������ �ִ� �Ѱ� �ӵ�
			if (v_dec > v_out[k]) {			// �ִ� �Ѱ� �ӵ��� k��°������ ���� �ӵ� ��ɺ��� Ŭ ���
				v_out[k] = v_dec;			// k��°������ �ӵ� ����� �ִ� �Ѱ� �ӵ��� �����Ѵ�.
			}
		}
	}

	// velocity command�� ������ ���� ���� �� ��� ���� �ӵ��� ����
	for (int i = 0; i < WP_size; i++) {
		if (v_out[i] > 100 || v_out[i] < 0) v_out[i] = lowest_velocity_command;
	}

	// velocity command output
	RTT_output->velocity_command = new double[WP_size];
	for (int i = 0; i < WP_size; i++) RTT_output->velocity_command[i] = v_out[i];


	////////// Section 3. �Ÿ��� ���� �ӵ� ����� �ð��� ���� �ӵ� ������� ��ȯ //////////

	// v(x) -> v(t)�� ��ȯ	
	// Cubic spline interpolation (�� �Ÿ� ������ �ӵ� �׷����� 3�� ���׽��� ��� a_i, b_i, c_i, d_i�� ���ϴ� ����)
	// Tridiagonal matrix�� ���, 1�� ���� �밢���, 2�� ���� �밢 ���, 3�� ���� �밢 ��� �����
	for (i = 0; i < WP_size - 2; i++)
	{
		element_1[i] = h_in[i];						// ���� �밢 ���
		element_2[i] = 2.*(h_in[i] + h_in[i + 1]);	// ���� �밢 ���
		element_3[i] = h_in[i + 1];					// ���� �밢 ���
	}
	// Tridiagonal matrix �Ŀ����� �캯 ����
	for (i = 0; i < WP_size - 2; i++) {
		element_4[i] = 6.0*((v_out[i + 2] - v_out[i + 1]) / h_in[i + 1] - (v_out[i + 1] - v_out[i]) / h_in[i]);	// �캯 ����
	}
	// ���ߴ밢��� solve
	solve_tridiagonal_in_place_destructive(element_4, WP_size - 2, element_1, element_2, element_3);	//���ߴ밢����� Ǯ�� (Ǯ�� �� element_4�� ��(S_i)�� �������ִ� �Լ�)

	// S_i: second derivate ��, �ʱ� S_i ���� ������ S_i ���� 0���� ���� --> natural cubic spline
	S_i[0] = S_i[WP_size - 1] = 0.0;

	// S_i �ʱⰪ�� ���������� �����Ͽ� ���ߴ밢����� ��� ��(S_i == element_4)�� ����
	for (i = 1; i < WP_size - 1; i++) {
		S_i[i] = element_4[i - 1];	// S_i�� �ش� ��ġ������ ���ӵ��� �ǹ���
	}

	// �� ���׽��� ��� ����
	for (i = 0; i < WP_size - 1; i++)
	{
		a_i[i] = (S_i[i + 1] - S_i[i]) / (6.0*h_in[i]);														// Cubic spline ���׽��� ���
		b_i[i] = S_i[i] / 2.0;																				// Cubic spline ���׽��� ���
		c_i[i] = (v_out[i + 1] - v_out[i]) / h_in[i] - (2.0*h_in[i] * S_i[i] + h_in[i] * S_i[i + 1]) / 6.0;	// Cubic spline ���׽��� ���
		d_i[i] = v_out[i];																					// Cubic spline ���׽��� ���
	}

	// ���б� (v(x) -> v(t)�� ��ȯ)
	velocity_indx = 0;		// ���� �Ÿ� �迭(input_distance)�� ���Ǵ� �ε����� 0���� �ʱ�ȭ
	t_current = 0;			// ���� ���� �ð� �ʱⰪ
	num_state = 2;			// absh3_for_control ���б� ����� ���� ���� ������ ����
	intcount = 1;			// absh3_for_control ���б� ����� ���� �ʱ� ����

	define_Y_vector_for_control(input_distance, v_out);		// ���� ���� �ʱⰪ ����

// 	// �ð��� ���� �Ÿ� �� �ӵ��� �ʱ� �� ����
// 	distance_by_time[0] = input_distance[0];	// �ð��� ���� �̵� �Ÿ�
// 	velocity_by_time[0] = v_out[0];				// �ð��� ���� ���� �ӵ�
	RTT_output->distance_by_time.assign(1,0);
	RTT_output->velocity_by_time.assign(1,0);
	RTT_output->distance_by_time[0] = input_distance[0];
	RTT_output->velocity_by_time[0] = v_out[0];

	while (velocity_indx < WP_size - 2)		// indx�� 0���� �Ÿ��� ���� �ӵ� ����� ���κп� ������ ������ �ݺ�
	{
		dqddq2Yp_for_control(input_distance);	// ���� �� ���� (�̺й����� ����)
		absh3_for_control(h, num_state, &intcount, AW1, AW, &t_current, Y, Yp, Y_next, &t_next);	// ���б� ���� (�ʱ� ���� ������ h���� ������ ����)
		for (i = 0; i < num_state; i++)
		{
			Y[i] = Y_next[i];	// state update
		}

		t_current = t_next;		// ���� �ð� ����

		if (Y[0] > input_distance[velocity_indx + 1]) {		// ���� ���� ������ �ش� indx ������ �ʰ��� ��� (���⼭ indx�� �Ÿ��� ���� �ӵ� �迭�� indx�� ��Ÿ��)
			velocity_indx++;									// indx�� 1 �߰���Ŵ
		}

		if (intcount == 5 || intcount >= 7) {	// absh3�� ���� step�� �ʹݿ� 0.01�� �������� �ʱ� ������ 0.01�ʰ� �� ������ ���е� ����� ������Ʈ ������ (2016.11.18 ȫȿ��)
// 			distance_by_time[time_indx] = Y[0];	// ���� ����� �ð��� ���� �Ÿ� �迭�� ����
// 			velocity_by_time[time_indx] = Y[1];	// ���� ����� �ð��� ���� �ӵ� �迭�� ����
			RTT_output->distance_by_time.push_back(Y[0]);
			RTT_output->velocity_by_time.push_back(Y[1]);
			time_indx++;		// �ð��� ���� �ӵ��� ��� ������ ������ �������� �������Ƿ� ������ �ε����� ����Ѵ�.
		}
		if (time_indx > time_indx_MAX)	// ���� �ð��� �ʹ� ������� time_indx_MAX ������ ����ϰ� ����
			break;
	}
	RTT_output->DataLength = RTT_output->distance_by_time.size();

#if !OptimalCPUTimeCheck
	FILE *fp;					// ��� ���� �ؽ�Ʈ ���Ͽ� �����ϱ� ���� ���� ������
	char buf[255];				// �ؽ�Ʈ ������ ��ο� ���ϸ��� �����ϱ� ���� ����
	sprintf_s(buf,sizeof(char)*255, "simulation_results/%d/data_sim_%d/_sim_flag4_velocity_profile_Cpp.dat",RTT_output->sim_count, data_indx);		// �ؽ�Ʈ ���� ��� �� �̸� ����

	for (k = 0; k < WP_size; k++) {
		if (k == 0) {
			//fp = fopen(buf, "w+");
			fopen_s(&fp, buf, "w+");
		}
		else {
			//fp = fopen(buf, "a+");
			fopen_s(&fp, buf, "a+");
		}

		fprintf_s(fp, "%d\t%10.10f\t%10.10f\t%10.10f\n", k, input_distance[k], v_in_temp[k], v_out[k]);	// ��� �����͸� �ؽ�Ʈ ���Ͽ� ����
		fclose(fp);
	}
#endif
	delete[] v_in_temp;
	opt_end = omp_get_wtime();	// CPU time ���� ����
}

void optimal_velocity_planning::absh3_for_control(double step_size, const int n, int *intcount, double AW1[2][2], double AW[2][2], double *t_current, double Y[2], double Yp[2], double Y_next[2], double *t_next)
{
	int i;

	/*  ABSH3  : constant step Adams Bashforth 3rd order formulation.
	written by Sung-Soo Kim
	Date: Oct. 19, 1998
	copyright reserved by Sung-Soo Kim

	input variables
	t_current: current time
	Y : current state
	Yp : current derivative of state
	step_size: integration step_size

	output variables
	Y_next : state at next time step
	t_next : Next time

	STARTER:  upto 2h, i.e., derivatives are stored for the initial  time steps at 0, h, 2h, to form
	3rd order  Adams Bashforth formula */

	switch (*intcount)
	{
		case 1:
			// Forward Euler method with 0.25 step_size for initial step
			// use derivative information at 0 step
			// y=y+step_size*Yp/4.0;
			for (i = 0; i<n; i++)
			{
				Y_next[i] = Y[i] + step_size*Yp[i] / 4;
			}
			// w(:,2) = Yp;
			for (i = 0; i<n; i++)
			{
				AW[i][1] = Yp[i];
			}
			// w1(:,2) = Yp;
			for (i = 0; i<n; i++)
			{
				AW1[i][1] = Yp[i];
			}
			*intcount = *intcount + 1;
			*t_next = *t_current + step_size / 4;
			break;
		case 2:
			// Adams Bashforth 2nd order method with 0.25 step_size for 2nd step
			// use derivative inforamtion at 0, h/4
			// y = y + step_size_h * ( 3.0*Yp - w1(:,2))/8.0;
			for (i = 0; i<n; i++)
			{
				Y_next[i] = Y[i] + step_size*(3 * Yp[i] - AW1[i][1]) / 8;
			}
			// w1(:,1) = Yp;
			for (i = 0; i<n; i++)
			{
				AW1[i][0] = Yp[i];
			}
			*intcount = *intcount + 1;
			*t_next = *t_current + step_size / 4;
			break;
		case 3:
			// Adams Bashforth 3rd order method with 0.25 step_size for 3rd step
			// use derivative information at 0, h/4, h/2
			// y = y + step_size * ( 23.0*Yp - 16.0*w1(:,1) + 5.0*w1(:,2))/48.0;
			for (i = 0; i<n; i++)
			{
				Y_next[i] = Y[i] + step_size*(23.0*Yp[i] - 16.0*AW1[i][0] + 5.0*AW1[i][1]) / 48.0;
			}
			// w1(:,2) = w1(:,1);
			for (i = 0; i<n; i++)
			{
				AW1[i][1] = AW1[i][0];
			}
			// w1(:,1) = Yp;
			for (i = 0; i<n; i++)
			{
				AW1[i][0] = Yp[i];
			}
			*intcount = *intcount + 1;
			*t_next = *t_current + step_size / 4.0;
			break;
		case 4:
			// Adams Bashforth 3rd order method with 0.25 step_size for 4th step
			// use derivative information at h/4, h/2, 3h/4
			// y = y + step_size * ( 23.0*Yp - 16.0*w1(:,1) + 5.0*w1(:,2))/48.0;
			for (i = 0; i<n; i++)
			{
				Y_next[i] = Y[i] + step_size*(23.0*Yp[i] - 16.0*AW1[i][0] + 5.0*AW1[i][1]) / 48.0;
			}
			// w1(:,2) = w(:,2);
			for (i = 0; i<n; i++)
			{
				AW1[i][1] = AW[i][1];
			}
			*intcount = *intcount + 1;
			*t_next = *t_current + step_size / 4.0;
			break;
		case 5:
			// Adams Bashforth 3rd order method with 0.5 step_size for 5th step
			// use derivative information at 0, h/2, h
			// y = y + step_size * ( 23.0*Yp - 16.0*w1(:,1) + 5.0*w1(:,2))/24.0;
			for (i = 0; i<n; i++)
			{
				Y_next[i] = Y[i] + step_size*(23.0*Yp[i] - 16.0*AW1[i][0] + 5.0*AW1[i][1]) / 24.0;
			}
			// w(:,1) = Yp;
			for (i = 0; i<n; i++)
			{
				AW[i][0] = Yp[i];
			}
			// w1(:,2) = w1(:,1);
			for (i = 0; i<n; i++)
			{
				AW1[i][1] = AW1[i][0];
			}
			// w1(:,1) = Yp;
			for (i = 0; i<n; i++)
			{
				AW1[i][0] = Yp[i];
			}
			*intcount = *intcount + 1;
			*t_next = *t_current + step_size / 2.0;
			break;
		case 6:
			// Adams Bashforth 3rd order method with 0.5 step_size for 6th step
			// use derivative information at h/2, h,  3h/2
			// y = y + step_size * ( 23.0*Yp - 16.0*w1(:,1) + 5.0*w1(:,2))/24.0;
			for (i = 0; i<n; i++)
			{
				Y_next[i] = Y[i] + step_size*(23.0*Yp[i] - 16.0*AW1[i][0] + 5.0*AW1[i][1]) / 24.0;
			}
			// w1(:,2) = w1(:,1);
			for (i = 0; i<n; i++)
			{
				AW1[i][1] = AW1[i][0];
			}
			// w1(:,1) = Yp;
			for (i = 0; i<n; i++)
			{
				AW1[i][0] = Yp[i];
			}
			*intcount = *intcount + 1;
			*t_next = *t_current + step_size / 2.0;
			break;
		case 7:
			// Adams Bashforth 3rd order method with step_size for 7th step
			// use derivative information at 0,  h,  2h
			// y = y + step_size * ( 23.0*Yp - 16.0*w(:,1) + 5.0*w(:,2))/12.0;
			for (i = 0; i<n; i++)
			{
				Y_next[i] = Y[i] + step_size*(23.0*Yp[i] - 16.0*AW[i][0] + 5.0*AW[i][1]) / 12.0;
			}
			// w(:,2) = w(:,1);
			for (i = 0; i<n; i++)
			{
				AW[i][1] = AW[i][0];
			}
			// w(:,1) = Yp;
			for (i = 0; i<n; i++)
			{
				AW[i][0] = Yp[i];
			}
			*intcount = *intcount + 1;
			*t_next = *t_current + step_size;
			break;
		default:
			// Adams Bashforth 3rd order method with step_size for more than 8th step
			// use derivative information t_current-2h, t_current-h, t_current
			// y = y + step_size * ( 23.0*Yp - 16.0*w(:,1) + 5.0*w(:,2))/12.0;
			for (i = 0; i<n; i++)
			{
				Y_next[i] = Y[i] + step_size*(23.0*Yp[i] - 16.0*AW[i][0] + 5.0*AW[i][1]) / 12.0;
			}
			// w(:,2) = w(:,1);
			for (i = 0; i<n; i++)
			{
				AW[i][1] = AW[i][0];
			}
			// w(:,1) = Yp;
			for (i = 0; i<n; i++)
			{
				AW[i][0] = Yp[i];
			}
			*intcount = *intcount + 1;
			*t_next = *t_current + step_size;
			break;
	}
}

void optimal_velocity_planning::define_Y_vector_for_control(double* input_distance, double* input_velocity)
{
	Y[0] = input_distance[0];
	Y[1] = input_velocity[0];
}

void optimal_velocity_planning::dqddq2Yp_for_control(double *input_distance)
{
	Yp[0] = Y[1];
	Yp[1] = (c_i[velocity_indx] + (Y[0] - input_distance[velocity_indx])*(2 * b_i[velocity_indx] + (3 * a_i[velocity_indx] * (Y[0] - input_distance[velocity_indx]))))*Y[1];
}

void optimal_velocity_planning::solve_tridiagonal_in_place_destructive(double *x, int X, double * a, double * b, double * c)
{
	/*
	solves Ax = v where A is a tridiagonal matrix consisting of vectors a, b, c
	x - initially contains the input vector v, and returns the solution x. indxed from 0 to X - 1 inclusive
	X - number of equations (length of vector x)
	a - subdiagonal (means it is the diagonal below the main diagonal), indxed from 1 to X - 1 inclusive
	b - the main diagonal, indxed from 0 to X - 1 inclusive
	c - superdiagonal (means it is the diagonal above the main diagonal), indxed from 0 to X - 2 inclusive

	Note: contents of input vector c will be modified, making this a one-time-use function (scratch space can be allocated instead for this purpose to make it reusable)
	Note 2: We don't check for diagonal dominance, etc.; this is not guaranteed stable
	*/

	/* indx variable is an unsigned integer of same size as pointer */
	int ix;

	c[0] = c[0] / b[0];
	x[0] = x[0] / b[0];

	/* loop from 1 to X - 1 inclusive, performing the forward sweep */
	for (ix = 1; ix < X; ix++) {
		double m = 1.0f / (b[ix] - a[ix] * c[ix - 1]);
		c[ix] = c[ix] * m;
		x[ix] = (x[ix] - a[ix] * x[ix - 1]) * m;
	}

	/* loop from X - 2 to 0 inclusive (safely testing loop condition for an unsigned integer), to perform the back substitution */
	for (ix = X - 1; ix-- > 0;)
		x[ix] = x[ix] - c[ix] * x[ix + 1];
}