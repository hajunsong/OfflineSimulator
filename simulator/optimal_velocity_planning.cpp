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
	WP_size = OID->WaypointSize;		// waypoint 개수 자동 입력
	
	input_distance = new double[WP_size];	// input_distance 배열 메모리 할당 input_distance[WP_size]
	v_in = new double[WP_size];				// v_in 배열 메모리 할당 v_in[WP_size]
	v_out = new double[WP_size];			// v_out 배열 메모리 할당 v_in[WP_size]
	h_in = new double[WP_size - 1];			// h_in 배열 메모리 할당 v_in[WP_size - 1]	

	v_x = OID->CurrentVelocity;			// 차량의 현재 주행 속도 입력

	element_1 = new double[WP_size - 2];	// Cubic spline 풀이에 사용하기 위한 배열 메모리 할당 (삼중 대각 행렬에서의 하위 대각 행렬)
	element_2 = new double[WP_size - 2];	// Cubic spline 풀이에 사용하기 위한 배열 메모리 할당 (삼중 대각 행렬에서의 메인 대각 행렬)
	element_3 = new double[WP_size - 2];	// Cubic spline 풀이에 사용하기 위한 배열 메모리 할당 (삼중 대각 행렬에서의 상위 대각 행렬)
	element_4 = new double[WP_size - 2];	// Cubic spline 풀이에 사용하기 위한 배열 메모리 할당 (우변 벡터)
	S_i = new double[WP_size];				// Cubic spline 풀이에 사용하기 위한 배열 메모리 할당 (가속도)
	a_i = new double[WP_size - 1];			// Cubic spline 풀이에 사용하기 위한 배열 메모리 할당 (다항식의 계수)
	b_i = new double[WP_size - 1];			// Cubic spline 풀이에 사용하기 위한 배열 메모리 할당 (다항식의 계수)
	c_i = new double[WP_size - 1];			// Cubic spline 풀이에 사용하기 위한 배열 메모리 할당 (다항식의 계수)
	d_i = new double[WP_size - 1];			// Cubic spline 풀이에 사용하기 위한 배열 메모리 할당 (다항식의 계수)

	h = OID->step_size;	// 스텝 사이즈. 동역학 모델에서 입력으로 받아와야 하지만 임시로 그냥 사용, 추후 수정 할 것(2017-02-28 song hajun)
}

void optimal_velocity_planning::run(ROUTPUTDATA *RTT_output) {

	////////// 입력 받아야할 차량 모델 파라미터 //////////
	double a_acc_max = 4;			// 최대 가속도
	double a_dec_max = -2;			// 최대 감속도
	double lowest_velocity_command = 1.0;	// 모든 thread의 주행 안정성이 fail 되었을 때 설정할 최저 속도

											////////// Sectoin 0. 지역 변수 설정 //////////
	int WP_thread, WP_indx;	// CPU 쓰레드 번호, waypoint의 순번(인덱스)
	int i, k;					// for 반복문 사용을 위한 변수
	int time_indx = 1;			// 시간에 대한 속도의 경우 데이터 개수가 기존보다 많아지므로 별도의 인덱스를 사용한다. time_indx = 1부터 새로 적분된 데이터 저장 (2016.11.18 홍효성)
	int check_flag = 1;			// 가감속 적용하여 속도 최적화 할 때 속도가 기존과 달라졌는지 여부를 확인하기 위한 flag
	double v_acc, v_dec;		// 최대 가속도가 적용된 속도, 최대 감속도가 적용된 속도	
	double WP_previous[4], WP_next[4];		// 차량이 현재 소속된 waypoint 및 그 다음 waypoint에서 위치, 속도, 주행성지표결과 저장을 위한 배열

	double *v_in_temp = new double[WP_size];

	opt_start = omp_get_wtime();		// CPU time 측정 시작

	////////// Section 1. 주행 안정성 지표를 바탕으로 waypoint 별 최대 주행 허용 속도 계산 //////////

// 	// 7가지 thread들에 있는 속도 중에서 stability indx가 pass 처리된 가장 빠른 속도를 해당 WP의 속도로 지정
 	for (WP_indx = 0; WP_indx < WP_size; WP_indx++) {
 		for (WP_thread = 0; WP_thread < 7; WP_thread++) {	// thread 순서는 속도가 느린 것(0번 thread)부터 빠른 것(6번 thread)까지
 			//if (WP[WP_indx][3][WP_thread] == 0) {			// 가장 낮은 속도 후보군부터 stability indx가 0 (fail)인 부분 찾기
 			if(OID->StabilityIndex[WP_indx][WP_thread] == 0) {
 				if (WP_thread == 0) {
 					//WP[WP_indx][2][7] = lowest_velocity_command;		// 최저 속도 thread에서 주행 안정성 지표 fail일 경우 최저 속도 강제 입력 (2016.12.22)
 					v_in_temp[WP_indx] = lowest_velocity_command;
 				}
 				else {
 					//WP[WP_indx][2][7] = WP[WP_indx][2][WP_thread - 1];	// 지표가 fail 되기 이전의 thread의 속도 후보군을 해당 WP에서의 최대 허용 속도로 할당
 					v_in_temp[WP_indx] = OID->VehicleVelocity[WP_indx][WP_thread - 1];
 				}
 				break;
 			}
 			//WP[WP_indx][2][7] = WP[WP_indx][2][WP_thread];	// 모든 thread의 지표가 만족한 경우 최대 속도인 6번 thread의 속도를 할당
 			v_in_temp[WP_indx] = OID->VehicleVelocity[WP_indx][WP_thread];
 		}
 	}
// 7가지 thread들에 있는 속도 중에서 stability index가 pass 처리된 가장 빠른 속도를 해당 WP의 속도로 지정
//	for (WP_indx = 0; WP_indx < WP_size; WP_indx++) {
//		for (WP_thread = 6; WP_thread >= 0; WP_thread--) {	// thread 순서는 속도가 빠른 것(6번 thread)부터 느린 것(0번 thread)까지
//// 			if (WP[WP_indx][3][WP_thread] == 1) {			// 가장 빠른 속도 후보군부터 stability index가 1 (pass)인 부분 찾기
//// 				WP[WP_indx][2][7] = WP[WP_indx][2][WP_thread];	// 지표가 fail 되기 이전의 thread의 속도 후보군을 해당 WP에서의 최대 허용 속도로 할당
//// 				break;
//// 			}
//// 			WP[WP_indx][2][7] = lowest_velocity_command;	// 최저 속도 thread에서 주행 안정성 지표 fail일 경우 최저 속도 1.5 강제 입력 (2016.12.22)
//			if (OID->StabilityIndex[WP_indx][WP_thread] == 1) {	// 가장 빠른 속도 후보군부터 stability index가 1 (pass)인 부분 찾기
//				v_in_temp[WP_indx] = OID->VehicleVelocity[WP_indx][WP_thread];	// 지표가 fail 되기 이전의 thread의 속도 후보군을 해당 WP에서의 최대 허용 속도로 할당
//				break;
//			}
//			v_in_temp[WP_indx] = lowest_velocity_command;	// 최저 속도 thread에서 주행 안정성 지표 fail일 경우 최저 속도 1.5 강제 입력 (2016.12.22)
//		}
//	}

	//WP[WP_size - 1][2][7] = WP[WP_size - 2][2][7];			// 마지막 waypoint의 속도 후보군은 직전 waypoint의 속도 후보군과 같도록 설정
	v_in_temp[WP_size - 1] = v_in_temp[WP_size - 2];


	////////// Section 2. 가감속 한계를 고려한 속도 명령 계산 //////////

	// 위치(x축)와 속도(y축) 입력
	for (WP_indx = 0; WP_indx < WP_size; WP_indx++) {
		if (WP_indx == 0) {
			input_distance[0] = 0;		// 초기 위치 0
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

			input_distance[WP_indx] = input_distance[WP_indx - 1] + sqrt((WP_next[0] - WP_previous[0])*(WP_next[0] - WP_previous[0]) + (WP_next[1] - WP_previous[1])*(WP_next[1] - WP_previous[1]));	// 누적 거리
			v_in[WP_indx] = WP_next[2];	// 해당 waypoint에서의 최대 주행 허용 속도 입력
		}
	}

	// 각 waypoint 사이의 거리 간격 계산
	for (i = 0; i < WP_size - 1; i++) {
		h_in[i] = input_distance[i + 1] - input_distance[i];		// h_i = x(i + 1) - x(i); x축(여기서는 거리) 간격
	}

	//v_in[0] = v_in_temp[0];		// 초기 속도 명령을 현재속도로 설정하는 코드 무시

	for (i = 0; i < WP_size; i++) {
		v_out[i] = v_in[i];          // 먼저 v_out을 v_in과 동기화 시키고, 다음 코드에서 v_out을 가감속 한계에 맞추어 변경한다.
	}

	// 가감속 한계를 고려한 속도 명령 계산
	while (check_flag == 1) {			//check_flag == 1인 경우 v_out 배열의 속도 명령 중 어느 하나라도 변경 되었음을 의미한다.
		check_flag = 0;
		for (k = WP_size - 1; k > 1; k--) {
			if (v_out[k] < v_out[k - 1]) {	// 감속이 되는 구간 검사
				v_dec = sqrt(v_out[k] * v_out[k] - 2 * a_dec_max*h_in[k - 1]);	// 감속 구간에서 k번째 위치의 속도 명령을 만족시키기 위해 k-1번째에서 가질 수 있는 최대 한계 속도
				if (v_dec < v_out[k - 1]) { // 최대 한계 속도가 k-1번째에서의 기존 속도 명령보다 작을 경우
					v_out[k - 1] = v_dec;	// k-1번째에서의 속도 명령을 최대 한계 속도로 변경한다.
					check_flag = 1;			// 속도가 변경 되었음
				}
			}
		}

		for (k = 1; k < WP_size; k++) {
			if (v_out[k] > v_out[k - 1]) {	// 가속이 되는 구간 검사
				v_acc = sqrt(v_out[k - 1] * v_out[k - 1] + 2 * a_acc_max*h_in[k - 1]);	// k-1번째에서 최대 가속하여 k번째까지 도달할 수 있는 최대 한계 속도
				if (v_acc < v_out[k]) {		// 최대 한계 속도가 k번째에서의 기존 속도 명령보다 작을 경우
					v_out[k] = v_acc;		// k번째에서의 속도 명령을 최대 한계 속도로 변경한다.
					check_flag = 1;			// 속도가 변경 되었음
				}
			}
		}
	}

	// 초기 속도를 고려한 감속 구간 검사
	for (k = 1; k < WP_size; k++) {
		if (v_out[k] < v_out[k - 1]) {		// 감속이 되는 구간 검사
			v_dec = sqrt(v_out[k - 1] * v_out[k - 1] + 2 * a_dec_max*h_in[k - 1]);	// k-1에서 k까지 감속 가능한 최대 한계 속도
			if (v_dec > v_out[k]) {			// 최대 한계 속도가 k번째에서의 기존 속도 명령보다 클 경우
				v_out[k] = v_dec;			// k번째에서의 속도 명령을 최대 한계 속도로 변경한다.
			}
		}
	}

	// velocity command에 쓰레기 값이 포함 될 경우 최저 속도로 보정
	for (int i = 0; i < WP_size; i++) {
		if (v_out[i] > 100 || v_out[i] < 0) v_out[i] = lowest_velocity_command;
	}

	// velocity command output
	RTT_output->velocity_command = new double[WP_size];
	for (int i = 0; i < WP_size; i++) RTT_output->velocity_command[i] = v_out[i];


	////////// Section 3. 거리에 대한 속도 명령을 시간에 대한 속도 명령으로 변환 //////////

	// v(x) -> v(t)로 변환	
	// Cubic spline interpolation (각 거리 구간별 속도 그래프의 3차 다항식의 계수 a_i, b_i, c_i, d_i를 구하는 과정)
	// Tridiagonal matrix의 요소, 1은 하위 대각행렬, 2는 메인 대각 행렬, 3은 상위 대각 행렬 요소임
	for (i = 0; i < WP_size - 2; i++)
	{
		element_1[i] = h_in[i];						// 하위 대각 행렬
		element_2[i] = 2.*(h_in[i] + h_in[i + 1]);	// 메인 대각 행렬
		element_3[i] = h_in[i + 1];					// 상위 대각 행렬
	}
	// Tridiagonal matrix 식에서의 우변 벡터
	for (i = 0; i < WP_size - 2; i++) {
		element_4[i] = 6.0*((v_out[i + 2] - v_out[i + 1]) / h_in[i + 1] - (v_out[i + 1] - v_out[i]) / h_in[i]);	// 우변 벡터
	}
	// 삼중대각행렬 solve
	solve_tridiagonal_in_place_destructive(element_4, WP_size - 2, element_1, element_2, element_3);	//삼중대각행렬을 풀이 (풀이 후 element_4에 해(S_i)를 저장해주는 함수)

	// S_i: second derivate 값, 초기 S_i 값과 마지막 S_i 값을 0으로 설정 --> natural cubic spline
	S_i[0] = S_i[WP_size - 1] = 0.0;

	// S_i 초기값과 마지막값을 포함하여 삼중대각행렬의 결과 값(S_i == element_4)을 대입
	for (i = 1; i < WP_size - 1; i++) {
		S_i[i] = element_4[i - 1];	// S_i는 해당 위치에서의 가속도를 의미함
	}

	// 각 다항식의 계수 저장
	for (i = 0; i < WP_size - 1; i++)
	{
		a_i[i] = (S_i[i + 1] - S_i[i]) / (6.0*h_in[i]);														// Cubic spline 다항식의 계수
		b_i[i] = S_i[i] / 2.0;																				// Cubic spline 다항식의 계수
		c_i[i] = (v_out[i + 1] - v_out[i]) / h_in[i] - (2.0*h_in[i] * S_i[i] + h_in[i] * S_i[i + 1]) / 6.0;	// Cubic spline 다항식의 계수
		d_i[i] = v_out[i];																					// Cubic spline 다항식의 계수
	}

	// 적분기 (v(x) -> v(t)로 변환)
	velocity_indx = 0;		// 누적 거리 배열(input_distance)에 사용되는 인덱스를 0으로 초기화
	t_current = 0;			// 적분 누적 시간 초기값
	num_state = 2;			// absh3_for_control 적분기 사용을 위한 상태 변수의 개수
	intcount = 1;			// absh3_for_control 적분기 사용을 위한 초기 변수

	define_Y_vector_for_control(input_distance, v_out);		// 상태 변수 초기값 정의

// 	// 시간에 대한 거리 및 속도의 초기 값 저장
// 	distance_by_time[0] = input_distance[0];	// 시간에 따른 이동 거리
// 	velocity_by_time[0] = v_out[0];				// 시간에 따른 주행 속도
	RTT_output->distance_by_time.assign(1,0);
	RTT_output->velocity_by_time.assign(1,0);
	RTT_output->distance_by_time[0] = input_distance[0];
	RTT_output->velocity_by_time[0] = v_out[0];

	while (velocity_indx < WP_size - 2)		// indx가 0부터 거리에 대한 속도 명령의 끝부분에 도달할 때까지 반복
	{
		dqddq2Yp_for_control(input_distance);	// 상태 식 대입 (미분방정식 형태)
		absh3_for_control(h, num_state, &intcount, AW1, AW, &t_current, Y, Yp, Y_next, &t_next);	// 적분기 실행 (초기 적분 간격은 h보다 작음에 주의)
		for (i = 0; i < num_state; i++)
		{
			Y[i] = Y_next[i];	// state update
		}

		t_current = t_next;		// 적분 시간 누적

		if (Y[0] > input_distance[velocity_indx + 1]) {		// 적분 누적 구간이 해당 indx 범위를 초과할 경우 (여기서 indx는 거리에 대한 속도 배열의 indx를 나타냄)
			velocity_indx++;									// indx를 1 추가시킴
		}

		if (intcount == 5 || intcount >= 7) {	// absh3의 적분 step이 초반에 0.01씩 증가하지 않기 때문에 0.01초가 될 때마다 적분된 결과를 업데이트 시켜줌 (2016.11.18 홍효성)
// 			distance_by_time[time_indx] = Y[0];	// 적분 결과를 시간에 따른 거리 배열에 저장
// 			velocity_by_time[time_indx] = Y[1];	// 적분 결과를 시간에 따른 속도 배열에 저장
			RTT_output->distance_by_time.push_back(Y[0]);
			RTT_output->velocity_by_time.push_back(Y[1]);
			time_indx++;		// 시간에 대한 속도의 경우 데이터 개수가 기존보다 많아지므로 별도의 인덱스를 사용한다.
		}
		if (time_indx > time_indx_MAX)	// 적분 시간이 너무 길어지면 time_indx_MAX 까지만 계산하고 끊음
			break;
	}
	RTT_output->DataLength = RTT_output->distance_by_time.size();

#if !OptimalCPUTimeCheck
	FILE *fp;					// 결과 값을 텍스트 파일에 저장하기 위한 파일 포인터
	char buf[255];				// 텍스트 파일의 경로와 파일명을 저장하기 위한 변수
	sprintf_s(buf,sizeof(char)*255, "simulation_results/%d/data_sim_%d/_sim_flag4_velocity_profile_Cpp.dat",RTT_output->sim_count, data_indx);		// 텍스트 파일 경로 및 이름 정의

	for (k = 0; k < WP_size; k++) {
		if (k == 0) {
			//fp = fopen(buf, "w+");
			fopen_s(&fp, buf, "w+");
		}
		else {
			//fp = fopen(buf, "a+");
			fopen_s(&fp, buf, "a+");
		}

		fprintf_s(fp, "%d\t%10.10f\t%10.10f\t%10.10f\n", k, input_distance[k], v_in_temp[k], v_out[k]);	// 결과 데이터를 텍스트 파일에 저장
		fclose(fp);
	}
#endif
	delete[] v_in_temp;
	opt_end = omp_get_wtime();	// CPU time 측정 종료
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