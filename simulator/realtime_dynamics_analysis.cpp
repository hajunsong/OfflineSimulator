#include "realtime_dynamics_analysis.h"
#include <memory.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <omp.h>

#include <Windows.h>
#include <process.h>
#include <conio.h>
#include <map>

#include <direct.h>

#include "pnu_map_tire.h"

#define TINY 1.0e-20;

using namespace std;

uint	Obj_Map::_num_x; //포인터로 수정
uint	Obj_Map::_num_y;
double	Obj_Map::_gab;
double	Obj_Map::_x_c;
double	Obj_Map::_y_c;
double	*Obj_Map::_x_data = NULL; //갯수가 런타임 결정이기 때문에 배열 사용곤란
double	*Obj_Map::_y_data = NULL;
const double*const*	Obj_Map::_map_data = NULL; //이중포인터 지시내용의 변경 방지

Obj_Map::~Obj_Map()
{
	delete[] _x_data;
	_x_data = NULL;

	delete[] _y_data;
	_y_data = NULL;
}

void Obj_Map::PNU_set_map(double **map_data, uint rows, uint cols, double gab, double x_c, double y_c)
{
	double length_x, length_y;
	_map_data = map_data;
	//_map_data[0][0]=0.0;
	_num_x = cols;
	_num_y = rows;
	_gab = gab;
	_x_c = x_c;
	_y_c = y_c;

	length_x = (double)(_num_x - 1)*gab;
	length_y = (double)(_num_y - 1)*gab;

	delete[] _x_data;
	delete[] _y_data;

	_x_data = new double[_num_x];
	_y_data = new double[_num_y];

	fn_linspace(_x_data, length_x, _num_x, _x_c);
	fn_linspace(_y_data, length_y, _num_y, _y_c);
}

void Obj_Map::PNU_fn_map(double x_p, double y_p, double *z)
{
	double MBMt[4][4], Ut[4], V[4];
	double u, v;
	uint i;
	if (_map_data == NULL || _x_data == NULL || _y_data == NULL) {
		printf("\n 노면이 초기화 되지 않음! 높이값을 0으로 임시출력!! \n");
		*z = 0;
	}
	else {
		fn_MBMt_uv(MBMt, &u, &v, x_p, y_p);

		for (i = 0; i < 4; ++i) {
			Ut[i] = pow(u, (int)(i));
			V[i] = pow(v, (int)(i));
		}
		fn_mat1441_multi(V, MBMt, Ut, z);
	}
}

void Obj_Map::fn_linspace(double *data, double L, uint num_data, double c)
{
	double gab = L / (double)(num_data - 1);
	uint i;
	for (i = 0; i < num_data; ++i) {
		data[i] = gab*(double)(i)-c;
	}
}

void Obj_Map::fn_find_patch_uv(uint* patch_r, uint* patch_c, double* u, double* v, double x_p, double y_p)
{
	double gab_th = 0.00000001; //magic, 경계값 처리를 위함
	double x_p_m, y_p_m, x_len, y_len;
	uint sw_error;

	x_p_m = x_p;
	y_p_m = y_p;

	sw_error = 0;
	if (x_p_m < _x_data[0]) {
		x_p_m = _x_data[0];
		sw_error = 1;
	} //x축 시작
	if (x_p_m >= _x_data[_num_x - 1]) {
		x_p_m = _x_data[_num_x - 1] - gab_th;
		sw_error = 1;
	} //x축 종료, 끝값 문제로 >에서 >=으로 수정
	if (y_p_m < _y_data[0]) {
		y_p_m = _y_data[0];
		sw_error = 1;
	} //y축 시작
	if (y_p_m >= _y_data[_num_y - 1]) {
		y_p_m = _y_data[_num_y - 1] - gab_th;
		sw_error = 1;
	} //y축 종료, 끝값 문제로 >에서 >=으로 수정

	  // 에러처리 메세지 표시여부 확인
	if (sw_error) {
		//printf("xy_position이 data범위를 넘어섬 \n");
	}

	////////// patch 계산
	x_len = x_p_m + _x_c; //내부적으로는 기준 위치가 row, column의 첫 번째이기 때문
	y_len = y_p_m + _y_c;

	*patch_r = (uint)floor(y_len / _gab); // 끝값 문제로 (y_len/_gab+gab_th)에서 (y_len/_gab) 으로 수정
	*patch_c = (uint)floor(x_len / _gab);

	////////// uv계산
	*u = fmod(x_len, _gab) / _gab;
	*v = fmod(y_len, _gab) / _gab;
}

void Obj_Map::fn_MBMt_uv(double(*MBMt)[4], double* u, double* v, double x_p, double y_p)
{
	double M[4][4] = { { 1, 0, 0, 0 },{ 0, 0, 1, 0 },{ -3, 3, -2, -1 },{ 2, -2, 1, 1 } };
	double B[4][4], MB[4][4], Mt[4][4];
	uint patch_r, patch_c;
	uint i, j;
	fn_find_patch_uv(&patch_r, &patch_c, u, v, x_p, y_p);
	fn_B_mat(B, patch_r, patch_c);

	for (j = 0; j < 4; ++j) {
		for (i = 0; i < 4; ++i) {
			Mt[i][j] = M[j][i];
		}
	}

	fn_mat444_multi(MB, M, B);
	fn_mat444_multi(MBMt, MB, Mt);
}

void Obj_Map::fn_patch_idxs(int row_idx[4], int col_idx[4], uint patch_r, uint patch_c)
{
	int idx_set[4] = { -1, 0, 1, 2 };
	uint i;
	for (i = 0; i < 4; ++i) {
		col_idx[i] = idx_set[i] + (int)(patch_c);
		row_idx[i] = idx_set[i] + (int)(patch_r);

		if (col_idx[i] >= (int)(_num_x)) //uword인데 col_idx가 ivec으로 표현하기가 쉬워서...
			col_idx[i] = (int)(_num_x - 1);
		else if (col_idx[i] < 0)
			col_idx[i] = 0;

		if (row_idx[i] >= (int)(_num_y))
			row_idx[i] = (int)(_num_y - 1);
		else if (row_idx[i] < 0)
			row_idx[i] = 0;
	}
}

void Obj_Map::fn_B_mat(double(*B)[4], uint patch_r, uint patch_c)
{
	int row_idx[4];
	int col_idx[4];
	double tmp;
	double R[2][2], Rv[2][2], Ru[2][2], Ruv[2][2];
	uint i, j;

	fn_patch_idxs(row_idx, col_idx, patch_r, patch_c);

	R[0][0] = _map_data[row_idx[1]][col_idx[1]];
	R[0][1] = _map_data[row_idx[1]][col_idx[2]];
	R[1][0] = _map_data[row_idx[2]][col_idx[1]];
	R[1][1] = _map_data[row_idx[2]][col_idx[2]];

	tmp = 1 / (2.0*_gab);
	Rv[0][0] = (_map_data[row_idx[2]][col_idx[1]] - _map_data[row_idx[0]][col_idx[1]])*tmp;
	Rv[1][0] = (_map_data[row_idx[3]][col_idx[1]] - _map_data[row_idx[0]][col_idx[1]])*tmp;
	Rv[0][1] = (_map_data[row_idx[2]][col_idx[2]] - _map_data[row_idx[0]][col_idx[2]])*tmp;
	Rv[1][1] = (_map_data[row_idx[3]][col_idx[2]] - _map_data[row_idx[0]][col_idx[2]])*tmp;

	Ru[0][0] = (_map_data[row_idx[1]][col_idx[2]] - _map_data[row_idx[1]][col_idx[0]])*tmp;
	Ru[1][0] = (_map_data[row_idx[2]][col_idx[2]] - _map_data[row_idx[2]][col_idx[0]])*tmp;
	Ru[0][1] = (_map_data[row_idx[1]][col_idx[3]] - _map_data[row_idx[1]][col_idx[1]])*tmp;
	Ru[1][1] = (_map_data[row_idx[2]][col_idx[3]] - _map_data[row_idx[2]][col_idx[1]])*tmp;

	Ruv[0][0] = (Ru[0][0] + Rv[0][0])*0.5;
	Ruv[0][1] = (Ru[0][1] + Rv[0][1])*0.5;
	Ruv[1][0] = (Ru[1][0] + Rv[1][0])*0.5;
	Ruv[1][1] = (Ru[1][1] + Rv[1][1])*0.5;

	for (j = 0; j < 2; ++j) {
		for (i = 0; i < 2; ++i) {
			B[i][j] = R[i][j];
			B[i][j + 2] = Ru[i][j];
			B[i + 2][j] = Rv[i][j];
			B[i + 2][j + 2] = Ruv[i][j];
		}
	}
}

void Obj_Map::fn_mat444_multi(double(*C)[4], double(*A)[4], double(*B)[4])
{
	uint i, j, k;
	for (j = 0; j < 4; ++j) {
		for (i = 0; i < 4; ++i) {
			C[i][j] = 0.0;
			for (k = 0; k < 4; ++k) {
				C[i][j] = C[i][j] + A[i][k] * B[k][j];
			}
		}
	}
}

void Obj_Map::fn_mat1441_multi(double *A, double(*B)[4], double *C, double *output)
{
	double sol = C[0] * (A[0] * B[0][0] + A[1] * B[1][0] + A[2] * B[2][0] + A[3] * B[3][0]) + C[1] * (A[0] * B[0][1] + A[1] * B[1][1] + A[2] * B[2][1] + A[3] * B[3][1]) + C[2] * (A[0] * B[0][2] + A[1] * B[1][2] + A[2] * B[2][2] + A[3] * B[3][2]) + C[3] * (A[0] * B[0][3] + A[1] * B[1][3] + A[2] * B[2][3] + A[3] * B[3][3]);
	*output = sol;
}

double Obj_Nan::z_min = -100.0; // nan 판정기준 최소 높이

void Obj_Nan::fn_nan_interp(double **Map_data, uint n_rows, uint n_cols)
{
	uint i = 0, j = 0;
	for (i = 0; i < n_rows; ++i) {
		for (j = 0; j < n_cols; ++j) {

			fn_nan_patch(Map_data, n_rows, n_cols, i, j);

		}
	}
}


double Obj_Nan::fn_nan_patch(double **Map_data, uint n_rows, uint n_cols, uint i, uint j)
{
	//const volatile uint num_p=4; //C89 문법은 const 붙여도 배열초기화 불가
	enum array_size { num_p = 4 }; //탐색 방향의 수, 4 또는 8로 지정  // 4는 대각선방향, 8은 대각선 + 직선방향
								   //array 사이즈 지정을 위해 enum을 사용, define 사용 자제 위함
	double z, z_p;
	uint k, p;
	int s;
	int i_p, j_p;
	double pr, num, den;
	double k_end, z_end;

	uint sw_set[num_p] = { 0 }; // 경계값 도달 여부, 0 : nan, 1 : 데이터 도달, 2 : 경계 도달(nan인체로)
	double z_set[num_p] = { 0 };  // 경계값 도달시 높이값
	uint k_set[num_p] = { 1 };  // 경게값 도달시 인덱스거리

	z = Map_data[i][j];
	if (is_nan(z) != 1) {
		return z;
	}

	// 경계값 탐색 시작
	i_p = 0;
	j_p = 0;
	for (p = 0; p < num_p; ++p) {
		if (sw_set[p] != 0) { continue; }

		k = 0;
		while (1) {
			k++;
			fn_ij_p(&i_p, &j_p, p, k, i, j);
			if (i_p < 0 || j_p < 0 || i_p >= (int)(n_rows) || j_p >= (int)(n_cols)) {
				sw_set[p] = 2;
				break;
			}
			else {
				z_p = Map_data[i_p][j_p];
				if (is_nan(z_p)) {
					sw_set[p] = 0;
				}
				else {
					sw_set[p] = 1;
					z_set[p] = z_p;
					k_set[p] = k;
					break;
				}
			}
		}
	}

	// Nan 거리가중 평균
	num = 0.0;
	den = 0.0;
	s = 0;
	for (p = 0; p < num_p; p++) {
		if (sw_set[p] == 1) {
			s++;
			pr = 1.0 / k_set[p];
			num = num + z_set[p] * pr;
			den = den + pr;
		}
	}

	if (s != 0) {		//s가 0이면 보간할 경계값이 없기때문
		z = num / den;
		Map_data[i][j] = z;
	}


	// 탐색 경로 보간
	for (p = 0; p < num_p; ++p) {
		if (sw_set[p] == 2) { continue; }

		k_end = k_set[p];
		z_end = z_set[p];
		for (k = 0; k < k_end; ++k) { //중요 : matlab은 k<k_end-1 임, 차이점 이해하기
			fn_ij_p(&i_p, &j_p, p, k, i, j);
			z_p = (z*(k_end - k) + z_end*k) / k_end;
			Map_data[i_p][j_p] = z_p;
		}
	}

	return z;
}

int Obj_Nan::is_nan(double z)
{
	if (z < z_min) {
		return 1;
	}
	else {
		return 0;
	}
}

void Obj_Nan::fn_ij_p(int* i_p, int* j_p, uint p, uint k, uint i, uint j)
{
	switch (p) {
	case 0: *i_p = i + k;    *j_p = j + k; break;
	case 1: *i_p = i + k;    *j_p = j - k; break;
	case 2: *i_p = i - k;    *j_p = j - k; break;
	case 3: *i_p = i - k;    *j_p = j + k; break;
	case 4: *i_p = i + k;    *j_p = j;   break;
	case 5: *i_p = i - k;    *j_p = j;   break;
	case 6: *i_p = i;      *j_p = j + k; break;
	case 7: *i_p = i;      *j_p = j - k; break;
	}
}

Obj_Tire::Obj_Tire()
{
	_c_t = 0;
	_k_t = 300000;// 490000;
	_Iyy = 21;                               // moment of inertia
	_w = 0.38;                               // patch width
	_R_u = 0.548;                            // 타이어 반지름
	_t_s = 0.2;
	//_slope_rr = 1;                           // 구름저항 계산을 위한 0에서 omega의 최대기울기, 적분기와 연관

	_eps = 2.22044604925031e-16;
	_pi = 3.141592653589793;
	
	_P[0] = 1.5260080e+01;
	_P[1] = 1.4351226e+01;
	_P[2] = -1.7859839e+00;
	_P[3] = -6.7357190e-02;
	_P[4] = -1.7717619e+00;
	_P[5] = -2.9680650e+00;
	_P[6] = 7.2179859e-01;
	_P[7] = 2.6179968e+00;

	_c_x = exp(_P[0]);
	_c_y = exp(_P[1]);
}

void Obj_Tire::PNU_Pre_road(double R_c_RF[3], double R_c_RM[3], double R_c_RR[3], double R_c_LF[3], double R_c_LM[3], double R_c_LR[3]) {
	Obj_Map::PNU_fn_map(R_c_RF[0], R_c_RF[1], &_Road_h_i[0]);
	Obj_Map::PNU_fn_map(R_c_RM[0], R_c_RM[1], &_Road_h_i[1]);
	Obj_Map::PNU_fn_map(R_c_RR[0], R_c_RR[1], &_Road_h_i[2]);
	Obj_Map::PNU_fn_map(R_c_LF[0], R_c_LF[1], &_Road_h_i[3]);
	Obj_Map::PNU_fn_map(R_c_LM[0], R_c_LM[1], &_Road_h_i[4]);
	Obj_Map::PNU_fn_map(R_c_LR[0], R_c_LR[1], &_Road_h_i[5]);

	double P_RF[3] = { R_c_RF[0], R_c_RF[1], _Road_h_i[0] };
	double P_RR[3] = { R_c_RR[0], R_c_RR[1], _Road_h_i[2] };
	double P_LF[3] = { R_c_LF[0], R_c_LF[1], _Road_h_i[3] };
	double P_LR[3] = { R_c_LR[0], R_c_LR[1], _Road_h_i[5] };

	double T_vec1[3], T_vec2[3];
	for (int i = 0; i < 3; i++) {
		T_vec1[i] = P_RF[i] - P_LR[i];
		T_vec2[i] = P_LF[i] - P_RR[i];
	}

	fn_cross(_N_vec, T_vec1, T_vec2);

	fn_normalize(_N_vec, 3);
}

void Obj_Tire::PNU_Tire_force(double R_c[3], double dR_c[3], double u_vec[3], double omega, int id)
{
	// 수직력 계산
	double R_x = R_c[0];
	double R_y = R_c[1];
	//Obj_Map::PNU_fn_map(R_x, R_y, &_road_h);
	//road_h = -0.63; // PNU_fn_map(R_x,R_y);	//평지로 가정

	double R_z = R_c[2];
	double dR_z = dR_c[2];
	double point_h;
	point_h = R_z - _R_u;	//tire의 point follow점
	fn_max(_Road_h_i[id] - point_h, 0.0, &_pen);	//tire penetration 침투량
	//fn_max(_road_h - point_h, 0.0, &_pen);		//tire penetration 침투량

	double Fz;
	if (_pen > 0)
		Fz = _k_t * _pen - (_c_t * dR_z)*(1);
	else if (_pen < 0)
		Fz = _k_t * _pen - (_c_t * dR_z)*(-1);
	else
		Fz = _k_t * _pen - (_c_t * dR_z)*(0);

	// effective radius 계산
	double R_e, a;
	fn_max(_R_u - _pen, 0.0, &_R_d);
	fn_R_eff(&R_e, &a, _R_d);

	// 좌표계 생성

	//double N_vec[3] = { 0, 0, 1 };	// 법선 벡터 해제
	double A_road[3][3], A_tire[3][3];
	fn_road_tire_A(A_road, A_tire, _N_vec, u_vec);

	// 타이어 병진속도의 global->road(패치) 좌표변환

	double V_x, V_sy;
	double B_dR_c[3];
	fn_mat33T31(B_dR_c, A_road, dR_c); // A_tire에서 A_road로 변경, 이게 맞는것 같음
	V_x = B_dR_c[0];
	V_sy = B_dR_c[1];

	// slip_ratio, slip_angle 계산
	//fn_slip_ratio(&_slip, &_angle, V_x, omega, R_e, V_sy, Fz);
	fn_slip_ratio(&_slip_x, &_slip_y, V_x, omega, R_e, V_sy, Fz, a);

	// road(패치)기준 타이어 힘 계산

	double F_road[3], M_road[3];
	double vec_p2c[3], temp1[3];
	//fn_Fiala_tire(F_road, M_road, Fz, _slip, _angle, omega); //contact patch기준 iso

	_mu = 0.5;
	fn_Brush_tire(F_road, M_road, Fz, _slip_x, _slip_y, _mu, a, _pen, omega);
	
															 // 타이어 힘의 road->tire 좌표변환
	//vec_p2c[0] = 0; //tire의 point center에서 point follow점까지의 vector, 좌표계 변환때 사용
	//vec_p2c[1] = 0;
	//vec_p2c[2] = -_R_d;

	//vec_p2c[0] = -_N_vec[0]; //tire의 point center에서 point follow점까지의 vector, 좌표계 변환때 사용
	//vec_p2c[1] = -_N_vec[1];
	//vec_p2c[2] = -_N_vec[2]*_R_d;

	vec_p2c[0] = -_N_vec[0] * _R_d; //tire의 point center에서 point follow점까지의 vector, 좌표계 변환때 사용
	vec_p2c[1] = -_N_vec[1] * _R_d;
	vec_p2c[2] = -_N_vec[2] * _R_d;
			
	//_F_road : Contact patch에서 발생하는 road local force
	//_F_global : Wheel center에서 발생하는 global force
	//_F_tire : Wheel center에서 발생하는 tire local force

	fn_mat3331(_F_global, A_road, F_road);
	fn_mat3331(_M_global, A_road, M_road);

	fn_cross(temp1, vec_p2c, _F_global);
	for (uint i = 0; i < 3; i++) {
		_M_global[i] += temp1[i];
	}
	//fn_mat3331(_M_tire, A_road, M_road);
	fn_mat33T31(_F_tire, A_tire, _F_global);
	fn_mat33T31(_M_tire, A_tire, _M_global);

	_Fx = _F_tire[0];
	_Fy = _F_tire[1];
	_Fz = _F_tire[2];

	_My = _M_tire[1];
}

void Obj_Tire::PNU_get_data(double F[3], double M[3], double *slip, double *angle, double *road_h, double *pen, double *R_d, double *Fx, double *Fy, double *Fz, double *My)
{
	memcpy(F, _F_global, sizeof(double) * 3); // 	F = _F_tire;
	memcpy(M, _M_global, sizeof(double) * 3); // 	M = _M_tire;
	*slip = _slip_x;
	*angle = _slip_y;
	*road_h = _road_h;
	*pen = _pen;
	*R_d = _R_d;
	*Fx = _Fx;
	*Fy = _Fy;
	*Fz = _Fz;
	*My = _My;
}

void Obj_Tire::fn_Brush_tire(double F[3], double M[3], double Fz, double slip_x, double slip_y, double mu, double a, double pen, double omega)
{

	double s_x = fabs(slip_x);
	double s_y = fabs(slip_y);

	const double b = _w;

	a = a + _eps;

	double k_x = _c_x / (2 * b);
	double k_y = _c_y / (2 * b);

	double u_s_x = mu*(exp(_P[2]) + exp(pen*_P[3]));
	double	u_s_y = mu*(exp(_P[4]) + exp(pen*_P[5]));

	double u_d_x = u_s_x*(1 - s_x / (1 + exp(_P[6])));
	double u_d_y = u_s_y*(1 - s_y / (1 + exp(_P[7])));

	// Combined adhesion force

	double e = 3 * Fz / (4 * a*a*b);

	double S_s_x = u_s_x / k_x*e + _eps;
	double S_s_y = u_s_y / k_y*e + _eps;
	double S_s;
	fn_min(sqrt(pow(s_x / S_s_x, 2) + pow(s_y / S_s_y, 2)), 1.0, &S_s);
	double X_s;
	fn_max(2 * a*(1 - S_s), 0.0, &X_s);

	double F_ax = b*k_x*s_x*X_s*X_s;
	double F_ay = b*k_y*s_y*X_s*X_s;
	double M_az = (b*k_y*s_y*X_s*X_s*(3 * a - 2 * X_s)) / 3;

	// Combined s sliding force
	double den = sqrt(pow(u_d_x*s_x, 2) + pow(u_d_y*s_y, 2) + _eps);
	u_d_x = u_d_x*(u_d_x*s_x / den);
	u_d_y = u_d_y*(u_d_y*s_y / den);
	double F_sz = (b*e*(a + X_s)*pow(2 * a - X_s, 2)) / (3 * a);

	double F_sx = u_d_x*F_sz;
	double F_sy = u_d_y*F_sz;
	double M_sz = -u_d_y*(b*e*X_s*X_s*pow(2 * a - X_s, 2)) / (4 * a);

	// Sum. of adhesion and sliding forces
	// ISO 규정 따름

	double sign_x, sign_y;
	fn_sign(slip_x, &sign_x);
	fn_sign(slip_y, &sign_y);


	double Fx = sign_x*(F_ax + F_sx);
	double Fy = -sign_y*(F_ay + F_sy);
	double Mz = -sign_y*(M_az + M_sz);

	// Output
	F[0] = Fx;
	F[1] = Fy;
	F[2] = Fz;
	M[0] = 0;
	M[1] = 0;
	M[2] = Mz;
}

void Obj_Tire::fn_R_eff(double *R_e, double *a, double R_d) {

	double e = acos(R_d / _R_u);
	*a = _R_u * sin(e);
	*R_e = *a / (e + _eps);
}

void Obj_Tire::fn_road_tire_A(double road_A[3][3], double tire_A[3][3], double road_z_Vec[3], double tire_y_Vec[3]) {

	// tire_Ay와 road_Az가 나란하면 에러
	double tire_Ax[3], tire_Ay[3], tire_Az[3];
	double road_Ax[3], road_Ay[3], road_Az[3];
	int i;

	memcpy(road_Az, road_z_Vec, sizeof(double) * 3);
	memcpy(tire_Ay, tire_y_Vec, sizeof(double) * 3);

	// 단위 벡터를 잘 입력해준다면 필요없음
	//fn_normalize(road_Az,3);
	//fn_normalize(tire_Ay,3); 

	// 연산
	fn_cross(road_Ax, tire_Ay, road_Az);
	fn_cross(road_Ay, road_Az, road_Ax);

	memcpy(tire_Ax, road_Ax, sizeof(double) * 3);

	fn_cross(tire_Az, tire_Ax, tire_Ay);

	// output
	for (i = 0; i < 3; i++) {
		road_A[i][0] = road_Ax[i]; //road의 좌표계
		road_A[i][1] = road_Ay[i];
		road_A[i][2] = road_Az[i];
		tire_A[i][0] = tire_Ax[i]; //tire의 좌표계
		tire_A[i][1] = tire_Ay[i];
		tire_A[i][2] = tire_Az[i];
	}
}

//void Obj_Tire::fn_slip_ratio(double *slip, double *angle, double V_x, double omega, double R_e, double V_sy, double Fz) {
//	double speed_m_x, speed_x, slip_x, speed_m_y, speed_y, slip_y, temp1, temp2;
//
//	speed_m_x = _CSLIP * R_e * R_e * _t_s / (_Iyy);  // 원래 슬립식에서 Euler method 사용시 unstable 조건의 속도
//	fn_max(fabs(V_x), speed_m_x, &speed_x);
//	slip_x = (R_e*omega - V_x) / speed_x;
//
//	speed_m_y = _Ca*9.81*_t_s / (2 * Fz + _eps);
//	fn_max(fabs(V_x), speed_m_y, &speed_y);
//	slip_y = V_sy / speed_y;
//	fn_max(slip_y, -1.0, &temp1);
//	fn_min(temp1, 1.0, &slip_y);
//
//	fn_max(slip_x, -1.0, &temp2);
//	fn_min(temp2, 1.0, slip);        // combined slip의 의미를 고려해서 - 1, +1로 제한
//	*angle = atan(slip_y); // atan(1) = pi / 4 = 45deg 이기 때문에 45도가 최대
//}

void Obj_Tire::fn_slip_ratio(double *slip_x, double *slip_y, double V_x, double omega, double R_e, double V_sy, double Fz, double a) {
	double speed_m_x, speed_x, speed_m_y, speed_y, temp1, temp2;
		
	_CSLIP = 2 * a*a * _c_x; // N / rad : cornering stiffness
	_Ca = 2 * a*a * _c_y;    // N / m : longitudinal stiffness

	speed_m_x = _CSLIP * R_e * R_e * _t_s / (2*_Iyy);  // 원래 슬립식에서 Euler method 사용시 unstable 조건의 속도
	fn_max(fabs(V_x), speed_m_x, &speed_x);
	*slip_x = (R_e*omega - V_x) / speed_x;

	speed_m_y = _Ca*9.81*_t_s / (2 * Fz + _eps);
	fn_max(fabs(V_x), speed_m_y, &speed_y);
	*slip_y = V_sy / speed_y;


	fn_max(*slip_x, -1.0, &temp2);
	fn_min(temp2, 1.0, slip_x);        // combined slip의 의미를 고려해서 - 1, +1로 제한

	fn_max(*slip_y, -1.0, &temp1);
	fn_min(temp1, 1.0, slip_y);
}

void Obj_Tire::fn_Fiala_tire(double F[3], double M[3], double Fz, double slip, double alpha, double omega)
{
	//  _mus % mu max
	//  _mud % mu min
	//  _w   % 0.235;
	//  _rr  % m : rolling resistance
	//  _Ca  % N/rad : cornering stiffness
	//  _CSLIP % N/m : longitudinal stiffness
	double  Fx, Fy, Mzo, My, alpha_cr, sign1, sign2, sign3;
	double sl_a, mu, Scr, Fx1, Fx2, H;

	// 조건에 따른 Longitudinal force 계산 %%%%
	if (Fz > 0) {
		sl_a = sqrt(pow(slip, 2) + pow(tan(alpha), 2));
		mu = _mus - (_mus - _mud) * sl_a;

		Scr = fabs(mu*Fz / (2 * _CSLIP));

		if (fabs(slip) < Scr) {// '='조건은 Fz >0 이라는 조건이 있기때문에 필요없음
			Fx = _CSLIP * slip; // C_Fk(k = 0일 때 Fx - k curve의 기울기) * longitudinal slip ratio(k)
		}
		else {
			Fx1 = mu*Fz;
			Fx2 = fabs((mu*Fz) * (mu*Fz) / (4 * fabs(slip)*_CSLIP));
			if (slip > 0) {
				Fx = (1)*(Fx1 - Fx2);
			}
			else if (slip < 0) {
				Fx = (-1)*(Fx1 - Fx2);
			}
			else {
				Fx = 0;
			}
		}
		// 조건에 따른 Lateral force & Aligning Moment 계산 %%%%
		alpha_cr = atan2(3 * mu*fabs(Fz), _Ca);
		if (fabs(alpha) <= alpha_cr) {
			H = 1 - (_Ca*fabs(tan(alpha)) / (3 * mu*fabs(Fz)));
			fn_sign(alpha, &sign1);
			Fy = -mu*fabs(Fz)*(1 - pow(H, 3))*sign1;

			fn_sign(alpha, &sign2);
			Mzo = mu*fabs(Fz)*_w*(1 - H)*pow(H, 3) * sign2;
		}
		else {

			fn_sign(alpha, &sign3);
			Fy = -mu*fabs(Fz)*sign3;
			Mzo = 0;
		}

		// 조건에 따른 Rolling Resistance Moment 계산 %%%%
		// t_r = _rr*_R_u;
		// fn_step_s(omega, _slope_rr, &step);
		// My = -Fz*t_r*step; // 기존:My = -_rr * Fz, 이렇게 되면 가만히 있어도, 뒤로가도 역토크가 걸림
		My = 0.0;
	}
	else {        // -- - Vertical force == 0,
		Fx = 0; Fy = 0; Fz = 0;
		My = 0; Mzo = 0;
	}

	// output
	F[0] = Fx;
	F[1] = Fy;
	F[2] = Fz;
	M[0] = 0;
	M[1] = My;
	M[2] = Mzo;
}

void Obj_Tire::fn_step_s(double x, double slope, double *output) {
	// 기울기가 정해진 step, -1, 0, 1
	// slope은 부호의미 없음;
	double xc, y;

	slope = fabs(slope);
	xc = _pi / slope * 0.5;
	y = (1 - cos(slope*(x + xc))) - 1;

	if (x < -xc) { y = -1; }
	else if (x > xc) { y = 1; }
	else {}

	*output = y;
}

void Obj_Tire::fn_cross(double c[3], double a[3], double b[3])
{
	c[0] = a[1] * b[2] - a[2] * b[1];
	c[1] = a[2] * b[0] - a[0] * b[2];
	c[2] = a[0] * b[1] - a[1] * b[0];
}

void Obj_Tire::fn_mat3331(double c[3], double a[3][3], double b[3])
{
	//c(3x1)=A(3x3)*B(3x1)
	c[0] = a[0][0] * b[0] + a[0][1] * b[1] + a[0][2] * b[2];
	c[1] = a[1][0] * b[0] + a[1][1] * b[1] + a[1][2] * b[2];
	c[2] = a[2][0] * b[0] + a[2][1] * b[1] + a[2][2] * b[2];
}

void Obj_Tire::fn_mat33T31(double c[3], double a[3][3], double b[3])
{
	//c(3x1)=A(3x3)'*B(3x1)
	c[0] = a[0][0] * b[0] + a[1][0] * b[1] + a[2][0] * b[2];
	c[1] = a[0][1] * b[0] + a[1][1] * b[1] + a[2][1] * b[2];
	c[2] = a[0][2] * b[0] + a[1][2] * b[1] + a[2][2] * b[2];
}

void Obj_Tire::fn_sign(double a, double *output)
{
	double b;
	if (a > 0) { b = 1; }
	else if (a < 0) { b = -1; }
	else { b = 0; }

	*output = b;
}

void Obj_Tire::fn_max(double a, double b, double *output)
{
	double c;
	if (a >= b) { c = a; }
	else { c = b; }

	*output = c;
}

void Obj_Tire::fn_min(double a, double b, double *output)
{
	double c;
	if (a <= b) { c = a; }
	else { c = b; }

	*output = c;
}

 void Obj_Tire::fn_normalize(double a[], uint num)
 {
 	double tmp = 0;
 	uint i;
 	for (i = 0; i < num; ++i)
 	{
 		tmp = tmp + a[i] * a[i];
 	}
 
 	for (i = 0; i < num; ++i)
 	{
 		a[i] = a[i] / sqrt(tmp);
 	}
 }

realtime_dynamics_analysis::realtime_dynamics_analysis(int WaypointSize) {
	DID = new DynamicInputData;
	NumberofThreads = 7;

	DID->NumberofThreads = NumberofThreads;
	DID->WaypointSize = WaypointSize;
	DID->LocalPath = new double*[WaypointSize];
	DID->StabilityIndex = new int*[WaypointSize];
	DID->VehicleVelocity = new double*[WaypointSize];
	for (int i = 0; i < WaypointSize; i++) {
		DID->LocalPath[i] = new double[3];
		DID->StabilityIndex[i] = new int[NumberofThreads];
		DID->VehicleVelocity[i] = new double[NumberofThreads];
	}
	DID->RSM = new double*[WaypointSize];
	DID->PSM = new double*[WaypointSize];
	DID->LSM = new double*[WaypointSize];
	DID->VSM = new double*[WaypointSize];
	for (int i = 0; i < WaypointSize; i++) {
		DID->RSM[i] = new double[NumberofThreads];
		DID->PSM[i] = new double[NumberofThreads];
		DID->LSM[i] = new double[NumberofThreads];
		DID->VSM[i] = new double[NumberofThreads];
		for (int j = 0; j < NumberofThreads; j++) {
			DID->RSM[i][j] = 0;
			DID->PSM[i][j] = 0;
			DID->LSM[i][j] = 0;
			DID->VSM[i][j] = 0;
		}
	}

	sim_indx = 0;
	heading = 0;
	vx = 0;
	WP_size = WaypointSize;

	// cputime 측정 변수(실차 적용시 불필요)
	read_start1 = read_end1 = read_start2 = read_end2 = read_time = 0;
	equilibrium_start = equilibrium_end = equilibrium_time = 0;
	parallel_start = parallel_end = parallel_time = 0;
}

realtime_dynamics_analysis::~realtime_dynamics_analysis() {
	for (int i = 0; i < WP_size; i++) {
		delete[] DID->LocalPath[i], DID->StabilityIndex[i], DID->VehicleVelocity[i];
		delete[] DID->RSM[i];
		delete[] DID->PSM[i];
		delete[] DID->LSM[i];
		delete[] DID->VSM[i];
	}
	delete[] DID->LocalPath, DID->StabilityIndex, DID->VehicleVelocity;
	delete[] DID->RSM;
	delete[] DID->PSM;
	delete[] DID->LSM;
	delete[] DID->VSM;
	delete DID;
}

void realtime_dynamics_analysis::init(RINPUTDATA *RTT_input) {
	read_start1 = omp_get_wtime();


 //	FILE *path, *map;
 //	char map_data[256], path_data[256];
 //	sprintf(map_data, "simulation_results/data_sim_%d/elevation_map.csv", sim_indx);
	//fopen_s(&map, map_data, "w+");
 //	for (int i = 0; i < RTT_input->MapInfo_Row; i++) {
 //		for (int j = 0; j < RTT_input->MapInfo_Col; j++) {
 //			if (j == RTT_input->MapInfo_Col - 1) {
 //				fprintf(map, "%5.5f\n", RTT_input->ElevationData[RTT_input->MapInfo_Row - 1 - i][j]);
 //			}
 //			else {
 //				fprintf(map, "%5.5f,", RTT_input->ElevationData[RTT_input->MapInfo_Row - 1 - i][j]);
 //			}
 //		}
 //	}
	//fclose(map);

 //	sprintf(path_data, "simulation_results/data_sim_%d/path.csv", sim_indx);
	//fopen_s(&path, path_data, "w+");
 //	for (int i = 0; i < RTT_input->WaypointSize; i++) {
 //		fprintf(path, "%5.5f,%5.5f\n", RTT_input->LocalPath[i][0] - RTT_input->MapInfo_ReferenceNorth, RTT_input->LocalPath[i][1] - RTT_input->MapInfo_ReferenceEast);
 //	}
	//fclose(path);
		
	//////////// Map //////////
	//unsigned int n_rows = RTT_input->MapInfo_Row, n_cols = RTT_input->MapInfo_Col;
	//double gab, xy_c[2];

	//// 노면 데이터 간격, 노면 기준 위치 설정
	//gab = RTT_input->MapInfo_Resolution_Row;
	//xy_c[0] = 0;
	//xy_c[1] = RTT_input->MapInfo_Row*RTT_input->MapInfo_Resolution_Row;

	///////////// Map 전처리 //////
	//Obj_Map::PNU_set_map(RTT_input->ElevationData, n_rows, n_cols, gab, xy_c[0], xy_c[1]);

	///////////// Nan_interp //////
	//Obj_Nan::fn_nan_interp(RTT_input->ElevationData, n_rows, n_cols);

	//RTT_input->Velocity = sqrt(RTT_input->Velocity_East*RTT_input->Velocity_East + RTT_input->Velocity_North*RTT_input->Velocity_North + RTT_input->Velocity_Down*RTT_input->Velocity_Down);
	RTT_input->Velocity = RTT_input->Velocity_Longitude;

	/////////// Path //////////
	// WP 3중 배열 변수에서 구조체로 변경 - 2017.03.14 송하준
	DID->CurrentVelocity = RTT_input->Velocity;

	for (int i = 0; i < WP_size; i++) {
		DID->LocalPath[i][0] = RTT_input->LocalPath[i][1];// - RTT_input->MapInfo_ReferenceEast;
		DID->LocalPath[i][1] = RTT_input->LocalPath[i][0];// - RTT_input->MapInfo_ReferenceNorth;		
		DID->LocalPath[i][2] = RTT_input->LocalPath[i][2];

		for (int j = 0; j < NumberofThreads; j++) {
			DID->StabilityIndex[i][j] = 0;
			DID->VehicleVelocity[i][j] = RTT_input->LocalPath[i][2];
			DID->StabilityIndex[0][j] = 1;
			DID->StabilityIndex[1][j] = 1;
		}
	}

	//heading = atan2((WP[3][1][0] - WP[0][1][0]), (WP[3][0][0] - WP[0][0][0])); // waypoint 1, 3 사이의 벡터를 통해 heading direction 계산(임시 값)
	//heading = atan2((DID->LocalPath[1][1] - DID->LocalPath[0][1]), (DID->LocalPath[1][0] - DID->LocalPath[0][0]));
	heading = RTT_input->Attitude_Yaw;	// 입력 받은 heading direction을 RTT class에서 사용하는 변수에 저장
	vx = RTT_input->Velocity; // 입력 받은 차량의 현재 속도를 RTT class에서 사용하는 멤버 변수에 저장
	road_mu = RTT_input->road_mu;
	// 고도맵, 경로 데이터, 항법 정보의 시간 확인
	//double map_path_time = fabs(RTT_input->Map_Time - RTT_input->Path_Time);
	//double map_navi_time = fabs(RTT_input->Map_Time - RTT_input->Navi_Time);
	//double path_navi_time = fabs(RTT_input->Path_Time - RTT_input->Navi_Time);
	//double ref_time = 0.58;
	//if (map_path_time > ref_time || map_navi_time > ref_time || path_navi_time > ref_time) {
	//	sprintf(RTT_input->input_data_err, "RTT input data time not match!!\n");
	//	printf("%s", RTT_input->input_data_err);
	//	RTT_input->err_flag = 2;
	//	return;
	//}

	
	// 경로 데이터가 고도맵 데이터 범위 안에 있는지 확인
	//for (int i = 0; i < RTT_input->WaypointSize; i++) {
	//	if (RTT_input->LocalPath[i][1] < RTT_input->MapInfo_ReferenceEast || RTT_input->LocalPath[i][1] > RTT_input->MapInfo_ReferenceEast + RTT_input->MapInfo_Col*RTT_input->MapInfo_Resolution_Col) {
	//		sprintf(RTT_input->input_data_err, "The elevation map and path data do not match!!\n");
	//		printf("%s", RTT_input->input_data_err);
	//		RTT_input->err_flag = 1;
	//		return;
	//	}
	//	if (RTT_input->LocalPath[i][0] > RTT_input->MapInfo_ReferenceNorth || RTT_input->LocalPath[i][0] < RTT_input->MapInfo_ReferenceNorth - RTT_input->MapInfo_Row*RTT_input->MapInfo_Resolution_Row) {
	//		sprintf(RTT_input->input_data_err, "The elevation map and path data do not match!!\n");
	//		printf("%s", RTT_input->input_data_err);
	//		RTT_input->err_flag = 1;
	//		return;
	//	}
	//}
	RTT_input->err_flag = 0;

	read_end1 = omp_get_wtime();
}

// Local Path Planner
void realtime_dynamics_analysis::LPP(RINPUTDATA *RTT_input) {
	double init_position[2] = { RTT_input->Position_East , RTT_input->Position_North};		// RTT 차량 초기 위치
	double init_position_gap[2];		// 첫 번째 waypoint와 RTT 차량 현재 위치 사이의 거리 차
	double LocalPathFiltered[40][2] = { 0 };
	double WP_gap, position_ratio, x_d, y_d, e_l;
	int i2 = 2;

	for (int i = 0; i < WP_size; i++) {
		LocalPathFiltered[i][0] = DID->LocalPath[i][0];		// Filtering 시킬 original local path 복사
		LocalPathFiltered[i][1] = DID->LocalPath[i][1];
	}
	init_position_gap[0] = init_position[0] - DID->LocalPath[0][0];
	init_position_gap[1] = init_position[1] - DID->LocalPath[0][1];
	LocalPathFiltered[0][0] = init_position[0];
	LocalPathFiltered[0][1] = init_position[1];
	LocalPathFiltered[1][0] = DID->LocalPath[1][0] + init_position_gap[0];		// 두 번째 waypoint도 원래의 waypoint와 평행하도록 나타냄
	LocalPathFiltered[1][1] = DID->LocalPath[1][1] + init_position_gap[1];

	for (int i = 2; i < WP_size; i++) {	// IIR filter
		LocalPathFiltered[i][0] = 0.0077769953*DID->LocalPath[i][0] + 0.0155539906*DID->LocalPath[i - 1][0] + 0.0077769953*DID->LocalPath[i - 2][0] + 1.7355001605*LocalPathFiltered[i - 1][0] + -0.7666081417*LocalPathFiltered[i - 2][0];
		LocalPathFiltered[i][1] = 0.0077769953*DID->LocalPath[i][1] + 0.0155539906*DID->LocalPath[i - 1][1] + 0.0077769953*DID->LocalPath[i - 2][1] + 1.7355001605*LocalPathFiltered[i - 1][1] + -0.7666081417*LocalPathFiltered[i - 2][1];
	}
	
	for (int i = 1; i < WP_size; i++) {
		WP_gap = sqrt((DID->LocalPath[i2][0] - DID->LocalPath[i2 - 1][0])*(DID->LocalPath[i2][0] - DID->LocalPath[i2 - 1][0]) + (DID->LocalPath[i2][1] - DID->LocalPath[i2 - 1][1])*(DID->LocalPath[i2][1] - DID->LocalPath[i2 - 1][1]));
		position_ratio = (((LocalPathFiltered[i][0] - DID->LocalPath[i2 - 1][0])*(DID->LocalPath[i2][0] - DID->LocalPath[i2 - 1][0])) + ((LocalPathFiltered[i][1] - DID->LocalPath[i2 - 1][1])*(DID->LocalPath[i2][1] - DID->LocalPath[i2 - 1][1]))) / (WP_gap*WP_gap);
		x_d = DID->LocalPath[i2 - 1][0] + position_ratio*(DID->LocalPath[i2][0] - DID->LocalPath[i2 - 1][0]);    // desired x position(perpendicular point from the path)
		y_d = DID->LocalPath[i2 - 1][1] + position_ratio*(DID->LocalPath[i2][1] - DID->LocalPath[i2 - 1][1]);    // desired y position(perpendicular point from the path)
		e_l = sqrt((x_d - LocalPathFiltered[i][0])*(x_d - LocalPathFiltered[i][0]) + (y_d - LocalPathFiltered[i][1])*(y_d - LocalPathFiltered[i][1]));
		if (e_l < 0.05) {
			break;
		}
		while (position_ratio > 1){
			i2 = i2 + 1;
			WP_gap = sqrt((DID->LocalPath[i2][0] - DID->LocalPath[i2 - 1][0])*(DID->LocalPath[i2][0] - DID->LocalPath[i2 - 1][0]) + (DID->LocalPath[i2][1] - DID->LocalPath[i2 - 1][1])*(DID->LocalPath[i2][1] - DID->LocalPath[i2 - 1][1]));
			position_ratio = (((LocalPathFiltered[i][0] - DID->LocalPath[i2 - 1][0])*(DID->LocalPath[i2][0] - DID->LocalPath[i2 - 1][0])) + ((LocalPathFiltered[i][1] - DID->LocalPath[i2 - 1][1])*(DID->LocalPath[i2][1] - DID->LocalPath[i2 - 1][1]))) / (WP_gap*WP_gap);
		}	
		DID->LocalPath[i2][0] = LocalPathFiltered[i][0];
		DID->LocalPath[i2][1] = LocalPathFiltered[i][1];
	}
}

void realtime_dynamics_analysis::run(RINPUTDATA *RTT_input) {
	SIM sim; CHASSIS chassis; SUS sus[6]; CTRL ctrl; Obj_Tire tire;

	read_start2 = omp_get_wtime();
	//	입력 받은 heading direction & current velocity 데이터를 class 멤버 변수에 저장
	sim.heading = heading;
	sim.vx = vx;

	sim.main_count = RTT_input->sim_count;

	tire.set_mu(road_mu);

	//	input parameter 저장
	//	initial position & orientation 계산
	read_system(&sim); DID->step_size = sim.h;
	read_chassis(&sim, &chassis, sus, tire.PNU_get_unloaded_radius());

	//	제어기에서 사용하는 변수 초기값 설정
	read_control(&ctrl);

	//	State vector 초기 값 설정
	define_Y_vector(&sim, &chassis, sus);
	read_end2 = omp_get_wtime();

	// pre process - 평형 상태 해석 및 차량 초기 속도 계산
	equilibrium_start = omp_get_wtime();
	equilibrium_process(&sim, &chassis, sus, &ctrl, &tire);
	equilibrium_end = omp_get_wtime();

	// 평형 상태 해석에서 에러 발생 할 경우 주행 시뮬레이션을 수행하지 않음
	if (sim.singular_flag == 0 && sim.singular_flag2 == 0) {

		// parallel process - 7개 스레드에서 주행 시뮬레이션 수행(병렬처리)
		parallel_start = omp_get_wtime();
		omp_set_num_threads(7);
#pragma omp parallel
		{
#pragma omp for firstprivate(sim, chassis, sus, ctrl, tire)
			for (int i = 0; i < NumberofThreads; i++) {
				sim.thread_indx = i;	// Thread id
				parallel_process(&sim, &chassis, sus, &ctrl, &tire);
			}
		}

		// 특정 스레드만 디버깅 할 때 사용
// 	 	sim.thread_indx = 0;
// 	 	parallel_process(&sim, &chassis, sus, &ctrl, &tire);
	}
	else {
		RTT_input->err_flag = 3;
	}

	parallel_end = omp_get_wtime();

	read_time = ((read_end1 - read_start1) + (read_end2 - read_start2)) * 1000;
	equilibrium_time = (equilibrium_end - equilibrium_start) * 1000;
	parallel_time = (parallel_end - parallel_start) * 1000;
}

void realtime_dynamics_analysis::equilibrium_process(SIM *sim, CHASSIS *chassis, SUS sus[6], CTRL *ctrl, Obj_Tire *tire) {
	const int n = 31;
	double p_count, p_step_size;
	double tol_acc, acc_flag;

	sim->t_current = 0; // simulation time 초기화
	acc_flag = 1;		// 챠량 가속도 값 저장 변수 초기화
	tol_acc = 0.01;		// 평형 상태 판단 기준 값

	sim->intcount = 1;	// absh3 integrator variable

	// data print step size control을 위한 변수
	p_count = 0;			
	p_step_size = sim->h;

	sim->thread_indx = 3;		// 평형 상태 해석 결과를 저장할 Thread id. 0 ~ 6 사이 아무 값을 넣어도 무방
	sim->simulation_flag = 1;	// 평형 상태 해석 시뮬레이션을 구분하는 flag 변수

	// singularity error flag 초기화
	sim->singular_flag = 0;
	sim->singular_flag2 = 0;

	while (1) {
		// 동역학 해석 sub main function
		analysis(sim, chassis, sus, ctrl, tire);

		// system down을 방지하기 위해서 운동방정식의 singularity를 감지하고 해석을 종료
		if (sim->singular_flag || sim->singular_flag2) {
			printf("Singular error occured!\n");
			break;
		}

		// 차량 가속도 값 계산 (3축 선형 가속도 합)
		acc_flag = sqrt(sim->Yp[13] * sim->Yp[13] + sim->Yp[14] * sim->Yp[14] + sim->Yp[15] * sim->Yp[15]);

		// 차량 가속도가 기준 값 이하일 경우 시뮬레이션 종료
		if (acc_flag < tol_acc) break;

		// 해석이 너무 오래 걸릴 경우 강제 종료
		if (sim->t_current >= 5.0) {
			sim->singular_flag = 0;
			break;
		}

		// explicit integrator
		absh3(sim->h, n, &sim->intcount, sim->AW1, sim->AW, &sim->t_current, sim->Y, sim->Yp, sim->Y_next, &sim->t_next);

		// 적분기를 통해 계산한 next step 의 state 를 다음 시뮬레이션에서 사용하기 위해 값을 복사
		memcpy(sim->Y, sim->Y_next, sizeof(double)*n);

		// cmd 출력 및 data 저장, printf step size 설정에 따라 바뀜
#if !DynamicCPUTimeCheck
		if ((sim->t_current == 0) || (sim->t_current + (p_step_size*0.0001) >= p_step_size * p_count)) {
			p_count = p_count + 1.0;
			//save_data(sim, chassis, sus, ctrl);
			printf("Simulation time : %5.5f\tSimulation Index : %d\n", sim->t_current, sim_indx);
		}
#endif
		// simulation time update
		sim->t_current = sim->t_next;
	}

	if (sim->singular_flag == 0 && sim->singular_flag2 == 0) {

		// 평형 해석 종료 후 차량의 position 값을 저장
		memcpy(sim->Y_equil, sim->Y, sizeof(double) * 13);

		velocity_candidate(ctrl, sim, tire->PNU_get_mu());
	}
}

void realtime_dynamics_analysis::velocity_candidate(CTRL *ctrl, SIM *sim, double road_mu) {
	double a_acc_max = 4.0;			// 최대 가속도
	double a_dec_max = -4.0;		// 최대 감속도
	double a_dec_critical;			// 최대 임계 감속도: 평지 상태에서의 최대 감속 가능한 감속도
	double vd_max_critical;			// 최대 임계 속도: 차량이 주어진 길이의 구간에서 급정거하여 최소 주행 속도(vd_min_critical = 0.5 m/s)에 도달할 수 있는 속도
	double vd_min_critical = 0.5;	// 최소 임계 속도
	double vd_mission = 15;			// 임무 속도: 광역 경로 계획에서 주어진 최대 주행 허용 속도 명령
	DID->vd_mission = vd_mission;
	double vd_max_allowable;		// 속도 후보군에서 허용 가능한 최대 값 (최대 임계 속도 > 임무 속도 : 임무속도, 임무속도 > 최대 임계 속도 : 최대 임계 속도)
	double vd_min_allowable = 1.5;	// 속도 후보군에서 허용 가능한 최저 값	
	double vd_max_accessable;		// 주어진 시간 동안 최대 가속 했을 때 도달하는 속도
	double vd_min_accessable;		// 주어진 시간 동안 최대 감속 했을 때 도달하는 속도
	double vd_max;					// 속도 후보군의 최대 값
	double vd_min;					// 속도 후보군의 최소 값
	double vd_gap;					// 속도 후보군 간격
	double step_time_margin = 0.5;	// 최대 가,감속 시킬 시간 (0.5 = 노면 데이터가 최대 500ms 동안 안 들어왔을 때까지의 속도 프로파일을 커버하겠다는 의미임)
	double WP_distance = 0;			// waypoint 총 길이
	double *LocalPathHeight = new double[ctrl->WP_size];
	double marginal_distance = 0;
	double LocalPathSlopeAngle = 0;
	int marginal_distance_flag; 

	for (int i = 0; i < ctrl->WP_size; i++) {
		//Obj_Map::PNU_fn_map(DID->LocalPath[i][0], DID->LocalPath[i][1], &LocalPathHeight[i]);
		Obj_Map_Global::PNU_fn_map(DID->LocalPath[i][0], DID->LocalPath[i][1], &LocalPathHeight[i]);
	}
	
	// waypoint의 총 길이 계산
	marginal_distance_flag = 0;
	for (int i = 1; i < ctrl->WP_size; i++) {
		//WP_distance += sqrt((WP[i][0][0] - WP[i - 1][0][0])*(WP[i][0][0] - WP[i - 1][0][0]) + (WP[i][1][0] - WP[i - 1][1][0])*(WP[i][1][0] - WP[i - 1][1][0]));
		WP_distance += sqrt(pow(DID->LocalPath[i][0] - DID->LocalPath[i - 1][0], 2) + pow(DID->LocalPath[i][1] - DID->LocalPath[i - 1][1], 2));
		// 노면 경사도 계산
		if (WP_distance > sim->vx * step_time_margin && marginal_distance_flag == 0) {			
			marginal_distance = WP_distance;	// step_time_margin 시간 동안 등속 이동 가능한 거리
			LocalPathSlopeAngle = atan2((LocalPathHeight[i] - LocalPathHeight[0]), marginal_distance);
			marginal_distance_flag = 1;
		}
	}
	DID->LocalPathSlopeAngle = LocalPathSlopeAngle;
	

	// 최대 임계 속도 계산 : 주어진 waypoint 총 길이에서 차량이 풀 브레이킹을 했을 때 vd_min_critical에 도달 가능한 최대 제한 속도
	a_dec_critical = sim->g*(sin(LocalPathSlopeAngle) + road_mu*cos(LocalPathSlopeAngle));	// 경사가 적용된 상태에서의 최대 감속도
	//a_dec_critical = -sim->g*(sin(LocalPathSlopeAngle) - road_mu*cos(LocalPathSlopeAngle));
	//a_dec_critical = 2;
	if (vd_min_critical*vd_min_critical - 2 * a_dec_critical*WP_distance >= 0)
		vd_max_critical = sqrt(vd_min_critical*vd_min_critical - 2 * a_dec_critical*WP_distance);	// 최대 임계 속도
	else
		vd_max_critical = sim->vx + 2;

	// 속도 후보군에서 허용 가능한 최대 값	 계산(최대 임계 속도 > 임무 속도 : 임무속도, 임무속도 > 최대 임계 속도 : 최대 임계 속도)
	if (vd_max_critical > vd_mission)
		vd_max_allowable = vd_mission;
	else
		vd_max_allowable = vd_max_critical;

	// 차량 현재 속도 대비 500ms 동안 도달 가능한 최대 및 최소 속도를 속도 후보군으로 설정
	vd_max_accessable = sim->vx + a_acc_max * step_time_margin;	
	vd_min_accessable = sim->vx + a_dec_max * step_time_margin;

	// 속도 후보군 최대 값 설정
	if (vd_max_accessable > vd_max_allowable)
		vd_max = vd_max_allowable;
	else
		vd_max = vd_max_accessable;
	
	// 현재 속도가 임무 속도 보다 빠를 경우 최저 속도 후보군을 임무 속도 아래로 설정
	if (sim->vx > vd_max) {
		vd_min_accessable = vd_max + a_dec_max * step_time_margin;
	}
		
	// 속도 후보군 최소 값 설정
	if (vd_min_accessable < vd_min_allowable)
		vd_min = vd_min_allowable;
	else
		vd_min = vd_min_accessable;

	// 각 쓰레드에 속도 후보군 배정
	ctrl->v_d[6] = vd_max;				// 최대 속도 후보군
	vd_gap = (vd_max - vd_min) / (NumberofThreads-1);	// 속도 후보군 간격
	for (int i = 5; i >= 0; i--) {
		ctrl->v_d[i] = ctrl->v_d[i + 1] - vd_gap;	// 나머지 쓰레드에 속도 후보군 배정
	}
	delete[] LocalPathHeight;
}

void realtime_dynamics_analysis::parallel_process(SIM *sim, CHASSIS *chassis, SUS sus[6], CTRL *ctrl, Obj_Tire *tire) {
	const int n = 31;
	double p_count, p_step_size;
	double temp[3];

	sim->intcount = 1;	// absh3 integrator variable

	// data print step size control을 위한 변수
	p_count = 0;
	p_step_size = sim->h;

	sim->t_current = 0;	// simulation time 초기화

	sim->simulation_flag = 3;	// 주행 해석 시뮬레이션을 구분하는 flag 변수

	// 제어기 사용 변수 초기화
	ctrl->e_v_sum = 0;
	ctrl->M_d_hat = 0;
	ctrl->yaw_hat = 0;
	ctrl->WP_indx = 1;

	memset(sim->Y, 0, sizeof(double)*n);	// state vector 0으로 초기화
	
	// 주행 시뮬레이션에서 사용할 state vector의 position 부분을 평형 해석 결과로 구한 값을 대입
	memcpy(sim->Y, sim->Y_equil, sizeof(double) * 13);	

	// 각 스레드 별 차량의 desired velocity 와 그에 따른 wheel omega 초기값 계산
  	sim->Y[13] = ctrl->v_d[sim->thread_indx];
  	sim->Y[25] = ctrl->v_d[sim->thread_indx] / tire->PNU_get_unloaded_radius();
  	sim->Y[26] = ctrl->v_d[sim->thread_indx] / tire->PNU_get_unloaded_radius();
  	sim->Y[27] = ctrl->v_d[sim->thread_indx] / tire->PNU_get_unloaded_radius();
  	sim->Y[28] = ctrl->v_d[sim->thread_indx] / tire->PNU_get_unloaded_radius();
  	sim->Y[29] = ctrl->v_d[sim->thread_indx] / tire->PNU_get_unloaded_radius();
  	sim->Y[30] = ctrl->v_d[sim->thread_indx] / tire->PNU_get_unloaded_radius();

	//// 각 스레드 별 차량의 desired velocity 와 그에 따른 wheel omega 초기값을 차량의 현재 속도 기준으로 계산 (2017.03.26 홍효성 수정)
	//sim->Y[13] = sim->vx;
	//sim->Y[25] = sim->vx / tire->PNU_get_unloaded_radius();
	//sim->Y[26] = sim->vx / tire->PNU_get_unloaded_radius();
	//sim->Y[27] = sim->vx / tire->PNU_get_unloaded_radius();
	//sim->Y[28] = sim->vx / tire->PNU_get_unloaded_radius();
	//sim->Y[29] = sim->vx / tire->PNU_get_unloaded_radius();
	//sim->Y[30] = sim->vx / tire->PNU_get_unloaded_radius();

	// 지역 경로 데이터에 따른 차량 자세 변환으로 인한 속도 벡터의 변환 값 계산
	temp[0] = sim->Y[13]; temp[1] = sim->Y[14]; temp[2] = sim->Y[15];
	mat3331(chassis->A0, temp, sim->Y + 13);

	// 노면 에러 회피 알고리즘을 적용할 최대 횟수
	sim->peak_error_count_max = 50;
	sim->peak_error_count = 0;

	while (1) {
		// 동역학 해석 sub main function
		analysis(sim, chassis, sus, ctrl, tire);

		// system down을 방지하기 위해서 운동방정식의 singularity를 감지하고 해석을 종료
		if (sim->singular_flag || sim->singular_flag2) {
			printf("Singular error occured!\n");
			break;
		}

		// way point 주행 종료시 해석 종료 조건
		if (ctrl->WP_indx >= ctrl->WP_size - 1) break;

		// 주행 해석 중 차량 전복 시 해석 종료
		if (sim->overturn_flag == 1) break;

		// 노면 에러 횟수가 기준값 이상일 경우 주행 불가능 환경으로 인식하고 해석 종료
		if (sim->peak_error_count > 50) break;

		// 해석 시간이 너무 오래 걸릴 경우 해석 강제 종료
		if (sim->t_current >= 30.0) break;

		// explicit integrator
		absh3(sim->h, n, &sim->intcount, sim->AW1, sim->AW, &sim->t_current, sim->Y, sim->Yp, sim->Y_next, &sim->t_next);

		// 적분기를 통해 계산한 next step 의 state 를 다음 시뮬레이션에서 사용하기 위해 값을 복사
		memcpy(sim->Y, sim->Y_next, sizeof(double)*n);

		// cmd 출력 및 data 저장, printf step size 설정에 따라 바뀜
#if !DynamicCPUTimeCheck
		if ((sim->t_current == 0) || (sim->t_current + (p_step_size*0.0001) >= p_step_size * p_count)) {
			p_count = p_count + 1.0;
			save_data(sim, chassis, sus, ctrl);

			printf("%3.2f m/s\t Simulation time : %f\t WayPoint : %d\t Simulation Index : %d\n", ctrl->v_d[sim->thread_indx], sim->t_current, ctrl->WP_indx, sim_indx);
		}
#endif
		// simulation time update
		sim->t_current = sim->t_next;
	}

#if !DynamicCPUTimeCheck
	//for (int i = 0; i < 23; i++) {
	fclose(sim->fp[7]);
	fclose(sim->fp[16]);
	//}
#endif

}

void realtime_dynamics_analysis::analysis(SIM *sim, CHASSIS *chassis, SUS sus[6], CTRL *ctrl, Obj_Tire *tire) {

	Y2qdq(sim, chassis, sus);

	orientation_chassis(chassis);
	position_chassis(chassis);
	velocity_state_chassis(chassis);
	cartesian_velocity_chassis(chassis);
	mass_force_state_chassis(chassis, sim);

	for (int i = 0; i < 6; i++) {
		orientation_suspension(&sus[i], chassis);
		position_suspension(&sus[i], chassis);
		velocity_state_suspension(&sus[i], chassis);
		cartesian_velocity_suspension(&sus[i]);
	}

	// 경사도에 따른 법선 벡터 계산하는 루틴 - 부산대
	// Equilibrium state 해석의 경우 초기 노면의 높이 데이터를 계속 사용하므로 time = 0 일 때 한번만 수행
// 	if (sim->simulation_flag == 1) { 
// 		if (sim->t_current == 0) {
// 			tire->PNU_Pre_road(sus[RF].rw, sus[RM].rw, sus[RR].rw, sus[LF].rw, sus[LM].rw, sus[LR].rw);
// 		}
// 	}
// 	else {
// 		tire->PNU_Pre_road(sus[RF].rw, sus[RM].rw, sus[RR].rw, sus[LF].rw, sus[LM].rw, sus[LR].rw);
// 	}

	// 경사도에 따른 법선 벡터 계산하는 루틴 - 충남대
	// Equilibrium state 해석의 경우 초기 노면의 높이 데이터를 계속 사용하므로 time = 0 일 때 한번만 수행
	if (sim->simulation_flag == 1) {
		if (sim->t_current == 0) {
			pre_road(sus, sim, tire);
			double N_vec[3] = { 0,0,1 };
			tire->set_N_vec(N_vec);
		}
	}
	else {
		pre_road(sus, sim, tire);
	}

	for(int i = 0; i < 6; i++) {
		sus[i].id = i;
		mass_force_state_suspension(&sus[i], sim, tire, chassis);
		velocity_coupling(chassis, &sus[i]);
		effective_mass_force(&sus[i]);
	}
	
	if(sim->simulation_flag == 3) sim->overturn_flag = overturn_detect(chassis, sus);
	else sim->overturn_flag = 0;

	acceleration_state_chassis(chassis, sus, sim);

	// singularity error 발생 또는 차량의 전복 발생시 더이상 해석 하지 않음
	if (sim->singular_flag == 0 && sim->overturn_flag == 0) {
		cartesian_acceleration_chassis(chassis);

		for (int i = 0; i < 6; i++) {
			acceleration_state_suspension(chassis, &sus[i]);
			cartesian_acceleration_suspension(&sus[i]);
		}

		LP_control(chassis, sus, ctrl, sim);
		wheel_spin_dyn(sus, ctrl);

		mat33T31(chassis->A0, chassis->dr0c, chassis->dr0cp);
		mat33T31(chassis->A0, chassis->ddr0c, chassis->ddr0cp);
		mat33T31(chassis->A0, chassis->w0, chassis->w0p);

		dqddq2Yp(sim, chassis, sus);
	}
}

void realtime_dynamics_analysis::Y2qdq(SIM *sim, CHASSIS *chassis, SUS sus[6]) {
	memcpy(chassis->r0, sim->Y, sizeof(double) * 3);
	memcpy(chassis->e0, sim->Y + 3, sizeof(double) * 4);

	sus[RF].q1 = sim->Y[7];
	sus[RM].q1 = sim->Y[8];
	sus[RR].q1 = sim->Y[9];
	sus[LF].q1 = sim->Y[10];
	sus[LM].q1 = sim->Y[11];
	sus[LR].q1 = sim->Y[12];

	memcpy(chassis->dr0, sim->Y + 13, sizeof(double) * 3);
	memcpy(chassis->w0, sim->Y + 16, sizeof(double) * 3);

	sus[RF].dq1 = sim->Y[19];
	sus[RM].dq1 = sim->Y[20];
	sus[RR].dq1 = sim->Y[21];
	sus[LF].dq1 = sim->Y[22];
	sus[LM].dq1 = sim->Y[23];
	sus[LR].dq1 = sim->Y[24];

	sus[RF].w_wh = sim->Y[25];
	sus[RM].w_wh = sim->Y[26];
	sus[RR].w_wh = sim->Y[27];
	sus[LF].w_wh = sim->Y[28];
	sus[LM].w_wh = sim->Y[29];
	sus[LR].w_wh = sim->Y[30];
}

void realtime_dynamics_analysis::orientation_chassis(CHASSIS *chassis) {
	chassis->E0[0][0] = -chassis->e0[1];
	chassis->E0[0][1] = chassis->e0[0];
	chassis->E0[0][2] = -chassis->e0[3];
	chassis->E0[0][3] = chassis->e0[2];
	chassis->E0[1][0] = -chassis->e0[2];
	chassis->E0[1][1] = chassis->e0[3];
	chassis->E0[1][2] = chassis->e0[0];
	chassis->E0[1][3] = -chassis->e0[1];
	chassis->E0[2][0] = -chassis->e0[3];
	chassis->E0[2][1] = -chassis->e0[2];
	chassis->E0[2][2] = chassis->e0[1];
	chassis->E0[2][3] = chassis->e0[0];

	chassis->A0[0][0] = 2 * (chassis->e0[0] * chassis->e0[0] + chassis->e0[1] * chassis->e0[1] - 0.5);
	chassis->A0[0][1] = 2 * (chassis->e0[1] * chassis->e0[2] - chassis->e0[0] * chassis->e0[3]);
	chassis->A0[0][2] = 2 * (chassis->e0[1] * chassis->e0[3] + chassis->e0[0] * chassis->e0[2]);
	chassis->A0[1][0] = 2 * (chassis->e0[1] * chassis->e0[2] + chassis->e0[0] * chassis->e0[3]);
	chassis->A0[1][1] = 2 * (chassis->e0[0] * chassis->e0[0] + chassis->e0[2] * chassis->e0[2] - 0.5);
	chassis->A0[1][2] = 2 * (chassis->e0[2] * chassis->e0[3] - chassis->e0[0] * chassis->e0[1]);
	chassis->A0[2][0] = 2 * (chassis->e0[1] * chassis->e0[3] - chassis->e0[0] * chassis->e0[2]);
	chassis->A0[2][1] = 2 * (chassis->e0[2] * chassis->e0[3] + chassis->e0[0] * chassis->e0[1]);
	chassis->A0[2][2] = 2 * (chassis->e0[0] * chassis->e0[0] + chassis->e0[3] * chassis->e0[3] - 0.5);

	chassis->roll_ang = -atan2(-chassis->A0[2][1], chassis->A0[2][2]);
	chassis->pitch_ang = atan2(chassis->A0[2][0], sqrt(chassis->A0[2][1] * chassis->A0[2][1] + chassis->A0[2][2] * chassis->A0[2][2]));
	chassis->yaw_ang = atan2(chassis->A0[1][0], chassis->A0[0][0]);
}

void realtime_dynamics_analysis::position_chassis(CHASSIS *chassis) {
	// rho0 = A0*rho0p
	mat3331(chassis->A0, chassis->rho0p, chassis->rho0);

	// r0c = r0 + rho0
	for (int i = 0; i < 3; i++) {
		chassis->r0c[i] = chassis->r0[i] + chassis->rho0[i];
	}
}

void realtime_dynamics_analysis::velocity_state_chassis(CHASSIS *chassis) {
	double temp[3];

	tilde(chassis->r0, chassis->r0t);

	// Yh0 = [dr0 + r0t*w0; w0]
	mat3331(chassis->r0t, chassis->w0, temp);

	chassis->Yh0[0] = chassis->dr0[0] + temp[0];
	chassis->Yh0[1] = chassis->dr0[1] + temp[1];
	chassis->Yh0[2] = chassis->dr0[2] + temp[2];
	chassis->Yh0[3] = chassis->w0[0];
	chassis->Yh0[4] = chassis->w0[1];
	chassis->Yh0[5] = chassis->w0[2];
}

void realtime_dynamics_analysis::cartesian_velocity_chassis(CHASSIS *chassis) {
	int i, j;
	double temp1[4], temp2[3];

	// de0 = 0.5*E0*w0
	mat34T31(chassis->E0, chassis->w0, temp1);
	for (i = 0; i < 4; i++) {
		chassis->de0[i] = 0.5*temp1[i];
	}

	// T0 = [eye(3) -r0t; zeros(3) eye(3)]
	for (i = 0; i < 6; i++) {
		for (j = 0; j < 6; j++) {
			if (i == j) chassis->T0[i][j] = 1;
			else chassis->T0[i][j] = 0;
		}
	}
	for (i = 0; i < 3; i++) {
		for (j = 3; j < 6; j++) {
			chassis->T0[i][j] = -chassis->r0t[i][j - 3];
		}
	}

	// dr0c = dr0 + w0t*rho0
	tilde(chassis->w0, chassis->w0t);
	mat3331(chassis->w0t, chassis->rho0, temp2);
	for (i = 0; i < 3; i++) {
		chassis->dr0c[i] = chassis->dr0[i] + temp2[i];
	}
}

void realtime_dynamics_analysis::mass_force_state_chassis(CHASSIS *chassis, SIM *sim) {
	int i, j;
	double A0_C00[3][3], A0_C00_J0p[3][3];

	mat3333(chassis->A0, chassis->C00, A0_C00);
	mat3333(A0_C00, chassis->J0, A0_C00_J0p);
	mat3333T(A0_C00_J0p, A0_C00, chassis->J0c);

	tilde(chassis->r0c, chassis->r0ct);
	tilde(chassis->dr0c, chassis->dr0ct);
	tilde(chassis->w0, chassis->w0t);

	// Mh=[m*eye(3),-m*rct;m*rct,Jc-m*rct*rct]
	for (i = 0; i < 6; i++)
	{
		for (j = 0; j < 6; j++)
		{
			if (i < 3 && j < 3)
			{
				if (i == j)
					chassis->Mh0[i][j] = chassis->m0;
				else
					chassis->Mh0[i][j] = 0;
			}
			else if (i < 3 && j>2)
				chassis->Mh0[i][j] = -chassis->m0*chassis->r0ct[i][j - 3];
			else if (i > 2 && j < 3)
				chassis->Mh0[i][j] = chassis->m0*chassis->r0ct[i - 3][j];
			else
				chassis->Mh0[i][j] = chassis->J0c[i - 3][j - 3] - chassis->m0*(chassis->r0ct[i - 3][0] * chassis->r0ct[0][j - 3] + chassis->r0ct[i - 3][1] * chassis->r0ct[1][j - 3] + chassis->r0ct[i - 3][2] * chassis->r0ct[2][j - 3]);
		}
	}

	chassis->F0c[0] = 0;
	chassis->F0c[1] = 0;
	chassis->F0c[2] = chassis->m0 * sim->g;

	chassis->T0c[0] = 0;
	chassis->T0c[1] = 0;
	chassis->T0c[2] = 0;

	// Qh=[Fc+m*drct*w;Tc+rct*Fc+m*rct*drct*w-wt*Jc*w]
	for (i = 0; i < 6; i++)
	{
		if (i < 3)
			chassis->Qh0[i] = chassis->F0c[i] + chassis->m0*(chassis->dr0ct[i][0] * chassis->w0[0] + chassis->dr0ct[i][1] * chassis->w0[1] + chassis->dr0ct[i][2] * chassis->w0[2]);
		else
		{
			chassis->Qh0[i] = chassis->T0c[i - 3] + chassis->r0ct[i - 3][0] * chassis->F0c[0] + chassis->r0ct[i - 3][1] * chassis->F0c[1] + chassis->r0ct[i - 3][2] * chassis->F0c[2] +
				chassis->m0*((chassis->r0ct[i - 3][0] * chassis->dr0ct[0][0] + chassis->r0ct[i - 3][1] * chassis->dr0ct[1][0] + chassis->r0ct[i - 3][2] * chassis->dr0ct[2][0])*chassis->w0[0] +
				(chassis->r0ct[i - 3][0] * chassis->dr0ct[0][1] + chassis->r0ct[i - 3][1] * chassis->dr0ct[1][1] + chassis->r0ct[i - 3][2] * chassis->dr0ct[2][1])*chassis->w0[1] +
					(chassis->r0ct[i - 3][0] * chassis->dr0ct[0][2] + chassis->r0ct[i - 3][1] * chassis->dr0ct[1][2] + chassis->r0ct[i - 3][2] * chassis->dr0ct[2][2])*chassis->w0[2]) +
				-((chassis->w0t[i - 3][0] * chassis->J0c[0][0] + chassis->w0t[i - 3][1] * chassis->J0c[1][0] + chassis->w0t[i - 3][2] * chassis->J0c[2][0])*chassis->w0[0] +
				(chassis->w0t[i - 3][0] * chassis->J0c[0][1] + chassis->w0t[i - 3][1] * chassis->J0c[1][1] + chassis->w0t[i - 3][2] * chassis->J0c[2][1])*chassis->w0[1] +
					(chassis->w0t[i - 3][0] * chassis->J0c[0][2] + chassis->w0t[i - 3][1] * chassis->J0c[1][2] + chassis->w0t[i - 3][2] * chassis->J0c[2][2])*chassis->w0[2]);
		}
	}
}

void realtime_dynamics_analysis::orientation_suspension(SUS *sus, CHASSIS *chassis) {
	double A01pp[3][3] = { 0, }, u[3] = { 0, 0, 1 };

	A01pp[0][0] = cos(sus->q1);
	A01pp[0][1] = -sin(sus->q1);
	A01pp[0][2] = 0;
	A01pp[1][0] = sin(sus->q1);
	A01pp[1][1] = cos(sus->q1);
	A01pp[1][2] = 0;
	A01pp[2][0] = 0;
	A01pp[2][1] = 0;
	A01pp[2][2] = 1;

	mat333333(chassis->A0, sus->C01, A01pp, sus->A1);
	mat333331(chassis->A0, sus->C01, u, sus->H1);
}

void realtime_dynamics_analysis::position_suspension(SUS *sus, CHASSIS *chassis) {
	int i;
	double A0_s01p[3];

	// s01 = A0*s01p
	mat3331(chassis->A0, sus->s01p, A0_s01p);

	// r1 = r0 + A0*s01p
	for (i = 0; i < 3; i++) {
		sus->r1[i] = chassis->r0[i] + A0_s01p[i];
	}

	// s12 = A0*s12p
	mat3331(sus->A1, sus->s12p, sus->s12);

	// rw = r1 + s12
	for (i = 0; i < 3; i++) {
		sus->rw[i] = sus->r1[i] + sus->s12[i];
	}

	// rho1 = A1*rho1p
	mat3331(sus->A1, sus->rho1p, sus->rho1);

	// r1c = r1 + rho1
	for (i = 0; i < 3; i++) {
		sus->r1c[i] = sus->r1[i] + sus->rho1[i];
	}
}

void realtime_dynamics_analysis::velocity_state_suspension(SUS *sus, CHASSIS *chassis) {
	int i;
	double r1t_H1[3];

	tilde(sus->r1, sus->r1t);

	// B1 = [r1t*H1; H1]
	mat3331(sus->r1t, sus->H1, r1t_H1);
	for (i = 0; i < 3; i++) {
		sus->B1[i] = r1t_H1[i];
		sus->B1[i + 3] = sus->H1[i];
	}
	
	// Yh1 = Yh0 + B1*dq1
	for (i = 0; i < 6; i++) {
		sus->Yh1[i] = chassis->Yh0[i] + sus->B1[i] * sus->dq1;
	}
}

void realtime_dynamics_analysis::cartesian_velocity_suspension(SUS *sus) {
	double w1t_s12[3], w1t_rho1[3];
	int i, j;

	// T1 = [eye(3) -r1t; zeros(3) eye(3)]
	for (i = 0; i < 6; i++) {
		for (j = 0; j < 6; j++) {
			if (i == j) sus->T1[i][j] = 1;
			else sus->T1[i][j] = 0;
		}
	}
	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			sus->T1[i][j + 3] = -sus->r1t[i][j];
		}
	}

	// Yb1 = T1*Yh1
	for (i = 0; i < 6; i++) {
		sus->Yb1[i] = 0;
		for (j = 0; j < 6; j++) {
			sus->Yb1[i] += sus->T1[i][j] * sus->Yh1[j];
		}
	}

	// dr1 = Yb1(1:3)
	// w1 = Yb1(4:6)
	for (i = 0; i < 3; i++) {
		sus->dr1[i] = sus->Yb1[i];
		sus->w1[i] = sus->Yb1[i + 3];
	}

	tilde(sus->w1, sus->w1t);

	// dr1c = dr1 + w1t*rho1
	// drwc = dr1 + w1t*s12
	mat3331(sus->w1t, sus->rho1, w1t_rho1);
	mat3331(sus->w1t, sus->s12, w1t_s12);
	for (i = 0; i < 3; i++) {
		sus->dr1c[i] = sus->dr1[i] + w1t_rho1[i];
		sus->drwc[i] = sus->dr1[i] + w1t_s12[i];
	}
}

void realtime_dynamics_analysis::pre_road(SUS sus[6], SIM *sim, Obj_Tire *tire) {
	double road_h[6];

	for (int i = 0; i < 6; i++) {
		//Obj_Map::PNU_fn_map(sus[i].rw[0], sus[i].rw[1], &road_h[i]);
		Obj_Map_Global::PNU_fn_map(sus[i].rw[0], sus[i].rw[1], &road_h[i]);
	}

	double P_RF[3] = { sus[RF].rw[0], sus[RF].rw[1], road_h[0] };
	double P_RR[3] = { sus[RR].rw[0], sus[RR].rw[1], road_h[2] };
	double P_LF[3] = { sus[LF].rw[0], sus[LF].rw[1], road_h[3] };
	double P_LR[3] = { sus[LR].rw[0], sus[LR].rw[1], road_h[5] };

	double T_vec1[3], T_vec2[3], T_vec1t[3][3], N_vec[3];
	for (int i = 0; i < 3; i++) {
		T_vec1[i] = P_RF[i] - P_LR[i];
		T_vec2[i] = P_LF[i] - P_RR[i];
	}

	tilde(T_vec1, T_vec1t);
	mat3331(T_vec1t, T_vec2, N_vec);
	// vector normalization
	double tmp = 0;
	for (int i = 0; i < 3; ++i) tmp = tmp + N_vec[i] * N_vec[i];
	for (int i = 0; i < 3; ++i) N_vec[i] = N_vec[i] / sqrt(tmp);

	tire->set_Road_h(road_h);
	tire->set_N_vec(N_vec);

	memcpy(sim->old_road_h, sim->pre_road_h, sizeof(double) * 6);
	memcpy(sim->pre_road_h, road_h, sizeof(double) * 6);
}

int realtime_dynamics_analysis::peak_error_detect(double current_road[6], double previous_road[6], double old_road[6]) {
	// 한 스텝 전, 두 스텝 전 데이터 사용
	double peak_error_ref = 0.31;
	double pre_road_ave = 0, old_road_ave = 0;
	int peak_error_flag = 0;

	for (int i = 0; i < 6; i++) {
		pre_road_ave += previous_road[i];
		old_road_ave += old_road[i];
	}
	pre_road_ave /= 6.0;
	old_road_ave /= 6.0;

	for (int i = 0; i < 6; i++) {
		if (current_road[i] - previous_road[i] >= peak_error_ref) {
			current_road[i] = pre_road_ave + (pre_road_ave - old_road_ave);
			peak_error_flag = 1;
		}
	}

	return peak_error_flag;

	// 이전 스텝 데이터 사용
	// 	double peak_error_ref = 0.31;
	// 	double delta_road[6], delta_ave = 0;
	// 	double count = 0;
	// 	double pre_road_ave = 0;
	// 
	// 	for (int i = 0; i < 6; i++) {
	// 		delta_road[i] = current_road[i] - previous_road[i];
	// 		if (delta_road[i] < peak_error_ref) {
	// 			delta_ave += delta_road[i];
	// 			count = count + 1.0;
	// 		}
	// 	}
	// 	delta_ave /= count;
	// 	for (int i = 0; i < 6; i++) pre_road_ave += previous_road[i];
	// 	pre_road_ave /= 6.0;
	// 	for (int i = 0; i < 6; i++) {
	// 		if (delta_road[i] >= peak_error_ref) {
	// 			current_road[i] = pre_road_ave + delta_ave;
	// 		}
	// 	}
}

void realtime_dynamics_analysis::mass_force_state_suspension(SUS *sus, SIM *sim, Obj_Tire *tire, CHASSIS *chassis) {
	double A1C11[3][3], A1_C11_J1p[3][3], s12t_rho1t_tire[3], s12t_rho1t[3][3];
	int i, j;

	// J1c = A1*C11*J1*(A1*C11)'
	mat3333(sus->A1, sus->C11, A1C11);
	mat3333(A1C11, sus->J1, A1_C11_J1p);
	mat3333T(A1_C11_J1p, A1C11, sus->J1c);

	tilde(sus->r1c, sus->r1ct);
	tilde(sus->dr1c, sus->dr1ct);

	// Mh=[m*eye(3),-m*rct;m*rct,Jc-m*rct*rct]
	for (i = 0; i < 6; i++) {
		for (j = 0; j < 6; j++) {
			if (i < 3 && j < 3) {
				if (i == j)
					sus->Mh1[i][j] = sus->m1;
				else
					sus->Mh1[i][j] = 0;
			}
			else if (i < 3 && j>2)
				sus->Mh1[i][j] = -sus->m1*sus->r1ct[i][j - 3];
			else if (i > 2 && j < 3)
				sus->Mh1[i][j] = sus->m1*sus->r1ct[i - 3][j];
			else
				sus->Mh1[i][j] = sus->J1c[i - 3][j - 3] - sus->m1*(sus->r1ct[i - 3][0] * sus->r1ct[0][j - 3] + sus->r1ct[i - 3][1] * sus->r1ct[1][j - 3] + sus->r1ct[i - 3][2] * sus->r1ct[2][j - 3]);
		}
	}

	// TSDA force 계산
	CNU_TSDA(sus, chassis, sim);

	// Tire input parameter - wheel 회전 축 vector
	double temp[3][3] = { 0, };
	mat3333(sus->A1, sus->C12, temp);
	double u_vec[3] = { -temp[0][2], -temp[1][2], -temp[2][2] };

	// Tire force & moment calculate
	tire->PNU_Tire_force(sus->rw, sus->drwc, u_vec, sus->w_wh, sus->id);
	// Get result of tire analysis
	tire->PNU_get_data(sus->F_tire, sus->M_tire, &sus->slip, &sus->angle, &sus->road_h, &sus->pen, &sus->R_d, &sus->Fx, &sus->Fy, &sus->Fz, &sus->My);
	
	// 타이어 침투량 변화 비율, 타이어 침투량 변화량, tire의 point follow 점, 이전 스텝 타이어 침투량의 평균
//	double pen_ratio = 0, delta_pen = 0, point_h = 0, pen_ave = 0;
//
// 	switch (sim->simulation_flag) {
// 		case 1:
// 			// Tire force & moment calculate
// 			tire->PNU_Tire_force(sus->rw, sus->drwc, u_vec, sus->w_wh, sus->id);
// 			// Get result of tire analysis
// 			tire->PNU_get_data(sus->F_tire, sus->M_tire, &sus->slip, &sus->angle, &sus->road_h, &sus->pen, &sus->R_d, &sus->Fx, &sus->Fy, &sus->Fz, &sus->My);
// 			sim->pre_pen[sus->id] = sus->pen;
// 			break;
// 		case 3:
// 			// 타이어 침투량 계산
// 			Obj_Map::PNU_fn_map(sus->rw[0], sus->rw[1], &sus->road_h);
// 			point_h = sus->rw[2] - tire->PNU_get_unloaded_radius();
// 
// 			if (sus->road_h - point_h >= 0.0) {
// 				sus->pen = sus->road_h - point_h;
// 			}
// 			else {
// 				sus->pen = 0;
// 			}
// 			
// 			// 노면의 peak error를 감지하는 루틴
// 			// equilibrium state 해석시에는 사용하지 않고 주행 시뮬레이션시에만 사용
// 
// 			// 이전 스템의 타이어 침투량 평균값 계산
// 			for (int i = 0; i < 6; i++) pen_ave += sim->pre_pen[i];
// 			pen_ave /= 6;
// 
// 			// 타이어 침투량 변화량 및 침투량 변화율 계산
// 			delta_pen = sus->pen - pen_ave;
// 			pen_ratio = delta_pen / pen_ave;
// 
// 			// 노면 peak error일 경우 침투량을 이전 스텝 침투량의 평균값으로 변경
// 			if (pen_ratio > sim->peak_err_ref_ratio) {
// 				sus->pen = pen_ave;
// 			}
// 
// 			// Tire force & moment calculate
// 			tire->PNU_Tire_force(sus->rw, sus->drwc, u_vec, sus->w_wh, sus->id, sus->pen);
// 			// Get result of tire analysis
// 			tire->PNU_get_data(sus->F_tire, sus->M_tire, &sus->slip, &sus->angle, &sus->road_h, &sus->pen, &sus->R_d, &sus->Fx, &sus->Fy, &sus->Fz, &sus->My);
// 
// 			// 타이어의 index에 맞춰 현재 스텝의 타이어 침투량 저장
// 			sim->pre_pen[sus->id] = sus->pen;
// 
// 			break;
// 		default:
// 			break;
// 	}

	// 현가장치의 질량 중심점에 작용하는 힘. 타이어력에 의한 효과와 중력에 의한 효과의 합
	sus->F1c[0] = sus->F_tire[0];
	sus->F1c[1] = sus->F_tire[1];
	sus->F1c[2] = sus->F_tire[2] + sus->m1*sim->g;

	// T1c = (s12t - rho1t)*F_tire + M_tire
	// 현가장치의 질량 중심점에 작용하는 모멘트. 타이어력으로 인해 발생하는 모멘트와 중력으로 인해 발생하는 모멘트의 합.
	tilde(sus->s12, sus->s12t);
	tilde(sus->rho1, sus->rho1t);
	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			s12t_rho1t[i][j] = sus->s12t[i][j] - sus->rho1t[i][j];
		}
	}
	mat3331(s12t_rho1t, sus->F_tire, s12t_rho1t_tire);
	sus->T1c[0] = s12t_rho1t_tire[0] + sus->M_tire[0];
	sus->T1c[1] = s12t_rho1t_tire[1] + sus->M_tire[1];
	sus->T1c[2] = s12t_rho1t_tire[2] + sus->M_tire[2];

	// Qh=[Fc+m*drct*w;Tc+rct*Fc+m*rct*drct*w-wt*Jc*w]
	for (i = 0; i < 6; i++)
	{
		if (i < 3)
			sus->Qh1[i] = sus->F1c[i] + sus->m1*(sus->dr1ct[i][0] * sus->w1[0] + sus->dr1ct[i][1] * sus->w1[1] + sus->dr1ct[i][2] * sus->w1[2]);
		else
			sus->Qh1[i] = sus->T1c[i - 3] + sus->r1ct[i - 3][0] * sus->F1c[0] + sus->r1ct[i - 3][1] * sus->F1c[1] + sus->r1ct[i - 3][2] * sus->F1c[2] +
			sus->m1*((sus->r1ct[i - 3][0] * sus->dr1ct[0][0] + sus->r1ct[i - 3][1] * sus->dr1ct[1][0] + sus->r1ct[i - 3][2] * sus->dr1ct[2][0])*sus->w1[0] +
			(sus->r1ct[i - 3][0] * sus->dr1ct[0][1] + sus->r1ct[i - 3][1] * sus->dr1ct[1][1] + sus->r1ct[i - 3][2] * sus->dr1ct[2][1])*sus->w1[1] +
				(sus->r1ct[i - 3][0] * sus->dr1ct[0][2] + sus->r1ct[i - 3][1] * sus->dr1ct[1][2] + sus->r1ct[i - 3][2] * sus->dr1ct[2][2])*sus->w1[2]) +
			-((sus->w1t[i - 3][0] * sus->J1c[0][0] + sus->w1t[i - 3][1] * sus->J1c[1][0] + sus->w1t[i - 3][2] * sus->J1c[2][0])*sus->w1[0] +
			(sus->w1t[i - 3][0] * sus->J1c[0][1] + sus->w1t[i - 3][1] * sus->J1c[1][1] + sus->w1t[i - 3][2] * sus->J1c[2][1])*sus->w1[1] +
				(sus->w1t[i - 3][0] * sus->J1c[0][2] + sus->w1t[i - 3][1] * sus->J1c[1][2] + sus->w1t[i - 3][2] * sus->J1c[2][2])*sus->w1[2]);
	}
}

void realtime_dynamics_analysis::CNU_TSDA(SUS *sus, CHASSIS *chassis, SIM *sim) {
	double r1t_s1st_d01[3], r1t_s1st[3][3], r0t_s0st_d01[3], r0t_s0st[3][3], f_L, s1st[3][3], s0st[3][3];
	double dL_s, w0t_A0_s0sp[3], w1t_A1_s1sp[3], p8, p7, p6, p5, p4, p3, p2, p1, r1s[3], r0s[3], s1s[3], s0s[3], L_free;
	int i, j;
	double temp1[3], temp2[3];

	L_free = 0.8;	// Spring initial length

	mat3331(chassis->A0, sus->s0sp, s0s);	// chassis의 스프링 Global 위치 벡터 계산
	mat3331(sus->A1, sus->s1sp, s1s);		// 현가장치의 스프링 Global 위치 벡터 계산

	// d01 = r1 + s1s + r0 + s0s
	// 스프링 길이 계산
	for (i = 0; i < 3; i++) {
		r0s[i] = chassis->r0[i] + s0s[i];
		r1s[i] = sus->r1[i] + s1s[i];
		sus->d01[i] = r1s[i] - r0s[i];
	}
	sus->L_spring = sqrt(sus->d01[0] * sus->d01[0] + sus->d01[1] * sus->d01[1] + sus->d01[2] * sus->d01[2]);

	// 스프링 변위
	sus->defo = sus->L_spring - L_free; // deformation
	// Spring curve fitting parameter
	p1 = 4.631e+09;
	p2 = -1.492e-06;
	p3 = -2.063e+08;
	p4 = 6.644e-08;
	p5 = 3.638e+06;
	p6 = -6.597e-10;
	p7 = 2279;
	p8 = -3000;
	// Spring force
	sus->T_spring = p1*pow(sus->defo, 7) + p2*pow(sus->defo, 6) + p3*pow(sus->defo, 5) + p4*pow(sus->defo, 4) + p5*pow(sus->defo, 3) + p6*pow(sus->defo, 2) + p7*sus->defo + p8;

	// 스프링 속도
	temp1[0] = sus->d01[0] / sus->L_spring;
	temp1[1] = sus->d01[1] / sus->L_spring;
	temp1[2] = sus->d01[2] / sus->L_spring;
	mat333331(sus->w1t, sus->A1, sus->s1sp, w1t_A1_s1sp);
	mat333331(chassis->w0t, chassis->A0, sus->s0sp, w0t_A0_s0sp);
	for (i = 0; i < 3; i++) {
		temp2[i] = sus->dr1[i] + w1t_A1_s1sp[i] - chassis->dr0[i] - w0t_A0_s0sp[i];
	}
	dL_s = temp1[0] * temp2[0] + temp1[1] * temp2[1] + temp1[2] * temp2[2];

	// Equilibrium state 해석 시에는 linear damper 사용. 주행 시뮬레이션 시에는 non-linear damper 사용
	if (sim->simulation_flag == 1) {
		sus->T_damper = 31500 * 2.0 * dL_s; // 상수 댐핑력
	}
	else {
		// 비선형 1.5A 일 때 - 1차 방정식 + 탄젠트함수
		double b5 = 1.075e+04;
		double b6 = 14.37;
		sus->T_damper = (b5*dL_s + b6) + 9300 * (2 / M_PI)*atan2(10 * dL_s, 1); // 1.5A 일때
	}

	// Spring force 와 damping force 를 합한 TSDA force 계산
	sus->f = sus->T_spring + sus->T_damper;

	// TSDA force로 인해 발생하는 Qh 계산
	tilde(s0s, s0st);
	tilde(s1s, s1st);

	// Qh0_TSDA = (f/L_spring)*[d01; (r0st+s0st)*d01]
	// Qh1_TSDA = -(f/L_spring)*[d01; (r1st+s1st)*d01]
	// Qh0_TSDA와 Qh1_TSDA는 크기는 같고 방향만 반대
	f_L = sus->f / sus->L_spring;
	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			r0t_s0st[i][j] = chassis->r0t[i][j] + s0st[i][j];
			r1t_s1st[i][j] = sus->r1t[i][j] + s1st[i][j];
		}
	}
	mat3331(r0t_s0st, sus->d01, r0t_s0st_d01);
	mat3331(r1t_s1st, sus->d01, r1t_s1st_d01);

	sus->Qh0_TSDA[0] = f_L*sus->d01[0];
	sus->Qh0_TSDA[1] = f_L*sus->d01[1];
	sus->Qh0_TSDA[2] = f_L*sus->d01[2];
	sus->Qh0_TSDA[3] = f_L*r0t_s0st_d01[0];
	sus->Qh0_TSDA[4] = f_L*r0t_s0st_d01[1];
	sus->Qh0_TSDA[5] = f_L*r0t_s0st_d01[2];

	sus->Qh1_TSDA[0] = -f_L*sus->d01[0];
	sus->Qh1_TSDA[1] = -f_L*sus->d01[1];
	sus->Qh1_TSDA[2] = -f_L*sus->d01[2];
	sus->Qh1_TSDA[3] = -f_L*r1t_s1st_d01[0];
	sus->Qh1_TSDA[4] = -f_L*r1t_s1st_d01[1];
	sus->Qh1_TSDA[5] = -f_L*r1t_s1st_d01[2];
}

void realtime_dynamics_analysis::velocity_coupling(CHASSIS *chassis, SUS *sus) {
	// dH1 = w0t*H1
	tilde(chassis->dr0, chassis->dr0t);
	tilde(sus->dr1, sus->dr1t);
	mat3331(chassis->w0t, sus->H1, sus->dH1);

	// D1=[dr1t*H1+r1t*dH1;dH1]*dq1
	sus->D1[0] = (sus->dr1t[0][0] * sus->H1[0] + sus->dr1t[0][1] * sus->H1[1] + sus->dr1t[0][2] * sus->H1[2] + sus->r1t[0][0] * sus->dH1[0] + sus->r1t[0][1] * sus->dH1[1] + sus->r1t[0][2] * sus->dH1[2])*sus->dq1;
	sus->D1[1] = (sus->dr1t[1][0] * sus->H1[0] + sus->dr1t[1][1] * sus->H1[1] + sus->dr1t[1][2] * sus->H1[2] + sus->r1t[1][0] * sus->dH1[0] + sus->r1t[1][1] * sus->dH1[1] + sus->r1t[1][2] * sus->dH1[2])*sus->dq1;
	sus->D1[2] = (sus->dr1t[2][0] * sus->H1[0] + sus->dr1t[2][1] * sus->H1[1] + sus->dr1t[2][2] * sus->H1[2] + sus->r1t[2][0] * sus->dH1[0] + sus->r1t[2][1] * sus->dH1[1] + sus->r1t[2][2] * sus->dH1[2])*sus->dq1;
	sus->D1[3] = sus->dH1[0] * sus->dq1;
	sus->D1[4] = sus->dH1[1] * sus->dq1;
	sus->D1[5] = sus->dH1[2] * sus->dq1;
}

void realtime_dynamics_analysis::effective_mass_force(SUS *sus) {
	double temp[6], Mqq, Myy[6][6];
	int i, j;

	// Myy=Mh1
	for (i = 0; i < 6; i++) {
		for (j = 0; j < 6; j++) {
			Myy[i][j] = sus->Mh1[i][j];
		}
	}

	// Myq=Mh1'*B1
	for (i = 0; i < 6; i++) {
		sus->Myq[i] = 0;
		for (j = 0; j < 6; j++) {
			sus->Myq[i] += sus->Mh1[j][i] * sus->B1[j];
		}
	}

	// Mqq=B1'*Mh1*B1
	Mqq = sus->B1[0] * sus->Myq[0] + sus->B1[1] * sus->Myq[1] + sus->B1[2] * sus->Myq[2] + sus->B1[3] * sus->Myq[3] + sus->B1[4] * sus->Myq[4] + sus->B1[5] * sus->Myq[5];


	// Pq=B1'*(Qh1+Qh1_TSDA-Mh1*D1)
	sus->Pq = sus->B1[0] * (-sus->Mh1[0][0] * sus->D1[0] - sus->Mh1[0][1] * sus->D1[1] - sus->Mh1[0][2] * sus->D1[2] - sus->Mh1[0][3] * sus->D1[3] - sus->Mh1[0][4] * sus->D1[4] - sus->Mh1[0][5] * sus->D1[5] + sus->Qh1[0] + sus->Qh1_TSDA[0]) +
		sus->B1[1] * (-sus->Mh1[1][0] * sus->D1[0] - sus->Mh1[1][1] * sus->D1[1] - sus->Mh1[1][2] * sus->D1[2] - sus->Mh1[1][3] * sus->D1[3] - sus->Mh1[1][4] * sus->D1[4] - sus->Mh1[1][5] * sus->D1[5] + sus->Qh1[1] + sus->Qh1_TSDA[1]) +
		sus->B1[2] * (-sus->Mh1[2][0] * sus->D1[0] - sus->Mh1[2][1] * sus->D1[1] - sus->Mh1[2][2] * sus->D1[2] - sus->Mh1[2][3] * sus->D1[3] - sus->Mh1[2][4] * sus->D1[4] - sus->Mh1[2][5] * sus->D1[5] + sus->Qh1[2] + sus->Qh1_TSDA[2]) +
		sus->B1[3] * (-sus->Mh1[3][0] * sus->D1[0] - sus->Mh1[3][1] * sus->D1[1] - sus->Mh1[3][2] * sus->D1[2] - sus->Mh1[3][3] * sus->D1[3] - sus->Mh1[3][4] * sus->D1[4] - sus->Mh1[3][5] * sus->D1[5] + sus->Qh1[3] + sus->Qh1_TSDA[3]) +
		sus->B1[4] * (-sus->Mh1[4][0] * sus->D1[0] - sus->Mh1[4][1] * sus->D1[1] - sus->Mh1[4][2] * sus->D1[2] - sus->Mh1[4][3] * sus->D1[3] - sus->Mh1[4][4] * sus->D1[4] - sus->Mh1[4][5] * sus->D1[5] + sus->Qh1[4] + sus->Qh1_TSDA[4]) +
		sus->B1[5] * (-sus->Mh1[5][0] * sus->D1[0] - sus->Mh1[5][1] * sus->D1[1] - sus->Mh1[5][2] * sus->D1[2] - sus->Mh1[5][3] * sus->D1[3] - sus->Mh1[5][4] * sus->D1[4] - sus->Mh1[5][5] * sus->D1[5] + sus->Qh1[5] + sus->Qh1_TSDA[5]);

	sus->inv_Mqq = 1 / Mqq;

	// Mhc = Myy - Myq*inv_Mqq*Myq'
	// Phc = -Myq*inv_Mqq*Pq
	for (i = 0; i < 6; i++) {
		temp[i] = -sus->Myq[i] * sus->inv_Mqq;
	}

	for (i = 0; i < 6; i++) {
		for (j = 0; j < 6; j++) {
			sus->Mhc[i][j] = Myy[i][j] + temp[i] * sus->Myq[j];
		}

		sus->Phc[i] = temp[i] * sus->Pq;
	}
}

int realtime_dynamics_analysis::overturn_detect(CHASSIS *chassis, SUS sus[6]) {
	double roll_over_ref = 65;
	double pitch_over_ref = 61;

	// -roll over
	if (sus[0].pen == 0 && sus[1].pen == 0 && sus[2].pen == 0) {
		if (chassis->roll_ang * 180 / M_PI <= -roll_over_ref) 
			return 1;
	}
	// +roll over
	if (sus[3].pen == 0 && sus[4].pen == 0 && sus[5].pen == 0) {
		if (chassis->roll_ang * 180 / M_PI >= roll_over_ref) 
			return 1;
	}
	// -pitch over
	if (sus[0].pen == 0 && sus[1].pen == 0 && sus[3].pen == 0 && sus[4].pen == 0) {
		if (chassis->pitch_ang * 180 / M_PI >= pitch_over_ref) 
			return 1;
	}
	// +pitch over
	if (sus[1].pen == 0 && sus[2].pen == 0 && sus[4].pen == 0 && sus[5].pen == 0) {
		if (chassis->pitch_ang * 180 / M_PI <= -pitch_over_ref) 
			return 1;
	}

	return 0;
}

void realtime_dynamics_analysis::acceleration_state_chassis(CHASSIS *chassis, SUS sus[6], SIM *sim) {
	double Py[6], Mh1_D1_RF[6], Mh1_D1_RM[6], Mh1_D1_RR[6], Mh1_D1_LF[6], Mh1_D1_LM[6], Mh1_D1_LR[6], M[6][6], P[6];
	int i, j;
	int indx[6];
	double M_fac[6][6];
	double Qh0_TSDA[6], Qh1_TSDA[6];

	for (i = 0; i < 6; i++) {
		Qh0_TSDA[i] = sus[RF].Qh0_TSDA[i] + sus[RM].Qh0_TSDA[i] + sus[RR].Qh0_TSDA[i] + sus[LF].Qh0_TSDA[i] + sus[LM].Qh0_TSDA[i] + sus[LR].Qh0_TSDA[i];
		Qh1_TSDA[i] = sus[RF].Qh1_TSDA[i] + sus[RM].Qh1_TSDA[i] + sus[RR].Qh1_TSDA[i] + sus[LF].Qh1_TSDA[i] + sus[LM].Qh1_TSDA[i] + sus[LR].Qh1_TSDA[i];
	}

	mat6661(sus[RF].Mh1, sus[RF].D1, Mh1_D1_RF);
	mat6661(sus[RM].Mh1, sus[RM].D1, Mh1_D1_RM);
	mat6661(sus[RR].Mh1, sus[RR].D1, Mh1_D1_RR);
	mat6661(sus[LF].Mh1, sus[LF].D1, Mh1_D1_LF);
	mat6661(sus[LM].Mh1, sus[LM].D1, Mh1_D1_LM);
	mat6661(sus[LR].Mh1, sus[LR].D1, Mh1_D1_LR);

	// Py = Qh0 + Qh1(1 ~ 6)
	for (i = 0; i < 6; i++) {
		Py[i] = chassis->Qh0[i] + (sus[RF].Qh1[i] - Mh1_D1_RF[i]) + (sus[RM].Qh1[i] - Mh1_D1_RM[i]) + (sus[RR].Qh1[i] - Mh1_D1_RR[i]) +
			(sus[LF].Qh1[i] - Mh1_D1_LF[i]) + (sus[LM].Qh1[i] - Mh1_D1_LM[i]) + (sus[LR].Qh1[i] - Mh1_D1_LR[i]) + Qh0_TSDA[i] + Qh1_TSDA[i];
	}

	// M = Mh0 + Mhc(1 ~ 6)
	// P = Py + Phc(1 ~ 6)
	for (i = 0; i < 6; i++) {
		for (j = 0; j < 6; j++) {
			M[i][j] = chassis->Mh0[i][j] + sus[RF].Mhc[i][j] + sus[RM].Mhc[i][j] + sus[RR].Mhc[i][j] + sus[LF].Mhc[i][j] + sus[LM].Mhc[i][j] + sus[LR].Mhc[i][j];
		}
		P[i] = Py[i] + sus[RF].Phc[i] + sus[RM].Phc[i] + sus[RR].Phc[i] + sus[LF].Phc[i] + sus[LM].Phc[i] + sus[LR].Phc[i];
	}

	// Ax=b solver(LU decomposition), 행렬 A의 singularity를 감지하기 위한 flag 변수 사용
	sim->singular_flag = ludcmp6(M, 6, indx, 0.0, M_fac);
	// 행렬 A가 singular matrix가 아닐 때 back substitution을 이용해서 Ax=b에서 x를 계산
	if (!sim->singular_flag) 
		lubksb6(M_fac, 6, indx, P, chassis->dYh0);
}

void realtime_dynamics_analysis::cartesian_acceleration_chassis(CHASSIS *chassis) {
	double w0t_w0t_rho0[3], w0t_rho0[3], dw0t_rho0[3], dYb0[6], dT0[6][6] = { 0, }, T0_dYh0[6] = { 0, }, dT0_Yh0[6] = { 0, };
	int i, j;

	// dT0 = [zeros(3) -dr0t; zeros(3,6)]
	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			dT0[i][j + 3] = -chassis->dr0t[i][j];
		}
	}

	// dYb0 = T0*dYh0 + dT0*Yh0
	mat6661(chassis->T0, chassis->dYh0, T0_dYh0);
	mat6661(dT0, chassis->Yh0, dT0_Yh0);
	for (i = 0; i < 6; i++) {
		dYb0[i] = T0_dYh0[i] + dT0_Yh0[i];
	}

	// ddr0 = dYb0(1:3)
	// dw0 = dYb0(4:6)
	for (i = 0; i < 3; i++) {
		chassis->ddr0[i] = dYb0[i];
		chassis->dw0[i] = dYb0[i + 3];
	}

	// ddr0c = ddrc + dw0t*rho0 + w0t*w0t*rho0
	tilde(chassis->dw0, chassis->dw0t);

	mat3331(chassis->dw0t, chassis->rho0, dw0t_rho0);
	mat3331(chassis->w0t, chassis->rho0, w0t_rho0);
	mat3331(chassis->w0t, w0t_rho0, w0t_w0t_rho0);

	for (i = 0; i < 3; i++) {
		chassis->ddr0c[i] = chassis->ddr0[i] + dw0t_rho0[i] + w0t_w0t_rho0[i];
	}
}

void realtime_dynamics_analysis::acceleration_state_suspension(CHASSIS *chassis, SUS *sus) {
	int i;

	// ddq1 = inv_Mqq*(Pq - Myq'*dYh0)
	sus->ddq1 = sus->inv_Mqq*((-sus->Myq[0] * chassis->dYh0[0] - sus->Myq[1] * chassis->dYh0[1] - sus->Myq[2] * chassis->dYh0[2] - sus->Myq[3] * chassis->dYh0[3] - sus->Myq[4] * chassis->dYh0[4] - sus->Myq[5] * chassis->dYh0[5]) + sus->Pq);

	// dYh1 = dYh0 + B1*ddq1 + D1
	for (i = 0; i < 6; i++) {
		sus->dYh1[i] = chassis->dYh0[i] + sus->B1[i] * sus->ddq1 + sus->D1[i];
	}
}

void realtime_dynamics_analysis::cartesian_acceleration_suspension(SUS *sus) {
	double w1t_w1t_rho1[3], w1t_rho1[3], dw1t_rho1[3], dYb1[6], dT1_Yh1[6], T1_dYh1[6], dT1[6][6] = { 0, };
	int i, j;

	// dT1 = [zeros(3) -dr1t; zeros(3,6)]
	for (i = 0; i < 3; i++) {
		for (j = 3; j < 6; j++) {
			dT1[i][j] = -sus->dr1t[i][j - 3];
		}
	}
	
	// dYb1 = T1*dYh1 + dT1*Yh1
	mat6661(sus->T1, sus->dYh1, T1_dYh1);
	mat6661(dT1, sus->Yh1, dT1_Yh1);
	for (i = 0; i < 6; i++) {
		dYb1[i] = T1_dYh1[i] + dT1_Yh1[i];
	}

	// ddr1 = dYb1(1:3)
	// dw1 = dYb1(4:6)
	for (i = 0; i < 3; i++) {
		sus->ddr1[i] = dYb1[i];
		sus->dw1[i] = dYb1[i + 3];
	}

	// ddr1c = ddr1 + dw1t*rho1 + dw1t*dw1t*rho1
	tilde(sus->dw1, sus->dw1t);

	mat3331(sus->dw1t, sus->rho1, dw1t_rho1);
	mat3331(sus->w1t, sus->rho1, w1t_rho1);
	mat3331(sus->w1t, w1t_rho1, w1t_w1t_rho1);

	for (i = 0; i < 3; i++) {
		sus->ddr1c[i] = sus->ddr1[i] + dw1t_rho1[i] + w1t_w1t_rho1[i];
	}
}

void realtime_dynamics_analysis::LP_control(CHASSIS *chassis, SUS sus[6], CTRL *ctrl, SIM *sim) {
/*
LP_control() : 자율 주행 제어기
작성자: 홍효성
Date: 2017.02.17

input variables
v_di: desired velocity
WP : Waypoint data
Yp : current derivative of state
step_size: integration step_size

output variables
motor_torque: 각 휠에 전달되는 모터 토크 명령

*/	

	////////// Section 0. 차량 파라미터 및 지역 변수 초기화 //////////

	////////// 입력 받아야할 차량 모델 파라미터 //////////
	double m_vehicle = 6762;		// vehicle mass [kg]
	double t_w = 1.948;             // 좌, 우측 휠 간 거리(m) = tread	(변경 2016.11.30 홍효성)	
	double I_z = 13201;             // total moment of inertia of z - axis(kgm ^ 2)
	double tire_radius = 0.5226;	// 타이어 반지름
	double torque_limit = 3416;		//[Nm], 기어비: 17.08 기준

	////////// 지역 변수 //////////
	int j;						// 반복문에 사용하기 위한 변수
	int indx[4] = { 0 };			// lusolve4 (LU decomponent) 함수에서 사용하기 위한 변수
	double M_fac[4][4] = { 0 };		// lusolve4 (LU decomponent) 함수에서 사용하기 위한 변수
	double dr0c_p[3];				// 차체의 로컬 속도
	int WP_pv_indx;						// preview point 계산을 위한 waypoint 인덱스
	double x_c, y_c;				// 글로벌 좌표에서 본 차량 무게 중심의 위치
	double yaw, yaw_rate;			//요 각도, 요 각속도
	double L = 1;                   // Preview distance 거리 배율(m) : 1로 고정	
	double F_z[6], F_x[6];			// 수직 방향 타이어력, 종 방향 타이어력
	double W1, W2, W3, W4, W5, W6;	// 타이어력 분배에 사용되는 weight factor	
	double F_x_A[4][4], F_x_B[4], F_x_1234[4], F_xd[6];	// 타이어력 분배에 사용되는 행렬
	double K_vp, K_vi, K_yp, K_yd;	// 제어 게인
	double M_c, M_c_true;			// 요 모멘트 제어 명령 및 실제 타이어력에서 발생되는 요 모멘트
	double T;						// step time
	double l, eta, p;				// 외란 요 모멘트 관측기에 사용되는 게인
	double F_tire_limit;				// 종 방향 타이어력 명령 limit, 
	double limit_ratio;				//모터 토크 명령이 한계 값을 초과할 경우 토크를 비례 배분하여 제한 위한 비율
	double WP_previous[4], WP_next[4];			// WP_previous: 차량이 위치하고 있는 현재 waypoint, WP_next: WP_previous 다음의 waypoint
	double WP_pv_previous[4], WP_pv_next[4];	// preview point 계산을 위한 현재 waypoint 및 다음 waypoint
	double d_wp;						// 각 waypoint 사이의 거리 (WP_previous 부터 WP_next 까지)
	double yaw_wp;					// 현재 소속된 waypoint의 방향(각도)
	double u_wp;					// 현재 소속된 waypoint 내에서 차량의 위치를 상대적으로 나타내기 위한 비율 (0~1까지 범위, u_wp > 1 인 경우 차량이 현재 waypoint를 지나쳤다는 의미)
	double x_n, y_n;				// 현재 소속된 waypoint 벡터 위에서 차량 무게중심에 대한 법선 벡터의 교차점
	double L_pv; // preview distance 
	double d_rest, L_pv_rest;	// preview distance에 의해 계산된 지향점까지의 거리가 해당 waypoint 범위보다 클 경우 초과된 만큼의 거리 (preview waypoint를 update 할 때 사용)
	double x_pv, y_pv, yaw_wp_pv; // preview distance에 의한 지향점 좌표(x_pv, y_pv) 및 차량 무게중심으로부터의 지향 각
	double de_psi;	// 지향각속도 오차 (desired yaw_rate - current yaw_rate)
	
	// 속도 및 지향각 제어 게인
	K_vp = 30;
	K_vi = 0;
	K_yd = 10;
	K_yp = 30;

	// 글로벌 속도 및 가속도를 rotation matrix(A0)를 이용하여 로컬 속도 및 가속도로 변환
	mat33T31(chassis->A0, chassis->ddr0c, ctrl->ddr0c_p);	// 가속도 ddr0c_p = A0 * ddr0c
	mat33T31(chassis->A0, chassis->dr0c, dr0c_p);			// 속도 dr0c_p = A0 * dr0c


	////////// 경로 추종 제어에 사용하기 위한 차량 상태 변수 입력 //////////

	yaw = chassis->yaw_ang;     // 차량의 지향각 (rad)
	yaw_rate = chassis->w0[2];  // 지향각속도 (rad/s)	
	x_c = chassis->r0c[0];		// 글로벌 좌표에서의 차량 무게중심의 x 좌표
	y_c = chassis->r0c[1];		// 글로벌 좌표에서의 차량 무게중심의 y 좌표
	ctrl->v_x = dr0c_p[0];		// 차량의 로컬 x 방향 속도(LOCAL)

	// 수직 방향 타이어력 (타이어 모델의 변수로부터 입력)
// 	F_z[0] = sus[LF].Fz;		// LF (Left, Front)
// 	F_z[1] = sus[RF].Fz;		// FR (Right, Front)
// 	F_z[2] = sus[LM].Fz;		// LM (Left, Middle)
// 	F_z[3] = sus[RM].Fz;		// RM (Right, Middle)
// 	F_z[4] = sus[LR].Fz;		// LR (Left, Right)
// 	F_z[5] = sus[RR].Fz;		// RR (Right, Rear)
	F_z[0] = m_vehicle / 6.0;   // LF (Left, Front)
	F_z[1] = m_vehicle / 6.0;   // FR (Right, Front)
	F_z[2] = m_vehicle / 6.0;   // LM (Left, Middle)
	F_z[3] = m_vehicle / 6.0;   // RM (Right, Middle)
	F_z[4] = m_vehicle / 6.0;   // LR (Left, Right)
	F_z[5] = m_vehicle / 6.0;   // RR (Right, Rear)

	// 종 방향 타이어력 (타이어 모델의 변수로부터 입력)
	F_x[0] = sus[LF].Fx;		// LF (Left, Front)
	F_x[1] = sus[RF].Fx;		// FR (Right, Front)
	F_x[2] = sus[LM].Fx;		// LM (Left, Middle)
	F_x[3] = sus[RM].Fx;		// RM (Right, Middle)
	F_x[4] = sus[LR].Fx;		// LR (Left, Right)
	F_x[5] = sus[RR].Fx;		// RR (Right, Rear)

	

	// 타이어력 분배에 사용되는 weight factor	정의
	W1 = 1; W2 = 1; W3 = 1; W4 = 1; W5 = 1; W6 = 1;

	// 외란 요 모멘트 관측기 게인
	T = sim->h;					// step time (동역학 모델과 동기화)
	l = 10 * I_z;				// observer gain
	eta = 25 * I_z;				// observer gain
	p = l / I_z;				// observer gain

	// 종 방향 타이어력 제한
	F_tire_limit = torque_limit / tire_radius;

	

	////////// Section 1. 지향 점 계산 ////////// 
	////////// Input: 차량 현재 위치 (x_c, y_c), Waypoint data (WP)
	////////// Output: 지향각 명령 (yaw_d)

	if (sim->simulation_flag == 3) {										// 병렬 주행 (RTT) 모드인 경우에만 계산 (Equilibrium 모드에서는 계산 안함)
		if (sim->t_current >= 0.0 && sim->t_current <= 0.0) {
			// waypoint 시작 점에 차량 현재 주행 속도 입력 (optimal_velocity_profile에 사용하기 위해 초기화 해줌)
			DID->VehicleVelocity[0][sim->thread_indx] = ctrl->v_x;
		}
		WP_previous[0] = DID->LocalPath[ctrl->WP_indx - 1][0];
		WP_previous[1] = DID->LocalPath[ctrl->WP_indx - 1][1];
		WP_previous[2] = DID->VehicleVelocity[ctrl->WP_indx - 1][sim->thread_indx];
		WP_previous[3] = DID->StabilityIndex[ctrl->WP_indx - 1][sim->thread_indx];
		WP_next[0] = DID->LocalPath[ctrl->WP_indx][0];
		WP_next[1] = DID->LocalPath[ctrl->WP_indx][1];
		WP_next[2] = DID->VehicleVelocity[ctrl->WP_indx][sim->thread_indx];
		WP_next[3] = DID->StabilityIndex[ctrl->WP_indx][sim->thread_indx];

		d_wp = sqrt((WP_next[0] - WP_previous[0])*(WP_next[0] - WP_previous[0]) + (WP_next[1] - WP_previous[1])*(WP_next[1] - WP_previous[1]));		// WP_previous과 WP_next 사이의 거리 계산
		yaw_wp = atan2((WP_next[1] - WP_previous[1]), (WP_next[0] - WP_previous[0]));										// waypoint의 direction(rad) 계산
		u_wp = ((x_c - WP_previous[0])*(WP_next[0] - WP_previous[0]) + (y_c - WP_previous[1])*(WP_next[1] - WP_previous[1])) / (d_wp*d_wp);	// WP_previous과 WP_next 사이에서 차량의 상대 위치를 비율로 계산 (0 ~ 1 범위)

		while (u_wp > 1)													// u_wp > 1 : 차량의 위치가 WP_next 범위를 초과한 경우 -> 새로운 waypoint 할당
		{
			//WP[ctrl->WP_indx][2][sim->thread_indx] = ctrl->v_x;			// 해당 WP_indx 앞에서의 차량 현재 속도 입력
			if (ctrl->v_x > DID->vd_mission) {
				DID->VehicleVelocity[ctrl->WP_indx][sim->thread_indx] = DID->vd_mission;
			}
			else {
				DID->VehicleVelocity[ctrl->WP_indx][sim->thread_indx] = ctrl->v_x;
			}

			if (ctrl->WP_indx >= ctrl->WP_size - 1) {						// WP_indx가 마지막 값에 도달하면 더 이상 WP_indx를 업데이트 하지 않음
				break;
			}

			ctrl->WP_indx++;												// WP_previous과 WP_next 좌표를 새로 업데이트하기 위해 WP_indx를 +1 시킨다.
			ctrl->WP_indx_update_flag = 1;									// WP_indx가 업데이트 된 경우 LP_stability_metric()가 실행됐을 때 stability indx 항목을 1로 초기화 시켜준다.

			// waypoint 새로 할당
// 			for (i = 0; i < 4; i++) {
// 				WP_previous[i] = WP[ctrl->WP_indx - 1][i][sim->thread_indx];		// 현재 차량이 소속된 좌표인 WP_previous에 waypoint의 X, Y 좌표와 주행 속도 명령 그리고 stability indx = 0 할당
// 				WP_next[i] = WP[ctrl->WP_indx][i][sim->thread_indx];			// 현재 차량이 소속된 좌표인 WP_next에 waypoint의 X, Y 좌표와 주행 속도 명령 그리고 stability indx = 0 할당
// 			}
			WP_previous[0] = DID->LocalPath[ctrl->WP_indx - 1][0];
			WP_previous[1] = DID->LocalPath[ctrl->WP_indx - 1][1];
			WP_previous[2] = DID->VehicleVelocity[ctrl->WP_indx - 1][sim->thread_indx];
			WP_previous[3] = DID->StabilityIndex[ctrl->WP_indx - 1][sim->thread_indx];
			WP_next[0] = DID->LocalPath[ctrl->WP_indx][0];
			WP_next[1] = DID->LocalPath[ctrl->WP_indx][1];
			WP_next[2] = DID->VehicleVelocity[ctrl->WP_indx][sim->thread_indx];
			WP_next[3] = DID->StabilityIndex[ctrl->WP_indx][sim->thread_indx];

			d_wp = sqrt((WP_next[0] - WP_previous[0])*(WP_next[0] - WP_previous[0]) + (WP_next[1] - WP_previous[1])*(WP_next[1] - WP_previous[1]));		// WP_previous과 WP_next 사이의 거리 계산
			yaw_wp = atan2((WP_next[1] - WP_previous[1]), (WP_next[0] - WP_previous[0]));										// waypoint의 direction(rad) 계산
			u_wp = ((x_c - WP_previous[0])*(WP_next[0] - WP_previous[0]) + (y_c - WP_previous[1])*(WP_next[1] - WP_previous[1])) / (d_wp*d_wp);	// WP_previous과 WP_next 사이에서 차량의 상대 위치를 비율로 계산 (0 ~ 1 범위)
		}

		x_n = WP_previous[0] + u_wp*(WP_next[0] - WP_previous[0]);		// 현재 소속된 waypoint 벡터 위에서 차량 무게중심에 대한 법선 벡터의 교차점 X 좌표
		y_n = WP_previous[1] + u_wp*(WP_next[1] - WP_previous[1]);		// 현재 소속된 waypoint 벡터 위에서 차량 무게중심에 대한 법선 벡터의 교차점 Y 좌표

		// 횡 방향 위치 오차 계산 (lateral position error)
		if (y_n > y_c) {							// 원래 LPE는 항상 양수로 계산되지만 경로에 대하여 차량 Y 좌표의 상대 위치에 따라 부호가 바뀌어 보이도록 임시 조치함
			ctrl->e_l = sqrt((x_n - x_c)*(x_n - x_c) + (y_n - y_c)*(y_n - y_c));		// e_l : 횡 방향 위치 오차
		}
		else if (y_n < y_c) {
			ctrl->e_l = sqrt((x_n - x_c)*(x_n - x_c) + (y_n - y_c)*(y_n - y_c)) * (-1);
		}
		else {
			ctrl->e_l = sqrt((x_n - x_c)*(x_n - x_c) + (y_n - y_c)*(y_n - y_c)) * 0;
		}

		// Preview distance
		if (fabs(ctrl->v_x) > L) {		// Preview distance는 차량 주행 속도에 비례하도록 설정함
			L_pv = fabs(ctrl->v_x)*L;
		}
		else {
			L_pv = L;					// 차량 주행 속도가 1 m/s 보다 작을 때는 L_pv = 1 (m) 로 지정함. 여기서 L = 1 이다.
		}

		if (ctrl->WP_indx < ctrl->WP_size) {	// WP_indx가 WP_indx의 마지막 지점에 도달하지 않았을 경우
			WP_pv_indx = ctrl->WP_indx;			// Preview waypoint의 indx를 현재 waypoint의 indx로 동기화
			//d_rest = sqrt((WP[WP_pv_indx][0][sim->thread_indx] - x_n)*(WP[WP_pv_indx][0][sim->thread_indx] - x_n) + (WP[WP_pv_indx][1][sim->thread_indx] - y_n)*(WP[WP_pv_indx][1][sim->thread_indx] - y_n));	// 현재 waypoint의 말단인 WP_next까지 남은 거리 (distance_rest)
			d_rest = sqrt((DID->LocalPath[WP_pv_indx][0] - x_n) * (DID->LocalPath[WP_pv_indx][0] - x_n) + (DID->LocalPath[WP_pv_indx][1] - y_n)*(DID->LocalPath[WP_pv_indx][1] - y_n));

			if (L_pv > d_rest) {			// Preview distance가 현재 waypoint의 범위를 초과할 경우
				L_pv_rest = L_pv - d_rest;	// 초과된 만큼을 계산하여 L_pv_rest에 저장
				while (L_pv_rest > 0) {		// 초과된 거리가 양수이면					
					WP_pv_indx = WP_pv_indx + 1;		// Preview waypoint의 indx을 +1 업데이트 시킴 (Preview waypoint 업데이트를 위함)
					if (WP_pv_indx >= ctrl->WP_size) {	// 만약 Preview waypoint의 indx가 전체 waypoint indx 범위를 초과할 경우 마지막 waypoint를 preview point로 할당
// 						for (i = 0; i < 4; i++) {
// 							WP_pv_previous[i] = WP[ctrl->WP_size - 2][i][sim->thread_indx];		// 마지막 waypoint를 preview waypoint로 할당 (X 좌표)
// 							WP_pv_next[i] = WP[ctrl->WP_size - 1][i][sim->thread_indx];		// 마지막 waypoint를 preview waypoint로 할당 (Y 좌표)
// 						}
						WP_previous[0] = DID->LocalPath[ctrl->WP_size - 2][0];
						WP_previous[1] = DID->LocalPath[ctrl->WP_size - 2][1];
						WP_previous[2] = DID->VehicleVelocity[ctrl->WP_size - 2][sim->thread_indx];
						WP_previous[3] = DID->StabilityIndex[ctrl->WP_size - 2][sim->thread_indx];
						WP_next[0] = DID->LocalPath[ctrl->WP_size - 1][0];
						WP_next[1] = DID->LocalPath[ctrl->WP_size - 1][1];
						WP_next[2] = DID->VehicleVelocity[ctrl->WP_size - 1][sim->thread_indx];
						WP_next[3] = DID->StabilityIndex[ctrl->WP_size - 1][sim->thread_indx];
						break;		
					}
					
// 					for (i = 0; i < 4; i++) {
// 						WP_pv_previous[i] = WP[WP_pv_indx - 1][i][sim->thread_indx];				// 다음 waypoint를 preview waypoint로 할당 (X 좌표)
// 						WP_pv_next[i] = WP[WP_pv_indx][i][sim->thread_indx];					// 다음 waypoint를 preview waypoint로 할당 (X 좌표)
// 					}
					WP_pv_previous[0] = DID->LocalPath[WP_pv_indx - 1][0];
					WP_pv_previous[1] = DID->LocalPath[WP_pv_indx - 1][1];
					WP_pv_previous[2] = DID->VehicleVelocity[WP_pv_indx - 1][sim->thread_indx];
					WP_pv_previous[3] = DID->StabilityIndex[WP_pv_indx - 1][sim->thread_indx];
					WP_pv_next[0] = DID->LocalPath[WP_pv_indx][0];
					WP_pv_next[1] = DID->LocalPath[WP_pv_indx][1];
					WP_pv_next[2] = DID->VehicleVelocity[WP_pv_indx][sim->thread_indx];
					WP_pv_next[3] = DID->StabilityIndex[WP_pv_indx][sim->thread_indx];

					d_rest = sqrt((WP_pv_next[0] - WP_pv_previous[0])*(WP_pv_next[0] - WP_pv_previous[0]) + (WP_pv_next[1] - WP_pv_previous[1])*(WP_pv_next[1] - WP_pv_previous[1]));	// 다음 waypoint의 WP_next까지 남은 거리 다시 계산
					L_pv_rest = L_pv_rest - d_rest;		// 초과된 만큼을 계산하여 L_pv_rest에 저장 (초과될 경우 양수, 초과되지 않을 경우 음수를 가진다.)
				}
				L_pv_rest = d_rest + L_pv_rest;			// L_pv_rest가 더이상 초과되지 않는 음수를 가진 경우 d_rest와 더해주어 해당 preview waypoint에서 남은 거리를 다시 계산한다.
				yaw_wp_pv = atan2((WP_pv_next[1] - WP_pv_previous[1]), (WP_pv_next[0] - WP_pv_previous[0]));	// Preview waypoint의 방향 (각도)
				x_pv = WP_pv_previous[0] + L_pv_rest*cos(yaw_wp_pv);	// Preview waypoint의 방향과 남은 거리를 이용하여 preview waypoint 벡터 내에서의 지향점(x_pv)를 계산한다.
				y_pv = WP_pv_previous[1] + L_pv_rest*sin(yaw_wp_pv);	// Preview waypoint의 방향과 남은 거리를 이용하여 preview waypoint 벡터 내에서의 지향점(y_pv)를 계산한다.
			}
			else {		// Preview distance가 현재 waypoint의 범위를 초과되지 않은 경우
				x_pv = x_n + L_pv*cos(yaw_wp);	// 현재 waypoint의 방향과 남은 거리를 이용하여 waypoint 벡터 내에서의 지향점(x_pv)를 계산한다.
				y_pv = y_n + L_pv*sin(yaw_wp);	// 현재 waypoint의 방향과 남은 거리를 이용하여 waypoint 벡터 내에서의 지향점(y_pv)를 계산한다.
			}
		}
		else {			// WP_indx가 WP_indx의 마지막 지점에 도달한 경우
			x_pv = x_n + L_pv*cos(yaw_wp);	// 현재 waypoint, (=마지막 waypoint)에서 waypoint의 방향과 preview distance를 이용하여 지향점 할당 (x_pv)
			y_pv = y_n + L_pv*sin(yaw_wp);	// 현재 waypoint, (=마지막 waypoint)에서 waypoint의 방향과 preview distance를 이용하여 지향점 할당 (y_pv)
		}

		ctrl->yaw_d = atan2((y_pv - y_c), (x_pv - x_c));    // 차량 무게중심으로부터 지향점(preview point)까지의 벡터의 방향을 지향각 명령(desired yaw angle)로 계산한다.
	}

	if (sim->simulation_flag == 1) {		// Equilibrium 상태에서는 최초에 한번 지향 각 명령을 현재 지향각으로 설정함
		if (sim->t_current >= 0.0 && sim->t_current <= 0.0) {
			ctrl->yaw_d = yaw;
		}
		ctrl->e_l = 0;						// save_data()에 사용하기 위한 입력
	}




	////////// Section 2. 속도 및 지향각 명령 계산 ////////// 
	////////// Input: 차량 주행 속도 (v_x), 지향각 (yaw), 지향각속도 (yaw_rate), 종 방향 타이어력(F_x)
	////////// Output: 종 방향 전체 힘 (F_xd_total), 요 모멘트 제어 명령 (M_c)
	
	// 외란 요 모멘트 추정 관측기
	M_c_true = t_w / 2 * (F_x[1] + F_x[3] + F_x[5] - F_x[0] - F_x[2] - F_x[4]);     // 실제 타이어력으로부터 발생하는 모멘트 계산
	ctrl->yaw_hat = ctrl->yaw_hat + T / I_z*ctrl->M_d_hat + T / I_z*M_c_true + T*p*(yaw_rate - ctrl->yaw_hat);	// yaw rate estimate
	ctrl->M_d_hat = ctrl->M_d_hat + T*eta*(yaw_rate - ctrl->yaw_hat);	// disturbance moment estimate	
	
	// 지향각 오차 계산
	if (yaw >= 0) {										// 지향각이 1,2 사분면에 위치한 경우 (양수)
		if (ctrl->yaw_d < (yaw - M_PI)) {					// 지향각 명령이 (지향각-pi)보다 작은 경우
			ctrl->e_psi = ctrl->yaw_d + 2 * M_PI - yaw;
		}
		else {												// 지향각 명령이 (지향각-pi)보다 큰 경우
			ctrl->e_psi = ctrl->yaw_d - yaw;
		}
	}
	else {												// 지향각이 3,4 사분면에 위치한 경우 (음수)
		if (ctrl->yaw_d > (yaw + M_PI)) {					// 지향각 명령이 (지향각+pi)보다 큰 경우
			ctrl->e_psi = (ctrl->yaw_d - 2 * M_PI) - yaw;
		}
		else {												// 지향각 명령이 (지향각+pi)보다 작은 경우
			ctrl->e_psi = ctrl->yaw_d - yaw;
		}
	}
	
	// 지향 각속도 오차 계산
	de_psi = -yaw_rate;	// de_psi = desired_yaw_rate - yaw_rate (여기서 desired_yaw_rate는 0으로 놓고 계산함)

	// desired velocity
	if (sim->simulation_flag == 1) {	// Equilibrium 상태인 경우
		ctrl->v_di = 0;					// Equilibrium 상태일 때는 차량 속도를 0으로 제어
	}
	else {
		ctrl->v_di = ctrl->v_d[sim->thread_indx];		// 각 CPU 쓰레드에 할당된 속도 명령을 입력
	}

	// 주행 속도 오차 계산
	ctrl->e_v = ctrl->v_di - ctrl->v_x;				// 속도 오차
	ctrl->e_v_sum = ctrl->e_v_sum + ctrl->e_v*T;	// 속도 오차 적분


	// 종방향 힘 및 제어 모멘트 명령
	ctrl->F_xd_total = m_vehicle*(K_vp*ctrl->e_v + K_vi*ctrl->e_v_sum);		// 종 방향 전체 힘 명령
	M_c = I_z*(K_yd*de_psi + K_yp*ctrl->e_psi) - ctrl->M_d_hat;					// 요 모멘트 제어 명령
	


	////////// Section 3. 인-휠 모터 토크 명령 계산 //////////
	////////// Output: 각 휠의 모터 토크 (motoreque[6])
	

	// 종 방향 타이어력 각 휠에 분배
	F_x_A[0][0] = 2 * (W1 / ((F_z[0] + 1)*(F_z[0] + 1)) + W5 / ((F_z[4] + 1)*(F_z[4] + 1)));
	F_x_A[0][1] = 0;
	F_x_A[0][2] = 2 * W5 / ((F_z[4] + 1)*(F_z[4] + 1));
	F_x_A[0][3] = 0;
	F_x_A[1][0] = 0;
	F_x_A[1][1] = 2 * (W2 / ((F_z[1] + 1)*(F_z[1] + 1)) + W6 / ((F_z[5] + 1)*(F_z[5] + 1)));
	F_x_A[1][2] = 0;
	F_x_A[1][3] = 2 * W6 / ((F_z[5] + 1)*(F_z[5] + 1));
	F_x_A[2][0] = 2 * W5 / ((F_z[4] + 1)*(F_z[4] + 1));
	F_x_A[2][1] = 0;
	F_x_A[2][2] = 2 * (W3 / ((F_z[2] + 1)*(F_z[2] + 1)) + W5 / ((F_z[4] + 1)*(F_z[4] + 1)));
	F_x_A[2][3] = 0;
	F_x_A[3][0] = 0;
	F_x_A[3][1] = 2 * W6 / ((F_z[5] + 1)*(F_z[5] + 1));
	F_x_A[3][2] = 0;
	F_x_A[3][3] = 2 * (W4 / ((F_z[3] + 1)*(F_z[3] + 1)) + W6 / ((F_z[5] + 1)*(F_z[5] + 1)));

	F_x_B[0] = W5 / ((F_z[4] + 1)*(F_z[4] + 1))*ctrl->F_xd_total - 2 * W5 / ((F_z[4] + 1)*(F_z[4] + 1))*M_c / t_w;
	F_x_B[1] = W6 / ((F_z[5] + 1)*(F_z[5] + 1))*ctrl->F_xd_total + 2 * W6 / ((F_z[5] + 1)*(F_z[5] + 1))*M_c / t_w;
	F_x_B[2] = W5 / ((F_z[4] + 1)*(F_z[4] + 1))*ctrl->F_xd_total - 2 * W5 / ((F_z[4] + 1)*(F_z[4] + 1))*M_c / t_w;
	F_x_B[3] = W6 / ((F_z[5] + 1)*(F_z[5] + 1))*ctrl->F_xd_total + 2 * W6 / ((F_z[5] + 1)*(F_z[5] + 1))*M_c / t_w;

	// LU decomposition 방법으로 F_x_1234 = inv(F_x_A)*F_x_B 계산 (2016.11.29 홍효성)
	sim->singular_flag2 = ludcmp4(F_x_A, 4, indx, 0.0, M_fac);
	if (sim->singular_flag2 == 0) {
		lusolve4(M_fac, 4, indx, F_x_B, F_x_1234);

		// 각 휠에서의 종 방향 타이어력 제어 명령
		F_xd[0] = F_x_1234[0];
		F_xd[1] = F_x_1234[1];
		F_xd[2] = F_x_1234[2];
		F_xd[3] = F_x_1234[3];
		F_xd[4] = ctrl->F_xd_total / 2 - M_c / t_w - F_xd[0] - F_xd[2];
		F_xd[5] = ctrl->F_xd_total / 2 + M_c / t_w - F_xd[1] - F_xd[3];

		// 종 방향 타이어력 명령이 한계를 초과할 경우 한계치 만큼만 제어 명령을 가지도록 각 휠에 비례 분배시킨다.
		if (fabs(F_xd[0]) > F_tire_limit || fabs(F_xd[1]) > F_tire_limit || fabs(F_xd[2]) > F_tire_limit || fabs(F_xd[3]) > F_tire_limit || fabs(F_xd[4]) > F_tire_limit || fabs(F_xd[5]) > F_tire_limit)
		{
			limit_ratio = F_tire_limit / Find_Max(F_xd, 6);	// 최대 타이어력 제어 명령을 타이어력 한계치 수준에 맞추기 위한 비율 계산
			for (j = 0; j < 6; j++) {
				F_xd[j] = F_xd[j] * limit_ratio;
			}
		}

		// 인-휠 모터 토크 명령 계산
		for (j = 0; j < 6; j++) {
			ctrl->motor_torque[j] = tire_radius*F_xd[j];		// 모터 토크 명령 = 타이어 반지름 * 종 방향 힘 명령
			if (F_z[j] >= 0 && F_z[j] <= 0)
				ctrl->motor_torque[j] = 0;						// 수직 타이어력이 0 인 경우 모터 토크도 0으로 설정
		}


		////////// Section 4. 주행안정성지표 계산 함수 실행 //////////
		if (sim->simulation_flag == 3) {	// 병렬 주행(RTT) 모드에서만 작동
			LP_stability_metric(chassis, sus, sim, ctrl, sim->thread_indx);
		}
	}
}

void realtime_dynamics_analysis::LP_stability_metric(CHASSIS *chassis, SUS sus[6], SIM *sim, CTRL *ctrl, int thread_indx) {		// 주행 안정성 지표 계산
	/*
	LP_stability_metric() : 주행 안정성 지표 계산
	작성자: 홍효성
	Date: 2017.02.20

	input variables
	차량 동역학 상태 변수

	output variables
	4-DOF 주행 안정성 지표 (Roll, Pitch, Lateral, Vertical acceleration)
	*/

	////////// Section 0. 차량 파라미터 및 지역 변수 초기화 //////////

	////////// 입력 받아야할 차량 모델 파라미터 //////////	
	#define MAF_TAP 10
	double VF_side_max = 32450;		// 차량이 평지 위에 가만히 놓여 있을 때 좌측 또는 우측 세 개의 수직 타이어력 합
	double VF_FM_max = 40888;		// 차량이 평지 위에 가만히 놓여 있을 때 전방과 중간 휠 네 개의 수직 타이어력 합
	double VF_RM_max = 45538;		// 차량이 평지 위에 가만히 놓여 있을 때 후방과 중간 휠 네 개의 수직 타이어력 합
	double lateral_position_max = 0.3 + 0.3*0;			// 최대 허용 가능한 횡 방향 위치 오차
	double a_z_max = 1.5;				// 최대 허용 수직 가속도. 단위: G
	double RSM_max = 0.65;				// 롤 안정성 한계값
	double PSM_max = 0.65;				// 피치 안정성 한계값
	double LSM_max = 0.9;//0.25;				// 횡 방향 안정성 한계값
	double VSM_max = 0.5;//0.53;				// 수직 가속도 안정성 한계값 0.53    0.5: 35.39 kph ,  0.56: 31.26 kph
	double rate_epsilon = 0.01;			// RSM 및 PSM이 작동하기 위한 최소 angular rate margin

	// RTT 속도를 줄이기 위한 주행 안정성 지표 별 성능 저하 파라미터 (1: 100% 성능, 0: 0% 성능 = 무조건 fail이므로 최저 속도)
	double RSM_degrade = 1.0;
	double PSM_degrade = 1.0;
	double LSM_degrade = 1.0;
	double VSM_degrade = 1.0;

	if (RSM_degrade > 1.0) RSM_degrade = 1.0;
	if (RSM_degrade < 0.0) RSM_degrade = 0.0;
	if (PSM_degrade > 1.0) PSM_degrade = 1.0;
	if (PSM_degrade < 0.0) PSM_degrade = 0.0;
	if (LSM_degrade > 1.0) LSM_degrade = 1.0;
	if (LSM_degrade < 0.0) LSM_degrade = 0.0;
	if (VSM_degrade > 1.0) VSM_degrade = 1.0;
	if (VSM_degrade < 0.0) VSM_degrade = 0.0;

	RSM_max = RSM_max + (1.0 - RSM_max)*(1.0 - RSM_degrade);
	PSM_max = PSM_max + (1.0 - PSM_max)*(1.0 - PSM_degrade);
	LSM_max = LSM_max + (1.0 - LSM_max)*(1.0 - LSM_degrade);
	VSM_max = VSM_max + (1.0 - VSM_max)*(1.0 - VSM_degrade);

	////////// 지역 변수 //////////
	int j, k;			// 반복문에 사용하기 위한 변수
	double roll_angle, roll_rate, pitch_angle, pitch_rate, yaw_angle, yaw_rate;		// 롤, 피치, 요 각 및 각속도
	double a_z;	// 로컬 수직 가속도
	double F_z[6], VF_left, VF_right;	// 수직 방향 타이어력, 좌우 측 휠의 수직 방향 타이어력
	double VF_RSM_mean[2];		// 롤 안정성 지표 계산에 사용할 수직 방향 타이어력의 평균 값 (0: 좌측, 1: 우측)
	double VF_PSM_mean[2];		// 피치 안정성 지표 계산에 사용할 수직 방향 타이어력의 평균 값 (0: 전방+중간, 1: 후방+중간)
	double VF_FM, VF_RM;		// 전방+중간 수직 타이어력, 후방+중간 수직 타이어력
	double RSM_threshold = 0.38;	// 약 21.7도
	double PSM_threshold = 0.78;	// 약 45도
	double Roll_rate_lambda = 0.5;
	double Pitch_rate_lambda = 0.1;
	double RSM_margin = 0.0;
	double PSM_margin = 0.0;
	static double VF_RSM_buffer[MAF_TAP][2], VF_PSM_buffer[MAF_TAP][2];	// 롤 및 피치 안정성 지표에서 수직 타이어력의 이동 평균 필터에 사용하기 위한 버퍼

	////////// 주행 안정성 지표에 사용하기 위한 차량 상태 변수 입력 //////////
	// 수직 방향 타이어력
	F_z[0] = sus[LF].Fz;		// Left Front
	F_z[1] = sus[RF].Fz;		// Right Front
	F_z[2] = sus[LM].Fz;		// Left Middle
	F_z[3] = sus[RM].Fz;		// Right Middle
	F_z[4] = sus[LR].Fz;		// Left Rear
	F_z[5] = sus[RR].Fz;		// Right Rear

	// 롤, 피치, 요 각도 및 각속도
	roll_angle = chassis->roll_ang;     // 롤 각
	roll_rate = chassis->w0p[0];			// 롤 각속도
	pitch_angle = chassis->pitch_ang;   // 피치 각
	pitch_rate = chassis->w0p[1];		// 피치 각속도
	yaw_angle = chassis->yaw_ang;       // 요 각
	yaw_rate = chassis->w0p[2];			// 요 각속도

										////////// 롤 안정성 지표 //////////

										// Roll Stability Metric (ver 5): 각도 + 각속도와 수직 타이어력 사용
	RSM_margin = fabs(roll_angle + Roll_rate_lambda*roll_rate) / RSM_threshold;		// 롤 각도와 롤 각속도를 결합한 에너지 컨셉 마진
	
	VF_left = F_z[0] + F_z[2] + F_z[4];		// 좌측 수직 타이어력의 합
	VF_right = F_z[1] + F_z[3] + F_z[5];	// 우측 수직 타이어력의 합

	for (k = 1; k < MAF_TAP; k++) {
		for (j = 0; j < 2; j++) {	// j = 0 : left, j = 1 : right
			VF_RSM_buffer[k - 1][j] = VF_RSM_buffer[k][j];		// Moving average filter 계산을 위해 버퍼 내부 데이터 이동
		}
	}
	VF_RSM_buffer[MAF_TAP - 1][0] = VF_left;	// 버퍼의 마지막 배열에 수직 타이어력 합 추가 (좌측 타이어)
	VF_RSM_buffer[MAF_TAP - 1][1] = VF_right;	// 버퍼의 마지막 배열에 수직 타이어력 합 추가 (우측 타이어)

												// 타이어력 평균 값 계산 (Moving average filter)
	VF_RSM_mean[0] = 0;
	VF_RSM_mean[1] = 0;

	for (k = 0; k < MAF_TAP; k++) {
		for (j = 0; j < 2; j++) {
			VF_RSM_mean[j] += VF_RSM_buffer[k][j];	// 타이어력 평균 계산을 위해 각 버퍼에 저장된 데이터를 모두 더해준다.
		}
	}

	VF_RSM_mean[0] = VF_RSM_mean[0] / MAF_TAP;		// 타이어력의 평균값 = moving average filter, left
	VF_RSM_mean[1] = VF_RSM_mean[1] / MAF_TAP;		// 타이어력의 평균값 = moving average filter, right
	if (RSM_margin > 0.5) {			// 롤 각도와 각속도가 일정 값을 초과할 경우 롤 모션 불안정 상태
		if (roll_rate > 0) {
			ctrl->RSM = VF_RSM_mean[0] / VF_side_max;	// 좌측 타이어력의 비율을 RSM으로 계산
		}
		else {
			ctrl->RSM = VF_RSM_mean[1] / VF_side_max;	// 좌측 타이어력의 비율을 RSM으로 계산
		}

		if (ctrl->RSM > 1.0) {
			ctrl->RSM = 1.0;
		}
	}
	else {
		ctrl->RSM = 1.0 - RSM_margin;
	}


	////////// 피치 안정성 지표 //////////

	// Pitch Stability Metric (ver 5): 각도 + 각속도와 수직 타이어력 사용
	PSM_margin = fabs(pitch_angle + Pitch_rate_lambda*pitch_rate) / PSM_threshold;

	VF_FM = F_z[0] + F_z[1] + F_z[2] + F_z[3];	// Front + Middle 수직 타이어력의 합
	VF_RM = F_z[2] + F_z[3] + F_z[4] + F_z[5];	//Rear + Middle 수직 타이어력의 합

	for (k = 1; k < MAF_TAP; k++) {
		for (j = 0; j < 2; j++) {	// j = 0 : Front + Middle, j = 1 : Rear + Middle
			VF_PSM_buffer[k - 1][j] = VF_PSM_buffer[k][j];	// Moving average filter 계산을 위해 버퍼 내부 데이터 이동
		}
	}
	VF_PSM_buffer[MAF_TAP - 1][0] = VF_FM;	// 버퍼의 마지막 배열에 수직 타이어력 합 추가 (전방+중간 타이어)
	VF_PSM_buffer[MAF_TAP - 1][1] = VF_RM;	// 버퍼의 마지막 배열에 수직 타이어력 합 추가 (후방+중간 타이어)

											// 타이어력 평균 값 계산 (Moving average filter)
	VF_PSM_mean[0] = 0;
	VF_PSM_mean[1] = 0;

	for (k = 0; k < MAF_TAP; k++) {
		for (j = 0; j < 2; j++) {
			VF_PSM_mean[j] += VF_PSM_buffer[k][j];	// 타이어력 평균 계산을 위해 각 버퍼에 저장된 데이터를 모두 더해준다.
		}
	}

	VF_PSM_mean[0] = VF_PSM_mean[0] / MAF_TAP;		// 타이어력의 평균값 = moving average filter, FM(uphill)
	VF_PSM_mean[1] = VF_PSM_mean[1] / MAF_TAP;		// 타이어력의 평균값 = moving average filter, RM(downhill)
		
	if (PSM_margin > 0.5) {
		if (pitch_rate > 0) {							// FM(uphill)
			ctrl->PSM = VF_PSM_mean[0] / VF_FM_max;		// 전방+중간 수직 타이어력의 비율을 PSM으로 계산
		}
		else {											// RM(downhill)
			ctrl->PSM = VF_PSM_mean[1] / VF_RM_max;		// 후방+중간 수직 타이어력의 비율을 PSM으로 계산
		}

		if (ctrl->PSM > 1.0) {
			ctrl->PSM = 1.0;
		}
	}
	else {
		ctrl->PSM = 1.0 - PSM_margin;
	}
	//////////// 롤 안정성 지표 //////////

	//// Roll Stability Metric (ver 4): 각속도와 수직 타이어력 사용

	//VF_left = F_z[0] + F_z[2] + F_z[4];		// 좌측 수직 타이어력의 합
	//VF_right = F_z[1] + F_z[3] + F_z[5];	// 우측 수직 타이어력의 합

	//for (k = 1; k < MAF_TAP; k++) {
	//	for (j = 0; j < 2; j++) {	// j = 0 : left, j = 1 : right
	//		VF_RSM_buffer[k - 1][j] = VF_RSM_buffer[k][j];		// Moving average filter 계산을 위해 버퍼 내부 데이터 이동
	//	}
	//}
	//VF_RSM_buffer[MAF_TAP - 1][0] = VF_left;	// 버퍼의 마지막 배열에 수직 타이어력 합 추가 (좌측 타이어)
	//VF_RSM_buffer[MAF_TAP - 1][1] = VF_right;	// 버퍼의 마지막 배열에 수직 타이어력 합 추가 (우측 타이어)

	//// 타이어력 평균 값 계산 (Moving average filter)
	//VF_RSM_mean[0] = 0;
	//VF_RSM_mean[1] = 0;

	//for (k = 0; k < MAF_TAP; k++) {
	//	for (j = 0; j < 2; j++) {
	//		VF_RSM_mean[j] += VF_RSM_buffer[k][j];	// 타이어력 평균 계산을 위해 각 버퍼에 저장된 데이터를 모두 더해준다.
	//	}
	//}

	//VF_RSM_mean[0] = VF_RSM_mean[0] / MAF_TAP;		// 타이어력의 평균값 = moving average filter, left
	//VF_RSM_mean[1] = VF_RSM_mean[1] / MAF_TAP;		// 타이어력의 평균값 = moving average filter, right

	//// RSM 계산
	//if (roll_rate > rate_epsilon) {							// 차량이 오른쪽으로 기울어지면
	//	ctrl->RSM = VF_RSM_mean[0] / VF_side_max;	// 좌측 타이어력의 비율을 RSM으로 계산
	//	if (VF_RSM_mean[0] > VF_side_max)			// 만약 타이어력의 크기가 평상시보다 클 경우 (분자가 분모보다 클 경우)
	//		ctrl->RSM = 1;							// RSM = 1로 고정
	//}
	//else if (roll_rate < -rate_epsilon) {						// 차량이 왼쪽으로 기울어지면
	//	ctrl->RSM = VF_RSM_mean[1] / VF_side_max;	// 우측 타이어력의 비율을 RSM으로 계산
	//	if (VF_RSM_mean[1] > VF_side_max)			// 만약 타이어력의 크기가 평상시보다 클 경우 (분자가 분모보다 클 경우)
	//		ctrl->RSM = 1;							// RSM = 1로 고정
	//}
	//else{
	//	ctrl->RSM = 1;								// roll_rate가 0인경우 롤 지표는 안정하다고 가정하고 RSM = 1 설정
	//}



	//////////// 피치 안정성 지표 //////////
	//
	//// Pitch Stability Metric (ver 4): 각속도와 수직 타이어력 사용

	//VF_FM = F_z[0] + F_z[1] + F_z[2] + F_z[3];	// Front + Middle 수직 타이어력의 합
	//VF_RM = F_z[2] + F_z[3] + F_z[4] + F_z[5];	//Rear + Middle 수직 타이어력의 합

	//for (k = 1; k < MAF_TAP; k++) {
	//	for (j = 0; j < 2; j++) {	// j = 0 : Front + Middle, j = 1 : Rear + Middle
	//		VF_PSM_buffer[k - 1][j] = VF_PSM_buffer[k][j];	// Moving average filter 계산을 위해 버퍼 내부 데이터 이동
	//	}
	//}
	//VF_PSM_buffer[MAF_TAP - 1][0] = VF_FM;	// 버퍼의 마지막 배열에 수직 타이어력 합 추가 (전방+중간 타이어)
	//VF_PSM_buffer[MAF_TAP - 1][1] = VF_RM;	// 버퍼의 마지막 배열에 수직 타이어력 합 추가 (후방+중간 타이어)

	//// 타이어력 평균 값 계산 (Moving average filter)
	//VF_PSM_mean[0] = 0;
	//VF_PSM_mean[1] = 0;

	//for (k = 0; k < MAF_TAP; k++) {
	//	for (j = 0; j < 2; j++) {
	//		VF_PSM_mean[j] += VF_PSM_buffer[k][j];	// 타이어력 평균 계산을 위해 각 버퍼에 저장된 데이터를 모두 더해준다.
	//	}
	//}

	//VF_PSM_mean[0] = VF_PSM_mean[0] / MAF_TAP;		// 타이어력의 평균값 = moving average filter, FM(uphill)
	//VF_PSM_mean[1] = VF_PSM_mean[1] / MAF_TAP;		// 타이어력의 평균값 = moving average filter, RM(downhill)

	//if (pitch_rate > rate_epsilon) {							// FM(uphill)
	//	ctrl->PSM = VF_PSM_mean[0] / VF_FM_max;		// 전방+중간 수직 타이어력의 비율을 PSM으로 계산
	//	if (VF_PSM_mean[0] > VF_FM_max)				// 만약 타이어력의 크기가 평상시보다 클 경우 (분자가 분모보다 클 경우)
	//		ctrl->PSM = 1;							// PSM = 1로 고정
	//}
	//else if (pitch_rate < -rate_epsilon) {						// RM(downhill)
	//	ctrl->PSM = VF_PSM_mean[1] / VF_RM_max;		// 후방+중간 수직 타이어력의 비율을 PSM으로 계산
	//	if (VF_PSM_mean[1] > VF_RM_max)				// 만약 타이어력의 크기가 평상시보다 클 경우 (분자가 분모보다 클 경우)
	//		ctrl->PSM = 1;							// PSM = 1로 고정
	//}
	//else {
	//	ctrl->PSM = 1;								// pitch_rate가 0인경우 롤 지표는 안정하다고 가정하고 RSM = 1 설정
	//}
	
	////////// 횡 방향 안정성 지표 //////////

	// Lateral Stability Metric
	ctrl->LSM = (lateral_position_max - fabs(ctrl->e_l)) / lateral_position_max;		// Lateral Position Error (normalized)


	////////// 수직 방향 안정성 지표 //////////
	
	// Vertical Acceleration Metric
	a_z = ctrl->ddr0c_p[2];	// local Z acceleration	
	ctrl->VSM = (a_z_max - fabs(a_z / 9.81)) / a_z_max;     // VSM이 1을 넘으면 fail
	


	////////// 안정성 지표 Pass/Fail 계산 //////////

	// stability_indx
	if (ctrl->RSM <= RSM_max)
		ctrl->RSM_indx = 0;    // RSM indx
	else
		ctrl->RSM_indx = 1;

	if (ctrl->PSM <= PSM_max)
		ctrl->PSM_indx = 0;    // PSM indx
	else
		ctrl->PSM_indx = 1;

	if (ctrl->LSM <= LSM_max)
		ctrl->LSM_indx = 0;    // LSM indx
	else
		ctrl->LSM_indx = 1;

	if (ctrl->VSM <= VSM_max)
		ctrl->VSM_indx = 0;    // VSM indx
	else
		ctrl->VSM_indx = 1;

	if (ctrl->WP_indx_update_flag == 1) {								// 차량이 새로운 Waypoint indx에 진입했을 때
		DID->RSM[ctrl->WP_indx - 1][thread_indx] = 1;
		DID->PSM[ctrl->WP_indx - 1][thread_indx] = 1;
		DID->LSM[ctrl->WP_indx - 1][thread_indx] = 1;
		DID->VSM[ctrl->WP_indx - 1][thread_indx] = 1;
	}

	// 주행 안정성 지표 결과 (RSM_indx, PSM_indx, LSM_indx, VSM_indx)를 WP 배열에 pass:1 / fail : 0 형태로 저장
	if (ctrl->RSM_indx && ctrl->PSM_indx && ctrl->LSM_indx && ctrl->VSM_indx == 1) {	// 4가지 인덱스가 모두 1(=Pass)인 경우
		if (ctrl->WP_indx_update_flag == 1) {								// 차량이 새로운 Waypoint indx에 진입했을 때
			//WP[ctrl->WP_indx - 1][3][thread_indx] = 1;						// 해당 Waypoint에서의 stability metric indx를 1(Pass)로 변경해줌 (맨 처음에는 모든 구간에서 0으로 초기화 되어 시뮬레이션 중단 시 stability metric을 fail 상태로 고정 시킴)
			DID->StabilityIndex[ctrl->WP_indx - 1][thread_indx] = 1;
			ctrl->WP_indx_update_flag = 0;									// 새로운 Waypoint indx 진입 flag 해제
		}
		//WP[ctrl->WP_indx - 1][3][thread_indx] = WP[ctrl->WP_indx - 1][3][thread_indx] * 1;    // i-1번째 waypoint에서의 주행 안정성 지표 모두 pass (1을 곱해줌)
		DID->StabilityIndex[ctrl->WP_indx - 1][thread_indx] *= 1;
	}
	else {																	// 4가지 인덱스 중 어느 하나라도 0(=Fail)인 경우
		//WP[ctrl->WP_indx - 1][3][thread_indx] = WP[ctrl->WP_indx - 1][3][thread_indx] * 0;    // fail	(0을 곱해줌)		
		DID->StabilityIndex[ctrl->WP_indx - 1][thread_indx] *= 0;
	}
	if (DID->RSM[ctrl->WP_indx - 1][thread_indx] >= ctrl->RSM)
		DID->RSM[ctrl->WP_indx - 1][thread_indx] = ctrl->RSM;
	if (DID->PSM[ctrl->WP_indx - 1][thread_indx] >= ctrl->PSM)
		DID->PSM[ctrl->WP_indx - 1][thread_indx] = ctrl->PSM;
	if (DID->LSM[ctrl->WP_indx - 1][thread_indx] >= ctrl->LSM)
		DID->LSM[ctrl->WP_indx - 1][thread_indx] = ctrl->LSM;
	if (DID->VSM[ctrl->WP_indx - 1][thread_indx] >= ctrl->VSM)
		DID->VSM[ctrl->WP_indx - 1][thread_indx] = ctrl->VSM;
}

void realtime_dynamics_analysis::wheel_spin_dyn(SUS *sus, CTRL *ctrl) {
	// 제어기 계산 결과인 토크를 index를 정렬해서 저장
	sus[RF].T_in = ctrl->motor_torque[1];
	sus[RM].T_in = ctrl->motor_torque[3];
	sus[RR].T_in = ctrl->motor_torque[5];
	sus[LF].T_in = ctrl->motor_torque[0];
	sus[LM].T_in = ctrl->motor_torque[2];
	sus[LR].T_in = ctrl->motor_torque[4];

	double Iyy = 21;   // 회전 축 관성 모멘트

	// dw = (T + My)/Iyy
	// My는 휠의 회전을 발생시키는 축 방향의 moment
	//    sus[RF].dw_wh = (sus[RF].T_in + sus[RF].My) / Iyy;
	//    sus[RM].dw_wh = (sus[RM].T_in + sus[RM].My) / Iyy;
	//    sus[RR].dw_wh = (sus[RR].T_in + sus[RR].My) / Iyy;
	//    sus[LF].dw_wh = (sus[LF].T_in + sus[LF].My) / Iyy;
	//    sus[LM].dw_wh = (sus[LM].T_in + sus[LM].My) / Iyy;
	//    sus[LR].dw_wh = (sus[LR].T_in + sus[LR].My) / Iyy;
	sus[RF].dw_wh = (sus[RF].T_in - sus[RF].Fx*sus[RF].R_d) / Iyy;
	sus[RM].dw_wh = (sus[RM].T_in - sus[RM].Fx*sus[RM].R_d) / Iyy;
	sus[RR].dw_wh = (sus[RR].T_in - sus[RR].Fx*sus[RR].R_d) / Iyy;
	sus[LF].dw_wh = (sus[LF].T_in - sus[LF].Fx*sus[LF].R_d) / Iyy;
	sus[LM].dw_wh = (sus[LM].T_in - sus[LM].Fx*sus[LM].R_d) / Iyy;
	sus[LR].dw_wh = (sus[LR].T_in - sus[LR].Fx*sus[LR].R_d) / Iyy;
}

void realtime_dynamics_analysis::dqddq2Yp(SIM *sim, CHASSIS *chassis, SUS sus[6]) {
	memcpy(sim->Yp, chassis->dr0, sizeof(double) * 3);
	memcpy(sim->Yp + 3, chassis->de0, sizeof(double) * 4);

	sim->Yp[7] = sus[RF].dq1;
	sim->Yp[8] = sus[RM].dq1;
	sim->Yp[9] = sus[RR].dq1;
	sim->Yp[10] = sus[LF].dq1;
	sim->Yp[11] = sus[LM].dq1;
	sim->Yp[12] = sus[LR].dq1;

	memcpy(sim->Yp + 13, chassis->ddr0, sizeof(double) * 3);
	memcpy(sim->Yp + 16, chassis->dw0, sizeof(double) * 3);

	sim->Yp[19] = sus[RF].ddq1;
	sim->Yp[20] = sus[RM].ddq1;
	sim->Yp[21] = sus[RR].ddq1;
	sim->Yp[22] = sus[LF].ddq1;
	sim->Yp[23] = sus[LM].ddq1;
	sim->Yp[24] = sus[LR].ddq1;

	sim->Yp[25] = sus[RF].dw_wh;
	sim->Yp[26] = sus[RM].dw_wh;
	sim->Yp[27] = sus[RR].dw_wh;
	sim->Yp[28] = sus[LF].dw_wh;
	sim->Yp[29] = sus[LM].dw_wh;
	sim->Yp[30] = sus[LR].dw_wh;
}

void realtime_dynamics_analysis::save_data(SIM *sim, CHASSIS *chassis, SUS sus[6], CTRL *ctrl) {

	int i;

	switch (sim->simulation_flag) {
		case 1:
			sprintf_s(sim->buf[0], sizeof(char)*255,"simulation_results/%d/data_sim_%d/_sim_flag1_Vertical_chassis_Cpp.dat", sim->main_count, sim_indx);
			sprintf_s(sim->buf[1], sizeof(char)*255,"simulation_results/%d/data_sim_%d/_sim_flag1_Longitudinal_chassis_Cpp.dat", sim->main_count, sim_indx);
			sprintf_s(sim->buf[2], sizeof(char)*255,"simulation_results/%d/data_sim_%d/_sim_flag1_Lateral_chassis_Cpp.dat", sim->main_count, sim_indx);
			sprintf_s(sim->buf[3], sizeof(char)*255,"simulation_results/%d/data_sim_%d/_sim_flag1_TSDA_force_Cpp.dat", sim->main_count, sim_indx);
			sprintf_s(sim->buf[4], sizeof(char)*255,"simulation_results/%d/data_sim_%d/_sim_flag1_R_P_Y_Cpp.dat", sim->main_count, sim_indx);
			sprintf_s(sim->buf[5], sizeof(char)*255,"simulation_results/%d/data_sim_%d/_sim_flag1_CHASSIS_X_Y_Cpp.dat", sim->main_count, sim_indx);
			sprintf_s(sim->buf[6], sizeof(char)*255,"simulation_results/%d/data_sim_%d/_sim_flag1_Motor_torque_Cpp.dat", sim->main_count, sim_indx);
			sprintf_s(sim->buf[7], sizeof(char)*255,"simulation_results/%d/data_sim_%d/_sim_flag1_LP_control_results_Cpp.dat", sim->main_count, sim_indx);
			sprintf_s(sim->buf[8], sizeof(char)*255,"simulation_results/%d/data_sim_%d/_sim_flag1_Tire_force_results_Cpp.dat", sim->main_count, sim_indx);
			sprintf_s(sim->buf[9], sizeof(char)*255, "simulation_results/%d/data_sim_%d/_sim_flag1_TSDA_K_C_Cpp.dat", sim->main_count, sim_indx);
			sprintf_s(sim->buf[10],sizeof(char)*255, "simulation_results/%d/data_sim_%d/_sim_flag1_VideoData_Cpp.dat", sim->main_count, sim_indx);
			sprintf_s(sim->buf[11],sizeof(char)*255, "simulation_results/%d/data_sim_%d/_sim_flag1_X_Y_Z_Cpp.dat", sim->main_count, sim_indx);
			sprintf_s(sim->buf[12],sizeof(char)*255, "simulation_results/%d/data_sim_%d/_sim_flag1_TIRE_alpha_Cpp.dat", sim->main_count, sim_indx);
			sprintf_s(sim->buf[13],sizeof(char)*255, "simulation_results/%d/data_sim_%d/_sim_flag1_TIRE_slip_Cpp.dat", sim->main_count, sim_indx);
			sprintf_s(sim->buf[14],sizeof(char)*255, "simulation_results/%d/data_sim_%d/_sim_flag1_Tire_X_force_results_Cpp.dat", sim->main_count, sim_indx);
			sprintf_s(sim->buf[15],sizeof(char)*255, "simulation_results/%d/data_sim_%d/_sim_flag1_Tire_pen_Cpp.dat", sim->main_count, sim_indx);
			sprintf_s(sim->buf[16],sizeof(char)*255, "simulation_results/%d/data_sim_%d/_sim_flag1_3D_Animation_Data_Cpp.dat", sim->main_count, sim_indx);
			sprintf_s(sim->buf[17],sizeof(char)*255, "simulation_results/%d/data_sim_%d/_sim_flag1_omega_Cpp.dat", sim->main_count, sim_indx);
			sprintf_s(sim->buf[18],sizeof(char)*255, "simulation_results/%d/data_sim_%d/_sim_flag1_Tire_Y_force_results_Cpp.dat", sim->main_count, sim_indx);
			sprintf_s(sim->buf[19],sizeof(char)*255, "simulation_results/%d/data_sim_%d/_sim_flag1_Tire_X_Moment_Cpp.dat", sim->main_count, sim_indx);
			sprintf_s(sim->buf[20],sizeof(char)*255, "simulation_results/%d/data_sim_%d/_sim_flag1_Tire_Y_Moment_Cpp.dat", sim->main_count, sim_indx);
			sprintf_s(sim->buf[21],sizeof(char)*255, "simulation_results/%d/data_sim_%d/_sim_flag1_Tire_Z_Moment_Cpp.dat", sim->main_count, sim_indx);

			if (sim->t_current == 0) {
				for (i = 0; i < 22; i++) {
					sim->fp[i] = fopen(sim->buf[i], "w+");
				}
			}
// 			else {
// 				for (i = 0; i < 22; i++) {
// 					sim->fp[i] = fopen(sim->buf[i], "a+");
// 				}
// 			}
			break;
		case 3:
			//sprintf_s(sim->buf[0], sizeof(char)*255,"simulation_results/%d/data_sim_%d/_sim_flag3_thread_%d_Vertical_chassis_Cpp.dat", sim->main_count, sim_indx, sim->thread_indx);
			//sprintf_s(sim->buf[1], sizeof(char)*255,"simulation_results/%d/data_sim_%d/_sim_flag3_thread_%d_Longitudinal_chassis_Cpp.dat", sim->main_count, sim_indx, sim->thread_indx);
			//sprintf_s(sim->buf[2], sizeof(char)*255,"simulation_results/%d/data_sim_%d/_sim_flag3_thread_%d_Lateral_chassis_Cpp.dat", sim->main_count, sim_indx, sim->thread_indx);
			//sprintf_s(sim->buf[3], sizeof(char)*255,"simulation_results/%d/data_sim_%d/_sim_flag3_thread_%d_TSDA_force_Cpp.dat", sim->main_count, sim_indx, sim->thread_indx);
			//sprintf_s(sim->buf[4], sizeof(char)*255,"simulation_results/%d/data_sim_%d/_sim_flag3_thread_%d_R_P_Y_Cpp.dat", sim->main_count, sim_indx, sim->thread_indx);
			//sprintf_s(sim->buf[5], sizeof(char)*255,"simulation_results/%d/data_sim_%d/_sim_flag3_thread_%d_CHASSIS_X_Y_Cpp.dat", sim->main_count, sim_indx, sim->thread_indx);
			//sprintf_s(sim->buf[6], sizeof(char)*255,"simulation_results/%d/data_sim_%d/_sim_flag3_thread_%d_Motor_torque_Cpp.dat", sim->main_count, sim_indx, sim->thread_indx);
			sprintf_s(sim->buf[7], sizeof(char)*255,"simulation_results/%d/data_sim_%d/_sim_flag3_thread_%d_LP_control_results_Cpp.dat", sim->main_count, sim_indx, sim->thread_indx);
			//sprintf_s(sim->buf[8], sizeof(char)*255,"simulation_results/%d/data_sim_%d/_sim_flag3_thread_%d_Tire_force_results_Cpp.dat", sim->main_count, sim_indx, sim->thread_indx);
			//sprintf_s(sim->buf[9], sizeof(char)*255,"simulation_results/%d/data_sim_%d/_sim_flag3_thread_%d_TSDA_K_C_Cpp.dat", sim->main_count, sim_indx, sim->thread_indx);
			//sprintf_s(sim->buf[10],sizeof(char)*255, "simulation_results/%d/data_sim_%d/_sim_flag3_thread_%d_VideoData_Cpp.dat", sim->main_count, sim_indx, sim->thread_indx);
			//sprintf_s(sim->buf[11],sizeof(char)*255, "simulation_results/%d/data_sim_%d/_sim_flag3_thread_%d_X_Y_Z_Cpp.dat", sim->main_count, sim_indx, sim->thread_indx);
			//sprintf_s(sim->buf[12],sizeof(char)*255, "simulation_results/%d/data_sim_%d/_sim_flag3_thread_%d_TIRE_alpha_Cpp.dat", sim->main_count, sim_indx, sim->thread_indx);
			//sprintf_s(sim->buf[13],sizeof(char)*255, "simulation_results/%d/data_sim_%d/_sim_flag3_thread_%d_TIRE_slip_Cpp.dat", sim->main_count, sim_indx, sim->thread_indx);
			//sprintf_s(sim->buf[14],sizeof(char)*255, "simulation_results/%d/data_sim_%d/_sim_flag3_thread_%d_Tire_X_force_results_Cpp.dat", sim->main_count, sim_indx, sim->thread_indx);
			//sprintf_s(sim->buf[15],sizeof(char)*255, "simulation_results/%d/data_sim_%d/_sim_flag3_thread_%d_omega_Cpp.dat", sim->main_count, sim_indx, sim->thread_indx);
			sprintf_s(sim->buf[16],sizeof(char)*255, "simulation_results/%d/data_sim_%d/_sim_flag3_thread_%d_Stability_index_Cpp.dat", sim->main_count, sim_indx, sim->thread_indx);
			//sprintf_s(sim->buf[17],sizeof(char)*255, "simulation_results/%d/data_sim_%d/_sim_flag3_thread_%d_3D_Animation_Data_Cpp.dat", sim->main_count, sim_indx, sim->thread_indx);
			//sprintf_s(sim->buf[18],sizeof(char)*255, "simulation_results/%d/data_sim_%d/_sim_flag3_thread_%d_Tire_Y_force_results_Cpp.dat", sim->main_count, sim_indx, sim->thread_indx);
			//sprintf_s(sim->buf[19],sizeof(char)*255, "simulation_results/%d/data_sim_%d/_sim_flag3_thread_%d_Tire_X_Moment_Cpp.dat", sim->main_count, sim_indx, sim->thread_indx);
			//sprintf_s(sim->buf[20],sizeof(char)*255, "simulation_results/%d/data_sim_%d/_sim_flag3_thread_%d_Tire_Y_Moment_Cpp.dat", sim->main_count, sim_indx, sim->thread_indx);
			//sprintf_s(sim->buf[21],sizeof(char)*255, "simulation_results/%d/data_sim_%d/_sim_flag3_thread_%d_Tire_Z_Moment_Cpp.dat", sim->main_count, sim_indx, sim->thread_indx);
			//sprintf_s(sim->buf[22],sizeof(char)*255, "simulation_results/%d/data_sim_%d/_sim_flag3_thread_%d_Tire_pen_Cpp.dat", sim->main_count, sim_indx, sim->thread_indx);

			if (sim->t_current == 0) {
				//for (i = 0; i < 23; i++) {
					//sim->fp[i] = fopen(sim->buf[i], "w+");
				fopen_s(&sim->fp[7], sim->buf[7], "w+");
				fopen_s(&sim->fp[16], sim->buf[16], "w+");

				//}
			}
// 			else {
// 				for (i = 0; i < 23; i++) {
// 					sim->fp[i] = fopen(sim->buf[i], "a+");
// 				}
// 			}
			break;
		
		default: break;
	}
	
	// 애니메이션 상에서 타이어의 회전을 표현하기 위해 angle 로 적분하는 부분
	if (sim->t_current == 0) {
		sus[RF].theta_wh = 0;
		sus[RM].theta_wh = 0;
		sus[RR].theta_wh = 0;
		sus[LF].theta_wh = 0;
		sus[LM].theta_wh = 0;
		sus[LR].theta_wh = 0;
	}
	else {
		sus[RF].theta_wh = sus[RF].theta_wh + sim->h*sus[RF].w_wh + 0.5*sim->h*sim->h*sus[RF].dw_wh;
		sus[RM].theta_wh = sus[RM].theta_wh + sim->h*sus[RM].w_wh + 0.5*sim->h*sim->h*sus[RM].dw_wh;
		sus[RR].theta_wh = sus[RR].theta_wh + sim->h*sus[RR].w_wh + 0.5*sim->h*sim->h*sus[RR].dw_wh;
		sus[LF].theta_wh = sus[LF].theta_wh + sim->h*sus[LF].w_wh + 0.5*sim->h*sim->h*sus[LF].dw_wh;
		sus[LM].theta_wh = sus[LM].theta_wh + sim->h*sus[LM].w_wh + 0.5*sim->h*sim->h*sus[LM].dw_wh;
		sus[LR].theta_wh = sus[LR].theta_wh + sim->h*sus[LR].w_wh + 0.5*sim->h*sim->h*sus[LR].dw_wh;
	}

	//fprintf_s(sim->fp[0], "%10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f\n", sim->t_current, chassis->r0c[2], chassis->dr0cp[2], chassis->ddr0cp[2], sus[RF].q1, sus[RM].q1, sus[RR].q1, sus[LF].q1, sus[LM].q1, sus[LR].q1);
	//fprintf_s(sim->fp[1], "%10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f\n", sim->t_current, chassis->r0c[0], chassis->dr0cp[0], chassis->ddr0cp[0], sus[RF].w_wh, sus[RM].w_wh, sus[RR].w_wh, sus[LF].w_wh, sus[LM].w_wh, sus[LR].w_wh);
	//fprintf_s(sim->fp[2], "%10.15f %10.15f %10.15f %10.15f %10.15f\n", sim->t_current, chassis->r0c[1], chassis->dr0cp[1], chassis->ddr0cp[1], ctrl->ddr0c_p[1]);
	//fprintf_s(sim->fp[3], "%10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f\n", sim->t_current, sus[RF].f, sus[RM].f, sus[RR].f, sus[LF].f, sus[LM].f, sus[LR].f);
	//fprintf_s(sim->fp[4], "%10.15f %10.15f %10.15f %10.15f \n", sim->t_current, chassis->roll_ang, chassis->pitch_ang, chassis->yaw_ang);
	//fprintf_s(sim->fp[5], "%10.15f %10.15f %10.15f \n", sim->t_current, chassis->r0c[0], chassis->r0c[1]);
	//fprintf_s(sim->fp[6], "%10.15f %10.15f %10.15f %10.15f %10.5f %10.15f %10.15f \n", sim->t_current, ctrl->motor_torque[1], ctrl->motor_torque[3], ctrl->motor_torque[5], ctrl->motor_torque[0], ctrl->motor_torque[2], ctrl->motor_torque[4]);
	fprintf_s(sim->fp[7], "%10.15f %10.15f %10.15f %10.15f %10.15f %4d %10.15f %10.15f\n", sim->t_current, ctrl->e_l, ctrl->e_psi, ctrl->v_di, ctrl->v_x, sim_indx, ctrl->F_xd_total, ctrl->v_di);
	//fprintf_s(sim->fp[8], "%10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f\n", sim->t_current, sus[RF].Fz, sus[RM].Fz, sus[RR].Fz, sus[LF].Fz, sus[LM].Fz, sus[LR].Fz);
	//fprintf_s(sim->fp[9], "%10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f\n", sim->t_current, sus[RF].T_spring, sus[RM].T_spring, sus[RR].T_spring, sus[LF].T_spring, sus[LM].T_spring, sus[LR].T_spring, sus[RF].T_damper, sus[RM].T_damper, sus[RR].T_damper, sus[LF].T_damper, sus[LM].T_damper, sus[LR].T_damper);
	//fprintf_s(sim->fp[10], "%10.15f %10.15f %10.15f %10.15f\n", sim->t_current, chassis->r0c[0], chassis->r0c[1], chassis->yaw_ang);
	//fprintf_s(sim->fp[11], "%10.15f %10.15f %10.15f %10.15f\n", sim->t_current, chassis->r0c[0], chassis->r0c[1], chassis->r0c[2]);
	//fprintf_s(sim->fp[12], "%10.15f %10.15f %10.15f %10.15f %10.5f %10.15f %10.15f \n", sim->t_current, sus[RF].angle, sus[RM].angle, sus[RR].angle, sus[LF].angle, sus[LM].angle, sus[LR].angle);
	//fprintf_s(sim->fp[13], "%10.15f %10.15f %10.15f %10.15f %10.5f %10.15f %10.15f \n", sim->t_current, sus[RF].slip, sus[RM].slip, sus[RR].slip, sus[LF].slip, sus[LM].slip, sus[LR].slip);
	//fprintf_s(sim->fp[14], "%10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f\n", sim->t_current, sus[RF].Fx, sus[RM].Fx, sus[RR].Fx, sus[LF].Fx, sus[LM].Fx, sus[LR].Fx);
	if (sim->simulation_flag == 3 || sim->simulation_flag == 5) {
		//fprintf_s(sim->fp[15], "%10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f\n", sim->t_current, sim->Y[25], sim->Y[26], sim->Y[27], sim->Y[28], sim->Y[29], sim->Y[30]);
		fprintf_s(sim->fp[16], "%10.15f %d %10.15f %d %10.15f %d %10.5f %d %10.15f %d \n", sim->t_current, ctrl->WP_indx, ctrl->RSM, ctrl->RSM_indx, ctrl->PSM, ctrl->PSM_indx, ctrl->LSM, ctrl->LSM_indx, ctrl->VSM, ctrl->VSM_indx);
		//fprintf_s(sim->fp[17], "%10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f\n",
			//sim->t_current, chassis->r0[0], chassis->r0[1], chassis->r0[2], chassis->A0[0][0], chassis->A0[0][1], chassis->A0[0][2], chassis->A0[1][0], chassis->A0[1][1], chassis->A0[1][2], chassis->A0[2][0], chassis->A0[2][1], chassis->A0[2][2], sus[RF].q1, sus[RM].q1, sus[RR].q1, sus[LF].q1, sus[LM].q1, sus[LR].q1, sus[RF].theta_wh, sus[RM].theta_wh, sus[RR].theta_wh, sus[LF].theta_wh, sus[LM].theta_wh, sus[LR].theta_wh, chassis->roll_ang, chassis->pitch_ang, chassis->yaw_ang);
		//fprintf_s(sim->fp[18], "%10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f\n", sim->t_current, sus[RF].Fy, sus[RM].Fy, sus[RR].Fy, sus[LF].Fy, sus[LM].Fy, sus[LR].Fy);
		//fprintf_s(sim->fp[19], "%10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f\n", sim->t_current, sus[RF].M_tire[0], sus[RM].M_tire[0], sus[RR].M_tire[0], sus[LF].M_tire[0], sus[LM].M_tire[0], sus[LR].M_tire[0]);
		//fprintf_s(sim->fp[20], "%10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f\n", sim->t_current, sus[RF].M_tire[1], sus[RM].M_tire[1], sus[RR].M_tire[1], sus[LF].M_tire[1], sus[LM].M_tire[1], sus[LR].M_tire[1]);
		//fprintf_s(sim->fp[21], "%10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f\n", sim->t_current, sus[RF].M_tire[2], sus[RM].M_tire[2], sus[RR].M_tire[2], sus[LF].M_tire[2], sus[LM].M_tire[2], sus[LR].M_tire[2]);
		//fprintf_s(sim->fp[22], "%10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f\n", sim->t_current, sus[RF].pen, sus[RM].pen, sus[RR].pen, sus[LF].pen, sus[LM].pen, sus[LR].pen);
// 		fclose(fp[15]); //omega.dat
// 		fclose(fp[16]); //Stability_index.dat
// 		fclose(fp[17]); //3D_Animation_Data.dat
// 		fclose(fp[18]); //Tire_Y_force_results.dat
// 		fclose(fp[19]); //Tire_X_Moment.dat
// 		fclose(fp[20]); //Tire_Y_Moment.dat
// 		fclose(fp[21]); //Tire_Z_Moment.dat
// 		fclose(fp[22]); //Tire pen.dat
	}
	else {
		fprintf_s(sim->fp[15], "%10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f\n", sim->t_current, sus[RF].pen, sus[RM].pen, sus[RR].pen, sus[LF].pen, sus[LM].pen, sus[LR].pen);
		fprintf_s(sim->fp[16], "%10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f\n",
			sim->t_current, chassis->r0[0], chassis->r0[1], chassis->r0[2], chassis->A0[0][0], chassis->A0[0][1], chassis->A0[0][2], chassis->A0[1][0], chassis->A0[1][1], chassis->A0[1][2], chassis->A0[2][0], chassis->A0[2][1], chassis->A0[2][2], sus[RF].q1, sus[RM].q1, sus[RR].q1, sus[LF].q1, sus[LM].q1, sus[LR].q1, sus[RF].theta_wh, sus[RM].theta_wh, sus[RR].theta_wh, sus[LF].theta_wh, sus[LM].theta_wh, sus[LR].theta_wh, chassis->roll_ang, chassis->pitch_ang, chassis->yaw_ang);
		fprintf_s(sim->fp[17], "%10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f\n", sim->t_current, sim->Y[25], sim->Y[26], sim->Y[27], sim->Y[28], sim->Y[29], sim->Y[30]);
		fprintf_s(sim->fp[18], "%10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f\n", sim->t_current, sus[RF].Fy, sus[RM].Fy, sus[RR].Fy, sus[LF].Fy, sus[LM].Fy, sus[LR].Fy);
		fprintf_s(sim->fp[19], "%10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f\n", sim->t_current, sus[RF].M_tire[0], sus[RM].M_tire[0], sus[RR].M_tire[0], sus[LF].M_tire[0], sus[LM].M_tire[0], sus[LR].M_tire[0]);
		fprintf_s(sim->fp[20], "%10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f\n", sim->t_current, sus[RF].M_tire[1], sus[RM].M_tire[1], sus[RR].M_tire[1], sus[LF].M_tire[1], sus[LM].M_tire[1], sus[LR].M_tire[1]);
		fprintf_s(sim->fp[21], "%10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f\n", sim->t_current, sus[RF].M_tire[2], sus[RM].M_tire[2], sus[RR].M_tire[2], sus[LF].M_tire[2], sus[LM].M_tire[2], sus[LR].M_tire[2]);
// 		fclose(fp[15]); //Tire pen.dat
// 		fclose(fp[16]); //3D_Animation_Data.dat
// 		fclose(fp[17]); //omega.dat
// 		fclose(fp[18]); //Tire_Y_force_results.dat
// 		fclose(fp[19]); //Tire_X_Moment.dat
// 		fclose(fp[20]); //Tire_Y_Moment.dat
// 		fclose(fp[21]); //Tire_Z_Moment.dat
	}

// 	fclose(fp[0]);	//Vertical_chassis.dat
// 	fclose(fp[1]);	//Longitudinal_chassis.dat
// 	fclose(fp[2]);	//Lateral_chassis.dat
// 	fclose(fp[3]);	//TSDA_force.dat
// 	fclose(fp[4]);	//R_P_Y.dat
// 	fclose(fp[5]);	//CHASSIS_X_Y.dat
// 	fclose(fp[6]);	//Motor_torque.dat
// 	fclose(fp[7]);	//LP_control_results.dat
// 	fclose(fp[8]);	//Tire_force_results.dat
// 	fclose(fp[9]);	//TSDA_K_C.dat
// 	fclose(fp[10]);	//VideoData.dat
// 	fclose(fp[11]);	//X_Y_Z.dat
// 	fclose(fp[12]);	//TIRE_alpha.dat
// 	fclose(fp[13]);	//TIRE_slip.dat
// 	fclose(fp[14]);	//Tire_X_force_results.dat
}

void realtime_dynamics_analysis::absh3(double step_size, const int n, int *intcount, double AW1[31][2], double AW[31][2], double *t_current, double Y[31], double Yp[31], double Y_next[31], double *t_next)
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

void realtime_dynamics_analysis::read_system(SIM *sim) {
	// start_time, end_time 은 시뮬레이션 상에서 사용하지 않으므로 없음.
	sim->h = 0.01;
	sim->g = -9.80665;
}

void realtime_dynamics_analysis::read_chassis(SIM *sim, CHASSIS *chassis, SUS sus[6], double R_u) {
	const double CH[35] = { 
		0, 0, 0, /*r0*/
		1, 0, 0, 0, /*p0*/
		0, 0, 0, /*dr0*/
		0, 0, 0, /*w0*/
		0, 0, 0, /*rho0p*/
		1, 0, 0, 0, 1, 0, 0, 0, 1, /*C00*/
		4362, /*m0*/
		2046.8, 0, 0, 0, 7113.2, 0, 0, 0, 7474.2 /*J0p*/
	};

	chassis->r0[0] = DID->LocalPath[0][0];//WP[0][0][0]; // 초기 글로벌 X;
	chassis->r0[1] = DID->LocalPath[0][1];//WP[0][1][0]; // 초기 글로벌 Y;
	chassis->r0[2] = 0;

	chassis->e0[0] = sqrt(2 * cos(sim->heading) + 2) / 2;
	chassis->e0[1] = 0;
	chassis->e0[2] = 0;
	chassis->e0[3] = (sin(sim->heading) * 2) / (4 * chassis->e0[0]);

	chassis->A0[0][0] = 2 * (chassis->e0[0] * chassis->e0[0] + chassis->e0[1] * chassis->e0[1] - 0.5);
	chassis->A0[0][1] = 2 * (chassis->e0[1] * chassis->e0[2] - chassis->e0[0] * chassis->e0[3]);
	chassis->A0[0][2] = 2 * (chassis->e0[1] * chassis->e0[3] + chassis->e0[0] * chassis->e0[2]);
	chassis->A0[1][0] = 2 * (chassis->e0[1] * chassis->e0[2] + chassis->e0[0] * chassis->e0[3]);
	chassis->A0[1][1] = 2 * (chassis->e0[0] * chassis->e0[0] + chassis->e0[2] * chassis->e0[2] - 0.5);
	chassis->A0[1][2] = 2 * (chassis->e0[2] * chassis->e0[3] - chassis->e0[0] * chassis->e0[1]);
	chassis->A0[2][0] = 2 * (chassis->e0[1] * chassis->e0[3] - chassis->e0[0] * chassis->e0[2]);
	chassis->A0[2][1] = 2 * (chassis->e0[2] * chassis->e0[3] + chassis->e0[0] * chassis->e0[1]);
	chassis->A0[2][2] = 2 * (chassis->e0[0] * chassis->e0[0] + chassis->e0[3] * chassis->e0[3] - 0.5);

	read_RF(&sus[RF], chassis);
	read_RM(&sus[RM], chassis);
	read_RR(&sus[RR], chassis);
	read_LF(&sus[LF], chassis);
	read_LM(&sus[LM], chassis);
	read_LR(&sus[LR], chassis);

	// 차량 초기 자세에서 각 타이어 별 침투량 계산
	for (int i = 0; i < 6; i++) {
		//Obj_Map::PNU_fn_map(sus[i].rw[0], sus[i].rw[1], &sus[i].road_h);
		Obj_Map_Global::PNU_fn_map(sus[i].rw[0], sus[i].rw[1], &sus[i].road_h);
		sus[i].pen = sus[i].road_h - (sus->rw[2] - R_u);
	}

	// 차체의 높이를 설정하기 위해 각 타이어 침투량 중 최대값을 찾는 routine
	double max_pen = 0;
	for (int i = 0; i < 6; i++) {
		if (max_pen < sus[i].pen) max_pen = sus[i].pen;
	}
	chassis->r0[2] = max_pen;

	chassis->dr0[0] = CH[7]; chassis->dr0[1] = CH[8]; chassis->dr0[2] = CH[9];

	chassis->w0[0] = CH[10]; chassis->w0[1] = CH[11]; chassis->w0[2] = CH[12];

	chassis->rho0p[0] = CH[13]; chassis->rho0p[1] = CH[14]; chassis->rho0p[2] = CH[15];

	chassis->C00[0][0] = CH[16]; chassis->C00[0][1] = CH[17]; chassis->C00[0][2] = CH[18];
	chassis->C00[1][0] = CH[19]; chassis->C00[1][1] = CH[20]; chassis->C00[1][2] = CH[21];
	chassis->C00[2][0] = CH[22]; chassis->C00[2][1] = CH[23]; chassis->C00[2][2] = CH[24];

	chassis->m0 = CH[25];

	chassis->J0[0][0] = CH[26]; chassis->J0[0][1] = CH[27]; chassis->J0[0][2] = CH[28];
	chassis->J0[1][0] = CH[29]; chassis->J0[1][1] = CH[30]; chassis->J0[1][2] = CH[31];
	chassis->J0[2][0] = CH[32]; chassis->J0[2][1] = CH[33]; chassis->J0[2][2] = CH[34];
}

void realtime_dynamics_analysis::read_RF(SUS *sus, CHASSIS *chassis) {
	const double RF_in[54] = { 
		0, /*q1*/
		0, /*dq1*/
		0.586234043, 0, 9.32E-02, /*rho1p*/
		-6.48E-12, 0.865112866, -0.501577241, 1, 3.04E-12, -7.68E-12, -5.12E-12, -0.501577241, -0.865112866, /*C11*/
		376, /*m1*/
		31.25706362, 0, 0, 0, 25.98627239, 0, 0, 0, 20.83704122, /*J1p*/
		1.13, -0.755, -0.108, /*s01p*/
		0.65, 0, 0.219, /*s12p*/
		1, 0, 0, 0, -5.10E-12, -1, 0, 1, -5.10E-12, /*C01*/
		1, 0, 0, 0, 1, 0, 0, 0, 1, /*C12*/
		1.276, -0.736, 0.658, /*s0sp*/
		0.281, 0.133, -1.90E-02 /*s1sp*/
	};

	sus->q1 = RF_in[0];
	sus->dq1 = RF_in[1];

	sus->rho1p[0] = RF_in[2]; sus->rho1p[1] = RF_in[3]; sus->rho1p[2] = RF_in[4];

	sus->C11[0][0] = RF_in[5]; sus->C11[0][1] = RF_in[6]; sus->C11[0][2] = RF_in[7];
	sus->C11[1][0] = RF_in[8]; sus->C11[1][1] = RF_in[9]; sus->C11[1][2] = RF_in[10];
	sus->C11[2][0] = RF_in[11]; sus->C11[2][1] = RF_in[12]; sus->C11[2][2] = RF_in[13];

	sus->m1 = RF_in[14];

	sus->J1[0][0] = RF_in[15]; sus->J1[0][1] = RF_in[16]; sus->J1[0][2] = RF_in[17];
	sus->J1[1][0] = RF_in[18]; sus->J1[1][1] = RF_in[19]; sus->J1[1][2] = RF_in[20];
	sus->J1[2][0] = RF_in[21]; sus->J1[2][1] = RF_in[22]; sus->J1[2][2] = RF_in[23];

	sus->s01p[0] = RF_in[24]; sus->s01p[1] = RF_in[25]; sus->s01p[2] = RF_in[26];

	sus->s12p[0] = RF_in[27]; sus->s12p[1] = RF_in[28]; sus->s12p[2] = RF_in[29];

	sus->C01[0][0] = RF_in[30]; sus->C01[0][1] = RF_in[31]; sus->C01[0][2] = RF_in[32];
	sus->C01[1][0] = RF_in[33]; sus->C01[1][1] = RF_in[34]; sus->C01[1][2] = RF_in[35];
	sus->C01[2][0] = RF_in[36]; sus->C01[2][1] = RF_in[37]; sus->C01[2][2] = RF_in[38];

	sus->C12[0][0] = RF_in[39]; sus->C12[0][1] = RF_in[40]; sus->C12[0][2] = RF_in[41];
	sus->C12[1][0] = RF_in[42]; sus->C12[1][1] = RF_in[43]; sus->C12[1][2] = RF_in[44];
	sus->C12[2][0] = RF_in[45]; sus->C12[2][1] = RF_in[46]; sus->C12[2][2] = RF_in[47];

	sus->s0sp[0] = RF_in[48]; sus->s0sp[1] = RF_in[49]; sus->s0sp[2] = RF_in[50];

	sus->s1sp[0] = RF_in[51]; sus->s1sp[1] = RF_in[52]; sus->s1sp[2] = RF_in[53];

	sus->theta_wh = 0; sus->w_wh = 0; sus->dw_wh = 0;

	sus->rw[0] = chassis->A0[0][0] * sus->s01p[0] + chassis->A0[0][1] * sus->s01p[1] + chassis->A0[0][2] * sus->s01p[2] + chassis->r0[0] + (chassis->A0[0][0] * sus->C01[0][0] + chassis->A0[0][1] * sus->C01[1][0] + chassis->A0[0][2] * sus->C01[2][0]) * sus->s12p[0] + (chassis->A0[0][0] * sus->C01[0][1] + chassis->A0[0][1] * sus->C01[1][1] + chassis->A0[0][2] * sus->C01[2][1]) * sus->s12p[1] + (chassis->A0[0][0] * sus->C01[0][2] + chassis->A0[0][1] * sus->C01[1][2] + chassis->A0[0][2] * sus->C01[2][2]) * sus->s12p[2];
	sus->rw[1] = chassis->A0[1][0] * sus->s01p[0] + chassis->A0[1][1] * sus->s01p[1] + chassis->A0[1][2] * sus->s01p[2] + chassis->r0[1] + (chassis->A0[1][0] * sus->C01[0][0] + chassis->A0[1][1] * sus->C01[1][0] + chassis->A0[1][2] * sus->C01[2][0]) * sus->s12p[0] + (chassis->A0[1][0] * sus->C01[0][1] + chassis->A0[1][1] * sus->C01[1][1] + chassis->A0[1][2] * sus->C01[2][1]) * sus->s12p[1] + (chassis->A0[1][0] * sus->C01[0][2] + chassis->A0[1][1] * sus->C01[1][2] + chassis->A0[1][2] * sus->C01[2][2]) * sus->s12p[2];
	sus->rw[2] = chassis->A0[2][0] * sus->s01p[0] + chassis->A0[2][1] * sus->s01p[1] + chassis->A0[2][2] * sus->s01p[2] + chassis->r0[2] + (chassis->A0[2][0] * sus->C01[0][0] + chassis->A0[2][1] * sus->C01[1][0] + chassis->A0[2][2] * sus->C01[2][0]) * sus->s12p[0] + (chassis->A0[2][0] * sus->C01[0][1] + chassis->A0[2][1] * sus->C01[1][1] + chassis->A0[2][2] * sus->C01[2][1]) * sus->s12p[1] + (chassis->A0[2][0] * sus->C01[0][2] + chassis->A0[2][1] * sus->C01[1][2] + chassis->A0[2][2] * sus->C01[2][2]) * sus->s12p[2];
}

void realtime_dynamics_analysis::read_RM(SUS *sus, CHASSIS *chassis) {
	const double RM_in[54] = { 0, 0, 0.585382979, 0, 9.32E-02, 1.53E-11, -0.868782601, -0.495193692, -1, -1.33E-11, -7.58E-12, 0, 0.495193692, -0.868782601, 376, 31.21662106, 0, 0, 0, 25.94143616, 0, 0, 0, 20.8414349, 0.83, -0.755, -0.108, 0.648, 0, 0.219, -1, 1.02E-11, 0, 5.21E-23, 5.10E-12, -1, -1.02E-11, -1, -5.10E-12, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0.683, -0.736, 0.658, 0.281, -0.133, -1.90E-02 };

	sus->q1 = RM_in[0];
	sus->dq1 = RM_in[1];

	sus->rho1p[0] = RM_in[2]; sus->rho1p[1] = RM_in[3]; sus->rho1p[2] = RM_in[4];

	sus->C11[0][0] = RM_in[5]; sus->C11[0][1] = RM_in[6]; sus->C11[0][2] = RM_in[7];
	sus->C11[1][0] = RM_in[8]; sus->C11[1][1] = RM_in[9]; sus->C11[1][2] = RM_in[10];
	sus->C11[2][0] = RM_in[11]; sus->C11[2][1] = RM_in[12]; sus->C11[2][2] = RM_in[13];

	sus->m1 = RM_in[14];

	sus->J1[0][0] = RM_in[15]; sus->J1[0][1] = RM_in[16]; sus->J1[0][2] = RM_in[17];
	sus->J1[1][0] = RM_in[18]; sus->J1[1][1] = RM_in[19]; sus->J1[1][2] = RM_in[20];
	sus->J1[2][0] = RM_in[21]; sus->J1[2][1] = RM_in[22]; sus->J1[2][2] = RM_in[23];

	sus->s01p[0] = RM_in[24]; sus->s01p[1] = RM_in[25]; sus->s01p[2] = RM_in[26];

	sus->s12p[0] = RM_in[27]; sus->s12p[1] = RM_in[28]; sus->s12p[2] = RM_in[29];

	sus->C01[0][0] = RM_in[30]; sus->C01[0][1] = RM_in[31]; sus->C01[0][2] = RM_in[32];
	sus->C01[1][0] = RM_in[33]; sus->C01[1][1] = RM_in[34]; sus->C01[1][2] = RM_in[35];
	sus->C01[2][0] = RM_in[36]; sus->C01[2][1] = RM_in[37]; sus->C01[2][2] = RM_in[38];

	sus->C12[0][0] = RM_in[39]; sus->C12[0][1] = RM_in[40]; sus->C12[0][2] = RM_in[41];
	sus->C12[1][0] = RM_in[42]; sus->C12[1][1] = RM_in[43]; sus->C12[1][2] = RM_in[44];
	sus->C12[2][0] = RM_in[45]; sus->C12[2][1] = RM_in[46]; sus->C12[2][2] = RM_in[47];

	sus->s0sp[0] = RM_in[48]; sus->s0sp[1] = RM_in[49]; sus->s0sp[2] = RM_in[50];

	sus->s1sp[0] = RM_in[51]; sus->s1sp[1] = RM_in[52]; sus->s1sp[2] = RM_in[53];

	sus->theta_wh = 0; sus->w_wh = 0; sus->dw_wh = 0;

	sus->rw[0] = chassis->A0[0][0] * sus->s01p[0] + chassis->A0[0][1] * sus->s01p[1] + chassis->A0[0][2] * sus->s01p[2] + chassis->r0[0] + (chassis->A0[0][0] * sus->C01[0][0] + chassis->A0[0][1] * sus->C01[1][0] + chassis->A0[0][2] * sus->C01[2][0]) * sus->s12p[0] + (chassis->A0[0][0] * sus->C01[0][1] + chassis->A0[0][1] * sus->C01[1][1] + chassis->A0[0][2] * sus->C01[2][1]) * sus->s12p[1] + (chassis->A0[0][0] * sus->C01[0][2] + chassis->A0[0][1] * sus->C01[1][2] + chassis->A0[0][2] * sus->C01[2][2]) * sus->s12p[2];
	sus->rw[1] = chassis->A0[1][0] * sus->s01p[0] + chassis->A0[1][1] * sus->s01p[1] + chassis->A0[1][2] * sus->s01p[2] + chassis->r0[1] + (chassis->A0[1][0] * sus->C01[0][0] + chassis->A0[1][1] * sus->C01[1][0] + chassis->A0[1][2] * sus->C01[2][0]) * sus->s12p[0] + (chassis->A0[1][0] * sus->C01[0][1] + chassis->A0[1][1] * sus->C01[1][1] + chassis->A0[1][2] * sus->C01[2][1]) * sus->s12p[1] + (chassis->A0[1][0] * sus->C01[0][2] + chassis->A0[1][1] * sus->C01[1][2] + chassis->A0[1][2] * sus->C01[2][2]) * sus->s12p[2];
	sus->rw[2] = chassis->A0[2][0] * sus->s01p[0] + chassis->A0[2][1] * sus->s01p[1] + chassis->A0[2][2] * sus->s01p[2] + chassis->r0[2] + (chassis->A0[2][0] * sus->C01[0][0] + chassis->A0[2][1] * sus->C01[1][0] + chassis->A0[2][2] * sus->C01[2][0]) * sus->s12p[0] + (chassis->A0[2][0] * sus->C01[0][1] + chassis->A0[2][1] * sus->C01[1][1] + chassis->A0[2][2] * sus->C01[2][1]) * sus->s12p[1] + (chassis->A0[2][0] * sus->C01[0][2] + chassis->A0[2][1] * sus->C01[1][2] + chassis->A0[2][2] * sus->C01[2][2]) * sus->s12p[2];
}

void realtime_dynamics_analysis::read_RR(SUS *sus, CHASSIS *chassis) {
	const double RR_in[54] = { 0, 0, 0.585382979, 0, 0.0932, 0.0000000000153, -0.868782601, -0.495193692, -1, -0.0000000000133, -0.00000000000758, 0, 0.495193692, -0.868782601, 376, 31.21662106, 0, 0, 0, 25.94143616, 0, 0, 0, 20.8414349, -0.77, -0.755, -0.108, 0.648, 0, 0.219, -1, 0.0000000000102, 0, 5.21E-23, 0.0000000000051, -1, -0.0000000000102, -1, -0.0000000000051, 1, 0, 0, 0, 1, 0, 0, 0, 1, -0.917, -0.736, 0.658, 0.281, -0.133, -0.019};

	sus->q1 = RR_in[0];
	sus->dq1 = RR_in[1];

	sus->rho1p[0] = RR_in[2]; sus->rho1p[1] = RR_in[3]; sus->rho1p[2] = RR_in[4];

	sus->C11[0][0] = RR_in[5]; sus->C11[0][1] = RR_in[6]; sus->C11[0][2] = RR_in[7];
	sus->C11[1][0] = RR_in[8]; sus->C11[1][1] = RR_in[9]; sus->C11[1][2] = RR_in[10];
	sus->C11[2][0] = RR_in[11]; sus->C11[2][1] = RR_in[12]; sus->C11[2][2] = RR_in[13];

	sus->m1 = RR_in[14];

	sus->J1[0][0] = RR_in[15]; sus->J1[0][1] = RR_in[16]; sus->J1[0][2] = RR_in[17];
	sus->J1[1][0] = RR_in[18]; sus->J1[1][1] = RR_in[19]; sus->J1[1][2] = RR_in[20];
	sus->J1[2][0] = RR_in[21]; sus->J1[2][1] = RR_in[22]; sus->J1[2][2] = RR_in[23];

	sus->s01p[0] = RR_in[24]; sus->s01p[1] = RR_in[25]; sus->s01p[2] = RR_in[26];

	sus->s12p[0] = RR_in[27]; sus->s12p[1] = RR_in[28]; sus->s12p[2] = RR_in[29];

	sus->C01[0][0] = RR_in[30]; sus->C01[0][1] = RR_in[31]; sus->C01[0][2] = RR_in[32];
	sus->C01[1][0] = RR_in[33]; sus->C01[1][1] = RR_in[34]; sus->C01[1][2] = RR_in[35];
	sus->C01[2][0] = RR_in[36]; sus->C01[2][1] = RR_in[37]; sus->C01[2][2] = RR_in[38];

	sus->C12[0][0] = RR_in[39]; sus->C12[0][1] = RR_in[40]; sus->C12[0][2] = RR_in[41];
	sus->C12[1][0] = RR_in[42]; sus->C12[1][1] = RR_in[43]; sus->C12[1][2] = RR_in[44];
	sus->C12[2][0] = RR_in[45]; sus->C12[2][1] = RR_in[46]; sus->C12[2][2] = RR_in[47];

	sus->s0sp[0] = RR_in[48]; sus->s0sp[1] = RR_in[49]; sus->s0sp[2] = RR_in[50];

	sus->s1sp[0] = RR_in[51]; sus->s1sp[1] = RR_in[52]; sus->s1sp[2] = RR_in[53];

	sus->theta_wh = 0; sus->w_wh = 0; sus->dw_wh = 0;

	sus->rw[0] = chassis->A0[0][0] * sus->s01p[0] + chassis->A0[0][1] * sus->s01p[1] + chassis->A0[0][2] * sus->s01p[2] + chassis->r0[0] + (chassis->A0[0][0] * sus->C01[0][0] + chassis->A0[0][1] * sus->C01[1][0] + chassis->A0[0][2] * sus->C01[2][0]) * sus->s12p[0] + (chassis->A0[0][0] * sus->C01[0][1] + chassis->A0[0][1] * sus->C01[1][1] + chassis->A0[0][2] * sus->C01[2][1]) * sus->s12p[1] + (chassis->A0[0][0] * sus->C01[0][2] + chassis->A0[0][1] * sus->C01[1][2] + chassis->A0[0][2] * sus->C01[2][2]) * sus->s12p[2];
	sus->rw[1] = chassis->A0[1][0] * sus->s01p[0] + chassis->A0[1][1] * sus->s01p[1] + chassis->A0[1][2] * sus->s01p[2] + chassis->r0[1] + (chassis->A0[1][0] * sus->C01[0][0] + chassis->A0[1][1] * sus->C01[1][0] + chassis->A0[1][2] * sus->C01[2][0]) * sus->s12p[0] + (chassis->A0[1][0] * sus->C01[0][1] + chassis->A0[1][1] * sus->C01[1][1] + chassis->A0[1][2] * sus->C01[2][1]) * sus->s12p[1] + (chassis->A0[1][0] * sus->C01[0][2] + chassis->A0[1][1] * sus->C01[1][2] + chassis->A0[1][2] * sus->C01[2][2]) * sus->s12p[2];
	sus->rw[2] = chassis->A0[2][0] * sus->s01p[0] + chassis->A0[2][1] * sus->s01p[1] + chassis->A0[2][2] * sus->s01p[2] + chassis->r0[2] + (chassis->A0[2][0] * sus->C01[0][0] + chassis->A0[2][1] * sus->C01[1][0] + chassis->A0[2][2] * sus->C01[2][0]) * sus->s12p[0] + (chassis->A0[2][0] * sus->C01[0][1] + chassis->A0[2][1] * sus->C01[1][1] + chassis->A0[2][2] * sus->C01[2][1]) * sus->s12p[1] + (chassis->A0[2][0] * sus->C01[0][2] + chassis->A0[2][1] * sus->C01[1][2] + chassis->A0[2][2] * sus->C01[2][2]) * sus->s12p[2];
}

void realtime_dynamics_analysis::read_LF(SUS *sus, CHASSIS *chassis) {
	const double LF_in[54] = { 0, 0, 0.585864863, 0, -0.0932, -0.0000000000051, 0.866708022, 0.498815802, 1, 0.00000000000442, 0.00000000000255, 0, 0.498815802, -0.866708022, 376, 31.23950966, 0, 0, 0, 25.96680966, 0, 0, 0, 20.83895001, 1.13, 0.755, -0.108, 0.65, 0, -0.219, 1, 0, 0, 0, -0.0000000000051, -1, 0, 1, -0.0000000000051, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1.276, 0.736, 0.658, 0.281, 0.133, 0.019};

	sus->q1 = LF_in[0];
	sus->dq1 = LF_in[1];

	sus->rho1p[0] = LF_in[2]; sus->rho1p[1] = LF_in[3]; sus->rho1p[2] = LF_in[4];

	sus->C11[0][0] = LF_in[5]; sus->C11[0][1] = LF_in[6]; sus->C11[0][2] = LF_in[7];
	sus->C11[1][0] = LF_in[8]; sus->C11[1][1] = LF_in[9]; sus->C11[1][2] = LF_in[10];
	sus->C11[2][0] = LF_in[11]; sus->C11[2][1] = LF_in[12]; sus->C11[2][2] = LF_in[13];

	sus->m1 = LF_in[14];

	sus->J1[0][0] = LF_in[15]; sus->J1[0][1] = LF_in[16]; sus->J1[0][2] = LF_in[17];
	sus->J1[1][0] = LF_in[18]; sus->J1[1][1] = LF_in[19]; sus->J1[1][2] = LF_in[20];
	sus->J1[2][0] = LF_in[21]; sus->J1[2][1] = LF_in[22]; sus->J1[2][2] = LF_in[23];

	sus->s01p[0] = LF_in[24]; sus->s01p[1] = LF_in[25]; sus->s01p[2] = LF_in[26];

	sus->s12p[0] = LF_in[27]; sus->s12p[1] = LF_in[28]; sus->s12p[2] = LF_in[29];

	sus->C01[0][0] = LF_in[30]; sus->C01[0][1] = LF_in[31]; sus->C01[0][2] = LF_in[32];
	sus->C01[1][0] = LF_in[33]; sus->C01[1][1] = LF_in[34]; sus->C01[1][2] = LF_in[35];
	sus->C01[2][0] = LF_in[36]; sus->C01[2][1] = LF_in[37]; sus->C01[2][2] = LF_in[38];

	sus->C12[0][0] = LF_in[39]; sus->C12[0][1] = LF_in[40]; sus->C12[0][2] = LF_in[41];
	sus->C12[1][0] = LF_in[42]; sus->C12[1][1] = LF_in[43]; sus->C12[1][2] = LF_in[44];
	sus->C12[2][0] = LF_in[45]; sus->C12[2][1] = LF_in[46]; sus->C12[2][2] = LF_in[47];

	sus->s0sp[0] = LF_in[48]; sus->s0sp[1] = LF_in[49]; sus->s0sp[2] = LF_in[50];

	sus->s1sp[0] = LF_in[51]; sus->s1sp[1] = LF_in[52]; sus->s1sp[2] = LF_in[53];

	sus->theta_wh = 0; sus->w_wh = 0; sus->dw_wh = 0;

	sus->rw[0] = chassis->A0[0][0] * sus->s01p[0] + chassis->A0[0][1] * sus->s01p[1] + chassis->A0[0][2] * sus->s01p[2] + chassis->r0[0] + (chassis->A0[0][0] * sus->C01[0][0] + chassis->A0[0][1] * sus->C01[1][0] + chassis->A0[0][2] * sus->C01[2][0]) * sus->s12p[0] + (chassis->A0[0][0] * sus->C01[0][1] + chassis->A0[0][1] * sus->C01[1][1] + chassis->A0[0][2] * sus->C01[2][1]) * sus->s12p[1] + (chassis->A0[0][0] * sus->C01[0][2] + chassis->A0[0][1] * sus->C01[1][2] + chassis->A0[0][2] * sus->C01[2][2]) * sus->s12p[2];
	sus->rw[1] = chassis->A0[1][0] * sus->s01p[0] + chassis->A0[1][1] * sus->s01p[1] + chassis->A0[1][2] * sus->s01p[2] + chassis->r0[1] + (chassis->A0[1][0] * sus->C01[0][0] + chassis->A0[1][1] * sus->C01[1][0] + chassis->A0[1][2] * sus->C01[2][0]) * sus->s12p[0] + (chassis->A0[1][0] * sus->C01[0][1] + chassis->A0[1][1] * sus->C01[1][1] + chassis->A0[1][2] * sus->C01[2][1]) * sus->s12p[1] + (chassis->A0[1][0] * sus->C01[0][2] + chassis->A0[1][1] * sus->C01[1][2] + chassis->A0[1][2] * sus->C01[2][2]) * sus->s12p[2];
	sus->rw[2] = chassis->A0[2][0] * sus->s01p[0] + chassis->A0[2][1] * sus->s01p[1] + chassis->A0[2][2] * sus->s01p[2] + chassis->r0[2] + (chassis->A0[2][0] * sus->C01[0][0] + chassis->A0[2][1] * sus->C01[1][0] + chassis->A0[2][2] * sus->C01[2][0]) * sus->s12p[0] + (chassis->A0[2][0] * sus->C01[0][1] + chassis->A0[2][1] * sus->C01[1][1] + chassis->A0[2][2] * sus->C01[2][1]) * sus->s12p[1] + (chassis->A0[2][0] * sus->C01[0][2] + chassis->A0[2][1] * sus->C01[1][2] + chassis->A0[2][2] * sus->C01[2][2]) * sus->s12p[2];
}

void realtime_dynamics_analysis::read_LM(SUS *sus, CHASSIS *chassis) {
	const double LM_in[54] = { 0, 0, 0.585752158, 0, -0.0932, -0.00000000000375, -0.867194003, 0.497970441, -1, 0.00000000000578, 0.00000000000254, -0.00000000000508, -0.497970441, -0.867194003, 376, 31.23415389, 0, 0, 0, 25.96087203, 0, 0, 0, 20.83953186, 0.83, 0.755, -0.108, 0.648, 0, -0.219, -1, 0.0000000000102, 0, 5.21E-23, 0.0000000000051, -1, -0.0000000000102, -1, -0.0000000000051, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0.683, 0.736, 0.658, 0.281, -0.133, 0.019};

	sus->q1 = LM_in[0];
	sus->dq1 = LM_in[1];

	sus->rho1p[0] = LM_in[2]; sus->rho1p[1] = LM_in[3]; sus->rho1p[2] = LM_in[4];

	sus->C11[0][0] = LM_in[5]; sus->C11[0][1] = LM_in[6]; sus->C11[0][2] = LM_in[7];
	sus->C11[1][0] = LM_in[8]; sus->C11[1][1] = LM_in[9]; sus->C11[1][2] = LM_in[10];
	sus->C11[2][0] = LM_in[11]; sus->C11[2][1] = LM_in[12]; sus->C11[2][2] = LM_in[13];

	sus->m1 = LM_in[14];

	sus->J1[0][0] = LM_in[15]; sus->J1[0][1] = LM_in[16]; sus->J1[0][2] = LM_in[17];
	sus->J1[1][0] = LM_in[18]; sus->J1[1][1] = LM_in[19]; sus->J1[1][2] = LM_in[20];
	sus->J1[2][0] = LM_in[21]; sus->J1[2][1] = LM_in[22]; sus->J1[2][2] = LM_in[23];

	sus->s01p[0] = LM_in[24]; sus->s01p[1] = LM_in[25]; sus->s01p[2] = LM_in[26];

	sus->s12p[0] = LM_in[27]; sus->s12p[1] = LM_in[28]; sus->s12p[2] = LM_in[29];

	sus->C01[0][0] = LM_in[30]; sus->C01[0][1] = LM_in[31]; sus->C01[0][2] = LM_in[32];
	sus->C01[1][0] = LM_in[33]; sus->C01[1][1] = LM_in[34]; sus->C01[1][2] = LM_in[35];
	sus->C01[2][0] = LM_in[36]; sus->C01[2][1] = LM_in[37]; sus->C01[2][2] = LM_in[38];

	sus->C12[0][0] = LM_in[39]; sus->C12[0][1] = LM_in[40]; sus->C12[0][2] = LM_in[41];
	sus->C12[1][0] = LM_in[42]; sus->C12[1][1] = LM_in[43]; sus->C12[1][2] = LM_in[44];
	sus->C12[2][0] = LM_in[45]; sus->C12[2][1] = LM_in[46]; sus->C12[2][2] = LM_in[47];

	sus->s0sp[0] = LM_in[48]; sus->s0sp[1] = LM_in[49]; sus->s0sp[2] = LM_in[50];

	sus->s1sp[0] = LM_in[51]; sus->s1sp[1] = LM_in[52]; sus->s1sp[2] = LM_in[53];

	sus->theta_wh = 0; sus->w_wh = 0; sus->dw_wh = 0;

	sus->rw[0] = chassis->A0[0][0] * sus->s01p[0] + chassis->A0[0][1] * sus->s01p[1] + chassis->A0[0][2] * sus->s01p[2] + chassis->r0[0] + (chassis->A0[0][0] * sus->C01[0][0] + chassis->A0[0][1] * sus->C01[1][0] + chassis->A0[0][2] * sus->C01[2][0]) * sus->s12p[0] + (chassis->A0[0][0] * sus->C01[0][1] + chassis->A0[0][1] * sus->C01[1][1] + chassis->A0[0][2] * sus->C01[2][1]) * sus->s12p[1] + (chassis->A0[0][0] * sus->C01[0][2] + chassis->A0[0][1] * sus->C01[1][2] + chassis->A0[0][2] * sus->C01[2][2]) * sus->s12p[2];
	sus->rw[1] = chassis->A0[1][0] * sus->s01p[0] + chassis->A0[1][1] * sus->s01p[1] + chassis->A0[1][2] * sus->s01p[2] + chassis->r0[1] + (chassis->A0[1][0] * sus->C01[0][0] + chassis->A0[1][1] * sus->C01[1][0] + chassis->A0[1][2] * sus->C01[2][0]) * sus->s12p[0] + (chassis->A0[1][0] * sus->C01[0][1] + chassis->A0[1][1] * sus->C01[1][1] + chassis->A0[1][2] * sus->C01[2][1]) * sus->s12p[1] + (chassis->A0[1][0] * sus->C01[0][2] + chassis->A0[1][1] * sus->C01[1][2] + chassis->A0[1][2] * sus->C01[2][2]) * sus->s12p[2];
	sus->rw[2] = chassis->A0[2][0] * sus->s01p[0] + chassis->A0[2][1] * sus->s01p[1] + chassis->A0[2][2] * sus->s01p[2] + chassis->r0[2] + (chassis->A0[2][0] * sus->C01[0][0] + chassis->A0[2][1] * sus->C01[1][0] + chassis->A0[2][2] * sus->C01[2][0]) * sus->s12p[0] + (chassis->A0[2][0] * sus->C01[0][1] + chassis->A0[2][1] * sus->C01[1][1] + chassis->A0[2][2] * sus->C01[2][1]) * sus->s12p[1] + (chassis->A0[2][0] * sus->C01[0][2] + chassis->A0[2][1] * sus->C01[1][2] + chassis->A0[2][2] * sus->C01[2][2]) * sus->s12p[2];
}

void realtime_dynamics_analysis::read_LR(SUS *sus, CHASSIS *chassis) {
	const double LR_in[54] = { 0, 0, 0.585382979, 0, -0.0932, -0.00000000000376, -0.868782601, 0.495193692, -1, 0.00000000000577, 0.00000000000253, -0.00000000000505, -0.495193692, -0.868782601, 376, 31.21662106, 0, 0, 0, 25.94143616, 0, 0, 0, 20.8414349, -0.77, 0.755, -0.108, 0.648, 0, -0.219, -1, 0.0000000000102, 0, 5.21E-23, 0.0000000000051, -1, -0.0000000000102, -1, -0.0000000000051, 1, 0, 0, 0, 1, 0, 0, 0, 1, -0.917, 0.736, 0.658, 0.281, -0.133, 0.019};

	sus->q1 = LR_in[0];
	sus->dq1 = LR_in[1];

	sus->rho1p[0] = LR_in[2]; sus->rho1p[1] = LR_in[3]; sus->rho1p[2] = LR_in[4];

	sus->C11[0][0] = LR_in[5]; sus->C11[0][1] = LR_in[6]; sus->C11[0][2] = LR_in[7];
	sus->C11[1][0] = LR_in[8]; sus->C11[1][1] = LR_in[9]; sus->C11[1][2] = LR_in[10];
	sus->C11[2][0] = LR_in[11]; sus->C11[2][1] = LR_in[12]; sus->C11[2][2] = LR_in[13];

	sus->m1 = LR_in[14];

	sus->J1[0][0] = LR_in[15]; sus->J1[0][1] = LR_in[16]; sus->J1[0][2] = LR_in[17];
	sus->J1[1][0] = LR_in[18]; sus->J1[1][1] = LR_in[19]; sus->J1[1][2] = LR_in[20];
	sus->J1[2][0] = LR_in[21]; sus->J1[2][1] = LR_in[22]; sus->J1[2][2] = LR_in[23];

	sus->s01p[0] = LR_in[24]; sus->s01p[1] = LR_in[25]; sus->s01p[2] = LR_in[26];

	sus->s12p[0] = LR_in[27]; sus->s12p[1] = LR_in[28]; sus->s12p[2] = LR_in[29];

	sus->C01[0][0] = LR_in[30]; sus->C01[0][1] = LR_in[31]; sus->C01[0][2] = LR_in[32];
	sus->C01[1][0] = LR_in[33]; sus->C01[1][1] = LR_in[34]; sus->C01[1][2] = LR_in[35];
	sus->C01[2][0] = LR_in[36]; sus->C01[2][1] = LR_in[37]; sus->C01[2][2] = LR_in[38];

	sus->C12[0][0] = LR_in[39]; sus->C12[0][1] = LR_in[40]; sus->C12[0][2] = LR_in[41];
	sus->C12[1][0] = LR_in[42]; sus->C12[1][1] = LR_in[43]; sus->C12[1][2] = LR_in[44];
	sus->C12[2][0] = LR_in[45]; sus->C12[2][1] = LR_in[46]; sus->C12[2][2] = LR_in[47];

	sus->s0sp[0] = LR_in[48]; sus->s0sp[1] = LR_in[49]; sus->s0sp[2] = LR_in[50];

	sus->s1sp[0] = LR_in[51]; sus->s1sp[1] = LR_in[52]; sus->s1sp[2] = LR_in[53];

	sus->theta_wh = 0; sus->w_wh = 0; sus->dw_wh = 0;

	sus->rw[0] = chassis->A0[0][0] * sus->s01p[0] + chassis->A0[0][1] * sus->s01p[1] + chassis->A0[0][2] * sus->s01p[2] + chassis->r0[0] + (chassis->A0[0][0] * sus->C01[0][0] + chassis->A0[0][1] * sus->C01[1][0] + chassis->A0[0][2] * sus->C01[2][0]) * sus->s12p[0] + (chassis->A0[0][0] * sus->C01[0][1] + chassis->A0[0][1] * sus->C01[1][1] + chassis->A0[0][2] * sus->C01[2][1]) * sus->s12p[1] + (chassis->A0[0][0] * sus->C01[0][2] + chassis->A0[0][1] * sus->C01[1][2] + chassis->A0[0][2] * sus->C01[2][2]) * sus->s12p[2];
	sus->rw[1] = chassis->A0[1][0] * sus->s01p[0] + chassis->A0[1][1] * sus->s01p[1] + chassis->A0[1][2] * sus->s01p[2] + chassis->r0[1] + (chassis->A0[1][0] * sus->C01[0][0] + chassis->A0[1][1] * sus->C01[1][0] + chassis->A0[1][2] * sus->C01[2][0]) * sus->s12p[0] + (chassis->A0[1][0] * sus->C01[0][1] + chassis->A0[1][1] * sus->C01[1][1] + chassis->A0[1][2] * sus->C01[2][1]) * sus->s12p[1] + (chassis->A0[1][0] * sus->C01[0][2] + chassis->A0[1][1] * sus->C01[1][2] + chassis->A0[1][2] * sus->C01[2][2]) * sus->s12p[2];
	sus->rw[2] = chassis->A0[2][0] * sus->s01p[0] + chassis->A0[2][1] * sus->s01p[1] + chassis->A0[2][2] * sus->s01p[2] + chassis->r0[2] + (chassis->A0[2][0] * sus->C01[0][0] + chassis->A0[2][1] * sus->C01[1][0] + chassis->A0[2][2] * sus->C01[2][0]) * sus->s12p[0] + (chassis->A0[2][0] * sus->C01[0][1] + chassis->A0[2][1] * sus->C01[1][1] + chassis->A0[2][2] * sus->C01[2][1]) * sus->s12p[1] + (chassis->A0[2][0] * sus->C01[0][2] + chassis->A0[2][1] * sus->C01[1][2] + chassis->A0[2][2] * sus->C01[2][2]) * sus->s12p[2];
}

void realtime_dynamics_analysis::read_control(CTRL *ctrl) {
	// 최초 1회만 실행되어 각 변수들을 초기화 시켜줌
	ctrl->e_v_sum = 0;
	ctrl->M_d_hat = 0;
	ctrl->yaw_hat = 0;
	ctrl->WP_size = WP_size;
}

void realtime_dynamics_analysis::define_Y_vector(SIM *sim, CHASSIS *chassis, SUS sus[6]) {
	memcpy(sim->Y, chassis->r0, sizeof(double) * 3);
	memcpy(sim->Y + 3, chassis->e0, sizeof(double) * 4);

	sim->Y[7] = sus[RF].q1;
	sim->Y[8] = sus[RM].q1;
	sim->Y[9] = sus[RR].q1;
	sim->Y[10] = sus[LF].q1;
	sim->Y[11] = sus[LM].q1;
	sim->Y[12] = sus[LR].q1;

	memcpy(sim->Y + 13, chassis->dr0, sizeof(double) * 3);
	memcpy(sim->Y + 16, chassis->w0, sizeof(double) * 3);

	sim->Y[19] = sus[RF].dq1;
	sim->Y[20] = sus[RM].dq1;
	sim->Y[21] = sus[RR].dq1;
	sim->Y[22] = sus[LF].dq1;
	sim->Y[23] = sus[LM].dq1;
	sim->Y[24] = sus[LR].dq1;

	sim->Y[25] = sus[RF].w_wh;
	sim->Y[26] = sus[RM].w_wh;
	sim->Y[27] = sus[RR].w_wh;
	sim->Y[28] = sus[LF].w_wh;
	sim->Y[29] = sus[LM].w_wh;
	sim->Y[30] = sus[LR].w_wh;
}

// LU dcomposition
int realtime_dynamics_analysis::ludcmp6(double a[6][6], int n, int indx[6], double d, double fac[6][6])
{
	int i, imax, j, k;
	double big, dum, sum, temp;
	//double *vv = (double*)malloc(sizeof(double)*n);
	double *vv = new double[n];

	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			fac[i][j] = a[i][j];

	d = 1.0;
	for (i = 0; i < n; i++) {
		big = 0.0;
		for (j = 0; j < n; j++)
			if ((temp = fabs(fac[i][j])) > big) big = temp;
		if (big == 0.0) {
			printf("Singular matrix in routine LUdcmp - dynamics");
			return 1;
		}
		vv[i] = 1.0 / big;
	}
	for (j = 0; j < n; j++) {
		for (i = 0; i < j; i++) {
			sum = fac[i][j];
			for (k = 0; k < i; k++) sum -= fac[i][k] * fac[k][j];
			fac[i][j] = sum;
		}
		big = 0.0;
		for (i = j; i < n; i++) {
			sum = fac[i][j];
			for (k = 0; k < j; k++) sum -= fac[i][k] * fac[k][j];
			fac[i][j] = sum;
			if ((dum = vv[i] * fabs(sum)) >= big) {
				big = dum;
				imax = i;
			}
		}
		if (j != imax) {
			for (k = 0; k < n; k++) {
				dum = fac[imax][k];
				fac[imax][k] = fac[j][k];
				fac[j][k] = dum;
			}
			d = -d;
			vv[imax] = vv[j];
		}
		indx[j] = imax;
		if (fac[j][j] == 0.0) fac[j][j] = TINY;
		if (j != n - 1) {
			dum = 1.0 / (fac[j][j]);
			for (i = j + 1; i < n; i++) fac[i][j] *= dum;
		}
	}

	//free(vv);
	delete[] vv;

	return 0;
}

// Back substitution
void realtime_dynamics_analysis::lubksb6(double a[6][6], int n, int indx[6], double b[6], double x[6])
{
	int i, ii = 0, ip, j;
	double sum;

	for (i = 0; i < n; i++) x[i] = b[i];
	for (i = 0; i < n; i++) {
		ip = indx[i];
		sum = x[ip];
		x[ip] = x[i];
		if (ii)
			for (j = ii - 1; j < i; j++) sum -= a[i][j] * x[j];
		else if (sum) ii = i + 1;
		x[i] = sum;
	}
	for (i = n - 1; i >= 0; i--) {
		sum = x[i];
		for (j = i + 1; j < n; j++) sum -= a[i][j] * x[j];
		x[i] = sum / a[i][i];
	}
}

/***********************************************************************
subroutine skew_symmetric(a,b)
B(3x3)=A(3x1)(~)
***********************************************************************/
void realtime_dynamics_analysis::tilde(double a[3], double b[3][3])
{
	b[0][0] = 0;
	b[0][1] = -a[2];
	b[0][2] = a[1];
	b[1][0] = a[2];
	b[1][1] = 0;
	b[1][2] = -a[0];
	b[2][0] = -a[1];
	b[2][1] = a[0];
	b[2][2] = 0;
}

/***********************************************************************
subroutine mat3333(a,b,c)
C(3x3)=A(3x3)*B(3x3)
***********************************************************************/
void realtime_dynamics_analysis::mat3333(double a[3][3], double b[3][3], double c[3][3])
{
	c[0][0] = a[0][0] * b[0][0] + a[0][1] * b[1][0] + a[0][2] * b[2][0];
	c[0][1] = a[0][0] * b[0][1] + a[0][1] * b[1][1] + a[0][2] * b[2][1];
	c[0][2] = a[0][0] * b[0][2] + a[0][1] * b[1][2] + a[0][2] * b[2][2];
	c[1][0] = a[1][0] * b[0][0] + a[1][1] * b[1][0] + a[1][2] * b[2][0];
	c[1][1] = a[1][0] * b[0][1] + a[1][1] * b[1][1] + a[1][2] * b[2][1];
	c[1][2] = a[1][0] * b[0][2] + a[1][1] * b[1][2] + a[1][2] * b[2][2];
	c[2][0] = a[2][0] * b[0][0] + a[2][1] * b[1][0] + a[2][2] * b[2][0];
	c[2][1] = a[2][0] * b[0][1] + a[2][1] * b[1][1] + a[2][2] * b[2][1];
	c[2][2] = a[2][0] * b[0][2] + a[2][1] * b[1][2] + a[2][2] * b[2][2];
}

/***********************************************************************
subroutine mat3331(a,b,c)
c(3x1)=A(3x3)*B(3x1)
***********************************************************************/
void realtime_dynamics_analysis::mat3331(double a[3][3], double b[3], double c[3])
{
	c[0] = a[0][0] * b[0] + a[0][1] * b[1] + a[0][2] * b[2];
	c[1] = a[1][0] * b[0] + a[1][1] * b[1] + a[1][2] * b[2];
	c[2] = a[2][0] * b[0] + a[2][1] * b[1] + a[2][2] * b[2];
}

/***********************************************************************
subroutine mat333333(a,b,c,d)
D(3x3)=A(3x3)*B(3x3)*C(3x3)
***********************************************************************/
void realtime_dynamics_analysis::mat333333(double a[3][3], double b[3][3], double c[3][3], double d[3][3])
{
	d[0][0] = (a[0][0] * b[0][0] + a[0][1] * b[1][0] + a[0][2] * b[2][0])*c[0][0] + (a[0][0] * b[0][1] + a[0][1] * b[1][1] + a[0][2] * b[2][1])*c[1][0] + (a[0][0] * b[0][2] + a[0][1] * b[1][2] + a[0][2] * b[2][2])*c[2][0];
	d[0][1] = (a[0][0] * b[0][0] + a[0][1] * b[1][0] + a[0][2] * b[2][0])*c[0][1] + (a[0][0] * b[0][1] + a[0][1] * b[1][1] + a[0][2] * b[2][1])*c[1][1] + (a[0][0] * b[0][2] + a[0][1] * b[1][2] + a[0][2] * b[2][2])*c[2][1];
	d[0][2] = (a[0][0] * b[0][0] + a[0][1] * b[1][0] + a[0][2] * b[2][0])*c[0][2] + (a[0][0] * b[0][1] + a[0][1] * b[1][1] + a[0][2] * b[2][1])*c[1][2] + (a[0][0] * b[0][2] + a[0][1] * b[1][2] + a[0][2] * b[2][2])*c[2][2];
	d[1][0] = (a[1][0] * b[0][0] + a[1][1] * b[1][0] + a[1][2] * b[2][0])*c[0][0] + (a[1][0] * b[0][1] + a[1][1] * b[1][1] + a[1][2] * b[2][1])*c[1][0] + (a[1][0] * b[0][2] + a[1][1] * b[1][2] + a[1][2] * b[2][2])*c[2][0];
	d[1][1] = (a[1][0] * b[0][0] + a[1][1] * b[1][0] + a[1][2] * b[2][0])*c[0][1] + (a[1][0] * b[0][1] + a[1][1] * b[1][1] + a[1][2] * b[2][1])*c[1][1] + (a[1][0] * b[0][2] + a[1][1] * b[1][2] + a[1][2] * b[2][2])*c[2][1];
	d[1][2] = (a[1][0] * b[0][0] + a[1][1] * b[1][0] + a[1][2] * b[2][0])*c[0][2] + (a[1][0] * b[0][1] + a[1][1] * b[1][1] + a[1][2] * b[2][1])*c[1][2] + (a[1][0] * b[0][2] + a[1][1] * b[1][2] + a[1][2] * b[2][2])*c[2][2];
	d[2][0] = (a[2][0] * b[0][0] + a[2][1] * b[1][0] + a[2][2] * b[2][0])*c[0][0] + (a[2][0] * b[0][1] + a[2][1] * b[1][1] + a[2][2] * b[2][1])*c[1][0] + (a[2][0] * b[0][2] + a[2][1] * b[1][2] + a[2][2] * b[2][2])*c[2][0];
	d[2][1] = (a[2][0] * b[0][0] + a[2][1] * b[1][0] + a[2][2] * b[2][0])*c[0][1] + (a[2][0] * b[0][1] + a[2][1] * b[1][1] + a[2][2] * b[2][1])*c[1][1] + (a[2][0] * b[0][2] + a[2][1] * b[1][2] + a[2][2] * b[2][2])*c[2][1];
	d[2][2] = (a[2][0] * b[0][0] + a[2][1] * b[1][0] + a[2][2] * b[2][0])*c[0][2] + (a[2][0] * b[0][1] + a[2][1] * b[1][1] + a[2][2] * b[2][1])*c[1][2] + (a[2][0] * b[0][2] + a[2][1] * b[1][2] + a[2][2] * b[2][2])*c[2][2];
}

/***********************************************************************
subroutine mat33T31(a,b,c)
C(3x1)=A(3x3)'*B(3x1)
***********************************************************************/
void realtime_dynamics_analysis::mat33T31(double a[3][3], double b[3], double c[3])
{
	c[0] = a[0][0] * b[0] + a[1][0] * b[1] + a[2][0] * b[2];
	c[1] = a[0][1] * b[0] + a[1][1] * b[1] + a[2][1] * b[2];
	c[2] = a[0][2] * b[0] + a[1][2] * b[1] + a[2][2] * b[2];
}

/***********************************************************************
subroutine mat34T31(a,b,c)
C(4x1)=A(3x4)'*B(3x1)
***********************************************************************/
void realtime_dynamics_analysis::mat34T31(double a[3][4], double b[3], double c[4])
{
	c[0] = a[0][0] * b[0] + a[1][0] * b[1] + a[2][0] * b[2];
	c[1] = a[0][1] * b[0] + a[1][1] * b[1] + a[2][1] * b[2];
	c[2] = a[0][2] * b[0] + a[1][2] * b[1] + a[2][2] * b[2];
	c[3] = a[0][3] * b[0] + a[1][3] * b[1] + a[2][3] * b[2];
}

/***********************************************************************
subroutine mat3333T(a,b,c)
C(3x3)=A(3x3)*B(3x3)'
***********************************************************************/
void realtime_dynamics_analysis::mat3333T(double a[3][3], double b[3][3], double c[3][3])
{
	c[0][0] = a[0][0] * b[0][0] + a[0][1] * b[0][1] + a[0][2] * b[0][2];
	c[0][1] = a[0][0] * b[1][0] + a[0][1] * b[1][1] + a[0][2] * b[1][2];
	c[0][2] = a[0][0] * b[2][0] + a[0][1] * b[2][1] + a[0][2] * b[2][2];
	c[1][0] = a[1][0] * b[0][0] + a[1][1] * b[0][1] + a[1][2] * b[0][2];
	c[1][1] = a[1][0] * b[1][0] + a[1][1] * b[1][1] + a[1][2] * b[1][2];
	c[1][2] = a[1][0] * b[2][0] + a[1][1] * b[2][1] + a[1][2] * b[2][2];
	c[2][0] = a[2][0] * b[0][0] + a[2][1] * b[0][1] + a[2][2] * b[0][2];
	c[2][1] = a[2][0] * b[1][0] + a[2][1] * b[1][1] + a[2][2] * b[1][2];
	c[2][2] = a[2][0] * b[2][0] + a[2][1] * b[2][1] + a[2][2] * b[2][2];
}

/***********************************************************************
subroutine mat333331(a,b,c,d)
D(3x1)=A(3x3)*B(3x3)*C(3x1)
***********************************************************************/
void realtime_dynamics_analysis::mat333331(double a[3][3], double b[3][3], double c[3], double d[3])
{
	d[0] = (a[0][0] * b[0][0] + a[0][1] * b[1][0] + a[0][2] * b[2][0]) * c[0] + (a[0][0] * b[0][1] + a[0][1] * b[1][1] + a[0][2] * b[2][1]) * c[1] + (a[0][0] * b[0][2] + a[0][1] * b[1][2] + a[0][2] * b[2][2]) * c[2];
	d[1] = (a[1][0] * b[0][0] + a[1][1] * b[1][0] + a[1][2] * b[2][0]) * c[0] + (a[1][0] * b[0][1] + a[1][1] * b[1][1] + a[1][2] * b[2][1]) * c[1] + (a[1][0] * b[0][2] + a[1][1] * b[1][2] + a[1][2] * b[2][2]) * c[2];
	d[2] = (a[2][0] * b[0][0] + a[2][1] * b[1][0] + a[2][2] * b[2][0]) * c[0] + (a[2][0] * b[0][1] + a[2][1] * b[1][1] + a[2][2] * b[2][1]) * c[1] + (a[2][0] * b[0][2] + a[2][1] * b[1][2] + a[2][2] * b[2][2]) * c[2];
}

/***********************************************************************
subroutine mat6661(a,b,c)
C(6x1)=A(6x6)*B(6x1)
***********************************************************************/
void realtime_dynamics_analysis::mat6661(double a[6][6], double b[6], double c[6])
{
	c[0] = a[0][0] * b[0] + a[0][1] * b[1] + a[0][2] * b[2] + a[0][3] * b[3] + a[0][4] * b[4] + a[0][5] * b[5];
	c[1] = a[1][0] * b[0] + a[1][1] * b[1] + a[1][2] * b[2] + a[1][3] * b[3] + a[1][4] * b[4] + a[1][5] * b[5];
	c[2] = a[2][0] * b[0] + a[2][1] * b[1] + a[2][2] * b[2] + a[2][3] * b[3] + a[2][4] * b[4] + a[2][5] * b[5];
	c[3] = a[3][0] * b[0] + a[3][1] * b[1] + a[3][2] * b[2] + a[3][3] * b[3] + a[3][4] * b[4] + a[3][5] * b[5];
	c[4] = a[4][0] * b[0] + a[4][1] * b[1] + a[4][2] * b[2] + a[4][3] * b[3] + a[4][4] * b[4] + a[4][5] * b[5];
	c[5] = a[5][0] * b[0] + a[5][1] * b[1] + a[5][2] * b[2] + a[5][3] * b[3] + a[5][4] * b[4] + a[5][5] * b[5];
}

double realtime_dynamics_analysis::Find_Max(double Array[], int length)
{
	double MAX = fabs(Array[0]);
	int i;
	for (i = 1; i<length; i++)
	{
		if (fabs(Array[i]) > MAX)
		{
			MAX = fabs(Array[i]);
		}
	}
	return MAX;
}

double realtime_dynamics_analysis::Find_min(double data, double min_val)
{
	if (data < min_val)
		min_val = data;

	return min_val;
}

double realtime_dynamics_analysis::fsign(double data)
{
	if (data > 0)
		return 1.0;
	else if (data < 0)
		return -1.0;
	else
		return 0;
}
int realtime_dynamics_analysis::ludcmp4(double a[4][4], int n, int indx[4], double d, double fac[4][4])
{
	int i, imax, j, k;
	double big, temp;
	double vv[4];
	for (i = 0; i < n; i++) {
		big = 0.0;
		for (j = 0; j < n; j++)
			if ((temp = fabs(a[i][j])) > big) big = temp;
		if (big == 0.0) {
			printf("Singular matrix in LUdcmp - controller");
			return 1;
		}
		vv[i] = 1.0 / big;
	}
	for (k = 0; k < n; k++) {
		big = 0.0;
		for (i = k; i<n; i++) {
			temp = vv[i] * fabs(a[i][k]);
			if (temp > big) {
				big = temp;
				imax = i;
			}
		}
		if (k != imax) {
			for (j = 0; j < n; j++) {
				temp = a[imax][j];
				a[imax][j] = a[k][j];
				a[k][j] = temp;
			}
			d = -d;
			vv[imax] = vv[k];
		}
		indx[k] = imax;
		if (a[k][k] == 0.0) a[k][k] = TINY;
		for (i = k + 1; i < n; i++) {
			temp = a[i][k] /= a[k][k];
			for (j = k + 1; j < n; j++)
				a[i][j] -= temp*a[k][j];
		}
	}
	//////////////////
	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
			fac[i][j] = a[i][j];
	//////////////////

	return 0;
}

void realtime_dynamics_analysis::lusolve4(double a[4][4], int n, int indx[4], double b[4], double x[4])
{
	int i, ii = 0, ip, j;
	double sum;
	for (i = 0; i < n; i++) x[i] = b[i];
	for (i = 0; i < n; i++) {
		ip = indx[i];
		sum = x[ip];
		x[ip] = x[i];
		if (ii != 0)
			for (j = ii - 1; j < i; j++) sum -= a[i][j] * x[j];
		else if (sum != 0.0)
			ii = i + 1;
		x[i] = sum;
	}
	for (i = n - 1; i >= 0; i--) {
		sum = x[i];
		for (j = i + 1; j < n; j++) sum -= a[i][j] * x[j];
		x[i] = sum / a[i][i];
	}
}