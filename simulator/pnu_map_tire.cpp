#include "pnu_map_tire.h"

#pragma warning(disable:4996)

using namespace std;

uint	Obj_Map_Global::_num_x; //포인터로 수정
uint	Obj_Map_Global::_num_y;
double	Obj_Map_Global::_gab;
double	Obj_Map_Global::_x_c;
double	Obj_Map_Global::_y_c;
double	*Obj_Map_Global::_x_data = NULL; //갯수가 런타임 결정이기 때문에 배열 사용곤란
double	*Obj_Map_Global::_y_data = NULL;
const double*const*	Obj_Map_Global::_map_data = NULL; //이중포인터 지시내용의 변경 방지

Obj_Map_Global::~Obj_Map_Global()
{
	delete[] _x_data;
	_x_data = NULL;

	delete[] _y_data;
	_y_data = NULL;
}

void Obj_Map_Global::PNU_set_map(double **map_data, uint rows, uint cols, double gab, double x_c, double y_c)
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

void Obj_Map_Global::PNU_fn_map(double x_p, double y_p, double *z)
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

void Obj_Map_Global::fn_linspace(double *data, double L, uint num_data, double c)
{
	double gab = L / (double)(num_data - 1);
	uint i;
	for (i = 0; i < num_data; ++i) {
		data[i] = gab*(double)(i)-c;
	}
}

void Obj_Map_Global::fn_find_patch_uv(uint* patch_r, uint* patch_c, double* u, double* v, double x_p, double y_p)
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

void Obj_Map_Global::fn_MBMt_uv(double(*MBMt)[4], double* u, double* v, double x_p, double y_p)
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

void Obj_Map_Global::fn_patch_idxs(int row_idx[4], int col_idx[4], uint patch_r, uint patch_c)
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

void Obj_Map_Global::fn_B_mat(double(*B)[4], uint patch_r, uint patch_c)
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

void Obj_Map_Global::fn_mat444_multi(double(*C)[4], double(*A)[4], double(*B)[4])
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

void Obj_Map_Global::fn_mat1441_multi(double *A, double(*B)[4], double *C, double *output)
{
	double sol = C[0] * (A[0] * B[0][0] + A[1] * B[1][0] + A[2] * B[2][0] + A[3] * B[3][0]) + C[1] * (A[0] * B[0][1] + A[1] * B[1][1] + A[2] * B[2][1] + A[3] * B[3][1]) + C[2] * (A[0] * B[0][2] + A[1] * B[1][2] + A[2] * B[2][2] + A[3] * B[3][2]) + C[3] * (A[0] * B[0][3] + A[1] * B[1][3] + A[2] * B[2][3] + A[3] * B[3][3]);
	*output = sol;
}

double Obj_Nan_Global::z_min = 100.0; // nan 판정기준 최소 높이

void Obj_Nan_Global::fn_nan_interp(double **Map_data, uint n_rows, uint n_cols)
{
	int i = 0, j = 0;
	for (i = 0; i < n_rows; ++i) {
		for (j = 0; j < n_cols; ++j) {
			fn_nan_patch(Map_data, n_rows, n_cols, i, j);
		}
	}
}


double Obj_Nan_Global::fn_nan_patch(double **Map_data, uint n_rows, uint n_cols, uint i, uint j)
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

int Obj_Nan_Global::is_nan(double z)
{
	if (z < z_min) {
		return 1;
	}
	else {
		return 0;
	}
}

void Obj_Nan_Global::fn_ij_p(int* i_p, int* j_p, uint p, uint k, uint i, uint j)
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

Obj_Tire_Global::Obj_Tire_Global()
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

void Obj_Tire_Global::PNU_Pre_road(double R_c_RF[3], double R_c_RM[3], double R_c_RR[3], double R_c_LF[3], double R_c_LM[3], double R_c_LR[3]) {
	Obj_Map_Global::PNU_fn_map(R_c_RF[0], R_c_RF[1], &_Road_h_i[0]);
	Obj_Map_Global::PNU_fn_map(R_c_RM[0], R_c_RM[1], &_Road_h_i[1]);
	Obj_Map_Global::PNU_fn_map(R_c_RR[0], R_c_RR[1], &_Road_h_i[2]);
	Obj_Map_Global::PNU_fn_map(R_c_LF[0], R_c_LF[1], &_Road_h_i[3]);
	Obj_Map_Global::PNU_fn_map(R_c_LM[0], R_c_LM[1], &_Road_h_i[4]);
	Obj_Map_Global::PNU_fn_map(R_c_LR[0], R_c_LR[1], &_Road_h_i[5]);

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

void Obj_Tire_Global::PNU_Tire_force(double R_c[3], double dR_c[3], double u_vec[3], double omega, int id)
{
	// 수직력 계산
	double R_x = R_c[0];
	double R_y = R_c[1];
	//Obj_Map::PNU_fn_map(R_x, R_y, &_road_h);
	//road_h = -0.63; // PNU_fn_map(R_x,R_y);	//평지로 가정
	
	//double _N_vec[3] = { 0,0,1 };	// 법선 벡터 해제
	double R_z = R_c[2];
	double dR_z = dR_c[2];
	double point_h;
	point_h = R_z - _R_u;	//tire의 point follow점
	fn_max(_Road_h_i[id] - point_h, 0.0, &_pen);	//tire penetration 침투량
	//fn_max(_road_h - point_h, 0.0, &_pen);		//tire penetration 침투량
	_road_h = _Road_h_i[id];

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

	//double N_vec[3] = { 0, 0, 1 };
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

	double mu = 0.8;
	fn_Brush_tire(F_road, M_road, Fz, _slip_x, _slip_y, mu, a, _pen, omega);

	// 타이어 힘의 road->tire 좌표변환
	//vec_p2c[0] = 0; //tire의 point center에서 point follow점까지의 vector, 좌표계 변환때 사용
	//vec_p2c[1] = 0;
	//vec_p2c[2] = -_R_d;

	//vec_p2c[0] = -_N_vec[0]; //tire의 point center에서 point follow점까지의 vector, 좌표계 변환때 사용
	//vec_p2c[1] = -_N_vec[1];
	//vec_p2c[2] = -_N_vec[2] * _R_d;

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

void Obj_Tire_Global::PNU_get_data(double F[3], double M[3], double *slip, double *angle, double *road_h, double *pen, double *R_d, double *Fx, double *Fy, double *Fz, double *My)
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

void Obj_Tire_Global::fn_Brush_tire(double F[3], double M[3], double Fz, double slip_x, double slip_y, double mu, double a, double pen, double omega)
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

void Obj_Tire_Global::fn_R_eff(double *R_e, double *a, double R_d) {

	double e = acos(R_d / _R_u);
	*a = _R_u * sin(e);
	*R_e = *a / (e + _eps);
}

void Obj_Tire_Global::fn_road_tire_A(double road_A[3][3], double tire_A[3][3], double road_z_Vec[3], double tire_y_Vec[3]) {

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

void Obj_Tire_Global::fn_slip_ratio(double *slip_x, double *slip_y, double V_x, double omega, double R_e, double V_sy, double Fz, double a) {
	double speed_m_x, speed_x, speed_m_y, speed_y, temp1, temp2;

	_CSLIP = 2 * a*a * _c_x; // N / rad : cornering stiffness
	_Ca = 2 * a*a * _c_y;    // N / m : longitudinal stiffness

	speed_m_x = _CSLIP * R_e * R_e * _t_s / (2 * _Iyy);  // 원래 슬립식에서 Euler method 사용시 unstable 조건의 속도
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

void Obj_Tire_Global::fn_Fiala_tire(double F[3], double M[3], double Fz, double slip, double alpha, double omega)
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

void Obj_Tire_Global::fn_step_s(double x, double slope, double *output) {
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

void Obj_Tire_Global::fn_cross(double c[3], double a[3], double b[3])
{
	c[0] = a[1] * b[2] - a[2] * b[1];
	c[1] = a[2] * b[0] - a[0] * b[2];
	c[2] = a[0] * b[1] - a[1] * b[0];
}

void Obj_Tire_Global::fn_mat3331(double c[3], double a[3][3], double b[3])
{
	//c(3x1)=A(3x3)*B(3x1)
	c[0] = a[0][0] * b[0] + a[0][1] * b[1] + a[0][2] * b[2];
	c[1] = a[1][0] * b[0] + a[1][1] * b[1] + a[1][2] * b[2];
	c[2] = a[2][0] * b[0] + a[2][1] * b[1] + a[2][2] * b[2];
}

void Obj_Tire_Global::fn_mat33T31(double c[3], double a[3][3], double b[3])
{
	//c(3x1)=A(3x3)'*B(3x1)
	c[0] = a[0][0] * b[0] + a[1][0] * b[1] + a[2][0] * b[2];
	c[1] = a[0][1] * b[0] + a[1][1] * b[1] + a[2][1] * b[2];
	c[2] = a[0][2] * b[0] + a[1][2] * b[1] + a[2][2] * b[2];
}

void Obj_Tire_Global::fn_sign(double a, double *output)
{
	double b;
	if (a > 0) { b = 1; }
	else if (a < 0) { b = -1; }
	else { b = 0; }

	*output = b;
}

void Obj_Tire_Global::fn_max(double a, double b, double *output)
{
	double c;
	if (a >= b) { c = a; }
	else { c = b; }

	*output = c;
}

void Obj_Tire_Global::fn_min(double a, double b, double *output)
{
	double c;
	if (a <= b) { c = a; }
	else { c = b; }

	*output = c;
}

void Obj_Tire_Global::fn_normalize(double a[], uint num)
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

