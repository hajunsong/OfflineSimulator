#include <stdio.h>
#include <direct.h>
#include "unmanned_ground_vehicle.h"

void GaussianBlur77(double** ElevationData, int row_min, int row_max, int col_min, int col_max, int row, int col, int iter=1) {
	const int box[7 * 7] = {
		1,  4,  7,  10, 7,  4,  1,
		4,  16, 26, 33, 26, 16, 4,
		7,  26, 55, 71, 55, 26, 7,
		10, 33, 71, 91, 71, 33, 10,
		7,  26, 55, 71, 55, 26, 7,
		4,  16, 26, 33, 26, 16, 4,
		1,  4,  7,  10, 7,  4,  1
	};

	int start_y = row_min;
	int start_x = col_min;
	int end_y = row_max;
	int end_x = col_max;
	double sum = 0;
	int height = row_max - row_min;
	int width = col_max - col_min;
	int ksize = 7;
	int offset = ksize / 2;
	double scale = 0;
	for (int i = 0; i < ksize*ksize; i++) {
		scale += box[i];
	}

	double **ElevationDataFiltered;
	ElevationDataFiltered = new double*[row];
	for (int i = 0; i < row; i++) {
		ElevationDataFiltered[i] = new double[col];
	}

	for (int count = 0; count < iter; count++) {
		for (int y = start_y; y < end_y; y++) {
			for (int x = start_x; x < end_x; x++) {
				sum = 0;
				for (int i = -offset, ik = 0; i <= offset; i++, ik++) {
					int ii = y + i < 0 ? 0 : y + i >= row ? row - 1 : y + i;
					for (int j = -offset, jk = 0; j <= offset; j++, jk++) {
						int jj = x + j < 0 ? 0 : x + j >= col ? col - 1 : x + j;
						sum += ElevationData[ii][jj] * box[ik*ksize + jk];
					}
				}
				ElevationDataFiltered[y][x] = sum / scale;
			}
		}

		for (int y = start_y; y < end_y; y++) {
			for (int x = start_x; x < end_x; x++) {
				ElevationData[y][x] = ElevationDataFiltered[y][x];
			}
		}
	}

	for (int i = 0; i < row; i++) {
		delete[] ElevationDataFiltered[i];
	}
	delete[] ElevationDataFiltered;
}

#pragma warning(disable:4996)
unmanned_ground_vehicle *UGV = new unmanned_ground_vehicle;

void UGV_run(GINPUT *GInput) {

	int i = 0;
	for (int i = 0; i < 30; i++) {
		char dir[255];
		sprintf_s(dir, sizeof(char) * 255, "simulation_results/%d", i + 1);
		mkdir(dir);
		UGV->sim_count = i + 1;
		UGV->init(GInput);
		UGV->run();
	}
}

int main() {

	mkdir("simulation_results");

	GINPUT *GInput = new GINPUT;
	int i, j;

	UGV->HilbertMap = false;

	FILE *fp_map;
	FILE *east_path, *north_path;
	char *ptr_map, *basic_map, *token_map;
	int buffer_map = 10000000;
	basic_map = new char[buffer_map];
	//fopen_s(&fp_map, "hilbertmapelevation.csv", "r");
	//fopen_s(&fp_map, "animation_map.csv", "r");

	fopen_s(&fp_map, "map.csv", "r");
	fopen_s(&east_path, "east_temp.csv", "r");
	fopen_s(&north_path, "north_temp.csv", "r");
	GInput->MapInfo_Resolution_Row = 0.1;
	GInput->MapInfo_Resolution_Col = 0.1;
	GInput->MapInfo_Row = 11700;
	GInput->MapInfo_Col = 10800;

	GInput->ElevationData = new double*[GInput->MapInfo_Row];
	for (int k = 0; k < GInput->MapInfo_Row; k++) {
		GInput->ElevationData[k] = new double[GInput->MapInfo_Col];
	}

	i = 0;
	while (fgets(basic_map, buffer_map, fp_map) != NULL) {
		j = 0;
		ptr_map = strtok_s(basic_map, ",", &token_map);
		while (ptr_map != NULL) {
			GInput->ElevationData[i][j] = atof(ptr_map);
			ptr_map = strtok_s(NULL, ",", &token_map);
			j++;
		}
		i++;
	}
	delete[] basic_map;
	fclose(fp_map);

	//////////////////////////////////////////////////////////////////////////
	/////////////////////////// 맵 데이터 수정 ////////////////////////////////
	int row_min, row_max, col_min, col_max;
	GaussianBlur77(GInput->ElevationData, 0, GInput->MapInfo_Row, 0, GInput->MapInfo_Col, GInput->MapInfo_Row, GInput->MapInfo_Col);

	GaussianBlur77(GInput->ElevationData, 10400, 10900, 850, 1250, GInput->MapInfo_Row, GInput->MapInfo_Col, 20);

	row_min = 8660; row_max = 8690; col_min = 2595; col_max = 2610;
	for (int i = 8660; i < 8690; i++) {
		for (int j = 2595; j < 2610; j++) {
			GInput->ElevationData[i][j] = GInput->ElevationData[i][2615];
		}
	}
	GaussianBlur77(GInput->ElevationData, row_min - 10, row_max + 10, col_min - 10, col_max + 10, GInput->MapInfo_Row, GInput->MapInfo_Col, 3);

	row_min = 10410; row_max = 10430; col_min = 1120; col_max = 1140;
	for (int i = 10410; i < 10430; i++) {
		for (int j = 1120; j < 1140; j++) {
			GInput->ElevationData[i][j] = GInput->ElevationData[10390][1121];
		}
	}
	GaussianBlur77(GInput->ElevationData, row_min - 10, row_max + 10, col_min - 10, col_max + 10, GInput->MapInfo_Row, GInput->MapInfo_Col, 3);

	row_min = 6850; row_max = 6860; col_min = 3740; col_max = 3750;
	for (int i = 6850; i < 6860; i++) {
		for (int j = 3740; j < 3750; j++) {
			GInput->ElevationData[i][j] = GInput->ElevationData[i][3735];
		}
	}
	GaussianBlur77(GInput->ElevationData, row_min - 10, row_max + 10, col_min - 10, col_max + 10, GInput->MapInfo_Row, GInput->MapInfo_Col, 3);

	row_min = 4390; row_max = 4410; col_min = 5815; col_max = 5825;
	for (int i = 4390; i < 4410; i++) {
		for (int j = 5815; j < 5825; j++) {
			GInput->ElevationData[i][j] = GInput->ElevationData[i][5830];
		}
	}
	GaussianBlur77(GInput->ElevationData, row_min - 10, row_max + 10, col_min - 10, col_max + 10, GInput->MapInfo_Row, GInput->MapInfo_Col, 3);

	row_min = 3760; row_max = 3790; col_min = 6560; col_max = 6585;
	for (int i = 3760; i < 3790; i++) {
		for (int j = 6560; j < 6585; j++) {
			GInput->ElevationData[i][j] = GInput->ElevationData[i][6560];
		}
	}
	GaussianBlur77(GInput->ElevationData, row_min - 10, row_max + 10, col_min - 10, col_max + 10, GInput->MapInfo_Row, GInput->MapInfo_Col, 3);

	row_min = 2540; row_max = 2560; col_min = 8190; col_max = 8210;
	for (int i = 2540; i < 2560; i++) {
		for (int j = 8190; j < 8210; j++) {
			GInput->ElevationData[i][j] = GInput->ElevationData[i][8180];
		}
	}
	GaussianBlur77(GInput->ElevationData, row_min - 10, row_max + 10, col_min - 10, col_max + 10, GInput->MapInfo_Row, GInput->MapInfo_Col, 3);

	row_min = 2360; row_max = 2385; col_min = 8660; col_max = 8680;
	for (int i = 2360; i < 2385; i++) {
		for (int j = 8660; j < 8680; j++) {
			GInput->ElevationData[i][j] = GInput->ElevationData[i][8650];
		}
	}
	GaussianBlur77(GInput->ElevationData, row_min - 10, row_max + 10, col_min - 10, col_max + 10, GInput->MapInfo_Row, GInput->MapInfo_Col, 3);

	row_min = 1790; row_max = 1810; col_min = 9530; col_max = 9545;
	for (int i = 1790; i < 1810; i++) {
		for (int j = 9530; j < 9545; j++) {
			GInput->ElevationData[i][j] = GInput->ElevationData[i][9520];
		}
	}
	GaussianBlur77(GInput->ElevationData, row_min - 10, row_max + 10, col_min - 10, col_max + 10, GInput->MapInfo_Row, GInput->MapInfo_Col, 3);

	GaussianBlur77(GInput->ElevationData, 680, 740, 9900, 9980, GInput->MapInfo_Row, GInput->MapInfo_Col, 3);
	
	row_min = 3140; row_max = 3165; col_min = 6200; col_max = 6230;
	for (int i = 3140; i < 3165; i++) {
		for (int j = 6200; j < 6230; j++) {
			GInput->ElevationData[i][j] = GInput->ElevationData[3140][j];
		}
	}
	GaussianBlur77(GInput->ElevationData, row_min - 10, row_max + 10, col_min - 10, col_max + 10, GInput->MapInfo_Row, GInput->MapInfo_Col, 3);

	row_min = 4220; row_max = 4290; col_min = 5120; col_max = 5170;
	for (int i = 4220; i < 4290; i++) {
		for (int j = 5120; j < 5170; j++) {
			GInput->ElevationData[i][j] = GInput->ElevationData[4250][5169];
		}
	}
	GaussianBlur77(GInput->ElevationData, row_min - 10, row_max + 10, col_min - 10, col_max + 10, GInput->MapInfo_Row, GInput->MapInfo_Col, 3);

	row_min = 4400; row_max = 4470; col_min = 4685; col_max = 4730;
	for (int i = 4400; i < 4470; i++) {
		for (int j = 4685; j < 4730; j++) {
			GInput->ElevationData[i][j] = GInput->ElevationData[i][4680];
		}
	}
	GaussianBlur77(GInput->ElevationData, row_min - 10, row_max + 10, col_min - 10, col_max + 10, GInput->MapInfo_Row, GInput->MapInfo_Col, 3);

	GaussianBlur77(GInput->ElevationData, 4390, 4480, 4675, 4740, GInput->MapInfo_Row, GInput->MapInfo_Col, 3);
	
	row_min = 9810; row_max = 9840; col_min = 860; col_max = 890;
	for (int i = 9810; i < 9840; i++) {
		for (int j = 860; j < 890; j++) {
			GInput->ElevationData[i][j] = GInput->ElevationData[9810][j];
		}
	}
	GaussianBlur77(GInput->ElevationData, row_min - 10, row_max + 10, col_min - 10, col_max + 10, GInput->MapInfo_Row, GInput->MapInfo_Col, 3);

	row_min = 650; row_max = 670; col_min = 9655; col_max = 9680;
	for (int i = 650; i < 670; i++) {
		for (int j = 9655; j < 9680; j++) {
			GInput->ElevationData[i][j] = GInput->ElevationData[645][j];
		}
	}
	GaussianBlur77(GInput->ElevationData, row_min - 10, row_max + 10, col_min - 10, col_max + 10, GInput->MapInfo_Row, GInput->MapInfo_Col, 3);

	row_min = 4190; row_max = 4210; col_min = 5300; col_max = 5330;
	for (int i = 4190; i < 4210; i++) {
		for (int j = 5300; j < 5330; j++) {
			GInput->ElevationData[i][j] = GInput->ElevationData[4180][j];
		}
	}
	GaussianBlur77(GInput->ElevationData, row_min - 10, row_max + 10, col_min - 10, col_max + 10, GInput->MapInfo_Row, GInput->MapInfo_Col, 3);
	//////////////////////////////////////////////////////////////////////////

	GInput->WaypointSize = 2518;
	GInput->LocalPath = new double*[GInput->WaypointSize];
	for (int k = 0; k < GInput->WaypointSize; k++) {
		GInput->LocalPath[k] = new double[2];
	}

	char *ptr, *basic, *token;
	int buffer = 500000;
	basic = new char[buffer];

	while (fgets(basic, buffer, east_path) != NULL) {
		j = 0;
		ptr = strtok_s(basic, ",", &token);
		while (ptr != NULL) {
			GInput->LocalPath[j][0] = atof(ptr);
			ptr = strtok_s(NULL, ",", &token);
			j++;
		}
		break;
	}
	fclose(east_path);

	while (fgets(basic, buffer, north_path) != NULL) {
		j = 0;
		ptr = strtok_s(basic, ",", &token);
		while (ptr != NULL) {
			GInput->LocalPath[j][1] = atof(ptr);
			ptr = strtok_s(NULL, ",", &token);
			j++;
		}
		break;
	}
	delete[] basic;
	fclose(north_path);

	UGV_run(GInput);

	return 0;
}