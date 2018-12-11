#include <stdio.h>
#include <direct.h>
#include "unmanned_ground_vehicle.h"
void GaussianBlur(double** ElevationData, int row, int col) {
	const int box[7 * 7] = {
		1,  4,  7,  10, 7,  4,  1,
		4,  16, 26, 33, 26, 16, 4,
		7,  26, 55, 71, 55, 26, 7,
		10, 33, 71, 91, 71, 33, 10,
		7,  26, 55, 71, 55, 26, 7,
		4,  16, 26, 33, 26, 16, 4,
		1,  4,  7,  10, 7,  4,  1
	};

	int start_y = 0;
	int start_x = 0;
	int end_y = row;
	int end_x = col;
	double sum = 0;
	int height = row;
	int width = col;
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

	for (int y = start_y; y < end_y; y++) {
		for (int x = start_x; x < end_x; x++) {
			sum = 0;
			for (int i = -offset, ik = 0; i <= offset; i++, ik++) {
				int ii = y + i < 0 ? 0 : y + i >= height ? height - 1 : y + i;
				for (int j = -offset, jk = 0; j <= offset; j++, jk++) {
					int jj = x + j < 0 ? 0 : x + j >= width ? width - 1 : x + j;
					sum += ElevationData[ii][jj] * box[ik*ksize + jk];
				}
			}
			ElevationDataFiltered[y][x] = sum / scale;
		}
	}

	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			ElevationData[i][j] = ElevationDataFiltered[i][j];
		}
	}

	for (int i = 0; i < row; i++) {
		delete[] ElevationDataFiltered[i];
	}
	delete[] ElevationDataFiltered;
}
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
					for (int j = -offset, jk = 0; j <= offset; j++, jk++) {
						sum += ElevationData[y + i][x + j] * box[ik*ksize + jk];
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
	//for (int i = 28; i < 30; i++) {
	char dir[255];
	sprintf_s(dir, sizeof(char) * 255, "simulation_results/%d", i + 1);
	mkdir(dir);
	UGV->sim_count = i + 1;
	UGV->init(GInput);
	UGV->run();
	//}
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

	GInput->flag = 1; // 1 : 10cm, 2 : 50cm, 3 : 1m

	if (GInput->flag == 1) {
		fopen_s(&fp_map, "map.csv", "r");
		fopen_s(&east_path, "east_temp.csv", "r");
		fopen_s(&north_path, "north_temp.csv", "r");
		GInput->MapInfo_Resolution_Row = 0.1;
		GInput->MapInfo_Resolution_Col = 0.1;
		GInput->MapInfo_Row = 11700;
		GInput->MapInfo_Col = 10800;
	}
	else if (GInput->flag == 2) {
		fopen_s(&fp_map, "map_50cm.csv", "r");
		fopen_s(&east_path, "east_50cm.csv", "r");
		fopen_s(&north_path, "north_50cm.csv", "r");
		GInput->MapInfo_Resolution_Row = 0.5;
		GInput->MapInfo_Resolution_Col = 0.5;
		GInput->MapInfo_Row = 2340;
		GInput->MapInfo_Col = 2161;
	}
	else if (GInput->flag == 3) {
		fopen_s(&fp_map, "map_1m.csv", "r");
		fopen_s(&east_path, "east.csv", "r");
		fopen_s(&north_path, "north.csv", "r");
		GInput->MapInfo_Resolution_Row = 1;
		GInput->MapInfo_Resolution_Col = 1;
		GInput->MapInfo_Row = 1170;
		GInput->MapInfo_Col = 1080;
	}
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
	GaussianBlur(GInput->ElevationData, GInput->MapInfo_Row, GInput->MapInfo_Col);
	GaussianBlur77(GInput->ElevationData, 10400, 10900, 850, 1250, GInput->MapInfo_Row, GInput->MapInfo_Col, 20);
	for (int i = 8660; i < 8690; i++) {
		for (int j = 2595; j < 2610; j++) {
			GInput->ElevationData[i][j] = GInput->ElevationData[i][2615];
		}
	}
	for (int i = 10410; i < 10430; i++) {
		for (int j = 1120; j < 1140; j++) {
			GInput->ElevationData[i][j] = GInput->ElevationData[10390][1121];
		}
	}
	for (int i = 6850; i < 6860; i++) {
		for (int j = 3740; j < 3750; j++) {
			GInput->ElevationData[i][j] = GInput->ElevationData[i][3735];
		}
	}
	for (int i = 4390; i < 4410; i++) {
		for (int j = 5815; j < 5825; j++) {
			GInput->ElevationData[i][j] = GInput->ElevationData[i][5830];
		}
	}
	for (int i = 3760; i < 3790; i++) {
		for (int j = 6560; j < 6585; j++) {
			GInput->ElevationData[i][j] = GInput->ElevationData[i][6560];
		}
	}
	for (int i = 2540; i < 2560; i++) {
		for (int j = 8190; j < 8210; j++) {
			GInput->ElevationData[i][j] = GInput->ElevationData[i][8180];
		}
	}
	for (int i = 2360; i < 2385; i++) {
		for (int j = 8660; j < 8680; j++) {
			GInput->ElevationData[i][j] = GInput->ElevationData[i][8650];
		}
	}
	for (int i = 1790; i < 1810; i++) {
		for (int j = 9530; j < 9545; j++) {
			GInput->ElevationData[i][j] = GInput->ElevationData[i][9520];
		}
	}
	GaussianBlur77(GInput->ElevationData, 680, 740, 9900, 9980, GInput->MapInfo_Row, GInput->MapInfo_Col);
	for (int i = 3140; i < 3165; i++) {
		for (int j = 6200; j < 6230; j++) {
			GInput->ElevationData[i][j] = GInput->ElevationData[3140][j];
		}
	}
	for (int i = 4220; i < 4290; i++) {
		for (int j = 5120; j < 5170; j++) {
			GInput->ElevationData[i][j] = GInput->ElevationData[4250][5169];
		}
	}
	for (int i = 4400; i < 4470; i++) {
		for (int j = 4685; j < 4730; j++) {
			GInput->ElevationData[i][j] = GInput->ElevationData[i][4680];
		}
	}
	GaussianBlur77(GInput->ElevationData, 4390, 4480, 4675, 4740, GInput->MapInfo_Row, GInput->MapInfo_Col, 3);
	for (int i = 9810; i < 9840; i++) {
		for (int j = 860; j < 890; j++) {
			GInput->ElevationData[i][j] = GInput->ElevationData[9810][j];
		}
	}
	for (int i = 650; i < 670; i++) {
		for (int j = 9655; j < 9680; j++) {
			GInput->ElevationData[i][j] = GInput->ElevationData[645][j];
		}
	}
	for (int i = 4190; i < 4210; i++) {
		for (int j = 5300; j < 53330; j++) {
			GInput->ElevationData[i][j] = GInput->ElevationData[4180][j];
		}
	}
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