clc; clear all; close all;

simulation_indx = 0;
addpath(sprintf('%d',simulation_indx));

map = load('../map_10cm_anim.csv'); % 애니메이션을 위해 1m 간격으로 별도로 저장한 맵 데이터
path = load('../3d_path_10cm.csv'); % 애니메니션을 위해 x,y,z 데이터로 저장한 경로 데이터

data = csvread('SimulationResult.csv',1,0); % 시뮬레이션 결과 데이터
rtt_data = csvread('TraversabilityAnalysisResult.csv',1,0); % 주행성 분석 결과 데이터

data_size = size(data,1);

wp_size = rtt_data(1,2);
num_threads = rtt_data(1,3);
wp = 1:wp_size;

RSM_max = 0;
PSM_max = 0;
LSM_max = 0;
VSM_max = 0;

indx = 0;

start_wp = 1;

% Video file open
makeVideo = VideoWriter('TraversabilityAnalysisResultVideo');
% Frame Rate - 숫자가 클 수록 재생 속도가 빠름
makeVideo.FrameRate = 10;
% Quality - 용량과 관련 됨 (0 ~ 100)
makeVideo.Quality = 80;
open(makeVideo);

%% plot
figure('Position',[100,100,1720,880])
set(gcf,'Color',[1,1,1])

%% ugv animation
subplot(4,4,[1,2,5,6])
x = linspace(0,1080,10800/2);
y = linspace(0,1170,11700/2);
[X,Y] = meshgrid(x,y);
mesh(X,Y,map);
hold on
grid on
plot3(path(:,1),path(:,2),path(:,3),'or','MarkerSize',5,'MarkerFaceColor','r')
xlabel('X [m]','FontSize',15), ylabel('Y [m]','FontSize',15), zlabel('Z [m]','FontSize',15)

chassis_length = 3.2; chassis_width = 1.1996;
chassis_depth_1 = 1.0; chassis_depth_2 = 0.3;
sus_length_1 = 0.7; sus_length_2 = 0.11; sus_width = 0.3*0.8; sus_depth = 0.12*2;
tire_radius = 0.548; tire_width = 0.372; tire_offset = 0.219;

[tire_x,tire_y,tire_z] = cylinder(tire_radius);

RF = [0	0	0.586234043	0	0.0932	-6.48E-12	0.865112866	-0.501577241	1	3.04E-12	-7.68E-12	-5.12E-12	-0.501577241	-0.865112866	376	31.25706362	0	0	0	25.98627239	0	0	0	20.83704122	1.13	-0.755	-0.108	0.65	0	0.219	1	0	0	0	-5.10E-12	-1	0	1	-5.10E-12	1	0	0	0	1	0	0	0	1	1.276	-0.736	0.658	0.281	0.133	-0.019];
RM = [0	0	0.585382979	0	0.0932	1.53E-11	-0.868782601	-0.495193692	-1	-1.33E-11	-7.58E-12	0.00E+00	0.495193692	-0.868782601	376	31.21662106	0	0	0	25.94143616	0	0	0	20.8414349	0.83	-0.755	-0.108	0.648	0	0.219	-1	1.02E-11	0	5.21E-23	5.10E-12	-1	-1.02E-11	-1	-5.10E-12	1	0	0	0	1	0	0	0	1	0.683	-0.736	0.658	0.281	-0.133	-1.90E-02];
RR = [0	0	0.585382979	0	0.0932	1.53E-11	-0.868782601	-0.495193692	-1	-1.33E-11	-7.58E-12	0.00E+00	0.495193692	-0.868782601	376	31.21662106	0	0	0	25.94143616	0	0	0	20.8414349	-0.77	-0.755	-0.108	0.648	0	0.219	-1	1.02E-11	0	5.21E-23	5.10E-12	-1	-1.02E-11	-1	-5.10E-12	1	0	0	0	1	0	0	0	1	-0.917	-0.736	0.658	0.281	-0.133	-1.90E-02];
LF = [0	0	0.585864863	0	-0.0932	-5.10E-12	0.866708022	0.498815802	1	4.42E-12	2.55E-12	0.00E+00	0.498815802	-0.866708022	376	31.23950966	0	0	0	25.96680966	0	0	0	20.83895001	1.13	0.755	-0.108	0.65	0	-0.219	1	0.00E+00	0	0.00E+00	-5.10E-12	-1	0.00E+00	1	-5.10E-12	1	0	0	0	1	0	0	0	1	1.276	0.736	0.658	0.281	0.133	1.90E-02];
LM = [0	0	0.585752158	0	-0.0932	-3.75E-12	-0.867194003	0.497970441	-1	5.78E-12	2.54E-12	-5.08E-12	-0.497970441	-0.867194003	376	31.23415389	0	0	0	25.96087203	0	0	0	20.83953186	0.83	0.755	-0.108	0.648	0	-0.219	-1	1.02E-11	0	5.21E-23	5.10E-12	-1	-1.02E-11	-1	-5.10E-12	1	0	0	0	1	0	0	0	1	0.683	0.736	0.658	0.281	-0.133	1.90E-02];
LR = [0	0	0.585382979	0	-0.0932	-3.76E-12	-0.868782601	0.495193692	-1	5.77E-12	2.53E-12	-5.05E-12	-0.495193692	-0.868782601	376	31.21662106	0	0	0	25.94143616	0	0	0	20.8414349	-0.77	0.755	-0.108	0.648	0	-0.219	-1	1.02E-11	0	5.21E-23	5.10E-12	-1	-1.02E-11	-1	-5.10E-12	1	0	0	0	1	0	0	0	1	-0.917	0.736	0.658	0.281	-0.133	1.90E-02];

[sus_RF_s01p, sus_RF_s12p, sus_RF_C01, sus_RF_C12, sus_RF_s0sp, sus_RF_s1sp] = ugv.sus_data_input(RF');
[sus_RM_s01p, sus_RM_s12p, sus_RM_C01, sus_RM_C12, sus_RM_s0sp, sus_RM_s1sp] = ugv.sus_data_input(RM');
[sus_RR_s01p, sus_RR_s12p, sus_RR_C01, sus_RR_C12, sus_RR_s0sp, sus_RR_s1sp] = ugv.sus_data_input(RR');
[sus_LF_s01p, sus_LF_s12p, sus_LF_C01, sus_LF_C12, sus_LF_s0sp, sus_LF_s1sp] = ugv.sus_data_input(LF');
[sus_LM_s01p, sus_LM_s12p, sus_LM_C01, sus_LM_C12, sus_LM_s0sp, sus_LM_s1sp] = ugv.sus_data_input(LM');
[sus_LR_s01p, sus_LR_s12p, sus_LR_C01, sus_LR_C12, sus_LR_s0sp, sus_LR_s1sp] = ugv.sus_data_input(LR');

init_flag = 1;
for i = 1 : 20 : data_size
    if (data(i,30) >= start_wp)
        subplot(4,4,[1,2,5,6])
        chassis_r0 = data(i,2:4)';
        chassis_A0 = [data(i,5:7);data(i,8:10);data(i,11:13)];

        chassis_RFT = chassis_r0 + chassis_A0*[chassis_length/2;-chassis_width/2;chassis_depth_1];
        chassis_LFT = chassis_r0 + chassis_A0*[chassis_length/2;chassis_width/2;chassis_depth_1];
        chassis_LBT = chassis_r0 + chassis_A0*[-chassis_length/2;chassis_width/2;chassis_depth_1];
        chassis_RBT = chassis_r0 + chassis_A0*[-chassis_length/2;-chassis_width/2;chassis_depth_1];
        chassis_RFD = chassis_r0 + chassis_A0*[chassis_length/2;-chassis_width/2;-chassis_depth_2];
        chassis_LFD = chassis_r0 + chassis_A0*[chassis_length/2;chassis_width/2;-chassis_depth_2];
        chassis_LBD = chassis_r0 + chassis_A0*[-chassis_length/2;chassis_width/2;-chassis_depth_2];
        chassis_RBD = chassis_r0 + chassis_A0*[-chassis_length/2;-chassis_width/2;-chassis_depth_2];

        sus_RF_q1 = data(i,14); sus_RF_angle = data(i,20);
        sus_RM_q1 = data(i,15); sus_RM_angle = data(i,21);
        sus_RR_q1 = data(i,16); sus_RR_angle = data(i,22);
        sus_LF_q1 = data(i,17); sus_LF_angle = data(i,23);
        sus_LM_q1 = data(i,18); sus_LM_angle = data(i,24);
        sus_LR_q1 = data(i,19); sus_LR_angle = data(i,25);

        [sus_RF_ori, sus_RF_wheel_ori] = ugv.sus_ori(sus_RF_q1, sus_RF_angle);
        [sus_RM_ori, sus_RM_wheel_ori] = ugv.sus_ori(sus_RM_q1, sus_RM_angle);
        [sus_RR_ori, sus_RR_wheel_ori] = ugv.sus_ori(sus_RR_q1, sus_RR_angle);
        [sus_LF_ori, sus_LF_wheel_ori] = ugv.sus_ori(sus_LF_q1, sus_LF_angle);
        [sus_LM_ori, sus_LM_wheel_ori] = ugv.sus_ori(sus_LM_q1, sus_LM_angle);
        [sus_LR_ori, sus_LR_wheel_ori] = ugv.sus_ori(sus_LR_q1, sus_LR_angle);

        [sus_RF_RFT, sus_RF_LFT, sus_RF_LBT, sus_RF_RBT, sus_RF_RFD, sus_RF_LFD, sus_RF_LBD, sus_RF_RBD] = ugv.sus_geo(chassis_r0, chassis_A0, sus_RF_s01p, sus_RF_C01, sus_RF_ori,sus_length_1, sus_length_2, sus_depth, sus_width);
        [sus_RM_RFT, sus_RM_LFT, sus_RM_LBT, sus_RM_RBT, sus_RM_RFD, sus_RM_LFD, sus_RM_LBD, sus_RM_RBD] = ugv.sus_geo(chassis_r0, chassis_A0, sus_RM_s01p, sus_RM_C01, sus_RM_ori,sus_length_1, sus_length_2, sus_depth, sus_width);
        [sus_RR_RFT, sus_RR_LFT, sus_RR_LBT, sus_RR_RBT, sus_RR_RFD, sus_RR_LFD, sus_RR_LBD, sus_RR_RBD] = ugv.sus_geo(chassis_r0, chassis_A0, sus_RR_s01p, sus_RR_C01, sus_RR_ori,sus_length_1, sus_length_2, sus_depth, sus_width);
        [sus_LF_RFT, sus_LF_LFT, sus_LF_LBT, sus_LF_RBT, sus_LF_RFD, sus_LF_LFD, sus_LF_LBD, sus_LF_RBD] = ugv.sus_geo(chassis_r0, chassis_A0, sus_LF_s01p, sus_LF_C01, sus_LF_ori,sus_length_1, sus_length_2, sus_depth, -sus_width);
        [sus_LM_RFT, sus_LM_LFT, sus_LM_LBT, sus_LM_RBT, sus_LM_RFD, sus_LM_LFD, sus_LM_LBD, sus_LM_RBD] = ugv.sus_geo(chassis_r0, chassis_A0, sus_LM_s01p, sus_LM_C01, sus_LM_ori,sus_length_1, sus_length_2, sus_depth, -sus_width);
        [sus_LR_RFT, sus_LR_LFT, sus_LR_LBT, sus_LR_RBT, sus_LR_RFD, sus_LR_LFD, sus_LR_LBD, sus_LR_RBD] = ugv.sus_geo(chassis_r0, chassis_A0, sus_LR_s01p, sus_LR_C01, sus_LR_ori,sus_length_1, sus_length_2, sus_depth, -sus_width);

        [chassis_top, chassis_bot, chassis_right, chassis_left, chassis_front, chassis_back] = ugv.body_vertices(chassis_RFT, chassis_LFT, chassis_LBT, chassis_RBT, chassis_RFD, chassis_LFD, chassis_LBD, chassis_RBD);
        [sus_RF_top, sus_RF_bot, sus_RF_right, sus_RF_left, sus_RF_front, sus_RF_back] = ugv.body_vertices(sus_RF_RFT, sus_RF_LFT, sus_RF_LBT, sus_RF_RBT, sus_RF_RFD, sus_RF_LFD, sus_RF_LBD, sus_RF_RBD);
        [sus_RM_top, sus_RM_bot, sus_RM_right, sus_RM_left, sus_RM_front, sus_RM_back] = ugv.body_vertices(sus_RM_RFT, sus_RM_LFT, sus_RM_LBT, sus_RM_RBT, sus_RM_RFD, sus_RM_LFD, sus_RM_LBD, sus_RM_RBD);
        [sus_RR_top, sus_RR_bot, sus_RR_right, sus_RR_left, sus_RR_front, sus_RR_back] = ugv.body_vertices(sus_RR_RFT, sus_RR_LFT, sus_RR_LBT, sus_RR_RBT, sus_RR_RFD, sus_RR_LFD, sus_RR_LBD, sus_RR_RBD);
        [sus_LF_top, sus_LF_bot, sus_LF_right, sus_LF_left, sus_LF_front, sus_LF_back] = ugv.body_vertices(sus_LF_RFT, sus_LF_LFT, sus_LF_LBT, sus_LF_RBT, sus_LF_RFD, sus_LF_LFD, sus_LF_LBD, sus_LF_RBD);
        [sus_LM_top, sus_LM_bot, sus_LM_right, sus_LM_left, sus_LM_front, sus_LM_back] = ugv.body_vertices(sus_LM_RFT, sus_LM_LFT, sus_LM_LBT, sus_LM_RBT, sus_LM_RFD, sus_LM_LFD, sus_LM_LBD, sus_LM_RBD);
        [sus_LR_top, sus_LR_bot, sus_LR_right, sus_LR_left, sus_LR_front, sus_LR_back] = ugv.body_vertices(sus_LR_RFT, sus_LR_LFT, sus_LR_LBT, sus_LR_RBT, sus_LR_RFD, sus_LR_LFD, sus_LR_LBD, sus_LR_RBD);

        [tire_RF_center, tire_RF2_center, tire_RF_vertices, tire_RF2_vertices] = ugv.tire_geo(chassis_r0, chassis_A0, sus_RF_s01p, sus_RF_C01, sus_RF_ori, sus_RF_wheel_ori, sus_RF_s12p, tire_offset, tire_width, tire_x, tire_y, tire_z);
        [tire_RM_center, tire_RM2_center, tire_RM_vertices, tire_RM2_vertices] = ugv.tire_geo(chassis_r0, chassis_A0, sus_RM_s01p, sus_RM_C01, sus_RM_ori, sus_RM_wheel_ori, sus_RM_s12p, tire_offset, tire_width, tire_x, tire_y, tire_z);
        [tire_RR_center, tire_RR2_center, tire_RR_vertices, tire_RR2_vertices] = ugv.tire_geo(chassis_r0, chassis_A0, sus_RR_s01p, sus_RR_C01, sus_RR_ori, sus_RR_wheel_ori, sus_RR_s12p, tire_offset, tire_width, tire_x, tire_y, tire_z);
        [tire_LF_center, tire_LF2_center, tire_LF_vertices, tire_LF2_vertices] = ugv.tire_geo(chassis_r0, chassis_A0, sus_LF_s01p, sus_LF_C01, sus_LF_ori, sus_LF_wheel_ori, sus_LF_s12p, -tire_offset, -tire_width, tire_x, tire_y, tire_z);
        [tire_LM_center, tire_LM2_center, tire_LM_vertices, tire_LM2_vertices] = ugv.tire_geo(chassis_r0, chassis_A0, sus_LM_s01p, sus_LM_C01, sus_LM_ori, sus_LM_wheel_ori, sus_LM_s12p, -tire_offset, -tire_width, tire_x, tire_y, tire_z);
        [tire_LR_center, tire_LR2_center, tire_LR_vertices, tire_LR2_vertices] = ugv.tire_geo(chassis_r0, chassis_A0, sus_LR_s01p, sus_LR_C01, sus_LR_ori, sus_LR_wheel_ori, sus_LR_s12p, -tire_offset, -tire_width, tire_x, tire_y, tire_z);

        if init_flag == 1
            box_chassis_top = patch('XData',chassis_top(:,1)','YData',chassis_top(:,2)','ZData',chassis_top(:,3)','FaceColor','y');
            box_chassis_bot = patch('XData',chassis_bot(:,1)','YData',chassis_bot(:,2)','ZData',chassis_bot(:,3)','FaceColor','y');
            box_chassis_right = patch('XData',chassis_right(:,1),'YData',chassis_right(:,2),'ZData',chassis_right(:,3),'FaceColor','c');
            box_chassis_left = patch('XData',chassis_left(:,1),'YData',chassis_left(:,2),'ZData',chassis_left(:,3),'FaceColor','c');
            box_chassis_front = patch('XData',chassis_front(:,1),'YData',chassis_front(:,2),'ZData',chassis_front(:,3),'FaceColor','g');
            box_chassis_back = patch('XData',chassis_back(:,1),'YData',chassis_back(:,2),'ZData',chassis_back(:,3),'FaceColor','g');

            [box_sus_RF_top, box_sus_RF_bot, box_sus_RF_right, box_sus_RF_left, box_sus_RF_front, box_sus_RF_back] = ugv.sus_display(sus_RF_top, sus_RF_bot, sus_RF_right, sus_RF_left, sus_RF_front, sus_RF_back);
            [box_sus_RM_top, box_sus_RM_bot, box_sus_RM_right, box_sus_RM_left, box_sus_RM_front, box_sus_RM_back] = ugv.sus_display(sus_RM_top, sus_RM_bot, sus_RM_right, sus_RM_left, sus_RM_front, sus_RM_back);
            [box_sus_RR_top, box_sus_RR_bot, box_sus_RR_right, box_sus_RR_left, box_sus_RR_front, box_sus_RR_back] = ugv.sus_display(sus_RR_top, sus_RR_bot, sus_RR_right, sus_RR_left, sus_RR_front, sus_RR_back);
            [box_sus_LF_top, box_sus_LF_bot, box_sus_LF_right, box_sus_LF_left, box_sus_LF_front, box_sus_LF_back] = ugv.sus_display(sus_LF_top, sus_LF_bot, sus_LF_right, sus_LF_left, sus_LF_front, sus_LF_back);
            [box_sus_LM_top, box_sus_LM_bot, box_sus_LM_right, box_sus_LM_left, box_sus_LM_front, box_sus_LM_back] = ugv.sus_display(sus_LM_top, sus_LM_bot, sus_LM_right, sus_LM_left, sus_LM_front, sus_LM_back);
            [box_sus_LR_top, box_sus_LR_bot, box_sus_LR_right, box_sus_LR_left, box_sus_LR_front, box_sus_LR_back] = ugv.sus_display(sus_LR_top, sus_LR_bot, sus_LR_right, sus_LR_left, sus_LR_front, sus_LR_back);

            [tire_RF, wheel_RF, wheel_RF2] = ugv.tire_display(tire_RF_vertices, tire_RF2_vertices, tire_RF_center, tire_RF2_center);
            [tire_RM, wheel_RM, wheel_RM2] = ugv.tire_display(tire_RM_vertices, tire_RM2_vertices, tire_RM_center, tire_RM2_center);
            [tire_RR, wheel_RR, wheel_RR2] = ugv.tire_display(tire_RR_vertices, tire_RR2_vertices, tire_RR_center, tire_RR2_center);
            [tire_LF, wheel_LF, wheel_LF2] = ugv.tire_display(tire_LF_vertices, tire_LF2_vertices, tire_LF_center, tire_LF2_center);
            [tire_LM, wheel_LM, wheel_LM2] = ugv.tire_display(tire_LM_vertices, tire_LM2_vertices, tire_LM_center, tire_LM2_center);
            [tire_LR, wheel_LR, wheel_LR2] = ugv.tire_display(tire_LR_vertices, tire_LR2_vertices, tire_LR_center, tire_LR2_center);
            init_flag = 0;
        else
            ugv.body_move(chassis_top, chassis_bot, chassis_right, chassis_left, chassis_front, chassis_back, box_chassis_top, box_chassis_bot, box_chassis_right, box_chassis_left, box_chassis_front, box_chassis_back);

            ugv.body_move(sus_RF_top, sus_RF_bot, sus_RF_right, sus_RF_left, sus_RF_front, sus_RF_back, box_sus_RF_top, box_sus_RF_bot, box_sus_RF_right, box_sus_RF_left, box_sus_RF_front, box_sus_RF_back);
            ugv.body_move(sus_RM_top, sus_RM_bot, sus_RM_right, sus_RM_left, sus_RM_front, sus_RM_back, box_sus_RM_top, box_sus_RM_bot, box_sus_RM_right, box_sus_RM_left, box_sus_RM_front, box_sus_RM_back);
            ugv.body_move(sus_RR_top, sus_RR_bot, sus_RR_right, sus_RR_left, sus_RR_front, sus_RR_back, box_sus_RR_top, box_sus_RR_bot, box_sus_RR_right, box_sus_RR_left, box_sus_RR_front, box_sus_RR_back);
            ugv.body_move(sus_LF_top, sus_LF_bot, sus_LF_right, sus_LF_left, sus_LF_front, sus_LF_back, box_sus_LF_top, box_sus_LF_bot, box_sus_LF_right, box_sus_LF_left, box_sus_LF_front, box_sus_LF_back);
            ugv.body_move(sus_LM_top, sus_LM_bot, sus_LM_right, sus_LM_left, sus_LM_front, sus_LM_back, box_sus_LM_top, box_sus_LM_bot, box_sus_LM_right, box_sus_LM_left, box_sus_LM_front, box_sus_LM_back);
            ugv.body_move(sus_LR_top, sus_LR_bot, sus_LR_right, sus_LR_left, sus_LR_front, sus_LR_back, box_sus_LR_top, box_sus_LR_bot, box_sus_LR_right, box_sus_LR_left, box_sus_LR_front, box_sus_LR_back);

            ugv.tire_move(tire_RF_vertices, tire_RF2_vertices, tire_RF_center, tire_RF2_center, tire_RF, wheel_RF, wheel_RF2);
            ugv.tire_move(tire_RM_vertices, tire_RM2_vertices, tire_RM_center, tire_RM2_center, tire_RM, wheel_RM, wheel_RM2);
            ugv.tire_move(tire_RR_vertices, tire_RR2_vertices, tire_RR_center, tire_RR2_center, tire_RR, wheel_RR, wheel_RR2);
            ugv.tire_move(tire_LF_vertices, tire_LF2_vertices, tire_LF_center, tire_LF2_center, tire_LF, wheel_LF, wheel_LF2);
            ugv.tire_move(tire_LM_vertices, tire_LM2_vertices, tire_LM_center, tire_LM2_center, tire_LM, wheel_LM, wheel_LM2);
            ugv.tire_move(tire_LR_vertices, tire_LR2_vertices, tire_LR_center, tire_LR2_center, tire_LR, wheel_LR, wheel_LR2);
        end

        view_size = 7;
        axis([chassis_r0(1)-view_size chassis_r0(1)+view_size chassis_r0(2)-view_size chassis_r0(2)+view_size chassis_r0(3)-view_size/2 chassis_r0(3)+view_size/2]);
        sim_status = sprintf('Simulation time: %5.2f s,  Velocity: %2.2f m/s',i*0.01, data(i,29));
        title(sim_status,'FontSize',15)
%         view(data(i,28)*180/pi+90, 45)

        %% global path & ugv trajectory
        subplot(4,4,[3,4,7,8])
        plot(path(:,1),path(:,2),'b','LineWidth',3)
        hold on
        plot(path(1:data(i,30),1), path(1:data(i,30),2),'c','LineWidth',3);
        plot(path(data(i,30),1), path(data(i,30),2),'ko','Markersize',7,'MarkerFaceColor','k');
        grid on
        xlabel('X [m]','FontSize',15), ylabel('Y [m]','FontSize',15)
        axis([0 1080 0 1170])
        hold off

        %% RTT analysis results
    %     if rtt_data(indx+1,4) <= data(i,1)
    %         indx = indx + 1;
    %     end
    %     subplot(4,4,[9,10])
    %     velocity = rtt_data(indx,6:6+wp_size-1);
    %     plot(wp',velocity','k','LineWidth',3)
    %     grid on
    %     ylabel('Velocity Profile [m/s]','FontSize',15)
    %     xlim([wp(1), wp(end)])
    %     set(gca,'XtickLabel','{}')

    %     subplot(4,4,[13,14])
    %     for j = 1 : num_threads
    %         thread(j).RSM = rtt_data(indx,46+j-1:28:1138+j-1);
    %         thread(j).PSM = rtt_data(indx,46+j-1+num_threads:28:1138+j-1+num_threads);
    %         thread(j).LSM = rtt_data(indx,46+j-1+num_threads*2:28:1138+j-1+num_threads*2);
    %         thread(j).VSM = rtt_data(indx,46+j-1+num_threads*3:28:1138+j-1+num_threads*3);
    % 
    %         plot(wp, thread(j).RSM,'ko','MarkerSize',5)
    %         hold on
    %         plot(wp, thread(j).PSM,'bo','MarkerSize',5)
    %         plot(wp, thread(j).LSM,'go','MarkerSize',5)
    %         plot(wp, thread(j).VSM,'mo','MarkerSize',5)
    %         plot([wp(1),wp(end)],[RSM_max, RSM_max],'k--','LineWidth',1.5)
    %         plot([wp(1),wp(end)],[PSM_max, PSM_max],'b--','LineWidth',1.5)
    %         plot([wp(1),wp(end)],[LSM_max, LSM_max],'g--','LineWidth',1.5)
    %         plot([wp(1),wp(end)],[VSM_max, VSM_max],'m--','LineWidth',1.5)
    %         grid on
    %         ylim([-1.5 1.5])
    %         xlim([wp(1), wp(end)])
    %         ylabel('Stability','FontSize',15)
    %         legend('RSM','PSM','LSM','VSM','RSM limit','PSM limit','LSM limit','VSM limit','Location','east')
    %     end
    %     hold off

        %% velocity command & response
        subplot(4,4,[11,12,15,16])
        plot(data(1:i,1),data(1:i,31),'b','LineWidth',3)
        hold on
        plot(data(1:i,1),data(1:i,32),'r:','LineWidth',3)
        grid on
        xlabel('Time [s]'), ylabel('Velocity [m/s]')
        legend('Velocity Command','Velocity Response','Location','NorthWest')
        set(gca,'FontSize',14)
        hold off

        pause(0.00001);
        disp([i, data(i,30)]);

        frame = getframe(gcf);
        writeVideo(makeVideo,frame);
    end
end

close(makeVideo);