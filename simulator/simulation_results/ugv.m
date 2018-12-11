classdef ugv
    methods (Static)        
        function [s01p, s12p, C01, C12, s0sp, s1sp] = sus_data_input(sus)
            s01p = sus(25:27);
            s12p = sus(28:30);
            C01 = [sus(31:33)';sus(34:36)'; sus(37:39)'];
            C12 = [sus(40:42)';sus(43:45)'; sus(46:48)'];
            s0sp = sus(49:51);
            s1sp = sus(52:54);
        end

        function [box_top, box_bot, box_right, box_left, box_front, box_back] = sus_display(top, bot, right, left, front, back)
            box_top = patch('XData',top(:,1),'YData',top(:,2),'ZData',top(:,3),'FaceColor',[0.5 0.5 0.5]);
            box_bot = patch('XData',bot(:,1),'YData',bot(:,2),'ZData',bot(:,3),'FaceColor',[0.5 0.5 0.5]);
            box_right = patch('XData',right(:,1),'YData',right(:,2),'ZData',right(:,3),'FaceColor',[0.5 0.5 0.5]);
            box_left = patch('XData',left(:,1),'YData',left(:,2),'ZData',left(:,3),'FaceColor',[0.5 0.5 0.5]);
            box_front = patch('XData',front(:,1),'YData',front(:,2),'ZData',front(:,3),'FaceColor',[0.5 0.5 0.5]);
            box_back = patch('XData',back(:,1),'YData',back(:,2),'ZData',back(:,3),'FaceColor',[0.5 0.5 0.5]);
        end

        function [RFT, LFT, LBT, RBT, RFD, LFD, LBD, RBD] = sus_geo(chassis_r0, chassis_A0, s01p, C01, ori,length_1, length_2, depth, width)
            RFT = chassis_r0 + chassis_A0*s01p + chassis_A0*C01*ori*[length_1; depth / 2; width];
            LFT = chassis_r0 + chassis_A0*s01p + chassis_A0*C01*ori*[length_1; depth / 2; 0];
            LBT = chassis_r0 + chassis_A0*s01p + chassis_A0*C01*ori*[-length_2; depth / 2; 0];
            RBT = chassis_r0 + chassis_A0*s01p + chassis_A0*C01*ori*[-length_2; depth / 2; width];
            RFD = chassis_r0 + chassis_A0*s01p + chassis_A0*C01*ori*[length_1; -depth / 2; width];
            LFD = chassis_r0 + chassis_A0*s01p + chassis_A0*C01*ori*[length_1; -depth / 2; 0];
            LBD = chassis_r0 + chassis_A0*s01p + chassis_A0*C01*ori*[-length_2; -depth / 2; 0];
            RBD = chassis_r0 + chassis_A0*s01p + chassis_A0*C01*ori*[-length_2; -depth / 2; width];
        end

        function [ori, wheel_ori] = sus_ori(q1, angle)
            ori = [cos(q1) sin(q1) 0; -sin(q1) cos(q1) 0; 0 0 1];
            wheel_ori = [cos(angle) sin(angle) 0; -sin(angle) cos(angle) 0; 0 0 1];
        end

        function [tire, wheel, wheel2] = tire_display(tire_vertices, tire_vertices2, tire_center, tire_center2)
            for j = 1 : size(tire_vertices,2)-1
                tire(j) = patch('XData',[tire_vertices(1,j:j+1)';flipud(tire_vertices2(1,j:j+1)')],...
                    'YData',[tire_vertices(2,j:j+1)';flipud(tire_vertices2(2,j:j+1)')],...
                    'ZData',[tire_vertices(3,j:j+1)';flipud(tire_vertices2(3,j:j+1)')],...
                    'FaceColor',[0.4,0.4,0.4],'EdgeColor','k');
                wheel(j) = patch('xdata',[tire_center(1),tire_vertices(1,j),tire_vertices(1,j+1)]',...
                    'ydata',[tire_center(2),tire_vertices(2,j),tire_vertices(2,j+1)]',...
                    'zdata',[tire_center(3),tire_vertices(3,j),tire_vertices(3,j+1)]',...
                'FaceColor',[0.4,0.4,0.4],'EdgeColor','k');
                wheel2(j) = patch('xdata',[tire_center2(1),tire_vertices2(1,j),tire_vertices2(1,j+1)]',...
                    'ydata',[tire_center2(2),tire_vertices2(2,j),tire_vertices2(2,j+1)]',...
                    'zdata',[tire_center2(3),tire_vertices2(3,j),tire_vertices2(3,j+1)]',...
                'FaceColor',[0.4,0.4,0.4],'EdgeColor','k');
                if j == 1
                    set(tire(j),'FaceColor','r')
                    set(wheel(j),'FaceColor','r')
                    set(wheel2(j),'FaceColor','r')
                end
            end
        end

        function [center, center2, vertices, vertices2] = tire_geo(r0, A0, s01p, C01, ori, wheel_ori, s12p, tire_offset, tire_width, tire_x, tire_y, tire_z)
            center = r0 + A0*s01p + A0*C01*ori*(s12p + [0;0;tire_offset]);
            center2 = r0 + A0*s01p + A0*C01*ori*(s12p + [0;0;tire_offset+tire_width]);
            for j = 1 : size(tire_x,2)
                vertices(:,j) = center + A0*C01*ori*wheel_ori*[tire_x(1,j);tire_y(1,j);tire_z(1,j)*tire_width];
                vertices2(:,j) = center + A0*C01*ori*wheel_ori*[tire_x(1,j);tire_y(1,j);tire_z(2,j)*tire_width];
            end
        end

        function tire_move(vertices, vertices2, center, center2, tire, wheel, wheel2)
            for j = 1 : size(vertices,2)-1
                set(tire(j), 'XData',[vertices(1,j:j+1)';flipud(vertices2(1,j:j+1)')],...
                    'YData',[vertices(2,j:j+1)';flipud(vertices2(2,j:j+1)')],...
                    'ZData',[vertices(3,j:j+1)';flipud(vertices2(3,j:j+1)')],...
                    'FaceColor',[0.4,0.4,0.4],'EdgeColor','k');
                set(wheel(j), 'xdata',[center(1),vertices(1,j),vertices(1,j+1)]',...
                    'ydata',[center(2),vertices(2,j),vertices(2,j+1)]',...
                    'zdata',[center(3),vertices(3,j),vertices(3,j+1)]',...
                'FaceColor',[0.4,0.4,0.4],'EdgeColor','k');
                set(wheel2(j), 'xdata',[center2(1),vertices2(1,j),vertices2(1,j+1)]',...
                    'ydata',[center2(2),vertices2(2,j),vertices2(2,j+1)]',...
                    'zdata',[center2(3),vertices2(3,j),vertices2(3,j+1)]',...
                'FaceColor',[0.4,0.4,0.4],'EdgeColor','k');
                if j == 1
                    set(tire(j),'FaceColor','r')
                    set(wheel(j),'FaceColor','r')
                    set(wheel2(j),'FaceColor','r')
                end
            end
        end

        function body_move(top, bot, right, left, front, back, box_top, box_bot, box_right, box_left, box_front, box_back)
            set(box_top,'XData',top(:,1)','YData',top(:,2)','ZData',top(:,3)')
            set(box_bot,'XData',bot(:,1)','YData',bot(:,2)','ZData',bot(:,3)')
            set(box_right,'XData',right(:,1)','YData',right(:,2)','ZData',right(:,3)')
            set(box_left,'XData',left(:,1)','YData',left(:,2)','ZData',left(:,3)')
            set(box_front,'XData',front(:,1)','YData',front(:,2)','ZData',front(:,3)')
            set(box_back,'XData',back(:,1)','YData',back(:,2)','ZData',back(:,3)')
        end
        
        function [top, bot, right, left, front, back] = body_vertices(RFT, LFT, LBT, RBT, RFD, LFD, LBD, RBD)
            top = [RFT';LFT';LBT';RBT'];
            bot = [RFD';LFD';LBD';RBD'];
            right = [RFT';RFD';RBD';RBT'];
            left = [LFT';LFD';LBD';LBT'];
            front = [RFT';RFD';LFD';LFT'];
            back = [RBT';RBD';LBD';LBT'];
        end
    end
end