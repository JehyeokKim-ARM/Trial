clear all
close all
clc

%% Geometric information (m)
% Center leg #2
l2_1 = 0.25;
l2_2 = 0.2;
l2_3 = 0.4;
% Side leg #1
l1_1 = 0.2;
l1_2 = 0.2;
l1_3 = 0.38;
% Side leg #3
l3_1 = 0.2;
l3_2 = 0.2;
l3_3 = 0.38;

L = 0.05;
d = 0.05;

%% Kinematics using product of exponential formula
syms q1 q2 q3 q1_1 q1_2 q1_3 q3_1 q3_2 q3_3 
syms dq1 dq2 dq3 
syms ddq1 ddq2 ddq3
%% leg 1
        omega1_1 = [1;0;0];
        omega1_2 = [0;0;1];
        omega1_3 = [1;0;0];

        p1_1 = [0; 0; 0;];
        p1_2 = [L; 0; 0;];
        p1_3 = [L; 0; l1_1 + l1_2;];

        S1_1 = [omega1_1;-cross(omega1_1,p1_1);];
        S1_2 = [omega1_2;-cross(omega1_2,p1_2);];
        S1_3 = [omega1_3;-cross(omega1_3,p1_3);];

        r1_1 = [L; 0; 0;];
        r1_2 = [L; 0; l1_1;];
        r1_3 = [L; 0; l1_1+l1_2;];
        r1_s = [L; l1_3; l1_1+l1_2;];

        M1_01 = [[0,0,1;0,1,0;-1,0,0],r1_1; 0,0,0,1;];
        M1_02 = [[1,0,0;0,1,0;0,0,1],r1_2; 0,0,0,1;];
        M1_03 = [[0,0,1;0,1,0;-1,0,0],r1_3; 0,0,0,1;];
        M1_0s = [[1,0,0;0,1,0;0,0,1],r1_s; 0,0,0,1;];

        Exp1_1 = screw2matrix(S1_1,q1_1);
        Exp1_12 = screw2matrix(S1_1,q1_1)*screw2matrix(S1_2,q1_2);
        Exp1_123 = screw2matrix(S1_1,q1_1)*screw2matrix(S1_2,q1_2)*screw2matrix(S1_3,q1_3);


        Exp1_1 = simplify(Exp1_1);
        Exp1_12 = simplify(Exp1_12);
        Exp1_123 = simplify(Exp1_123);

        T1_01 = Exp1_1*M1_01;
        T1_02 = Exp1_12*M1_02;
        T1_03 = Exp1_123*M1_03;
        T1_0s = Exp1_123*M1_0s;

        T1_01 = simplify(T1_01);
        T1_02 = simplify(T1_02);
        T1_03 = simplify(T1_03);
        T1_0s = simplify(T1_0s);

        %% Symbolic function
        Q1_1 = symfun(T1_01(1:3,4),[q1_1,q1_2,q1_3]);
        Q1_2 = symfun(T1_02(1:3,4),[q1_1,q1_2,q1_3]);
        Q1_3 = symfun(T1_03(1:3,4),[q1_1,q1_2,q1_3]);
        Q1_s = symfun(T1_0s(1:3,4),[q1_1,q1_2,q1_3]);

        
%% leg 2
        omega2_1 = [1;0;0];
        omega2_2 = [0;0;1];
        omega2_3 = [1;0;0];

        p2_1 = [0; 0; 0;];
        p2_2 = [0; 0; 0;];
        p2_3 = [0; 0; l2_1 + l2_2;];

        S2_1 = [omega2_1;-cross(omega2_1,p2_1);];
        S2_2 = [omega2_2;-cross(omega2_2,p2_2);];
        S2_3 = [omega2_3;-cross(omega2_3,p2_3);];

        r2_1 = [0; 0; 0;];
        r2_2 = [0; 0; l2_1;];
        r2_3 = [0; 0; l2_1+l2_2;];
        r2_s = [0; l2_3; l2_1+l2_2;];

        M2_01 = [[0,0,1;0,1,0;-1,0,0],r2_1; 0,0,0,1;];
        M2_02 = [[1,0,0;0,1,0;0,0,1],r2_2; 0,0,0,1;];
        M2_03 = [[0,0,1;0,1,0;-1,0,0],r2_3; 0,0,0,1;];
        M2_0s = [[1,0,0;0,1,0;0,0,1],r2_s; 0,0,0,1;];

        Exp2_1 = screw2matrix(S2_1,q1);
        Exp2_12 = screw2matrix(S2_1,q1)*screw2matrix(S2_2,q2);
        Exp2_123 = screw2matrix(S2_1,q1)*screw2matrix(S2_2,q2)*screw2matrix(S2_3,q3);


        Exp2_1 = simplify(Exp2_1);
        Exp2_12 = simplify(Exp2_12);
        Exp2_123 = simplify(Exp2_123);

        T2_01 = Exp2_1*M2_01;
        T2_02 = Exp2_12*M2_02;
        T2_03 = Exp2_123*M2_03;
        T2_0s = Exp2_123*M2_0s;

        T2_01 = simplify(T2_01);
        T2_02 = simplify(T2_02);
        T2_03 = simplify(T2_03);
        T2_0s = simplify(T2_0s);

        %% Symbolic function
        Q2_1 = symfun(T2_01(1:3,4),[q1,q2,q3]);
        Q2_2 = symfun(T2_02(1:3,4),[q1,q2,q3]);
        Q2_3 = symfun(T2_03(1:3,4),[q1,q2,q3]);
        Q2_s = symfun(T2_0s(1:3,4),[q1,q2,q3]);

%% leg 3
        omega3_1 = [1;0;0];
        omega3_2 = [0;0;1];
        omega3_3 = [1;0;0];

        p3_1 = [0; 0; 0;];
        p3_2 = [-L; 0; 0;];
        p3_3 = [-L; 0; l1_1 + l1_2;];

        S3_1 = [omega3_1;-cross(omega3_1,p3_1);];
        S3_2 = [omega3_2;-cross(omega3_2,p3_2);];
        S3_3 = [omega3_3;-cross(omega3_3,p3_3);];

        r3_1 = [-L; 0; 0;];
        r3_2 = [-L; 0; l1_1;];
        r3_3 = [-L; 0; l1_1+l1_2;];
        r3_s = [-L; l1_3; l1_1+l1_2;];

        M3_01 = [[0,0,1;0,1,0;-1,0,0],r3_1; 0,0,0,1;];
        M3_02 = [[1,0,0;0,1,0;0,0,1],r3_2; 0,0,0,1;];
        M3_03 = [[0,0,1;0,1,0;-1,0,0],r3_3; 0,0,0,1;];
        M3_0s = [[1,0,0;0,1,0;0,0,1],r3_s; 0,0,0,1;];

        Exp3_1 = screw2matrix(S3_1,q3_1);
        Exp3_12 = screw2matrix(S3_1,q3_1)*screw2matrix(S3_2,q3_2);
        Exp3_123 = screw2matrix(S3_1,q3_1)*screw2matrix(S3_2,q3_2)*screw2matrix(S3_3,q3_3);


        Exp3_1 = simplify(Exp3_1);
        Exp3_12 = simplify(Exp3_12);
        Exp3_123 = simplify(Exp3_123);

        T3_01 = Exp3_1*M3_01;
        T3_02 = Exp3_12*M3_02;
        T3_03 = Exp3_123*M3_03;
        T3_0s = Exp3_123*M3_0s;

        T3_01 = simplify(T3_01);
        T3_02 = simplify(T3_02);
        T3_03 = simplify(T3_03);
        T3_0s = simplify(T3_0s);

        %% Symbolic function
        Q3_1 = symfun(T3_01(1:3,4),[q3_1,q3_2,q3_3]);
        Q3_2 = symfun(T3_02(1:3,4),[q3_1,q3_2,q3_3]);
        Q3_3 = symfun(T3_03(1:3,4),[q3_1,q3_2,q3_3]);
        Q3_s = symfun(T3_0s(1:3,4),[q3_1,q3_2,q3_3]);
        
%% Moving platform kinematics 
syms theta 
        omega_mp = [1/sqrt(3);1/sqrt(3);1/sqrt(3)];
        v = [0;0.4;0.3];
        S_mp = [omega_mp;-cross(omega_mp,v);];
        rmp_1 = v + [d*cos(pi/6); 0; -d*sin(pi/6);];
        rmp_2 = v + [0; 0; d;];
        rmp_3 = v + [-d*cos(pi/6); 0; -d*sin(pi/6);];
        
        RV_mp_1 = [[1,0,0;0,1,0;0,0,1],rmp_1; 0,0,0,1;];
        RV_mp_2 = [[1,0,0;0,1,0;0,0,1],rmp_2; 0,0,0,1;];
        RV_mp_3 = [[1,0,0;0,1,0;0,0,1],rmp_3; 0,0,0,1;];   
        
        gamma1 = 10*pi/180;
        gamma2 = 10*pi/180;
        gamma3 = 10*pi/180;
        
        d4 = 0.02;
        
        RRmp_1 = rmp_1 + [d4*cos(7*pi/6 + gamma1); 0; d4*sin(7*pi/6 + gamma1);];
        RRmp_2 = rmp_2 + [d4*cos(pi/2 + gamma2); 0; d4*sin(pi/2 + gamma2);];
        RRmp_3 = rmp_3 + [d4*cos(11*pi/6 + gamma3); 0; d4*sin(11*pi/6 + gamma3);];

        M_mp_1 = [[1,0,0;0,1,0;0,0,1],RRmp_1; 0,0,0,1;];
        M_mp_2 = [[1,0,0;0,1,0;0,0,1],RRmp_2; 0,0,0,1;];
        M_mp_3 = [[1,0,0;0,1,0;0,0,1],RRmp_3; 0,0,0,1;];

        Exp_mp = screw2matrix(S_mp,theta);
        Exp_mp = simplify(Exp_mp);
        
        % Sphrical joint 
        T_mp_1 = Exp_mp*M_mp_1;
        T_mp_2 = Exp_mp*M_mp_2;
        T_mp_3 = Exp_mp*M_mp_3;
        
        T_mp_1 = simplify(T_mp_1);
        T_mp_2 = simplify(T_mp_2);
        T_mp_3 = simplify(T_mp_3);
        
        Q_mp_1 = symfun(T_mp_1(1:3,4),theta);
        Q_mp_2 = symfun(T_mp_2(1:3,4),theta);
        Q_mp_3 = symfun(T_mp_3(1:3,4),theta);
        
        % Revolute joint at moving platform  
        RT_mp_1 = Exp_mp*RV_mp_1;
        RT_mp_2 = Exp_mp*RV_mp_2;
        RT_mp_3 = Exp_mp*RV_mp_3;
        
        RT_mp_1 = simplify(RT_mp_1);
        RT_mp_2 = simplify(RT_mp_2);
        RT_mp_3 = simplify(RT_mp_3);
        
        RQ_mp_1 = symfun(RT_mp_1(1:3,4),theta);
        RQ_mp_2 = symfun(RT_mp_2(1:3,4),theta);
        RQ_mp_3 = symfun(RT_mp_3(1:3,4),theta);       
        
thetat = 0*pi/180;
Q_mp_1t = Q_mp_1(thetat);
Q_mp_2t = Q_mp_2(thetat);
Q_mp_3t = Q_mp_3(thetat);
     
%% INVERSE KINEMATICS
    %leg 1
        Pend_1 = [Q_mp_1t(1), Q_mp_1t(2), Q_mp_1t(3)];
        Pend_star_1 = [L, Q_mp_1t(2), Q_mp_1t(3)];
        OPend_1 = distance3d(L,0,0,Pend_1(1),Pend_1(2),Pend_1(3));
        OPend_1star = distance3d(L,0,0,Pend_star_1(1),Pend_star_1(2),Pend_star_1(3));
        l1_3_star = sqrt(l1_3^2 - (L - Pend_1(1))^2);
        alpha_1 = acos(((l1_1+l1_2)^2 + OPend_1star^2 - l1_3_star^2)/(2*(l1_1+l1_2)*OPend_1star));

        %angle 3
            angle1_3 = -(pi/2 - acos(((l1_1+l1_2)^2 + l1_3^2 - OPend_1^2)/(2*(l1_1+l1_2)*l1_3)));
        %angle 1 & 2
        if Pend_1(2) >= 0
            angle1_1 = -((pi/2 - asin(Q_mp_1t(3)/OPend_1star)) - alpha_1);
            angle1_2 = -atan((Q_mp_1t(1)-L)/(OPend_1star*sin(alpha_1)));
        else
            angle1_1 = acos(Q_mp_1t(3)/OPend_1star)+alpha_1;
            angle1_2 = -atan((Q_mp_1t(1)-L)/(OPend_1star*sin(alpha_1)));
        end

	%leg 2
        Pend = [Q_mp_2t(1), Q_mp_2t(2), Q_mp_2t(3)];
        Pend_star = [0, Q_mp_2t(2), Q_mp_2t(3)];
        OPend = distance3d(0,0,0,Pend(1),Pend(2),Pend(3));
        OPend_star = distance3d(0,0,0,Pend_star(1),Pend_star(2),Pend_star(3));
        l3_2_star = sqrt(l2_3^2 - Pend(1)^2);
        alpha = acos(((l2_1+l2_2)^2 + OPend_star^2 - l3_2_star^2)/(2*(l2_1+l2_2)*OPend_star));

        %angle 3
            angle2_3 = -(pi/2 - acos(((l2_1+l2_2)^2 + l2_3^2 - OPend^2)/(2*(l2_1+l2_2)*l2_3)));
        %angle 1 & 2
        if Pend(2) >= 0
            angle2_1 = -((pi/2 - asin(Q_mp_2t(3)/OPend_star)) - alpha);
            angle2_2 = -atan(Q_mp_2t(1)/(OPend_star*sin(alpha)));
        else
            angle2_1 = acos(Q_mp_2t(3)/OPend_star)+alpha;
            angle2_2 = -atan(Q_mp_2t(1)/(OPend_star*sin(alpha)));
        end

    %leg 3
        Pend_3 = [Q_mp_3t(1), Q_mp_3t(2), Q_mp_3t(3)];
        Pend_star_3 = [-L, Q_mp_3t(2), Q_mp_3t(3)];
        OPend_3 = distance3d(-L,0,0,Pend_3(1),Pend_3(2),Pend_3(3));
        OPend_3star = distance3d(-L,0,0,Pend_star_3(1),Pend_star_3(2),Pend_star_3(3));
        l3_3_star = sqrt(l3_3^2 - (-L - Pend_3(1))^2);
        alpha_3 = acos(((l3_1+l3_2)^2 + OPend_3star^2 - l3_3_star^2)/(2*(l3_1+l3_2)*OPend_3star));

        %angle 3
            angle3_3 = -(pi/2 - acos(((l3_1+l3_2)^2 + l3_3^2 - OPend_3^2)/(2*(l3_1+l3_2)*l3_3)));
        %angle 1 & 2
        if Pend_3(2) >= 0
            angle3_1 = -((pi/2 - asin(Q_mp_3t(3)/OPend_3star)) - alpha_3);
            angle3_2 = -atan((Q_mp_3t(1)-(-L))/(OPend_3star*sin(alpha_3)));
        else
            angle3_1 = acos(Q_mp_3t(3)/OPend_3star)+alpha_3;
            angle3_2 = -atan((Q_mp_3t(1)-(-L))/(OPend_3star*sin(alpha_3)));
        end



q1_1t = angle1_1;
q1_2t = angle1_2;
q1_3t = angle1_3;

Q1_1t = Q1_1(q1_1t,q1_2t,q1_3t);
Q1_2t = Q1_2(q1_1t,q1_2t,q1_3t);
Q1_3t = Q1_3(q1_1t,q1_2t,q1_3t);
Q1_st = Q1_s(q1_1t,q1_2t,q1_3t);

q2_1t = angle2_1;
q2_2t = angle2_2;
q2_3t = angle2_3;

Q2_1t = Q2_1(q2_1t,q2_2t,q2_3t);
Q2_2t = Q2_2(q2_1t,q2_2t,q2_3t);
Q2_3t = Q2_3(q2_1t,q2_2t,q2_3t);
Q2_st = Q2_s(q2_1t,q2_2t,q2_3t);

q3_1t = angle3_1;
q3_2t = angle3_2;
q3_3t = angle3_3;

Q3_1t = Q3_1(q3_1t,q3_2t,q3_3t);
Q3_2t = Q3_2(q3_1t,q3_2t,q3_3t);
Q3_3t = Q3_3(q3_1t,q3_2t,q3_3t);
Q3_st = Q3_s(q3_1t,q3_2t,q3_3t);

figure(1)

%leg1
plot3([0,Q1_1t(1)], [0,Q1_1t(2)], [0,Q1_1t(3)], 'k','Linewidth',2)
hold on;
plot3([Q1_1t(1),Q1_2t(1)], [Q1_1t(2),Q1_2t(2)], [Q1_1t(3),Q1_2t(3)], 'b','Linewidth',2)
plot3([Q1_2t(1),Q1_3t(1)], [Q1_2t(2),Q1_3t(2)], [Q1_2t(3),Q1_3t(3)], 'r','Linewidth',2)
plot3([Q1_3t(1),Q1_st(1)], [Q1_3t(2),Q1_st(2)], [Q1_3t(3),Q1_st(3)], 'g','Linewidth',2)
plot3(Q1_st(1), Q1_st(2),Q1_st(3), 'ko')

%leg2
plot3([0,Q2_1t(1)], [0,Q2_1t(2)], [0,Q2_1t(3)], 'k','Linewidth',3)
plot3([Q2_1t(1),Q2_2t(1)], [Q2_1t(2),Q2_2t(2)], [Q2_1t(3),Q2_2t(3)], 'b','Linewidth',2)
plot3([Q2_2t(1),Q2_3t(1)], [Q2_2t(2),Q2_3t(2)], [Q2_2t(3),Q2_3t(3)], 'r','Linewidth',2)
plot3([Q2_3t(1),Q2_st(1)], [Q2_3t(2),Q2_st(2)], [Q2_3t(3),Q2_st(3)], 'g','Linewidth',2)
plot3(Q2_st(1), Q2_st(2),Q2_st(3), 'ko')

%leg3
plot3([0,Q3_1t(1)], [0,Q3_1t(2)], [0,Q3_1t(3)], 'k','Linewidth',3)
plot3([Q3_1t(1),Q3_2t(1)], [Q3_1t(2),Q3_2t(2)], [Q3_1t(3),Q3_2t(3)], 'b','Linewidth',2)
plot3([Q3_2t(1),Q3_3t(1)], [Q3_2t(2),Q3_3t(2)], [Q3_2t(3),Q3_3t(3)], 'r','Linewidth',2)
plot3([Q3_3t(1),Q3_st(1)], [Q3_3t(2),Q3_st(2)], [Q3_3t(3),Q3_st(3)], 'g','Linewidth',2)
plot3(Q3_st(1), Q3_st(2),Q3_st(3), 'ko')
%moving platform 
plot3([v(1),Q_mp_1t(1)], [v(2),Q_mp_1t(2)], [v(3),Q_mp_1t(3)], 'k','Linewidth',3)
plot3([v(1),Q_mp_2t(1)], [v(2),Q_mp_2t(2)], [v(3),Q_mp_2t(3)], 'k','Linewidth',3)
plot3([v(1),Q_mp_3t(1)], [v(2),Q_mp_3t(2)], [v(3),Q_mp_3t(3)], 'k','Linewidth',3)

grid on 
axis equal 
xlabel('X(m)')
ylabel('Y(m)')
zlabel('Z(m)')
axis([-0.8 0.8 -0.5 1 -0.2 1]);
view(145,20);


%% Profile generation (Workspace)

% Define operation conditions
Cycle_Time = 3;                   %% Cycle time [sec]
t_acc = 0.5;                      %% Acc time [sec]
t_dcc = 0.5;                      %% Dcc time [sec]
% v = [0;0.4;0.3];
% Xs1 = 0.45; Ys1 =  W1+0.38; Zs1 = H1+0.38;           %% Starting point of Path 1 [m]
% Xe1 = 0.45; Ye1 =  W1-0.38; Ze1 = H1+0.38;           %% Ending point of Path 1 [m]
% 
% Xs1 = 0.45; Ys1 =  W1-0.38; Zs1 = H1+0.38;           %% Starting point of Path 1 [m]
% Xe1 = 0.45; Ye1 =  W1-0.38; Ze1 = H1-0.38;           %% Ending point of Path 1 [m]
% 
% Xs1 = 0.45; Ys1 =  W1-0.38; Zs1 = H1-0.38;           %% Starting point of Path 1 [m]
% Xe1 = 0.45; Ye1 =  W1+0.38; Ze1 = H1-0.38;           %% Ending point of Path 1 [m]

Xs1 = -0; Ys1 =  0.3; Zs1 = 0.3;           %% Starting point of Path 1 [m]
Xe1 = 0; Ye1 =  0.3; Ze1 = -0.3;           %% Ending point of Path 1 [m]

% profile generation 
[t1,X1,dX1,ddX1,dddX1] = Profile_XYZ (Xs1,Xe1,Cycle_Time,t_acc);
[t1,Y1,dY1,ddY1,dddY1] = Profile_XYZ (Ys1,Ye1,Cycle_Time,t_acc);
[t1,Z1,dZ1,ddZ1,dddZ1] = Profile_XYZ (Zs1,Ze1,Cycle_Time,t_acc);
N1 = max(size(t1)); 

%% Profile generation (Joint space)

for i = 1 : 10 : N1
%moving platform
        omega_mp = [1/sqrt(3);1/sqrt(3);1/sqrt(3)];
        v = [X1(i);Y1(i);Z1(i)];
        S_mp = [omega_mp;-cross(omega_mp,v);];
        rmp_1 = v + [d*cos(pi/6); 0; -d*sin(pi/6);];
        rmp_2 = v + [0; 0; d;];
        rmp_3 = v + [-d*cos(pi/6); 0; -d*sin(pi/6);];

        M_mp_1 = [[1,0,0;0,1,0;0,0,1],rmp_1; 0,0,0,1;];
        M_mp_2 = [[1,0,0;0,1,0;0,0,1],rmp_2; 0,0,0,1;];
        M_mp_3 = [[1,0,0;0,1,0;0,0,1],rmp_3; 0,0,0,1;];

        Exp_mp = screw2matrix(S_mp,theta);
        Exp_mp = simplify(Exp_mp);
        
        T_mp_1 = Exp_mp*M_mp_1;
        T_mp_2 = Exp_mp*M_mp_2;
        T_mp_3 = Exp_mp*M_mp_3;
        
        T_mp_1 = simplify(T_mp_1);
        T_mp_2 = simplify(T_mp_2);
        T_mp_3 = simplify(T_mp_3);
        
        Q_mp_1 = symfun(T_mp_1(1:3,4),theta);
        Q_mp_2 = symfun(T_mp_2(1:3,4),theta);
        Q_mp_3 = symfun(T_mp_3(1:3,4),theta);
        
        thetat = 45*pi/180;
        Q_mp_1t(:,i) = Q_mp_1(thetat);
        Q_mp_2t(:,i) = Q_mp_2(thetat);
        Q_mp_3t(:,i) = Q_mp_3(thetat);
    
%leg1
        Pend_1 = [Q_mp_1t(1,i), Q_mp_1t(2,i), Q_mp_1t(3,i)];
        Pend_star_1 = [L, Q_mp_1t(2,i), Q_mp_1t(3,i)];
        OPend_1 = distance3d(L,0,0,Pend_1(1),Pend_1(2),Pend_1(3));
        OPend_1star = distance3d(L,0,0,Pend_star_1(1),Pend_star_1(2),Pend_star_1(3));
        l1_3_star = sqrt(l1_3^2 - (L - Pend_1(1))^2);
        alpha_1 = acos(((l1_1+l1_2)^2 + OPend_1star^2 - l1_3_star^2)/(2*(l1_1+l1_2)*OPend_1star));

        %angle 3
            angle1_3(i) = -(pi/2 - acos(((l1_1+l1_2)^2 + l1_3^2 - OPend_1^2)/(2*(l1_1+l1_2)*l1_3)));
        %angle 1 & 2
        if Pend_1(2) >= 0
            angle1_1(i) = -((pi/2 - asin(Q_mp_1t(3,i)/OPend_1star)) - alpha_1);
            angle1_2(i) = -atan((Q_mp_1t(1,i)-L)/(OPend_1star*sin(alpha_1)));
        else
            angle1_1(i) = acos(Q_mp_1t(3,i)/OPend_1star)+alpha_1;
            angle1_2(i) = -atan((Q_mp_1t(1,i)-L)/(OPend_1star*sin(alpha_1)));
        end


%leg2
        Pend = [Q_mp_2t(1,i), Q_mp_2t(2,i), Q_mp_2t(3,i)];
        Pend_star = [0, Q_mp_2t(2,i), Q_mp_2t(3,i)];
        OPend = distance3d(0,0,0,Pend(1),Pend(2),Pend(3));
        OPend_star = distance3d(0,0,0,Pend_star(1),Pend_star(2),Pend_star(3));
        l3_2_star = sqrt(l2_3^2 - Pend(1)^2);
        alpha = acos(((l2_1+l2_2)^2 + OPend_star^2 - l3_2_star^2)/(2*(l2_1+l2_2)*OPend_star));

        %angle 3
            angle2_3(i) = -(pi/2 - acos(((l2_1+l2_2)^2 + l2_3^2 - OPend^2)/(2*(l2_1+l2_2)*l2_3)));
        %angle 1 & 2
        if Pend(2) >= 0
            angle2_1(i) = -((pi/2 - asin(Q_mp_2t(3,i)/OPend_star)) - alpha);
            angle2_2(i) = -atan(Q_mp_2t(1,i)/(OPend_star*sin(alpha)));
        else
            angle2_1(i) = acos(Q_mp_2t(3,i)/OPend_star)+alpha;
            angle2_2(i) = -atan(Q_mp_2t(1,i)/(OPend_star*sin(alpha)));
        end
%leg3
        Pend_3 = [Q_mp_3t(1,i), Q_mp_3t(2,i), Q_mp_3t(3,i)];
        Pend_star_3 = [-L, Q_mp_3t(2,i), Q_mp_3t(3,i)];
        OPend_3 = distance3d(-L,0,0,Pend_3(1),Pend_3(2),Pend_3(3));
        OPend_3star = distance3d(-L,0,0,Pend_star_3(1),Pend_star_3(2),Pend_star_3(3));
        l3_3_star = sqrt(l3_3^2 - (-L - Pend_3(1))^2);
        alpha_3 = acos(((l3_1+l3_2)^2 + OPend_3star^2 - l3_3_star^2)/(2*(l3_1+l3_2)*OPend_3star));

        %angle 3
            angle3_3(i) = -(pi/2 - acos(((l3_1+l3_2)^2 + l3_3^2 - OPend_3^2)/(2*(l3_1+l3_2)*l3_3)));
        %angle 1 & 2
        if Pend_3(2) >= 0
            angle3_1(i) = -((pi/2 - asin(Q_mp_3t(3,i)/OPend_3star)) - alpha_3);
            angle3_2(i) = -atan((Q_mp_3t(1,i)-(-L))/(OPend_3star*sin(alpha_3)));
        else
            angle3_1(i) = acos(Q_mp_3t(3,i)/OPend_3star)+alpha_3;
            angle3_2(i) = -atan((Q_mp_3t(1,i)-(-L))/(OPend_3star*sin(alpha_3)));
        end
end

%     %% Right arm new profile 
%     for i = 1 : 1 : N1
%     M0EE = [[1,0,0;0,0,1;0,-1,0],[X1(i); Y1(i); Z1(i);]; 0,0,0,1;];
%     T0EE = Exp1*M0EE;
%     Q0EE = symfun(T0EE(1:3,4),[q1]); 
%     Q0EEt = Q0EE(-q1t);
%     P = eval(Q0EEt);
%     X1(i) = P(1); Y1(i) = P(2); Z1(i) = P(3);
%     end
%     plot3(X1,Y1,Z1);
%     axis equal
%     grid on 
% %% Profile generation (Joint space)
% 
% for i = 1 : 1 : N1
%     [q2t_rad(i), q3t_rad(i), q4t_rad(i)] = inv_ra(X1(i), Y1(i), Z1(i), L01, L1n, Ln2, L23, L34, L45,q1t); %%Inv_4Relbow is a fuction of inverse kinematics 
% end
% 
% for i = 2 : 1 : N1
%     dq2t_rad(i) = (q2t_rad(i) - q2t_rad(i-1)) / 0.01; 
%     dq3t_rad(i) = (q3t_rad(i) - q3t_rad(i-1)) / 0.01;
%     dq4t_rad(i) = (q4t_rad(i) - q4t_rad(i-1)) / 0.01;
%  end
% dq2t_rad(1) = 0; dq3t_rad(1) = 0; dq4t_rad(1) = 0; %% Initial condition of velocity
% 
% for i = 3 : 1 : N1
%     ddq2t_rad(i) = (dq2t_rad(i) - dq2t_rad(i-1)) / 0.01;
%     ddq3t_rad(i) = (dq3t_rad(i) - dq3t_rad(i-1)) / 0.01;
%     ddq4t_rad(i) = (dq4t_rad(i) - dq4t_rad(i-1)) / 0.01;     
%  end
% ddq2t_rad(1) = 0; ddq3t_rad(1) = 0; ddq4t_rad(1) = 0;
% ddq2t_rad(2) = 0; ddq3t_rad(2) = 0; ddq4t_rad(2) = 0;  %% Initial condition of acceleration 
% 
% q2t_deg = q2t_rad * 180 / pi(); q3t_deg = q3t_rad * 180 / pi(); q4t_deg = q4t_rad * 180 / pi(); 
% dq2t_deg = dq2t_rad * 180 / pi(); dq3t_deg = dq3t_rad * 180 / pi(); dq4t_deg = dq4t_rad * 180 / pi(); 
% ddq2t_deg = ddq2t_rad * 180 / pi(); ddq3t_deg = ddq3t_rad * 180 / pi(); ddq4t_deg = ddq4t_rad * 180 / pi(); 
% 
% figure(7)
% subplot(3,1,1)
% plot(t1, q2t_deg, 'r'); hold on; plot(t1,q3t_deg, 'b'); hold on; plot(t1,q4t_deg, 'k'); 
% title('Joint Space Profile (deg)')
% legend('q2','q3','q4' )
% xlabel('Time(sec)')
% ylabel('Angle(deg)')
% 
% subplot(3,1,2)
% plot(t1, dq2t_deg, 'r'); hold on; plot(t1,dq3t_deg, 'b'); hold on; plot(t1,dq4t_deg, 'k'); 
% title('Joint angular velocity Profile (rad)')
% legend('q2','q3','q4' )
% xlabel('Time(sec)')
% ylabel('Angular veolciy(deg/s)')
% 
% subplot(3,1,3)
% plot(t1, ddq2t_deg, 'r'); hold on; plot(t1,ddq3t_deg, 'b'); hold on; plot(t1,ddq4t_deg, 'k'); 
% title('Joint Angular acceleration (rad/s)')
% legend('q2','q3','q4' )
% xlabel('Time(sec)')
% ylabel('Angular acceleration (deg/s^2)')

% figure(8)
% subplot(3,1,1)
% plot(t1, q1t_rad, 'r'); hold on; plot(t1,q2t_rad, 'b'); hold on; plot(t1,q3t_rad, 'k'); 
% title('Joint angle Profile (rad)')
% legend('q1','q2','q3','q4' )
% xlabel('Time(sec)')
% ylabel('Angle(rad)')
% 
% subplot(3,1,2)
% plot(t1, dq1t_rad, 'r'); hold on; plot(t1,dq2t_rad, 'b'); hold on; plot(t1,dq3t_rad, 'k'); 
% title('Joint angular velocity Profile (rad)')
% legend('q1','q2','q3','q4' )
% xlabel('Time(sec)')
% ylabel('Angular veolciy (rad/s)')
% 
% subplot(3,1,3)
% plot(t1, ddq1t_rad, 'r'); hold on; plot(t1,ddq2t_rad, 'b'); hold on; plot(t1,ddq3t_rad, 'k'); 
% title('Joint angular acceleration Profile (rad)')
% legend('q1','q2','q3','q4' )
% xlabel('Time(sec)')
% ylabel('Angular acceleration (rad/s^2)')

%% Animation of movement 


VV = VideoWriter('Zmotion_arbi_0');
open(VV)
figure(1)
for i = 1 : 10 : N1

Q1_1t = Q1_1(angle1_1(i),angle1_2(i),angle1_3(i));
Q1_2t = Q1_2(angle1_1(i),angle1_2(i),angle1_3(i));
Q1_3t = Q1_3(angle1_1(i),angle1_2(i),angle1_3(i));
Q1_st = Q1_s(angle1_1(i),angle1_2(i),angle1_3(i));

Q2_1t = Q2_1(angle2_1(i),angle2_2(i),angle2_3(i));
Q2_2t = Q2_2(angle2_1(i),angle2_2(i),angle2_3(i));
Q2_3t = Q2_3(angle2_1(i),angle2_2(i),angle2_3(i));
Q2_st = Q2_s(angle2_1(i),angle2_2(i),angle2_3(i));

Q3_1t = Q3_1(angle3_1(i),angle3_2(i),angle3_3(i));
Q3_2t = Q3_2(angle3_1(i),angle3_2(i),angle3_3(i));
Q3_3t = Q3_3(angle3_1(i),angle3_2(i),angle3_3(i));
Q3_st = Q3_s(angle3_1(i),angle3_2(i),angle3_3(i));


%leg1
plot3([0,Q1_1t(1)], [0,Q1_1t(2)], [0,Q1_1t(3)], 'k','Linewidth',2)
hold on;
plot3([Q1_1t(1),Q1_2t(1)], [Q1_1t(2),Q1_2t(2)], [Q1_1t(3),Q1_2t(3)], 'b','Linewidth',2)
plot3([Q1_2t(1),Q1_3t(1)], [Q1_2t(2),Q1_3t(2)], [Q1_2t(3),Q1_3t(3)], 'r','Linewidth',2)
plot3([Q1_3t(1),Q1_st(1)], [Q1_3t(2),Q1_st(2)], [Q1_3t(3),Q1_st(3)], 'g','Linewidth',2)
plot3(Q1_st(1), Q1_st(2),Q1_st(3), 'ko')

%leg2
plot3([0,Q2_1t(1)], [0,Q2_1t(2)], [0,Q2_1t(3)], 'k','Linewidth',3)
plot3([Q2_1t(1),Q2_2t(1)], [Q2_1t(2),Q2_2t(2)], [Q2_1t(3),Q2_2t(3)], 'b','Linewidth',2)
plot3([Q2_2t(1),Q2_3t(1)], [Q2_2t(2),Q2_3t(2)], [Q2_2t(3),Q2_3t(3)], 'r','Linewidth',2)
plot3([Q2_3t(1),Q2_st(1)], [Q2_3t(2),Q2_st(2)], [Q2_3t(3),Q2_st(3)], 'g','Linewidth',2)
plot3(Q2_st(1), Q2_st(2),Q2_st(3), 'ko')

%leg3
plot3([0,Q3_1t(1)], [0,Q3_1t(2)], [0,Q3_1t(3)], 'k','Linewidth',3)
plot3([Q3_1t(1),Q3_2t(1)], [Q3_1t(2),Q3_2t(2)], [Q3_1t(3),Q3_2t(3)], 'b','Linewidth',2)
plot3([Q3_2t(1),Q3_3t(1)], [Q3_2t(2),Q3_3t(2)], [Q3_2t(3),Q3_3t(3)], 'r','Linewidth',2)
plot3([Q3_3t(1),Q3_st(1)], [Q3_3t(2),Q3_st(2)], [Q3_3t(3),Q3_st(3)], 'g','Linewidth',2)
plot3(Q3_st(1), Q3_st(2),Q3_st(3), 'ko')

%moving platform 
plot3([X1(i),Q_mp_1t(1,i)], [Y1(i),Q_mp_1t(2,i)], [Z1(i),Q_mp_1t(3,i)], 'k','Linewidth',3)
plot3([X1(i),Q_mp_2t(1,i)], [Y1(i),Q_mp_2t(2,i)], [Z1(i),Q_mp_2t(3,i)], 'k','Linewidth',3)
plot3([X1(i),Q_mp_3t(1,i)], [Y1(i),Q_mp_3t(2,i)], [Z1(i),Q_mp_3t(3,i)], 'k','Linewidth',3)
grid on 
axis equal 
xlabel('X(m)')
ylabel('Y(m)')
zlabel('Z(m)')
axis([-0.8 0.8 -0.5 1 -0.4 1]);
view(145,20);

hold off;



 frame = getframe(figure(1));
 writeVideo(VV,frame);
end
close(VV)


%% Plotting 1
% q1_1t = -20*pi/180;
% q1_2t = -20*pi/180;
% q1_3t = -20*pi/180;
% 
% q2_1t = 0*pi/180;
% q2_2t = 0*pi/180;
% q2_3t = 0*pi/180;
% 
% q3_1t = 0*pi/180;
% q3_2t = 0*pi/180;
% q3_3t = 0*pi/180;
% 
% 
% Q1_1t = Q1_1(q1_1t,q1_2t,q1_3t);
% Q1_2t = Q1_2(q1_1t,q1_2t,q1_3t);
% Q1_3t = Q1_3(q1_1t,q1_2t,q1_3t);
% Q1_st = Q1_s(q1_1t,q1_2t,q1_3t);
% 
% Q2_1t = Q2_1(q2_1t,q2_2t,q2_3t);
% Q2_2t = Q2_2(q2_1t,q2_2t,q2_3t);
% Q2_3t = Q2_3(q2_1t,q2_2t,q2_3t);
% Q2_st = Q2_s(q2_1t,q2_2t,q2_3t);
% 
% Q3_1t = Q3_1(q3_1t,q3_2t,q3_3t);
% Q3_2t = Q3_2(q3_1t,q3_2t,q3_3t);
% Q3_3t = Q3_3(q3_1t,q3_2t,q3_3t);
% Q3_st = Q3_s(q3_1t,q3_2t,q3_3t);
% 
% figure(1)
% %leg1
% plot3([0,Q1_1t(1)], [0,Q1_1t(2)], [0,Q1_1t(3)], 'k','Linewidth',2)
% hold on;
% plot3([Q1_1t(1),Q1_2t(1)], [Q1_1t(2),Q1_2t(2)], [Q1_1t(3),Q1_2t(3)], 'b','Linewidth',2)
% plot3([Q1_2t(1),Q1_3t(1)], [Q1_2t(2),Q1_3t(2)], [Q1_2t(3),Q1_3t(3)], 'r','Linewidth',2)
% plot3([Q1_3t(1),Q1_st(1)], [Q1_3t(2),Q1_st(2)], [Q1_3t(3),Q1_st(3)], 'g','Linewidth',2)
% plot3(Q1_st(1), Q1_st(2),Q1_st(3), 'ko')
% %leg2
% plot3([0,Q2_1t(1)], [0,Q2_1t(2)], [0,Q2_1t(3)], 'k','Linewidth',3)
% plot3([Q2_1t(1),Q2_2t(1)], [Q2_1t(2),Q2_2t(2)], [Q2_1t(3),Q2_2t(3)], 'b','Linewidth',2)
% plot3([Q2_2t(1),Q2_3t(1)], [Q2_2t(2),Q2_3t(2)], [Q2_2t(3),Q2_3t(3)], 'r','Linewidth',2)
% plot3([Q2_3t(1),Q2_st(1)], [Q2_3t(2),Q2_st(2)], [Q2_3t(3),Q2_st(3)], 'g','Linewidth',2)
% plot3(Q2_st(1), Q2_st(2),Q2_st(3), 'ko')
% %leg3
% plot3([0,Q3_1t(1)], [0,Q3_1t(2)], [0,Q3_1t(3)], 'k','Linewidth',3)
% plot3([Q3_1t(1),Q3_2t(1)], [Q3_1t(2),Q3_2t(2)], [Q3_1t(3),Q3_2t(3)], 'b','Linewidth',2)
% plot3([Q3_2t(1),Q3_3t(1)], [Q3_2t(2),Q3_3t(2)], [Q3_2t(3),Q3_3t(3)], 'r','Linewidth',2)
% plot3([Q3_3t(1),Q3_st(1)], [Q3_3t(2),Q3_st(2)], [Q3_3t(3),Q3_st(3)], 'g','Linewidth',2)
% plot3(Q3_st(1), Q3_st(2),Q3_st(3), 'ko')
% %moving platform 
% plot3([v(1),Q_mp_1t(1)], [v(2),Q_mp_1t(2)], [v(3),Q_mp_1t(3)], 'k','Linewidth',3)
% plot3([v(1),Q_mp_2t(1)], [v(2),Q_mp_2t(2)], [v(3),Q_mp_2t(3)], 'k','Linewidth',3)
% plot3([v(1),Q_mp_3t(1)], [v(2),Q_mp_3t(2)], [v(3),Q_mp_3t(3)], 'k','Linewidth',3)
% 
% grid on 
% axis equal 
% xlabel('X(m)')
% ylabel('Y(m)')
% zlabel('Z(m)')
% axis([-0.5 0.5 -0.5 1 -0.2 1]);
% view(145,20);









% 
% %% Dynamics using Lagrangian formulation 
% 
% %%%%%----- Inertia calculration -----%%%%%
% % Inertia of Link for body coordinate 
% Ixx23_2 = mL23*R23^2*1/4 + mL23*L23^2*1/12;
% Iyy23_2 = mL23*R23^2*1/2;
% Izz23_2 = mL23*R23^2*1/4 + mL23*L23^2*1/12;
% 
% Ixx34_3 = mL34*R34^2*1/4 + mL34*L34^2*1/12;
% Iyy34_3 = mL34*R34^2*1/2;
% Izz34_3 = mL34*R34^2*1/4 + mL34*L34^2*1/12;
% 
% Ixx45_4 = mL45*R45^2*1/4 + mL45*L45^2*1/12;
% Iyy45_4 = mL45*R45^2*1/2;
% Izz45_4 = mL45*R45^2*1/4 + mL45*L45^2*1/12;
% 
% for k = 1 : 1 : 3
%     for s = 1 : 1 : 3
%         R01(k,s) = T01(k,s);
%         R02(k,s) = T02(k,s);
%         R03(k,s) = T03(k,s);
%         R04(k,s) = T04(k,s);      
%     end
% end
% 
% % Inertia matrix for body coordinate 
% I23_2 = [Ixx23_2, 0, 0;
%          0, Iyy23_2, 0;
%          0, 0, Izz23_2];
% I34_3 = [Ixx34_3, 0, 0;
%          0, Iyy34_3, 0;
%          0, 0, Izz34_3];
% I45_4 = [Ixx45_4, 0, 0;
%          0, Iyy45_4, 0;
%          0, 0, Izz45_4];  
%  % Inertia matrix for fixed coordinate 
%  I23_0 = R02*I23_2*transpose(R02);
%  I34_0 = R03*I34_3*transpose(R03);
%  I45_0 = R04*I45_4*transpose(R04);
% 
%  %% Mass Matrix %%
%  %% Linear kinetic mass of links %%
%  
%  % Position of link's center for body coordinates
%  pLc23_2 = [0; -L23/2; 0; 1];
%  pLc34_3 = [0; -L34*43/100; 0; 1];
%  pLc45_4 = [0; -L45*43/100; 0; 1];
%  
%  % Position of link's center for fixed coordinates
%  pLc23_0 = T02*pLc23_2;
%  pLc34_0 = T03*pLc34_3;
%  pLc45_0 = T04*pLc45_4;
%  
%  % Jacobian for Link's linear motion 
%  tempJv23_2 = diff(pLc23_0, q2);
%  tempJv34_2 = diff(pLc34_0, q2);
%  tempJv34_3 = diff(pLc34_0, q3);
%  tempJv45_2 = diff(pLc45_0, q2);
%  tempJv45_3 = diff(pLc45_0, q3);
%  tempJv45_4 = diff(pLc45_0, q4);
% 
% Jv23 = [tempJv23_2(1,1),0,0;
%         tempJv23_2(2,1),0,0;
%         tempJv23_2(3,1),0,0];
% 
% Jv34 = [tempJv34_2(1,1),tempJv34_3(1,1),0;
%         tempJv34_2(2,1),tempJv34_3(2,1),0;
%         tempJv34_2(3,1),tempJv34_3(3,1),0];
% 
% Jv45 = [tempJv45_2(1,1),tempJv45_3(1,1),tempJv45_4(1,1);
%         tempJv45_2(2,1),tempJv45_3(2,1),tempJv45_4(2,1);
%         tempJv45_2(3,1),tempJv45_3(3,1),tempJv45_4(3,1)];
%  
% Kl_links = mL23*transpose(Jv23)*Jv23 + mL34*transpose(Jv34)*Jv34 + mL45*transpose(Jv45)*Jv45;
% 
% %% Linear kinetic mass of joints %%
%  % position of joints for body coordinates
%  pJc3_2 = [0; -L23; 0; 1];
%  pJc4_3 = [0; -L34; 0; 1];
%  pJc5_4 = [0; -L45; 0; 1];
%   
%  % position of joints for fixed coordinates
%  pJc3_0 = T02*pJc3_2;
%  pJc4_0 = T03*pJc4_3;
%  pJc5_0 = T04*pJc5_4;
%  
%  % Jacobian for joint's linear motion 
%  tempJvj3_2 = diff(pJc3_0, q2);
%  tempJvj4_2 = diff(pJc4_0, q2);
%  tempJvj4_3 = diff(pJc4_0, q3);
%  tempJvj5_2 = diff(pJc5_0, q2);
%  tempJvj5_3 = diff(pJc5_0, q3);
%  tempJvj5_4 = diff(pJc5_0, q4);
%  
%  Jvj3 = [tempJvj3_2(1,1),0,0;
%         tempJvj3_2(2,1),0,0;
%         tempJvj3_2(3,1),0,0];
% 
% Jvj4 = [tempJvj4_2(1,1),tempJvj4_3(1,1),0;
%         tempJvj4_2(2,1),tempJvj4_3(2,1),0;
%         tempJvj4_2(3,1),tempJvj4_3(3,1),0];
% 
% Jvj5 = [tempJvj5_2(1,1),tempJvj5_3(1,1),tempJvj5_4(1,1);
%         tempJvj5_2(2,1),tempJvj5_3(2,1),tempJvj5_4(2,1);
%         tempJvj5_2(3,1),tempJvj5_3(3,1),tempJvj5_4(3,1)];
% 
%  Kl_joints = mA3*transpose(Jvj3)*Jvj3 + mA4*transpose(Jvj4)*Jvj4 + (mEE+mPAY)*transpose(Jvj5)*Jvj5;
%  
% %% angular kinetic mass of links %%
% 
% omega2_0 = R02*omega2;
% omega3_0 = R03*omega3;
% omega4_0 = R04*omega4;
% 
% Jw23 = [[omega2_0(1:3,1)],zeros(3,2)];
% Jw34 = [zeros(3,1),[omega3_0(1:3,1)],zeros(3,1)];
% Jw45 = [zeros(3,2),[omega4_0(1:3,1)]];
% 
% Jw45 = Jw45 + Jw34 + Jw23;
% Jw34 = Jw34 + Jw23;
% 
% Ka_links = transpose(Jw23)*I23_0*Jw23 + transpose(Jw34)*I34_0*Jw34 + transpose(Jw45)*I45_0*Jw45;
% 
% %% Mass matrix 
% M = Kl_links + Kl_joints + Ka_links; 
% M = simplify(M);
% Mf = symfun(M,[q1,q2,q3,q4]);   
% T_mass = M*[ddq2; ddq3; ddq4];
%    
% %% Gravity matrix 
% 
% gm = [0;0;g];
% G =  -(transpose(Jv23)*mL23*gm + transpose(Jv34)*mL34*gm + transpose(Jv45)*mL45*gm + transpose(Jvj3)*mA3*gm + transpose(Jvj4)*mA4*gm + transpose(Jvj5)*(mEE+mPAY)*gm); % original 
% G = simplify(G);
% Gf = symfun(G,[q1,q2,q3,q4]);
% T_gravity = G;
% 
% %% C matrix (Centrifigual force) 
% 
% C(1,1) = (diff(M(1,1),q2) + diff(M(1,1),q2) - diff(M(1,1),q2))/2;
% C(2,1) = (diff(M(2,1),q2) + diff(M(2,1),q2) - diff(M(1,1),q3))/2;
% C(3,1) = (diff(M(3,1),q2) + diff(M(3,1),q2) - diff(M(1,1),q4))/2;
% 
% C(1,2) = (diff(M(1,2),q3) + diff(M(1,2),q3) - diff(M(2,2),q2))/2;
% C(2,2) = (diff(M(2,2),q3) + diff(M(2,2),q3) - diff(M(2,2),q3))/2;
% C(3,2) = (diff(M(3,2),q3) + diff(M(3,2),q3) - diff(M(2,2),q4))/2;
% 
% C(1,3) = (diff(M(1,3),q4) + diff(M(1,3),q4) - diff(M(3,3),q2))/2;
% C(2,3) = (diff(M(2,3),q4) + diff(M(2,3),q4) - diff(M(3,3),q3))/2;
% C(3,3) = (diff(M(3,3),q4) + diff(M(3,3),q4) - diff(M(3,3),q4))/2;
% 
% C = simplify(C);
% 
% T_centri = C*[dq2^2; dq3^2;dq4^2];
% 
% %% B matrix (Coliolis force)
% 
% B(1,1) = (diff(M(1,1),q3) + diff(M(1,2),q2)) - diff(M(1,2),q2);
% B(1,2) = (diff(M(1,1),q4) + diff(M(1,3),q2)) - diff(M(1,3),q2);
% B(1,3) = (diff(M(1,2),q4) + diff(M(1,3),q3)) - diff(M(2,3),q2);
% 
% B(2,1) = (diff(M(2,1),q3) + diff(M(2,2),q2)) - diff(M(1,2),q3);
% B(2,2) = (diff(M(2,1),q4) + diff(M(2,3),q2)) - diff(M(1,3),q3);
% B(2,3) = (diff(M(2,2),q4) + diff(M(2,3),q3)) - diff(M(2,3),q3);
% 
% B(3,1) = (diff(M(3,1),q3) + diff(M(3,2),q2)) - diff(M(1,2),q4);
% B(3,2) = (diff(M(3,1),q4) + diff(M(3,3),q2)) - diff(M(1,3),q4);
% B(3,3) = (diff(M(3,2),q4) + diff(M(3,3),q3)) - diff(M(3,3),q3);
% 
% B = simplify(B);
% 
% T_coli = B*[dq2*dq3; dq2*dq4; dq3*dq4];
% 
% %% B + C (centrifigual + coliolis)
% 
% T_BC = T_centri + T_coli;
% 
% %% M + B + C + G (Total torque)
% % T_total = T_mass + T_gravity + T_BC;
% T_total = T_gravity; % GRAVITY ONLY
% Tf_total = symfun(T_total, [q1,q2,q3,q4, dq2,dq3,dq4, ddq2,ddq3,ddq4]);
% 
% 
% %% Torque calculation [Nm] (All considered Mass, Centrifigual, Coliolis, Gravity)
% for i = 1 : 1 : N1
%     Total_torque(:,i) = Tf_total(q1t, q2t_rad(i), q3t_rad(i), q4t_rad(i), dq2t_rad(i), dq3t_rad(i), dq4t_rad(i), ddq2t_rad(i), ddq3t_rad(i), ddq4t_rad(i));
% end
% 
% Total_torque = simplify(Total_torque); 
% figure(2)
% plot(t1, Total_torque(1,:), 'r'); hold on; plot(t1,Total_torque(2,:), 'b'); hold on; plot(t1,Total_torque(3,:), 'k'); 
% 
% figure(100)
% plot(q4t_deg,Total_torque(3,:), 'b');
% 
% 
% 
% figure(101)
% plot(q3t_deg,Total_torque(2,:), 'b');
% grid on;
% 
% %% Animation of movement 
% 
% for i = 1 : 1 : N1
%     q3t_deg(i) = -90 + q3t_deg(i);
% end
% for i = 1 : 1 : N1
%     q4t_deg(i) = -90 + q4t_deg(i);
% end
% 
%  subplot(2,1,1)
%  plot(q3t_deg,Total_torque(2,:), 'b');
% % axis([-90 90 0 50]);
% grid on 
%  hold on;
%  subplot(2,1,2)
%  plot(q4t_deg,Total_torque(3,:), 'k');
%  legend('Torque4' )
%  xlabel('angle_bodyframe (deg)')
% ylabel('Torque (Nm)')
% title('Aassistance diagram joint 4');
%  
% v = VideoWriter('Inversetest4');
% open(v)
% % figure(200)
% for i = 1 : 10 : N1
%  figure(200)
%  subplot(2,1,1)
%  plot(q3t_deg(i),Total_torque(2,i), 'bo');
%  hold on;
%  subplot(2,1,2)
%  plot(q4t_deg(i),Total_torque(3,i), 'ko');
%  
%  frame = getframe(figure(1));
%  writeVideo(v,frame);
% end
% close(v)
% 
% 
