clear all
close all
clc

%% Geometric information (m)
% Center leg #2
L = 0.3;


%% Kinematics using product of exponential formula
syms q1 q2 q3 q1_1 q1_2 q1_3 q1_4 q1_5 

%% leg 1
        omega1_1 = [0;0;1];
        omega1_2 = [0;1;0];
        omega1_3 = [0;1;0];
        omega1_4 = [0;1;0];
        omega1_5 = [0;0;1];

        p1_1 = [0; 0; 0;];
        p1_2 = [0; 0; 0;];
        p1_3 = [0; 0; 2*L];
        p1_4 = [2*L; 0; 2*L];
        p1_5 = [2.1*L; 0; 2*L];

        S1_1 = [omega1_1;-cross(omega1_1,p1_1);];
        S1_2 = [omega1_2;-cross(omega1_2,p1_2);];
        S1_3 = [omega1_3;-cross(omega1_3,p1_3);];
        S1_4 = [omega1_4;-cross(omega1_4,p1_4);];
        S1_5 = [omega1_5;-cross(omega1_5,p1_5);];
        
        r1_1 = [0; 0; 0;];
        r1_2 = [0; 0; 0;];
        r1_3 = [0; 0; 2*L;];
        r1_4 = [2*L; 0; 2*L;];
        r1_5 = [2.1*L; 0; 2*L;];
        
        M1_01 = [[1,0,0;0,1,0;0,0,1],r1_1; 0,0,0,1;];
        M1_02 = [[1,0,0;0,0,-1;0,1,0],r1_2; 0,0,0,1;];
        M1_03 = [[1,0,0;0,0,-1;0,1,0],r1_3; 0,0,0,1;];
        M1_04 = [[1,0,0;0,0,-1;0,1,0],r1_4; 0,0,0,1;];
        M1_05 = [[1,0,0;0,1,0;0,0,1],r1_5; 0,0,0,1;];
        
        Exp1_1 = screw2matrix(S1_1,q1_1);
        Exp1_12 = screw2matrix(S1_1,q1_1)*screw2matrix(S1_2,q1_2);
        Exp1_123 = screw2matrix(S1_1,q1_1)*screw2matrix(S1_2,q1_2)*screw2matrix(S1_3,q1_3);
        Exp1_1234 = screw2matrix(S1_1,q1_1)*screw2matrix(S1_2,q1_2)*screw2matrix(S1_3,q1_3)*screw2matrix(S1_4,q1_4);
        Exp1_12345 = screw2matrix(S1_1,q1_1)*screw2matrix(S1_2,q1_2)*screw2matrix(S1_3,q1_3)*screw2matrix(S1_4,q1_4)*screw2matrix(S1_5,q1_5);

        Exp1_1 = simplify(Exp1_1);
        Exp1_12 = simplify(Exp1_12);
        Exp1_123 = simplify(Exp1_123);
        Exp1_1234 = simplify(Exp1_1234);
        Exp1_12345 = simplify(Exp1_12345);
        
        Exp1_45 = screw2matrix(S1_4,q1_4)*screw2matrix(S1_5,q1_5);
        Exp1_45 = simplify(Exp1_45);
        Mexp45 = Exp1_45(1:3,1:3);
        Inverse45 = inv(Mexp45);
        Zrequired = Inverse45*[1,0,0;0,1,0;0,0,1];
        Zrequired = simplify(Zrequired);
        
        T1_01 = Exp1_1*M1_01;
        T1_02 = Exp1_12*M1_02;
        T1_03 = Exp1_123*M1_03;
        T1_04 = Exp1_1234*M1_04;
        T1_05 = Exp1_12345*M1_05;

        T1_01 = simplify(T1_01);
        T1_02 = simplify(T1_02);
        T1_03 = simplify(T1_03);
        T1_04 = simplify(T1_04);
        T1_05 = simplify(T1_05);
        
        %% Symbolic function
        Q1_1 = symfun(T1_01(1:3,4),[q1_1,q1_2,q1_3,q1_4,q1_5]);
        Q1_2 = symfun(T1_02(1:3,4),[q1_1,q1_2,q1_3,q1_4,q1_5]);
        Q1_3 = symfun(T1_03(1:3,4),[q1_1,q1_2,q1_3,q1_4,q1_5]);
        Q1_4 = symfun(T1_04(1:3,4),[q1_1,q1_2,q1_3,q1_4,q1_5]);
        Q1_5 = symfun(T1_05(1:3,4),[q1_1,q1_2,q1_3,q1_4,q1_5]);
        
q1_1t = 0*pi/180;
q1_2t = 0*pi/180;
q1_3t = 0*pi/180;
q1_4t = 0*pi/180;
q1_5t = 0*pi/180;

Q1_1t = Q1_1(q1_1t,q1_2t,q1_3t,q1_4t,q1_5t);
Q1_2t = Q1_2(q1_1t,q1_2t,q1_3t,q1_4t,q1_5t);
Q1_3t = Q1_3(q1_1t,q1_2t,q1_3t,q1_4t,q1_5t);
Q1_4t = Q1_4(q1_1t,q1_2t,q1_3t,q1_4t,q1_5t);
Q1_5t = Q1_5(q1_1t,q1_2t,q1_3t,q1_4t,q1_5t);

figure(1)

%leg1
plot3([0,Q1_1t(1)], [0,Q1_1t(2)], [0,Q1_1t(3)], 'k','Linewidth',2)
hold on;
plot3([Q1_1t(1),Q1_2t(1)], [Q1_1t(2),Q1_2t(2)], [Q1_1t(3),Q1_2t(3)], 'b','Linewidth',2)
plot3([Q1_2t(1),Q1_3t(1)], [Q1_2t(2),Q1_3t(2)], [Q1_2t(3),Q1_3t(3)], 'r','Linewidth',2)
plot3([Q1_3t(1),Q1_4t(1)], [Q1_3t(2),Q1_4t(2)], [Q1_3t(3),Q1_4t(3)], 'g','Linewidth',2)
plot3([Q1_5t(1),Q1_4t(1)], [Q1_5t(2),Q1_4t(2)], [Q1_5t(3),Q1_4t(3)], 'm','Linewidth',2)


grid on 
axis equal 
xlabel('X(m)')
ylabel('Y(m)')
zlabel('Z(m)')
% axis([-0.8 0.8 -0.5 1 -0.8 1]);
view(145,20);


VV = VideoWriter('3Rtest1');
open(VV)
figure(1)

for i = 1 : 1 : 180
            
q1_1t = 10*pi/180;
q1_2t = i*pi/180;
q1_3t = 10*pi/180;
        

Q1_1t = Q1_1(q1_1t,q1_2t,q1_3t);
Q1_2t = Q1_2(q1_1t,q1_2t,q1_3t);
Q1_3t = Q1_3(q1_1t,q1_2t,q1_3t);
Q1_st = Q1_s(q1_1t,q1_2t,q1_3t);

%leg1
plot3([0,Q1_1t(1)], [0,Q1_1t(2)], [0,Q1_1t(3)], 'k','Linewidth',2)
hold on;
plot3([Q1_1t(1),Q1_2t(1)], [Q1_1t(2),Q1_2t(2)], [Q1_1t(3),Q1_2t(3)], 'b','Linewidth',2)
plot3([Q1_2t(1),Q1_3t(1)], [Q1_2t(2),Q1_3t(2)], [Q1_2t(3),Q1_3t(3)], 'r','Linewidth',2)
plot3([Q1_3t(1),Q1_st(1)], [Q1_3t(2),Q1_st(2)], [Q1_3t(3),Q1_st(3)], 'g','Linewidth',2)
plot3(Q1_st(1), Q1_st(2),Q1_st(3), 'ko')

grid on 
axis equal 
xlabel('X(m)')
ylabel('Y(m)')
zlabel('Z(m)')
axis([-1 1 -0.5 1 -0.8 1]);
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
