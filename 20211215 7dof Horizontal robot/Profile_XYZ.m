function [t,x,dx,ddx,dddx] = Profile_XYZ (X1,X2,Cycle_Time,T_Acc)
index = 0;
%% ���� ���� (����� ����)
q_init = X1;
q_final = X2;
vel_max = 5;  % �ִ� �ӵ� ���� : 1/2*jerk*TS^2 ���� Ŀ�� ����.
TA = T_Acc;  % ������ �ð�
TS = T_Acc * 0.3;  % S-curve �ð� (2�� �Լ� ����)
T = Cycle_Time;  % �� �̵� �ð�
dt = 0.01; % ������ ����


%% ���� ���� ���
if TA < TS*2
    disp('TA value error! Data is not reliable!')
end
dq = q_final - q_init;
vel = dq / ( T - TA );
acc = vel / ( TA - TS );
jerk = acc / TS;

T1 = TA - 2*TS; % ��ӱ��� �ð�
T2 = T - 2*TA; % ��ӱ��� �ð�

%% ��� ����

%% ��ٸ��� ��������

if T>2*TA  % ��ٸ��� �������� ����
fprintf('Trapezoidal profile\n')    

% ��� ���� �� ���� 

t0 = 0;
v0 = 0;
q0 = q_init;

t1 = TS; % ����ũ ������ ������ �ð�
v1 = (1/2)*jerk*TS^2;
q1 = (1/6)*jerk*TS^3 + q_init;

t2 = TS + T1; % ��� ������ ������ �ð�
v2 = acc*T1 + v1;
q2 = (1/2)*acc*T1^2 + v1*T1 + q1;

t3 = TA; % ����ũ ������ ������ �ð�
v3 = -(1/2)*jerk*TS^2 + acc*TS + v2;
v_max = v3;
q3 = -(1/6)*jerk*TS^3 + (1/2)*acc*TS^2 + v2*TS + q2;

t4 = TA + T2; % ��ӵ� ������ ������ �ð�
v4 = v_max;
q4 = v_max*T2 + q3;

t5 = TA + T2 + TS; % ����ũ ������ ������ �ð�
v5 = v2; % ��Ī�� ����
q5 = q4 + (q3 - q2); % ��Ī�� ����

t6 = TA + T2 + TS + T1; % ��� ������ ������ �ð�
v6 = v1; % ��Ī�� ����
q6 = q5 + (q2 - q1); % ��Ī�� ����

t7 = T;
v7 = 0;
q7 = q_final;


for j = 0:floor(t1/dt) % floor : ����
    index = index + 1;
    t = j*dt; % ���ؽð� 0.01�� ����
    q_data(index) = (1/6)*jerk*t^3 + q_init; 
end

for j = floor(t1/dt)+1:floor(t2/dt)
    index = index + 1;
    t = j*dt - t1;
    q_data(index) = (1/2)*acc*t^2 + v1*t + q1;
end

for j = floor(t2/dt)+1:floor(t3/dt)
    index = index + 1;
    t = j*dt - t2;
    q_data(index) = -(1/6)*jerk*t^3 + (1/2)*acc*t^2 + v2*t + q2;
end

for j = floor(t3/dt)+1:floor(t4/dt)
    index = index + 1;
    t = j*dt - t3;
    q_data(index) = v_max*t + q3;
end

for j = floor(t4/dt)+1:floor(t5/dt)
    index = index + 1;
    t = j*dt - t4;
    q_data(index) = -(1/6)*jerk*t^3 + v_max*t + q4;
end

for j = floor(t5/dt)+1:floor(t6/dt)
    index = index + 1;
    t = j*dt - t5;
    q_data(index) = -(1/2)*acc*t^2 + v5*t + q5;
end

for j = floor(t6/dt)+1:round(t7/dt)-1  % 100*t7  (or 100*T) �� ������ �ȶ������� ������ �־� round �� ��.
    index = index + 1;
    t = j*dt - t6;
    q_data(index) = (1/6)*jerk*t^3 - (1/2)*acc*t^2 + v6*t + q6;
end

%% �ﰢ�� ��������

elseif T <= 2*TA && T > 4*TS  % �ﰢ�� �������� ���� (�̵� �Ÿ��� 2*jerk*TS^3 ���� ���� �� )
fprintf('Triangular profile\n')  

TA = T/2;
T1 = TA - 2*TS; % ��ӱ��� �ð�
jerk = dq / (TS^2 + TS*T1) / TA;
acc = jerk*TS;

t1 = TS; % ����ũ ������ ������ �ð�
v1 = (1/2)*jerk*TS^2;
q1 = (1/6)*jerk*TS^3 + q_init;

t2 = TS + T1; % ��� ������ ������ �ð�
v2 = acc*T1 + v1;
q2 = (1/2)*acc*T1^2 + v1*T1 + q1;

t3 = TA; % ����ũ ������ ������ �ð�
v3 = -(1/2)*jerk*TS^2 + acc*TS + v2;
v_max = v3;
q3 = -(1/6)*jerk*TS^3 + (1/2)*acc*TS^2 + v2*TS + q2;

% t4 ��� ������ ����.

t5 = TA + TS; % ����ũ ������ ������ �ð�
v5 = v2; % ��Ī�� ����
q5 = q3 + (q3 - q2); % ��Ī�� ����

t6 = TA + TS + T1; % ��� ������ ������ �ð�
v6 = v1; % ��Ī�� ����
q6 = q5 + (q2 - q1); % ��Ī�� ����

t7 = T;
v7 = 0;
q7 = q_final;

for j = 0:floor(t1/dt) % floor : ����
    index = index + 1;
    t = j*dt; % ���ؽð� 0.01�� ����
    q_data(index) = (1/6)*jerk*t^3 + q_init; 
end

for j = floor(t1/dt)+1:floor(t2/dt)
    index = index + 1;
    t = j*dt - t1;
    q_data(index) = (1/2)*acc*t^2 + v1*t + q1;
end

for j = floor(t2/dt)+1:floor(t3/dt)
    index = index + 1;
    t = j*dt - t2;
    q_data(index) = -(1/6)*jerk*t^3 + (1/2)*acc*t^2 + v2*t + q2;
end

for j = floor(t3/dt)+1:floor(t5/dt)
    index = index + 1;
    t = j*dt - t3;
    q_data(index) = -(1/6)*jerk*t^3 + v_max*t + q3;
end

for j = floor(t5/dt)+1:floor(t6/dt)
    index = index + 1;
    t = j*dt - t5;
    q_data(index) = -(1/2)*acc*t^2 + v5*t + q5;
end

for j = floor(t6/dt)+1:round(t7/dt)-1  % 100*t7  (or 100*T) �� ������ �ȶ������� ������ �־� round �� ��.
    index = index + 1;
    t = j*dt - t6;
    q_data(index) = (1/6)*jerk*t^3 - (1/2)*acc*t^2 + v6*t + q6;
end

%% ���� �������� (����ũ ������ �ִ� ���)

else % ���� �������� ���� (��ٸ���, �ﰢ�� ���������� �ƴ� ���)
fprintf('Bell profile\n')  
    
TS = T/4;
jerk = dq/(2*TS^3);
TS = (dq/(2*jerk))^(1/3);
acc = jerk*TS;

T = 4*TS;

t1 = TS; % ����ũ ������ ������ �ð�
v1 = (1/2)*jerk*TS^2;
q1 = (1/6)*jerk*TS^3 + q_init;

% t2 ��� ������ ����.

t3 = 2*TS; % ����ũ ������ ������ �ð�
v3 = -(1/2)*jerk*TS^2 + acc*TS + v1;
v_max = v3;
q3 = -(1/6)*jerk*TS^3 + (1/2)*acc*TS^2 + v1*TS + q1;

% t4 ��� ������ ����.

t5 = 3*TS; % ����ũ ������ ������ �ð�
v5 = v1; % ��Ī�� ����
q5 = q3 + (q3 - q1); % ��Ī�� ����

% t6 ��� ������ ����.

t7 = T;
v7 = 0;
q7 = q_final;

for j = 0:floor(t1/dt) % floor : ����
    index = index + 1;
    t = j*dt; % ���ؽð� 0.01�� ����
    q_data(index) = (1/6)*jerk*t^3 + q_init; 
end

for j = floor(t1/dt)+1:floor(t3/dt)
    index = index + 1;
    t = j*dt - t1;
    q_data(index) = -(1/6)*jerk*t^3 + (1/2)*acc*t^2 + v1*t + q1;
end

for j = floor(t3/dt)+1:floor(t5/dt)
    index = index + 1;
    t = j*dt - t3;
    q_data(index) = -(1/6)*jerk*t^3 + v_max*t + q3;
end

for j = floor(t5/dt)+1:round(t7/dt)-1  % 100*t7  (or 100*T) �� ������ �ȶ������� ������ �־� round �� ��.
    index = index + 1;
    t = j*dt - t5;
    q_data(index) = (1/6)*jerk*t^3 - (1/2)*acc*t^2 + v5*t + q5;
end

end % �ﰢ�� or ��ٸ��� �������� if end

q_data = [q_init q_data]';
dq_data = diff(q_data)/dt;
dq_data = [0 ; dq_data];
ddq_data = diff(diff(q_data))/dt/dt;
ddq_data = [0 ; 0 ; ddq_data];
dddq_data = diff(diff(diff(q_data)))/dt/dt/dt;
dddq_data = [0 ; 0 ; 0 ; dddq_data];
% q_data plot
t_data = [0:dt:T];


t = t_data;
x = q_data;
dx = dq_data;
ddx = ddq_data;
dddx = dddq_data;


% figure(100)
% subplot(4,1,1)
% plot(t_data,q_data) 
% title('q data')
% xlim([0 T]);
% subplot(4,1,2)
% plot(t_data,dq_data) 
% title('dq data') 
% xlim([0 T]);
% subplot(4,1,3)
% plot(t_data,ddq_data) 
% title('ddq data')
% xlim([0 T]);
% subplot(4,1,4)
% plot(t_data,dddq_data) 
% title('dddq data')
% xlim([0 T]);
% 
% 
% if v_max > vel_max
%     fprintf('Velocity exceeded maximum velocity.\n',vel)
% end  
% fprintf('Max Velocity: %5.3f/s\n',v_max)
% fprintf('Max Acceleration: %5.3f/s\n',acc)
% fprintf('Jerk: %5.3f/s\n',jerk)
%         
