function [t,x,dx,ddx,dddx] = Profile_XYZ (X1,X2,Cycle_Time,T_Acc)
index = 0;
%% 변수 설정 (사용자 설정)
q_init = X1;
q_final = X2;
vel_max = 5;  % 최대 속도 설정 : 1/2*jerk*TS^2 보다 커야 좋음.
TA = T_Acc;  % 가감속 시간
TS = T_Acc * 0.3;  % S-curve 시간 (2차 함수 형태)
T = Cycle_Time;  % 총 이동 시간
dt = 0.01; % 데이터 간격


%% 종속 변수 계산
if TA < TS*2
    disp('TA value error! Data is not reliable!')
end
dq = q_final - q_init;
vel = dq / ( T - TA );
acc = vel / ( TA - TS );
jerk = acc / TS;

T1 = TA - 2*TS; % 등가속구간 시간
T2 = T - 2*TA; % 등속구간 시간

%% 경로 생성

%% 사다리꼴 프로파일

if T>2*TA  % 사다리꼴 프로파일 조건
fprintf('Trapezoidal profile\n')    

% 경계 지점 값 정의 

t0 = 0;
v0 = 0;
q0 = q_init;

t1 = TS; % 등저크 구간이 끝나는 시간
v1 = (1/2)*jerk*TS^2;
q1 = (1/6)*jerk*TS^3 + q_init;

t2 = TS + T1; % 등가속 구간이 끝나는 시간
v2 = acc*T1 + v1;
q2 = (1/2)*acc*T1^2 + v1*T1 + q1;

t3 = TA; % 등저크 구간이 끝나는 시간
v3 = -(1/2)*jerk*TS^2 + acc*TS + v2;
v_max = v3;
q3 = -(1/6)*jerk*TS^3 + (1/2)*acc*TS^2 + v2*TS + q2;

t4 = TA + T2; % 등속도 구간이 끝나는 시간
v4 = v_max;
q4 = v_max*T2 + q3;

t5 = TA + T2 + TS; % 등저크 구간이 끝나는 시간
v5 = v2; % 대칭에 의해
q5 = q4 + (q3 - q2); % 대칭에 의해

t6 = TA + T2 + TS + T1; % 등가속 구간이 끝나는 시간
v6 = v1; % 대칭에 의해
q6 = q5 + (q2 - q1); % 대칭에 의해

t7 = T;
v7 = 0;
q7 = q_final;


for j = 0:floor(t1/dt) % floor : 내림
    index = index + 1;
    t = j*dt; % 기준시간 0.01초 단위
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

for j = floor(t6/dt)+1:round(t7/dt)-1  % 100*t7  (or 100*T) 가 정수로 안떨어지는 에러가 있어 round 로 씀.
    index = index + 1;
    t = j*dt - t6;
    q_data(index) = (1/6)*jerk*t^3 - (1/2)*acc*t^2 + v6*t + q6;
end

%% 삼각형 프로파일

elseif T <= 2*TA && T > 4*TS  % 삼각형 프로파일 조건 (이동 거리가 2*jerk*TS^3 보단 길어야 함 )
fprintf('Triangular profile\n')  

TA = T/2;
T1 = TA - 2*TS; % 등가속구간 시간
jerk = dq / (TS^2 + TS*T1) / TA;
acc = jerk*TS;

t1 = TS; % 등저크 구간이 끝나는 시간
v1 = (1/2)*jerk*TS^2;
q1 = (1/6)*jerk*TS^3 + q_init;

t2 = TS + T1; % 등가속 구간이 끝나는 시간
v2 = acc*T1 + v1;
q2 = (1/2)*acc*T1^2 + v1*T1 + q1;

t3 = TA; % 등저크 구간이 끝나는 시간
v3 = -(1/2)*jerk*TS^2 + acc*TS + v2;
v_max = v3;
q3 = -(1/6)*jerk*TS^3 + (1/2)*acc*TS^2 + v2*TS + q2;

% t4 등속 구간은 없음.

t5 = TA + TS; % 등저크 구간이 끝나는 시간
v5 = v2; % 대칭에 의해
q5 = q3 + (q3 - q2); % 대칭에 의해

t6 = TA + TS + T1; % 등가속 구간이 끝나는 시간
v6 = v1; % 대칭에 의해
q6 = q5 + (q2 - q1); % 대칭에 의해

t7 = T;
v7 = 0;
q7 = q_final;

for j = 0:floor(t1/dt) % floor : 내림
    index = index + 1;
    t = j*dt; % 기준시간 0.01초 단위
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

for j = floor(t6/dt)+1:round(t7/dt)-1  % 100*t7  (or 100*T) 가 정수로 안떨어지는 에러가 있어 round 로 씀.
    index = index + 1;
    t = j*dt - t6;
    q_data(index) = (1/6)*jerk*t^3 - (1/2)*acc*t^2 + v6*t + q6;
end

%% 종형 프로파일 (등저크 구간만 있는 경우)

else % 종형 프로파일 조건 (사다리꼴, 삼각형 프로파일이 아닌 경우)
fprintf('Bell profile\n')  
    
TS = T/4;
jerk = dq/(2*TS^3);
TS = (dq/(2*jerk))^(1/3);
acc = jerk*TS;

T = 4*TS;

t1 = TS; % 등저크 구간이 끝나는 시간
v1 = (1/2)*jerk*TS^2;
q1 = (1/6)*jerk*TS^3 + q_init;

% t2 등가속 구간은 없음.

t3 = 2*TS; % 등저크 구간이 끝나는 시간
v3 = -(1/2)*jerk*TS^2 + acc*TS + v1;
v_max = v3;
q3 = -(1/6)*jerk*TS^3 + (1/2)*acc*TS^2 + v1*TS + q1;

% t4 등속 구간은 없음.

t5 = 3*TS; % 등저크 구간이 끝나는 시간
v5 = v1; % 대칭에 의해
q5 = q3 + (q3 - q1); % 대칭에 의해

% t6 등가속 구간은 없음.

t7 = T;
v7 = 0;
q7 = q_final;

for j = 0:floor(t1/dt) % floor : 내림
    index = index + 1;
    t = j*dt; % 기준시간 0.01초 단위
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

for j = floor(t5/dt)+1:round(t7/dt)-1  % 100*t7  (or 100*T) 가 정수로 안떨어지는 에러가 있어 round 로 씀.
    index = index + 1;
    t = j*dt - t5;
    q_data(index) = (1/6)*jerk*t^3 - (1/2)*acc*t^2 + v5*t + q5;
end

end % 삼각형 or 사다리꼴 프로파일 if end

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
