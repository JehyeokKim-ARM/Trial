function Mindistance = mindistance3d(P0x,P0y,P0z,P1x,P1y,P1z, Q0x,Q0y,Q0z, Q1x,Q1y,Q1z)
u = [P1x,P1y,P1z] - [P0x,P0y,P0z];
v = [Q1x,Q1y,Q1z] - [Q0x,Q0y,Q0z];
W0 = [P0x,P0y,P0z] - [Q0x,Q0y,Q0z];

a = dot(u,u);
b = dot(u,v);
c = dot(v,v);
d = dot(u,W0);
e = dot(v,W0);

det = b^2  - a*c;
if det ==0
    sc = 0;
    tc = 0;
else
    sc = (b*e - c*d)/det;
    tc = (a*e - b*d)/det;
end
%% 현재는 시간상 약식으로 계산 
if sc >= 0 && sc < 0
    P = [P0x,P0y,P0z] + sc*u;
    Q = [Q0x,Q0y,Q0z] + tc*v;
    Mindistance = distance3d(P(1),P(2),P(3),Q(1),Q(2),Q(3));
elseif sc < 0
    Mindistance = distance3d(P0x,P0y,P0z,Q0x,Q0y,Q0z);
else 
    Mindistance = distance3d(P1x,P1y,P1z,Q1x,Q1y,Q1z);
end

end
