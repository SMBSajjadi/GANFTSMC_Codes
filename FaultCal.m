function D = FaultCal (num,U1,U2,U3,U4,phi,teta,sai,alpha,gama)

m=2;
l=0.2;
I_x=48*10^-4;
I_y=48*10^-4;
I_z=81*10^-4;
kt=4*10^-5;
kd=3*10^-6;

 A = [kt kt kt kt ; 0 -kt*l 0 kt*l ; -kt*l 0 kt*l 0 ; -kd kd -kd kd];
 B = [U1;U2;U3;U4];
 O = (A)^-1*B;

Rx = [1 0 0 ; 0 cos(phi) -sin(phi) ; 0 sin(phi) cos(phi) ];
Ry = [cos(teta) 0 sin(teta) ; 0 1 0 ; -sin(teta) 0 cos(teta) ];
Rz = [cos(sai) -sin(sai) 0 ; sin(sai) cos(sai) 0 ; 0 0 1 ];
r = Rz*Ry*Rx;

I = [1 0 0;0 1 0;0 0 1];
z = zeros(3,3);
R = [r z ; z I];

E = [kt 0 0 0 0 0 ; 0 kt 0 0 0 0 ; 0 0 kt kt kt kt ; -kd 0 0 -kt*l 0 kt*l ; 0 -kd -kt*l 0 kt*l 0 ; kt*l 0 -kd kd -kd kd];

x = 20*pi/180;
y = 25*pi/180;
% 
% x = 0.002;
% y = 0.003;

if num == 1
    F = [sin(gama)*cos(alpha) 0 0 0; sin(gama)*sin(alpha) 0 0 0 ; cos(gama) 0 0 0 ; 0 1 0 0 ; 0 0 1 0 ; 0 0 0 1 ];
else if num == 2
    F = [0 sin(gama)*cos(alpha) 0 0 ; 0 sin(gama)*sin(alpha) 0 0 ; 1 0 0 0 ; 0 cos(gama) 0 0 ; 0 0 1 0 ; 0 0 0 1 ];
    else if num == 3
    F = [0 0 sin(gama)*cos(alpha) 0 ; 0 0  sin(gama)*sin(alpha) 0 ; 1 0 0 0; 0 1 0 0 ; 0 0 cos(gama) 0 ; 0 0 0 1 ];
        else if num == 4
    F = [0 0 0 sin(gama)*cos(alpha) ; 0 0 0 sin(gama)*sin(alpha) ; 1 0 0 0; 0 1 0 0 ; 0 0 1 0 ; 0 0 0 cos(gama)  ];
            else if num == 12
    F = [sin(gama)*cos(alpha) sin(x)*cos(y) 0 0; sin(gama)*sin(alpha) sin(x)*sin(y) 0 0 ; cos(gama) 0 0 0 ; 0 cos(x) 0 0 ; 0 0 1 0 ; 0 0 0 1 ];
                else if num == 13 
    F = [sin(gama)*cos(alpha) 0 sin(x)*cos(y) 0 ; sin(gama)*sin(alpha) 0 sin(x)*sin(y) 0 ; cos(gama) 0 0 0 ; 0 1 0 0 ; 0 0 cos(x) 0 ; 0 0 0 1 ];
end
end
end
end
end
end




% F = [sin(gama)*cos(alpha) sin(x)*cos(y) 0 0; sin(gama)*sin(alpha) sin(x)*sin(y) 0 0 ; cos(gama) 0 0 0 ; 0 cos(x) 0 0 ; 0 0 1 0 ; 0 0 0 1 ];
%  F = [sin(gama)*cos(alpha) 0 sin(x)*cos(y) 0 ; sin(gama)*sin(alpha) 0 sin(x)*sin(y) 0 ; cos(gama) 0 0 0 ; 0 1 0 0 ; 0 0 cos(x) 0 ; 0 0 0 1 ];

prefault =  (R*E*F*O);


 Dyn = [(sin(phi)*sin(sai) + cos(phi)*cos(sai)*sin(teta))*U1 ; (-cos(sai)*sin(phi) + cos(phi)*sin(teta)*sin(sai))*U1 ; (cos(phi)*cos(teta))*U1 ; U2 ; U3 ; U4 ];
    
D =  [ (prefault(1)-Dyn(1))/m  (prefault(2)-Dyn(2))/m  (prefault(3)-Dyn(3))/m  (prefault(4)-Dyn(4))/I_x  (prefault(5)-Dyn(5))/I_y  (prefault(6)-Dyn(6))/I_z ] ;   


end