function Matrix = screw2matrix(S,theta)
omega = S(1:3);
v = S(4:6);
omega_hat = [0, -omega(3), omega(2);
             omega(3), 0 , -omega(1);
             -omega(2), omega(1), 0];
screw_hat = [omega_hat, v; 0, 0, 0, 0;];
Matrix = expm(screw_hat*theta);
Matrix = simplify(Matrix,'Steps', 50);