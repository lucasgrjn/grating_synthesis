% testing some tensor/gridding operations

clear; close all;

a = [ 1, 2 ];
b = [ 3, 4 ];
c = [ 5, 6, 7 ];

[A,B,C] = ndgrid( a, b, c );

n_abc = size(A);

A_vec = A(:);
B_vec = B(:);
C_vec = C(:);

% try reshaping to see if i can recover the original grids
A_reshaped = reshape(A, n_abc);
B_reshaped = reshape(B, n_abc);
C_reshaped = reshape(C, n_abc);

A
A_reshaped

B
B_reshaped

C
C_reshaped

% seems to work!