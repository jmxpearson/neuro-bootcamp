%%
a = 5
a = 10;

vec = [1 , 2 , 18 , 36]
size(vec)
vec'
size(vec')

%%
A = [ [1 5 7] ; [2 18 10] ];
size(A)
B = [ [8 ; 9 ; 16] [7 ; -2 ; -3] ];
size(B)

% or

A = vertcat([1 5 7] , [2 18 10]);
B = horzcat([8 ; 9 ; 16] , [7 ; -2 ; -3]);

%% operations
5 * 2
5 - 3
7 / 4
8 + 7

2 * vec
vec + [ 1 1 7 8 ]
vec + [ 1 ; 1 ; 7 ; 8]
vec + [ 1 ; 1 ; 7 ; 8]'

vec * [ 1 1 7 8 ]
vec .* [ 1 1 7 8 ]
vec * [ 1 1 7 8 ]'

A * B
B * A

%% getting data from arrays
A
A(1, 2)
A(2, 2)
A(3)
A(5)

%% slicing/indexing
A(:, 1)
A(2, :)
A(:, :)
A(:)
A(:, end)
A(end)
A(end - 1, :)

%% sequences
1:10
1:2:10
10:-2:-5

B(1:2:end, :)