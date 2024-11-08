% Define symbolic variables for theta values
syms th1 th2 th3 th4 th5 th6 th7
syms d1 d3 d5 
syms L4 L7 
% Define DH parameters [a, d, alpha, theta]
dh_matrix = [0, d1, 0, th1;          % theta1
             0, 0, -pi/2, th2;          % theta2
             0, d3, pi/2, th3;          % theta3
             L4, 0, pi/2, th4;     % theta4
             -L4, d5, -pi/2, th5;   % theta5
             0, 0, pi/2, th6;           % theta6
             L7, 0, pi/2, th7];         % theta7

% DH transformation function with Craig's method
function A = dh_transform(a, d, alpha, theta)
    % Compute each transformation step according to Craig's method
    Rx_alpha = [1, 0, 0, 0;
                0, cos(alpha), -sin(alpha), 0;
                0, sin(alpha), cos(alpha), 0;
                0, 0, 0, 1];
                
    Dx_a = [1, 0, 0, a;
            0, 1, 0, 0;
            0, 0, 1, 0;
            0, 0, 0, 1];
            
    Rz_theta = [cos(theta), -sin(theta), 0, 0;
                sin(theta), cos(theta), 0, 0;
                0, 0, 1, 0;
                0, 0, 0, 1];
                
    Dz_d = [1, 0, 0, 0;
            0, 1, 0, 0;
            0, 0, 1, d;
            0, 0, 0, 1];
    
    % Combine transformations in the order specified by Craig's method
    A = Rx_alpha * Dx_a * Rz_theta * Dz_d;
end

% Compute individual transformation matrices A_i from frame i-1 to frame i
A_matrices = cell(1, 7); % Cell array to store each A_i matrix

for i = 1:7
    a = dh_matrix(i, 1);
    d = dh_matrix(i, 2);
    alpha = dh_matrix(i, 3);
    theta = dh_matrix(i, 4);
    A_matrices{i} = dh_transform(a, d, alpha, theta); % Compute A_i matrix
end

% Display each A_i transformation matrix (T_i^i-1)
for i = 1:7
    fprintf('Transformation Matrix A_%d^%d (from frame %d to frame %d):\n', i, i-1, i-1, i);
    disp(A_matrices{i});
    disp('-------------------------------------');
end

% Compute cumulative transformation matrices step by step (from base to frame i)
T_matrices = cell(1, 7); % Cell array to store cumulative T_0^i matrices

% Start by calculating T_2^0 = A_1 * A_2
T_matrices{2} = A_matrices{1} * A_matrices{2};  % Multiply first two matrices

% Now loop through and calculate the cumulative transformations for the remaining frames
for j = 3:7
    % Initialize T_j^0 with the first matrix A_1
    T_matrices{j} = A_matrices{1};
    
    % Multiply progressively from A_1 to A_j
    for k = 2:j
        T_matrices{j} = T_matrices{j} * A_matrices{k}; % Keep multiplying in sequence
    end
end

% Display cumulative transformation matrices T_2^0 to T_7^0
for i = 2:7
    fprintf('Cumulative Transformation Matrix T_%d^0 (from base to frame %d):\n', i, i);
    disp(T_matrices{i});
    disp('-------------------------------------');
end

% Compute the end-effector transformation directly as T_7^0 by multiplying all matrices
T_endeff = A_matrices{1} * A_matrices{2} * A_matrices{3} * A_matrices{4} * A_matrices{5} * A_matrices{6} * A_matrices{7};

% Display the end-effector transformation matrix T_7^0
fprintf('End-Effector Transformation Matrix T_7^0:\n');
disp(T_endeff);

% Check if the end-effector transformation matches the last T matrix
if simplify(T_endeff - T_matrices{7}) == sym(zeros(4))  % Compare the symbolic difference
    disp('Correct! The transformation matches.');
else
    disp('Fail! The transformation does not match.');
end
