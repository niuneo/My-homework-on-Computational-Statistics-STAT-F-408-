
function X = randMCbinary(m,alfa,beta);

% ThreshLab/Tests/randMCbinary - binary Markov Chain
%  Usage
%    X = randMCbinary(m,alfa,beta);
%  Inputs
%    m    integer, size of output vector
%    alfa real number in [0,1]; alfa = P(X(n) = 1| X(n-1) = -1)
%    beta real number in [0,1]; beta = P(X(n) = -1| X(n-1) = 1)
%  Outputs
%    X    binary vector
%  Description
%  See also

U = rand(1,m);
X = zeros(size(U));
X(1) = (U(1)<alfa/(alfa+beta));
X(1) = 2*X(1)-1;
U(1) = Inf; % X(1) is fixed now, no changes needed
changes = zeros(size(U));
changeplus = find(U<alfa); 
changemin = find(U<beta);
changes(changeplus) = changes(changeplus)+1;
changes(changemin) = changes(changemin)-1;
changes(1) = X(1);
changelocations = [1 union(changeplus,changemin)];
changes = changes(changelocations);
% changes(i) = 1 means apply the following function:
% -1 -> 1; 1 -> 1
% changes(i) = -1 means apply the following function:
% -1 -> -1; 1 -> -1
% changes(i) = 0 means apply the following function:
% -1 -> 1; 1 -> -1
% So, after a "-1", we are certainly in state -1, whether the changes was real
% or not. That means that a 0 after -1 means actually a 1 and vice versa

while any(changes==0),
   z = find(changes==0); % note that z cannot contain 1
   changes(z) = -changes(z-1);
end

X(changelocations) = changes;
while any(X==0),
   z = find(X==0);
   X(z) = X(z-1);
end
