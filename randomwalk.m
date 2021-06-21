function P = randomwalk(N, D, makeplot)
% randomwalk - create random walk in arbitrary dimension
%
%   P  = randomwalk(N, D) produces a N-by-D matrix P with the positions
%   of a random walk of N steps in D dimensions, starting from (0,..,0). 
%   P(K,:) holds the coordinates of the position at step K. Each step
%   from one position to the next is +1 or -1 along a random dimension.
%
%   With a third logical argument, randomwalk(N, D, true) creates a plot
%   of the random walk in a new figure.
%
%   Examples:
%      P = randomwalk(2, 10) % -> for example
%      % P = [0 0 ; 1 0 ; 1 1 ; 2 1 ; 2 0 ; 2 -1 ; 3 -1 ; 3 -2 ; 2 -2 ; 2 -1] 
%      randomwalk(1000,4,true) ; % create a plot
%
%   See also: RANDI, RANDPERM

% tested in Matlab 2017b
% version 1.0 (may 2018)
% (c) Jos van der Geest
% Matlab File Exchange Author ID: 10584
% email: samelinoa@gmail.com

% 1.0 (may 2018) - created after post on Matlab Answers

narginchk(2, 3) ;

if ~(isnumeric(N) && N==fix(N) && N>0)
    error('The number of steps, N, should be a positive integer.') ;
end
if ~(isnumeric(D) && D==fix(D) && D>0)
    error('The number of dimensions, D, should be a positive integer.') ;
end

% Basically the code can be written as a one-liner. I have split it here
% for clarity:
% 1) % function to get K random values of -1 and +1 (steps in either direction)
direction = 2*randi([0 1],1,N-1)-1 ; 
% 2) in which dimension should we move
whichDim = randi(D,1,N-1) ;
% 3) use sparse to create a step matrix: each row now specifies a single step
% we start at (0,0,...0)
P = full(sparse(2:N, whichDim, direction, N, D)) ;
% 4) use cumsum to transform steps into positions
P = cumsum(P) ;

if nargin==3  && ~isequal(makeplot, false)
    % visualisation
    figure ;
    hold on ;
    for k=1:D
        % plot each dimension as a separate line
        h = plot(1:N, P(:,k),'.-') ;
        text(N, P(N, k), sprintf(' Dim %d',k),'color',get(h,'color')) ;
    end
    xlabel('Step #') ;
    ylabel('Position') ;
    if D==1 
        title('Random walk in 1 dimension') ;
    else
        title(sprintf('Random walk in %d dimensions', D)) ;
    end
    hold off ;
end