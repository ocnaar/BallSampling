%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Octavio Narvaez-Aroche                                                  %
% ocnaar@berkeley.edu                                                     %
% Berkeley Center for Control and Identification                          %
% Fall 2017                                                               %
%                                                                         %
% Function to sample m points from an n-dimensional L-p ball of radius r  %
% with center at c.                                                       %
%                                                                         %
% Input                                                                   %
% m: number of samples.                                                   %
% r: ball radius in L-p.                                                  %
% c: n by 1 array with the center of the ball.                            %
% p: norm defining the L-p ball. Only p={1, 2, Inf} are supported.        %
% sd: seed for random number generation.                                  %
%                                                                         %
% Output                                                                  %
% X: n by m array of samples.                                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function X = nBallSampling(m,r,c,p,sd)

% Dimension of the L-p ball.
n = numel(c);

% Array for storing points within the n-dimensional L-p ball.
X = zeros(n,m);

% Control random number generation for Latin Hypercube Sampling (LHS).
rng(sd);

if n==1
    % Obtain m samples from a Latin Hypercube (LH).
    LH = lhsdesign(m,n);
    X = repmat(-r,[1,m])+2*r*LH';
else
    switch p
        case 1
            % Use LHS for obtaining random values between 0 and 1.
            mext = 10*m;
            LH = lhsdesign(mext,n);
            
            % Transform values from LHS into n-dimensional spherical 
			% coordinates.
            SC = repmat([r pi*ones(1,n-2) 2*pi],[mext,1]).*LH;
            
            % Compute trigonometric functions of angular coordinates.
            ss = [ones(mext,1), sin(SC(:,2:end))];
            cs = cos(SC(:,2:end));
            
            % Array for storing cartesian coordinates of the points 
			% within the n-dimensional L-2 ball.
            Xl2 = zeros(n,mext);
            
            % Rejection method to retain points inside the L-1 ball.  
            k = 1;
            for j=1:mext
                % Map spherical coordinates into cartesian coordinates.
                rj = SC(j,1);
                for i=1:n-1
                    Xl2(i,j) = rj*prod(ss(j,1:i))*cs(j,i);
                end
                Xl2(n,j) = rj*prod(ss(j,:));
				
                % Check 1-norm of sampled point. 
                l1 = norm(Xl2(:,j),1);
                if l1 <= r
                    % Retain point if inside the n-dimensional L-1 ball.
                    X(:,k) = Xl2(:,j);
                    k = k+1;
                end
                if k > m
                    fprintf('\nFinished obtaining points inside the L-1 ball with radius %d. \nThe ratio of accepted over rejected samples from the L-2 ball with the same radius is %d/%d.\n %.3d%% of the samples from the L-2 ball were rejected.\n',r,k-1,j-k-1,(j-k-1)/(j-2)*100);
                    break
                end
            end
            
        case 2
            % Use LHS for obtaining random values between 0 and 1.
            LH = lhsdesign(m,n);
            
            % Transform values from LHS into n-dimensional spherical 
			% coordinates.
            SC = repmat([r pi*ones(1,n-2) 2*pi],[m,1]).*LH;
            
            % Compute trigonometric functions of angular variables.
            ss = [ones(m,1), sin(SC(:,2:end))];
            cs = cos(SC(:,2:end));
            
            % Map spherical coordinates into cartesian coordinates.
            for j=1:m
                rj = SC(j,1);
                for i=1:n-1
                    X(i,j) = rj*prod(ss(j,1:i))*cs(j,i);
                end
                X(n,j) = rj*prod(ss(j,:));
            end
            
        case Inf
            % Obtain the m values from LHS.
            LH = lhsdesign(m,n);
            X = repmat(-r,[n,m])+2*r*LH';
            
        otherwise
            error('\nOnly L-p balls for p={1, 2, Inf} are supported.\n')
    end
end

% Translate the center of the ball from the origin to point c.
X = X + repmat(c,[1,m]);