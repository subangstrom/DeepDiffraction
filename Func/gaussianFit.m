function [P, PE] = gaussianFit( X, Y, opt )
%gaussianFit Gaussian fit
%   Makes a 1-D Gaussian fit of a distribution of points
%
%   [P,PE] = gaussianFit(X,Y,opt)
%
%   X is the x-axis input array (n x 1). 
%   Y is the y-axis (n x 1) input array. Y can be a (n x m) matrix.
%   In this case each column must be an array of points 
%   of the same length of X and the fit is done for each column. 
%   opt is the optional input argument:
%       1. Logical scalar (false/true). If false (default), no weights and
%       identity covariance matrix are used. If true, weights are
%       calculated using an algorithm that creates a gaussian distribution
%       around the peak with sigma equals to 1 r.m.s. of the distribution.
%       2. Vector of weights, must be a (n x 1) array
%       3. Covariance matrix, must be a (n x n) matrix.
%       
%   P is the output (3 x m) matrix of the fitted parametes (mu,s,A). 
%   PE is the output (3 x m) matrix of the errors on the fitted parameters.
%
%   Author: Nicola Galante 09/10/2015
%           <galante.nico@gmail.com>

switch nargin
    case 2
        if size(X,1)==1
            X = X';
        end

        if size(Y,1)==1
            Y = Y';
        end

        w = eye(length(X));
    case 3
        if size(X,1)==1
            X = X';
        end

        if size(Y,1)==1
            Y = Y';
        end

        if isscalar(opt)
            if islogical(opt)
                if opt==true
                    [pk,idx] = max(sum(Y,2));
                    m = X(idx);
                    s = sqrt(sum(sum(Y,2).*((X-m).^2))./sum(X.*sum(Y,2)));
                    w = pk*exp(-((m-X).^2)./(2.*s.^2));
                else
                    w = ones(size(X));
                end
            else
                error('localfunctions:gaussianFit:IllegalOption',...
                    'option must be either a logical statement or a numerical matrix');
            end
        elseif isvector(opt)
            if length(opt)==length(X)
                w = opt;
            else
                error('localfunctions:gaussianFit:IllegalVector',...
                    'weight array must be the same size as input X array');
            end
        elseif size(opt,1)==length(X) && size(opt,2)==length(X)
            w = opt;
        else
            error('localfunctions:gaussianFit:IllegalMatrix',...
                    'covariance matrix must be square and the same size as input X array');
        end
    otherwise
        error('localfunctions:gaussianFit:WrongNumberOfInputs',...
                    'function requires 2 or 3 arguments');
end


A = [X.*X, X, ones(size(X))];
[H,HE] = lscov(A,log(Y),w);
P = zeros(3,size(Y,2));
PE = zeros(3,size(Y,2));

P(1,:) = -H(2,:)./(2.*H(1,:));
P(2,:) = sqrt(-1./(2.*H(1,:)));
P(3,:) = exp(H(3,:)-H(2,:).^2./(4.*H(1,:)));

PE(1,:) = sqrt((H(2,:).^2 ./ H(1,:).^4) .* HE(1,:).^2 - (1 ./ (2.*H(1,:))) .* HE(2,:).^2);
PE(2,:) = sqrt(HE(1,:).^2 ./ H(1,:).^4);
PE(3,:) = sqrt(exp(H(3,:)+H(2,:).^2./2).^2 .* (H(2,:).^2 .* HE(2,:).^2 + HE(3,:).^2) );

end