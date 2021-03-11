function [density,densityCutting,densityHO,boundarydistances] = BoundaryKDE(data,bandwidth,dimension)
%BoundaryKDE : Estimate the dimension of the underlying manifold and
%       optimal bandwidth parameter for kernel based learning algorithms
%   INPUTS:  distances  - n-by-N matrix with N data points in R^n
%            kNN        - Number of nearest neighbors (saves memory)
%   OUTPUTS: density            - Consistent density estimator for manifolds with boundary
%            densityCutting     - Consistent density estimator based on the cutting method
%            densityHO          - Higher order density estimator based on Richardson extrapolation
%            boundarydistances  - distance to the boundary for each point (only
%                                  accurate for distances less than the bandwidth)

    [n,N] = size(data);
    
    density = zeros(N,1);
    densityCutting = zeros(N,1);
    densityHO = zeros(N,1);
    boundarydistances = zeros(N,1);
    
    if (nargin < 2)
        kNN = min(N,5*ceil(sqrt(N)));
        [distances,~] = pdist2(data',data','euclidean','smallest',kNN);
        [bandwidth,dimension] = DimensionEstimate(distances);
        bandwidth = (bandwidth)^2; %%% convert units to distance squared
    end

    m0  = (2*pi*bandwidth)^(dimension/2);
    
    for i=1:N
        ds = sum((data - repmat(data(:,i),1,N)).^2);
        weights = exp(-ds/bandwidth/2);
        vecs = (data - repmat(data(:,i),1,N)); %%% vectors pointing to neighboring data points

        %%% Standard Kernel Density Estimator (KDE) (not consistent at the boundary)
        qest = (1/m0)*mean(weights);
        
        %%% Boundary Direction Estimator (BDE)
        muvec = (2*pi)^(-(dimension-1)/2)*bandwidth^(-(dimension+1)/2)*(mean(repmat(weights,n,1).*vecs,2));

        mu = norm(muvec); %%% Norm of the boundary direction estimator
        muvec = muvec/mu; %%% normalize the BDE to get a unit vector in the direction of the boundary
        dotprods = sum(vecs.*repmat(muvec,1,N)); %%% dot product of each neighbor with the boundary direction

        
        %%% Use Newton's method to estimate the distance to the boundary
        %dMguess = sqrt(2*bandwidth/pi)*2*pi*abs(mu(i)-2*qest(i));
        dMguess = sqrt(max(-2*bandwidth*(log(mu/qest)-.4),0));
        contin=1;count = 0;
        while (contin)
            f = exp(dMguess^2/2/bandwidth)*(1+erf(dMguess/sqrt(2*bandwidth)))/2 - qest/mu;
            fp = 1/sqrt(2*pi*bandwidth) + exp(dMguess^2/2/bandwidth)*erf(dMguess/sqrt(2*bandwidth))*dMguess/2/bandwidth;
            delta = f/fp;
            if ((abs(delta) < bandwidth*1e-3)||(count>20)) contin=0; end
            dMguess = dMguess - delta;
            count = count+1;
        end
        bdist=max(min(dMguess,sqrt(max(ds))),0);
        boundarydistances(i) = bdist;
        
        
        %%% Consistent density estimator using Newton's estimate of distance to the boundary
        m0b = (1+erf(bdist/sqrt(2*bandwidth)))/2;
        density(i) = qest/m0b;

        
        %%% Cutting Method
        densityCutting(i) = (1/m0)*sum(weights(dotprods>-bdist))/N/m0b;
   
        
        %%% Introduce a second bandwidth for Richardson extrapoloation
        bandwidth2 = 2*bandwidth;
        weights2 = exp(-ds./bandwidth2/2);
        m02  = (2*pi*bandwidth2)^(dimension/2);
        m0b2 = (1+erf(bdist/sqrt(2*bandwidth2)))/2;
        qest21 = (1/m0)*sum(weights(dotprods>-bdist))/N;
        qest22 = (1/m02)*sum(weights2(dotprods>-bdist))/N;
        densityHO(i) = sqrt(2)*qest21*exp(bdist^2/2/bandwidth) - qest22*exp(bdist^2/2/bandwidth2);
        densityHO(i) = densityHO(i)/(sqrt(2)*m0b*exp(bdist^2/2/bandwidth)-m0b2*exp(bdist^2/2/bandwidth2));

    end


end

