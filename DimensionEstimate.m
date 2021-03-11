function [bandwidth,dimension] = DimensionEstimate(distances)
%DimensionEstimate : Estimate the dimension of the underlying manifold and
%       optimal bandwidth parameter for kernel based learning algorithms
%   INPUTS:  distances  - N-by-Knn matrix of distances between each of N
%                         data points and the corresponding Knn nearest neighbors
%   OUTPUTS: bandwidth  - Kernel parameter, units of distance
%            dimension  - Estimate of the intrinsic dimension of the data
    
    bandwidths = mean(distances(:))*2.^(-20:.1:10);
    
    %%% Tune epsilon and estimate dimension
    D=zeros(1,length(bandwidths));
    for i=1:length(bandwidths)
        D(i) = sum(sum(exp(-distances.^2./(2*bandwidths(i)^2))));
    end

    [maxval,maxind] = max(diff(log(D))./diff(log(bandwidths)));
    
    dimension=maxval;
    bandwidth=bandwidths(maxind);
    
    figure;
    semilogx(bandwidths(2:end),diff(log(D))./diff(log(bandwidths)));
    hold on;
    semilogx(bandwidth,dimension,'ok');
    ylabel('Dimension','fontsize',22);
    xlabel('Bandwidth','fontsize',22);
    set(gca,'fontsize',18);
    
end

