clear;close all;


%%% Generate data on a Hemisphere with an oscillating boundary
N = 50000;
N2 = floor(sqrt(N*4/pi));
[X,Y] = meshgrid(2/N2:2/N2:2);
X = [X(:)-1 Y(:)-1];
X = rand(size(X))*2-1;
thetas = atan2(X(:,2),X(:,1));
X = X(sqrt(sum(X.^2,2))<= (sin(6*(thetas-pi/12))/8 + 3/4),:);
X = X';
X(3,:) = sqrt(1-X(1,:).^2-X(2,:).^2);

%%% Analytic formula for the true sampling density, for details see 
%%% T. Berry, T. Sauer, Density Estimation on Manifolds with Boundary, 2016
maxr = .875;
M = 1/(2.81893*(1-maxr^2))/(128/73/pi);
X = X(:,1./(2.81893*(1-sum(X(1:2,:).^2)))/(M*128/73/pi) >= rand(1,size(X,2)));
N=size(X,2);
n=size(X,1);
r = sqrt(sum(X(1:2,:).^2));
TrueDensity = 1./sqrt(1-r.^2)/2.81893;




%%% Estimate the density (see function for input/output descriptions)
%%% Optionally, you can omit bandwidth/dimension and it will fit them for you

bandwidth=5e-3;
dimension = 2;
[density,densityCutting,densityHO] = BoundaryKDE(X,bandwidth,dimension);





%%% Compare the true density to the various estimators

cr = [.2 .75];

figure;
subplot(1,4,1);
scatter3(X(1,:),X(2,:),X(3,:),10,TrueDensity,'filled');caxis(cr);title('True Density')
subplot(1,4,2);
scatter3(X(1,:),X(2,:),X(3,:),10,density,'filled');caxis(cr);title('Consistent Estimator')
subplot(1,4,3);
scatter3(X(1,:),X(2,:),X(3,:),10,densityCutting,'filled');caxis(cr);title('Cut Estimator')
subplot(1,4,4);
scatter3(X(1,:),X(2,:),X(3,:),10,densityHO,'filled');caxis(cr);title('Higher Order Estimator')

figure;
subplot(1,4,1);
scatter(X(1,:),X(2,:),10,TrueDensity,'filled');caxis(cr);title('True Density')
subplot(1,4,2);
scatter(X(1,:),X(2,:),10,density,'filled');caxis(cr);title('Consistent Estimator')
subplot(1,4,3);
scatter(X(1,:),X(2,:),10,densityCutting,'filled');caxis(cr);title('Cut Estimator')
subplot(1,4,4);
scatter(X(1,:),X(2,:),10,densityHO,'filled');caxis(cr);title('Higher Order Estimator')


