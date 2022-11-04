function KMN = K_PER(XM,XN,theta)

%   KMN = mykernel(XM,XN,theta) takes a M-by-D matrix XM, a N-by-D matrix
%   XN and computes a M-by-N matrix KMN of kernel products such that
%   KMN(i,j) is the kernel product between XM(i,:) and XN(j,:). theta is
%   the R-by-1 unconstrained parameter vector for the kernel.

% Periodic Kernel

D = size(XM,2);

params  = exp(theta);
sigmaL1 = params(1:D,1);
sigmaL2 = params(D+1,1);
sigmaF1 = params(D+2,1);
sigmaF2 = params(D+3,1);
sigmaP1 = params(D+4,1);

KMN = (sigmaF1^2)*exp(-2*sin(pi*(pdist2(XN,XM))/sigmaP1).^2)/(sigmaL1^2);