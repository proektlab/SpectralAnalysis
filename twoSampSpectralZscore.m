function [dz,pdz] = twoSampSpectralZscore(sample1,dof1,sample2,dof2)
% function [dz,pdz] = twoSampSpectralZscore(sample1,dof1,sample2,dof2)

if any(iscomplex(sample1))||any(iscomplex(sample2))
    sample1 = abs(sample1);
    sample2 = abs(sample2);
end;

bias1=psi(dof1)-log(dof1); bias2=psi(dof2)-log(dof2); % bias from Thomson & Chave
var1=psi(1,dof1); var2=psi(1,dof2); % variance from Thomson & Chave
z1=log(sample1)-bias1; % Bias-corrected Fisher z, condition 1
z2=log(sample2)-bias2; % Bias-corrected Fisher z, condition 2
dz=(z1-z2)/sqrt(var1+var2); % z statistic
if nargout>1
    pdz=normpdf(dz,0,1); % probability of observing value dz
end