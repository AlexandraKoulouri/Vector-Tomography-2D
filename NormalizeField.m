function [nE] = NormalizeField(E)

amp = sqrt(E(:,1).^2 + E(:,2).^2);
nEx = E(:,1)./amp;
nEy = E(:,2)./amp;

nE = [nEx nEy];