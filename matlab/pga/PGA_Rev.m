function iM = PGA_Rev(M)
%PGA_REV Reversion of an eight-coefficient even PGA multivector.
% Ordered basis: [1,e23,e31,e12,e01,e02,e03,e0123].
% Reversion equals the inverse only for a suitably normalized versor.
% Adapted utility; SPDX-License-Identifier: Apache-2.0.
% Modified for the support sample: size validation and column output added;
% original reversion signs retained.
validateattributes(M,{'numeric'},{'vector','numel',8,'real','finite'});
iM=double(M(:));
iM(2:7)=-iM(2:7);
end
