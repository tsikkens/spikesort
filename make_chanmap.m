function chanmap = make_chanmap(shank1,nshanks)
% 
if nshanks == 1
    chanmap = shank1;
else
    chanmap = [make_chanmap(shank1,nshanks-1), nshanks*shank1];
end
% chanmap = zeros(numel(shank1),nshanks);
% for iShank = 1:nshanks
%    chanmap(:,iShank) = iShank*shank1;
% end