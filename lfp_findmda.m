function [channels, varargout] = lfp_findmda(datafolder)
% Use as: channels = nlx_findncs(datafolder). Returns a 1xnChannels cell
% that contains the filenames (including extension) of all 
% Neuralynx .NCS files in the datafolder

%TS Last edit 23-01-2018

% Read filenames in datafolder
contents = dir(datafolder);

%Find Filenames that end with *.ncs
channels = cellfun(@(x) any(strfind(x, '.mda')),{contents.name});
channels = {contents(channels).name};

%Find channels numbers for sorting (windows sorts numbers as strings so 10
%is before 2 etc, this is not what we want)
channum = regexp(channels,'|\d+|[^.mda}','match');
channum(cellfun(@(x) isempty(x),channum)) = [];
channum = cellfun(@(x) sort(str2double(cell2mat(x))),channum);

%Place channels in correct order
[~,chansort] = sort(channum);
channels = channels(chansort);
channum = channum(chansort);
if nargout > 1
varargout{1} = channum;
end