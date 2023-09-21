function Data_norm = TimeNormalize(data)
%% input: 
    %data is a matrix consisting of the pattern trace that you want to time-normalize.
    %(i.e. joint angle/moment, contact force), and will normalize each collumn.
    nf = 101; %amount of frames that you want to time-normalize to
% output:
    % Data_norm: data normalized matrix. 
%%
for i = 1 : size(data,2)
temp = data(:,i); 
Data_norm(:,i) = spline([1:size(temp,1)]',temp,[1:((size(temp,1))/nf):size(temp,1)]');
end 