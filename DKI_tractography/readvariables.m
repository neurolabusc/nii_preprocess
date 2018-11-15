function S = readvariables(fn)
%READVARIABLES(fn) reads all variables defined in the text file fn and out
%puts a structure, S, with the variables af field names

%Read File
fid=fopen(fn); %EM
file = textscan(fid, '%s', 'delimiter', '\n', 'whitespace', ''); %EM
for i = 1:length(file{1}) %EM
    eval(file{1}{i})%EM
end
fclose(fid);

%Remove empty strings and comments
x = cellfun(@(x)x(regexpi(x,'\w|[%]','once'):min([regexp(x,'[%]')-1,length(x)])),file{1},'uniformoutput',0);
xi = cellfun(@(x)~isempty(x),x); %Index of non-empty / non-comment entries
C = cell(sum(xi),1);             %Cell of non-empty / non-comment entries
[C{:}] = deal(x{xi});

%Create output structure 
for i = 1:length(C)
    eval(['S.' C{i}]); 
end


