function params = read_parameter_file(param_file,out_dir)

% Copy param file to output directory
copyfile(param_file,[out_dir '/params.csv']);

% Read param file
P = readtable([out_dir '/params.csv'], ...
	'ReadVariableNames',false);
P.Properties.VariableNames = {'Var','Value'};

% Convert to struct
params = struct();
for p = 1:height(P)
	
	val = str2double(P.Value(p));
	
	% If it's not a number don't convert. If it is, do.
	if isnan(val)
		params.(P.Var{p}) = P.Value{p};
	else
		params.(P.Var{p}) = val;
	end
	
end


