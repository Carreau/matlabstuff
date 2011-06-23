function m=wmean(data,weight)
	m=  sum(data .* weight / sum(weight));
end
