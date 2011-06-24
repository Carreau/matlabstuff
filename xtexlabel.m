function xtexlabel(str)
	str = regexprep(str, 'é', '\''e');
    str = regexprep(str, 'è', '\`e');
    str = regexprep(str, 'à', '\`a');
    str = regexprep(str, '(\\[A-Za-z]+)', '$$1$');
	set(0,'defaulttextinterpreter','none')
	xlabel(str)
end
