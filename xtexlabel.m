function xtexlabel(str)
	str = regexprep(str, '�', '\''e');
    str = regexprep(str, '�', '\`e');
    str = regexprep(str, '�', '\`a');
    str = regexprep(str, '(\\[A-Za-z]+)', '$$1$');
	set(0,'defaulttextinterpreter','none')
	xlabel(str)
end
