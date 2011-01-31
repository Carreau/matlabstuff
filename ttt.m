function ttt()
	files = liveviscosity();
	power = cellfun(@(x) x.Traps(5),files);
	rate = m_extractfield(files,'rate');
	etax = m_extractfield(files,'eta_x');
	etay = m_extractfield(files,'eta_y');
	figure(3);
	plot(rate,etax,'+');
	hold on ;
	plot(rate,etay,'ro');
	xlabel('rate');
	ylabel('\eta');
end
