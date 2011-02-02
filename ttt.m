function ttt()
	files = liveviscosity();
	power = cellfun(@(x) x.Traps(5),files);
	rate = m_extractfield(files,'rate');
	avg = (1e6./rate-4) ./ 2

	etax = m_extractfield(files,'eta_x');
	etay = m_extractfield(files,'eta_y');
	figure(3);
	plot(avg,etax,'+');
	hold on ;
	plot(avg,etay,'ro');
	xlabel('avg');
	ylabel('\eta');
	hold off;
	figure(4);
	plot(rate,etax,'+');
	hold on ;
	plot(rate,etay,'ro');
	xlabel('rate');
	ylabel('\eta');
	hold off;
end
