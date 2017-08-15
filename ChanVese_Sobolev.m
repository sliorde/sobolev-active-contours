%% Marina Kokotov, 305891350
%% Lior Deutsch, 200700987
function ChanVese_Sobolev
	close all;
	
	% load and save parameters - see documentation 
	file_name = 'test.bmp';
	load_initial_contour = 0;
	initial_contour_file_name = 'pic1_init_snake.dat';
	record_movie = 0;
	movie_file_name = '';
	save_frames = 0;
	frames_file_names = '';
	
	% algorithm parameters - see documentation
	use_Sobolev = 1;
	lambda = 0.242;	
	steps = 1000;
 	dt = 0.045;
	out_size = 30;
 	min_dist = 1;
 	max_dist = 2;
	smooth_length = 0.095;	
	epsilon = 10^(-9);

	if (record_movie == 1)
		aviobj=avifile(movie_file_name,'fps',10,'compression','None');
	end
	
	% read image and convert to gray-double
	im = imread(file_name);
	if (size(im,3)==3)
		im = rgb2gray(im);
	end
	im = im2double(im);

	figure;
	imshow(im);

	% get initial contour
	if (load_initial_contour == 0)
		[BW, x, y] = roipoly;
		snake = [x, y];
		save(initial_contour_file_name,'snake','-ascii','-double','-tabs');
	else
		snake = load(initial_contour_file_name);
	end

	imshow(im);
	hold on;
	plot([snake(:,1);snake(1,1)],[snake(:,2);snake(1,2)],'y','LineWidth',2);
	drawnow;
	if (save_frames == 1)
		saveas(gcf,[frames_file_names, '_init.bmp'],'bmp');
	end
	hold off;	
	
	% fix initial snake by adding points and smoothing
	snake=addpoints(snake,min_dist,max_dist);
	if (snake(end, :) == snake(1, :))
		snake(end, :) = [];
	end
   	snake(:,1)=cconv(snake(:,1),(1/floor(1+size(snake,1)*smooth_length))*ones(floor(1+size(snake,1)*smooth_length),1),length(snake(:,1)));
   	snake(:,2)=cconv(snake(:,2),(1/floor(1+size(snake,1)*smooth_length))*ones(floor(1+size(snake,1)*smooth_length),1),length(snake(:,2)));

	imshow(im);
	hold on;
	plot([snake(:,1);snake(1,1)],[snake(:,2);snake(1,2)],'y','LineWidth',2);
	drawnow;
	if (save_frames == 1)
		saveas(gcf,[frames_file_names, '_smoothed.bmp'],'bmp');
	end
	hold off;	
	
	cont_flag = 1;
	k = 1;
	while k < steps
		BW = roipoly(im, snake(:,1), snake(:,2));
		
		if (strcmp(out_size,'all'))
			BW_out = ones(size(im))-BW;
		else
			BW_out = imdilate(BW,strel('disk',out_size))-BW;
		end
		
		mean_in = sum(sum(BW.*im))/sum(BW(:));
		mean_out = sum(sum(BW_out.*im))/sum(BW_out(:));

		% compute arclenth parameterization
		s = (snake([2:end,1],:)-snake);
		s = sqrt(s(:,1).^2+s(:,2).^2);
		L = sum(s);

		% find unit normal 
		normal = snake([2:end,1], 2) - snake([end,1:(end-1)], 2);
		normal(:, 2) = snake([end,1:(end-1)], 1) - snake([2:end,1], 1);
		normal = normal./repmat(sqrt(normal(:,1).^2+normal(:,2).^2),[1,2]);
		
		% make normal point inwards
		in_snake = round(snake + normal/1.5);
		in_snake(in_snake < 1) = 1;
		in_snake(in_snake(:,1)>size(im,2),1) = size(im,2);
		in_snake(in_snake(:,2)>size(im,1),2) = size(im,1);
		normal_orient = 2*BW(in_snake(:,2)+size(im,1)*(in_snake(:,1)-1))-1;		
		normal = normal.*(2*(sum(normal_orient) > 0)-1);
		
		% find normals that point outwards and delete corresponding points
		% (self-crossing elimination)
		bad_points = find(normal_orient == (1 - 2*(sum(normal_orient) > 0)));
		if ((length(bad_points) > 0) && (cont_flag == 1))		
			snake(bad_points,:)=[];
			snake = addpoints(snake,min_dist,max_dist);
			cont_flag = 1 - cont_flag;
			continue;
		end
		cont_flag = 1;
		
		% nearest neighbor interpolation of image
		snake_interp = round(snake);
		snake_interp(snake_interp < 1) = 1;
		snake_interp(snake_interp(:,1)>size(im,2),1) = size(im,2);
		snake_interp(snake_interp(:,2)>size(im,1),2) = size(im,1);

		% calculate L2 force
		F_L2 = (mean_in - mean_out).*(mean_in + mean_out - 2 * im(snake_interp(:,2)+size(im,1)*(snake_interp(:,1)-1)));
		F_L2 = normal.*repmat(F_L2,[1,2]);
		
		if (use_Sobolev == 0)
			F_L2 = L*F_L2/(sqrt(sum(F_L2(:,1).^2+F_L2(:,2).^2))+epsilon);
			snake = snake + F_L2*dt;
		else
			% calculate Sobolev force
			K = (1+(cumsum(s).^2-cumsum(s)*L+L^2/6)/(2*lambda*L^2))/L;
			F_Sobolev = cconv(F_L2(:,1).*s,K,length(F_L2(:,1)));
			F_Sobolev(:,2) = cconv(F_L2(:,2).*s,K,length(F_L2(:,2)));

			F_Sobolev = L*F_Sobolev/(sqrt(sum(F_Sobolev(:,1).^2+F_Sobolev(:,2).^2))+epsilon);
			snake = snake + F_Sobolev*dt;
		end
		
		% show snake in figure
		imshow(im);
		hold on;
 		plot([snake(:,1);snake(1,1)],[snake(:,2);snake(1,2)],'y','LineWidth',2);
		drawnow;
		
		if (save_frames == 1)
			saveas(gcf,[frames_file_names '_frame' num2str(k) '.bmp'],'bmp');
		end
		
		hold off;

		if (record_movie == 1)
			frame = getframe(gca);
			aviobj = addframe(aviobj,frame);
		end
		
		k = k+1;
	end

	if (record_movie == 1)
		aviobj = close(aviobj);
	end
end



function new_snake = addpoints(snake,min_dist,max_dist)
% add to snake points such that distance between succesive points is in
% range [min_dist, max_dist]
	x = [snake(:,1);snake(1,1)];
	y = [snake(:,2);snake(1,2)];
	dist = sqrt(diff(x).^2+diff(y).^2);
	x(dist<min_dist)=[];
	y(dist<min_dist)=[];

	ind=0;
	while ~isempty(ind)
		dist = sqrt(diff(x).^2+diff(y).^2);
		ind = find(dist>max_dist);
		ind = sort(ind,'descend')';
		for t = ind
			x = [x(1:t);(x(t)+x(t+1))/2; x(t+1:end)];
			y = [y(1:t);(y(t)+y(t+1))/2; y(t+1:end)];
		end
	end
	x = x(1:(end-1),1);
	y = y(1:(end-1),1);
	new_snake = [x,y];
end