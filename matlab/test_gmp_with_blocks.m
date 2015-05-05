
%float
function test_gmp_with_blocks
	codebook_file = '/home/ntrang/project/output/hmdb51/feature/bow.codebook.devel/imprvdensetraj.hoghof/data/codebook.gmm.256.128.mat';
	load(codebook_file);
	
	low_proj_ = load('/home/ntrang/project/output/hmdb51/feature/bow.codebook.devel/imprvdensetraj.hoghof/data/lowproj.128.204.mat', 'low_proj');
	low_proj = low_proj_.low_proj;
	
	fisher_params.grad_weights = false;		% "soft" BOW
	fisher_params.grad_means = true;		% 1st order
	fisher_params.grad_variances = true;	% 2nd order
	fisher_params.alpha = single(1.0);		% power normalization (set to 1 to disable)
	fisher_params.pnorm = single(0.0);		% norm regularisation (set to 0 to disable)
	
	F = load('/home/ntrang/data/F.mat');
	F = F.X;

	tic;

	X = zeros(65536, size(F, 2));
	for i = 1:size(F, 2),
		cpp_handle = mexFisherEncodeHelperSP('init', codebook, fisher_params);
		mexFisherEncodeHelperSP('accumulate', cpp_handle, single(low_proj * F(:, i)));
		X(:, i) = mexFisherEncodeHelperSP('getfk', cpp_handle);
		mexFisherEncodeHelperSP('clear', cpp_handle);
		X(:, i) = X(:, i)/norm(X(:, i));
	end
	toc;
	
	tic;
	
	X_c = cell(1,256);
	for i = 1:size(F, 2),
		cpp_handle = mexFisherEncodeHelperSP('init', codebook, fisher_params);
		k = mexFisherEncodeHelperSP('accumulate', cpp_handle, single(low_proj * F(:, i))) + 1;
		code = mexFisherEncodeHelperSP('getfk', cpp_handle);
		mexFisherEncodeHelperSP('clear', cpp_handle);
		code = code/norm(code);
		X_c{k}(:, end+1) = code;
	end
	X_c = X_c(~cellfun(@isempty,X_c));
	%X_c = cell2mat(X_c);
	toc;
	
	tic;	
	X_c2 = cell(1,256);
	X2 = zeros(65536, size(F, 2));
	for i = 1:size(F, 2),
		cpp_handle = mexFisherEncodeHelperSP('init', codebook, fisher_params);
		[code, stats] = mexFisherEncodeHelperSP('encodestats', cpp_handle, single(low_proj * F(:, i)));
		mexFisherEncodeHelperSP('clear', cpp_handle);
		encode_entry = find(stats(2:257)>0);
		k = encode_entry(1);
		T = X_c2{k};
		s = size(T,2) + 1;
		T(:, s) = code;
		T(:, s) = T(:, s)/norm(T(:, s));
		X2(:, i) = T(:, s);
		X_c2{k} = T;
	end
	X_c2 = X_c2(~cellfun(@isempty,X_c2));
	%X_c2 = cell2mat(X_c2);
	toc;
	
	%if sum(~eq(X_c, X_c2)) == 0,
	%	fprintf('1-2 Equal!\n');
	%end
	
	tic;
	alpha = solve_gmp(10, X2.');
	code = X2 * alpha';
	toc;
	
	tic;
	alpha2 = solve_gmp_cell(10, X_c);	
	code2 = cell2mat(X_c) * alpha2.';
	toc;
	
	if sum(~eq(code, code2)) == 0,
		fprintf('1-2 Equal!\n');
	end
	
end

function w = solve_gmp(lambda, x)
	[N,~] = size(x);
	k	 = x*x';	
	nse = lambda*eye(N);
	e   = ones(N,1);
    b   = (k*e).^0;
	% solve gmp:
	w = (k + nse)\b;
	w = w';%/sum(w);
end

function w = solve_gmp_cell(lambda, x)
	k = cell(1,256);
	for i = 1:size(x,2),
		xi = x{i};
		k{i} = xi'*xi;
	end
	k = blkdiag(k{:});
	[N,~] = size(k);
	nse = lambda*eye(N);
	e   = ones(N,1);
	b   = (k*e).^0;
	% solve gmp:
	w = (k + nse)\b;
	w = w';%/sum(w);
end
