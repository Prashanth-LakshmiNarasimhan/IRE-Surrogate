% ========================================================================
% Class definition with various metamodel settings and error computation
% ========================================================================

classdef Models
    methods(Static)
	
		% Function to build a metamodel for a given training data 
		% --------------------------------------------------------
		function [mdl, valerr, looerr, btime] = buildModel(Xtrain, Ytrain, Xval, Yval, mdlOpt)
			MetaOpts.Type = 'Metamodel';
			MetaOpts.ExpDesign.X = Xtrain;
			MetaOpts.ExpDesign.Y = Ytrain;
			MetaOpts.ValidationSet.X = Xval;
			MetaOpts.ValidationSet.Y = Yval;
			
			if contains(mdlOpt,'GPR') % Kriging/GPR
				MetaOpts.MetaType = 'Kriging';
				MetaOpts.EstimMethod  = 'CV';
                % MetaOpts.Optim.Display = 'iter';
				MetaOpts.Optim.Maxiter = 60;
            elseif contains(mdlOpt,'PCE') % Polynomial chaos expansion
				MetaOpts.MetaType = 'PCE';
				MetaOpts.TruncOptions.qNorm = 0.7:0.1:1;
                % MetaOpts.Display = 'verbose';
            elseif contains(mdlOpt,'PCK') % Polynomial chaos kriging
				MetaOpts.MetaType = 'PCK';
				MetaOpts.PCE.Method = 'LARS';
				MetaOpts.PCE.Degree = 3:15;
				MetaOpts.PCE.TruncOptions.qNorm = 0.7:0.1:1;
				MetaOpts.Kriging.Corr.Family = 'matern-5_2';
				MetaOpts.Optim.Maxiter = 60;
                % MetaOpts.Optim.Display = 'iter';
			end
			
			switch mdlOpt % Model choice
			
				case "GPR-linear" % Linear kernel
					MetaOpts.Corr.Family = 'linear';
				
				case "GPR-exponential" % Exponential kernel
					MetaOpts.Corr.Family = 'exponential';
				
				case "GPR-Gaussian" % Gaussian kernel
					MetaOpts.Corr.Family = 'Gaussian';
				
				case "GPR-Matern5/2" % Matern-5_2 kernel
					MetaOpts.Corr.Family = 'Matern-5_2';
				
				case "PCE-OLS" % Least-square calculation
					MetaOpts.Method = 'OLS';
					MetaOpts.Degree = 3:15;
						
				case "PCE-LARS" % Sparse Least-Angle-Regression-based
					MetaOpts.Method = 'LARS';
					MetaOpts.Degree = 3:15;

				case "PCE-OMP" % Orthogonal Matching Pursuit
					MetaOpts.Method = 'OMP';
					MetaOpts.Degree = 3:8;

				case "PCE-SP" % Subspace Pursuit
					MetaOpts.Method = 'SP';
					MetaOpts.Degree = 3:15;
					MetaOpts.SP.CVMethod = 'kfold';
				
				case "PCE-BCS" % Bayesian Compressive Sensing
					MetaOpts.Method = 'BCS';
					MetaOpts.Degree = 3:15;
				
				case "PCK-Seq" % Sequential polynomial chaos kriging
					MetaOpts.Mode = 'sequential';
				
				case "PCK-Opt" % Optimal polynomial chaos kriging - Sequential
					MetaOpts.Mode = 'optimal';  % slower but more accurate
				
				case "GPR-Cust"
					MetaOpts.Scaling = true;
					[Xtrainc, Ytrainc] = Models.SVMdata(Xtrain, Ytrain);
                    [~, Yvalc] = Models.SVMdata(Xval, Yval);
                    muX = mean(Xtrain);
                    StdX = std(Xtrain);
                    Xvalc = bsxfun(@rdivide, (bsxfun(@minus, Xval, muX)), StdX);
                    global SVC;
                    SVC = Models.buildModel(Xtrainc, Ytrainc, Xvalc, Yvalc, "SVMC");
					MetaOpts.Corr.Handle = @Models.my_eval_R;
					nOpars = 3;
					BoundsL = [ones(1,nOpars)*1e-2] ;
					BoundsU = [ones(1,nOpars)*10^3] ;
					MetaOpts.Optim.Bounds =[BoundsL ;BoundsU];
				
				case "SVMC" % Support vector machine classification
					MetaOpts.MetaType = 'SVC';
					MetaOpts.Penalization = 'linear';
					MetaOpts.QPSolver = 'SMO';
					MetaOpts.Kernel.Family = 'Gaussian';
					MetaOpts.Optim.Method = 'GS';
					MetaOpts.EstimMethod = 'SpanLOO';
				
				otherwise
					error("Invalid model choice");
			
			end
			
			tic;
			mdl = uq_createModel(MetaOpts);
			mdl.Internal.buildTime = toc;
			uq_print(mdl)
			valerr = mdl.Error.Val;
			looerr = mdl.Error.LOO;
			btime = mdl.Internal.buildTime;
            
            if mdlOpt == "GPR-Cust"
                mdl.Internal.SVM = SVC;
            end

		end
		
		
		% Custom kernel definition based on Gibbs kernel and SVM classifier 
		% ------------------------------------------------------------------
		function R = my_eval_R(x1,x2,theta, ~)

			% Initialize R matrix
			R = zeros(size(x1,1), size(x2,1));

			% Evaluating SVM classifier values
            % !!! If SVC is empty, kindly set the value of the global 
            % variable SVC using the following code. !!!
            %
            %   global SVC;
            %   SVC = GPRcustom.Internal.SVM;
            % 
            % !!! Run the above code before evaluating GPR custom model !!!
            
			global SVC; 
            
            if isempty(SVC)
                error('Set _GPRmodelVariable_.Internal.SVM as global SVC');
            end
            
			[~,d1] = uq_evalModel(SVC,x1);
			[~,d2] = uq_evalModel(SVC,x2);
			
			maxlength = 1e5; % Lengthscale value at AR
			nug = 1e-5; % Nugget term for the kernel
			
			% R matrix computation
			for i=1:size(x1,1)
				for j=1:size(x2,1)
				
					% Length scale value
                    if (d1(i) > 0) % x is in AR
                        l1 = maxlength;
                    else % x is in UR  ==>  l(x) = lambda_1 + lambda_2*abs(SVM(x)) + lambda_3*SVM(x)^2
                        l1 = theta(1) + theta(2) * abs(d1(i)) + theta(3) * (d1(i))^2 ;
                    end
                    
                    if (d2(j) > 0) % x' is in AR
                        l2 = maxlength;
                    else % x' is in UR  ==>  l(x') = lambda_1 + lambda_2*abs(SVM(x')) + lambda_3*SVM(x')^2
                        l2 = theta(1) + theta(2) * abs(d2(j)) + theta(3) * (d2(j))^2 ;
                    end
					
					% Compute the R values
					den = l1^2 + l2^2;
					fs = 2*l1*l2/den;                 
					ss=0;
					for k = 1:size(x1,2)
						ss = ss+((x1(i,k)-x2(j,k))^2);
					end
					R(i,j) = (fs^(k/2) * exp(-ss/den));
					
					if(i == j && size(x1,1) == size(x2,1))
						R(i,j) = R(i,j) + nug; % Adding nugget to diagonal terms
					end
				end
			end

        end
        
		% Preping data for SVM classifier 
		% -----------------------------------
		function [X, Y] = SVMdata(Xval,Yval)
			X = normalize(Xval);
			Y = (abs(Yval-100)<1e-8)*2-1;
		end
		
		% Postprocessing metamodel value to be bounded to phyiscal values 
		% ------------------------------------------------------------------
		function out = postpro(MModel, x)
            out = uq_evalModel(MModel,x);
            out(out > 100) = 100;
            out(out < 0) = 0;
        end
	
		% Root-mean-square error (RMSE) calculation for a metamodel at different regions
		% -------------------------------------------------------------------------------
		function errs = errorBins(MMs, X, Y, bins, isPostpro)
            if nargin < 5
                isPostpro = false;
            end
            
            nerrs = length(MMs);
            errs = ones(nerrs,length(bins)-1)*-1;
			
            for i=1:nerrs
                if(nerrs==1)
                    myMdl = MMs;
                else
                    myMdl = MMs{i};
                end
                
                for j = 1:length(bins)-1
                    sel = (Y<bins(j+1))&(Y>=bins(j));
                    if isPostpro
                        errs(i,j) = sqrt(mean((Y(sel) - Models.postpro(myMdl,X(sel,:))).^2));
                    else
                        errs(i,j) = sqrt(mean((Y(sel) - uq_evalModel(myMdl,X(sel,:))).^2));
                    end
                end
            end
            
        end
        
		% Root-mean-square error (RMSE) calculation for a metamodel
		% ----------------------------------------------------------
        function errs = RMSEcalc(MMs, X, Y, isPostpro)
            if nargin < 4
                isPostpro = false;
            end
            
            nerrs = length(MMs);
            errs = ones(nerrs,1)*-1;
			
            for i=1:nerrs
                if(nerrs==1)
                    myMdl = MMs;
                else
                    myMdl = MMs{i};
                end
				
                if isPostpro
                    errs(i) = sqrt(mean((Y - Models.postpro(myMdl,X)).^2));
                else
                    errs(i) = sqrt(mean((Y - uq_evalModel(myMdl,X)).^2));
                end
            end
        end
        
        % Function to calculate the lower bound of nugget parameter
		% ----------------------------------------------------------
        function nu = nuggetLB(GP)
            a = 20;
            lamb = max(eig(GP.Internal.Kriging.GP.R));
            C = cond(GP.Internal.Kriging.GP.R);
            nu = max([lamb*(C-exp(a))/C/(exp(a)-1),0]);
        end
	
    end
end
        