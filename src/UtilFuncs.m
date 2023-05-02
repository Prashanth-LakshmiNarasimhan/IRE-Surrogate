% ================================================================================================
% Class definition with functions used in sampling and plotting biophysical and metamodel results
% ================================================================================================

classdef UtilFuncs

    methods(Static)
        
		% Function to start the COMSOL livelink server and load  a COMSOL model
		% -----------------------------------------------------------------------
        function model = loadCOMSOLfile(modelname)
			format long
			addpath('C:\Program Files\COMSOL\COMSOL56\Multiphysics\mli')
			system("C:\Program Files\COMSOL\COMSOL56\Multiphysics\bin\win64\comsolmphserver.exe &")

			mphstart

			import com.comsol.model.*
			import com.comsol.model.util.*

			model = mphload(modelname);

		end

		% Function to end the COMSOL server (by closing the Command prompt in which the COMSOL server was started)
		% ---------------------------------------------------------------------------------------------------------
		function endComsolmphserver()
			system('Taskkill/IM cmd.exe');
		end
	
		% Function to run the COMSOL model for a set of input parameters
		% -----------------------------------------------------------------
		function QoI = solve_COMSOL_EPlist(model,params)
			nParamSize = size(params,1);

			% Empty containers for QoIs
			RelTA= zeros(1,nParamSize);
			TVol= zeros(1,nParamSize);
			TAbl= zeros(1,nParamSize);
			AZne= zeros(1,nParamSize);
			HtAbl = zeros(1,nParamSize);

			for i=1:nParamSize
				fprintf(['iter ', int2str(i), '  -->   ', num2str(params(i,:)), '\n'])
				try
					
					[RelTA(i), TVol(i), TAbl(i), AZne(i), HtAbl(i)] = UtilFuncs.solve_COMSOL_EP(model, params(i,:));
					
					fprintf(['Results  -->   ', num2str([RelTA(i), TVol(i), TAbl(i), AZne(i), HtAbl(i)]), '\n\n\n'])
					
				catch exception
					fprintf(["Exception. COMSOL model did not run \n"])
					continue
				end
			end

			QoI = [RelTA', TVol', TAbl', AZne', HtAbl'];
		end

		% Function to compute the QoI from a COMSOL model for a given input
		% ------------------------------------------------------------------
		function [RelTA, TVol, TAbl, AZne, HtAbl] = solve_COMSOL_EP(model,param)
			
			% update parameters in COMSOL model
			model.param.set('Tmr_E_Hlt', [num2str(param(1)) '[S/m]']);
			model.param.set('Lvr_E_Hlt', [num2str(param(2)) '[S/m]']);
			model.param.set('Tmr_E_Abl_Fac', num2str(param(3)));
			model.param.set('Lvr_E_Abl_Fac', num2str(param(4)));

			if(length(param)>4)
				model.param.set('V0', [num2str(param(5)) '[V]']);
			end

			mphgetexpressions(model.param)

			% update problem solution
			model.sol('sol28').runAll;
			model.result('pg32').run;
			model.result('pg33').run;
			model.result('pg34').run;
			model.result.evaluationGroup('eg1').run;

			% get solution
			RelTA = mphglobal(model,'comp1.RelTA');
			TVol = mphglobal(model,'comp1.TVol','unit','cm^3');
			TAbl = mphglobal(model,'comp1.TAbl','unit','cm^3');
			AZne = mphglobal(model,'comp1.AZne','unit','cm^3');
			HtAbl = mphglobal(model,'comp1.HtAbl','unit','cm^3');

		end
        
		% Function to change a given input parameter
		% -------------------------------------------
        function params = valReplace(params, ind, newvals)
            params(ind) = newvals;
        end
        
		% Function to define the input parameter distributions
		% ------------------------------------------------------
        function InputOpts = createInputs(noInps)
            
            % define the parameter pdf
            CDP = UtilFuncs.combinedDist(...
                [4.11, 4.25, 4.19, 4.30]/10, [2.56, 2.58, 2.35, 2.28]/10);
            InputOpts.Marginals(1).Name = 'TmrEHlt';
            InputOpts.Marginals(1).Type = 'Uniform';
            InputOpts.Marginals(1).Parameters = [CDP(1)-CDP(2) CDP(1)+CDP(2)];
            
            CDP = UtilFuncs.combinedDist(...
                [0.75, 0.74, 0.74, 0.92]/10, [0.28, 0.24, 0.20, 0.23]/10);
            InputOpts.Marginals(2).Name = 'LvrEHlt';
            InputOpts.Marginals(2).Type = 'Uniform' ;
            InputOpts.Marginals(2).Parameters = [CDP(1)-CDP(2) CDP(1)+CDP(2)];
            
            CDP = UtilFuncs.combinedDist(...
                [1.2826, 1.2826, 1.1923, 1.0517],...
                [0.25*1.2826, 0.25*1.2826, 0.25*1.1923, 0.25*1.0517]);
            InputOpts.Marginals(3).Name = 'TmrEAblFac';
            InputOpts.Marginals(3).Type = 'Uniform';
            InputOpts.Marginals(3).Parameters = [CDP(1)-CDP(2) CDP(1)+CDP(2)];
            
            CDP = UtilFuncs.combinedDist(...
                [3.3684, 3.3684, 3.5, 2.9583],...
                [0.25*3.3684, 0.25*3.3684, 0.25*3.5, 0.25*2.9583]);
            InputOpts.Marginals(4).Name = 'LvrEAblFac';
            InputOpts.Marginals(4).Type = 'Uniform' ;
            InputOpts.Marginals(4).Parameters = [CDP(1)-CDP(2) CDP(1)+CDP(2)];
            
            if(noInps > 4)
                InputOpts.Marginals(5).Name = 'V0';
                InputOpts.Marginals(5).Type = 'Uniform';
                InputOpts.Marginals(5).Parameters = [700 1500];
            end
        end
        
		% Function to combine multiple Normal distributions
		% --------------------------------------------------
        function out = combinedDist(mus, sigmas)
            noOfPoints = 10000;
            rng default % for reproducibility
            ys = zeros(1,length(mus)*noOfPoints);
            for i = 1:length(mus)
                ys(1,(i-1)*noOfPoints+1:(i)*noOfPoints) = normrnd(mus(i),sigmas(i),[1,noOfPoints]);
            end
            
            newpd = fitdist(ys','Normal');
            
            newMu = mean(newpd);
            newSD = std(newpd);
            
            out = [newMu newSD];
        end
		
		% Function to evaluate SVM classifier for a given input
		% ------------------------------------------------------
		function PV = getProbVal(mySVC, X)
            [~,PV] = uq_evalModel(mySVC, X);
        end
        
		% Function to plot a 3D slice of the QoI against two input parameters
		% ---------------------------------------------------------------------
        function plotMM2D(X, MM, paramList, paramInd, nsp)
            [Xm,Ym] = meshgrid(linspace(min(X(:,paramInd(1))),max(X(:,paramInd(1))),nsp),linspace(min(X(:,paramInd(2))),max(X(:,paramInd(2))),nsp));
            
            Zm = Xm;
            for i = 1:size(Xm,1)
                for j = 1:size(Xm,2)
                    Zm(i,j) = uq_evalModel(MM, UtilFuncs.valReplace(paramList, paramInd, [Xm(i,j),Ym(i,j)]));
                end
            end
            
            surf(Xm,Ym,Zm)
        end
        
		% Function to visualize the accuracy of model prediction against validation data
		% -------------------------------------------------------------------------------
        function plotValSet(Ytrue,YMM, titlename, Xlbl, Ylbl)
            uq_plot(Ytrue, YMM, '+')
            hold on
            uq_plot([min([Ytrue,YMM]) max([Ytrue,YMM])],...
                [min([Ytrue,YMM]) max([Ytrue,YMM])], 'k')
            hold off
            axis([min(Ytrue) max(Ytrue) min(YMM) max(YMM)*1.001])
            
            title(titlename)
            xlabel(Xlbl)
            ylabel(Ylbl)
        end
        
		
		% Function to plot metamodel prediction for given data (with variance and postprocessing)
		% ----------------------------------------------------------------------------------------
        function plotMM1D(varargin)
            Defaults = {[],[],[],[1],1,['X','Y'],0,0,0};
            idx = ~cellfun('isempty',varargin);
            Defaults(idx) = varargin(idx);
            
            X = Defaults{1};
            MM = Defaults{2};
            paramList = Defaults{3};
            paramInd = Defaults{4};
            nOut = Defaults{5};
            labels = Defaults{6};
            isHold = Defaults{7};
            isVar = Defaults{8};
            isPostpro = Defaults{9};
            
            if(isHold == 1)
                hold on;
            else
                figure;
            end
            
            if(nOut == 2 || isVar)
                [myMMmean, myMMvar] = arrayfun(@(x) uq_evalModel(MM,UtilFuncs.valReplace(paramList, paramInd, [x])),X);
                if(isPostpro == 1)
                    myMMmean(myMMmean > 100) = 100;
                    myMMmean(myMMmean < 0) = 0;
                end
                if(isVar)
                    yyaxis right
                    plot(X, sqrt(myMMvar))
                    hold on;
                    yyaxis left
                    plot(X, myMMmean)
                else
                    UtilFuncs.plotStochastic1D(X,myMMmean,sqrt(myMMvar))
                end
            else
                if(isPostpro == 1)
                    myMMmean = arrayfun(@(x) Models.postpro(MM,UtilFuncs.valReplace(paramList, paramInd, [x])),X);
                else
                    myMMmean = arrayfun(@(x) uq_evalModel(MM,UtilFuncs.valReplace(paramList, paramInd, [x])),X);
                end
                plot(X, myMMmean)
            end
            
            xlabel(labels(1))
            ylabel(labels(2))
        end
        
		% Function to plot mean and one standard deviation of given data 
		% ---------------------------------------------------------------
        function plotStochastic1D(X,myMMmean,myMMvar)
            plot(X, myMMmean,'LineWidth',1)
            hold on;
            curve1 = myMMmean + myMMvar;
            curve2 = myMMmean - myMMvar;
            x2 = [X, fliplr(X)];
            inBetween = [curve1, fliplr(curve2)];
            h = fill(x2, inBetween, '');
            set(h,'facealpha',.1,'HandleVisibility','off')
        end
		
		
		% Function to plot SVM prediction against applied voltage for given Econd parameters 
		% -----------------------------------------------------------------------------------
		function plotSVMprobValCurve(mySVC, Vs, ConParams)
            Xnew = ones(length(Vs),4).*ConParams;
            Xnew = [Xnew Vs'];
            cutoff = 0;
            probYhat_val = UtilFuncs.getProbVal(mySVC, Xnew);
            Yhat2 = zeros(size(probYhat_val));
            Yhat2(probYhat_val < cutoff) = -1;
            Yhat2(probYhat_val >= cutoff) = 1;
            plot(Vs, probYhat_val)
            hold on;
            plot(Vs, Yhat2)            
        end
        
		
		% Function to plot SVM prediction against a given data from file
		% ----------------------------------------------------------------
        function SVMvalidationplotwithData(mySVC, conVal, filename, X)
            Xnew = ones(size(X,2),4).*conVal;
            Xnew = [Xnew X'];
            
            load(filename);
            figure;
            yyaxis left
            plot(X, evals(1:length(X),1))
            hold on;
            yyaxis right
            UtilFuncs.plotSVMprobValCurve(mySVC, X, conVal)
            grid on;
            legend("True","ProbVal","Classification");
        end
        
    end
end