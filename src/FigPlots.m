% ======================================================================================
% MATLAB file to reproduce the figure plots and tables in the work
% "Surrogate Modeling in Irreversible Electroporation for Real-time Treatment Planning"
% ======================================================================================

function FigPlots(plots)

uqlab

if exist('metaModels','var') ~= 1 && exist('methodNames','var') ~= 1
    MMs = load('MetaModels');
    methodNames = MMs.methodNames;
    metaModels = MMs.metaModels;
end

% methodNames = ["GPR-linear", "GPR-exponential", "GPR-Gaussian", "GPR-Matern5/2", ...
% 	"PCE-OLS", "PCE-LARS", "PCE-OMP", "PCE-SP", "PCE-BCS", ...
% 	"PCK-Seq", "PCK-Opt", "GPR-Cust"];
% metaModels = [myKrigingLin myKrigingExp myKrigingGaus myKrigingMat52 ...
%         myPCE_OLS myPCE_LARS myPCE_OMP myPCE_SP myPCE_BCS ...
%         mySeqPCK mySeqPCKopt myKrigingCust];

if exist('plots','var') ~=1
    plots = 7;
end


% Two regions UR and AR based on different FEM data configurations (Figure - 3)
% ------------------------------------------------------------------------------
if ismember(1,plots)
    Nodatas = [6 3 4 8 0 7];
    Xplot = 1000:10:1800;
    N = length(Xplot);
    plotData = zeros(length(Nodatas)*N,3);
    k=0;
    for i=Nodatas
        k=k+1;
        load(join(['Config',num2str(i)]));
        plotData((k-1)*N+1:k*N,2) = Xplot';
        plotData((k-1)*N+1:k*N,3) = evals(1:N,1)';
        plotData((k-1)*N+1:k*N,1) = k*ones(N,1)';
    end
    
    colors = ['r','g','b','m','k','y','c'];
    figure(),
    for ii=plotData'
        hold on;
        if (ii(3) == 100)
            scatter3(ii(1),ii(2),ii(3),50,colors(ii(1)),'+')
        else
            scatter3(ii(1),ii(2),ii(3),40,colors(ii(1)),'o')
        end
    end
    
    set(gca,'TickLabelInterpreter','latex');
    set(gca,'FontSize',20)
    xlabel('Conductivity parameter configuration')
    ylabel("V_0")
    zlabel("y_{RelTA}")
    
    xticks(1:6)
    xticklabels({'$\mathrm{\mathbf{C}}_1$','$\mathrm{\mathbf{C}}_2$','$\mathrm{\mathbf{C}}_3$', '$\mathrm{\mathbf{C}}_4$','$\mathrm{\mathbf{C}}_5$','$\mathrm{\mathbf{C}}_6$'})
    grid on;
    grid minor;
    set(gca,'XMinorGrid','off')
    view(-50,40);
end


% Plot MCS data (Figure - 4)
% ---------------------------
if ismember(2,plots)
    load('MCSdata');
    Vopts = [600, 700, 800, 900, 1000, 1079, 1200, 1400, 1600, 1800];
    mus = -1*ones(1,10);
    stds = -1*ones(1,10);
    nSolsperV = 1550;
    
    for i = 1:10
        inds = (i-1)*nSolsperV+1:(i)*nSolsperV;
        mus(i) = mean(Ymcs(inds,1));
        stds(i) = std(Ymcs(inds,1));
    end
    figure(),
    UtilFuncs.plotStochastic1D(Vopts,mus,stds)
    xlabel("V_0")
    ylabel("y_{RelTA}")
    grid on; set(gca, 'FontSize', 20)
    hold on;
    grid minor;
    plot(Vopts,mus,'kx')
    xlim([600, 1800]);
end


% Calculate RMSE for GP, PCK and GP-custom (Table - 4 & C.6)
% -----------------------------------------------------------
if ismember(3,plots)
    modlchoice = 1:12; % Choice based on the variable "metaModels"
    Ebin = [0 85 100 101];
    errorBinorg = -1*ones(length(modlchoice),length(Ebin));
    errorBin = -1*ones(length(modlchoice),length(Ebin));
    
    for i = modlchoice
        errorBinorg(i,:)=[Models.errorBins(metaModels{i}, metaModels{i,1}.Internal.ValidationSet.X, metaModels{i,1}.Internal.ValidationSet.Y, Ebin), ...
            Models.RMSEcalc(metaModels{i}, metaModels{i,1}.Internal.ValidationSet.X, metaModels{i,1}.Internal.ValidationSet.Y)];
        errorBin(i,:)=[Models.errorBins(metaModels{i}, metaModels{i,1}.Internal.ValidationSet.X, metaModels{i,1}.Internal.ValidationSet.Y,Ebin,1), ...
            Models.RMSEcalc(metaModels{i}, metaModels{i,1}.Internal.ValidationSet.X, metaModels{i,1}.Internal.ValidationSet.Y,1)];
    end
end


% Compare GP, PCK and GP-custom against Validation data points (Figure - 5 & C.8)
% --------------------------------------------------------------------------------
if ismember(4,plots)
    modlchoice = 1:9; % Choice based on the variable "metaModels"
    % Set modlchoice as 12 for Figure 5  and  1:9 for Figure C.8 when all
    % 12 metamodels are built
    
    for i = modlchoice
        if(strcmp(metaModels{i}.Options.MetaType, 'Kriging'))
            [GPmean,GPstd] = uq_evalModel(metaModels{i}, metaModels{i,1}.Internal.ValidationSet.X);
            
            if i ==12 % (Figure - 5)
                figure,
                uq_plot(metaModels{i,1}.Internal.ValidationSet.Y, GPstd, 'r+')
                xlabel("$y_{RelTA}$")
                ylabel("Variance of predicted $y_{RelTA}$")
                set(gca,'FontSize',20)
                grid minor;
                title(methodNames(i))
            end
        else
            GPmean = uq_evalModel(metaModels{i}, metaModels{i,1}.Internal.ValidationSet.X);
        end
        figure,
        UtilFuncs.plotValSet(metaModels{i,1}.Internal.ValidationSet.Y, GPmean, '', 'Ytrue', 'YGP');
        xlabel("$y_{RelTA}$")
        ylabel("Predicted $y_{RelTA}$")
        set(gca, 'FontSize', 20)
        grid minor;
        title(methodNames(i))
    end
end

% Compare GP, PCK and GP-custom for different FEM data configurations (Figure - 7)
% ---------------------------------------------------------------------------------
if ismember(5,plots)
    Nodatas = [6 3 4 8 0 7]; % Choice based on the data in "Config" files
    Xplot = 1000:10:1800;
    modlchoice = [4 11 12]; % Choice based on the variable "metaModels"
    legendInfos = arrayfun(@(x) methodNames(x), modlchoice);
    legendInfos = [legendInfos, "True"];
    set(groot,'defaultLineLineWidth',1.5)
    
    for ii=Nodatas
        load(join(['Config',num2str(ii)]));
        figure;
        for i = modlchoice
            UtilFuncs.plotMM1D(Xplot, metaModels{i}, sample(1,:), 5, 1, ['V', 'RelTA'], 1, 0);
        end
        hold on;
        plot(Xplot, evals(1:length(Xplot),1),'--')
        legend(legendInfos, 'Location', 'southeast');
        xlabel("V_0")
        ylabel("y_{RelTA}")
        grid on;set(gca, 'FontSize', 20)
        title(join(string(sample(1,1:4))," , "))
    end
end


% 3D plot of metamodel against two input parameters
% --------------------------------------------------
if ismember(6,plots)
    refVals = [0.411, 0.075, 1.2826, 3.3684, 1079];
    modlchoice = [4 11 12]; % Choice based on the variable "metaModels"
    paramchoice = [1 5]; % \sigma_{T0} and V_0
    
    for i = modlchoice
        figure,
        UtilFuncs.plotMM2D(metaModels{i}.Internal.ValidationSet.X ,metaModels{i}, refVals, paramchoice, 30);
        zlabel('RelTA')
        xlabel('\sigma_{T0}')
        ylabel('V_0')
        title(methodNames(i))
    end
end

% Results of treatment plannong (Figure - 6)
% ---------------------------------------------
if ismember(7,plots)
    if exist('resultsCOMSOL','var') ~=1
        TPopt = load('data\TPopt');
    end
    
    set(groot,'defaultLineLineWidth',1.5, 'defaultAxesFontSize', 20)
    
    InputOpts = UtilFuncs.createInputs(4, true);
    
    subplot(2,4,1)
    yyaxis left
    histogram(TPopt.paramconfigs(:,1),20,'Normalization','probability')
    xlabel('\sigma_{T0}')
    ylabel('Probability')
    hold on;
    UtilFuncs.plotdist(InputOpts.Marginals(1).Parameters(1), InputOpts.Marginals(1).Parameters(2), 'Normal')
    grid minor;
    
    subplot(2,4,2)
    yyaxis left
    histogram(TPopt.paramconfigs(:,2),20,'Normalization','probability')
    xlabel('\sigma_{N0}')
    hold on;
    UtilFuncs.plotdist(InputOpts.Marginals(2).Parameters(1), InputOpts.Marginals(2).Parameters(2), 'Normal')
    grid minor;
    
    subplot(2,4,3)
    yyaxis left
    histogram(TPopt.paramconfigs(:,3),20,'Normalization','probability')
    xlabel('\theta_T')
    hold on;
    UtilFuncs.plotdist(InputOpts.Marginals(3).Parameters(1), InputOpts.Marginals(3).Parameters(2), 'Normal')
    grid minor;
    
    subplot(2,4,4)
    yyaxis left
    histogram(TPopt.paramconfigs(:,4),20,'Normalization','probability')
    xlabel('\theta_N')
    hold on;
    UtilFuncs.plotdist(InputOpts.Marginals(4).Parameters(1), InputOpts.Marginals(4).Parameters(2), 'Normal')
    grid minor;
    
    subplot(2,4,[5 6])
    yyaxis left
    histogram(TPopt.resultsCOMSOL,20,'Normalization','probability')
    xlabel('Optimum V_0 with FEM model')
    ylabel('Probability')
    hold on;
    pdCOMSOL = fitdist(TPopt.resultsCOMSOL','LogNormal');
    UtilFuncs.plotdist(exp(pdCOMSOL.mu), exp(pdCOMSOL.sigma), pdCOMSOL.DistributionName)
    grid minor;
    
    subplot(2,4,[7 8])
    yyaxis left
    histogram(TPopt.resultsGPR,20,'Normalization','probability')
    xlabel('Optimum V_0 with GPR custom model')
    hold on;
    pdGPR = fitdist(TPopt.resultsGPR','LogNormal');
    UtilFuncs.plotdist(exp(pdGPR.mu), exp(pdGPR.sigma), pdGPR.DistributionName)
    grid minor;
    
end

end