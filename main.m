% ======================================================================================
% MATLAB file to build the metamodels shown in the work
% "Surrogate Modeling in Irreversible Electroporation for Real-time Treatment Planning"
% ======================================================================================

clearvars
clc
uqlab
format long
rng default % for reproducibility

Markers = {'+','o','*','x','v','d','^','s','>','<'};
figTitles = ["RelTA", "TVol", "TAbl", "AblZne", "HtAbl"];
VarI = 1;
modelname = 'EP2Dnew.mph';


% Setup the input parameters
% ---------------------------
InputOpts = UtilFuncs.createInputs(5);
myInput = uq_createInput(InputOpts);

% Choose the phases.
% --------------------------------------------------------------------
% 1-MCS, 2-Load/train data, 3-SVM, 4-Build metamodels, 5-Plot results
% --------------------------------------------------------------------
phases = 1:5;

% Monte-Carlo Simulation (MCS)
% =============================
if ismember(1, phases)
    
    Vopts = [600, 700, 800, 900, 1000, 1079, 1200, 1400, 1600, 1800];
    
    if exist('MCSdata.mat','file') ~= 2
        
        % Choice of applied voltages for MCS
        % -----------------------------------
        SamplesList = [50,100,200,400,800];
        
        samEval = cell(length(Vopts), length(SamplesList));
        
        % Run the MCS for the chosen applied voltages
        % --------------------------------------------
        model = UtilFuncs.loadCOMSOLfile(modelname);
        i=0;
        for V = Vopts
            i = i + 1;
            j = 0;
            for noOfSamples = SamplesList
                j = j + 1;
                sample = uq_getSample(noOfSamples,'MC');
                sample(:,end) = ones(noOfSamples,1)*V;
                evals = UtilFuncs.solve_COMSOL_EPlist(model,sample);
                samEval{i,j} = {sample,evals};
            end
        end
        UtilFuncs.endComsolmphserver();
        
    else
        load('MCSdata');
    end
    
    % Calculate mean and standard deviation for each voltage value
    % -------------------------------------------------------------
    newSamEval = cell(size(samEval,1),size(samEval,2)*2-1);
    for i = 1:size(samEval,1)
        for j = 1:size(samEval,2)*2-1
            if (j == 1 || j==2)
                newSamEval{i,j} = samEval{i,j};
            elseif ((-1)^j == 1)
                newSamEval{i,j} = samEval{i,j-(j/2)+1};
            else
                newSamEval{i,j}{1, 1}  = [newSamEval{i,j-2}{1, 1} ;newSamEval{i,j-1}{1, 1} ];
                newSamEval{i,j}{1, 2} = [newSamEval{i,j-2}{1, 2} ;newSamEval{i,j-1}{1, 2} ];
            end
            
            mus(i,j) = mean(newSamEval{i,j}{1, 2}(:,VarI));
            stds(i,j) = std(newSamEval{i,j}{1, 2}(:,VarI));
        end
    end
    
    % Plot the convergence of the mean at each voltage value
    % -------------------------------------------------------
    figure,
    semilogy(abs(mus'-mus(:,end)')),legend(split(num2str(Vopts)))
    figure,
    UtilFuncs.plotStochastic1D(Vopts,mus(:,end)',stds(:,end)')
    
    if exist('MCSdata.mat','file') ~= 2
        
        % Save the MCS data to use as validation data
        % --------------------------------------------
        Xmcs = [];
        Ymcs = [];
        for i = 1:size(newSamEval,1)
            Xmcs = [Xmcs; newSamEval{i,end}{1, 1}];
            Ymcs = [Ymcs; newSamEval{i,end}{1, 2}];
        end
        
        save('MCSdata','samEval','Xmcs','Ymcs');
        
    end
    
end


% Loading the training and validation set
% =======================================
if ismember(2, phases)
    
    % Training data
    % -------------
    if exist('trainData.mat','file') == 2
        evals = load('trainData','params','QoI');
        sample = evals.params;
        evals = evals.QoI;
    else
        Npts = 1000;
        params = uq_getSample(Npts,'LHS');
        model = UtilFuncs.loadCOMSOLfile(modelname);
        QoI = UtilFuncs.solve_COMSOL_EPlist(model,params(1:end,:));
        UtilFuncs.endComsolmphserver();
        save('trainData','params','QoI');
        sample = params;
        evals = QoI;
    end
    
    % Clean the data
    % ---------------
    nullValues = find(~evals(:,1));
    sample(nullValues,:) = [];
    evals(nullValues,:) = [];
    
    % Shuffle the training data
    % --------------------------
    noOfTrainsets = length(evals);
    orgInd = 1:length(sample);
    newInd = orgInd(randperm(length(orgInd)));
    sample = sample(newInd,:);
    evals = evals(newInd,:);
    Xtrain = sample(1:noOfTrainsets,:);
    Ytrain = evals(1:noOfTrainsets,VarI);
    
    % Validation data
    % ----------------
    if ismember(1, phases)
        Xval = Xmcs(1551:12400,:); % A subset of data with V0 in parameter range was chosen
        Yval = Ymcs(1551:12400,VarI);
    else
        mcsdata = load('MCSdata');
        Xval = mcsdata.Xmcs(1551:12400,:); % A subset of data with V0 in parameter range was chosen
        Yval = mcsdata.Ymcs(1551:12400,VarI);
    end
end


% SVM metamodel training
% =======================
if ismember(3, phases)
    
    % Build SVM classifier
    % ---------------------
    [Xtrainc, Ytrainc] = Models.SVMdata(Xtrain, Ytrain);
    [~, Yvalc] = Models.SVMdata(Xval, Yval);
    
    muX = mean(Xtrain);
    StdX = std(Xtrain);
    %     Xtrainc = bsxfun(@rdivide, (bsxfun(@minus, Xtrain, muX)), StdX);
    Xvalc = bsxfun(@rdivide, (bsxfun(@minus, Xval, muX)), StdX);
    
    mySVC = Models.buildModel(Xtrainc, Ytrainc, Xvalc, Yvalc, "SVMC");
    
    
    % Evaluate the SVC metamodel at the validation set points
    % ---------------------------------------------------------
    [Yhat,probVal] = uq_evalModel(mySVC,Xvalc);
    
    % Model accuracy
    % ---------------
    accSVC = sum(Yhat == Yvalc)*100 /length(Yvalc); % Percentage of the SVC metamodel accuracy
    sensSVC = sum(Yhat == 1 & Yvalc == 1)*100/sum(Yvalc == 1); % Percentage of true positives: Senstivity
    specSVC = sum(Yhat == -1 & Yvalc == -1)*100/sum(Yvalc == -1); % Percentage of true negatives: Specificity
    falsePosSVC = sum(Yhat == 1 & Yvalc == -1)*100/sum(Yvalc == -1); % Percentage of false positives
    falseNegSVC = sum(Yhat == -1 & Yvalc == 1)*100/sum(Yvalc == 1); % Percentage of false negatives
    
end


% Kriging/GPR, PCE and PCK metamodel training
% ============================================
if ismember(4, phases)
    
    allmethodNames = ["GPR-linear", "GPR-exponential", "GPR-Gaussian", "GPR-Matern5/2", ...
        "PCE-OLS", "PCE-LARS", "PCE-OMP", "PCE-SP", "PCE-BCS", ...
        "PCK-Seq", "PCK-Opt", "GPR-Cust"];
    modlchoice = 1:12;
    metaModels = cell(length(modlchoice),1);
    
    % Empty containers for errors and build time
    % -------------------------------------------
    valerrs = nan(1,length(modlchoice));
    looerrs = nan(1,length(modlchoice));
    times = nan(1,length(modlchoice));
    methodNames = strings(1,length(modlchoice));
    
    % Build different metamodels with training data
    % ----------------------------------------------
    for ii = 1:length(modlchoice)
        [metaModels{ii}, valerrs(ii), looerrs(ii), times(ii)] = Models.buildModel(Xtrain, Ytrain, Xval, Yval, allmethodNames(modlchoice(ii)));
        methodNames(ii) = allmethodNames(modlchoice(ii));
    end
    
    save('MetaModels', 'metaModels', 'valerrs', 'looerrs', 'times', 'methodNames')
end


% Plot metamodel results
% =======================
if ismember(5, phases)
    
    if exist('times','var') ~= 1 && exist('valerrs','var') ~= 1 && exist('looerrs','var') ~= 1
        load('MetaModels');
    end
    
    % Error and Performance plots for the metamodels
    % -----------------------------------------------
    titleNames = ["Time", "Val Err", "LOO Err"];
    ylabels = ["sec","%","%"];
    plotdata = [times;valerrs;looerrs];
    
    for k = 1:3
        figure(),
        bar(plotdata(k,:))
        xticklabels(methodNames), xtickangle(45),
        title(titleNames(k)), ylabel(ylabels(k)),
        set(gca,'FontSize',16), set(gcf, 'Position', [100 100 700 500])
        grid on;
    end
    
    
    % Comparsion of validation and cross validation errors
    % -----------------------------------------------------
    figure(),
    bar(plotdata(2:3,:)')
    xticklabels(methodNames), xtickangle(45),
    title("Error Comparsion"), ylabel(ylabels(k)),
    set(gca,'FontSize',16), set(gcf, 'Position', [100 100 700 500])
    grid on;
    legend(["ValErr","LOOErr"], 'Location', 'best')
    
    
    % Figure plots from the paper
    % -----------------------------
    FigPlots
    
end
