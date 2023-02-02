function doClassification( measF, clasF, classifiers, folds, nrrep, statDisp)
%
% Function for general classification
% Input parameters:
%       measF - name of the text file containing measurements, e.g., 'featureVectors.txt'
%       clasF - name of the text file containing reference class for each
%         measurement in a measF file, e.g., 'referenceClass.txt'
%       classifiers - cell array of integere values (0/1) indicating which
%         of the available classifiers to use, e.g.,  {1,1}
%         available classifiers:
%            LDA - Linear Discriminant Analysis classifier 
%            QDA - Quadratic Discriminant Analysis classifier
%       folds    - number of folds to use when performing cross-validation
%       nrrep    - the number of overall repetitions
%       statDisp - whether to display statistics (0/1)
% Example of a command: 
%       doClassification('featureVectors.txt', 'referenceClass.txt', {1, 1}, 5, 30, 0)
% This function presumes, that the class of the first measurement represents 
% the "EVENT" class, while the second class represents the "NON-EVENT" class 
%

    meas=mytextread(measF);
    clas=textread(clasF,'%s');

    unq=unique(clas);
    [~, idxs, ~] = unique(clas, 'first');
    unq = clas(sort(idxs));
    if (size(unq,1)<2)
        return;
    end

    clas1=unq{1,1};
    clas2=unq{2,1};
    fprintf('Size: %dx%d; EVENT: %s; NON-EVENT: %s\n',size(meas,1), size(meas,2), clas1, clas2);

    if (statDisp==1)

        % Font size
        fSize=12;

        idx1=find(strcmp(clas(:), clas1));
        idx2=find(strcmp(clas(:), clas2));

        for cIdx=1:size(meas,2)

            % p - Student's t test
            [h, p] = ttest2(meas(idx1,cIdx), meas(idx2,cIdx));

            fprintf('Difference between the group %s and %s, column %2d\n', clas1, clas2,cIdx);
            fprintf('  p-values: %f (student), %f (anova)\n',p,p); 

            % Boxplots

            fig1=figure; % No saving, OR
            %picName=sprintf('ANOVA-%s-Col%02d.eps',regexprep(measF, '\..*',''), cIdx);
            %fig1=figure('Name', picName);

            set(gca,'FontSize',fSize);
            aLimit=[floor(min(meas(:,cIdx))) ceil(max(meas(:,cIdx)))];
            ylim(aLimit);

            boxplot([meas(idx1,cIdx),meas(idx2,cIdx)],'Labels',{clas1,clas2});
            hold all;

            text(gca,1.25,(double(aLimit(2))-double(aLimit(1)))/2.0, strcat('p=',num2str(p)));
            text1=sprintf('Boxplots, feature %d', cIdx);
            title(gca,text1);
            ylabel(gca,'Feature value');
            set(findall(gcf,'-property','FontSize'),'FontSize',fSize)
            %print('-depsc',picName);

            % Scatterplot

            fig2=figure; % No saving, OR
            %picName=sprintf('Scatter-%s-Col%02d.eps',regexprep(measF, '\..*',''), cIdx);
            %fig2=figure('Name', picName);

            set(gca,'FontSize',fSize);
            aLimit=[floor(min(meas(:,cIdx))) ceil(max(meas(:,cIdx)))];
            ylim(aLimit);

            scatter(idx1, meas(idx1,cIdx), 'fill');
            hold all;
            scatter(idx2, meas(idx2,cIdx));

            text(gca,double(length(meas))/2.0,(double(aLimit(2))-double(aLimit(1)))/2.0, strcat('p=',num2str(p)));
            text1=sprintf('Scatterplot, feature %d', cIdx);
            title(gca,text1);
            ylabel(gca,'Feature value');
            xlabel(gca,'Feature index');
            set(findall(gcf,'-property','FontSize'),'FontSize',fSize) 
            set(gca,'box','on'); % display the top and right border of the window as well  
            %print('-depsc',picName);
         
        pause

        end
    end
    classifyMeasurements(meas, clas, clas1, clas2, classifiers, folds, nrrep);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  classifyMeasurements( meas, clas, clas1, clas2, classifiers, folds, nrrep)
%
% The function that performs the actual classification
%   meas  - measurements (one line per measurement)
%   clas  - true class designations for each measurement
%   clas1 - class of the "EVENT", e.g., imagining of left hand movement
%   clas2 - class of the "NON-EVENT", e.g., imagining of right hand movement
%   classifiers - which classifiers to use (LDA, QDA)
%   folds - number of folds when doing cross-validation
%   nrepp - the number of averall repetitions

    warning('off', 'stats:obsolete:ReplaceThisWithMethodOfObjectReturnedBy');
    warning('off', 'stats:obsolete:ReplaceThisWith');

    if (size(classifiers,2)<size(classifiers,1))
        classifiers=classifiers';
    end
    
    % Partitioning the measurements to train and test sets - holdout
    % Holdout set to 20%
    P = cvpartition(clas, 'Holdout', 0.20);
    P0 = P;
    
    % Partitioning the measurements into folds - cross-validation
    cvFolds = crossvalind('Kfold',clas, folds);
    cvFolds0 = cvFolds;

    for idx=1:size(classifiers,2)

         % Classification 
         if (classifiers{idx}~=0)

                if idx==1
                    fprintf('LDA\n');
                    classifier = ClassificationDiscriminant.fit(meas, clas, 'discrimType', 'pseudoLinear');
                elseif idx==2
                    fprintf('QDA\n');
                    classifier = ClassificationDiscriminant.fit(meas, clas, 'discrimType', 'pseudoQuadratic');
                end

                [C, score, ~] = classifier.predict(meas);
                [~, score2, ~] = resubPredict(classifier);
                probsCl = score(:,1); % probsCl = score(:,2);

                [TP, FN, FP, TN] = evaluateClassifier(C, clas, clas1, clas2);
                %fprintf('\t all records Se: %5.2f Sp: %5.2f\n', TP*100.0/(TP+FN), TN*100.0/(TN+FP));
                %probsTr=[];
                %probsCl=[];
                %for cIdx=1:size(clas)
                %    if (strcmp(clas{cIdx}, C{cIdx})==1) probsCl(cIdx)=1;
                %    else probsCl(cIdx)=0;
                %    end
                %end
                %score'
                [aucX,aucY, ~, auc] = perfcurve(clas, probsCl, clas1);
                [aucX1,aucY1, ~, auc1] = perfcurve(clas, score2(:,1), clas1); %  ,score2(:,2),
                %display(auc);
                %figure;
                %plot (aucX,aucY);
                displayResults('All records-Test  ', TP, FN, FP, TN, auc);
                displayResults('All records-Learn ', TP, FN, FP, TN, auc1);


         %Holdout
         TN2=0;
         TP2=0;
         FN2=0;
         FP2=0; 
         auc2=0;
         auc12=0;

         for myidx=1:nrrep
            if myidx==1
                P = P0;
            else
                P = cvpartition(clas, 'Holdout', 0.20);
            end

                if idx==1
                    classifier = ClassificationDiscriminant.fit(meas(P.training,:), clas(P.training,:), 'discrimType', 'pseudoLinear');
                elseif idx==2
                    classifier = ClassificationDiscriminant.fit(meas(P.training,:), clas(P.training,:), 'discrimType', 'pseudoQuadratic');
                end

                [C, score,~] = classifier.predict(meas(P.test,:));
                [~, score2] = resubPredict(classifier);
                probsCl = score(:,1); % probsCl = score(:,2);

                [TP, FN, FP, TN] = evaluateClassifier(C, clas(P.test), clas1, clas2);
                % fprintf('\t 20%% holdout Se: %5.2f Sp: %5.2f\n', TP*100.0/(TP+FN), TN*100.0/(TN+FP));  
                % displayResults('all records', TP, FN, FP, TN);

                TN2=TN2+TN;
                TP2=TP2+TP;
                FN2=FN2+FN;
                FP2=FP2+FP;
 
                [aucX,aucY, ~, auc] = perfcurve(clas(P.test,:), probsCl, clas1);
                [aucX1,aucY1, ~, auc1] = perfcurve(clas(P.training,:), score2(:,1), clas1); %  ,score2(:,2),

                auc2=auc2+auc;
                auc12=auc12+auc1;

                if myidx==1                
                  %display(auc);
                  %figure;
                  %plot (aucX,aucY);
                  displayResults('20% holdout-Test  Single ', TP, FN, FP, TN, auc);
                  displayResults('20% holdout-Learn Single ', TP, FN, FP, TN, auc1);
                elseif myidx==nrrep
                  TP2=round(TP2*1.0/nrrep);
                  TN2=round(TN2*1.0/nrrep);
                  FP2=round(FP2*1.0/nrrep);
                  FN2=round(FN2*1.0/nrrep);
                  auc2=auc2/nrrep;
                  auc12=auc12/nrrep;
                  displayResults('  20% holdout-Test  Nrep ', TP2, FN2, FP2, TN2, auc2);
                  displayResults('  20% holdout-Learn Nrep ', TP2, FN2, FP2, TN2, auc12);
                end
         end


         %Crossvalidation
         TN2=0;
         TP2=0;
         FN2=0;
         FP2=0; 
         auc2=0;
         auc12=0;

         for myidx=1:nrrep
            if myidx==1
               cvFolds = cvFolds0;
            else
               cvFolds = crossvalind('Kfold',clas, folds);
            end

                TN=0;
                TP=0;
                FN=0;
                FP=0;
                probsCl=[];
                probsCl2=[];
                probsTr=[];
                auc1=[];
                for fdx=1:folds
                    testIdx = (cvFolds == fdx);

                    trainIdx = ~testIdx;
                    if idx==1                  
                        classifier = ClassificationDiscriminant.fit(meas(trainIdx,:), clas(trainIdx,:), 'discrimType', 'pseudoLinear');
                    elseif idx==2
                        classifier = ClassificationDiscriminant.fit(meas(trainIdx,:), clas(trainIdx,:), 'discrimType', 'pseudoQuadratic');
                    end

                    [C, score, ~] = classifier.predict(meas(testIdx,:));
                    [C2,score2] = resubPredict(classifier);
                    probsCl=[probsCl;score(:,1)];   %  probsCl=[probsCl;score(:,2)];                   
                    [X1,Y1,~, tauc1] = perfcurve(clas(trainIdx), score2(:,1), clas1);  %  ,score2(:,2), 
                    auc1=[auc1 tauc1];

                    [TP1, FN1, FP1, TN1] = evaluateClassifier(C, clas(testIdx), clas1, clas2);

                    TP=TP+TP1;
                    FN=FN+FN1;
                    FP=FP+FP1;
                    TN=TN+TN1;

                    TP2=TP2+TP1;
                    FN2=FN2+FN1;
                    FP2=FP2+FP1;
                    TN2=TN2+TN1;

                    probsTr=[probsTr;clas(testIdx)];
                end
                %fprintf('\t cross-valid Se: %5.2f Sp: %5.2f\n', TP*100.0/(TP+FN), TN*100.0/(TN+FP));
                [aucX,aucY, ~, auc] = perfcurve(probsTr, probsCl, clas1);
                %display(auc);
                %figure;
                %plot(aucX,aucY);

                auc2=auc2+auc;
                auc12=auc12+mean(auc1);

                if myidx==1                
                   %display(auc);
                   %figure;
                   %plot(aucX,aucY);
                   displayResults('Cross-valid-Test  Single ', TP, FN, FP, TN, auc);
                   displayResults('Cross-valid-Learn Single ', TP, FN, FP, TN, mean(auc1));
                elseif myidx==nrrep                
                %display(auc);
                %%
                figure;
                plot(aucX,aucY);
                if (idx == 1)
                    text1=sprintf('ROC curve, %s, cross validation, %d folds','LDA',folds);
                elseif (idx == 2)
                    text1=sprintf('ROC curve, %s, cross validation, %d folds','QDA',folds);
                end 
                title(text1);
                xlabel('1 - Specificity');
                ylabel('Sensitivity');
                %%

                   TP2=round(TP2*1.0/nrrep);
                   TN2=round(TN2*1.0/nrrep);
                   FP2=round(FP2*1.0/nrrep);
                   FN2=round(FN2*1.0/nrrep);
                   auc2=auc2/nrrep;
                   auc12=auc12/nrrep;
                   displayResults('  Cross-valid-Test  Nrep ', TP2, FN2, FP2, TN2, auc2);
                   displayResults('  Cross-valid-Learn Nrep ', TP2, FN2, FP2, TN2, auc12);
                end
         end

         %Leave-One-Out
                clasC={};      %zeros(size(meas,1),1);
                clasTrth={};
                probsCl=[];
                aucAgg=[];
                for k=1:size(meas,1)
                    measSingle = meas(k,:);  % pick the t-th measurement
                    measTmp = meas;          % copy measurement 
                    clasTmp = clas;          % copy classes
                    clasTrth{end+1}=clas{k}; % remember the class of current measurement
                    measTmp(k,:)=[];         % erase measurement from the training set
                    clasTmp(k)=[];           % erase class from the training set
                    if idx==1
                        classifier = ClassificationDiscriminant.fit(measTmp, clasTmp, 'discrimType', 'pseudoLinear');
                    elseif idx==2
                        classifier = ClassificationDiscriminant.fit(measTmp, clasTmp, 'discrimType', 'pseudoQuadratic');
                    end
                    
                    [C, scoreS, ~] = classifier.predict(measSingle);
                    probsCl=[probsCl;scoreS(:,1)];   %    probsCl=[probsCl;scoreS(:,2)];
                    [C2, score2] = resubPredict(classifier);
       
                    clasC{end+1}=C{1};
                    [aX1, aY1, ~, tauc] = perfcurve(clasTmp, score2(:,1), clas1); %  ,score2(:,2),
                    aucAgg=[aucAgg tauc];
                end
                clasC=clasC';
                clasTrth=clasTrth';                
                %[TP, FN, FP, TN] = evaluateClassifier(clasC, clasTrth, clas1, clas2);
                %displayResults('leave one out', TP, FN, FP, TN);
                [TP, FN, FP, TN] = evaluateClassifier(clasC, clas, clas1, clas2);
                [aucX,aucY, ~, auc] = perfcurve(clasTrth, probsCl, clas1);
                %display(auc);
                %figure;
                %plot (aucX,aucY);
                displayResults('Leave one o-Test  ', TP, FN, FP, TN,auc);
                displayResults('Leave one o-Learn ', TP, FN, FP, TN, mean(aucAgg));

         end  % if (classifiers{idx}~=0)
    end     % for idx=1:size(classifiers,2)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ TP, FN, FP, TN ] = evaluateClassifier( pred, orig, cls1, cls2 )
%
% The function to evaluate performance 
%   pred - predicted class of measurements
%   orig - true class of measurements
%   cls1 - class of the "EVENT"     (Class_1)
%   cls2 - class of the "NON-EVENT" (Class_2)
%

    TP=0; FN=0; TN=0; FP=0;
   
    for i = 1:size(pred)
        if     (strcmp(pred{i},orig{i})==1) & (strcmp(pred{i},cls1)==1)
            TP=TP+1;
        elseif (strcmp(pred{i},orig{i})==0) & (strcmp(pred{i},cls2)==1)
            FN=FN+1;
        elseif (strcmp(pred{i},orig{i})==1) & (strcmp(pred{i},cls2)==1)
            TN=TN+1;
        elseif (strcmp(pred{i},orig{i})==0) & (strcmp(pred{i},cls1)==1)
            FP=FP+1;
        else
            fprintf('Error *%s +%s =%s =>%s\n',pred{i},orig{i},cls1,cls2);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function displayResults(aText, TP, FN, FP, TN, AUC)
%
% The function for displaying the results
%   aText - description of a result
%   TP    - true positives
%   FN    - false negative
%   FP    - false positives
%   TN    - true negatives
%
% Mode of display, short=0; - long; short=1; - long in line; short=2; - brief;
%

    short=1;
    if (short==0)
        fprintf('  %s   TP=%3d FN=%3d\n',aText, TP, FN);
        fprintf('                   FP=%3d TN=%3d\n', FP, TN);
        fprintf('                   Se: %6.2f Sp: %6.2f\n',TP*100.0/(TP+FN), TN*100.0/(TN+FP));
    elseif (short==1)
        if (nargin==5)
            fprintf('   %s TP=%3d FN=%3d; FP=%3d TN=%3d; Se: %6.2f Sp: %6.2f\n',aText,TP,FN,FP,TN,TP*100.0/(TP+FN), TN*100.0/(TN+FP));
        else
            fprintf('   %s TP=%3d FN=%3d; FP=%3d TN=%3d; Se: %6.2f Sp: %6.2f CA: %6.2f AUC: %6.2f\n',aText,TP,FN,FP,TN,TP*100.0/(TP+FN),TN*100.0/(TN+FP), (TP+TN)*100.0/(TP+FN+FP+TN), AUC*100);
        end
    else
        fprintf('   %s Se: %6.2f Sp: %6.2f\n',aText, TP*100.0/(TP+FN), TN*100.0/(TN+FP));
    end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ meas ] = mytextread( fileName )
%
%  My 'textread' function, which ignores leading, in-between and trailing
%  whitespaces, intended to replace original 'textread' function. This
%  version of the 'textread' function reads a file with n rows containing 
%  m measurements separated with arbitrary whitespaces. It is not intended
%  for reading text files containing non-numerical values e.g. strings. 
%      Input arguments:
%               fileName - name of the text file containging measurements
%      Output arguments: 
%               meas     - a nxm size matrix of numerical values
%

mfid=fopen(fileName, 'rt');
  meas=[];
  while ~feof(mfid)
      meas=[meas; str2num(fgetl(mfid))];
  end
  fclose(mfid);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


