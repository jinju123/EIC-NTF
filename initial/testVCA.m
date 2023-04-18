%% compare between MDC and MVC algorithm
    clear all
    load('newdata2.mat')
    load('Y.mat')
    % generate true W, H, V
    sampleNum = 900;
    noiseLevel = 0.2;
    bandNum = 5;
    endNum = 9;
    
    % compare several times
    HTrue = E;
    V = Y;
    
    % run VCA
    HVca = hyperVca(Y, endNum);
    
    % visualize result
    