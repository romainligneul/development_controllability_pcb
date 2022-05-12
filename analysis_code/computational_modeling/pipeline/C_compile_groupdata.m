%%%%% Concatenate into group data %%%%%
% only required if the fits are re-run
% this script produces the fitted_model.mat file used by D-I
% romain.ligneul@gmail.com

clear all

fitdir = 'fit_data_rerun/'
dirlist = dir(fitdir)
dirlist(1:2)=[];
dirlist=dirlist([dirlist.isdir])

for m=1:length(dirlist)
    
    slist = dir([fitdir dirlist(m).name '/*.mat']);
    
    removeind=strmatch('fitted_model.mat', {slist.name}');
    slist(removeind)=[];
    
    clear GoF phiFitted thetaFitted u muX options phiRaw thetaRaw
    
    for s=1:length(slist)
        sdat = load([fitdir dirlist(m).name '/' slist(s).name]);
        GoF(s,:)=sdat.GoF;
        
        try
            phiFitted(s,:)=sdat.phiFitted;
        end
        
        try
            thetaFitted(s,:)=sdat.thetaFitted;
        end
        
        try
            u{s}=sdat.keep_u;
        catch
            u{s}=sdat.u;
        end
        
        muX{s}=sdat.muX;
        options = sdat.options;
        
        try
            phiNames=sdat.phiNames;
            thetaNames=sdat.thetaNames;
        end
        
        Traw= @(x) x;
        InvMin1to1= @(x) -log((2./(x+1))-1);
        
        if isfield(sdat', 'phiRaw')
            phiRaw(s,:)=sdat.phiRaw;
        else
            for ppp=1:length(sdat.phiFitted)
                if strcmp(char(options.inG.param_transform{ppp}),'@(x)exp(x)*5')
                    phiRaw(s,ppp) = log(sdat.phiFitted(ppp)/5);
                elseif strcmp(char(options.inG.param_transform{ppp}),'@(x)-1+2*VBA_sigmoid(x)')
                    phiRaw(s,ppp) = InvMin1to1(sdat.phiFitted(ppp));
                elseif strcmp(char(options.inG.param_transform{ppp}),'@(x)x');
                    phiRaw(s,ppp) = sdat.phiFitted(ppp);
                else
                    error('unknown transform')
                end
            end
        end
        
        if isfield(sdat', 'thetaRaw')
            thetaRaw(s,:)=sdat.thetaRaw;
        else
            for ppp=1:length(sdat.thetaFitted)
                if strcmp(char(options.inF.param_transform{ppp}),'@(x) VBA_sigmoid(x)') || strcmp(char(options.inF.param_transform{ppp}),'@(x)VBA_sigmoid(x)')
                    thetaRaw(s,ppp) = VBA_sigmoid(sdat.thetaFitted(ppp),'inverse', true);
                elseif strcmp(char(options.inF.param_transform{ppp}),'@(x)exp(x)*5')
                    thetaRaw(s,ppp) = log(sdat.thetaFitted(ppp)/5);
                elseif strcmp(char(options.inF.param_transform{ppp}),'@(x)-1+2*VBA_sigmoid(x)')
                    thetaRaw(s,ppp) = InvMin1to1(sdat.thetaFitted(ppp));
                elseif strcmp(char(options.inF.param_transform{ppp}),'@(x)x');
                else
                    error('unknown transform')
                end
            end
        end
        
    end
    
    clear Traw InvMin1to1;
    
    
    save([fitdir dirlist(m).name '/fitted_model.mat'])
    
end