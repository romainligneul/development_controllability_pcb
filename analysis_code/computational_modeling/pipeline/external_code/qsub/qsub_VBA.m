function qsub_VBA(y, u, evof, obsf, dim, options, output_file)
%QSUB_VBA Summary of this function goes here
%   Detailed explanation goes here

[posterior, out] = VBA_NLStateSpaceModel(y, u, evof, obsf, dim, options);

GoF =  out.F;
GoF =  out.fit.BIC;
GoF =  out.fit.AIC;
if ~isempty(posterior.muPhi)
 phiRaw=posterior.muPhi;
 for pp=1:length(posterior.muPhi)
     phiFitted(1,pp)=options.inG.param_transform{pp}(posterior.muPhi(pp));
 end
end
if ~isempty(posterior.muPhi)
 thetaRaw=posterior.muTheta;
 for pp=1:length(posterior.muTheta)
     thetaFitted(1,pp)=options.inF.param_transform{pp}(posterior.muTheta(pp));
 end
end
muX = out.suffStat.muX;
suffStat = out.suffStat;

keep_u = out.u;

save([output_file],'suffStat', 'muX', 'phiFitted','thetaFitted', 'thetaRaw', 'phiRaw', 'GoF', 'keep_u', 'options');


end

