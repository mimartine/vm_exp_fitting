function [rmsd] = rmsd(preds, data)

sd = (preds-data).^2;
rmsd = sqrt(sum(sd)/numel(sd));

end

