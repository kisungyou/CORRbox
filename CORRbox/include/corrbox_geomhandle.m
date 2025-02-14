function obj = corrbox_geomhandle(geometry)

% this is an internal function to create a struct of function handles so
% that it controls multiple geometries at once. In this handler, the
% following functions are supported. Geometries are care-insensitive.
%   
%   dist  : pairwise distance
%   fmean : Frechet mean and variation given a 3d array

obj = struct;
glc = lower(geometry);

switch glc
    case "euclidean"
        obj.pga   = @(x3,ndim) euclid_pga(x3, ndim);
        obj.dist  = @(x,y) euclid_dist(x,y);
        obj.pdist = @(x3) euclid_pdist(x3);
        obj.fmean = @(x3,maxiter,thr) euclid_mean(x3);
        obj.fmedian = @(x3,maxiter,thr) euclid_median(x3,maxiter,thr);
        obj.localise = @(x3, xref) euclid_localise(xref, x3);
        obj.localcoord = @(x3) euclid_coord_near_mean(x3);
        obj.pred_kernreg = @(x3new, xloc, model, funcname) euclid_pred_kernreg(x3new, xloc, model, funcname);
    case "qam"
        obj.pga   = @(x3,ndim) qam_pga(x3,ndim);
        obj.dist  = @(x,y) qam_dist(x,y);
        obj.pdist = @(x3) qam_pdist(x3);
        obj.fmean = @(x3,maxiter,thr) qam_mean(x3,maxiter,thr);
        obj.fmedian = @(x3,maxiter,thr) qam_median(x3,maxiter,thr);
        obj.localise = @(x3, xref) qam_localise(xref, x3);
        obj.localcoord = @(x3) qam_coord_near_mean(x3);
        obj.pred_kernreg = @(x3new, xloc, model, funcname) qam_pred_kernreg(x3new, xloc, model, funcname);
    case "lec"
        obj.pga   = @(x3,ndim) lec_pga(x3,ndim);
        obj.dist  = @(x,y) lec_dist(x,y);
        obj.pdist = @(x3) lec_pdist(x3);
        obj.fmean = @(x3,maxiter,thr) lec_mean(x3);
        obj.fmedian = @(x3,maxiter,thr) lec_median(x3,maxiter,thr);
        obj.localise = @(x3, xref) lec_localise(xref, x3);
        obj.localcoord = @(x3) lec_coord_near_mean(x3);
        obj.pred_kernreg = @(x3new, xloc, model, funcname) lec_pred_kernreg(x3new, xloc, model, funcname);
    case "ecm"
        obj.pga   = @(x3,ndim) ecm_pga(x3, ndim);
        obj.dist  = @(x,y) ecm_dist(x,y);
        obj.pdist = @(x3) ecm_pdist(x3);
        obj.fmean = @(x3,maxiter,thr) ecm_mean(x3);
        obj.fmedian = @(x3,maxiter,thr) ecm_median(x3,maxiter,thr);
        obj.localise = @(x3,xref) ecm_localise(xref,x3);
        obj.localcoord = @(x3) ecm_coord_near_mean(x3);
        obj.pred_kernreg = @(x3new, xloc, model, funcname) ecm_pred_kernreg(x3new, xloc, model, funcname);
    otherwise
        warning("* CORRbox : the provided 'geometry' is not supported. Return the 'ECM' geometry.")
        obj = corrbox_geomhandle("ECM");
end



end