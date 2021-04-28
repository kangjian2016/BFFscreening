



#' Find the neighbors for each voxel in images
#' @param img an nifti object.
#' @param mask an nifti object.
#' \code{mask > 0} specifies which voxels are on the mask.
#' Default value is NULL indicating all voxels are considered.
#' @param radius the size (voxel) of neighborhood.
#' @return the neighbor indices of each voxel in each row.
#' @author Jian Kang <jiankang@umich.edu>
#' @examples
#' maskfile <- file.path(system.file("nifti", package="BFFscreening"),"AAL_MNI_2mm.nii")
#' mask <- oro.nifti::readNIfTI(maskfile)
#' imgfile <- file.path(system.file("nifti", package="BFFscreening"),"VBM_example.nii.gz")
#' img <- oro.nifti::readNIfTI(imgfile)
#' nb <- find_brain_image_neighbors(img, mask,radius=1)
#' @export
#'
find_brain_image_neighbors <- function(img, mask=NULL, radius = 1){


	grids <- list(X = 1:img@dim_[2],
								Y = 1:img@dim_[3],
								Z = 1:img@dim_[4])

	d = length(grids)
	dim = sapply(1:d,function(i) length(grids[[i]]))
	coords = expand.grid(grids)
	if(!is.null(mask)){
		maskidx <- which(mask>0)
		num_voxels = length(maskidx)
	} else{
		num_voxels = prod(dim)
	}


	nb_idx = seq(-radius,radius,by=1)
	nb_idx_list = list()
	for(i in 1:d){
		nb_idx_list[[i]] = nb_idx
	}
	idx_patterns = expand.grid(nb_idx_list)
	zero_idx = which(apply(idx_patterns==0,1,all))
	idx_patterns = idx_patterns[-zero_idx,]
	nb = array(NA,dim=c(num_voxels,(2*radius+1)^d))
	#nb_dist = array(0, dim=c(num_voxels,(2*radius+1)^d))
	if(!is.null(mask)){
		nb[,1] = maskidx
	} else{
		nb[,1] = 1:num_voxels
	}
	pb = txtProgressBar(style=3)
	for(l in 1:nrow(idx_patterns)){
		img_arr_idx = list()
		nb_arr_idx = list()
		for(i in 1:d){
			img_arr_idx[[i]] = max(1-idx_patterns[l,i],1):min(dim[i] - idx_patterns[l,i],dim[i])
			nb_arr_idx[[i]] = max(1+idx_patterns[l,i],1):min(dim[i] + idx_patterns[l,i],dim[i])
		}
		img_arr = expand.grid(img_arr_idx)
		img_vec_idx = 0
		for(i in d:1){
			img_vec_idx = dim[i]*img_vec_idx + (img_arr[,i]-1)
		}
		img_vec_idx = img_vec_idx + 1

		nb_arr = expand.grid(nb_arr_idx)
		nb_vec_idx = 0
		for(i in d:1){
			nb_vec_idx = dim[i]*nb_vec_idx + (nb_arr[,i]-1)
		}
		nb_vec_idx = nb_vec_idx + 1
		if(!is.null(mask)){
			nb_vec_idx_0 = intersect(nb_vec_idx,maskidx)
			subidx = match(nb_vec_idx_0,nb_vec_idx)
			nb_vec_idx = match(nb_vec_idx_0,maskidx)
			img_vec_idx = img_vec_idx[subidx]
		}
		nb[nb_vec_idx,l+1] = img_vec_idx
		#for(s in 1:d){
		#  nb_dist[,l+1] = nb_dist[,l+1] + (coords[nb[,l+1],s] - coords[,s])^2
		#}
		setTxtProgressBar(pb,l/nrow(idx_patterns))
	}

	if(!is.null(mask)){
		for(l in 2:ncol(nb)){
			nb[,l] = ifelse(is.element(nb[,l],maskidx),nb[,l],NA)
		}
	}

	close(pb)
	return(nb)
}
