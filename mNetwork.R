mnetwork = function(dat, d, e) {
  dat = unname(as.matrix(dat))
  n = dim(dat)[1]
  p_est_mat = matrix(rep(0,d^2), nrow=d)
  cov_inv_est_ecov = matrix(rep(0,d^2), nrow=d)

  cb = combn(d,2)
  for (i in 1:dim(cb)[2]) {
    idx1 = cb[,i][1]
    idx2 = cb[,i][2]
    Y = dat[,idx1]
    Z = dat[,idx2]
    X = as.data.frame(dat[,-cb[,i]])

    Hat.muY = ranger(Y~., data=cbind(X,Y), num.trees = 1000)$predictions
    Hat.muZ = ranger(Z~., data=cbind(X,Z), num.trees = 1000)$predictions

    ECov2 =  mean((Y-Hat.muY)*(Z-Hat.muZ))
    EVarY2 = mean((Y-Hat.muY)^2)
    EVarZ2 = mean((Z-Hat.muZ)^2)
    onestep = ECov2/sqrt(EVarY2*EVarZ2)

    cov_inv_est_ecov[idx1,idx2] = cov_inv_est_ecov[idx2,idx1] = onestep

    Dp = (Y-Hat.muY)*(Z-Hat.muZ)/sqrt(EVarY2*EVarZ2) - ECov2/(2*sqrt(EVarY2*EVarZ2))*((Y-Hat.muY)^2/EVarY2 + (Z-Hat.muZ)^2/EVarZ2)
    sd = sqrt(mean((Dp)^2)/n)

    p_est_mat[idx1,idx2] = p_est_mat[idx2,idx1] = 2*pnorm(-abs(onestep/sd))
  }

  SEcov_adj_obj = make_adj(pv_mat = p_est_mat, idx = e)
  SEcov_adj = SEcov_adj_obj$est_adj
  SEcov_corr = cov_inv_est_ecov

  res = list()
  res$p_mat = p_est_mat
  res$adj = SEcov_adj
  res$partial_corr = SEcov_corr

  return(res)
}








