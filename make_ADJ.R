make_adj = function(pv_mat, idx) {
  pv_vec = pv_mat[upper.tri(pv_mat, diag = FALSE)]
  cut_off = sort(pv_vec)[idx]
  est_adj = pv_mat <= cut_off
  diag(est_adj) = 0

  return(list(est_adj = est_adj, cut = cut_off))
}
