ntrait <- ncol(fito$B_hat);
if(type == "limit"){
    fct_ix <- which(fito$pve > 1/ntrait);
    s <- fito$fit;
    u <- s$u[,fct_ix, drop = FALSE];
    v <- s$v[,fct_ix, drop = FALSE];
    d <- s$d[fct_ix];
    pve <- fito$pve[fct_ix];
    fit <- list("F_hat" = v, "L_hat" = u %*% diag(d), "B_hat" = u %*% diag(d) %*% t(v), "fit" = s, "ix" = fito$ix);
}else{
    fit <- fito;
};
